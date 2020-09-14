#!/usr/bin/env python3
"""
Using ntJoin concept with long reads as input reference
"""
__author__ = 'laurencoombe'

import argparse
import datetime
import multiprocessing
import re
from collections import defaultdict
from collections import namedtuple
import sys
import igraph as ig
from ntjoin_assemble import PathNode, Ntjoin

__author__ = 'laurencoombe'

EdgeInfo = namedtuple("EdgeInfo", ["gap_estimate", "num_links"])

class NtjoinLongScaffold():
    "Scaffold layout using a directed scaffold graph"

    def read_scaffold_graph(self, in_graph_file):
        "Reads in a scaffold graph in dot format"

        graph = ig.Graph(directed=True)

        vertices = set()
        edges = defaultdict(dict) # source -> target -> EdgeInfo

        node_re = re.compile(r'\"(\S+[+-])\"\s+\[l\=\d+\]')
        edge_re = re.compile(r'\"(\S+[+-])\"\s+\-\>\s+\"(\S+[+-])\"\s+\[d\=(\-?\d+)\s+e\=\d+\s+n\=(\d+)\]')

        past_header = False

        with open(in_graph_file, 'r') as in_graph:
            for line in in_graph:
                line = line.strip()
                if not past_header:
                    past_header = True
                    continue
                node_match = re.search(node_re, line)
                if node_match:
                    vertices.add(node_match.group(1))
                    continue

                edge_match = re.search(edge_re, line)
                if edge_match:
                    source, target, gap_est, num_links = edge_match.group(1), edge_match.group(2), \
                                                         edge_match.group(3), edge_match.group(4)
                    if NtjoinLongScaffold.scaffolds[source[:-1]].length < self.args.z or \
                        NtjoinLongScaffold.scaffolds[target[:-1]].length < self.args.z:
                        continue
                    edges[source][target] = EdgeInfo(int(gap_est), int(num_links))
                elif line != "}":
                    print("Error! Unexpected line in input dot file:", line)
                    sys.exit(1)

        formatted_edges = [(s, t) for s in edges for t in edges[s]]
        graph.add_vertices(list(vertices))
        graph.add_edges(formatted_edges)

        edge_attributes = {Ntjoin.edge_index(graph, s, t): {'d': edges[s][t].gap_estimate,
                                                            "n": edges[s][t].num_links}
                           for s in edges for t in edges[s]}
        graph.es()["d"] = [edge_attributes[e]['d'] for e in sorted(edge_attributes.keys())]
        graph.es()["n"] = [edge_attributes[e]['n'] for e in sorted(edge_attributes.keys())]

        max_weight = max([edge_attributes[e]['n'] for e in edge_attributes.keys()])

        return graph, max_weight

    @staticmethod
    def is_graph_linear(graph):
        "Given a graph, return True if all the components are linear"
        for component in graph.components(mode="weak"):
            component_graph = graph.subgraph(component)
            if not all(u.indegree() < 2 for u in component_graph.vs()):
                return False
            if not all(u.outdegree() < 2 for u in component_graph.vs()):
                return False
        return True

    @staticmethod
    def filter_graph(graph, min_weight):
        "Filter the graph by edge weights on edges incident to branch nodes"
        branch_nodes_in = [node.index for node in graph.vs() if node.indegree() > 1]
        to_remove_edges_in = [edge for node in branch_nodes_in for edge in graph.incident(node, mode="in")
                              if graph.es()[edge]['n'] < min_weight]

        branch_nodes_out = [node.index for node in graph.vs() if node.outdegree() > 1]
        to_remove_edges_out = [edge for node in branch_nodes_out for edge in graph.incident(node, mode="out")
                               if graph.es()[edge]['n'] < min_weight]
        to_remove_edges = to_remove_edges_in + to_remove_edges_out

        new_graph = graph.copy()
        new_graph.delete_edges(to_remove_edges)
        return new_graph

    def format_path_contigs(self, path, component_graph):
        "Given a path (sequence of oriented contigs), format to a path of PathNode"
        return_path = []
        for ctga, ctgb in zip(path, path[1:]):
            ctga_name, ctga_ori = ctga[:-1], ctga[-1]
            edge_index = Ntjoin.edge_index(component_graph, ctga, ctgb)
            gap_estimate = component_graph.es()[edge_index]['d']
            return_path.append(PathNode(contig=ctga_name, ori=ctga_ori,
                                        start=0, end=NtjoinLongScaffold.scaffolds[ctga_name].length,
                                        contig_size=NtjoinLongScaffold.scaffolds[ctga_name].length,
                                        first_mx=None, terminal_mx=None,
                                        gap_size=max(gap_estimate, self.args.g)))
        last_ctg_name, last_ctg_ori = path[-1][:-1], path[-1][-1]
        return_path.append(PathNode(contig=last_ctg_name, ori=last_ctg_ori,
                                    start=0, end=NtjoinLongScaffold.scaffolds[last_ctg_name].length,
                                    contig_size=NtjoinLongScaffold.scaffolds[last_ctg_name].length,
                                    first_mx=None, terminal_mx=None,
                                    gap_size=0))
        return return_path


    def find_paths_process(self, component):
        "Find paths given a component of the graph"
        return_paths = []
        min_edge_weight = self.args.n
        max_edge_weight = NtjoinLongScaffold.max_weight
        component_graph = NtjoinLongScaffold.gin.subgraph(component)
        while not self.is_graph_linear(component_graph) and \
                min_edge_weight <= max_edge_weight:
            component_graph = self.filter_graph(component_graph, min_edge_weight)
            min_edge_weight += 1

        visited = set()

        for subcomponent in component_graph.components(mode="weak"):
            subcomponent_graph = component_graph.subgraph(subcomponent)
            source_nodes = [node.index for node in subcomponent_graph.vs() if node.indegree() == 0]
            if len(source_nodes) == 1:
                target_nodes = [node.index for node in subcomponent_graph.vs() if node.outdegree() == 0]
                assert len(target_nodes) == 1
                source, target = source_nodes[0], target_nodes[0]
                path = subcomponent_graph.get_shortest_paths(source, target)[0]
                if any([Ntjoin.vertex_name(subcomponent_graph, node)[:-1] in visited for node in path]):
                    continue
                num_edges = len(path) - 1
                if len(path) == len(subcomponent_graph.vs()) and \
                        num_edges == len(subcomponent_graph.es()) and len(path) == len(set(path)):
                    # All the nodes/edges from the graph are in the simple path, no repeated nodes
                    path = Ntjoin.convert_path_index_to_name(subcomponent_graph, path)
                    ctg_path = self.format_path_contigs(path, subcomponent_graph)
                    for ctg in ctg_path:
                        visited.add(ctg.contig)
                    return_paths.append(ctg_path)
            else:
                print(subcomponent_graph)
        return return_paths

    @staticmethod
    def remove_duplicate_paths(paths):
        "Removes paths that are reverse complements of each other or have dup contigs"
        visited = set()
        return_paths = []
        for path in paths:
            if not any([node.contig in visited for node in path]):
                return_paths.append(path)
            for node in path:
                visited.add(node.contig)
        return return_paths


    def find_paths(self, graph):
        "Finds paths through input scaffold graph"
        print(datetime.datetime.today(), ": Finding paths", file=sys.stdout)
        NtjoinLongScaffold.gin = graph
        components = graph.components(mode="weak")
        print("\nTotal number of components in graph:", len(components), "\n", sep=" ", file=sys.stdout)

        if self.args.t == 1:
            paths = [self.find_paths_process(component) for component in components]
        else:
            with multiprocessing.Pool(self.args.t) as pool:
                paths = pool.map(self.find_paths_process, components)
        paths_return = [path for path_list in paths for path in path_list]
        paths_return = self.remove_duplicate_paths(paths_return)
        return paths_return


    @staticmethod
    def parse_arguments():
        "Parse ntJoin arguments"
        parser = argparse.ArgumentParser(
            description="ntJoin: Scaffolding genome assemblies using reference assemblies and minimizer graphs",
            epilog="Note: Script expects that each input minimizer TSV file has a matching fasta file.\n"
                   "Example: myscaffolds.fa.k32.w1000.tsv - myscaffolds.fa is the expected matching fasta",
            formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("FASTA", nargs="+", help="Assembly to scaffold")
        parser.add_argument("-d", help="Directed graph (dot format)", required=True, type=str)
        parser.add_argument("-p", help="Output prefix [out]", default="out",
                            type=str, required=False)
        parser.add_argument("-n", help="Minimum edge weight [1]", default=1, type=int)
        parser.add_argument("-g", help="Minimum gap size (bp)", required=False, default=20, type=int)
        parser.add_argument("-z", help="Minimum size of contig to scaffold", required=False, default=500, type=int)
        parser.add_argument('-t', help="Number of threads [1]", default=1, type=int)
        parser.add_argument("-v", "--version", action='version', version='ntJoin v1.0.2')
        parser.add_argument("-k", help="Kmer size used for minimizer step", required=True, type=int)
        parser.add_argument("-w", help="NtJoin w value", required=True, type=int)
        parser.add_argument("--agp", help="Output AGP file describing scaffolds", action="store_true")
        return parser.parse_args()


    def print_parameters(self):
        "Print the set parameters for the ntJoin run"
        print("Assembly to be scaffolded:", self.args.FASTA)
        print("Directed graph:", self.args.d)
        print("Parameters:")
        print("\t-p ", self.args.p)
        print("\t-n ", self.args.n)
        print("\t-g ", self.args.g)
        print("\t-z ", self.args.z)
        print("\t-t ", self.args.t)
        if self.args.agp:
            print("\t--agp")

    def main(self):
        "Run ntJoin graph stage"
        print("Running ntJoin long scaffolding stage ...\n")
        self.print_parameters()

        # Load target scaffolds into memory
        if len(self.args.FASTA) != 1:
            print("ERROR: Only one input fasta file is expected,", len(self.args.FASTA), "found. Exiting...")
            sys.exit(1)
        scaffolds = Ntjoin.read_fasta_file(self.args.FASTA[0])  # scaffold_id -> Scaffold
        NtjoinLongScaffold.scaffolds = scaffolds
        Ntjoin.scaffolds = scaffolds

        # Read in the directed graph
        graph, max_weight = self.read_scaffold_graph(self.args.d)
        NtjoinLongScaffold.max_weight = max_weight

        # Find the paths through the graph
        paths = self.find_paths(graph)

        # Print the final scaffolds
        Ntjoin.print_scaffolds(self, paths, self.args.FASTA[0], "k{k}.w{w}.n{n}.z{z}".format(k=self.args.k,
                                                                                             w=self.args.w,
                                                                                             n=self.args.n,
                                                                                             z=self.args.z))

        print(datetime.datetime.today(), ": DONE!", file=sys.stdout)

    def __init__(self):
        "Create an ntJoin instance"
        self.args = self.parse_arguments()

def main():
    "Run ntJoin"
    NtjoinLongScaffold().main()

if __name__ == "__main__":
    main()
