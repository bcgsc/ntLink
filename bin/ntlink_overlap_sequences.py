#!/usr/bin/env python3
import argparse
import datetime
import igraph as ig
import numpy as np
import re
from collections import defaultdict
import sys

from ntjoin_assemble import *
from ntlink_stitch_paths import NtLinkPath

'''
Use minimizers to overlap the sequence pairs that are likely overlapping
'''

class Scaffold:
    "Defines a scaffold, and the cut points for that scaffold"

    def __init__(self, ctg_id, sequence):
        self.ctg_id = ctg_id
        self.sequence = sequence
        self.length = len(sequence)
        self._ori = None
        self._source_cut = len(sequence)
        self._target_cut = 0

    @property
    def ori(self, orientation):
        if self._ori is not None and self._ori != orientation:
            raise AssertionError("Ori is already set")
        if orientation not in ["+", "-"]:
            raise ValueError("Orientation must be + or -")
        self._ori = orientation

    @ori.setter
    def ori(self):
        return self._ori

    @property
    def source_cut(self, pos):
        if self._source_cut != len(self.sequence):
            raise AssertionError("Source cut is already set")
        self._source_cut = pos

    @source_cut.setter
    def source_cut(self):
        return self._source_cut

    @property
    def target_cut(self, pos):
        if self._target_cut != 0:
            raise AssertionError("Target cut is already set")
        self._target_cut = pos

    @target_cut.setter
    def target_cut(self):
        return self._target_cut


def edge_index(graph, source_name, target_name):
    "Returns graph edge index based on source/target names"
    return graph.get_eid(source_name, target_name)

def calc_total_weight(list_files, weights):
    "Calculate the total weight of an edge given the assembly support"
    return sum([weights[f] for f in list_files])

def set_edge_attributes(graph, edge_attributes):
    "Sets the edge attributes for a python-igraph graph"
    graph.es()["support"] = [edge_attributes[e]['support'] for e in sorted(edge_attributes.keys())]
    graph.es()["weight"] = [edge_attributes[e]['weight'] for e in sorted(edge_attributes.keys())]

def vertex_name(graph, index):
    "Returns vertex name based on vertex id"
    return graph.vs[index]['name']

def read_minimizers(tsv_filename):
    "Read the minimizers from a file, removing duplicate minimizers for a given contig"
    print(datetime.datetime.today(), ": Reading minimizers", tsv_filename, file=sys.stdout)
    mx_info = defaultdict(dict)  # contig -> mx -> (contig, position)
    mx_info_filt = defaultdict(dict)
    mxs = {}  # Contig -> [list of minimizers]
    with open(tsv_filename, 'r') as tsv:
        for line in tsv:
            dup_mxs = set()  # Set of minimizers identified as duplicates
            line = line.strip().split("\t")
            if len(line) > 1:
                name = line[0]
                mx_pos_split = line[1].split(" ")
                for mx_pos in mx_pos_split:
                    mx, pos = mx_pos.split(":")
                    if name in mx_info and mx in mx_info[name]:  # This is a duplicate, add to dup set, don't add to dict
                        dup_mxs.add(mx)
                    else:
                        mx_info[name][mx] = (name, int(pos))
                mxs[name] = [[mx_pos.split(":")[0] for mx_pos in mx_pos_split if mx_pos.split(":")[0] not in dup_mxs]]
                mx_info_filt[name] = {}
                mx_info_filt[name] = {mx: mx_info[name][mx] for mx in mx_info[name] if mx not in dup_mxs}

    return mx_info_filt, mxs

def read_fasta_file(filename):
    "Read a fasta file into memory. Returns dictionary of scafID -> Scaffold"
    print(datetime.datetime.today(), ": Reading fasta file", filename, file=sys.stdout)
    scaffolds = {}

    with open(filename, 'r') as fasta:
        for header, seq, _, _ in read_fasta(fasta):
            scaffolds[header] = Scaffold(ctg_id=header, sequence=seq)

    return scaffolds

def print_graph(graph, list_mx_info, prefix):
    "Prints the minimizer graph in dot format"
    out_graph = prefix + ".mx.dot"
    outfile = open(out_graph, 'a')
    print(datetime.datetime.today(), ": Printing graph", out_graph, sep=" ", file=sys.stdout)

    outfile.write("graph G {\n")

    colours = ["red", "green", "blue", "purple", "orange",
               "turquoise", "pink", "yellow", "orchid", "salmon"]
    list_files = list(list_mx_info.keys())
    if len(list_files) > len(colours):
        colours = ["red"]*len(list_files)

    for node in graph.vs():
        mx_ctg_pos_labels = "\n".join([str(list_mx_info[assembly][node['name']])
                                       for assembly in list_mx_info if node['name'] in list_mx_info[assembly]])
        node_label = "\"%s\" [label=\"%s\n%s\"]" % (node['name'], node['name'], mx_ctg_pos_labels)
        outfile.write("%s\n" % node_label)

    for edge in graph.es():
        outfile.write("\"%s\" -- \"%s\"" %
                      (vertex_name(graph, edge.source),
                       vertex_name(graph, edge.target)))
        weight = edge['weight']
        support = edge['support']
        if len(support) == 1:
            colour = colours[list_files.index(support[0])]
        elif len(support) == 2:
            colour = "lightgrey"
        else:
            colour = "black"
        outfile.write(" [weight=%s color=%s]\n" % (weight, colour))

    outfile.write("}\n")

    print("\nfile_name\tnumber\tcolour")
    for i, filename in enumerate(list_files):
        print(filename, i, colours[i], sep="\t")
    print("")

def build_graph(list_mxs, weights):
    "Builds an undirected graph: nodes=minimizers; edges=between adjacent minimizers"
    print(datetime.datetime.today(), ": Building graph", file=sys.stdout)
    graph = ig.Graph()

    vertices = set()
    edges = defaultdict(dict)  # source -> target -> [list assembly support]

    for assembly in list_mxs:
        for assembly_mx_list in list_mxs[assembly]:
            for i, j in zip(range(0, len(assembly_mx_list)),
                            range(1, len(assembly_mx_list))):
                if assembly_mx_list[i] in edges and \
                        assembly_mx_list[j] in edges[assembly_mx_list[i]]:
                    edges[assembly_mx_list[i]][assembly_mx_list[j]].append(assembly)
                elif assembly_mx_list[j] in edges and \
                        assembly_mx_list[i] in edges[assembly_mx_list[j]]:
                    edges[assembly_mx_list[j]][assembly_mx_list[i]].append(assembly)
                else:
                    edges[assembly_mx_list[i]][assembly_mx_list[j]] = [assembly]
                vertices.add(assembly_mx_list[i])
            if assembly_mx_list:
                vertices.add(assembly_mx_list[-1])

    formatted_edges = [(s, t) for s in edges for t in edges[s]]

    print(datetime.datetime.today(), ": Adding vertices", file=sys.stdout)
    graph.add_vertices(list(vertices))

    print(datetime.datetime.today(), ": Adding edges", file=sys.stdout)
    graph.add_edges(formatted_edges)

    print(datetime.datetime.today(), ": Adding attributes", file=sys.stdout)
    edge_attributes = {edge_index(graph, s, t): {"support": edges[s][t],
                                                      "weight": calc_total_weight(edges[s][t],
                                                                                       weights)}
                       for s in edges for t in edges[s]}
    set_edge_attributes(graph, edge_attributes)

    return graph

def is_graph_linear(graph):
    "Given a graph, return True if all the components are linear"
    for component in graph.components():
        component_graph = graph.subgraph(component)
        if not all(u.degree() < 3 for u in component_graph.vs()):
            return False
    return True

def filter_graph_global(graph, n):
    "Filter the graph globally based on minimum edge weight"
    print(datetime.datetime.today(), ": Filtering the graph", file=sys.stdout)
    to_remove_edges = [edge.index for edge in graph.es()
                       if edge['weight'] < n]
    new_graph = graph.copy()
    new_graph.delete_edges(to_remove_edges)
    return new_graph


def is_valid_pos(mx, mx_pos, start, end):
    if mx_pos[mx][1] >= start and mx_pos[mx][1] <= end:
        return True
    return False


def filter_minimizers_position(list_mxs_pair, source, target, overlap, scaffolds, list_mx_info, fudge_factor):
    "Filter to keep minimizers in particular positions"
    list_mxs_pair_return = {}
    source_noori, source_ori = source.strip("+-"), source[-1]
    target_noori, target_ori = target.strip("+-"), target[-1]

    if source_ori == "+":
        start, end = (scaffolds[source_noori].length - overlap*-1 - 15) - int(overlap*-1*fudge_factor), \
                     scaffolds[source_noori].length
    else:
        start, end = 0, int(overlap*-1*(fudge_factor+1))
    print(source, start, end)
    print(list_mxs_pair[source_noori])
    list_mxs_pair_return[source_noori] = [[mx for mx in list_mxs_pair[source_noori][0]
                                     if is_valid_pos(mx, list_mx_info[source_noori], start, end)]]

    if target_ori == "-":
        start, end = scaffolds[target_noori].length - overlap*-1 - 15 - int(overlap*-1*fudge_factor), \
                     scaffolds[target_noori].length
    else:
        start, end = 0, int(overlap*-1*(fudge_factor+1))
    print(target, start, end)
    list_mxs_pair_return[target_noori] = [[mx for mx in list_mxs_pair[target_noori][0]
                                     if is_valid_pos(mx, list_mx_info[target_noori], start, end)]]

    list_mxs_pair_return = Ntjoin.filter_minimizers(list_mxs_pair_return)

    return list_mxs_pair_return

def set_scaffold_info(ctg_ori, pos, scaffolds, cut_type):
    "Set the cut and orientation information about the scaffold"
    ctg = ctg_ori.strip("+-")
    ori = ctg_ori[-1]
    scaffolds[ctg].ori = ori
    if cut_type == "source":
        scaffolds[ctg].source_cut = pos
    elif cut_type == "target":
        scaffolds[ctg].target_cut = pos
    else:
        raise ValueError("cut_type must be set to source or target")

def merge_overlapping(list_mxs, list_mx_info, source, target, gap, scaffolds, args):
    source_noori = source.strip("+-")
    target_noori = target.strip("+-")

    weights = {source_noori: 1, target_noori: 1}

    list_mxs_pair = {name: list_mxs[name] for name in list_mxs if (name == source_noori or name == target_noori)}

    list_mxs_pair = Ntjoin.filter_minimizers(list_mxs_pair)

    list_mxs_pair = filter_minimizers_position(list_mxs_pair, source, target, gap, scaffolds, list_mx_info, args.f)

    graph = build_graph(list_mxs_pair, weights)

    graph = filter_graph_global(graph, 2)

    # Print the DOT graph
    if args.v:
        print_graph(graph, list_mx_info, args.p)

    paths_components = []
    for component in graph.components():
        component_graph = graph.subgraph(component)
        source_nodes = [node.index for node in component_graph.vs() if node.degree() == 1]
        singleton_node = [node.index for node in component_graph.vs() if node.degree() == 0]
        if len(source_nodes) == 2:
            source_comp, target_comp = source_nodes
            if vertex_name(component_graph, source_comp) > vertex_name(component_graph, target_comp):
                source_tmp = source_comp
                source_comp = target_comp
                target_comp = source_tmp
            paths = component_graph.get_shortest_paths(source_comp, target_comp)

            if len(paths) > 1:
                print("NOTE: more than one path found")
            path = paths[0]
            source_start, target_start = [list_mx_info[assembly][vertex_name(component_graph, path[0])][1] for assembly in [source_noori, target_noori]]
            source_end, target_end = [list_mx_info[assembly][vertex_name(component_graph, path[-1])][1] for assembly in [source_noori, target_noori]]
            source_align_len = abs(source_start - source_end)
            target_align_len = abs(target_start - target_end)

            path = [vertex_name(component_graph, mx) for mx in path]
            paths_components.append((path, np.median([source_align_len, target_align_len])))
        elif singleton_node:
            paths_components.append((vertex_name(component_graph, singleton_node), 1))
    if not paths_components:
        return
    print(sorted(paths_components, key=lambda x: x[1], reverse=True))
    path = sorted(paths_components, key=lambda x: x[1], reverse=True)[0][0]
    mx = path[int(len(path)/2)]
    cuts = {list_mx_info[assembly][mx][0]: list_mx_info[assembly][mx][1] for assembly in [source_noori, target_noori]}

    if not cuts:
        return
    source_cut = cuts[source.strip("+-")]
    target_cut = cuts[target.strip("+-")]

    set_scaffold_info(source, source_cut, scaffolds, "source")
    set_scaffold_info(target, target_cut, scaffolds, "target")


def main():
    print("Assessing putative overlaps...")

    parser = argparse.ArgumentParser(description="Find coordinates for combining overlapping sequences")
    parser.add_argument("-m", help="Minimizer TSV file", type=str, required=True)
    parser.add_argument("-f", help="Fudge factor for estimated overlap [0.5]", type=float, default=0.5)
    parser.add_argument("-a", help="Path file", required=True, type=str)
    parser.add_argument("-s", help="Scaffold sequences", required=True, type=str)
    parser.add_argument("-d", help="Scaffold dot file", required=True, type=str)
    parser.add_argument("-g", help="Minimum gap size (bp) [20]", default=20, type=int)
    parser.add_argument("-p", help="Output file prefix [ntlink_merge]", default="ntlink_merge", type=str)
    parser.add_argument("-v", help="Verbose output logging", action="store_true")

    args = parser.parse_args()

    scaffolds = Ntjoin.read_fasta_file(args.s)
    graph = NtLinkPath.read_scaffold_graph(args.d)

    # !! TODO only load minimizers into file that are useful

    mxs_info, mxs = read_minimizers(args.m)

    gap_re = re.compile(r'^(\d+)N$')

    out_pathfile = open(args.p + ".path", 'w')

    with open(args.a, 'r') as path_fin:
        for path in path_fin:
            new_path = []
            path_id, path_seq = path.strip().split("\t")
            path_seq = path_seq.split(" ")
            for i, j, k in zip(path_seq, path_seq[1:], path_seq[2:]):
                source, gap, target = i, j, k
                gap_match = re.search(gap_re, gap)
                if not gap_match:
                    continue
                if int(gap_match.group(1)) <= args.g + 1 and graph.es()[edge_index(graph, source, target)]["d"]  < 0:
                    gap = graph.es()[edge_index(graph, source, target)]["d"]
                    merge_overlapping(mxs, mxs_info, source, target, gap, scaffolds, args) #!! TODO: output file name
                    gap = "{}N".format(2) #!! TODO: make parameter
                if not new_path:
                    new_path.append(source)
                new_path.append(gap)
                new_path.append(target)
            out_pathfile.write("{path_id}\t{ctgs}\n".format(path_id=path_id, ctgs=" ".join(new_path)))
    out_pathfile.close()

    # Print out all scaffolds
    fasta_outfile = open(args.p + ".trimmed_scafs.fa", 'w')
    for out_scaffold in scaffolds:
        scaffold = scaffolds[out_scaffold]
        if scaffold.ori == "+":
            sequence = scaffold.sequence[scaffold.target_cut:scaffold.source_cut]
        else:
            sequence = scaffold.sequence[scaffold.source_cut + 15:scaffold.target_cut+15] # !! TODO: make k parameter
        if len(sequence) == 0:
            sequence = "N"
            print("HERE")
        fasta_outfile.write(">{}\n{}\n".format(scaffold.id, sequence))
    fasta_outfile.close()

    # if source[-1] == "+":
    #     source_piece = scaffolds[source.strip("+-")].sequence[:source_cut]
    # else:
    #     source_piece = scaffolds[source.strip("+-")].sequence[source_cut + 15:]
    #
    # if target[-1] == "+":
    #     target_piece = scaffolds[target.strip("+-")].sequence[target_cut:]
    # else:
    #     target_piece = scaffolds[target.strip("+-")].sequence[:target_cut + 15]
    #
    # print(source[:-1], source_piece)
    # print(target[:-1], target_piece)
    #
    # # if len(source_piece) == 0:
    # #     source_piece = "N"
    # #     print("HERE")
    # # if len(target_piece) == 0:
    # #     target_piece = "N"
    # #     print("HERE")
    # scaffolds[source.strip("+-")] = Scaffold(id=source.strip("+-"), sequence=source_piece, length=len(source_piece))
    # scaffolds[target.strip("+-")] = Scaffold(id=target.strip("+-"), sequence=target_piece, length=len(target_piece))



if __name__ == '__main__':
    main()

