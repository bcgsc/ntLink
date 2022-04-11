#!/usr/bin/env python3
"""
Use minimizers to overlap the sequence pairs that are likely overlapping
"""
import argparse
import datetime
import re
import sys
from collections import defaultdict
from collections import namedtuple
import igraph as ig
import numpy as np
from read_fasta import read_fasta


import ntjoin_utils
import ntlink_utils

MappedPathInfo = namedtuple("MappedPathInfo",
                            ["mapped_region_length", "mid_mx", "median_length_from_end"])

class ScaffoldCut:
    "Defines a scaffold, and the cut points for that scaffold"

    def __init__(self, ctg_id, sequence):
        self.ctg_id = ctg_id
        self.length = len(sequence)
        self._ori = None
        self._source_cut = None
        self._target_cut = None

    @property
    def ori(self):
        "Return orientation"
        return self._ori

    @ori.setter
    def ori(self, orientation):
        if self._ori is not None and self._ori != orientation:
            raise AssertionError("Ori is already set")
        if orientation not in ["+", "-"]:
            raise ValueError("Orientation must be + or -")
        if self._ori is None:
            if orientation == "+":
                self._target_cut, self._source_cut = 0, self.length
            elif orientation == "-":
                self._target_cut, self._source_cut = self.length, 0
        self._ori = orientation

    @property
    def source_cut(self):
        "Return source cut"
        return self._source_cut

    @source_cut.setter
    def source_cut(self, pos):
        if (self.ori == "+" and self._source_cut != self.length) or \
                (self.ori == "-" and self._source_cut != 0):
            raise AssertionError("Source cut is already set")
        self._source_cut = pos

    @property
    def target_cut(self):
        "Return target cut"
        return self._target_cut

    @target_cut.setter
    def target_cut(self, pos):
        if (self.ori == "+" and self._target_cut != 0) or \
                (self.ori == "-" and self._target_cut != self.length):
            raise AssertionError("Target cut is already set")
        self._target_cut = pos

    def get_trim_coordinates(self, k):
        if self.ori == "+":
            return self.target_cut, self.source_cut
        if self.ori == "-":
            return self.source_cut + k, self.target_cut + k
        if self.ori is None:
            return 0, self.length
        raise ValueError("Orientation should be +, - or None")


    def __str__(self):
        return f"{self.ctg_id}{self._ori} {self.length} - s:{self._source_cut} t:{self._target_cut}"

def calc_total_weight(list_files, weights):
    "Calculate the total weight of an edge given the assembly support"
    return sum([weights[f] for f in list_files])

def set_edge_attributes(graph, edge_attributes):
    "Sets the edge attributes for a python-igraph graph"
    graph.es()["support"] = [edge_attributes[e]['support'] for e in sorted(edge_attributes.keys())]
    graph.es()["weight"] = [edge_attributes[e]['weight'] for e in sorted(edge_attributes.keys())]

def is_in_valid_region(pos, valid_minimizer_positions):
    "Returns True if the minimizer is in a valid position for overlap detection, else False"
    for start, end in valid_minimizer_positions:
        if start <= pos <= end:
            return True
    return False

def read_minimizers(tsv_filename, valid_mx_positions):
    "Read the minimizers from a file, removing duplicate minimizers for a given contig"
    print(datetime.datetime.today(), ": Reading minimizers", tsv_filename, file=sys.stdout)
    mx_info = defaultdict(dict)  # contig -> mx -> (contig, position)
    mxs = {}  # Contig -> [list of minimizers]
    with open(tsv_filename, 'r') as tsv:
        for line in tsv:
            dup_mxs = set()  # Set of minimizers identified as duplicates
            line = line.strip().split("\t")
            if len(line) > 1:
                name = line[0]
                if name not in valid_mx_positions:
                    continue
                mx_pos_split = line[1].split(" ")
                for mx_pos in mx_pos_split:
                    mx, pos = mx_pos.split(":")
                    if not is_in_valid_region(int(pos), valid_mx_positions[name]):
                        continue
                    if name in mx_info and mx in mx_info[name]:  # This is a duplicate
                        dup_mxs.add(mx)
                    else:
                        mx_info[name][mx] = (name, int(pos))
                mx_info[name] = {mx: mx_info[name][mx] for mx in mx_info[name] if mx not in dup_mxs}
                mxs[name] = [[mx for mx, pos in (mx_pos.split(":") for mx_pos in mx_pos_split)
                              if mx not in dup_mxs and mx in mx_info[name] and
                              is_in_valid_region(int(pos), valid_mx_positions[name])]]

    return mx_info, mxs

def read_fasta_file_trim_prep(filename):
    "Read a fasta file into memory. Returns dictionary of scafID -> Scaffold"
    print(datetime.datetime.today(), ": Reading fasta file", filename, file=sys.stdout)
    scaffolds = {}

    with open(filename, 'r') as fasta:
        for header, seq, _, _ in read_fasta(fasta):
            scaffolds[header] = ScaffoldCut(ctg_id=header, sequence=seq)

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
                                       for assembly in list_mx_info
                                       if node['name'] in list_mx_info[assembly]])
        node_label = "\"%s\" [label=\"%s\n%s\"]" % (node['name'], node['name'], mx_ctg_pos_labels)
        outfile.write("%s\n" % node_label)

    for edge in graph.es():
        outfile.write("\"%s\" -- \"%s\"" %
                      (ntlink_utils.vertex_name(graph, edge.source),
                       ntlink_utils.vertex_name(graph, edge.target)))
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

    graph.add_vertices(list(vertices))
    graph.add_edges(formatted_edges)
    edge_attributes = {ntlink_utils.edge_index(graph, s, t): {"support": edges[s][t],
                                                  "weight": calc_total_weight(edges[s][t], weights)}
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
    to_remove_edges = [edge.index for edge in graph.es()
                       if edge['weight'] < n]
    new_graph = graph.copy()
    new_graph.delete_edges(to_remove_edges)
    return new_graph


def is_valid_pos(mx, mx_pos, start, end):
    "Return True if minimizer is within start/end coordinates, else False"
    if mx_pos[mx][1] >= start and mx_pos[mx][1] <= end:
        return True
    return False


def filter_minimizers_position(list_mxs_pair, source, target, overlap,
                               scaffolds, list_mx_info, args):
    "Filter to keep minimizers in particular positions"
    list_mxs_pair_return = {}
    source_noori, source_ori = source.strip("+-"), source[-1]
    target_noori, target_ori = target.strip("+-"), target[-1]

    start, end = ntlink_utils.find_valid_mx_region(source_noori, source_ori, scaffolds, overlap, args)
    list_mxs_pair_return[source_noori] = [[mx for mx in list_mxs_pair[source_noori][0]
                                           if is_valid_pos(mx, list_mx_info[source_noori], start, end)]]

    start, end = ntlink_utils.find_valid_mx_region(target_noori, target_ori, scaffolds, overlap, args,
                                                   source=False)
    list_mxs_pair_return[target_noori] = [[mx for mx in list_mxs_pair[target_noori][0]
                                           if is_valid_pos(mx, list_mx_info[target_noori], start, end)]]
    with ntlink_utils.HiddenPrints():
        list_mxs_pair_return = ntjoin_utils.filter_minimizers(list_mxs_pair_return)

    return list_mxs_pair_return

def set_scaffold_info(ctg_ori, pos, scaffolds, cut_type):
    "Set the cut and orientation information about the scaffold"
    ctg, orientation = ctg_ori.strip("+-"), ctg_ori[-1]
    scaffolds[ctg].ori = orientation
    if cut_type == "source":
        scaffolds[ctg].source_cut = pos
    elif cut_type == "target":
        scaffolds[ctg].target_cut = pos
    else:
        raise ValueError("cut_type must be set to source or target")

def get_dist_from_end(ori, pos, scaf_len, target=False):
    "Given the orientation, calculate the dist of the mx from the scaffold end (return -ve value)"
    if (ori == "+" and not target) or (ori == "-" and target):
        return (scaf_len - pos)*-1
    return pos*-1

def merge_overlapping(list_mxs, list_mx_info, source, target, scaffolds, args, gap):
    "Find the cut points for overlapping adjacent contigs"
    source_noori = source.strip("+-")
    target_noori = target.strip("+-")

    weights = {source_noori: 1, target_noori: 1}

    list_mxs_pair = {source_noori: list_mxs[source_noori], target_noori: list_mxs[target_noori]}
    list_mxs_pair = filter_minimizers_position(list_mxs_pair, source, target, gap, scaffolds, list_mx_info, args)
    with ntlink_utils.HiddenPrints():
        list_mxs_pair = ntjoin_utils.filter_minimizers(list_mxs_pair)

    graph = build_graph(list_mxs_pair, weights)
    graph = filter_graph_global(graph, 2)

    # Print the DOT graph if in verbose mode
    if args.v:
        print_graph(graph, list_mx_info, args.p)

    paths_components = []
    for component in graph.components():
        component_graph = graph.subgraph(component)
        source_nodes = [node.index for node in component_graph.vs() if node.degree() == 1]
        singleton_nodes = [node.index for node in component_graph.vs() if node.degree() == 0]
        if len(source_nodes) == 2:
            source_node, target_node = source_nodes
            if ntlink_utils.vertex_name(component_graph, source_node) > \
                    ntlink_utils.vertex_name(component_graph, target_node):
                source_node, target_node = target_node, source_node
            paths = component_graph.get_shortest_paths(source_node, target_node)
            assert len(paths) == 1
            path = [ntlink_utils.vertex_name(component_graph, mx) for mx in paths[0]]
            start_mx, end_mx = path[0], path[-1]
            source_start, target_start = [list_mx_info[assembly][start_mx][1]
                                          for assembly in [source_noori, target_noori]]
            source_end, target_end = [list_mx_info[assembly][end_mx][1]
                                      for assembly in [source_noori, target_noori]]
            source_align_len = abs(source_start - source_end)
            target_align_len = abs(target_start - target_end)

            mid_mx = path[int(len(path)/2)]
            mid_mx_dist_end_source = get_dist_from_end(source[-1],
                                                       list_mx_info[source_noori][mid_mx][1],
                                                       scaffolds[source_noori].length)
            mid_mx_dist_end_target = get_dist_from_end(target[-1],
                                                       list_mx_info[target_noori][mid_mx][1],
                                                       scaffolds[target_noori].length, target=True)
            paths_components.append(MappedPathInfo(mapped_region_length=np.median([source_align_len,
                                                                                   target_align_len]),
                                                   mid_mx=mid_mx,
                                                   median_length_from_end=np.median(
                                                       [mid_mx_dist_end_source, mid_mx_dist_end_target])))
        elif singleton_nodes:
            assert len(singleton_nodes) == 1
            mid_mx = ntlink_utils.vertex_name(component_graph, singleton_nodes[0])
            mid_mx_dist_end_source = get_dist_from_end(source[-1], list_mx_info[source_noori][mid_mx][1],
                                                       scaffolds[source_noori].length)
            mid_mx_dist_end_target = get_dist_from_end(target[-1], list_mx_info[target_noori][mid_mx][1],
                                                       scaffolds[target_noori].length, target=True)
            paths_components.append(MappedPathInfo(mapped_region_length=1, mid_mx=mid_mx,
                                                   median_length_from_end=np.median([mid_mx_dist_end_source,
                                                                                     mid_mx_dist_end_target])))
        else:
            print("NOTE: non-singleton, {} source nodes".format(len(source_nodes)))
    if not paths_components:
        return False
    path = sorted(paths_components, key=lambda x: (x.mapped_region_length, x.median_length_from_end,
                                                   x.mid_mx), reverse=True)[0]
    source_cut, target_cut = list_mx_info[source_noori][path.mid_mx][1], list_mx_info[target_noori][path.mid_mx][1]

    if source_cut is None or target_cut is None:
        return False

    set_scaffold_info(source, source_cut, scaffolds, "source")
    set_scaffold_info(target, target_cut, scaffolds, "target")

    return True

def merge_overlapping_pathfile(args, gap_re, graph, mxs, mxs_info, scaffolds):
    "Read through pathfile, and merge overlapping pieces, updating path file"
    print(datetime.datetime.today(), ": Finding scaffold overlaps", file=sys.stdout)
    out_pathfile = open(args.p + ".trimmed_scafs.path", 'w')
    with open(args.path, 'r') as path_fin:
        for path in path_fin:
            new_path = []
            path_id, path_seq = path.strip().split("\t")
            path_seq = path_seq.split(" ")
            path_seq = ntlink_utils.normalize_path(path_seq, gap_re)
            for source, gap, target in zip(path_seq, path_seq[1:], path_seq[2:]):
                gap_match = re.search(gap_re, gap)
                if not gap_match:
                    continue
                if int(gap_match.group(1)) <= args.g + 1 and \
                        ntlink_utils.has_estimated_overlap(graph, source, target):
                    cuts_found = merge_overlapping(mxs, mxs_info, source, target, scaffolds, args,
                                      graph.es()[ntlink_utils.edge_index(graph, source, target)]["d"])
                    if cuts_found:
                        gap = "{}N".format(args.outgap)
                if not new_path:
                    new_path.append(source)
                new_path.append(gap)
                new_path.append(target)
            out_pathfile.write("{path_id}\t{ctgs}\n".format(path_id=path_id, ctgs=" ".join(new_path)))
    out_pathfile.close()

def print_trimmed_scaffolds(args, scaffolds):
    "Print the trimmed scaffolds fasta to file"
    print(datetime.datetime.today(), ": Printing trimmed scaffolds", file=sys.stdout)
    fasta_outfile = open(args.p + ".trimmed_scafs.fa", 'w')
    with open(args.fasta, 'r') as fin:
        for name, seq, _, _ in read_fasta(fin):
            scaffold = scaffolds[name]
            if scaffold.ori == "+":
                sequence_out = seq[scaffold.target_cut:scaffold.source_cut]
            elif scaffold.ori == "-":
                sequence_out = seq[scaffold.source_cut + args.k:scaffold.target_cut + args.k]
            elif scaffold.ori is None:
                sequence_out = seq
            else:
                raise ValueError("Invalid orientation for Scaffold:", scaffold)
            if len(sequence_out) == 0:
                sequence_out = "N"
            fasta_outfile.write(
                ">{} {}-{}\n{}\n".format(scaffold.ctg_id, scaffold.source_cut, scaffold.target_cut, sequence_out))
    fasta_outfile.close()

def print_trim_coordinates(args, scaffolds):
    "Print coordinates of trimming done for scaffolds"
    with open(args.p + ".trimmed_scafs.tsv", 'w') as tsvfile:
        for scaffold in scaffolds:
            scaffold_entry = scaffolds[scaffold]
            start, end = scaffold_entry.get_trim_coordinates(args.k)
            out_str = "{}\t{}\t{}\n".format(scaffold_entry.ctg_id, start, end, sep="\t")
            tsvfile.write(out_str)

def parse_arguments():
    "Parse arguments for ntLink overlap"
    parser = argparse.ArgumentParser(description="Find coordinates for combining overlapping sequences")
    parser.add_argument("-m", help="Minimizer TSV file", type=str, required=True)
    parser.add_argument("-f", help="Fudge factor for estimated overlap [0.5]", type=float, default=0.5)
    parser.add_argument("-a", "--path", help="Path file", required=True, type=str)
    parser.add_argument("-s", "--fasta", help="Scaffold sequences", required=True, type=str)
    parser.add_argument("-k", help="Indexlr k", required=True, type=int)
    parser.add_argument("-d", help="Scaffold dot file", required=True, type=str)
    parser.add_argument("-g", help="Minimum gap size (bp) [20]", default=20, type=int)
    parser.add_argument("--outgap", help="Gap size between trimmed overlapping sequences (bp) [2]", default=2, type=int)
    parser.add_argument("-p", help="Output file prefix [ntlink_merge]", default="ntlink_merge", type=str)
    parser.add_argument("-v", help="Verbose output logging", action="store_true")
    parser.add_argument("--trim_info", help="Verbose log of trimming info", action="store_true")
    parser.add_argument("--version", action='version', version='ntLink v1.2.0')

    return parser.parse_args()

def print_args(args):
    "Print the parameters for the ntLink overlap step"
    print("Parameters for overlap stage:")
    print("\t-m", args.m)
    print("\t-f", args.f)
    print("\t-a", args.path)
    print("\t-s", args.fasta)
    print("\t-k", args.k)
    print("\t-d", args.d)
    print("\t-g", args.g)
    print("\t--outgap", args.outgap)
    print("\t-p", args.p)

def main():
    "Run overlap sequence detection and trimming step of ntLink"
    print("Assessing putative overlaps...")

    args = parse_arguments()
    print_args(args)

    gap_re = re.compile(r'^(\d+)N$')
    args.outgap = args.outgap + 1 # Accounting for abyss-scaffold incrementing by one in path file

    scaffolds = read_fasta_file_trim_prep(args.fasta)
    graph, _ = ntlink_utils.read_scaffold_graph(args.d)

    valid_mx_positions = ntlink_utils.find_valid_mx_regions(args, gap_re, graph, scaffolds)

    args.m = "/dev/stdin" if args.m == "-" else args.m
    mxs_info, mxs = read_minimizers(args.m, valid_mx_positions)

    merge_overlapping_pathfile(args, gap_re, graph, mxs, mxs_info, scaffolds)

    if args.trim_info:
        print_trim_coordinates(args, scaffolds)

    print_trimmed_scaffolds(args, scaffolds)

    print(datetime.datetime.today(), ": DONE!", file=sys.stdout)


if __name__ == '__main__':
    main()
