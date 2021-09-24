#!/usr/bin/env python3
"""
Use minimizers to overlap the sequence pairs that are likely overlapping
"""
import argparse
import datetime
import re
import sys
import os
from collections import defaultdict
from collections import namedtuple
import igraph as ig
import numpy as np

import ntjoin_utils
from ntlink_stitch_paths import NtLinkPath
import ntlink_utils
from read_fasta import read_fasta

MappedPathInfo = namedtuple("MappedPathInfo",
                            ["mapped_region_length", "mid_mx", "median_length_from_end"])

class Scaffold:
    "Defines a scaffold, and the cut points for that scaffold"

    def __init__(self, ctg_id, sequence):
        self.ctg_id = ctg_id
        self.sequence = sequence
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

    def __str__(self):
        return f"{self.ctg_id}{self._ori} {self.length} - s:{self._source_cut} t:{self._target_cut}"

class HiddenPrints:
    "Adapted from: https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print"
    def __init__(self):
        self._original_stdout = sys.stdout

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

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
                                       for assembly in list_mx_info
                                       if node['name'] in list_mx_info[assembly]])
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
    edge_attributes = {edge_index(graph, s, t): {"support": edges[s][t],
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

    start, end = find_valid_mx_region(source_noori, source_ori, scaffolds, overlap, args)
    list_mxs_pair_return[source_noori] = [[mx for mx in list_mxs_pair[source_noori][0]
                                     if is_valid_pos(mx, list_mx_info[source_noori], start, end)]]

    start, end = find_valid_mx_region(target_noori, target_ori, scaffolds, overlap, args,
                                      source=False)
    list_mxs_pair_return[target_noori] = [[mx for mx in list_mxs_pair[target_noori][0]
                                     if is_valid_pos(mx, list_mx_info[target_noori], start, end)]]
    with HiddenPrints():
        list_mxs_pair_return = ntjoin_utils.filter_minimizers(list_mxs_pair_return)

    return list_mxs_pair_return


def find_valid_mx_region(scaf_noori, scaf_ori, scaffolds, overlap, args, source=True):
    "Return start/end of valid minimizer region on the scaffold"
    if (scaf_ori == "+" and source) or (scaf_ori == "-" and not source):
        start, end = (scaffolds[scaf_noori].length - overlap * -1 - args.k) - \
                     int(overlap * -1 * args.f), scaffolds[scaf_noori].length
    else:
        start, end = 0, int(overlap * -1 * (args.f + 1))

    return start, end


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
    with HiddenPrints():
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
            if vertex_name(component_graph, source_node) > \
                    vertex_name(component_graph, target_node):
                source_node, target_node = target_node, source_node
            paths = component_graph.get_shortest_paths(source_node, target_node)
            assert len(paths) == 1
            path = [vertex_name(component_graph, mx) for mx in paths[0]]
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
            mid_mx = vertex_name(component_graph, singleton_nodes[0])
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

def normalize_path(path_sequence, gap_re):
    "Given a path, normalize it to ensure deterministic running"
    if path_sequence[0].strip("+-") < path_sequence[-1].strip("+-"):
        return path_sequence
    new_seq = []
    for node in reversed(path_sequence):
        if re.search(gap_re, node):
            new_seq.append(node)
        else:
            new_seq.append(ntlink_utils.reverse_scaf_ori(node))
    return new_seq

def find_valid_mx_regions(args, gap_re, graph, scaffolds):
    "Return a dictionary with scaffold -> [(start, end)]"
    valid_regions = {}

    with open(args.a, 'r') as path_fin:
        for path in path_fin:
            _, path_seq = path.strip().split("\t")
            path_seq = path_seq.split(" ")
            path_seq = normalize_path(path_seq, gap_re)
            for source, gap, target in zip(path_seq, path_seq[1:], path_seq[2:]):
                source_noori, target_noori = source.strip("+-"), target.strip("+-")
                gap_match = re.search(gap_re, gap)
                if not gap_match:
                    continue
                if int(gap_match.group(1)) <= args.g + 1 and graph.es()[edge_index(graph, source, target)]["d"] < 0:
                    gap = graph.es()[edge_index(graph, source, target)]["d"]
                    source_start, source_end = find_valid_mx_region(source_noori, source[-1],
                                                                    scaffolds, gap, args)
                    if source_noori not in valid_regions:
                        valid_regions[source_noori] = []
                    valid_regions[source_noori].append((source_start, source_end))

                    target_start, target_end = find_valid_mx_region(target_noori, target[-1],
                                                                    scaffolds, gap, args, source=False)
                    if target_noori not in valid_regions:
                        valid_regions[target_noori] = []
                    valid_regions[target_noori].append((target_start, target_end))
    return valid_regions

def merge_overlapping_pathfile(args, gap_re, graph, mxs, mxs_info, scaffolds):
    "Read through pathfile, and merge overlapping pieces, updating path file"
    print(datetime.datetime.today(), ": Finding scaffold overlaps", file=sys.stdout)
    out_pathfile = open(args.p + ".trimmed_scafs.path", 'w')
    with open(args.a, 'r') as path_fin:
        for path in path_fin:
            new_path = []
            path_id, path_seq = path.strip().split("\t")
            path_seq = path_seq.split(" ")
            path_seq = normalize_path(path_seq, gap_re)
            for source, gap, target in zip(path_seq, path_seq[1:], path_seq[2:]):
                gap_match = re.search(gap_re, gap)
                if not gap_match:
                    continue
                if int(gap_match.group(1)) <= args.g + 1 and graph.es()[edge_index(graph, source, target)]["d"] < 0:
                    cuts_found = merge_overlapping(mxs, mxs_info, source, target, scaffolds, args,
                                      graph.es()[edge_index(graph, source, target)]["d"])
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
    for out_scaffold in scaffolds:
        scaffold = scaffolds[out_scaffold]
        if scaffold.ori == "+":
            sequence = scaffold.sequence[scaffold.target_cut:scaffold.source_cut]
        elif scaffold.ori == "-":
            sequence = scaffold.sequence[scaffold.source_cut + args.k:scaffold.target_cut + args.k]
        elif scaffold.ori is None:
            sequence = scaffold.sequence
        else:
            raise ValueError("Invalid orientation for Scaffold:", scaffold)
        if len(sequence) == 0:
            sequence = "N"
        fasta_outfile.write(
            ">{} {}-{}\n{}\n".format(scaffold.ctg_id, scaffold.source_cut, scaffold.target_cut, sequence))
    fasta_outfile.close()

def parse_arguments():
    "Parse arguments for ntLink overlap"
    parser = argparse.ArgumentParser(description="Find coordinates for combining overlapping sequences")
    parser.add_argument("-m", help="Minimizer TSV file", type=str, required=True)
    parser.add_argument("-f", help="Fudge factor for estimated overlap [0.5]", type=float, default=0.5)
    parser.add_argument("-a", help="Path file", required=True, type=str)
    parser.add_argument("-s", help="Scaffold sequences", required=True, type=str)
    parser.add_argument("-k", help="Indexlr k", required=True, type=int)
    parser.add_argument("-d", help="Scaffold dot file", required=True, type=str)
    parser.add_argument("-g", help="Minimum gap size (bp) [20]", default=20, type=int)
    parser.add_argument("--outgap", help="Gap size between trimmed overlapping sequences (bp) [2]", default=2, type=int)
    parser.add_argument("-p", help="Output file prefix [ntlink_merge]", default="ntlink_merge", type=str)
    parser.add_argument("-v", help="Verbose output logging", action="store_true")

    return parser.parse_args()

def print_args(args):
    "Print the parameters for the ntLink overlap step"
    print("Parameters for overlap stage:")
    print("\t-m", args.m)
    print("\t-f", args.f)
    print("\t-a", args.a)
    print("\t-s", args.s)
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
    args.outgap = args.outgap + 1

    scaffolds = read_fasta_file(args.s)
    graph = NtLinkPath.read_scaffold_graph(args.d)

    valid_mx_positions = find_valid_mx_regions(args, gap_re, graph, scaffolds)

    args.m = "/dev/stdin" if args.m == "-" else args.m
    mxs_info, mxs = read_minimizers(args.m, valid_mx_positions)

    merge_overlapping_pathfile(args, gap_re, graph, mxs, mxs_info, scaffolds)

    print_trimmed_scaffolds(args, scaffolds)


if __name__ == '__main__':
    main()
