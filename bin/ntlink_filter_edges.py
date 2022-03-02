#!/usr/bin/env python3
"""
Use minimizers to overlap the sequence pairs that are likely overlapping
"""
import sys
import btllib
import argparse
import re
from collections import namedtuple, defaultdict
import ntlink_utils
import ntjoin_utils
import ntlink_overlap_sequences
import numpy as np
import igraph as ig
import ntlink_pair

Scaffold = namedtuple("Scaffold", ["seq", "length"])

MappedPathInfo = namedtuple("MappedPathInfo",
                            ["mapped_region_length", "mid_mx", "median_length_from_end", "num_nodes",
                             "source_range", "target_range"])

OverlapEdgeInfo = namedtuple("OverlapEdgeInfo", ["source", "target", "gap_estimate", "detected_overlap", "num_paths",
                                                 "verbose_info"])


def read_fasta(fasta_filename: str) -> dict:
    "Read sequences into memory"
    contigs = {}
    with btllib.SeqReader(fasta_filename, btllib.SeqReaderFlag.LONG_MODE) as reads:
        for read in reads:
            contigs[read.id] = Scaffold(read.seq, len(read.seq))
    return contigs

def assess_edge(source: str, target: str, sequences: dict, gap_estimate: int, args: argparse.Namespace, edge_info: dict) -> None:
    "Assess the edge using the overlap code"
    source_start, source_end = ntlink_utils.find_valid_mx_region(source.strip("+-"), source[-1],
                                                                 sequences, gap_estimate, args, source=True)
    target_start, target_end = ntlink_utils.find_valid_mx_region(target.strip("+-"), target[-1],
                                                                 sequences, gap_estimate, args, source=False)
    out_fasta = open("tmp.fa", 'w') # !! TODO: Right now, only printing the overlapping region segment
    out_fasta.write(f">{source.strip('+-')}\n{sequences[source.strip('+-')].seq[source_start:source_end]}\n")
    out_fasta.write(f">{target.strip('+-')}\n{sequences[target.strip('+-')].seq[target_start:target_end]}\n")
    out_fasta.close()

    mx_info = defaultdict(dict) # contig -> mx -> (contig, position)
    mxs = {} # Contig -> [List minimizers]

    with btllib.Indexlr("tmp.fa", args.k, 10, btllib.IndexlrFlag.LONG_MODE) as minimizers:
        for minimizer_line in minimizers:
            dup_mxs = set()
            name = minimizer_line.id
            for mx in minimizer_line.minimizers:
                mx_hash, pos = str(mx.out_hash), mx.pos
                if name in mx_info and mx_hash in mx_info[name]:
                    dup_mxs.add(mx_hash)
                else:
                    mx_info[name][mx_hash] = (name, pos)
            mx_info[name] = {mx: mx_info[name][mx] for mx in mx_info[name] if mx not in dup_mxs}
            mxs[name] = [[str(mx_pos.out_hash) for mx_pos in minimizer_line.minimizers
                          if str(mx_pos.out_hash) not in dup_mxs and str(mx_pos.out_hash) in mx_info[name]]]

    return_info = merge_overlapping(mxs, mx_info, source, target, sequences, args, 0)
    if not return_info:
        edge_info[(source, target)] = OverlapEdgeInfo(source, target, gap_estimate, 0, 0, return_info)
    else:
        edge_info[(source, target)] = OverlapEdgeInfo(source, target, gap_estimate,
                                                      return_info[0].mapped_region_length*-1, len(return_info), return_info)

def merge_overlapping(list_mxs, list_mx_info, source, target, scaffolds, args, gap):
    "Find the cut points for overlapping adjacent contigs"
    source_noori = source.strip("+-")
    target_noori = target.strip("+-")

    weights = {source_noori: 1, target_noori: 1}

    list_mxs_pair = {source_noori: list_mxs[source_noori], target_noori: list_mxs[target_noori]}
   # list_mxs_pair = filter_minimizers_position(list_mxs_pair, source, target, gap, scaffolds, list_mx_info, args)
    with ntlink_utils.HiddenPrints():
        list_mxs_pair = ntjoin_utils.filter_minimizers(list_mxs_pair)

    graph = ntlink_overlap_sequences.build_graph(list_mxs_pair, weights)
    graph = ntlink_overlap_sequences.filter_graph_global(graph, 2)

    # Print the DOT graph if in verbose mode
    if args.v:
        ntlink_overlap_sequences.print_graph(graph, list_mx_info, args.p)

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
            mid_mx_dist_end_source = ntlink_overlap_sequences.get_dist_from_end(source[-1],
                                                       list_mx_info[source_noori][mid_mx][1],
                                                       scaffolds[source_noori].length)
            mid_mx_dist_end_target = ntlink_overlap_sequences.get_dist_from_end(target[-1],
                                                       list_mx_info[target_noori][mid_mx][1],
                                                       scaffolds[target_noori].length, target=True)
            paths_components.append(MappedPathInfo(mapped_region_length=np.median([source_align_len,
                                                                                   target_align_len]),
                                                   mid_mx=mid_mx,
                                                   median_length_from_end=np.median(
                                                       [mid_mx_dist_end_source, mid_mx_dist_end_target]), num_nodes=len(path),
                                                   source_range=(source_start, source_end),
                                                   target_range=(target_start, target_end)))
        elif singleton_nodes:
            assert len(singleton_nodes) == 1
            mid_mx = ntlink_utils.vertex_name(component_graph, singleton_nodes[0])
            mid_mx_dist_end_source = ntlink_overlap_sequences.get_dist_from_end(source[-1], list_mx_info[source_noori][mid_mx][1],
                                                       scaffolds[source_noori].length)
            mid_mx_dist_end_target = ntlink_overlap_sequences.get_dist_from_end(target[-1], list_mx_info[target_noori][mid_mx][1],
                                                       scaffolds[target_noori].length, target=True)
            paths_components.append(MappedPathInfo(mapped_region_length=1, mid_mx=mid_mx,
                                                   median_length_from_end=np.median([mid_mx_dist_end_source,
                                                                                     mid_mx_dist_end_target]), num_nodes=1,
                                                   source_range=(list_mx_info[source_noori][mid_mx][1], list_mx_info[source_noori][mid_mx][1]),
                                                   target_range=(list_mx_info[target_noori][mid_mx][1], list_mx_info[target_noori][mid_mx][1])))
        else:
            print("NOTE: non-singleton, {} source nodes".format(len(source_nodes)))
    if not paths_components:
        return False
    paths = sorted(paths_components, key=lambda x: (x.num_nodes, x.mapped_region_length, x.median_length_from_end,
                                                   x.mid_mx), reverse=True)
    return paths


def assess_scaffold_graph_edges(args: argparse.Namespace, fasta_seqs: dict) -> dict:
    "Assess each overlapping edge based on minimizer overlap, find mapping region length (and #)"
    edge_re = re.compile(r'^\"(\S+)\" -> \"(\S+)\"\s+\[d=([+-]?\d+)\s+e=\d+\s+n=(\d+)')
    #print("source", "target", "gap_estimate", "detected_overlap", "num_paths", "verbose_info", sep="\t")
    visited = set()
    edge_info = {} # (source, target) -> EdgeInfo
    with open(args.g, 'r') as fin:
        for line in fin:
            line = line.strip()
            edge_match = re.search(edge_re, line)
            if edge_match:
                source, target, gap_est, n = edge_match.group(1), edge_match.group(2), \
                                             int(edge_match.group(3)), int(edge_match.group(4))
                if gap_est >= 0 or n < args.n:
                    continue
                if (reverse_orientation(target), reverse_orientation(source)) in visited:
                    continue
                assess_edge(source, target, fasta_seqs, gap_est, args, edge_info)
                visited.add((source, target))
    return edge_info

def reverse_orientation(node: str) -> str:
    "Reverse the orientation of the given path node"
    if node[-1] == "+":
        return node[:-1] + "-"
    if node[-1] == "-":
        return node[:-1] + "+"
    raise ValueError("Node must end in + or -")

def read_checkpoint_edges(checkpoint_filename: str) -> dict:
    "Read checkpoint file into dictionary"
    edge_info = {}
    with open(checkpoint_filename, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            #"source", "target", "gap_estimate", "detected_overlap", "num_paths", "verbose_info"
            source, target, gap_est, detect_overlap, num_paths, verbose_info = line
            if source == "source":
                continue
            edge_info[(source, target)] = OverlapEdgeInfo(source, target, float(gap_est), float(detect_overlap), int(num_paths), verbose_info)
    return edge_info

def get_overlap_diff(edge_info: dict, source: int, target:int) -> int:
    "Get overlap (detected, estimated) difference for the edge"
    if (source, target) in edge_info:
        return abs(edge_info[(source, target)].detected_overlap - edge_info[(source, target)].gap_estimate)
    elif (reverse_orientation(target), reverse_orientation(source)) in edge_info:
        return abs(edge_info[(reverse_orientation(target), reverse_orientation(source))].detected_overlap -
                   edge_info[(reverse_orientation(target), reverse_orientation(source))].gap_estimate)
    raise ValueError("{}, {} edge not found in edge_info dict".format(source, target))

def choose_best_edge(in_edge_list: list, edge_info: dict, graph: ig.Graph) -> ig.Edge:
    "Choose 'best' edge - difference of estimated to detected"
    diffs = sorted([(edge.index, get_overlap_diff(edge_info, ntlink_utils.vertex_name(graph, edge.source), ntlink_utils.vertex_name(graph, edge.target))) for edge in in_edge_list],
                   key=lambda x: x[1])
    return diffs[0][0]


def filter_scaffold_graph_edges(edge_info: dict, args: argparse.Namespace, sequences: dict) -> None:
    "Filter the scaffold graph edges"
    graph, _ = ntlink_utils.read_scaffold_graph(args.g)
    edges_to_remove = set()
    for node in graph.vs():
        in_edge_overlap = [edge for edge in node.in_edges() if edge['d'] < 0]
        if len(in_edge_overlap) > 1:
            best_edge = choose_best_edge(in_edge_overlap, edge_info, graph)
            for edge in in_edge_overlap:
                if edge.index != best_edge:
                    edges_to_remove.add(edge.index)

        out_edge_overlap = [edge for edge in node.out_edges() if edge['d'] < 0]
        if len(out_edge_overlap) > 1:
            best_edge = choose_best_edge(out_edge_overlap, edge_info, graph)
            for edge in out_edge_overlap:
                if edge.index != best_edge:
                    edges_to_remove.add(edge.index)
    graph_copy = graph.copy()
    graph_copy.delete_edges(edges_to_remove)

    ntlink_pair.NtLink.print_directed_graph(graph_copy, args.o, sequences)


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify edges based on minimzer overlap checks")
    parser.add_argument("-g", help="Scaffold graph from ntLink", required=True, type=str)
    parser.add_argument("--fasta", help="Fasta file with sequences", required=True, type=str)
    parser.add_argument("-n", help="Minimum edge weight to consider edge", required=False, default=2, type=int)
    parser.add_argument("-f", help="Fudge factor [0.5]", default=0.5, type=float, required=False)
    parser.add_argument("-k", help="K-mer size for indexlr", default=15, type=int, required=False)
    parser.add_argument("--checkpoint", help="Checkpoint file", required=False)
    parser.add_argument("-v", help="Verbose mode", action="store_true")
    parser.add_argument("-o", help="Output filename", default="ntlink_scaffold_graph.dot", required=False, type=str)
    args = parser.parse_args()

    sequences = read_fasta(args.fasta)

    if not args.checkpoint:
        edge_info = assess_scaffold_graph_edges(args, sequences)
    else:
        edge_info = read_checkpoint_edges(args.checkpoint)

    filter_scaffold_graph_edges(edge_info, args, sequences)

if __name__ == "__main__":
    main()
