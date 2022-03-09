'''
General utility functions for ntLink
'''

__author__ = "Lauren Coombe @lcoombe"

import datetime
from collections import namedtuple, defaultdict
import os
import re
import sys
import igraph as ig
import itertools
import numpy as np
import ntlink_pair

from read_fasta import read_fasta

Scaffold = namedtuple("Scaffold", ["id", "length"])

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

def vertex_name(graph, index):
    "Returns vertex name based on vertex id"
    return graph.vs()[index]['name']

def vertex_index(graph, name):
    "Returns vertex index based on vertex name"
    return graph.vs.find(name).index

def edge_index(graph, source_name, target_name):
    "Returns graph edge index based on source/target names"
    return graph.get_eid(source_name, target_name)

def has_vertex(graph, name):
    "Returns True if graph has vertex, else False"
    try:
        graph.vs().find(name)
    except ValueError:
        return False
    return True

def has_estimated_overlap(graph, source, target):
    "Returns True if the edge has an estimated overlap, else False"
    try:
        overlap = graph.es()[edge_index(graph, source, target)]["d"]
        return overlap < 0
    except ig.InternalError:
        return False


def read_fasta_file(filename):
    "Read a fasta file into memory. Returns dictionary of scafID -> Scaffold"
    print(datetime.datetime.today(), ": Reading fasta file", filename, file=sys.stdout)
    scaffolds = {}
    with open(filename, 'r') as fasta:
        for header, seq, _, _ in read_fasta(fasta):
            scaffolds[header] = Scaffold(id=header, length=len(seq))

    return scaffolds

def convert_path_index_to_name(graph, path):
    "Convert path of vertex indices to path of vertex names"
    return [vertex_name(graph, vs) for vs in path]

def reverse_orientation(orientation):
    "Flip the given orientation"
    assert orientation in ("+", "-")
    if orientation == "+":
        return "-"
    return "+"

def reverse_scaf_ori(scaffold):
    "Reverse orientation of scaffold"
    return scaffold[:-1] + reverse_orientation(scaffold[-1])

def read_scaffold_graph(in_graph_file):
    "Reads in a scaffold graph in dot format"
    print(datetime.datetime.today(), ": Reading scaffold file", in_graph_file, file=sys.stdout)

    graph = ig.Graph(directed=True)

    vertices = set()
    edges = defaultdict(dict)  # source -> target -> EdgeInfo

    scaf_num = None

    scaf_num_re = re.compile(r'graph \[scaf_num=(\S+)\]')
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
                edges[source][target] = (int(gap_est), int(num_links))
                continue
            scaf_num_match = re.search(scaf_num_re, line)
            if scaf_num_match:
                try:
                    scaf_num = int(scaf_num_match.group(1))
                except ValueError:
                    scaf_num = None
            elif line != "}":
                print("Error! Unexpected line in input dot file:", line)
                sys.exit(1)

    formatted_edges = [(s, t) for s in edges for t in edges[s]]
    graph.add_vertices(list(vertices))
    graph.add_edges(formatted_edges)

    edge_attributes = {edge_index(graph, s, t): {'d': edges[s][t][0],
                                                 "n": edges[s][t][1]}
                       for s in edges for t in edges[s]}
    graph.es()["d"] = [edge_attributes[e]['d'] for e in sorted(edge_attributes.keys())]
    graph.es()["n"] = [edge_attributes[e]['n'] for e in sorted(edge_attributes.keys())]

    return graph, scaf_num

def find_valid_mx_regions(args, gap_re, graph, scaffolds):
    "Return a dictionary with scaffold -> [(start, end)], marking the valid overlap positions for minimizers on contigs"
    print(datetime.datetime.today(), ": Finding valid minimizer regions", file=sys.stdout)

    valid_regions = {}

    with open(args.path, 'r') as path_fin:
        for path in path_fin:
            _, path_seq = path.strip().split("\t")
            path_seq = path_seq.split(" ")
            path_seq = normalize_path(path_seq, gap_re)
            for source, gap, target in zip(path_seq, path_seq[1:], path_seq[2:]):
                source_noori, target_noori = source.strip("+-"), target.strip("+-")
                gap_match = re.search(gap_re, gap)
                if not gap_match:
                    continue
                if int(gap_match.group(1)) <= args.g + 1 and has_estimated_overlap(graph, source, target):
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

def normalize_path(path_sequence, gap_re):
    "Given a path, normalize it to ensure deterministic running"
    if path_sequence[0].strip("+-") < path_sequence[-1].strip("+-"):
        return path_sequence
    new_seq = []
    for node in reversed(path_sequence):
        if re.search(gap_re, node):
            new_seq.append(node)
        else:
            new_seq.append(reverse_scaf_ori(node))
    return new_seq

def find_valid_mx_region(scaf_noori, scaf_ori, scaffolds, overlap, args, source=True):
    "Return start/end of valid minimizer region on the scaffold"
    if (scaf_ori == "+" and source) or (scaf_ori == "-" and not source):
        start, end = (scaffolds[scaf_noori].length - overlap * -1 - args.k) - \
                     int(overlap * -1 * args.f), scaffolds[scaf_noori].length
    else:
        start, end = 0, int(overlap * -1 * (args.f + 1))

    return start, end


def get_accepted_anchor_contigs(mx_list, read_length, scaffolds, list_mx_info, args):
    "Returns dictionary of contigs of appropriate length, mx hits, whether subsumed"
    MinimizerPositions = namedtuple("MinimizerPositions", ["ctg_pos", "ctg_strand", "read_pos", "read_strand"])
    contig_list = []
    contig_positions = {}  # contig -> [mx positions]
    for mx, pos, strand in mx_list:
        contig = list_mx_info[mx].contig
        if scaffolds[contig].length >= args.z:
            contig_list.append(contig)
            if contig not in contig_positions:
                contig_positions[contig] = []
            contig_positions[contig].append(MinimizerPositions(ctg_pos=list_mx_info[mx].position,
                                                            ctg_strand=list_mx_info[mx].strand,
                                                               read_pos=int(pos),
                                                               read_strand=strand))

    # Filter out hits where mapped length on contig is greater than the read length
    noisy_contigs = set()
    for contig in contig_positions:
        positions = contig_positions[contig]
        if len(positions) < 2:
            continue
        ctg_positions = [position.ctg_pos for position in positions]
        start_idx, end_idx = np.argmin(ctg_positions), np.argmax(ctg_positions)
        start_positions = positions[start_idx]
        end_positions = positions[end_idx]
        if args.x == 0:
            if abs(end_positions.ctg_pos - start_positions.ctg_pos) > read_length + args.k:
                noisy_contigs.add(contig)
        else:
            threshold = min(read_length + args.k,
                            (args.x * abs(end_positions.read_pos - start_positions.read_pos)) + args.k)
            if abs(end_positions.ctg_pos - start_positions.ctg_pos) > threshold:
                noisy_contigs.add(contig)
    contig_list = [contig for contig in contig_list if contig not in noisy_contigs]

    contig_runs = [(ctg, len(list(hits)), list(hits)) for ctg, hits in itertools.groupby(contig_list)]
    contigs_hits = {}
    for i, run_tup in enumerate(contig_runs):
        ctg, cnt, list_hits = run_tup
        if ctg in contigs_hits:
            for j in range(contigs_hits[ctg].index + 1, i):
                contigs_hits[contig_runs[j][0]].subsumed = True
            contigs_hits[ctg].hits.extend(list_hits)
        else:
            contigs_hits[ctg] = ntlink_pair.ContigRun(ctg, i, cnt)
            contigs_hits[ctg].hits.extend(list_hits)

    return_contigs_hits = {ctg: contigs_hits[ctg] for ctg in contigs_hits if not contigs_hits[ctg].subsumed}

    return_contig_runs_tmp = [ctg for ctg, hits, list_hits in contig_runs if not contigs_hits[ctg].subsumed]
    return_contig_runs = [ctg for ctg, hits in itertools.groupby(return_contig_runs_tmp)]

    for ctg in contigs_hits:
        contigs_hits[ctg].hits = contig_positions[ctg]

    return return_contigs_hits, return_contig_runs