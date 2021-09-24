'''
General utility functions for ntLink
'''

__author__ = "Lauren Coombe @lcoombe"

import datetime
from collections import namedtuple, defaultdict
import re
import sys
import igraph as ig

from read_fasta import read_fasta

Scaffold = namedtuple("Scaffold", ["id", "length"])

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

    return graph
