'''
General utility functions for ntLink
'''

__author__ = "Lauren Coombe @lcoombe"

import datetime
from collections import namedtuple
import sys
import igraph as ig

from read_fasta import read_fasta

Scaffold = namedtuple("Scaffold", ["id", "length", "sequence"])

def vertex_name(graph, index):
    "Returns vertex name based on vertex id"
    return graph.vs[index]['name']

def edge_index(graph, source_name, target_name):
    "Returns graph edge index based on source/target names"
    return graph.get_eid(source_name, target_name)

def read_fasta_file(filename):
    "Read a fasta file into memory. Returns dictionary of scafID -> Scaffold"
    print(datetime.datetime.today(), ": Reading fasta file", filename, file=sys.stdout)
    scaffolds = {}
    with open(filename, 'r') as fasta:
        for header, seq, _, _ in read_fasta(fasta):
            scaffolds[header] = Scaffold(id=header, length=len(seq), sequence=seq)

    return scaffolds

def filter_graph(graph, min_weight):
    "Filter the graph by edge weights on edges incident to branch nodes"
    in_branch_nodes = [node.index for node in graph.vs() if node.indegree() > 1]
    to_remove_in_edges = [edge for node in in_branch_nodes for edge in graph.incident(node, mode=ig.IN) # pylint: disable=no-member
                          if graph.es()[edge]['n'] < min_weight]
    new_graph = graph.copy()
    new_graph.delete_edges(to_remove_in_edges)

    out_branch_nodes = [node.index for node in new_graph.vs() if node.outdegree() > 1]
    to_remove_out_edges = [edge for node in out_branch_nodes for edge in new_graph.incident(node, mode=ig.OUT) # pylint: disable=no-member
                           if new_graph.es()[edge]['n'] < min_weight]

    return_graph = new_graph.copy()
    return_graph.delete_edges(to_remove_out_edges)

    return return_graph

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
