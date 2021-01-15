'''
General utility functions for ntLink
'''

__author__ = "Lauren Coombe @lcoombe"

import datetime
from collections import namedtuple
import sys

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
