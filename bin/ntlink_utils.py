'''
General utility functions for ntLink
'''

__author__ = "Lauren Coombe @lcoombe"

import igraph as ig
import datetime
from read_fasta import read_fasta
from collections import namedtuple
import sys

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