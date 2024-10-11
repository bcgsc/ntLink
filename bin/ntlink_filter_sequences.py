#!/usr/bin/env python3
"""
ntLink - Filter sequences input to indexlr.
Performance improvement - only compute minimizers on pieces which have a putative overlap
Written by Lauren Coombe @lcoombe
"""

import argparse
import re
import btllib

import ntlink_utils
import ntlink_overlap_sequences

gap_re = re.compile(r'^(\d+)N$')

def filter_sequences(args, scaffolds, valid_mx_regions):
    "Filter the input fasta file to only output those in valid_mx_regions, in order from the path file"
    with open(args.path, 'r') as path_fin:
        for path in path_fin:
            path_name, path_seq = path.strip().split("\t")
            path_seq = path_seq.split(" ")
            path_seq = ntlink_utils.normalize_path(path_seq, gap_re)
            for node in path_seq:
                if re.search(gap_re, node):
                    continue
                node_name = node.strip("+-")
                if node_name in valid_mx_regions:
                    print(">{}\n{}".format(scaffolds[node_name].ctg_id,
                                           scaffolds[node_name].sequence))
            print(">LAST{}\nN".format(path_name))

def read_fasta_file(filename, args):
    "Read a fasta file into memory. Returns dictionary of scafID -> Scaffold"
    scaffolds = {}

    with btllib.SeqReader(filename, btllib.SeqReaderFlag.LONG_MODE, args.t) as reader:
        for record in reader:
            scaffolds[record.id] = ntlink_overlap_sequences.ScaffoldCut(ctg_id=record.id,
                                                                        sequence=record.seq, store_seq=True)

    return scaffolds

def main():
    "Run filtering of input sequences, print to the command line"
    parser = argparse.ArgumentParser(description="Filter input fasta sequences based on ntLink paths")
    parser.add_argument("-s", "--fasta", help="Input fasta file", required=True, type=str)
    parser.add_argument("-d", "--dot", help="Input scaffold graph dot file", required=True, type=str)
    parser.add_argument("-a", "--path", help="Input path file", required=True, type=str)
    parser.add_argument("-k", help="K-mer size (bp)", required=True, type=int)
    parser.add_argument("-f", help="Fudge factor for estimated overlap [0.5]", type=float, default=0.5)
    parser.add_argument("-g", help="Minimum gap size (bp) [20]", required=False, type=int, default=20)
    parser.add_argument("-t", help="Number of threads", default=8, type=int)
    parser.add_argument("-v", "--version", action='version', version='ntLink v1.3.11')

    args = parser.parse_args()


    with ntlink_utils.HiddenPrints():
        graph, _ = ntlink_utils.read_scaffold_graph(args.dot)
        scaffolds = read_fasta_file(args.fasta, args)
        valid_mx_regions = ntlink_utils.find_valid_mx_regions(args, gap_re, graph, scaffolds)

    filter_sequences(args, scaffolds, valid_mx_regions)


if __name__ == "__main__":
    main()
