#!/usr/bin/env python3
"""
ntLink - Filter sequences input to indexlr.
Performance improvement - only compute minimizers on pieces which have a putative overlap
Written by Lauren Coombe @lcoombe
"""

import argparse
import re
from read_fasta import read_fasta
import ntlink_utils
import ntlink_overlap_sequences

def filter_sequences(fa_in, valid_mx_regions):
    "Filter the input fasta file to only output those in valid_mx_regions"
    with open(fa_in, 'r') as fin:
        for name, seq, _, _ in read_fasta(fin):
            if name in valid_mx_regions:
                print(">{}\n{}".format(name, seq))

def main():
    "Run filtering of input sequences, print to the command line"
    parser = argparse.ArgumentParser(description="Filter input fasta sequences based on ntLink paths")
    parser.add_argument("-s", "--fasta", help="Input fasta file", required=True, type=str)
    parser.add_argument("-d", "--dot", help="Input scaffold graph dot file", required=True, type=str)
    parser.add_argument("-a", "--path", help="Input path file", required=True, type=str)
    parser.add_argument("-k", help="K-mer size (bp)", required=True, type=int)
    parser.add_argument("-f", help="Fudge factor for estimated overlap [0.5]", type=float, default=0.5)
    parser.add_argument("-g", help="Minimum gap size (bp) [20]", required=False, type=int, default=20)
    parser.add_argument("-v", "--version", action='version', version='ntLink v1.2.0')

    args = parser.parse_args()

    gap_re = re.compile(r'^(\d+)N$')

    with ntlink_utils.HiddenPrints():
        graph, _ = ntlink_utils.read_scaffold_graph(args.dot)
        scaffolds = ntlink_overlap_sequences.read_fasta_file_trim_prep(args.fasta)
        valid_mx_regions = ntlink_utils.find_valid_mx_regions(args, gap_re, graph, scaffolds)

    filter_sequences(args.fasta, valid_mx_regions)


if __name__ == "__main__":
    main()
