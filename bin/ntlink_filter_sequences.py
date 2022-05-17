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

def filter_sequences(fa_in, valid_mx_regions, mask_internal=False):
    "Filter the input fasta file to only output those in valid_mx_regions"
    with open(fa_in, 'r') as fin:
        for name, seq, _, _ in read_fasta(fin):
            if name in valid_mx_regions:
                if mask_internal:
                    valid_regions = sorted(valid_mx_regions[name])
                    five_prime_coord, three_prime_coord = 0, len(seq)
                    for region in valid_regions:
                        if region[0] == 0: # 5' end
                            five_prime_coord = region[1]
                        if region[1] == len(seq): # 3' end
                            three_prime_coord = region[0]
                    assert len(valid_regions) > 0 and len(valid_regions) < 3
                    new_seq = seq[:five_prime_coord] + 'N'*(three_prime_coord - five_prime_coord) + seq[three_prime_coord:]
                    assert len(new_seq) == len(seq)
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
    parser.add_argument("--mask", help="Mask internal region of sequences", required=False, action="store_true")
    parser.add_argument("-v", "--version", action='version', version='ntLink v1.2.1')

    args = parser.parse_args()

    gap_re = re.compile(r'^(\d+)N$')

    with ntlink_utils.HiddenPrints():
        graph, _ = ntlink_utils.read_scaffold_graph(args.dot)
        scaffolds = ntlink_overlap_sequences.read_fasta_file_trim_prep(args.fasta)
        valid_mx_regions = ntlink_utils.find_valid_mx_regions(args, gap_re, graph, scaffolds)

    filter_sequences(args.fasta, valid_mx_regions, mask_internal=args.mask)


if __name__ == "__main__":
    main()
