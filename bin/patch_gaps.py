#!/usr/bin/env python3
'''
Use minimizer anchors to patch gaps in scaffolds
'''
import argparse
import re
from collections import namedtuple
import itertools
import btllib
import numpy
import sys
import ntlink_pair
import ntlink_utils

MinimizerEntry = namedtuple("MinimizerEntry", ["ctg_pos", "read_pos"])

class ScaffoldGaps:
    def __init__(self, seq):
        self.seq = seq
        self.length = len(seq)
        self.five_prime_cut = 0
        self.three_prime_cut = self.length

    def __str__(self):
        return f"Length:{self.length} 5'cut:{self.five_prime_cut} 3'cut:{self.three_prime_cut}"

    def get_cut_sequence(self):
        return self.seq[self.five_prime_cut: self.three_prime_cut+1]

class PairInfo:
    def __init__(self, gap_size):
        self.gap_size = int(gap_size)
        self.mapping_reads = set()
        self.chosen_read = None
        self.source_ctg_cut = None
        self.source_read_cut = None
        self.target_ctg_cut = None
        self.target_read_cut = None


    def __str__(self):
        return f"Gap: {self.gap_size}; Chosen read: {self.chosen_read}; source ctg/read cuts: {self.source_ctg_cut}/{self.source_read_cut}" \
               f"target ctg/read cuts: {self.target_ctg_cut}/{self.target_read_cut}"

    def get_cut_read_sequence(self, reads): # !!TODO consider orientation of read
        return reads[self.chosen_read][self.source_read_cut: self.target_read_cut]

def read_path_file_pairs(path_filename: str) -> dict:
    "Read through the path file, storing the found pairs: pair -> gap estimate"
    pairs = {}
    gap_re = re.compile('^(\d+)N$')
    with open(path_filename, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if len(line) < 2:
                continue
            ctg_id, path = line
            path = path.split(" ")
            for idx in range(len(path) - 2):
                i, j, k = path[idx:idx+3]
                gap_match = re.search(gap_re, j)
                if gap_match:
                    pairs[(i, k)] = PairInfo(gap_match.group(1))
    return pairs

def parse_minimizers(minimizer_positions: str) -> list:
    "Parse the minimizer positions string"
    mx_pos_re = re.compile(f'MinimizerPositions\(ctg_pos=(\d+), read_pos=(\d+)\)')
    return_mxs = []
    for match in re.findall(mx_pos_re, minimizer_positions):
        return_mxs.append(MinimizerEntry(ctg_pos=int(match[0]), read_pos=int(match[1])))
    return return_mxs


def find_orientation(mx_positions: list) -> str:
    "Infer orientation from minimizer positions"
    if all(x < y for x, y in zip(mx_positions, mx_positions[1:])):
        return "+"
    if all(x > y for x, y in zip(mx_positions, mx_positions[1:])):
        return "-"
    return "None"


def reverse_complement_pair(source: str, target: str) -> tuple:
    "Reverse complement the given contig pair"
    if source[-1] == "+":
        source_new_ori = "-"
    elif source[-1] == "-":
        source_new_ori = "+"
    else:
        raise(ValueError("+ or - needed for last character of node, found" + source[-1]))

    if target[-1] == "+":
        target_new_ori = "-"
    elif target[-1] == "-":
        target_new_ori = "+"
    else:
        raise(ValueError("+ or - needed for last character of node, found" + target[-1]))

    return target[:-1] + target_new_ori, source[:-1] + source_new_ori


def tally_contig_mapping_info(read_id: str, mappings: list, read_info: dict, pairs: dict) -> None:
    "Tally the contig mapping information for the given read"
    read_info[read_id] = {}
    mapping_order = []
    for _, ctg_id, anchors, minimizer_positions in mappings:
        minimizer_positions = parse_minimizers(minimizer_positions)
        orientation = find_orientation(minimizer_positions)
        if orientation not in ["+", "-"]:
            continue
        read_info[read_id][ctg_id] = (int(anchors), minimizer_positions, orientation)
        mapping_order.append(ctg_id + orientation)

    for i, j in itertools.combinations(mapping_order, 2):
        if (i, j) in pairs:
            pairs[(i, j)].mapping_reads.add(read_id)
        if reverse_complement_pair(i, j) in pairs:
            pairs[reverse_complement_pair(i, j)].mapping_reads.add(read_id)



def read_verbose_mappings(mappings_filename: str, pairs: dict) -> dict:
    "Read through the verbose mappings, tallying the information for each read"
    read_info = {}
    curr_read_id = None
    curr_mappings = []

    with open(mappings_filename, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if line[0] != curr_read_id and curr_read_id is not None:
                # Tally information for the contig's mappings
                tally_contig_mapping_info(curr_read_id, curr_mappings, read_info, pairs)
                curr_read_id = line[0]
                curr_mappings = [line]
            else:
                curr_read_id = line[0]
                curr_mappings.append(line)
    if curr_read_id is not None:
        tally_contig_mapping_info(curr_read_id, curr_mappings, read_info, pairs)


    return read_info

def read_scaffold_file(scaffold_filename: str) -> dict:
    "Read the scaffolds into memory"
    scaffolds = {}
    with btllib.SeqReader(scaffold_filename, btllib.SeqReaderFlag.LONG_MODE) as reads:
        for read in reads:
            scaffolds[read.id] = ScaffoldGaps(read.seq)
    return scaffolds

def choose_best_read_per_pair(pairs: dict, mappings: dict) -> None:
    "For each pair, choose the 'best' read to fill in the gap"
    for source, target in pairs:
        reads = [(read_id, mappings[read_id][source.strip("+-")][0],
                  mappings[read_id][target.strip("+-")][0])
                 for read_id in pairs[(source, target)].mapping_reads]
        sorted_reads = sorted(reads, key=lambda x: numpy.mean([x[1], x[2]]), reverse=True)
        pairs[(source, target)].chosen_read = sorted_reads[0][0]


def get_gap_fill_reads(reads_filename: str, pairs: dict) -> dict:
    "Collect the reads needed for gap-filling"
    reads = {} # read_id -> sequence
    target_reads = {pairs[pair].chosen_read for pair in pairs}
    with btllib.SeqReader(reads_filename, btllib.SeqReaderFlag.LONG_MODE) as reads_in:
        for read in reads_in:
            if read.id in target_reads:
                reads[read.id] = read.seq
    return reads

def find_cut_points(sequences: dict, pairs: dict, reads: dict, mappings: dict) -> None:
    "Find the points for cuts for gap-filling"
    for source, target in pairs:
        read_id = pairs[(source, target)].chosen_read
        print(read_id)
        source_read_mxs = mappings[read_id][source.strip("+-")][1]
        source_ctg_pos, source_read_pos = source_read_mxs[-1] #!!TODO: change with orientation considered?

        target_read_mxs = mappings[read_id][target.strip("+-")][1]
        target_ctg_pos, target_read_pos = target_read_mxs[0]

        pairs[(source, target)].source_ctg_cut = source_ctg_pos
        pairs[(source, target)].source_read_cut = source_read_pos
        pairs[(source, target)].target_ctg_cut = target_ctg_pos
        pairs[(source, target)].target_read_cut = target_read_pos

        print((source, target), pairs[(source, target)])

        #print(">source", sequences[source.strip("+-")].seq[:source_ctg_pos], file=sys.stderr, sep="\n")
        #print(">target", sequences[target.strip("+-")].seq[target_ctg_pos:], file=sys.stderr, sep="\n")
        #print(">read_fill", reads[read_id][source_read_pos:target_read_pos+1], file=sys.stderr, sep="\n")

def print_masked_sequences(scaffolds: dict, reads: dict, pairs: dict, args: argparse.Namespace):
    "Print the masked sequences to a temp file"
    out_scaffolds = open(args.s + ".masked_temp.fa", 'w')
    out_reads = open(args.reads + ".masked_temp.fa", 'w')

    for source, target in pairs:
        pair_info = pairs[(source, target)]
        ## !! TODO allow fudge factor with sequences printed out?
        source_cut = pair_info.source_ctg_cut
        source_noori = source.strip("+-")
        if source[-1] == "+":
            out_scaffolds.write(f">{source}_source\n"
                                f"{'N'*source_cut}"
                                f"{scaffolds[source_noori].seq[source_cut:]}\n")
        elif source[-1] == "-":
            out_scaffolds.write(f">{source}_source\n"
                                f"{scaffolds[source_noori].seq[:source_cut]}"
                                f"{'N'*(len(scaffolds[source_noori].seq) - source_cut)}\n")
        else:
            raise ValueError(f"{source[-1]} must be + or -")

        target_cut = pair_info.target_ctg_cut
        target_noori = target.strip("+-")
        if target[-1] == "+":
            out_scaffolds.write(f">{target}_target\n"
                                f"{scaffolds[target_noori].seq[:target_cut]}"
                                f"{'N'*(len(scaffolds[target_noori].seq) - target_cut)}\n")
        elif target[-1] == "-":
            out_scaffolds.write(f">{target}_target\n"
                                f"{'N'*target_cut}"
                                f"{scaffolds[target_noori].seq[target_cut:]}")
        else:
            raise ValueError(f"{target[-1]} must be + or -")

        read_start, read_end = min(pair_info.source_read_cut, pair_info.target_read_cut), \
                               max(pair_info.source_read_cut, pair_info.target_read_cut)
        out_reads.write(f">{pair_info.chosen_read}__{source}__{target}\n" #!! TODO: consider orientation of read?
                        f"{'N'*read_start}{reads[pair_info.chosen_read][read_start:read_end]}"
                        f"{'N'*(len(reads[pair_info.chosen_read]) - read_end)}\n")

    out_scaffolds.close()
    out_reads.close()

def convert_btllib_strand(strand: bool) -> str:
    "Convert btllib strand to string"
    if strand:
        return "+"
    return "-"

def read_btllib_minimizers(minimizer_entries: list) -> dict:
    "Read minimizers from list of Indexlr entries from btllib"
    mx_info = {}
    dup_mxs = set()
    for entry in minimizer_entries:
        for minimizer in entry.minimizers:
            mx_hash = str(minimizer.out_hash)
            if mx_hash in mx_info:
                dup_mxs.add(mx_hash)
            else:
                mx_info[mx_hash] = ntlink_pair.Minimizer(entry.id, minimizer.pos, convert_btllib_strand(minimizer.forward))
    mx_info = {mx: mx_info[mx] for mx in mx_info if mx not in dup_mxs}
    return mx_info

def map_long_reads(pairs: dict, scaffolds: dict, args: argparse.Namespace) -> None:
    "Map the long reads to the sequences, print verbose output (for now)"
    read_header_re = re.compile(r'^(\S+)__(\S+)__(\S+)$')
    scaffold_header_re = re.compile(r'^(\S+)_(source|target)$')

    with btllib.Indexlr(args.s + ".masked_temp.fa", 15, 10, btllib.IndexlrFlag.LONG_MODE) as scaffolds_btllib: # !!TODO magic numbers
        with btllib.Indexlr(args.reads + ".masked_temp.fa", 15, 10, btllib.IndexlrFlag.LONG_MODE) as reads:
            for chosen_read in reads:
                print(chosen_read.id)
                read_id, source, target = re.search(read_header_re, chosen_read.id).groups()

                # Read source scaffold
                source_scaf = scaffolds_btllib.read()
                source_id, label = re.search(scaffold_header_re, source_scaf.id).groups()
                source_scaf.id = source_id.strip("+-")
                assert source_id == source
                assert label == "source"

                # Read target scaffold
                target_scaf = scaffolds_btllib.read()
                target_id, label = re.search(scaffold_header_re, target_scaf.id).groups()
                target_scaf.id = target_id.strip("+-")
                assert target_id == target
                assert label == "target"

                mx_info = read_btllib_minimizers([source_scaf, target_scaf])

                mxs = [(str(mx.out_hash), mx.pos, convert_btllib_strand(mx.forward)) for mx in chosen_read.minimizers if str(mx.out_hash) in mx_info]
                accepted_anchor_contigs, contig_runs = ntlink_utils.get_accepted_anchor_contigs(mxs, chosen_read.readlen,
                                                                                    scaffolds, mx_info, args)
                assert len(accepted_anchor_contigs) == 2
                for ctg_run in accepted_anchor_contigs:
                    ctg_run_entry = accepted_anchor_contigs[ctg_run]
                    if ctg_run_entry.contig == source_scaf.id:
                        terminal_mx = ctg_run_entry.hits[-1]
                        pairs[(source, target)].source_ctg_cut = terminal_mx.ctg_pos
                        pairs[(source, target)].source_read_cut = terminal_mx.read_pos
                        if source[-1] == "+":
                            scaffolds[source_scaf.id].three_prime_cut = terminal_mx.ctg_pos
                        else:
                            scaffolds[source_scaf.id].five_prime_cut = terminal_mx.ctg_pos

                    if ctg_run_entry.contig == target_scaf.id:
                        terminal_mx = ctg_run_entry.hits[0]
                        pairs[(source, target)].target_ctg_cut = terminal_mx.ctg_pos
                        pairs[(source, target)].target_read_cut = terminal_mx.read_pos
                        if target[-1] == "+":
                            scaffolds[target_scaf.id].five_prime_cut = terminal_mx.ctg_pos
                        else:
                            scaffolds[target_scaf.id].three_prime_cut = terminal_mx.ctg_pos

    print([str(pairs[pair]) for pair in pairs])
    print(str(scaffolds[source.strip("+-")]))
    print(str(scaffolds[target.strip("+-")]))

def main() -> None:
    parser = argparse.ArgumentParser(description="Use minimizer mappings to fill gaps")
    parser.add_argument("--path", help="Input path file for gap patching", required=True, type=str)
    parser.add_argument("--mappings", help="ntLink verbose mapping TSV", required=True, type=str)
    parser.add_argument("-s", help="Input scaffolds", required=True, type=str)
    parser.add_argument("--reads", help="Input reads", required=True, type=str)
    parser.add_argument("-z", help="Minimum contig size (bp) [1000]", type=int, required=False, default=1000)
    parser.add_argument("-k", help="Kmer size used in minimizer step [15]", type=int, required=False, default=15)
    parser.add_argument("-x", help="Fudge factor", type=float, required=False, default=0)
    parser.add_argument("-m", help="PLACEHOLDER", type=str)
    args = parser.parse_args()

    # Read path file into pairs ((source, target) -> (gap_est, supporting reads)
    pairs = read_path_file_pairs(args.path)

    # Read through verbose mappings read_id -> sequence_ori -> (anchor, minimizer_pos_list)
    mappings = read_verbose_mappings(args.mappings, pairs)

    # Read scaffold sequences into memory sequence_id -> [sequence, start, end]
    print("Reading scaffolds..")
    sequences = read_scaffold_file(args.s)

    # Choose best read for patching each pair's gap (adjust pairs supporting reads)
    print("Choosing best read..")
    choose_best_read_per_pair(pairs, mappings)

    # Read in the reads that are needed
    reads = get_gap_fill_reads(args.reads, pairs)

    # Find cut points
    find_cut_points(sequences, pairs, reads, mappings)

    # Print masked sequences for assembly, reads for minimizer generation
    print_masked_sequences(sequences, reads, pairs, args)

    # Compute minimizers, and map the long read to the sequences at a lower k/w
    map_long_reads(pairs, sequences, args)

    for source, target in pairs:
        print(">source", sequences[source.strip("+-")].get_cut_sequence(), file=sys.stderr, sep="\n")
        print(">target", sequences[target.strip("+-")].get_cut_sequence(), file=sys.stderr, sep="\n")
        print(">read_fill", pairs[(source, target)].get_cut_read_sequence(reads), file=sys.stderr, sep="\n")


if __name__ == "__main__":
    main()
