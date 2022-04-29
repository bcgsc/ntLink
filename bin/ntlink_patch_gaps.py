#!/usr/bin/env python3
'''
Use minimizer anchors to patch gaps in scaffolds
'''
import argparse
import re
from collections import namedtuple, Counter
import itertools
import datetime
import shlex
import subprocess
import sys
import btllib
import numpy
import ntlink_pair
import ntlink_utils
from ntlink_utils import MinimizerPositions

MinimizerMapping = namedtuple("MinimizerMapping", ["anchors", "minimizer_positions", "orientation"])

class ScaffoldGaps:
    "Representing a scaffold, adjusting for gap-filling coordinates"
    def __init__(self, seq):
        self.seq = seq
        self.length = len(seq)
        self.five_prime_cut = 0
        self.three_prime_cut = self.length
        self.five_prime_trim = 0
        self.three_prime_trim = self.length

    def __str__(self):
        return f"Length:{self.length} 5'cut:{self.five_prime_cut} 3'cut:{self.three_prime_cut}"

    def get_cut_sequence(self, reverse_compl: str):
        "Get the cut sequence, reverse complementing if indicated"
        five_prime_cut_site = max(self.five_prime_trim, self.five_prime_cut)
        three_prime_cut_site = min(self.three_prime_trim, self.three_prime_cut)
        if reverse_compl == "-":
            return self.reverse_complement(self.seq[five_prime_cut_site: three_prime_cut_site])
        return self.seq[five_prime_cut_site: three_prime_cut_site]

    @staticmethod
    def reverse_complement(sequence):
        "Reverse complements a given sequence"
        translation_table = str.maketrans(
            "ACGTUNMRWSYKVHDBacgtunmrwsykvhdb",
            "TGCAANKYWSRMBDHVtgcaankywsrmbdhv")
        return sequence[::-1].translate(translation_table)

class PairInfo:
    "Information about a pair in the path file"
    def __init__(self, gap_size):
        self.gap_size = int(gap_size)
        self.mapping_reads = set()
        self.chosen_read = None
        self.source_ctg_cut = None
        self.source_read_cut = None
        self.target_ctg_cut = None
        self.target_read_cut = None
        self.old_anchor_used = False

    @staticmethod
    def reverse_complement(sequence):
        "Reverse complements a given sequence"
        translation_table = str.maketrans(
            "ACGTUNMRWSYKVHDBacgtunmrwsykvhdb",
            "TGCAANKYWSRMBDHVtgcaankywsrmbdhv")
        return sequence[::-1].translate(translation_table)


    def __str__(self):
        return f"Gap: {self.gap_size}; Chosen read: {self.chosen_read}; " \
               f"source ctg/read cuts: {self.source_ctg_cut}/{self.source_read_cut}" \
               f" target ctg/read cuts: {self.target_ctg_cut}/{self.target_read_cut}; " \
               f"Anchor used: {self.old_anchor_used}"

    def get_cut_read_sequence(self, reads, reverse_compl: str):
        "Get the cut sequence, reverse complementing if indicated"
        if reverse_compl == "-":
            return self.reverse_complement(reads[self.chosen_read][self.target_read_cut: self.source_read_cut])
        return reads[self.chosen_read][self.source_read_cut: self.target_read_cut]

def read_path_file_pairs(path_filename: str, min_gap_size: int) -> dict:
    "Read through the path file, storing the found pairs: pair -> gap estimate"
    pairs = {}
    gap_re = re.compile(r'^(\d+)N$')
    with open(path_filename, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if len(line) < 2:
                continue
            _, path = line
            path = path.split(" ")
            for idx in range(len(path) - 2):
                i, j, k = path[idx:idx+3]
                gap_match = re.search(gap_re, j)
                if gap_match and int(gap_match.group(1)) > min_gap_size:
                    # Accounting for abyss-scaffold adding 1 to gaps in path file
                    pairs[(i, k)] = PairInfo(int(gap_match.group(1)) - 1)
    return pairs

def parse_minimizers(minimizer_positions: str) -> list:
    "Parse the minimizer positions string"
    minimizer_positions = minimizer_positions.split(" ")
    return_mxs = []
    for mx_str in minimizer_positions:
        ctg, read = mx_str.split("_")
        ctg_pos, ctg_strand = ctg.split(":")
        read_pos, read_strand = read.split(":")
        return_mxs.append(MinimizerPositions(ctg_pos=int(ctg_pos), ctg_strand=ctg_strand,
                                             read_pos=int(read_pos), read_strand=read_strand))
    return return_mxs


def find_orientation(mx_positions: list):
    "Infer orientation relative to contig from minimizer positions"
    if all(x.ctg_strand == x.read_strand for x in mx_positions):
        return "+"
    if all(x.ctg_strand != x.read_strand for x in mx_positions):
        return "-"
    return None

def check_position_consistency(mx_positions: list):
    "Check that all positions in mapping are monotonically increasing or decreasing"
    if all(i.ctg_pos < j.ctg_pos for i, j in zip(mx_positions, mx_positions[1:])):
        return True
    if all(i.ctg_pos > j.ctg_pos for i, j in zip(mx_positions, mx_positions[1:])):
        return True
    return False


def reverse_complement_pair(source: str, target: str) -> tuple:
    "Reverse complement the given contig pair"
    if source[-1] == "+":
        source_new_ori = "-"
    elif source[-1] == "-":
        source_new_ori = "+"
    else:
        raise ValueError("+ or - needed for last character of node, found" + source[-1])

    if target[-1] == "+":
        target_new_ori = "-"
    elif target[-1] == "-":
        target_new_ori = "+"
    else:
        raise ValueError("+ or - needed for last character of node, found" + target[-1])

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
        if not check_position_consistency(minimizer_positions):
            continue
        read_info[read_id][ctg_id] = MinimizerMapping(anchors=int(anchors), minimizer_positions=minimizer_positions,
                                                      orientation=orientation)
        mapping_order.append(ctg_id + orientation)
        read_info[read_id]["length"] = minimizer_positions[-1].read_pos

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

def read_scaffold_file(args: argparse.Namespace) -> dict:
    "Read the scaffolds into memory"
    scaffolds = {}
    with btllib.SeqReader(args.s, btllib.SeqReaderFlag.LONG_MODE, args.t) as reads:
        for read in reads:
            scaffolds[read.id] = ScaffoldGaps(read.seq)
    return scaffolds

def calculate_est_gap_size(source_mx: MinimizerMapping, source: str,
                           target_mx: MinimizerMapping, target: str, sequences: dict, k: int) -> int:
    "Calculate the estimated gap size from the given minimizer coordinates"
    source_ctg, source_ori, target_ctg, target_ori = source[:-1], source[-1], target[:-1], target[-1]
    if source_ori == "+":
        a = sequences[source_ctg].length - source_mx.ctg_pos - k
    else:
        a = source_mx.ctg_pos
    if target_ori == "+":
        b = target_mx.ctg_pos
    else:
        b = sequences[target_ctg].length - target_mx.ctg_pos - k

    try:
        assert a >= 0
        assert b >= 0
    except ValueError as e:
        print(a, b, source, sequences[source_ctg].length, target, sequences[target_ctg].length, source_mx, target_mx)
        raise ValueError(e) from e

    gap_size = target_mx.read_pos - source_mx.read_pos - a - b
    return gap_size


def is_valid_supporting_read(source: str, target:str, read_id: str, mappings: dict,
                             sequences: dict, args: argparse.Namespace) -> bool:
    "Return true if the support read is valid. That is, infers a gao estimate < the read length"
    if source[-1] != mappings[read_id][source[:-1]].orientation:
        assert target[-1] != mappings[read_id][target[:-1]].orientation
        source, target = reverse_complement_pair(source, target)
    source_terminal_mx = mappings[read_id][source[:-1]].minimizer_positions[-1]
    target_terminal_mx = mappings[read_id][target[:-1]].minimizer_positions[0]

    gap_est = calculate_est_gap_size(source_terminal_mx, source, target_terminal_mx, target,
                                     sequences, args.large_k)
    if abs(gap_est) > mappings[read_id]["length"]:
        return False

    return True


def choose_best_read_per_pair(pairs: dict, mappings: dict, sequences: dict, args: argparse.Namespace) -> None:
    "For each pair, choose the 'best' read to fill in the gap - average anchors on each incident sequence"
    for source, target in pairs:
        reads = [(read_id, mappings[read_id][source.strip("+-")].anchors,
                  mappings[read_id][target.strip("+-")].anchors)
                 for read_id in pairs[(source, target)].mapping_reads]
        if not reads:
            continue
        sorted_reads = sorted(reads, key=lambda x: (numpy.mean([x[1], x[2]]), x[0]), reverse=True)
        for read_id, _, _ in sorted_reads:
            if is_valid_supporting_read(source, target, read_id, mappings, sequences, args):
                pairs[(source, target)].chosen_read = read_id
                break


def get_gap_fill_reads(reads_filename: str, pairs: dict, args: argparse.Namespace) -> dict:
    "Collect the reads needed for gap-filling"
    reads = {} # read_id -> sequence
    target_reads = {pairs[pair].chosen_read for pair in pairs if pairs[pair].chosen_read is not None}
    with btllib.SeqReader(reads_filename, btllib.SeqReaderFlag.LONG_MODE, args.t) as reads_in:
        for read in reads_in:
            if read.id in target_reads:
                reads[read.id] = read.seq
    return reads


# Cut positions explanations
# Situation A: ctg+ read+
#     Read orientation relative to ctg is +
#     No adjustments in cuts needed - positions in minimizers already line up
# Situation B: ctg+ read-
#     Read orientation relative to ctg is -
#     Therefore, make adjustment relative to + orientation, adjust read by k
# Situation C: ctg- read+
#     Read orientation relative to ctg is -. Given that ctg is -, means read in +ve orientation
#     Therefore, adjust relative to +ve ori read, adjust ctg cut by k
# Situation D: ctg- read-
#     Read orientation relative to ctg is +
#     No adjustments in cuts needed - positions in minimizers already line up


def assign_ctg_cut(position: int, read_est_orientation: str, ctg_orientation: str, k: int) -> int:
    "Determine adjustments needed for the cut position (if any)"
    if read_est_orientation == ctg_orientation:
        if ctg_orientation == "-":
            # Read orientation is relative to source - so read ori different from ctg, therefore "+"
            # Adjust relative to +ve orientation read
            return position + k
        return position
    return position

def assign_read_cut(position: int, read_est_orientation: str, ctg_orientation: str, k: int) -> int:
    "Determine adjustments needed for the read cut position (if any)"
    if read_est_orientation != ctg_orientation:
        if ctg_orientation == "+":
            # Read is in the "-" orientation
            return position + k
        return position
    return position


def find_masking_cut_points(pairs: dict, mappings: dict, args: argparse.Namespace) -> None:
    "Find the initial points for cuts for masking sequences for more precise cut point determination"
    for source, target in pairs:
        read_id = pairs[(source, target)].chosen_read
        if read_id is None:
            continue
        source_read_mxs = mappings[read_id][source.strip("+-")].minimizer_positions
        source_ori = source[-1]
        if mappings[read_id][source.strip("+-")].orientation == source_ori: # Read, ctg in same orientation
            source_ctg_pos, source_read_pos = source_read_mxs[-1].ctg_pos, source_read_mxs[-1].read_pos
        else:
            source_ctg_pos, source_read_pos = source_read_mxs[0].ctg_pos, source_read_mxs[0].read_pos

        target_read_mxs = mappings[read_id][target.strip("+-")].minimizer_positions
        target_ori = target[-1]
        if mappings[read_id][target.strip("+-")].orientation == target_ori:
            target_ctg_pos, target_read_pos = target_read_mxs[0].ctg_pos, target_read_mxs[0].read_pos
        else:
            target_ctg_pos, target_read_pos = target_read_mxs[-1].ctg_pos, target_read_mxs[-1].read_pos

        pairs[(source, target)].source_ctg_cut = assign_ctg_cut(source_ctg_pos,
                                                                mappings[read_id][source.strip("+-")].orientation,
                                                                source_ori, args.large_k)
        pairs[(source, target)].source_read_cut = assign_read_cut(source_read_pos,
                                                                  mappings[read_id][source.strip("+-")].orientation,
                                                                  source_ori, args.large_k)
        pairs[(source, target)].target_ctg_cut = assign_ctg_cut(target_ctg_pos,
                                                                mappings[read_id][target.strip("+-")].orientation,
                                                                target_ori, args.large_k)
        pairs[(source, target)].target_read_cut = assign_read_cut(target_read_pos,
                                                                  mappings[read_id][target.strip("+-")].orientation,
                                                                  target_ori, args.large_k)



def print_masked_sequences(scaffolds: dict, reads: dict, pairs: dict, args: argparse.Namespace):
    "Print the masked sequences to a temp file"
    out_scaffolds = open(args.s + ".masked_temp.fa", 'w')
    out_reads = open(args.reads + ".masked_temp.fa", 'w')

    for source, target in pairs:
        pair_info = pairs[(source, target)]
        if pair_info.chosen_read is None:
            continue
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
                                f"{scaffolds[target_noori].seq[target_cut:]}\n")
        else:
            raise ValueError(f"{target[-1]} must be + or -")

        read_start, read_end = min(pair_info.source_read_cut, pair_info.target_read_cut),\
                               max(pair_info.source_read_cut, pair_info.target_read_cut)

        out_reads.write(f">{pair_info.chosen_read}__{source}__{target}\n"
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
                mx_info[mx_hash] = ntlink_pair.Minimizer(entry.id, minimizer.pos,
                                                         convert_btllib_strand(minimizer.forward))
    mx_info = {mx: mx_info[mx] for mx in mx_info if mx not in dup_mxs}
    return mx_info

def map_long_reads(pairs: dict, scaffolds: dict, args: argparse.Namespace) -> None:
    "Map the long reads to the sequences, print verbose output (for now)"
    read_header_re = re.compile(r'^(\S+)__(\S+)__(\S+)$')
    scaffold_header_re = re.compile(r'^(\S+)_(source|target)$')

    with btllib.Indexlr(args.s + ".masked_temp.fa", args.k, args.w, btllib.IndexlrFlag.LONG_MODE,
                        args.t) as scaffolds_btllib:
        with btllib.Indexlr(args.reads + ".masked_temp.fa", args.k, args.w, btllib.IndexlrFlag.LONG_MODE,
                            args.t) as reads:
            for chosen_read in reads:
                _, source, target = re.search(read_header_re, chosen_read.id).groups()

                # Read source scaffold
                source_scaf = scaffolds_btllib.read()
                source_id, label = re.search(scaffold_header_re, source_scaf.id).groups()
                source_scaf.id, source_ori = source_id.strip("+-"), source_id[-1]
                assert source_id == source and label == "source"

                # Read target scaffold
                target_scaf = scaffolds_btllib.read()
                target_id, label = re.search(scaffold_header_re, target_scaf.id).groups()
                target_scaf.id, target_ori = target_id.strip("+-"), target_id[-1]
                assert target_id == target and label == "target"

                mx_info = read_btllib_minimizers([source_scaf, target_scaf])

                mxs = [(str(mx.out_hash), mx.pos, convert_btllib_strand(mx.forward))
                       for mx in chosen_read.minimizers if str(mx.out_hash) in mx_info]
                accepted_anchor_contigs, _ = \
                    ntlink_utils.get_accepted_anchor_contigs(mxs, chosen_read.readlen,
                                                             scaffolds, mx_info, args)
                if len(accepted_anchor_contigs) != 2: # Fall back on previous anchors, if option specified
                    if args.stringent:
                        pairs[(source, target)].source_read_cut = None
                        pairs[(source, target)].target_read_cut = None
                    else:
                        fallback_old_anchor_cuts(pairs, scaffolds, source, source_scaf, target, target_scaf)
                    continue

                assert len(accepted_anchor_contigs) == 2
                source_ctg_ori_read_based, source_pos_consistency, source_terminal_mx, \
                target_ctg_ori_read_based, target_pos_consistency, target_terminal_mx = \
                    assess_accepted_anchor_contigs(accepted_anchor_contigs, source_ori, source_scaf,
                                                   target_ori, target_scaf)
                if source_ctg_ori_read_based is None or target_ctg_ori_read_based is None or \
                        not source_pos_consistency or not target_pos_consistency:
                    if args.stringent:
                        pairs[(source, target)].source_read_cut = None
                        pairs[(source, target)].target_read_cut = None
                    else:
                        fallback_old_anchor_cuts(pairs, scaffolds, source, source_scaf, target, target_scaf)
                    continue

                pairs[(source, target)].source_ctg_cut = source_terminal_mx.ctg_pos
                pairs[(source, target)].source_read_cut = assign_read_cut(source_terminal_mx.read_pos,
                                                                          source_ctg_ori_read_based,
                                                                          source_ori, args.k)
                if source[-1] == "+":
                    scaffolds[source_scaf.id].three_prime_cut = assign_ctg_cut(source_terminal_mx.ctg_pos,
                                                                               source_ctg_ori_read_based,
                                                                               source_ori, args.k)
                else:
                    scaffolds[source_scaf.id].five_prime_cut = assign_ctg_cut(source_terminal_mx.ctg_pos,
                                                                              source_ctg_ori_read_based,
                                                                              source_ori, args.k)

                pairs[(source, target)].target_ctg_cut = target_terminal_mx.ctg_pos
                pairs[(source, target)].target_read_cut = assign_read_cut(target_terminal_mx.read_pos,
                                                                          target_ctg_ori_read_based,
                                                                          target_ori, args.k)
                if target[-1] == "+":
                    scaffolds[target_scaf.id].five_prime_cut = assign_ctg_cut(target_terminal_mx.ctg_pos,
                                                                              target_ctg_ori_read_based,
                                                                              target_ori, args.k)
                else:
                    scaffolds[target_scaf.id].three_prime_cut = assign_ctg_cut(target_terminal_mx.ctg_pos,
                                                                               target_ctg_ori_read_based,
                                                                               target_ori, args.k)


def assess_accepted_anchor_contigs(accepted_anchor_contigs, source_ori, source_scaf,
                                   target_ori, target_scaf):
    "Assess the accepted anchor contigs for orientation, consistency, and terminal mx"
    source_ctg_ori_read_based, source_terminal_mx = None, None
    target_terminal_mx, target_ctg_ori_read_based = None, None
    source_pos_consistency, target_pos_consistency = False, False
    for ctg_run in accepted_anchor_contigs:
        ctg_run_entry = accepted_anchor_contigs[ctg_run]
        ctg_ori_read_based = find_orientation(ctg_run_entry.hits)
        if ctg_run_entry.contig == source_scaf.id:
            if source_ori == ctg_ori_read_based:  # Read, source in same ori
                source_terminal_mx = ctg_run_entry.hits[-1]
            else:
                source_terminal_mx = ctg_run_entry.hits[0]
            source_ctg_ori_read_based = ctg_ori_read_based
            source_pos_consistency = check_position_consistency(ctg_run_entry.hits)

        if ctg_run_entry.contig == target_scaf.id:
            if target_ori == ctg_ori_read_based:  # Read, target in same ori
                target_terminal_mx = ctg_run_entry.hits[0]
            else:
                target_terminal_mx = ctg_run_entry.hits[-1]
            target_ctg_ori_read_based = ctg_ori_read_based
            target_pos_consistency = check_position_consistency(ctg_run_entry.hits)
    return source_ctg_ori_read_based, source_pos_consistency, source_terminal_mx, target_ctg_ori_read_based,\
           target_pos_consistency, target_terminal_mx


def fallback_old_anchor_cuts(pairs, scaffolds, source, source_scaf, target, target_scaf):
    "Adjust cuts to fallback to the higher k/w anchors"
    pairs[(source, target)].old_anchor_used = True
    if source[-1] == "+":
        scaffolds[source_scaf.id].three_prime_cut = pairs[(source, target)].source_ctg_cut
    else:
        scaffolds[source_scaf.id].five_prime_cut = pairs[(source, target)].source_ctg_cut
    if target[-1] == "+":
        scaffolds[target_scaf.id].five_prime_cut = pairs[(source, target)].target_ctg_cut
    else:
        scaffolds[target_scaf.id].three_prime_cut = pairs[(source, target)].target_ctg_cut


def print_gap_filled_sequences(pairs: dict, mappings: dict, sequences: dict, reads: dict,
                               args: argparse.Namespace) -> None:
    "Print out the gap-filled sequences"
    gap_re = re.compile(r'^(\d+)N$')
    outfile = open(args.o, 'w')

    gaps_counter = Counter()
    printed_scaffolds = set()
    overlap_gap = False

    with open(args.path, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if len(line) < 2:
                continue
            ctg_id, path = line
            sequence = ""
            path = path.split(" ")
            for idx, node in enumerate(path):
                gap_match = re.search(gap_re, node)
                if gap_match:
                    gap_size = int(gap_match.group(1))
                    gaps_counter["num_gaps"] += 1
                    overlap_gap = tally_small_gaps(args, gap_size, overlap_gap, gaps_counter)
                    source, target = path[idx-1], path[idx+1]
                    if (source, target) not in pairs:
                        # Accounting for gaps being one larger in abyss-scaffold path file
                        sequence += "N"*(gap_size - 1)
                        continue
                    gaps_counter["potential_fills"] += 1
                    pair_entry = pairs[(source, target)]

                    if pair_entry.source_read_cut is None or pair_entry.target_read_cut is None:
                        sequence += "N"*pair_entry.gap_size
                    elif mappings[pair_entry.chosen_read][source.strip("+-")].orientation != source[-1]:
                        if args.soft_mask:
                            sequence += pair_entry.get_cut_read_sequence(reads, "-").lower()
                        else:
                            sequence += pair_entry.get_cut_read_sequence(reads, "-")
                        gaps_counter["filled_gaps"] += 1
                        tally_anchors(pair_entry, gaps_counter)
                    else:
                        if args.soft_mask:
                            sequence += pair_entry.get_cut_read_sequence(reads, "+").lower()
                        else:
                            sequence += pair_entry.get_cut_read_sequence(reads, "+")
                        gaps_counter["filled_gaps"] += 1
                        tally_anchors(pair_entry, gaps_counter)
                else:
                    ctg = node
                    printed_scaffolds.add(ctg.strip("+-"))
                    new_sequence = sequences[ctg.strip("+-")].get_cut_sequence(ctg[-1])
                    if overlap_gap:
                        new_sequence = new_sequence[:1].lower() + new_sequence[1:]
                        overlap_gap = False
                    sequence += new_sequence
                    if args.verbose:
                        print(">{}\n{}".format(ctg, sequences[ctg.strip("+-")].get_cut_sequence(ctg[-1])),
                              file=sys.stderr)
            outfile.write(">{}\n{}\n".format(ctg_id, sequence))

    print_filling_stats(gaps_counter)

    print_unassigned_contigs(outfile, printed_scaffolds, sequences)

    outfile.close()


def print_filling_stats(counter: Counter) -> None:
    "Print statistics about gap filling"
    print("\nGap filling summary:")
    print("Number of detected sequence joins", counter["num_gaps"], sep="\t")
    print("Number of overlap sequence joins", counter["overlap_pts"], sep="\t")
    print("Number of gaps smaller than threshold", counter["small_gaps"], sep="\t")
    print("Number of potentially fillable gaps", counter["potential_fills"], sep="\t")
    print("Number of filled gaps", counter["filled_gaps"], sep="\t")
    print("Number of pass 2 anchors used", counter["new_anchor_used"], sep="\t")
    print("Number of pass 1 anchors used", counter["old_anchor_used"], sep="\t")
    print()


def print_unassigned_contigs(outfile, printed_scaffolds, sequences):
    "Print scaffolds NOT in paths"
    for ctg in sequences:
        if ctg in printed_scaffolds:
            continue
        outfile.write(">{}\n{}\n".format(ctg, sequences[ctg].seq))


def tally_anchors(pair_entry: PairInfo, counter: Counter) -> None:
    "Tally anchors used"
    if pair_entry.old_anchor_used:
        counter["old_anchor_used"] += 1
    else:
        counter["new_anchor_used"] += 1


def tally_small_gaps(args: argparse.Namespace, gap_size: int, overlap_gap: bool, counter: Counter) -> bool:
    "Tally info for small gaps"
    if gap_size == 1:  # Indicates gap size of 0 for path file
        overlap_gap = True
        counter["overlap_pts"] += 1
    if args.min_gap >= gap_size > 1:
        counter["small_gaps"] += 1
    return overlap_gap


def read_trim_coordinates(sequences: dict, args: argparse.Namespace) -> None:
    "Read in the trim coordinates from ntLink, tracking in the scaffold entry"
    with open(args.trims, 'r') as fin:
        for line in fin:
            ctg, start, end = line.strip().split("\t")
            sequences[ctg].five_prime_trim = int(start)
            sequences[ctg].three_prime_trim = int(end)

def remove_temp_masked_sequences(args: argparse.Namespace) -> None:
    "Remove temporary masked sequence files"
    out_scaffolds = args.s + ".masked_temp.fa"
    out_reads = args.reads + ".masked_temp.fa"

    shlex_cmd = shlex.split(f"rm {out_scaffolds}")
    return_code = subprocess.call(shlex_cmd)
    assert return_code == 0

    shlex_cmd = shlex.split(f"rm {out_reads}")
    return_code = subprocess.call(shlex_cmd)
    assert return_code == 0


def print_log_message(message: str) -> None:
    "Print given log message with time"
    print(datetime.datetime.today(), message, file=sys.stdout)

def print_parameters(args: argparse.Namespace) -> None:
    "Print set parameters for gap filling"
    print("Running ntLink gap-filling...\n")
    print("Parameters:")
    print("\t--path", args.path)
    print("\t--mappings", args.mappings)
    print("\t--trims", args.trims)
    print("\t-s", args.s)
    print("\t--reads", args.reads)
    print()
    print("\t-z", args.z)
    print("\t-k", args.k)
    print("\t-w", args.w)
    print("\t-t", args.t)
    print("\t--large_k", args.large_k)
    print("\t-x", args.x)
    print("\t--min_gap", args.min_gap)
    print("\t-o", args.o)
    if args.verbose:
        print("\t--verbose")
    if args.stringent:
        print("\t--stringent")
    if args.soft_mask:
        print("\t--soft_mask")

    print()

def main() -> None:
    "Run ntLink-style gap-filling"
    parser = argparse.ArgumentParser(description="Use minimizer mappings to fill gaps")
    parser.add_argument("--path", help="Input path file for gap patching", required=True, type=str)
    parser.add_argument("--mappings", help="ntLink verbose mapping TSV", required=True, type=str)
    parser.add_argument("--trims", help="ntLink file listing trims made", required=True, type=str)
    parser.add_argument("-s", help="Input scaffolds", required=True, type=str)
    parser.add_argument("--reads", help="Input reads", required=True, type=str)
    parser.add_argument("-z", help="Minimum contig size (bp) [1000]", type=int, required=False, default=1000)
    parser.add_argument("-k", help="Kmer size used in minimizer step [20]", type=int, required=False, default=20)
    parser.add_argument("-w", help="Window size used in minimizer step [10]", type=int, required=False, default=10)
    parser.add_argument("-t", help="Number of threads [4]", type=int, required=False, default=4)
    parser.add_argument("--large_k", help="K-mer size used in generating verbose mapping TSV", required=True, type=int)
    parser.add_argument("-x", help="Fudge factor allowed between mapping block lengths on read and assembly "
                                   "for re-mapping reads", type=float, required=False, default=0)
    parser.add_argument("--min_gap", help="Minimum gap size [20]", type=int, default=20)
    parser.add_argument("-o", help="Output file name", required=False, default="ntLink_patch_gaps_out.fa", type=str)
    parser.add_argument("--stringent", help="If specified, will only use lower k/w re-mapping for filling gaps,"
                                            " will not fall back on original anchors", action="store_true")
    parser.add_argument("--soft_mask", help="If specified, will soft mask the filled gap", action="store_true")
    parser.add_argument("--verbose", help="Verbose logging - print out trimmed scaffolds without gaps",
                        action="store_true")
    parser.add_argument("-v", "--version", action='version', version='ntLink v1.2.0')
    args = parser.parse_args()

    print_parameters(args)

    # Read path file into pairs ((source, target) -> (gap_est, supporting reads)
    args.min_gap = args.min_gap + 1
    pairs = read_path_file_pairs(args.path, args.min_gap)

    # Read through verbose mappings read_id -> sequence_ori -> (anchor, minimizer_pos_list)
    print_log_message("Reading ntLink read mappings..")
    mappings = read_verbose_mappings(args.mappings, pairs)

    # Read scaffold sequences into memory sequence_id -> ScaffoldGaps
    print_log_message("Reading scaffolds..")
    sequences = read_scaffold_file(args)

    # Read in trim coordinates - used when printing sequences for adjusting as needed
    print_log_message("Reading trim coordinates..")
    read_trim_coordinates(sequences, args)

    # Choose best read for patching each pair's gap (adjust pairs supporting reads)
    print_log_message("Choosing best read..")
    choose_best_read_per_pair(pairs, mappings, sequences, args)

    # Read in the reads that are needed
    print_log_message("Collecting reads...")
    reads = get_gap_fill_reads(args.reads, pairs, args)

    # Find cut points
    print_log_message("Finding masking points..")
    find_masking_cut_points(pairs, mappings, args)

    # Print masked sequences for assembly, reads for minimizer generation
    print_log_message("Printing masked sequences..")
    print_masked_sequences(sequences, reads, pairs, args)

    # Compute minimizers, and map the long read to the sequences at a lower k/w
    print_log_message("Mapping long reads..")
    map_long_reads(pairs, sequences, args)

    # Print out the sequences
    print_log_message("Printing output scaffolds..")
    print_gap_filled_sequences(pairs, mappings, sequences, reads, args)

    # Clean up masked tmp files
    print_log_message("Cleaning up..")
    remove_temp_masked_sequences(args)

    print_log_message("DONE!")


if __name__ == "__main__":
    main()
