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

MinimizerPositions = namedtuple("MinimizerEntry", ["ctg_pos", "ctg_strand", "read_pos", "read_strand"])

MinimizerMapping = namedtuple("MinimizerMapping", ["anchors", "minimizer_positions", "orientation"])

class ScaffoldGaps:
    def __init__(self, seq):
        self.seq = seq
        self.length = len(seq)
        self.five_prime_cut = 0
        self.three_prime_cut = self.length

    def __str__(self):
        return f"Length:{self.length} 5'cut:{self.five_prime_cut} 3'cut:{self.three_prime_cut}"

    def get_cut_sequence(self, reverse_compl: str):
        if reverse_compl == "-":
            return self.reverse_complement(self.seq[self.five_prime_cut: self.three_prime_cut])
        return self.seq[self.five_prime_cut: self.three_prime_cut]

    def reverse_complement(self, sequence):
        "Reverse complements a given sequence"
        translation_table = str.maketrans(
            "ACGTUNMRWSYKVHDBacgtunmrwsykvhdb",
            "TGCAANKYWSRMBDHVtgcaankywsrmbdhv")
        return sequence[::-1].translate(translation_table)

class PairInfo:
    def __init__(self, gap_size):
        self.gap_size = int(gap_size)
        self.mapping_reads = set()
        self.chosen_read = None
        self.source_ctg_cut = None
        self.source_read_cut = None
        self.target_ctg_cut = None
        self.target_read_cut = None

    def reverse_complement(self, sequence):
        "Reverse complements a given sequence"
        translation_table = str.maketrans(
            "ACGTUNMRWSYKVHDBacgtunmrwsykvhdb",
            "TGCAANKYWSRMBDHVtgcaankywsrmbdhv")
        return sequence[::-1].translate(translation_table)


    def __str__(self):
        return f"Gap: {self.gap_size}; Chosen read: {self.chosen_read}; source ctg/read cuts: {self.source_ctg_cut}/{self.source_read_cut}" \
               f"target ctg/read cuts: {self.target_ctg_cut}/{self.target_read_cut}"

    def get_cut_read_sequence(self, reads, reverse_compl: str): # !!TODO consider orientation of read
        if reverse_compl == "-":
            return self.reverse_complement(reads[self.chosen_read][self.target_read_cut: self.source_read_cut])
        return reads[self.chosen_read][self.source_read_cut: self.target_read_cut]

def read_path_file_pairs(path_filename: str, min_gap_size: int) -> dict:
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
                if gap_match and int(gap_match.group(1)) > min_gap_size:
                    pairs[(i, k)] = PairInfo(gap_match.group(1))
    return pairs

def parse_minimizers(minimizer_positions: str) -> list:
    "Parse the minimizer positions string"
    mx_pos_re = re.compile(f"MinimizerPositions\(ctg_pos=(\d+),\s+ctg_strand=\'([+-])\',\s+read_pos=(\d+),\s+"
                           f"read_strand=\'([+-])\'\)")
    return_mxs = []
    for match in re.findall(mx_pos_re, minimizer_positions):
        return_mxs.append(MinimizerPositions(ctg_pos=int(match[0]), ctg_strand=match[1],
                                             read_pos=int(match[2]), read_strand=match[3]))
    return return_mxs


def find_orientation(mx_positions: list) -> str:
    "Infer orientation from minimizer positions"
    if all(x.ctg_strand == x.read_strand for x in mx_positions):
        return "+"
    if all(x.ctg_strand != x.read_strand for x in mx_positions):
        return "-"
    return None


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
        read_info[read_id][ctg_id] = MinimizerMapping(anchors=int(anchors), minimizer_positions=minimizer_positions,
                                                      orientation=orientation)
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
    "For each pair, choose the 'best' read to fill in the gap - average anchors on each incident sequence"
    for source, target in pairs:
        reads = [(read_id, mappings[read_id][source.strip("+-")].anchors,
                  mappings[read_id][target.strip("+-")].anchors)
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

def assign_ctg_cut(position: int, read_est_orientation: str, ctg_orientation: str, k: int) -> int:
    "Determine adjustments needed for the cut position (if any)"
    if read_est_orientation == ctg_orientation:
        #Means that read orientation is "+"
        if ctg_orientation == "-":
            return position + k
        return position
    return position

def assign_read_cut(position: int, read_est_orientation: str, ctg_orientation: str, k: int) -> int:
    "Determine adjustments needed for the read cut position (if any)"
    if read_est_orientation != ctg_orientation:
        # Read is in the "-" orientation
        if ctg_orientation == "+":
            return position + k
        return position
    return position


def find_masking_cut_points(pairs: dict, mappings: dict, args: argparse.Namespace) -> None:
    "Find the initial points for cuts for masking sequences for more precise cut point determination"
    for source, target in pairs:
        read_id = pairs[(source, target)].chosen_read
        source_read_mxs = mappings[read_id][source.strip("+-")].minimizer_positions
        source_ori = source[-1]
        if mappings[read_id][source.strip("+-")].orientation == source_ori: # Read, ctg in same orientation
            source_ctg_pos, source_read_pos = source_read_mxs[-1].ctg_pos, source_read_mxs[-1].read_pos
        else:
            source_ctg_pos, source_read_pos = source_read_mxs[0].ctg_pos, source_read_mxs[0].read_pos

        target_read_mxs = mappings[read_id][target.strip("+-")].minimizer_positions
        target_ori = target[-1]
        if mappings[read_id][target.strip("+-")][2] == target_ori:
            target_ctg_pos, target_read_pos = target_read_mxs[0].ctg_pos, target_read_mxs[0].read_pos
        else:
            target_ctg_pos, target_read_pos = target_read_mxs[-1].ctg_pos, target_read_mxs[-1].read_pos

        pairs[(source, target)].source_ctg_cut = assign_ctg_cut(source_ctg_pos, mappings[read_id][source.strip("+-")].orientation,
                                                                source_ori, args.k)
        pairs[(source, target)].source_read_cut = assign_read_cut(source_read_pos, mappings[read_id][source.strip("+-")].orientation,
                                                                  source_ori, args.k)
        pairs[(source, target)].target_ctg_cut = assign_ctg_cut(target_ctg_pos, mappings[read_id][target.strip("+-")].orientation,
                                                                target_ori, args.k)
        pairs[(source, target)].target_read_cut = assign_read_cut(target_read_pos, mappings[read_id][target.strip("+-")].orientation,
                                                                  target_ori, args.k)



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
                mx_info[mx_hash] = ntlink_pair.Minimizer(entry.id, minimizer.pos, convert_btllib_strand(minimizer.forward))
    mx_info = {mx: mx_info[mx] for mx in mx_info if mx not in dup_mxs}
    return mx_info

def map_long_reads(pairs: dict, scaffolds: dict, args: argparse.Namespace) -> None:
    "Map the long reads to the sequences, print verbose output (for now)"
    read_header_re = re.compile(r'^(\S+)__(\S+)__(\S+)$')
    scaffold_header_re = re.compile(r'^(\S+)_(source|target)$')

    with btllib.Indexlr(args.s + ".masked_temp.fa", args.k, 10, btllib.IndexlrFlag.LONG_MODE) as scaffolds_btllib: # !!TODO magic numbers
        with btllib.Indexlr(args.reads + ".masked_temp.fa", args.k, 10, btllib.IndexlrFlag.LONG_MODE) as reads:
            for chosen_read in reads:
                read_id, source, target = re.search(read_header_re, chosen_read.id).groups()

                # Read source scaffold
                source_scaf = scaffolds_btllib.read()
                source_id, label = re.search(scaffold_header_re, source_scaf.id).groups()
                source_scaf.id = source_id.strip("+-")
                source_ori = source_id[-1]
                assert source_id == source
                assert label == "source"

                # Read target scaffold
                target_scaf = scaffolds_btllib.read()
                target_id, label = re.search(scaffold_header_re, target_scaf.id).groups()
                target_scaf.id = target_id.strip("+-")
                target_ori = target_id[-1]
                assert target_id == target
                assert label == "target"

                mx_info = read_btllib_minimizers([source_scaf, target_scaf])

                mxs = [(str(mx.out_hash), mx.pos, convert_btllib_strand(mx.forward)) for mx in chosen_read.minimizers if str(mx.out_hash) in mx_info]
                accepted_anchor_contigs, contig_runs = ntlink_utils.get_accepted_anchor_contigs(mxs, chosen_read.readlen,
                                                                                                scaffolds, mx_info, args)
                source_terminal_mx, source_ctg_ori_read_based = None, None
                target_terminal_mx, target_ctg_ori_read_based = None, None
                assert len(accepted_anchor_contigs) == 2
                for ctg_run in accepted_anchor_contigs:
                    ctg_run_entry = accepted_anchor_contigs[ctg_run]
                    ctg_ori_read_based = find_orientation(ctg_run_entry.hits)
                    if ctg_run_entry.contig == source_scaf.id:
                        if source_ori == ctg_ori_read_based:  # Read, source in same ori
                            source_terminal_mx = ctg_run_entry.hits[-1]
                        else:
                            source_terminal_mx = ctg_run_entry.hits[0]
                        source_ctg_ori_read_based = ctg_ori_read_based

                    if ctg_run_entry.contig == target_scaf.id:
                        if target_ori == ctg_ori_read_based:  # Read, target in same ori
                            target_terminal_mx = ctg_run_entry.hits[0]
                        else:
                            target_terminal_mx = ctg_run_entry.hits[-1]
                        target_ctg_ori_read_based = ctg_ori_read_based
                if source_ctg_ori_read_based is None or target_ctg_ori_read_based is None:
                    pairs[(source, target)].source_read_cut = None
                    pairs[(source, target)].target_read_cut = None
                    continue

                pairs[(source, target)].source_ctg_cut = source_terminal_mx.ctg_pos
                pairs[(source, target)].source_read_cut = assign_read_cut(source_terminal_mx.read_pos, source_ctg_ori_read_based,
                                                                          source_ori, args.k)
                if source[-1] == "+":
                    scaffolds[source_scaf.id].three_prime_cut = assign_ctg_cut(source_terminal_mx.ctg_pos, source_ctg_ori_read_based,
                                                                               source_ori, args.k)
                else:
                    scaffolds[source_scaf.id].five_prime_cut = assign_ctg_cut(source_terminal_mx.ctg_pos, source_ctg_ori_read_based,
                                                                              source_ori, args.k)

                pairs[(source, target)].target_ctg_cut = target_terminal_mx.ctg_pos
                pairs[(source, target)].target_read_cut = assign_read_cut(target_terminal_mx.read_pos, target_ctg_ori_read_based,
                                                                          target_ori, args.k)
                if target[-1] == "+":
                    scaffolds[target_scaf.id].five_prime_cut = assign_ctg_cut(target_terminal_mx.ctg_pos, target_ctg_ori_read_based,
                                                                              target_ori, args.k)
                else:
                    scaffolds[target_scaf.id].three_prime_cut = assign_ctg_cut(target_terminal_mx.ctg_pos, target_ctg_ori_read_based,
                                                                               target_ori, args.k)

def print_gap_filled_sequences(pairs: dict, mappings: dict, sequences: dict, reads: dict, args: argparse.Namespace) -> None:
    "Print out the gap-filled sequences"
    gap_re = re.compile('^(\d+)N$')
    outfile = open(args.o, 'w')

    num_gaps, filled_gaps = 0, 0

    with open(args.path, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if len(line) < 2:
                continue
            ctg_id, path = line
            sequence = ""
            path = path.split(" ")
            for idx in range(len(path)):
                gap_match = re.search(gap_re, path[idx])
                if gap_match:
                    num_gaps += 1
                    source, target = path[idx-1], path[idx+1]
                    if (source, target) not in pairs:
                        sequence += "N"*int(gap_match.group(1))
                        continue
                    pair_entry = pairs[(source, target)]
                    if pair_entry.source_read_cut is None or pair_entry.target_read_cut is None:
                        sequence += "N"*pair_entry.gap_size
                    elif mappings[pair_entry.chosen_read][source.strip("+-")].orientation != source[-1]:
                        sequence += pair_entry.get_cut_read_sequence(reads, "-")
                        filled_gaps += 1
                    else:
                        sequence += pair_entry.get_cut_read_sequence(reads, "+")
                        filled_gaps += 1
                else:
                    ctg = path[idx]
                    sequence += sequences[ctg.strip("+-")].get_cut_sequence(ctg[-1])
                    if args.verbose:
                        print(">{}\n{}".format(ctg, sequences[ctg.strip("+-")].get_cut_sequence(ctg[-1])), file=sys.stderr)
            outfile.write(">{}\n{}\n".format(ctg_id, sequence))
    outfile.close()

    print("Number of detected gaps", num_gaps, sep="\t")
    print("Number of filled gaps", filled_gaps, sep="\t")


def main() -> None:
    parser = argparse.ArgumentParser(description="Use minimizer mappings to fill gaps")
    parser.add_argument("--path", help="Input path file for gap patching", required=True, type=str)
    parser.add_argument("--mappings", help="ntLink verbose mapping TSV", required=True, type=str)
    parser.add_argument("-s", help="Input scaffolds", required=True, type=str)
    parser.add_argument("--reads", help="Input reads", required=True, type=str)
    parser.add_argument("-z", help="Minimum contig size (bp) [1000]", type=int, required=False, default=1000)
    parser.add_argument("-k", help="Kmer size used in minimizer step [15]", type=int, required=False, default=15)
    parser.add_argument("-x", help="Fudge factor", type=float, required=False, default=0)
    parser.add_argument("--min_gap", help="Minimum gap size [20]", type=int, default=20)
    parser.add_argument("-o", help="Output file name", required=False, default="ntLink_patch_gaps_out.fa", type=str)
    parser.add_argument("--verbose", help="Verbose logging - print out trimmed scaffolds without gaps", action="store_true")
    args = parser.parse_args()

    # Read path file into pairs ((source, target) -> (gap_est, supporting reads)
    args.min_gap = args.min_gap + 1
    pairs = read_path_file_pairs(args.path, args.min_gap)

    # Read through verbose mappings read_id -> sequence_ori -> (anchor, minimizer_pos_list)
    mappings = read_verbose_mappings(args.mappings, pairs)

    # Read scaffold sequences into memory sequence_id -> ScaffoldGaps
    print("Reading scaffolds..")
    sequences = read_scaffold_file(args.s)

    # Choose best read for patching each pair's gap (adjust pairs supporting reads)
    print("Choosing best read..")
    choose_best_read_per_pair(pairs, mappings)

    # Read in the reads that are needed
    print("Collecting reads...")
    reads = get_gap_fill_reads(args.reads, pairs)

    # Find cut points
    print("Finding cut points..")
    find_masking_cut_points(pairs, mappings, args)

    # Print masked sequences for assembly, reads for minimizer generation
    print("Printing masked sequences..")
    print_masked_sequences(sequences, reads, pairs, args)

    # Compute minimizers, and map the long read to the sequences at a lower k/w
    print("Mapping long reads..")
    map_long_reads(pairs, sequences, args)

    # Print out the sequences
    print("Printing output scaffolds..")
    print_gap_filled_sequences(pairs, mappings, sequences, reads, args)


if __name__ == "__main__":
    main()
