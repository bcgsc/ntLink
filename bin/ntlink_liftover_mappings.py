#!/usr/bin/env python3
'''
Takes the verbose mappings from ntLink, then lifts them over to the coordinate
system of the scaffolded sequences
'''
import argparse
import itertools
import re
from ntlink_utils import MinimizerPositions, reverse_orientation
from ntlink_pair import ContigRun, NtLink
import io

class AGP:
    def __init__(self, path_id, scaf_start, scaf_end, contig_id, orientation, ctg_start, ctg_end, component_id):
        self.path_id = path_id
        self.scaf_start = int(scaf_start)
        self.scaf_end = int(scaf_end)
        self.contig_id = contig_id
        self.orientation = orientation
        self.ctg_start = int(ctg_start)
        self.ctg_end = int(ctg_end)
        self.component_id = int(component_id)

    def get_ctg_length(self):
        return self.ctg_end - self.ctg_start + 1


def read_agp(agp_filename: str) -> dict:
    "Read the AGP file into a dictionary, key being contig ID and value being an AGP object"
    agp_dict = {}
    with open(agp_filename, 'r') as agp_file:
        for line in agp_file:
            line = line.strip().split('\t')
            path_id, scaf_start, scaf_end, component_id, gap_indic, ctg_id, ctg_start, ctg_end, orientation = line
            if gap_indic == "N":
                continue
            agp_dict[ctg_id] = AGP(path_id, scaf_start, scaf_end, ctg_id, orientation, ctg_start, ctg_end, component_id)
    return agp_dict

def try_int(s):
    "Try to convert a string to an int"
    try:
        return int(s)
    except ValueError:
        return s

def parse_mappings(mappings: str) -> list:
    "Parse the mappings string into a list of MinimizerPositions objects"
    mappings = mappings.split(" ")
    mappings = [MinimizerPositions(*list(map(try_int, re.split(r':|_', m)))) for m in mappings]
    return mappings

def liftover_ctg_mappings(mappings_list: list, agp_dict: dict, k: int) -> tuple:
    "Liftover the mappings for an entry"
    read_id, ctg, num_anchors, mappings = mappings_list
    mappings = parse_mappings(mappings)
    adjusted_mappings = []
    agp_entry = agp_dict[ctg]
    for m in mappings:
        if not (agp_entry.ctg_start - 1 <= m.ctg_pos <= agp_entry.ctg_end):
            continue # Mapping is outside of the assigned contig region
        adjust_pos = m.ctg_pos - (agp_entry.ctg_start - 1)
        offset = agp_entry.scaf_start - 1
        if agp_entry.orientation == '+' and agp_entry.path_id != ctg:
            new_pos = offset + adjust_pos
            adjusted_mappings.append(MinimizerPositions(new_pos, m.ctg_strand, m.read_pos, m.read_strand))
        elif agp_entry.orientation == "-" and agp_entry.path_id != ctg:
            new_pos = offset + (agp_entry.get_ctg_length() - adjust_pos) - k
            adjusted_mappings.append(MinimizerPositions(new_pos, reverse_orientation(m.ctg_strand),
                                                        m.read_pos, m.read_strand))
        else:
            adjusted_mappings.append(m)

    new_ctg_id = agp_entry.path_id
    return (read_id, new_ctg_id, num_anchors, adjusted_mappings, agp_entry)

def print_adjusted_mappings(read_id: str, mappings: list, outfile: io.TextIOWrapper) -> None:
    "Print the adjusted mapping, grouping sequences from the same path ID"
    # Group the mappings by contig, and mark subsumed
    contig_runs = [(ctg, list(tup)) for ctg, tup in itertools.groupby(mappings, lambda x: x[1])]
    contig_hits = {}
    for i, ctg_run in enumerate(contig_runs):
        ctg, list_mappings = ctg_run
        if ctg not in contig_hits:
            contig_hits[ctg] = ContigRun(ctg, i, len(list_mappings))
            contig_hits[ctg].hits = list_mappings
        else:
            for j in range(contig_hits[ctg].index + 1, i):
                contig_hits[contig_runs[j][0]].subsumed = True
            contig_hits[ctg] = ContigRun(ctg, i, len(list_mappings))
            contig_hits[ctg].hits = list_mappings

    contig_subsumed = [ctg for ctg in contig_hits if contig_hits[ctg].subsumed]
    filtered_mappings = [m for m in mappings if m[1] not in contig_subsumed]

    # Group the reads again, this time adjusting and printing out the mappings
    for ctg, tup in itertools.groupby(filtered_mappings, lambda x: x[1]):
        tup = list(tup)

        concat_mappings = [m for run in tup for m in run[3]]
        if not concat_mappings:
            continue # Don't print if empty list
        monotonic_increase = all(i.ctg_pos < j.ctg_pos for i, j in zip(concat_mappings, concat_mappings[1:]))
        monotonic_decrease = all(i.ctg_pos > j.ctg_pos for i, j in zip(concat_mappings, concat_mappings[1:]))
        if not monotonic_increase and not monotonic_decrease: #!!TODO look into this more?
            continue
#        assert monotonic_increase or monotonic_decrease, (ctg, read_id, concat_mappings)
        mx_string = NtLink.print_minimizer_positions(concat_mappings)
        outfile.write(f"{read_id}\t{ctg}\t{len(concat_mappings)}\t{mx_string}\n")



def liftover_mappings(mappings_filename: str, agp_dict: dict, prefix: str, k: int) -> None:
    "Liftover the verbose mappings file"
    cur_read_id = None
    cur_read_id_mappings = []
    with open(mappings_filename, 'r') as mappings_file:
        with open(prefix + '.liftover.mappings', 'w') as liftover_mappings_file:
            for line in mappings_file:
                line = line.strip().split('\t')
                mapping = liftover_ctg_mappings(line, agp_dict, k)
                if mapping[0] != cur_read_id:
                    # Deal with previous mappings
                    if cur_read_id is not None:
                        print_adjusted_mappings(cur_read_id, cur_read_id_mappings, liftover_mappings_file)
                    cur_read_id = mapping[0]
                    cur_read_id_mappings = [mapping]
                else:
                    cur_read_id_mappings.append(mapping)
            # Don't forget to deal with the last read
            print_adjusted_mappings(cur_read_id, cur_read_id_mappings, liftover_mappings_file)


def main() -> None:
    parser = argparse.ArgumentParser(description='Liftover of ntLink mappings')
    parser.add_argument("-m", "--mappings", help="Path to the verbose mappings file", required=True)
    parser.add_argument("-a", "--agp", help="Path to the AGP file", required=True)
    parser.add_argument("-p", "--prefix", help="Prefix for the output files", required=True)
    parser.add_argument("-k", "--kmer", help="Kmer size", required=True, type=int)
    args = parser.parse_args()

    # Read in the AGP file
    agp_dict = read_agp(args.agp)

    # Liftover the verbose mappings file
    liftover_mappings(args.mappings, agp_dict, args.prefix, args.kmer)

if __name__ == "__main__":
    main()
