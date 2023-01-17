#!/usr/bin/env python3
'''
Takes the verbose mappings from ntLink, then lifts them over to the coordinate
system of the scaffolded sequences
'''
import argparse
import itertools
import io
from collections import namedtuple
from ntlink_utils import MinimizerPositions, reverse_orientation
from ntlink_pair import ContigRun, NtLink

NewMinimizerMapping = namedtuple("NewMinimizerMapping", ["read_id", "new_ctg_id", "num_anchors",
                                                         "adjusted_mappings", "agp_entry"])

class AGP:
    "Represents an AGP file entry"
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
        "Returns the length of the contig"
        return self.ctg_end - self.ctg_start + 1

    def get_scaf_length(self):
        "Returns the length of the scaffold"
        return self.scaf_end - self.scaf_start + 1


def read_agp(agp_filename: str) -> dict:
    "Read the AGP file into a dictionary, key being contig ID and value being an AGP object"
    agp_dict = {}
    with open(agp_filename, 'r') as agp_file:
        for line in agp_file:
            line = line.strip().split('\t')
            path_id, scaf_start, scaf_end, component_id, component_type, ctg_id, ctg_start, ctg_end, orientation = line
            if component_type in ("N", "P"):
                continue
            agp_dict[ctg_id] = AGP(path_id, scaf_start, scaf_end, ctg_id, orientation, ctg_start, ctg_end, component_id)
    return agp_dict


def parse_mappings(mappings: str) -> list:
    "Parse the mappings string into a list of MinimizerPositions objects"
    for m in mappings.split(" "):
        ctg_maps, read_maps = m.split("_")
        ctg_pos, ctg_strand = ctg_maps.split(":")
        read_pos, read_strand = read_maps.split(":")
        yield MinimizerPositions(None, int(ctg_pos), ctg_strand, int(read_pos), read_strand)

def liftover_ctg_mappings(mappings_list: list, agp_dict: dict, k: int) -> NewMinimizerMapping:
    "Liftover the mappings for an entry"
    read_id, ctg, num_anchors, mappings = mappings_list
    adjusted_mappings = []
    if ctg not in agp_dict:
        return NewMinimizerMapping(read_id, ctg, 0, [], None)
    agp_entry = agp_dict[ctg]
    for m in parse_mappings(mappings):
        if not agp_entry.ctg_start - 1 <= m.ctg_pos <= (agp_entry.ctg_end - k):
            continue # Mapping is outside of the assigned contig region
        adjust_pos = m.ctg_pos - (agp_entry.ctg_start - 1)
        offset = agp_entry.scaf_start - 1
        if agp_entry.orientation == '+' and agp_entry.path_id != ctg:
            new_pos = offset + adjust_pos
            adjusted_mappings.append(MinimizerPositions(None, new_pos, m.ctg_strand, m.read_pos, m.read_strand))
        elif agp_entry.orientation == "-" and agp_entry.path_id != ctg:
            new_pos = offset + (agp_entry.get_ctg_length() - adjust_pos) - k
            adjusted_mappings.append(MinimizerPositions(None, new_pos, reverse_orientation(m.ctg_strand),
                                                        m.read_pos, m.read_strand))
        else:
            adjusted_mappings.append(m)

    new_ctg_id = agp_entry.path_id
    return NewMinimizerMapping(read_id, new_ctg_id, num_anchors, adjusted_mappings, agp_entry)

def print_adjusted_mappings(read_id: str, mappings: list, outfile: io.TextIOWrapper) -> None:
    "Print the adjusted mapping, grouping sequences from the same path ID"
    # Group the mappings by contig, and mark subsumed
    contig_runs = [(ctg, list(tup)) for ctg, tup in itertools.groupby(mappings, lambda x: x.new_ctg_id)]
    contig_hits = {}
    for i, ctg_run in enumerate(contig_runs):
        ctg, list_mappings = ctg_run
        if ctg not in contig_hits:
            contig_hits[ctg] = ContigRun(ctg, list_mappings)
            contig_hits[ctg].index = i
        else:
            for j in range(contig_hits[ctg].index + 1, i):
                contig_hits[contig_runs[j][0]].subsumed = True
            contig_hits[ctg].hits.extend(list_mappings)
            contig_hits[ctg].hit_count += len(list_mappings)

    filtered_mappings = [m for m in mappings if not contig_hits[m.new_ctg_id].subsumed]

    # Group the reads again, this time adjusting and printing out the mappings
    for ctg, tup in itertools.groupby(filtered_mappings, lambda x: x.new_ctg_id):
        tup = list(tup)

        concat_mappings = [m for run in tup for m in run.adjusted_mappings]
        if not concat_mappings:
            continue # Don't print if empty list
        monotonic_increase = all(i.ctg_pos < j.ctg_pos for i, j in zip(concat_mappings, concat_mappings[1:]))
        if not monotonic_increase:
            monotonic_decrease = all(i.ctg_pos > j.ctg_pos for i, j in zip(concat_mappings, concat_mappings[1:]))
            if not monotonic_decrease:
                continue # !! TODO: deal with these cases?
        mx_string = NtLink.print_minimizer_positions(concat_mappings)
        outfile.write(f"{read_id}\t{ctg}\t{len(concat_mappings)}\t{mx_string}\n")



def liftover_mappings(mappings_filename: str, agp_dict: dict, output: str, k: int) -> None:
    "Liftover the verbose mappings file"
    cur_read_id = None
    cur_read_id_mappings = []
    with open(mappings_filename, 'r') as mappings_file:
        with open(output, 'w') as liftover_mappings_file:
            for line in mappings_file:
                line = line.strip().split('\t')
                mapping = liftover_ctg_mappings(line, agp_dict, k)
                if mapping.read_id != cur_read_id:
                    # Deal with previous mappings
                    if cur_read_id is not None:
                        print_adjusted_mappings(cur_read_id, cur_read_id_mappings, liftover_mappings_file)
                    cur_read_id = mapping.read_id
                    cur_read_id_mappings = [mapping]
                else:
                    cur_read_id_mappings.append(mapping)
            # Don't forget to deal with the last read
            print_adjusted_mappings(cur_read_id, cur_read_id_mappings, liftover_mappings_file)


def main() -> None:
    "Liftover the ntLink verbose mappings file"
    parser = argparse.ArgumentParser(description='Liftover of ntLink mappings')
    parser.add_argument("-m", "--mappings", help="Path to the verbose mappings file", required=True)
    parser.add_argument("-a", "--agp", help="Path to the AGP file", required=True)
    parser.add_argument("-o", "--output", help="Output file name", required=True)
    parser.add_argument("-k", "--kmer", help="Kmer size", required=True, type=int)
    parser.add_argument("-v", "--version", action='version', version='ntLink v1.3.8')
    args = parser.parse_args()

    # Read in the AGP file
    agp_dict = read_agp(args.agp)

    # Liftover the verbose mappings file
    liftover_mappings(args.mappings, agp_dict, args.output, args.kmer)

if __name__ == "__main__":
    main()
