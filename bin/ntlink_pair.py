#!/usr/bin/env python3
"""
Finding contig pairs using lightweight long read mapping with minimizers
"""
__author__ = 'laurencoombe'

import argparse
import datetime
from collections import defaultdict, namedtuple, Counter
import itertools
import os
import random
import re
import shlex
import subprocess
import sys
import numpy as np
import igraph as ig
import ntlink_utils

MinimizerEdge = namedtuple("MinimizerEdge", ["mx_i", "mx_i_pos", "mx_i_strand",
                                             "mx_j", "mx_j_pos", "mx_j_strand"])
Minimizer = namedtuple("Minimizer", ["contig", "position", "strand"])
MinimizerWithHash = namedtuple("Minimizer_with_hash", ["mx_hash", "contig", "position", "strand"])
MappingEntry = namedtuple("MappingEntry", ["read_id", "contig_id", "num_hits", "list_hits"])

class NtlinkPairError(Exception):
    pass


class ScaffoldPair:
    "Object that represents a scaffold pair"
    def __init__(self, source_contig, source_ori, target_contig, target_ori):
        self.source_contig = source_contig
        self.source_ori = source_ori
        self.target_contig = target_contig
        self.target_ori = target_ori

    def format_source(self):
        "Format the source contig with orientation"
        return self.source_contig + self.source_ori

    def format_target(self):
        "Format the target contig with orientation"
        return self.target_contig + self.target_ori

    def __str__(self):
        return "Source: {source} -> Target: {target}".format(source=self.format_source(), target=self.format_target())

    def __eq__(self, other):
        return other.format_source() == self.format_source() and other.format_target() == self.format_target()

    def __hash__(self):
        return hash((self.source_contig, self.source_ori, self.target_contig, self.target_ori))

class PairInfo:
    "Represents more information for a scaffold pair"
    def __init__(self):
        self._gap_est = []
        self.anchor = 0
        self._gap_estimate = None

    def add_gap_estimate(self, gap_est):
        "Add a gap estimate"
        self._gap_est.append(gap_est)
        self._gap_estimate = None

    def get_gap_estimate(self):
        "Calculate gap estimate"
        if self._gap_estimate is None:
            self._gap_estimate = int(np.median(self._gap_est))
        return self._gap_estimate

    def n_supporting_reads(self):
        "Return the number of supporting reads for a pair"
        return len(self._gap_est)

    def __str__(self):
        return "n={n}, gap_estimates={gap_ests}, anchor={a}".format(n=self.n_supporting_reads(),
                                                                    gap_ests=self._gap_est,
                                                                    a=self.anchor)

class ContigRun:
    "Represents information about a contig run based on a long read"
    def __init__(self, contig, list_hits):
        self.contig = contig
        self.index = None
        self.hit_count = len(list_hits)
        self.hits = list_hits
        self._subsumed = False
        self.first_mx = None
        self.terminal_mx = None

    @property
    def subsumed(self):
        "Set the subsumed field"
        return self._subsumed

    @subsumed.setter
    def subsumed(self, bool_val):
        if isinstance(bool_val, bool):
            self._subsumed = bool_val
        else:
            raise ValueError("subsumed must be Boolean")

    def __str__(self):
        return "contig={contig}, hit_count={hit_count}, subsumed={subsumed}, " \
               "index={index}, hits={hits}, first_mx={first_mx}, terminal_mx={terminal_mx}".format(
            contig=self.contig, hit_count=self.hit_count,
            subsumed=self.subsumed, index=self.index, hits=self.hits, first_mx=self.first_mx,
            terminal_mx=self.terminal_mx)

class NtLink():
    "Represents an ntLink graph construction run"

    @staticmethod
    def get_largest_ntlink_scaffold_id(scaffolds):
        "Detect if any headers adhere to ntLink_{num} convention, and if so return the largest number"
        header_re = re.compile(r'^ntLink_(\d+)$')
        largest_num = None

        for scaffold in scaffolds:
            header_match = re.search(header_re, scaffold)
            if header_match:
                largest_num = int(header_match.group(1)) \
                    if largest_num is None or int(header_match.group(1)) > largest_num \
                    else largest_num

        return largest_num

    @staticmethod
    def print_directed_graph(graph, out_prefix, scaffolds):
        "Prints the directed scaffold graph in dot format"
        out_graph = out_prefix + ".scaffold.dot"
        outfile = open(out_graph, 'w')
        print(datetime.datetime.today(), ": Printing graph", out_graph, sep=" ", file=sys.stdout)

        outfile.write("digraph G {\n")
        outfile.write("graph [scaf_num={}]\n".format(NtLink.get_largest_ntlink_scaffold_id(scaffolds)))

        for node in graph.vs():
            node_label = "\"{scaffold}\" [l={length}]\n".\
                format(scaffold=node['name'], length=scaffolds[node['name'][:-1]].length)
            outfile.write(node_label)

        for edge in graph.es():
            edge_str = "\"{source}\" -> \"{target}\" [d={d} e={e} n={n}]\n".\
                format(source=ntlink_utils.vertex_name(graph, edge.source),
                       target=ntlink_utils.vertex_name(graph, edge.target),
                       d=int(edge['d']), e=edge['e'], n=edge['n'])
            outfile.write(edge_str)

        outfile.write("}\n")

    def calculate_gap_size(self, i_mx, i_ori, j_mx, j_ori, est_distance):
        "Calculates the estimated distance between two contigs"

        # Correct for the overhanging sequence before/after terminal minimizers
        if i_ori == "+":
            u_ctg = NtLink.list_mx_info[i_mx].contig
            u_ctglen = NtLink.scaffolds[u_ctg].length
            a = u_ctglen - NtLink.list_mx_info[i_mx].position - self.args.k
        else:
            a = NtLink.list_mx_info[i_mx].position
        if j_ori == "+":
            b = NtLink.list_mx_info[j_mx].position
        else:
            v_ctg = NtLink.list_mx_info[j_mx].contig
            v_ctglen = NtLink.scaffolds[v_ctg].length
            b = v_ctglen - NtLink.list_mx_info[j_mx].position - self.args.k
        try:
            assert a >= 0
            assert b >= 0
        except AssertionError as assert_error:
            print("ERROR: Gap distance estimation less than 0", "Vertex 1:", i_mx, "Vertex 2:", j_mx,
                  sep="\n")
            print("Minimizer positions:", NtLink.list_mx_info[i_mx].position,
                  NtLink.list_mx_info[j_mx].position)
            print("Estimated distance: ", est_distance)
            raise assert_error

        gap_size = est_distance - a - b
        return int(gap_size)

    def read_minimizers(self):
        "Read the minimizers from a file, removing duplicate minimizers"
        print(datetime.datetime.today(), ": Reading minimizers", self.args.s, file=sys.stdout)
        mx_info = {}  # mx -> Minimizer object
        dup_mxs = set()  # Set of minimizers identified as duplicates

        inmx_file = "/dev/stdin" if self.args.m == "-" else self.args.m

        with open(inmx_file, 'r') as fin:
            for record in fin:
                record = record.strip().split("\t")
                if len(record) > 1:
                    ctg_name = record[0]
                    for mx_pos_strand in record[1].split(" "):
                        mx, pos, strand = mx_pos_strand.split(":")
                        if mx in mx_info:  # This is a duplicate, add to dup set, don't add to dict
                            dup_mxs.add(mx)
                        else:
                            mx_info[mx] = Minimizer(ctg_name, int(pos), strand)

        mx_info = {mx: mx_info[mx] for mx in mx_info if mx not in dup_mxs}

        return mx_info

    @staticmethod
    def normalize_pair(source_ctg, source_ori, target_ctg, target_ori):
        "Normalize pairs. Lexicographically smallest contig first in pair"
        if source_ctg < target_ctg:
            return ScaffoldPair(source_ctg, source_ori, target_ctg, target_ori)
        return ScaffoldPair(target_ctg, ntlink_utils.reverse_orientation(target_ori),
                            source_ctg, ntlink_utils.reverse_orientation(source_ori))


    def calculate_pair_info(self, mx_edge):
        "Given a contig pair, normalizes, defines orientation and estimates the gap size"

        assert mx_edge.mx_i_pos < mx_edge.mx_j_pos
        source_ctg = NtLink.list_mx_info[mx_edge.mx_i].contig
        target_ctg = NtLink.list_mx_info[mx_edge.mx_j].contig
        if mx_edge.mx_i_strand == NtLink.list_mx_info[mx_edge.mx_i].strand:
            source_ori = "+"
        else:
            source_ori = "-"
        if mx_edge.mx_j_strand == NtLink.list_mx_info[mx_edge.mx_j].strand:
            target_ori = "+"
        else:
            target_ori = "-"
        new_pair = self.normalize_pair(source_ctg, source_ori, target_ctg, target_ori)
        gap_estimate = self.calculate_gap_size(mx_edge.mx_i, source_ori, mx_edge.mx_j, target_ori,
                                               mx_edge.mx_j_pos - mx_edge.mx_i_pos)
        return new_pair, gap_estimate

    def filter_weak_anchor_pairs(self, pairs):
        "Filter out edges where there isn't at least threshold # well-anchored reads"
        new_pair = {pair: pairs[pair] for pair in pairs if pairs[pair].anchor >= self.args.a}
        return new_pair

    @staticmethod
    def filter_pairs_distances(pairs):
        "Filter out edges where distance estimate is incongruous with either scaffold"
        new_pairs = {}
        for pair in pairs:
            if pairs[pair].get_gap_estimate() <= NtLink.scaffolds[pair.source_contig].length*-1 or \
                    pairs[pair].get_gap_estimate() <= NtLink.scaffolds[pair.target_contig].length*-1:
                continue
            new_pairs[pair] = pairs[pair]
        return new_pairs

    @staticmethod
    def reverse_complement_pair(pair):
        "Reverse complemented a directed pair of scaffolds"
        return ScaffoldPair(pair.target_contig, ntlink_utils.reverse_orientation(pair.target_ori),
                            pair.source_contig, ntlink_utils.reverse_orientation(pair.source_ori))

    def build_scaffold_graph(self, pairs):
        "Builds a scaffold graph given the pairs info"
        print(datetime.datetime.today(), ": Building scaffold graph", file=sys.stdout)

        graph = ig.Graph(directed=True)

        vertices = set()
        edges = defaultdict(dict) # source -> target -> [gap distance estimates]

        for pair in pairs:
            reversed_pair = self.reverse_complement_pair(pair)

            vertices.add(pair.format_source())
            vertices.add(pair.format_target())
            vertices.add(reversed_pair.format_source())
            vertices.add(reversed_pair.format_target())

            try:
                assert not (pair.format_source() in edges and pair.format_target() in edges[pair.format_source()])
                assert not (reversed_pair.format_source() in edges and
                            reversed_pair.format_target() in edges[reversed_pair.format_source()])
            except AssertionError:
                print(pair)
                print(reversed_pair)
                print(edges[pair.format_source()])
                sys.exit(1)

            edges[pair.format_source()][pair.format_target()] = pairs[pair]
            edges[reversed_pair.format_source()][reversed_pair.format_target()] = pairs[pair]

        formatted_edges = [(s, t) for s in edges for t in edges[s]]
        graph.add_vertices(list(vertices))
        graph.add_edges(formatted_edges)

        edge_attributes = {ntlink_utils.edge_index(graph, s, t): {'d': edges[s][t].get_gap_estimate(),
                                                                  "e": 100,
                                                                  "n": edges[s][t].n_supporting_reads()}
                           for s in edges for t in edges[s]}
        graph.es()["d"] = [edge_attributes[e]['d'] for e in sorted(edge_attributes.keys())]
        graph.es()["e"] = [edge_attributes[e]['e'] for e in sorted(edge_attributes.keys())]
        graph.es()["n"] = [edge_attributes[e]['n'] for e in sorted(edge_attributes.keys())]

        return graph

    @staticmethod
    def print_minimizer_positions(list_minimizers):
        "Print information about minimizer positions/strands on ctg/read in brief format"
        return_list = []
        for mx in list_minimizers:
            return_list.append(f"{mx.ctg_pos}:{mx.ctg_strand}_{mx.read_pos}:{mx.read_strand}")
        return " ".join(return_list)

    def add_pair(self, accepted_anchor_contigs, ctg_i, ctg_j, pairs, length_read, check_added=None):
        "Add pair to dictionary of pairs"
        mx_i = accepted_anchor_contigs[ctg_i].terminal_mx
        mx_j = accepted_anchor_contigs[ctg_j].first_mx
        pair, gap_est = self.calculate_pair_info(MinimizerEdge(mx_i.mx_hash, mx_i.position, mx_i.strand,
                                                               mx_j.mx_hash, mx_j.position, mx_j.strand))

        if abs(gap_est) > length_read:
            return None
        if check_added is not None and pair in check_added:
            return None

        if pair not in pairs:
            pairs[pair] = PairInfo()
        pairs[pair].add_gap_estimate(gap_est)
        if accepted_anchor_contigs[ctg_i].hit_count > 1 and \
                        accepted_anchor_contigs[ctg_j].hit_count > 1:
            pairs[pair].anchor += 1

        return pair

    def find_scaffold_pairs(self):
        "Builds up pairing information between scaffolds"
        print(datetime.datetime.today(), ": Finding pairs", file=sys.stdout)

        target_mxs = set(NtLink.list_mx_info.keys())

        pairs = {} # source -> target -> [gap estimate]

        # Open file for outputting verbose logging if option specified
        verbose_file = None
        if self.args.verbose:
            verbose_file = open(self.args.p + ".verbose_mapping.tsv", 'w')
        if self.args.paf:
            paf_file = open(self.args.p + ".paf", 'w')

        # Add the long read edges to the graph
        for mx_long_file in self.args.FILES:
            mx_long_filename = "/dev/stdin" if mx_long_file == "-" else mx_long_file
            with open(mx_long_filename, 'r') as long_mxs:
                for line in long_mxs:
                    line = line.strip().split("\t")
                    if len(line) < 3:
                        continue
                    mx_pos_split_tups = line[2].split(" ")
                    mx_pos_split = []
                    mx_seen = set()
                    mx_dups = set()

                    for mx_pos in mx_pos_split_tups:
                        mx, pos, strand = mx_pos.split(":")
                        if mx in target_mxs:
                            mx_pos_split.append((mx, int(pos), strand))
                            if self.args.repeat_filter:
                                if mx in mx_seen:
                                    mx_dups.add(mx)
                                else:
                                    mx_seen.add(mx)
                    if self.args.repeat_filter:
                        mx_pos_split = [(mx, pos, strand) for mx, pos, strand in mx_pos_split if mx not in mx_dups]
                    if not mx_pos_split:
                        continue
                    read_name = line[0]
                    length_long_read = int(line[1])
                    accepted_anchor_contigs, contig_runs = \
                        ntlink_utils.get_accepted_anchor_contigs(mx_pos_split,length_long_read,
                                                                 NtLink.scaffolds, NtLink.list_mx_info, self.args)
                    if self.args.verbose and accepted_anchor_contigs:
                        for ctg_run in accepted_anchor_contigs:
                            verbose_file.write("{}\t{}\t{}\t{}\n".
                                               format(read_name, accepted_anchor_contigs[ctg_run].contig,
                                                      accepted_anchor_contigs[ctg_run].hit_count,
                                                      self.print_minimizer_positions(
                                                          accepted_anchor_contigs[ctg_run].hits)))
                    if self.args.paf and accepted_anchor_contigs:
                        self.print_paf(paf_file, accepted_anchor_contigs, length_long_read, read_name)


                    # Set first and terminal minimizers for the hits
                    for contig in accepted_anchor_contigs:
                        contig_run = accepted_anchor_contigs[contig].hits
                        first_mx = contig_run[0]
                        accepted_anchor_contigs[contig].first_mx = MinimizerWithHash(first_mx.mx,
                                                                                     contig,
                                                                                     first_mx.read_pos,
                                                                                     first_mx.read_strand)
                        last_mx = contig_run[-1]
                        accepted_anchor_contigs[contig].terminal_mx = MinimizerWithHash(last_mx.mx,
                                                                                        contig,
                                                                                        last_mx.read_pos,
                                                                                        last_mx.read_strand)

                    self.tally_pairs_from_mappings(accepted_anchor_contigs, contig_runs, length_long_read, pairs)
        if self.args.verbose:
            verbose_file.close()
        if self.args.paf:
            paf_file.close()

        return pairs

    @staticmethod
    def is_consistent(ctg_pos: list, increasing: bool, i1: int, i2: int) -> bool:
        "Given the expected monotonicity and a list with contig positions of mappings, determine whether the given indices are consistent"
        if increasing:
            return ctg_pos[i1].read_pos <= ctg_pos[i2].read_pos
        else:
            return ctg_pos[i2].read_pos >= ctg_pos[i2].read_pos

    @staticmethod
    def break_alignment_blocks(sorted_ctg_pos, breaks, filters):
        "Break the alignment blocks at detected locations"
        return_alignment_blocks = []
        current_alignment_block = []
        for i, mapping in enumerate(sorted_ctg_pos):
            if i in filters:
                continue
            if i in breaks:
                return_alignment_blocks.append(current_alignment_block)
                current_alignment_block = [mapping]
            else:
                current_alignment_block.append(mapping)
        if current_alignment_block:
            return_alignment_blocks.append(current_alignment_block)

        return return_alignment_blocks


    @staticmethod
    def filter_and_break_alignment_blocks(transitions, sorted_ctg_pos, duplicate_positions, increasing=True):
        breaks = set()
        filters = set()
        for i, transition in enumerate(transitions):
            if not transition:
                if sorted_ctg_pos[i].ctg_pos in duplicate_positions or sorted_ctg_pos[i+1].ctg_pos in duplicate_positions:
                    continue # Just skip when transitions include duplicate contig positions
                if i + 2 >= len(transitions):
                    # This is an end minimizer that's an issue. Remove it
                    breaks.add(i + 1)
                elif NtLink.is_consistent(sorted_ctg_pos, increasing, i, i + 2):
                    # This is a single issue minimizer. Remove it
                    filters.add(i + 1)
                elif i > 0 and NtLink.is_consistent(sorted_ctg_pos, increasing, i - 1, i + 1):
                    # This is a single issue minimizer. Remove it
                    filters.add(i)
                else:
                    # This is a larger segment problem or problem minimizer at the beginning. Break the alignment block
                    breaks.add(i + 1)

        if not breaks and not filters:
            return [sorted_ctg_pos]
        return NtLink.break_alignment_blocks(sorted_ctg_pos, breaks, filters)

    @staticmethod
    def get_mapped_blocks(sorted_ctg_pos: list, min_consistent=0.75) -> list:
        "Go through a list of transitions, deciding whether to filter or cut the alignment block accordingly"
        ctg_positions = set()
        dup_positions = set()
        transitions_incr, transitions_decr = [], []

        for i, j in zip(sorted_ctg_pos, sorted_ctg_pos[1:]):
            transitions_incr.append(i.read_pos <= j.read_pos)
            transitions_decr.append(i.read_pos >= j.read_pos)
            if i.ctg_pos in ctg_positions:
                dup_positions.add(i.ctg_pos)
            else:
                ctg_positions.add(i.ctg_pos)
        if sorted_ctg_pos[-1].ctg_pos in ctg_positions:
            dup_positions.add(sorted_ctg_pos[-1].ctg_pos)

        if all(transitions_incr):
            return [sorted_ctg_pos]
        if all(transitions_decr):
            return [sorted_ctg_pos]

        transition_counts = Counter(transitions_incr)
        if (transition_counts[True]/len(transitions_incr)) >= min_consistent:
            return NtLink.filter_and_break_alignment_blocks(transitions_incr, sorted_ctg_pos, dup_positions, increasing=True)
        if (transition_counts[False]/len(transitions_incr)) >= min_consistent:
            return NtLink.filter_and_break_alignment_blocks(transitions_decr, sorted_ctg_pos, dup_positions, increasing=False)
        return []


    def print_paf(self, outfile, accepted_contigs, read_len, read_name):
        "Print the given read mappings in PAF-like format"
        for ctg in accepted_contigs:
            ctg_run = accepted_contigs[ctg]
            sorted_mx_positions = sorted(ctg_run.hits, key=lambda x:x.ctg_pos)
            mapped_blocks = self.get_mapped_blocks(sorted_mx_positions)
            for mapping in mapped_blocks:
                first_mx_mapping = mapping[0]
                last_mx_mapping = mapping[-1]
                strand_counter = Counter([hit.ctg_strand == hit.read_strand for hit in mapping])
                if strand_counter[True]/len(strand_counter)*100 >= 50:
                    strand = "+"
                else:
                    strand = "-"
                target_start, target_end = min(first_mx_mapping.ctg_pos, last_mx_mapping.ctg_pos), \
                                           max(first_mx_mapping.ctg_pos, last_mx_mapping.ctg_pos) + self.args.k
                query_start, query_end = min(first_mx_mapping.read_pos, last_mx_mapping.read_pos), \
                                         max(first_mx_mapping.read_pos, last_mx_mapping.read_pos) + self.args.k
                assert query_start < query_end
                assert query_start >= 0
                assert query_end <= read_len

                out_str = f"{read_name}\t{read_len}\t{query_start}\t{query_end}\t{strand}\t" \
                          f"{ctg_run.contig}\t{NtLink.scaffolds[ctg_run.contig].length}\t" \
                          f"{target_start}\t{target_end}\t{len(mapping)}\t" \
                          f"{target_end - target_start}\t255\n"
                outfile.write(out_str)


    def tally_pairs_from_mappings(self, accepted_anchor_contigs, contig_runs, length_long_read, pairs):
        "Tally the pairs from the given mappings"
        if len(contig_runs) <= self.args.f:
            # Add all transitive edges for pairs
            for ctg_pair in itertools.combinations(contig_runs, 2):
                ctg_i, ctg_j = ctg_pair
                self.add_pair(accepted_anchor_contigs, ctg_i, ctg_j, pairs, length_long_read)
        else:
            added_pairs = set()
            # Add adjacent pairs
            for ctg_i, ctg_j in zip(contig_runs, contig_runs[1:]):
                new_pair = self.add_pair(accepted_anchor_contigs, ctg_i, ctg_j, pairs, length_long_read)
                added_pairs.add(new_pair)

            # Add transitive edges over weakly supported contigs
            contig_runs_filter = [ctg for ctg in contig_runs
                                  if accepted_anchor_contigs[ctg].hit_count > 1]
            for ctg_i, ctg_j in zip(contig_runs_filter, contig_runs_filter[1:]):
                self.add_pair(accepted_anchor_contigs, ctg_i, ctg_j, pairs, length_long_read,
                              check_added=added_pairs)

    def find_scaffold_pairs_checkpoints(self):
        "Build up pairing information between scaffold - using checkpoint mapping file"
        print(datetime.datetime.today(), ": Finding pairs", file=sys.stdout)
        pairs = {}

        # Read through checkpoint file
        curr_read_id = None
        curr_mappings = []
        with open(self.args.checkpoint, 'r') as checkpoint_file:
            for line in checkpoint_file:
                read_id, contig_id, num_hits, mx_hits = line.strip().split('\t')
                if read_id != curr_read_id:
                    if curr_read_id is not None:
                        self.parse_verbose_entries(curr_mappings, pairs)
                    curr_read_id = read_id
                    curr_mappings = [MappingEntry(read_id, contig_id, num_hits, mx_hits)]
                else:
                    curr_mappings.append(MappingEntry(read_id, contig_id, num_hits, mx_hits))
        # Don't forget the last read mapping
        self.parse_verbose_entries(curr_mappings, pairs)

        return pairs

    def parse_verbose_entries(self, mappings, pairs):
        "Parse the verbose entries from the given read, convert to accepted contigs runs"
        accepted_anchor_contigs = {}
        contig_runs = []
        read_mapping_positions = []
        for i, m in enumerate(mappings):
            contig_runs.append(m.contig_id)
            accepted_anchor_contigs[m.contig_id] = ContigRun(m.contig_id, ntlink_utils.parse_minimizers(m.list_hits))
            accepted_anchor_contigs[m.contig_id].index = i
            # Generate random 64-bit integer to represent the minimizer hash
            first_hash = random.getrandbits(64)
            last_hash = random.getrandbits(64)
            # Easier access to the first and last minimizer hits
            first_mx_hit = accepted_anchor_contigs[m.contig_id].hits[0]
            last_mx_hit = accepted_anchor_contigs[m.contig_id].hits[-1]
            # Load minimizer positions for the read
            accepted_anchor_contigs[m.contig_id].first_mx = MinimizerWithHash(first_hash, m.contig_id,
                                                                              first_mx_hit.read_pos,
                                                                              first_mx_hit.read_strand)
            accepted_anchor_contigs[m.contig_id].terminal_mx = MinimizerWithHash(last_hash, m.contig_id,
                                                                                 last_mx_hit.read_pos,
                                                                                 last_mx_hit.read_strand)
            read_mapping_positions.append(first_mx_hit.read_pos)
            read_mapping_positions.append(last_mx_hit.read_pos)
            # Load the minimizer positions for the contig
            NtLink.list_mx_info[first_hash] = Minimizer(m.contig_id, first_mx_hit.ctg_pos, first_mx_hit.ctg_strand)
            NtLink.list_mx_info[last_hash] = Minimizer(m.contig_id, last_mx_hit.ctg_pos, last_mx_hit.ctg_strand)
        length_read = max(read_mapping_positions)
        self.tally_pairs_from_mappings(accepted_anchor_contigs, contig_runs, length_read, pairs)

    def write_pairs(self, pairs):
        "Write the scaffold pairs to file"
        pair_out = open(self.args.p + ".pairs.tsv", 'w')
        for pair in pairs:
            pair_out.write("\t".join((pair.format_source(), pair.format_target(),
                                      str(pairs[pair]))) + "\n")
        pair_out.close()

    @staticmethod
    def filter_graph_global(graph, min_weight):
        "Filter the graph globally based on minimum edge weight"
        print(datetime.datetime.today(), ": Filtering the graph", file=sys.stdout)
        to_remove_edges = [edge.index for edge in graph.es()
                           if edge['n'] < min_weight]
        new_graph = graph.copy()
        new_graph.delete_edges(to_remove_edges)
        return new_graph

    @staticmethod
    def parse_arguments():
        "Parse ntLink arguments"
        parser = argparse.ArgumentParser(description="ntLink: Scaffolding genome assemblies using long reads")
        parser.add_argument("FILES", nargs="+", help="Long read minimizer TSV files")
        parser.add_argument("-s", help="Target scaffolds fasta file", required=True)
        parser.add_argument("-m", help="Target scaffolds minimizer TSV file", required=True)
        parser.add_argument("-p", help="Output prefix [out]", default="out",
                            type=str, required=False)
        parser.add_argument("-n", help="Minimum edge weight [1]", default=1, type=int)
        parser.add_argument("-k", help="Kmer size used for minimizer step", required=True, type=int)
        parser.add_argument("-z", help="Minimum size of contig to scaffold", required=False, default=500, type=int)
        parser.add_argument("-a", help="Minimum number of anchoring long reads for an edge", required=False,
                            type=int, default=1)
        parser.add_argument("-f", help="Maximum number of contigs in a run for full transitive edge addition",
                            required=False, default=10, type=int)
        parser.add_argument("-x", help="Fudge factor allowed between mapping block lengths on read and assembly. "
                                       "Set to 0 to allow mapping block to be up to read length",
                            type=float, default=0)
        parser.add_argument("-c", "--checkpoint", help="Mappings checkpoint file", required=False)
        parser.add_argument("--pairs", help="Output pairs TSV file", action="store_true")
        parser.add_argument("--paf", help="Output mappings in PAF-like format", action="store_true")
        parser.add_argument("--sensitive", help="Run more sensitive read mapping", action="store_true")
        parser.add_argument("--repeat-filter", help="Remove repetitive minimizers within a long read's sketch",
                            action="store_true")
        parser.add_argument("-v", "--version", action='version', version='ntLink v1.3.5')
        parser.add_argument("--verbose", help="Verbose output logging", action='store_true')

        return parser.parse_args()

    def print_parameters(self):
        "Print the set parameters for the ntLink run"
        print("Parameters:")
        print("\tRead minimizer files: ", self.args.FILES)
        print("\t-s ", self.args.s)
        print("\t-m ", self.args.m)
        print("\t-p ", self.args.p)
        print("\t-n ", self.args.n)
        print("\t-k ", self.args.k)
        print("\t-a ", self.args.a)
        print("\t-z ", self.args.z)
        print("\t-f ", self.args.f)
        print("\t-x ", self.args.x)
        if self.args.checkpoint:
            print("\t-c ", self.args.checkpoint)
        if self.args.sensitive:
            print("\t--sensitive")
        if self.args.repeat_filter:
            print("\t--repeat-filter")

    def main(self):
        "Run ntLink graph stage"
        print("Running pairing stage of ntLink ...\n")
    
        try:
            # Check if the checkpoint mapping file exists
            if os.path.isfile(self.args.p + ".verbose_mapping.tsv"):
                self.args.checkpoint = self.args.p + ".verbose_mapping.tsv"
    
            self.print_parameters()
    
            if self.args.checkpoint:
                print("Found checkpoint file, bypassing read mapping...\n")
                print("Warning: --paf specified, but not compatible with checkpoint")
                NtLink.list_mx_info = {}
            else:
                # Read in the minimizers for target assembly
                mxs_info = self.read_minimizers()
                NtLink.list_mx_info = mxs_info
    
            # Load target scaffolds into memory
            scaffolds = ntlink_utils.read_fasta_file(self.args.s)  # scaffold_id -> Scaffold
            NtLink.scaffolds = scaffolds
    
            if self.args.checkpoint:
                pairs = self.find_scaffold_pairs_checkpoints()
            else:
                # Get directed scaffold pairs, gap estimates from long reads
                pairs = self.find_scaffold_pairs()
    
            pairs = self.filter_pairs_distances(pairs)
    
            pairs = self.filter_weak_anchor_pairs(pairs)
    
            if self.args.pairs:
                self.write_pairs(pairs)
    
            # Build directed graph
            graph = self.build_scaffold_graph(pairs)
    
            # Filter graph
            graph = self.filter_graph_global(graph, int(self.args.n))
    
            # Print out the directed graph
            self.print_directed_graph(graph, "{0}.n{1}".format(self.args.p, self.args.n), scaffolds)
    
            print(datetime.datetime.today(), ": DONE!", file=sys.stdout)
        except:
            if not self.args.checkpoint and self.args.verbose:
               subprocess.call(shlex.split(f"rm {self.args.p}.verbose_mapping.tsv")) 
            if not self.args.checkpoint and self.args.paf:
               subprocess.call(shlex.split(f"rm {self.args.p}.paf")) 
            raise NtlinkPairError("ntLink pairing stage encountered an error..") 

    def __init__(self):
        "Create an ntLink instance"
        self.args = self.parse_arguments()

def main():
    "Run ntLink"
    NtLink().main()

if __name__ == "__main__":
    main()
