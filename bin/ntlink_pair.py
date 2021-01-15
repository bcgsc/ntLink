#!/usr/bin/env python3
"""
Finding contig pairs using lightweight long read mapping with minimizers
"""
__author__ = 'laurencoombe'

import argparse
import datetime
from collections import defaultdict
from collections import namedtuple
import itertools
import sys
import numpy as np
import igraph as ig
import ntlink_utils

MinimizerEdge = namedtuple("MinimizerEdge", ["mx_i", "mx_i_pos", "mx_i_strand",
                                             "mx_j", "mx_j_pos", "mx_j_strand"])
Minimizer = namedtuple("Minimizer", ["contig", "position", "strand"])
MinimizerWithHash = namedtuple("Minimizer_with_hash", ["mx_hash", "contig", "position", "strand"])

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
    def __init__(self, contig, index, hit_count):
        self.contig = contig
        self.index = index
        self.hit_count = hit_count
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
               "index={index}".format(contig=self.contig, hit_count=self.hit_count,
                                      subsumed=self.subsumed, index=self.index)

class NtLink():
    "Represents an ntLink graph construction run"

    @staticmethod
    def print_directed_graph(graph, out_prefix):
        "Prints the directed scaffold graph in dot format"
        out_graph = out_prefix + ".scaffold.dot"
        outfile = open(out_graph, 'w')
        print(datetime.datetime.today(), ": Printing graph", out_graph, sep=" ", file=sys.stdout)

        outfile.write("digraph G {\n")

        for node in graph.vs():
            node_label = "\"{scaffold}\" [l={length}]\n".\
                format(scaffold=node['name'], length=NtLink.scaffolds[node['name'][:-1]].length)
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

    @staticmethod
    def read_minimizers(tsv_filename):
        "Read the minimizers from a file, removing duplicate minimizers"
        print(datetime.datetime.today(), ": Reading minimizers", tsv_filename, file=sys.stdout)
        mx_info = {}  # mx -> Minimizer object
        dup_mxs = set()  # Set of minimizers identified as duplicates
        with open(tsv_filename, 'r') as tsv:
            for line in tsv:
                line = line.strip().split("\t")
                if len(line) > 1:
                    ctg_name = line[0]
                    mx_pos_split = line[1].split(" ")
                    for mx_pos in mx_pos_split:
                        mx, pos, strand = mx_pos.split(":")
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

    def get_accepted_anchor_contigs(self, mx_list, read_length):
        "Returns dictionary of contigs of appropriate length, mx hits, whether subsumed"
        MinimizerPositions = namedtuple("MinimizerPositions", ["ctg_pos", "read_pos"])
        contig_list = []
        contig_positions = {} # contig -> [mx positions]
        for mx, pos, _ in mx_list:
            contig = NtLink.list_mx_info[mx].contig
            if NtLink.scaffolds[contig].length >= self.args.z:
                contig_list.append(contig)
                if contig not in contig_positions:
                    contig_positions[contig] = []
                contig_positions[contig].append(MinimizerPositions(ctg_pos=NtLink.list_mx_info[mx].position,
                                                                   read_pos=int(pos)))

        # Filter out hits where mapped length on contig is greater than the read length
        noisy_contigs = set()
        for contig in contig_positions:
            positions = contig_positions[contig]
            if len(positions) < 2:
                continue
            ctg_positions = [position.ctg_pos for position in positions]
            start_idx, end_idx = np.argmin(ctg_positions), np.argmax(ctg_positions)
            start_positions = positions[start_idx]
            end_positions = positions[end_idx]
            if self.args.x == 0:
                if abs(end_positions.ctg_pos - start_positions.ctg_pos) > read_length + self.args.k:
                    noisy_contigs.add(contig)
            else:
                threshold = min(read_length + self.args.k,
                                (self.args.x * abs(end_positions.read_pos - start_positions.read_pos)) + self.args.k)
                if abs(end_positions.ctg_pos - start_positions.ctg_pos) > threshold:
                    noisy_contigs.add(contig)
        contig_list = [contig for contig in contig_list if contig not in noisy_contigs]

        contig_runs = [(ctg, len(list(hits))) for ctg, hits in itertools.groupby(contig_list)]
        contigs_hits = {}
        for i, run_tup in enumerate(contig_runs):
            ctg, cnt = run_tup
            if ctg in contigs_hits:
                for j in range(contigs_hits[ctg].index + 1, i):
                    contigs_hits[contig_runs[j][0]].subsumed = True
                contigs_hits[ctg].hit_count += cnt
            else:
                contigs_hits[ctg] = ContigRun(ctg, i, cnt)

        return_contigs_hits = {ctg: contigs_hits[ctg] for ctg in contigs_hits if not contigs_hits[ctg].subsumed}

        return_contig_runs_tmp = [ctg for ctg, hits in contig_runs if not contigs_hits[ctg].subsumed]
        return_contig_runs = [ctg for ctg, hits in itertools.groupby(return_contig_runs_tmp)]

        return return_contigs_hits, return_contig_runs

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

        # Add the long read edges to the graph
        for mx_long_file in self.args.FILES:
            mx_long_filename = "/dev/stdin" if mx_long_file == "-" else mx_long_file
            with open(mx_long_filename, 'r') as long_mxs:
                for line in long_mxs:
                    line = line.strip().split("\t")
                    if len(line) > 1:
                        mx_pos_split_tups = line[1].split(" ")
                        mx_pos_split = []
                        for mx_pos in mx_pos_split_tups:
                            mx, pos, strand = mx_pos.split(":")
                            if mx in target_mxs:
                                mx_pos_split.append((mx, pos, strand))
                        if not mx_pos_split:
                            continue
                        length_long_read = int(mx_pos_split[-1][1])
                        accepted_anchor_contigs, contig_runs = self.get_accepted_anchor_contigs(mx_pos_split,
                                                                                                length_long_read)
                        if self.args.verbose and accepted_anchor_contigs and len(accepted_anchor_contigs) > 1:
                            print(line[0], [str(accepted_anchor_contigs[ctg_run])
                                            for ctg_run in accepted_anchor_contigs], file=sys.stdout)

                        # Filter ordered minimizer list for accepted contigs, keep track of hashes for gap sizes
                        mx_pos_split = [mx_tup for mx_tup in mx_pos_split
                                        if NtLink.list_mx_info[mx_tup[0]].contig in
                                        accepted_anchor_contigs]
                        for mx, pos, strand in mx_pos_split:
                            mx_contig = NtLink.list_mx_info[mx].contig
                            if accepted_anchor_contigs[mx_contig].first_mx is None:
                                accepted_anchor_contigs[mx_contig].first_mx = MinimizerWithHash(mx, mx_contig,
                                                                                                int(pos), strand)
                            accepted_anchor_contigs[mx_contig].terminal_mx = MinimizerWithHash(mx, mx_contig,
                                                                                               int(pos), strand)

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

        return pairs

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
        parser.add_argument("FILES", nargs="+", help="Minimizer TSV files of long reads")
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
        parser.add_argument("-v", "--version", action='version', version='ntLink v0.0.1')
        parser.add_argument("--verbose", help="Verbose output logging", action='store_true')

        return parser.parse_args()

    def print_parameters(self):
        "Print the set parameters for the ntLink run"
        print("Parameters:")
        print("\tReads TSV files: ", self.args.FILES)
        print("\t-s ", self.args.s)
        print("\t-m ", self.args.m)
        print("\t-p ", self.args.p)
        print("\t-n ", self.args.n)
        print("\t-k ", self.args.k)
        print("\t-a ", self.args.a)
        print("\t-z ", self.args.z)
        print("\t-f ", self.args.f)
        print("\t-x ", self.args.x)

    def main(self):
        "Run ntLink graph stage"
        print("Running pairing stage of ntLink ...\n")
        self.print_parameters()

        # Read in the minimizers for target assembly
        mxs_info = self.read_minimizers(self.args.m)
        NtLink.list_mx_info = mxs_info

        # Load target scaffolds into memory
        scaffolds = ntlink_utils.read_fasta_file(self.args.s)  # scaffold_id -> Scaffold
        NtLink.scaffolds = scaffolds

        # Get directed scaffold pairs, gap estimates from long reads
        pairs = self.find_scaffold_pairs()

        pairs = self.filter_pairs_distances(pairs)

        pairs = self.filter_weak_anchor_pairs(pairs)

        self.write_pairs(pairs)

        # Build directed graph
        graph = self.build_scaffold_graph(pairs)

        # Filter graph
        graph = self.filter_graph_global(graph, int(self.args.n))

        # Print out the directed graph
        self.print_directed_graph(graph, "{0}.n{1}".format(self.args.p, self.args.n))

        print(datetime.datetime.today(), ": DONE!", file=sys.stdout)

    def __init__(self):
        "Create an ntLink instance"
        self.args = self.parse_arguments()

def main():
    "Run ntLink"
    NtLink().main()

if __name__ == "__main__":
    main()
