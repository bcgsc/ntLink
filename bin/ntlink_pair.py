#!/usr/bin/env python3
"""
Using ntJoin concept with long reads as input reference
"""
__author__ = 'laurencoombe'

import argparse
import datetime
import re
from collections import defaultdict
from collections import namedtuple
import itertools
import sys
import numpy as np
import igraph as ig
from ntjoin_assemble import Ntjoin


MinimizerEdge = namedtuple("MinimizerEdge", ["mx_i", "mx_i_pos", "mx_i_strand",
                                             "mx_j", "mx_j_pos", "mx_j_strand",
                                             "ONT_name"])
Minimizer = namedtuple("Minimizer", ["contig", "position", "strand"])
Minimizer_with_hash = namedtuple("Minimizer_with_hash", ["mx_hash", "contig", "position", "strand"])
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

    def add_gap_estimate(self, gap_est):
        "Add a gap estimate"
        self._gap_est.append(gap_est)

    def get_gap_estimate(self):
        "Calculate gap estimate"
        return int(np.median(self._gap_est))

    def n_supporting_reads(self):
        "Return the number of supporting reads for a pair"
        return len(self._gap_est)

    def __str__(self):
        return "n={n}, gap_estimates={gap_ests}, anchor={a}".format(n=self.n_supporting_reads(),
                                                                    gap_ests=self._gap_est,
                                                                    a=self.anchor)

class ContigRun:
    "Represents information about a contig run based on a long reaad"
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

class NtjoinLong():
    "Represents an Ntjoin graph construction run"

    @staticmethod
    def print_directed_graph(graph, out_prefix):
        "Prints the directed scaffold graph in dot format"
        out_graph = out_prefix + ".scaffold.dot"
        outfile = open(out_graph, 'w')
        print(datetime.datetime.today(), ": Printing graph", out_graph, sep=" ", file=sys.stdout)

        outfile.write("digraph G {\n")

        for node in graph.vs():
            node_label = "\"{scaffold}\" [l={length}]\n".\
                format(scaffold=node['name'], length=NtjoinLong.scaffolds[node['name'][:-1]].length)
            outfile.write(node_label)

        for edge in graph.es():
            edge_str = "\"{source}\" -> \"{target}\" [d={d} e={e} n={n}]\n".\
                format(source=Ntjoin.vertex_name(graph, edge.source),
                       target=Ntjoin.vertex_name(graph, edge.target),
                       d=int(edge['d']), e=edge['e'], n=edge['n'])
            outfile.write(edge_str)

        outfile.write("}\n")

    @staticmethod
    def convert_links_scaf_id(vertex_name, links_numbering):
        "Convert abyss-scaffold style naming to LINKS-style"
        if vertex_name[-1] == "+":
            return "{ori}{vertex}".format(ori="f", vertex=links_numbering[vertex_name[:-1]])
        return "{ori}{vertex}".format(ori="r", vertex=links_numbering[vertex_name[:-1]])

    def print_links_tigpair(self, graph, prefix):
        "Prints pairs in LINKS format"
        min_match = re.search(r'^(\S+).k\d+.w\d+\.tsv', self.args.s)
        assembly_fa = min_match.group(1)
        links_numbering = {}
        header_re = re.compile(r'>(\S+)')
        counter = 1
        with open(assembly_fa, 'r') as assembly:
            for line in assembly:
                line = line.strip()
                header_match = re.search(header_re, line)
                if header_match:
                    links_numbering[header_match.group(1)] = counter
                    counter += 1

        out_file_name = prefix + ".tigpair_checkpoint.tsv"
        out_file = open(out_file_name, 'w')
        for edge in graph.es():
            source = self.convert_links_scaf_id(Ntjoin.vertex_name(graph, edge.source), links_numbering)
            target = self.convert_links_scaf_id(Ntjoin.vertex_name(graph, edge.target), links_numbering)
            distance_est = int(edge['d'])
            links = int(edge['n'])

            if distance_est < 0:
                dist_category = -1
            elif distance_est <= 10:
                dist_category = 10
            elif distance_est <= 500:
                dist_category = 500
            elif distance_est <= 1000:
                dist_category = 1000
            elif distance_est <= 5000:
                dist_category = 5000
            else:
                dist_category = 10000
            out_str = "{dist_category}\t{source}\t{target}\t{links}\t{gap_links}\n".format(dist_category=dist_category,
                                                                                           source=source, target=target,
                                                                                           links=links,
                                                                                           gap_links=links*distance_est)
            out_file.write(out_str)


    def calculate_gap_size_check(self, u_mx, u_ori, v_mx, v_ori, est_distance):
        "Calculates the estimated distance between two contigs, assuming full contig"

        # Correct for the overhanging sequence before/after terminal minimizers
        if u_ori == "+":
            u_ctg = NtjoinLong.list_mx_info[self.args.s][u_mx].contig
            u_ctglen = NtjoinLong.scaffolds[u_ctg].length
            a = u_ctglen - NtjoinLong.list_mx_info[self.args.s][u_mx].position - self.args.k
        else:
            a = NtjoinLong.list_mx_info[self.args.s][u_mx].position
        if v_ori == "+":
            b = NtjoinLong.list_mx_info[self.args.s][v_mx].position
        else:
            v_ctg = NtjoinLong.list_mx_info[self.args.s][v_mx].contig
            v_ctglen = NtjoinLong.scaffolds[v_ctg].length
            b = v_ctglen - NtjoinLong.list_mx_info[self.args.s][v_mx].position - self.args.k

        try:
            assert a >= 0
            assert b >= 0
        except:
            print("ERROR: Gap distance estimation less than 0", "Vertex 1:", u_mx, "Vertex 2:", v_mx,
                  sep="\n")
            print("Minimizer positions:", NtjoinLong.list_mx_info[self.args.s][u_mx].position,
                  NtjoinLong.list_mx_info[self.args.s][v_mx].position)
            print("Estimated distance: ", est_distance)
            raise AssertionError

        gap_size = est_distance - a - b
        return int(gap_size)

    @staticmethod
    def read_minimizers(tsv_filename):
        "Read the minimizers from a file, removing duplicate minimizers"
        print(datetime.datetime.today(), ": Reading minimizers", tsv_filename, file=sys.stdout)
        mx_info = {}  # mx -> (contig, position)
        mxs = []  # List of lists of minimizers
        dup_mxs = set()  # Set of minimizers identified as duplicates
        with open(tsv_filename, 'r') as tsv:
            for line in tsv:
                line = line.strip().split("\t")
                if len(line) > 1:
                    mx_pos_split = line[1].split(" ")
                    mxs.append([mx_pos.split(":")[0] for mx_pos in mx_pos_split])
                    for mx_pos in mx_pos_split:
                        mx, pos, strand = mx_pos.split(":")
                        if mx in mx_info:  # This is a duplicate, add to dup set, don't add to dict
                            dup_mxs.add(mx)
                        else:
                            mx_info[mx] = Minimizer(line[0], int(pos), strand)

        mx_info = {mx: mx_info[mx] for mx in mx_info if mx not in dup_mxs}
        mxs_filt = []
        for mx_list in mxs:
            mx_list_filt = [mx for mx in mx_list if mx not in dup_mxs]
            mxs_filt.append(mx_list_filt)
        return mx_info, mxs_filt

    @staticmethod
    def reverse_orientation(orientation):
        "Flip the given orientation"
        assert orientation in ("+", "-")
        if orientation == "+":
            return "-"
        return "+"

    def normalize_pair(self, source_ctg, source_ori, target_ctg, target_ori):
        "Normalize pairs. Lexicographically smallest contig first in pair"
        if source_ctg < target_ctg:
            return ScaffoldPair(source_ctg, source_ori, target_ctg, target_ori)
        return ScaffoldPair(target_ctg, self.reverse_orientation(target_ori),
                            source_ctg, self.reverse_orientation(source_ori))


    def find_scaffolding_edges(self, graph):
        "Iterates over edges in a graph, and tallies connections between different scaffolds"
        pairs = {} # (source + ori, target + ori)
        for edge in graph.es():
            if "target" in graph.es()[edge.index]["support"]:
                continue
            for mx_edge in graph.es()[edge.index]["support"]:
                assert mx_edge.mx_i_pos < mx_edge.mx_j_pos
                assert mx_edge.mx_i == Ntjoin.vertex_name(graph, edge.source) or \
                       mx_edge.mx_i == Ntjoin.vertex_name(graph, edge.target)
                assert mx_edge.mx_j == Ntjoin.vertex_name(graph, edge.source) or \
                       mx_edge.mx_j == Ntjoin.vertex_name(graph, edge.target)
                if mx_edge.mx_i_strand == NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_i].strand:
                    source_ctg, source_ori = NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_i].contig, "+"
                else:
                    source_ctg, source_ori = NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_i].contig, "-"
                if mx_edge.mx_j_strand == NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_j].strand:
                    target_ctg, target_ori = NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_j].contig, "+"
                else:
                    target_ctg, target_ori = NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_j].contig, "-"
                new_pair = self.normalize_pair(source_ctg, source_ori, target_ctg, target_ori)
                gap_estimate = self.calculate_gap_size_check(mx_edge.mx_i, source_ori, mx_edge.mx_j, target_ori,
                                                             mx_edge.mx_j_pos - mx_edge.mx_i_pos)
                if new_pair not in pairs:
                    pairs[new_pair] = []
                pairs[new_pair].append((gap_estimate, mx_edge.ONT_name))
        return pairs

    def calculate_pair_info(self, mx_edge):
        "Given a contig pair, normalizes, defines orientation and estimates the gap size"

        assert mx_edge.mx_i_pos < mx_edge.mx_j_pos
        if mx_edge.mx_i_strand == NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_i].strand:
            source_ctg, source_ori = NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_i].contig, "+"
        else:
            source_ctg, source_ori = NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_i].contig, "-"
        if mx_edge.mx_j_strand == NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_j].strand:
            target_ctg, target_ori = NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_j].contig, "+"
        else:
            target_ctg, target_ori = NtjoinLong.list_mx_info[self.args.s][mx_edge.mx_j].contig, "-"
        new_pair = self.normalize_pair(source_ctg, source_ori, target_ctg, target_ori)
        gap_estimate = self.calculate_gap_size_check(mx_edge.mx_i, source_ori, mx_edge.mx_j, target_ori,
                                                     mx_edge.mx_j_pos - mx_edge.mx_i_pos)
        return (new_pair, gap_estimate)

    def filter_weak_anchor_pairs(self, pairs):
        "Filter out edges where there isn't at least threshold # well-anchored read"
        new_pair = {pair: pairs[pair] for pair in pairs if pairs[pair].anchor >= self.args.a}
        return new_pair

    @staticmethod
    def filter_pairs_distances(pairs):
        "Filter out edges where distance estimate is incongruous with either scaffold"
        new_pairs = {}
        for pair in pairs:
            if pairs[pair].get_gap_estimate() <= NtjoinLong.scaffolds[pair.source_contig].length*-1 or \
                pairs[pair].get_gap_estimate() <= NtjoinLong.scaffolds[pair.target_contig].length*-1:
                continue
            new_pairs[pair] = pairs[pair]
        return new_pairs


    def reverse_complement_pair(self, pair):
        "Reverse complemented a directed pair of scaffolds"
        return ScaffoldPair(pair.target_contig, self.reverse_orientation(pair.target_ori),
                            pair.source_contig, self.reverse_orientation(pair.source_ori))

    def build_scaffold_graph(self, pairs):
        "Builds a scaffold graph given the pairs info"

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

        edge_attributes = {Ntjoin.edge_index(graph, s, t): {'d': edges[s][t].get_gap_estimate(),
                                                            "e": 100,
                                                            "n": edges[s][t].n_supporting_reads()}
                           for s in edges for t in edges[s]}
        graph.es()["d"] = [edge_attributes[e]['d'] for e in sorted(edge_attributes.keys())]
        graph.es()["e"] = [edge_attributes[e]['e'] for e in sorted(edge_attributes.keys())]
        graph.es()["n"] = [edge_attributes[e]['n'] for e in sorted(edge_attributes.keys())]

        max_weight = max([edge_attributes[e]['n'] for e in edge_attributes.keys()])

        return graph, max_weight

    def get_accepted_anchor_contigs(self, mx_list):
        "Returns dictionary of contigs of appropriate length, mx hits, whether subsumed"
        contig_list = []
        for mx, _, _ in mx_list:
            contig = NtjoinLong.list_mx_info[self.args.s][mx].contig
            if NtjoinLong.scaffolds[contig].length >= self.args.z:
                contig_list.append(contig)

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

    def add_pair(self, accepted_anchor_contigs, ctg_i, ctg_j, ont_read, pairs, check_added=None):
        "Add pair to dictionary of pairs"
        mx_i = accepted_anchor_contigs[ctg_i].terminal_mx
        mx_j = accepted_anchor_contigs[ctg_j].first_mx
        pair, gap_est = self.calculate_pair_info(MinimizerEdge(mx_i.mx_hash, mx_i.position, mx_i.strand,
                                                               mx_j.mx_hash, mx_j.position, mx_j.strand,
                                                               ont_read))
        if check_added is not None and pair in check_added:
            return

        if pair not in pairs:
            pairs[pair] = PairInfo()
        pairs[pair].add_gap_estimate(gap_est)
        if accepted_anchor_contigs[ctg_i].hit_count > 1 and \
                        accepted_anchor_contigs[ctg_j].hit_count > 1:
            pairs[pair].anchor += 1

    def find_scaffold_pairs(self, list_mxs_target):
        "Builds up pairing information between scaffolds"
        print(datetime.datetime.today(), ": Finding pairs", file=sys.stdout)

        target_mxs = {mx for mx_list in list_mxs_target for mx in mx_list}

        pairs = {} # source -> target -> [gap estimate]

        # Add the long read edges to the graph
        file_count = 0
        for mx_long_file in self.args.FILES:
            file_count += 1
            with open(mx_long_file, 'r') as long_mxs:
                for line in long_mxs:
                    line = line.strip().split("\t")
                    if len(line) > 1:
                        mx_pos_split_tups = line[1].split(" ")
                        mx_pos_split = []
                        for mx_pos in mx_pos_split_tups:
                            mx, pos, strand = mx_pos.split(":")
                            if mx in target_mxs:
                                mx_pos_split.append((mx, pos, strand))
                        accepted_anchor_contigs, contig_runs = self.get_accepted_anchor_contigs(mx_pos_split)
                        if self.args.verbose and accepted_anchor_contigs and len(accepted_anchor_contigs) > 1:
                            print(line[0], [str(accepted_anchor_contigs[ctg_run])
                                            for ctg_run in accepted_anchor_contigs])
                        mx_pos_split = [mx_tup for mx_tup in mx_pos_split
                                        if NtjoinLong.list_mx_info[self.args.s][mx_tup[0]].contig in
                                        accepted_anchor_contigs]
                        for mx, pos, strand in mx_pos_split:
                            mx_contig = NtjoinLong.list_mx_info[self.args.s][mx].contig
                            if mx_contig not in accepted_anchor_contigs:
                                continue
                            if accepted_anchor_contigs[mx_contig].first_mx is None:
                                accepted_anchor_contigs[mx_contig].first_mx = Minimizer_with_hash(mx, mx_contig, int(pos), strand)
                            accepted_anchor_contigs[mx_contig].terminal_mx = Minimizer_with_hash(mx, mx_contig, int(pos), strand)

                        if len(contig_runs) <= self.args.m:
                            # Add all transitive edges for pairs
                            for ctg_pair in itertools.combinations(contig_runs, 2):
                                ctg_i, ctg_j = ctg_pair
                                self.add_pair(accepted_anchor_contigs, ctg_i, ctg_j, line[0], pairs)
                        else:
                            added_pairs = set()
                            # Add adjacent pairs
                            for ctg_i, ctg_j in zip(contig_runs, contig_runs[1:]):
                                new_pair = self.add_pair(accepted_anchor_contigs, ctg_i, ctg_j, line[0], pairs)
                                added_pairs.add(new_pair)

                            # Add transitive edges over weakly supported contigs
                            contig_runs_filter = [ctg for ctg in contig_runs if accepted_anchor_contigs[ctg].hit_count > 1]
                            for ctg_i, ctg_j in zip(contig_runs_filter, contig_runs_filter[1:]):
                                self.add_pair(accepted_anchor_contigs, ctg_i, ctg_j, line[0], pairs, check_added=added_pairs)

        return pairs

    def filter_graph_global(self, graph):
        "Filter the graph globally based on minimum edge weight"
        print(datetime.datetime.today(), ": Filtering the graph", file=sys.stdout)
        to_remove_edges = [edge.index for edge in graph.es()
                           if edge['n'] < self.args.n]
        new_graph = graph.copy()
        new_graph.delete_edges(to_remove_edges)
        return new_graph

    @staticmethod
    def parse_arguments():
        "Parse ntJoin arguments"
        parser = argparse.ArgumentParser(
            description="ntJoin: Scaffolding genome assemblies using reference assemblies and minimizer graphs",
            epilog="Note: Script expects that each input minimizer TSV file has a matching fasta file.\n"
                   "Example: myscaffolds.fa.k32.w1000.tsv - myscaffolds.fa is the expected matching fasta",
            formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("FILES", nargs="+", help="Minimizer TSV files of references")
        parser.add_argument("-a", help="Minimum number of anchoring ONT for an edge", required=False,
                            type=int, default=1)
        parser.add_argument("-s", help="Target scaffolds minimizer TSV file", required=True)
        parser.add_argument("-p", help="Output prefix [out]", default="out",
                            type=str, required=False)
        parser.add_argument("-n", help="Minimum edge weight [1]", default=1, type=int)
        parser.add_argument("-k", help="Kmer size used for minimizer step", required=True, type=int)
        parser.add_argument("-m", help="Maximum number of contigs in a run for full transitive edge addition",
                            required=False, default=10, type=int)
        parser.add_argument("-z", help="Minimum size of contig to scaffold", required=False, default=500, type=int)
        parser.add_argument("-v", "--version", action='version', version='ntJoin v1.0.2')
        parser.add_argument("--verbose", help="Verbose output logging", action='store_true')
        parser.add_argument("--agp", help="Output AGP file describing scaffolds", action="store_true")

        return parser.parse_args()

    def print_parameters(self):
        "Print the set parameters for the ntJoin run"
        print("Parameters:")
        print("\tReference TSV files: ", self.args.FILES)
        print("\t-s ", self.args.s)
        print("\t-p ", self.args.p)
        print("\t-n ", self.args.n)
        print("\t-k ", self.args.k)
        print("\t-a ", self.args.a)
        print("\t-z ", self.args.z)
        print("\t-m ", self.args.m)
        if self.args.agp:
            print("\t--agp")

    def main(self):
        "Run ntJoin graph stage"
        print("Running ntJoin long ...\n")
        self.print_parameters()

        # Read in the minimizers for target assembly
        mxs_info, mxs = self.read_minimizers(self.args.s)
        NtjoinLong.list_mx_info = {self.args.s: mxs_info}
        Ntjoin.list_mx_info = {self.args.s: mxs_info}

        # Load target scaffolds into memory
        min_match = re.search(r'^(\S+).k\d+.w\d+\.tsv', self.args.s)
        if not min_match:
            print("ERROR: Target assembly minimizer TSV file must follow the naming convention:")
            print("\ttarget_assembly.fa.k<k>.w<w>.tsv, where <k> and <w> are parameters used for minimizering")
            sys.exit(1)
        assembly_fa = min_match.group(1)
        scaffolds = Ntjoin.read_fasta_file(assembly_fa)  # scaffold_id -> Scaffold
        NtjoinLong.scaffolds = scaffolds
        Ntjoin.scaffolds = scaffolds

        # Get directed scaffold pairs, gap estimates from long reads
        pairs = self.find_scaffold_pairs(mxs)

        pairs = self.filter_pairs_distances(pairs)

        pairs = self.filter_weak_anchor_pairs(pairs)

        pair_out = open(self.args.p + ".pairs.tsv", 'w')
        for pair in pairs:
            pair_out.write("\t".join([pair.source_contig + pair.source_ori,
                                      pair.target_contig + pair.target_ori, str(pairs[pair])]) + "\n")
        pair_out.close()

        # Build directed graph
        graph, max_weight = self.build_scaffold_graph(pairs)
        graph = self.filter_graph_global(graph)

        NtjoinLong.max_weight = max_weight

        # Print out the directed graph
        self.print_directed_graph(graph, self.args.p)

        # Print LINKS-formatted graph
        self.print_links_tigpair(graph, self.args.p)

        print(datetime.datetime.today(), ": DONE!", file=sys.stdout)

    def __init__(self):
        "Create an ntJoin instance"
        self.args = self.parse_arguments()

def main():
    "Run ntJoin"
    NtjoinLong().main()

if __name__ == "__main__":
    main()
