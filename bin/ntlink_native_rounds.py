#!/usr/bin/env python
"""
Adjust scaffold graph for native iterations of ntLink
"""
import argparse
import btllib
import igraph as ig
import ntlink_utils
from collections import namedtuple, defaultdict
import re
import ntlink_pair

PathNodeInfo = namedtuple('PathNodeInfo', ['path_id', "orientation", "five_prime_terminal", "three_prime_terminal"])

def read_scaffolds(fasta: str, threads: argparse.Namespace) -> dict:
    "Read the scaffolds from the fasta file, storing as a dictonary with the scaffold name \
    as the key and an instance of the Scaffold namedtuple as the value"
    scaffolds = {}
    with btllib.SeqReader(fasta, btllib.SeqReaderFlag.LONG_MODE, threads) as seqs:
        for seq in seqs:
            scaffolds[seq.id] = ntlink_utils.Scaffold(seq.id, len(seq.seq))
    return scaffolds

def read_path_file(path_file: str) -> tuple:
    "Read the path file and return a dictionary with the constituent scaffold as the key, and PathNodeInfo as the value"
    path_graphs = defaultdict(dict)
    scaf_paths = {}
    with open(path_file, 'r') as path_file:
        for line in path_file:
            path_id, nodes = line.strip().split('\t')
            nodes = nodes.split(" ")
            nodes = [node for node in nodes if not node[-1] == "N"]
            rev_nodes = list(reversed([ntlink_utils.reverse_scaf_ori(node) for node in nodes]))

            vertices = nodes + rev_nodes
            fwd_edges = []
            for i, j in zip(nodes, nodes[1:]):
                fwd_edges.append((i, j))
            rev_edges = []
            for i, j in zip(rev_nodes, rev_nodes[1:]):
                rev_edges.append((i, j))
            path_graphs[path_id]["+"] = ig.Graph(directed=True)
            path_graphs[path_id]["+"].add_vertices(nodes)
            path_graphs[path_id]["+"].add_edges(fwd_edges)

            path_graphs[path_id]["-"] = ig.Graph(directed=True)
            path_graphs[path_id]["-"].add_vertices(rev_nodes)
            path_graphs[path_id]["-"].add_edges(rev_edges)

            for node in nodes:
                scaf_paths[node.strip("+-")] = path_id

    return path_graphs, scaf_paths

def adjust_node(node: str, path_graphs: dict, scaf_path: dict, source: bool = False) -> str:
    "Return True if the node is terminal compatible with the paths, False otherwise"
    node_name, node_ori = node[:-1], node[-1]
    if node_name not in scaf_path:
        return node # Node is not in the path file, so it is compatible
    path_id = scaf_path[node_name]
    for ori in path_graphs[path_id]:
        graph = path_graphs[path_id][ori]
        if source and ntlink_utils.has_vertex(graph, node) and graph.outdegree(node) == 0:
            return path_id + ori
        if not source and ntlink_utils.has_vertex(graph, node) and graph.indegree(node) == 0:
            return path_id + ori
    return "None"



def adjust_edge(source: str, target: str, path_graphs: dict, scaf_path: dict) -> tuple:
    "Adjust the edge for the paths"
    source_adjust = adjust_node(source, path_graphs, scaf_path, source=True)
    target_adjust = adjust_node(target, path_graphs, scaf_path, source=False)
    return (source_adjust, target_adjust)


def adjust_scaffold_graph(graph_file: str, scaffold_lengths: dict, path_graphs: dict, scaf_paths: dict) -> ig.Graph:
    "Read the scaffold graph and adjust it for the paths"
    g = ig.Graph(directed=True)

    # Add the scaffolds to the graph in both orientations
    g.add_vertices(["{}+".format(s) for s in scaffold_lengths.keys()])
    g.add_vertices(["{}-".format(s) for s in scaffold_lengths.keys()])

    # Step through the edges in the graph, adjusting as needed
    edge_re = re.compile(r'\"(\S+[+-])\"\s+\-\>\s+\"(\S+[+-])\"\s+\[d\=(\-?\d+)\s+e\=(\d+)\s+n\=(\d+)\]')
    edges = defaultdict(dict)
    with open(graph_file, 'r') as graph_file:
        for line in graph_file:
            edge_match = edge_re.match(line)
            if not edge_match:
                continue
            source, target, distance, edge_id, edge_weight = edge_match.groups()
            source, target = adjust_edge(source, target, path_graphs, scaf_paths)
            if source == "None" or target == "None":
                continue
            edges[source][target] = (distance, edge_id, edge_weight)
    formatted_edges = [(s, t) for s in edges for t in edges[s]]
    g.add_edges(formatted_edges)
    edge_attributes = {ntlink_utils.edge_index(g, s, t): {'d': edges[s][t][0],
                                                          "e": edges[s][t][1],
                                                          "n": edges[s][t][2]}
                       for s in edges for t in edges[s]}
    g.es()["d"] = [edge_attributes[e]['d'] for e in sorted(edge_attributes.keys())]
    g.es()["e"] = [edge_attributes[e]['e'] for e in sorted(edge_attributes.keys())]
    g.es()["n"] = [edge_attributes[e]['n'] for e in sorted(edge_attributes.keys())]
    return g


def main():
    parser = argparse.ArgumentParser(description='Adjust scaffold graph for native iterations of ntLink')
    parser.add_argument('-i', '--input', help='Input scaffold graph', required=True)
    parser.add_argument('-f', '--fasta', help='Input fasta file', required=True)
    parser.add_argument("--path", help="Path file from ntLink", required=True)
    parser.add_argument("-t", "--threads", help="Number of threads to use", type=int, default=1)
    parser.add_argument("--out", help="Output scaffold graph prefix", required=True)
    args = parser.parse_args()

    # Read in scaffold lengths
    scaffold_lengths = read_scaffolds(args.fasta, args.threads)

    # Read in the path file
    path_graphs, scaf_paths = read_path_file(args.path)

    # Read through the scaffold graph, adjusting the graph for the paths
    g = adjust_scaffold_graph(args.input, scaffold_lengths, path_graphs, scaf_paths)

    ntlink_pair.NtLink.print_directed_graph(g, args.out, scaffold_lengths)


if __name__ == "__main__":
    main()