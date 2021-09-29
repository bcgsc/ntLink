"""Tests for ntLink"""

import shlex
import subprocess
import re

def cleanup_files(target, prefix, k=32, w=100, n=2):
    "Remove all files in the input list"
    file_list = [f"{target}.k{k}.w{w}.z500.stitch.abyss-scaffold.fa", f"{target}.k{k}.w{w}.tsv",
                 f"{prefix}.pairs.tsv",
                 f"{prefix}.n{n}.scaffold.dot", f"{prefix}.stitch.path",
                 f"{prefix}.trimmed_scafs.fa", f"{prefix}.trimmed_scafs.path",
                 f"{target}.k{k}.w{w}.z500.ntLink.scaffolds.fa"]

    for out_file in file_list:
        command = "rm {0}".format(out_file)
        command_shlex = shlex.split(command)
        return_code = subprocess.call(command_shlex)
        assert return_code == 0

def run_ntLink(target, reads, prefix, k=32, w=100, n=2):
    "Run ntLink, return paths"
    command = "../ntLink scaffold -B target={target} reads={reads} prefix={prefix} k={k} w={w} z=500 n={n}".format(target=target,
                                                                                                                   reads=reads,
                                                                                                                   prefix=prefix,
                                                                                                                   k=k, w=w, n=n)
    command_shlex = shlex.split(command)
    return_code = subprocess.call(command_shlex)
    assert return_code == 0

    test_paths = []
    with open("{prefix}.trimmed_scafs.path".format(prefix=prefix), 'r') as test1_path:
        for line in test1_path:
            line = line.strip().split("\t")
            test_paths.append(line[1])

    return test_paths

def test_1():
    "Testing two sequences together, long reads in fasta format"
    test_paths = run_ntLink("scaffolds_1.fa", "long_reads_1.fa", "test1", w=250)

    expected_paths = ["188266+ 4529N 189231-"]
    for path in test_paths:
        assert path in expected_paths

    # Clean-up files
    cleanup_files("scaffolds_1.fa", "test1", w=250)


def test_2():
    "Testing 4 sequences together, long reads in gzipped fastq format"
    test_paths = run_ntLink("scaffolds_2.fa", "long_reads_2.fq.gz", "test2")

    expected_paths = ["189459+ 73N 183836- 448N 182169- 1311N 190964+"]
    for path in test_paths:
        assert path in expected_paths

    # Clean-up files
    cleanup_files("scaffolds_2.fa", "test2")

def test_3():
    "Testing multiple output paths, long reads in gzipped fasta format"
    test_paths = run_ntLink("scaffolds_3.fa", "long_reads_3.fa.gz", "test3", k=24, w=250)

    expected_paths = ["189459+ 77N 183836- 434N 182169- 1294N 190964+",
                      "188266+ 4566N 189231-"]
    for path in test_paths:
        assert path in expected_paths

    # Clean-up files
    cleanup_files("scaffolds_3.fa", "test3", k=24, w=250)
    