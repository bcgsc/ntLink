"""Tests for ntLink"""

import shlex
import subprocess
import pandas as pd

def run_abyssfac(scaffolds):
    "Run abyss-fac on scaffolds"
    cmd = "abyss-fac {}".format(scaffolds)
    cmd_shlex = shlex.split(cmd)
    with open(scaffolds + ".abyssfac.tsv", 'w') as outfile:
        return_code = subprocess.call(cmd_shlex, stdout=outfile)
    assert return_code == 0

def check_stats(abyssfac_filename):
    "Check stats of longstitch scaffolds"
    reference_stats = pd.read_csv("expected_outputs/{}".format(abyssfac_filename), sep="\t")
    ci_stats = pd.read_csv(abyssfac_filename, sep="\t")

    assert int(reference_stats["N50"]) == int(ci_stats["N50"])
    assert int(reference_stats["n"]) == int(ci_stats["n"])

def cleanup_files(target, prefix, k=32, w=100, n=2):
    "Remove all files in the input list"
    file_list = [f"{target}.k{k}.w{w}.z1000.stitch.abyss-scaffold.fa", f"{target}.k{k}.w{w}.tsv",
                 f"{prefix}.pairs.tsv",
                 f"{prefix}.n{n}.scaffold.dot", f"{prefix}.stitch.path",
                 f"{prefix}.trimmed_scafs.fa", f"{prefix}.trimmed_scafs.path",
                 f"{target}.k{k}.w{w}.z1000.ntLink.scaffolds.fa"]

    for out_file in file_list:
        command = "rm {0}".format(out_file)
        command_shlex = shlex.split(command)
        return_code = subprocess.call(command_shlex)
        assert return_code == 0

def run_ntLink(target, reads, prefix, k=32, w=100, n=2):
    "Run ntLink, return paths"
    command = "../ntLink scaffold -B target={target} reads={reads} prefix={prefix} k={k} w={w} z=1000 n={n}".format(target=target,
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

    # Test abyss-fac
    scaffolds = "{target}.k{k}.w{w}.z1000.ntLink.scaffolds.fa".format(target="scaffolds_1.fa", k=32, w=250)
    run_abyssfac(scaffolds)
    check_stats(scaffolds + ".abyssfac.tsv")

    # Clean-up files
    cleanup_files("scaffolds_1.fa", "test1", w=250)


def test_2():
    "Testing 4 sequences together, long reads in gzipped fastq format"
    test_paths = run_ntLink("scaffolds_2.fa", "long_reads_2.fq.gz", "test2")

    expected_paths = ["189459+ 73N 183836- 448N 182169- 1311N 190964+"]
    for path in test_paths:
        assert path in expected_paths

    # Test abyss-fac
    scaffolds = "{target}.k{k}.w{w}.z1000.ntLink.scaffolds.fa".format(target="scaffolds_2.fa", k=32, w=100)
    run_abyssfac(scaffolds)
    check_stats(scaffolds + ".abyssfac.tsv")

    # Clean-up files
    cleanup_files("scaffolds_2.fa", "test2")

def test_3():
    "Testing multiple output paths, long reads in gzipped fasta format"
    test_paths = run_ntLink("scaffolds_3.fa", "long_reads_3.fa.gz", "test3", k=24, w=250)

    expected_paths = ["189459+ 77N 183836- 434N 182169- 1294N 190964+",
                      "188266+ 4566N 189231-"]
    for path in test_paths:
        assert path in expected_paths

    # Test abyss-fac
    scaffolds = "{target}.k{k}.w{w}.z1000.ntLink.scaffolds.fa".format(target="scaffolds_3.fa", k=24, w=250)
    run_abyssfac(scaffolds)
    check_stats(scaffolds + ".abyssfac.tsv")

    # Clean-up files
    cleanup_files("scaffolds_3.fa", "test3", k=24, w=250)
