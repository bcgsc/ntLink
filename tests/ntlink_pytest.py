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

def check_trimmed_scaffolds(prefix):
    "Check that the trimmed scaffolds are identical to expected"
    cmd = "cmp {0}.trimmed_scafs.fa expected_outputs/{0}.trimmed_scafs.fa".format(prefix)
    cmd_shlex = shlex.split(cmd)

    return_code = subprocess.call(cmd_shlex)
    assert return_code == 0

def cleanup_files(target, prefix, k=32, w=100, n=2, **kwargs):
    "Remove all files in the input list"
    file_list = [f"{target}.k{k}.w{w}.z1000.stitch.abyss-scaffold.fa", f"{target}.k{k}.w{w}.tsv",
                 f"{prefix}.pairs.tsv",
                 f"{prefix}.n{n}.scaffold.dot", f"{prefix}.stitch.path",
                 f"{target}.k{k}.w{w}.z1000.ntLink.scaffolds.fa", f"{target}.k{k}.w{w}.z1000.ntLink.scaffolds.fa.abyssfac.tsv"]
    if "overlap" not in kwargs or kwargs["overlap"] is not False:
        file_list.extend([f"{prefix}.trimmed_scafs.fa", f"{prefix}.trimmed_scafs.path"])

    for out_file in file_list:
        command = "rm {0}".format(out_file)
        command_shlex = shlex.split(command)
        return_code = subprocess.call(command_shlex)


def run_ntLink(target, reads, prefix, k=32, w=100, n=1, gap_fill=False, **kwargs):
    "Run ntLink, return paths"
    args_str = " ".join(f"{k}={v}" for k, v in kwargs.items())
    if gap_fill:
        args_str += " gap_fill"
    command = "../ntLink scaffold -B target={target} reads={reads} prefix={prefix} k={k} w={w} z=1000 n={n} {extra_args}".format(target=target,
                                                                                                                   reads=reads,
                                                                                                                   prefix=prefix,
                                                                                                                   k=k, w=w, n=n, extra_args=args_str)
    command_shlex = shlex.split(command)
    return_code = subprocess.call(command_shlex)
    assert return_code == 0

    test_paths = []
    path_file = None
    if "overlap" in kwargs and kwargs["overlap"] is False:
        path_file = f"{prefix}.stitch.path"
    else:
        path_file = f"{prefix}.trimmed_scafs.path"
    with open(path_file.format(prefix=prefix), 'r') as test1_path:
        for line in test1_path:
            line = line.strip().split("\t")
            test_paths.append(line[1])

    return test_paths

def test_1():
    "Testing two sequences together, long reads in fasta format"
    test_paths = run_ntLink("scaffolds_1.fa", "long_reads_1.fa", "test1", w=250)

    expected_paths = ["188266+ 4542N 189231-"]
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
    test_paths = run_ntLink("scaffolds_2.fa", "long_reads_2.fq.gz", "test2", k=32, w=100, overlap=False)

    expected_paths = ["189459+ 90N 183836- 449N 182169- 1294N 190964+", '190964- 1294N 182169+ 449N 183836+ 90N 189459-']
    for path in test_paths:
        assert path in expected_paths

    # Test abyss-fac
    scaffolds = "{target}.k{k}.w{w}.z1000.ntLink.scaffolds.fa".format(target="scaffolds_2.fa", k=32, w=100)
    run_abyssfac(scaffolds)
    check_stats(scaffolds + ".abyssfac.tsv")

    # Clean-up files
    cleanup_files("scaffolds_2.fa", "test2", overlap=False)

def test_3():
    "Testing multiple output paths, long reads in gzipped fasta format"
    test_paths = run_ntLink("scaffolds_3.fa", "long_reads_3.fa.gz", "test3", k=24, w=250)

    expected_paths = ["189459+ 71N 183836- 433N 182169- 1315N 190964+",
                      "188266+ 4579N 189231-"]
    for path in test_paths:
        assert path in expected_paths

    # Test abyss-fac
    scaffolds = "{target}.k{k}.w{w}.z1000.ntLink.scaffolds.fa".format(target="scaffolds_3.fa", k=24, w=250)
    run_abyssfac(scaffolds)
    check_stats(scaffolds + ".abyssfac.tsv")

    # Clean-up files
    cleanup_files("scaffolds_3.fa", "test3", k=24, w=250)

def test_4():
    "Testing multiple output paths, long reads in gzipped fasta format, sequences overlap"
    test_paths = run_ntLink("scaffolds_4.fa", "long_reads_4.fa.gz", "scaffolds_4.fa.k40.w100.z1000", k=40, w=100, merge_gap=20)

    expected_paths = ["scaf1+ 21N scaf2+", "scaf3- 21N scaf4+"]

    for path in test_paths:
        assert path in expected_paths

    # Test abyss-fac
    scaffolds = "{target}.k{k}.w{w}.z1000.ntLink.scaffolds.fa".format(target="scaffolds_4.fa", k=40, w=100)
    run_abyssfac(scaffolds)
    check_stats(scaffolds + ".abyssfac.tsv")

    # Compare trimmed scaffolds to expected
    check_trimmed_scaffolds("scaffolds_4.fa.k40.w100.z1000")

    # Clean-up files
    cleanup_files("scaffolds_4.fa", "test4", k=40, w=100)

def test_5():
    "Testing gap-filling target"
    test_paths = run_ntLink("scaffolds_1.fa", "long_reads_1.fa", "test1", gap_fill=True, w=250, gap_k=35)

    # Compare with expected output
    cmd = "cmp {} {}".format("scaffolds_1.fa.k32.w250.z1000.ntLink.scaffolds.gap_fill.fa",
                             "expected_outputs/scaffolds_1.fa.k32.w250.z1000.ntLink.scaffolds.gap_fill.fa")
    cmd_shlex = shlex.split(cmd)
    return_code = subprocess.call(cmd_shlex)
    assert return_code == 0

    # Clean-up files
    cleanup_files("scaffolds_1.fa", "test1", w=250)
