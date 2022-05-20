#!/bin/bash

# Use this script to test your ntLink installation.
# This script expects you to have ntJoin in your PATH, and all dependencies installed (See README for details)

set -eux

ntLink scaffold -B target=scaffolds_1.fa reads=long_reads_1.fa  w=250 

ntLink scaffold -B target=scaffolds_2.fa reads=long_reads_2.fq.gz w=100

ntLink scaffold -B target=scaffolds_3.fa reads=long_reads_3.fa.gz k=24 w=250

ntLink scaffold -B target=scaffolds_4.fa reads=long_reads_4.fa.gz k=40 w=100

echo "Done tests! Compare your generated files with the files in the expected_outputs folder to ensure the tests were successful."
