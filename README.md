![GitHub release (latest by date)](https://img.shields.io/github/v/release/bcgsc/ntlink)
![Conda](https://img.shields.io/conda/dn/bioconda/ntlink?label=Conda)
[![Published in BMC Bioinformatics](https://img.shields.io/badge/Published%20in-BMC%20Bioinformatics-blue)](https://doi.org/10.1186/s12859-021-04451-7)

![Logo](https://github.com/bcgsc/ntLink/blob/master/ntlink-logo.png)

## Minimizer-based genome assembly scaffolding and mapping using long reads and minimizers

## Description of the algorithm
ntLink uses minimizers to perform a lightweight mapping between the input target assembly and the supplied long reads. These long-read mappings are then used as evidence to orient and order the output scaffolds.

### General steps in the algorithm:
1. Compute ordered minimizer sketches of the input target assembly and long reads
2. Use minimizers to map the long reads to the target assembly contigs
3. Find contig pairs, where joins are suggested by the long-read mapping evidence
4. Output a scaffold graph, where the nodes are oriented contigs and the edges are joins suggested by the long read data
5. Traverse the scaffold graph using `abyss-scaffold` to output the final scaffolds

## Credits
Concept: Rene Warren and Lauren Coombe

Design and implementation: Lauren Coombe

## Citing ntLink
If you use ntLink in your research, please cite:

Coombe L, Li JX, Lo T, Wong J, Nikolic V, Warren RL and Birol I. LongStitch: high-quality genome assembly correction and scaffolding using long reads. BMC Bioinformatics 22, 534 (2021). https://doi.org/10.1186/s12859-021-04451-7

Coombe L, Warren RL, Wong J, Nikolic V, and Birol I. ntLink: A toolkit for de novo genome assembly scaffolding and mapping using long reads. arXiv (2023). https://doi.org/10.48550/ARXIV.2301.08785


## Usage
```
ntLink: Scaffolding assemblies using long reads
Usage: ntLink scaffold target=<target scaffolds> reads='List of long read files'

To additionally run gap-filling (fill gap regions with raw read sequence):
Usage: ntLink scaffold gap_fill target=<target scaffolds> reads='List of long read files'

Options:
target			Target assembly to be scaffolded in fasta format
reads		        List of long read files (separated by a space)
prefix			Prefix of intermediate output files [<target>.k<k>.w<w>.n<n>]
t			Number of threads [4]
k			K-mer size for minimizers [32]
w			Window size for minimizers [100]
n			Minimum graph edge weight [1]
g			Minimum gap size (bp) [20]
f			Maximum number of contigs in a run for full transitive edge addition [10]
a                       Minimum number of anchored ONT reads required for an edge [1]
z			Minimum size of contig (bp) to scaffold [1000]
v                       If 1, track time and memory for each step of the pipeline [0]
paf			If True, outputs read to contig mappings in PAF-like format [False]
overlap			If True, runs extra step to attempt to identify and trim overlapping joined sequences [True]
sensitive	        If True, runs mapping in sensitive mode [False]
soft_mask		If True, gaps are filled with lowercase bases [False]

Note: 
	- Ensure all assembly and read files are in the current working directory, making soft links if necessary
```

Running `ntLink help` prints the help documentation.

* Input reads files can be gzipped (or not), and in either fastq or fasta format

### Example
Input files:
* target assembly `my_assembly.fa` 
* long read file `long_reads.fq.gz`

ntLink command:
```
ntLink scaffold target=my_assembly.fa reads=long_reads.fq.gz k=32 w=250
```

The post-ntLink scaffolds file will have the suffix `*ntLink.scaffolds.fa`

See our [wiki](https://github.com/bcgsc/ntLink/wiki) for more information about output file formats.


### Gap-filling
As of ntLink v1.2.0, ntLink can also run gap-filling after the scaffolding stage. This mode is enabled by adding the `gap_fill` target to the `ntLink` command. `overlap=True` is required when using the `gap_fill` feature.

Note that the gaps will be filled with raw read sequence, so subsequent polishing is a good idea. See the wiki page for more details.

### Rounds

To maximize the scaffolding gains, ntLink can be run iteratively in rounds. As of ntLink v1.3.0, these rounds can be launched using the `ntLink_rounds` Makefile, which uses mapping liftover to reduce the computational cost of additional ntLink rounds.

Example command without gap-filling (target `run_rounds_gaps` runs gap-filling, while `run_rounds` does not), running 5 rounds of ntLink:
```
ntLink_rounds run_rounds target=my_assembly.fa reads=long_reads.fq.gz k=24 w=250 rounds=5
```
See the wiki page for more details.

### Mapping only

To only run the pairing stage of `ntLink` (the stage where the long reads are mapped to the contigs), use the `pair` target for the `ntLink` command. The mappings can also be output in PAF-like format by specifying `paf=True`.


**For more information about the ntLink algorithm and tips for running ntLink see our [wiki](https://github.com/bcgsc/ntLink/wiki)**

 ## Installation
 ntLink is available from conda and homebrew package managers.
 
 Installing using conda:
 ```
 conda install -c bioconda -c conda-forge ntlink
 ```
 
 Installing using brew:
 ```
 brew install brewsci/bio/ntlink
 ```
 
 Installing from source code:
 ```
curl -L --output ntLink-1.3.8.tar.gz https://github.com/bcgsc/ntLink/releases/download/v1.3.8/ntLink-1.3.8.tar.gz && tar xvzf ntLink-1.3.8.tar.gz 
```

#### Testing your installation
To test your ntLink installation:
```
cd tests
./test_installation.sh
```
The expected output files can be found in: `tests/expected_outputs`

## Dependencies
* Python 3.7+ ([Numpy](https://numpy.org/), [Python-igraph](https://igraph.org/python/))
* [btllib](https://github.com/bcgsc/btllib)
* [ABySS v2.3.0+](https://github.com/bcgsc/abyss)
* GCC 6+ or Clang 5+ with OpenMP
* [zlib](https://zlib.net/)

Python dependencies can be installed with:
```
conda install -c bioconda --file requirements.txt
```

## License
ntLink Copyright (c) 2020 British Columbia Cancer Agency Branch. All rights reserved.

ntLink is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

For commercial licensing options, please contact Patrick Rebstein (prebstein@bccancer.bc.ca).

