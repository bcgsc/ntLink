![Logo](https://github.com/bcgsc/ntLink/blob/logo/ntlink-logo.png)

## Scaffolding genome assemblies using long reads and minimizers

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

## Usage
```
ntLink: Scaffolding assemblies using long reads
ntLink v1.0.0
Usage: ntLink scaffold target=<target scaffolds> reads='List of long read files'

Options:
target			Target assembly to be scaffolded in fasta format
reads		        List of long read files (separated by a space)
prefix			Prefix of intermediate output files [out.k<k>.w<w>.n<n>]
t			Number of threads [4]
k			K-mer size for minimizers [32]
w			Window size for minimizers [100]
n			Minimum graph edge weight [1]
g			Minimum gap size (bp) [20]
f			Maximum number of contigs in a run for full transitive edge addition [10]
a                       Minimum number of anchored ONT reads required for an edge [1]
z			Minimum size of contig (bp) to scaffold [1000]
v                       If 1, track time and memory for each step of the pipeline [0]
conservative		If False, runs ntLink in stitching mode [True]

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
ntLink scaffold target=my_assembly.fa reads=long_reads.fq.gz k=32 w=500
```

The post-ntLink scaffolds file will have the suffix `*ntLink.scaffolds.fa`

 ## Installation
 Installing from source code:
 ```
git clone https://github.com/bcgsc/ntLink.git
cd src
make
```

#### Testing your installation
To test your ntLink installation:
```
cd tests
./test_installation.sh
```
The expected output files can be found in: `tests/expected_outputs`

## Dependencies
* Python3 ([Numpy](https://numpy.org/), [Python-igraph](https://igraph.org/python/))
* [ABySS](https://github.com/bcgsc/abyss)
* [zlib](https://zlib.net/)

Python dependencies can be installed with:
```
pip3 install -r requirements.txt
```

## License
ntLink Copyright (c) 2020 British Columbia Cancer Agency Branch. All rights reserved.

ntLink is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

For commercial licensing options, please contact Patrick Rebstein (prebstein@bccancer.bc.ca).

