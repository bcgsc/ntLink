all: indexlr

indexlr: btllib/recipes/indexlr.cpp
	btllib/compile
	cp btllib/bin/indexlr .

clean:
	rm -f indexlr

# Check the C++ source code for errors.
lint: clang-format clang-tidy


# Fix white-space errors with clang-format.
fix:
	for i in *.cpp; do clang-format -i -style=file $$i; done

check: check-index-fasta check-index-pos-fasta check-index-pos-strand-fasta check-index-pos-strand-seq-fasta check-index-fastq

check-index-fasta: all
	./indexlr --long -t16 -k100 -w5 --id data/mt.fa | diff -q - data/mt.physlr.tsv.good

check-index-pos-fasta: all
	./indexlr --long -t16 -k100 -w5 --id --pos data/mt.fa | diff -q - data/mt.pos.physlr.tsv.good

check-index-pos-strand-fasta: all
	./indexlr --long -t16 -k100 -w5 --id --pos --strand data/mt.fa | diff -q - data/mt.pos.strand.physlr.tsv.good

check-index-pos-strand-seq-fasta: all
	./indexlr --long -t16 -k100 -w5 --id --pos --strand --seq data/mt.fa | diff -q - data/mt.pos.strand.seq.physlr.tsv.good

check-index-fastq: all
	./indexlr --long -t16 -k100 -w5 --id --bx data/tiny.fq | diff -q - data/tiny.physlr.tsv.good
