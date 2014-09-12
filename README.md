metafast
========

Fast metagenome analysis toolkit.

Print help:
~~~ sh
java -jar out/metafast.jar -h
~~~

Print available tools:
~~~ sh
java -jar out/metafast.jar -ts
~~~

Print options for selected tool (kmer-counter, for example):
~~~ sh
java -jar out/metafast.jar -t kmer-counter
~~~

### Authors

Toolkit developed by Vladimir Ulyantsev and Sergey Kazakov (both from ITMO University).
Inspired and supervised by Dmitry Alexeev (NII PCM). 
Experiments and testing - Alexander Tyakht, Veronika Golovanova (NII PCM).

### Examples

Print all options for tool matrix-builder.
~~~ sh
java -jar out/metafast.jar -t matrix-builder
~~~

Count and save binary file with 31-mers from SRR413558merged.fastq. 14 GB given to java virtual machine. 
Locate all files in tmp directory. -v to print all debug information. --force to allow files overwriting.
~~~ sh
java -Xmx14G -jar metafast.jar -t kmer-counter -w tmp -k 31 -i SRR413558merged.fastq -v --force
~~~

Build the distance matrix for input reads from two samples. Input reads are data/mytest_A.fastq and data/mytest_B.fastq.
Use the standard settings (k-mer size is 31 nucleotide, etc.)
~~~ sh
java -jar out/metafast.jar -t matrix-builder -i data/mytest_A.fastq data/mytest_B.fastq
~~~

Build the distance matrix for input reads from two samples. Input reads are data/mytest_A.fastq and data/mytest_B.fastq.
K-mer size is 7 nucleotides, maximal frequency for a k-mer to be assumed erroneous is 0 (all k-mers is good), 
minimum sequence length to be added to a component is 8 nucleotides, -b1 3 - minimum component size is 3 k-mers.
~~~ sh
java -jar out/metafast.jar -t matrix-builder -k 7 -b 0 -l 8 -b1 3 -i data/mytest_A.fastq data/mytest_B.fastq
~~~

### See also

* [khmer](https://github.com/ged-lab/khmer) - toolkit to split your reads
* [crAss](http://edwards.sdsu.edu/crass/) - Cross-Assembly of Metagenomes
* [MaryGold](http://sourceforge.net/projects/metavar/) - Variation analysis of metagenomic samples
