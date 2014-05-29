metafast
========

Fast metagenome analysis toolkit.

Print available tools:
~~~ sh
java -jar metafast.jar -ts
~~~

Print options for selected tool (kmer-counter, for example)
~~~ sh
java -jar metafast.jar -t kmer-counter
~~~

### Authors

Toolkit developed by Vladimir Ulyantsev (ITMO University); inspired and supervised by Dmitry Alexeev (NII PCM). Experiments and testing - Alexander Tyakht, Veronika Golovanova (NII PCM).

### Examples

Count and save binary file with 31-mers from SRR413558merged.fastq. 14 GB given to java virtual machine. Locate all files  in /tmpWD. -v to print all debug information. --force to allow files overwriting.
~~~ sh
java -Xmx14G -jar metafast.jar -t kmer-counter -k 31 -i SRR413558merged.fastq -v --force -w tmpWD
~~~

### See also

* [khmer](https://github.com/ged-lab/khmer) - toolkit to split your reads
* [crAss](http://edwards.sdsu.edu/crass/) - Cross-Assembly of Metagenomes
* [MaryGold](http://sourceforge.net/projects/metavar/) - Variation analysis of metagenomic samples
