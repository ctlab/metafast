metafast
========

***MetaFast*** (METAgenome FAST analysis toolkit) is a toolkit for calculating a number of statistics 
of metagenome sequences and building the distance matrix between them.

Authors:
* **Software:** *Sergey Kazakov* and *Vladimir Ulyantsev*, ITMO University, Saint-Petersburg, Russia.
* **Testing:** *Veronika Dubinkina* and *Alexandr Tyakht*, SRI of Physical-Chemical Medicine, Moscow, Russia.
* **Idea, supervisor:** *Dmitry Alexeev*, SRI of Physical-Chemical Medicine, Moscow, Russia.


### Installation

To run ***metafast*** only one script is required (`metafast.sh`, `metafast.bat` or `metafast.jar`). 
Also you need to have JRE 1.6 or higher installed.

You can download the ***metafast*** run script from the last stable release from <https://github.com/ctlab/metafast/releases>.

* For *Linux* and *Mac OS*: download `metafast.sh`, run the command `chmod a+x metafast.sh`, then run `./metafast.sh` from the command line.
* For *Windows*: download `metafast.bat` and run it from the command line.
* For other OS: download `metafast.jar` and run it via command `java -jar metafast.jar`.


Alternatively, you can build the newest version of ***metafast*** from the repository:
~~~
git clone https://github.com/ctlab/metafast.git
cd metafast 
ant
./out/metafast.sh --version
~~~


### Running ***metafast***

To run ***metafast*** use the following syntax:
* `metafast.sh [<Launch options>] [<Input parameters>]`
* `metafast.bat [<Launch options>] [<Input parameters>]`
* `java -jar metafast.jar [<Launch options>] [<Input parameters>]`

To view possible launch options and input parameters run `metafast.sh --help` or `metafast.sh --help-all`.

By running ***metafast*** a working directory is created (by default `./workDir/`). 
All intermidiate files, log file and final results are saved in it. 

File `output_description.txt` is created after every run in current and working directories. 
It contains the description of every output file produced by the ***metafast***.

Metafast script also allows you to run subtools of whole process or different tools, that was included in the package. 
To see the list of available additional tools, run `metafast.sh --tools`.


### Example

Download [tinytest_A.fastq](https://github.com/ctlab/metafast/raw/master/test_data/tinytest_A.fastq) and [tinytest_B.fastq](https://github.com/ctlab/metafast/raw/master/test_data/tinytest_B.fastq) and run the command:
~~~
./metafast.sh -k 7 -b 0 -l 8 -b1 3 -i tinytest_A.fastq tinytest_B.fastq
~~~

After it has finished, a distance matrix can be found in `workDir/matrices/dist_matrix_<date>_<time>.txt`:
~~~
#	tinytest_A.vec	tinytest_B.vec
tinytest_A.vec	0.0	0.09090909090909091
tinytest_B.vec	0.09090909090909091	0.0
~~~

The element `matrix[i][j]` is a distance between *sample i* and *sample j*.

K-mers frequency statistics is saved in `workDir/kmer-counter-many/stats/<in-file>.stat.txt`.


### Full documentation

To see the full documentation visit <https://github.com/ctlab/metafast/wiki>.


### Contact

Sergey Kazakov - researcher at Computer Technologies Laboratory, ITMO University
Email: <a href="mailto:svkazakov@rain.ifmo.ru">svkazakov@rain.ifmo.ru</a>.


### See also

* [khmer](https://github.com/ged-lab/khmer) - a toolkit to split reads.
* [crAss](http://edwards.sdsu.edu/crass/) - Cross-Assembly of Metagenomes.
* [MaryGold](http://sourceforge.net/projects/metavar/) - Variation analysis of metagenomic samples.

