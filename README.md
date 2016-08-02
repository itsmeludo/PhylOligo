# PhylOligo
Bioinformatics / Genomics 
Explore oligonucleotide composition similarity between assembly contigs or scaffolds to detect contaminant DNA.


Generate the all-by-all contig distance matrix
----------------------------------------------
```bash
phyloligo.py -d Eucl -i genome.fasta -o genome.Eucl.mat -u 64
```
Parameters:
* -i	--assembly    multifasta of the genome assembly
* -k	--lgMot       word length / kmer length / k
* -s	--strand      Strand used to compute microcomposition. "both", "plus" or "minus"
* -d	--distance    How to compute distance between two signatures "KL", "Eucl" or "JSD" KL: Kullback-Leibler, Eucl : Euclidean, JSD : Jensen-Shannon divergence
* -u	--cpu         How many parallel threads to use for microcomposition computation
* -g    --granularity Factor to refine the granularity of parallel threads. The higher the factor, the greater the number of smaller bins
* -o    --out         Output file
* -h    --help        Exactly what it says
 

Regroup contigs by compositional similarity on a tree and explore branching
---------------------------------------------------------------------------

```bash
phyloselect.R -d -v  -i genome.Eucl.mat -a genome.fasta -o genome_conta
```

Parameters:
* -i    --matrix            All-by-all contig distance matrix, tab separated (required)
* -a    --assembly          Multifasta file of the contigs (required)
* -t    --tree_method       Tree building method (optional), by default NJ. No other option implemented ATM
* -d    --dump_R_session    Should the R environment saved for later exploration? the filename will be generated from the outfile parameter or its default value
* -g    --max_per           Max edge assembly length percentage displayed (%)
* -l    --min_perc          Min edge assembly length percentage displayed (%)
* -k    --keep_perc         Ratio of out-of-range percentages to display (%)
* -o    --outfile           Outfile name, default:phyloligo.out
* -b    --branchlength      Use the branch length information  for tree display if invoked
* -v    --verbos            Say what the program do.
* -h    --help              Yep, does that.

note: PhyloSelect uses the library Ape and its interactive clade selection function on a tree plot with the mouse. X11 is therefore required. If the program has to run on a server -typically for memory reasons- please use the -X option of ssh to allow X11 forwarding.


Install
-------

* Dependencies:
    * Python 3.x
        * [BioPython](biopython.org)
        * [Numpy](numpy.org)
        * [Cython](http://cython.org/)
    * R 3.x
        * [ape](http://ape-package.ird.fr/)
        * [getopt](https://cran.r-project.org/web/packages/getopt/getopt.pdf)
        * [gtools](https://cran.r-project.org/web/packages/gtools/index.html)
    * [EMBOSS](http://emboss.sourceforge.net/download/)
    * [Samtools](http://www.htslib.org/)
    * X11 (only required to run phyloselect)

* Install python3, the latest R version, EMBOSS and Samtools [according to your system](https://xkcd.com/1654/) 

In the Bash/Shell, as root/admin if wanted installed globally.
```Bash
#Ubuntu/Debian-based
apt-get install python3-dev python3-setuptools r-base emboss samtools
easy_install3 -U setuptools
pip3 install biopython 
pip3 install cython
pip3 install numpy
```

in R, as root or user
```R
install.packages(c("ape","getopt","gtools"))
```

* clone the repo

```Bash
git clone https://github.com/itsmeludo/PhylOligo.git
```
or download it from https://github.com/itsmeludo/PhylOligo

* Link the programs into a directory listed in your $PATH

```Bash
cd PhylOligo
sudo ln -s `pwd`/{phyloligo.py,phyloselect.R} /usr/local/bin/
```
