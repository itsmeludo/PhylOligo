# PhylOligo
===========
Explore oligonucleotide composition similarity between assembly contigs or scaffolds to detect contaminant DNA.



Generate the all-by-all contig distance matrix
----------------------------------------------
```bash
phyloligo.py -d Eucl -i genome.fasta -o genome.Eucl.mat -u 64
```

Regroup contigs by compositional similarity on a tree and explore branching
---------------------------------------------------------------------------

```bash
phyloselect.R -d -v  -i genome.Eucl.mat -a genome.fasta -o genome_conta
```


Install
-------

Dependencies:
* Python 3.x
 * BioPython
 * numpy
 * Cython
* R 3.x
 * ape
 * getopt
 * gtools

In the Bash/Shell, as root/admin if wanted installed globally.
Install python3 and the latest R version [according to your system](https://xkcd.com/1654/) 
```Bash
apt-get install python3-dev python3-setuptools r-base
easy_install3 -U setuptools
pip3 install biopython 
pip3 install cython
pip3 install numpy
```


in R, as root or user
```R
install.packages(c("ape","getopt","gtools"))
```

