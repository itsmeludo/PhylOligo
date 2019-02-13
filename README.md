# PhylOligo
Bioinformatics / Genomics 
Explore oligonucleotide composition similarity between assembly contigs or scaffolds to detect contaminant DNA.


New install with conda:
-----------------------

```Bash
conda install phyloligo -c itsmeludo
```

You might want to install it in an environnement:
```bash
conda create -n phyloligo
conda activate phyloligo
conda install phyloligo -c itsmeludo
```



(Optional) Preprocess the original assembly/raw reads in order to filter out entries, reduce computational time and increase signal
-----------------------------------------------------------------------------------------------------------------------

Filter short sequences or highly conserved repeats.

* Reads an assembly or long sequencing reads multi-fasta file
* Output filtered dataset


```bash
phylopreprocess.py [-h] -i INPUTFASTA [-p PERCENTILE] [-m MIN_SEQSIZE]
                          [-c CUMULATED_SEQSIZE] [-g CUMULATED_PERCENTSIZE]
                          [-s SAMPLING] [-u SAMPLE_SIZE] [-r] [-o OUTPUTFASTA]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFASTA
  -p PERCENTILE         remove sequences of size not in Xth percentile
  -m MIN_SEQSIZE        remove sequences shorter than the provided minimal
                        size
  -c CUMULATED_SEQSIZE  select sequences until their cumulated size reach this
                        parameter. if -r is used, sequences are picked
                        randomly.
  -g CUMULATED_PERCENTSIZE
                        select sequences until their cumulated size reach a
                        percentage of the sequences in the entry file as a
                        whole. if -r is used, sequences are picked randomly.
  -s SAMPLING           percentage of reads to sample
  -u SAMPLE_SIZE        number of reads to sample
  -r                    the order of the sequences are randomized
  -o OUTPUTFASTA

```





Generate the all-by-all contig distance matrix
----------------------------------------------

* Load and index the genome assembly sequences.
* Compute the kmer/spaced-pattern composition profile of each sequence in the assembly.
* Compute a pairwise distance matrix for all sequences. 


```bash
phyloligo.py -d JSD -i genome.fasta -o genome.JSD.mat -c 64
```

Parameters:
```
phyloligo.py [-h] -i GENOME [-k PATTERN] [-s {both,plus,minus}]
                    [-d {Eucl,JSD,KT,BC,SC}] [--freq-chunk-size FREQCHUNKSIZE]
                    [--dist-chunk-size DISTCHUNKSIZE] --method {scoop,joblib}
                    [--large {None,memmap,h5py}] [-c THREADS_MAX]
                    [-o OUT_FILE] [-q OUT_FREQ_FILE] [-w WORKDIR] [-p PATTERN]

optional arguments:
  -h, --help            show this help message and exit
  -i GENOME, --assembly GENOME
                        multifasta of the genome assembly
  -k PATTERN, --lgMot PATTERN
                        word lenght / kmer length / k [default:4]
  -s {both,plus,minus}, --strand {both,plus,minus}
                        strand used to compute microcomposition.
                        [default:both]
  -d {Eucl,JSD,KT,BC,SC}, --distance {Eucl,JSD,KT,BC,SC}
                        how to compute distance between two signatures : Eucl
                        : Euclidean[default:Eucl], JSD : Jensen-Shannon
                        divergence, KT: Kendall's tau, BC: Bray-Curtis,
                        SC:Spearman Correlation
  --freq-chunk-size FREQCHUNKSIZE
                        the size of the chunk to use in scoop to compute
                        frequencies
  --dist-chunk-size DISTCHUNKSIZE
                        the size of the chunk to use in scoop to compute
                        distances
  --method {scoop,joblib}
                        don't use scoop to compute distances use joblib
  --large {None,memmap,h5py}
                        used in combination with joblib for large dataset
  -c THREADS_MAX, --cpu THREADS_MAX
                        how many threads to use for windows microcomposition
                        computation[default:4]
  -o OUT_FILE, --out OUT_FILE
                        output file[default:phyloligo.out]
  -q OUT_FREQ_FILE, --outfreq OUT_FREQ_FILE
                        kmer frequencies output file [default:infile_None]
  -w WORKDIR, --workdir WORKDIR
                        working directory
  -p PATTERN, --pattern PATTERN
                        spaced-word pattern string, only containing 1s and 0s,
                        i.e. '100101001', default='1111'

```



 

Regroup contigs by compositional similarity on a tree and explore branching
---------------------------------------------------------------------------

* Load the distance matrix produced by PhylOligo.
* Optionally create a hierarchically sorted distance matrix.
* Build a cladogram from the distance matrix.
* Interactively ask the user to explore the cladogram and select clads that might correspond to untargeted sequences based on the interpretation of the topology.
* Export clad-specific fasta files:
  * To inspect their potential origin for example with blast or GOHTAM Ménigaud /et al./ 2012
  * To use as learning material in ContaLocate


```bash
phyloselect.R -d -m -c 0.95 -s 4000 -t BIONJ -f c -w 20  -i genome.JSD.mat -a genome.fasta -o genome_conta 

```



Parameters:
```
    -i|--matrix                      
                    All-by-all contig distance matrix, tab separated (required)
    -a|--assembly                    
                    Multifasta file of the contigs (required)
    -f|--tree_draw_method            
                    Tree building type. [phylogram, cladogram,
                    fan, unrooted, radial] by default cladogram.
    -t|--tree_building_method        
                    Tree drawing type [NJ, UPGMA, BIONJ, wardD,
                    wardD2, Hsingle, Hcomplete, WPGMA, WPGMC, UPGMC] by default NJ.
    -m|--matrix_heatmap              
                    Should a matrix heatmap should be produced
    -c|--distance_clip_percentile    
                    Threshold to exclude very distant contigs based on the distance
                    distribution. Use if the tree is squashed by repeats or 
                    degenerated/uninformative contigs [0.97]
    -s|--contig_min_size             
                    Min length in bp of contigs to use in the matrix and tree.
                    Use if the tree is squashed by repeats or 
                    degenerated/uninformative contigs [4000]
    -d|--dump_R_session              
                    Should the R environment be saved for later exploration?
                    The filename will be generated from the outfile parameter
                    or its default value
    -g|--max_perc                    
                    Max edge assembly length percentage displayed (%)
    -l|--min_perc                    
                    Min edge assembly length percentage displayed (%)
    -k|--keep_perc                   
                    Ratio of out-of-range percentages to display (%)
    -o|--outfile                     
                    Outfile name, default:phyloligo.out
    -b|--branchlength                
                    Display branch length
    -w|--branchwidth                 
                    Branch width factor [40]
    -v|--verbose                     
                    Says what the program do.
    -h|--help                        
                    This help.
```

note: PhyloSelect uses the library Ape and its interactive clade selection function on a tree plot with the mouse. X11 is therefore required. If the program has to run on a server -typically for memory reasons- please use the -X option of ssh to allow X11 forwarding.


Regroup contigs by compositional similarity: hierarchical DBSCAN or K-medoids clustering and multidimensional scaling display with t-SNE
----------------------------------------------------------------------------------------------------------------------------------------

* Load the distance matrix produced by PhylOligo.
* Clusterize the sequences
* Export cluster-specific fasta files:
  * To inspect their potential origin for example with blast or GOHTAM Ménigaud /et al./ 2012
  * To use as learning material in ContaLocate

```bash
phyloselect.py -i genome.JSD.mat -t -m hdbscan --noX -o genome_conta


```
Parameters:
```
phyloselect.py [-h] -i DISTMAT [-t] [-p PERPLEXITY] -m
                      {hdbscan,kmedoids} [--minclustersize MIN_CLUSTER_SIZE]
                      [--minsamples MIN_SAMPLES] [-k NBK] [-f FASTAFILE]
                      [--interactive] [--large {memmap,h5py}] [--noX] -o
                      OUTPUTDIR [-q IN_FREQ_FILE]

optional arguments:
  -h, --help            show this help message and exit
  -i DISTMAT            The input matrix file
  -t                    Perform tsne for visualization and pre-clustering
  -p PERPLEXITY         Change the perplexity value
  -m {hdbscan,kmedoids}
                        Method to use to compute cluster on transformed
                        distance matrix
  --minclustersize MIN_CLUSTER_SIZE
                        Set the minimal cluster size of an HDBSCAN cluster
  --minsamples MIN_SAMPLES
                        Set the minimal sample size of an HDBSCAN cluster
  -k NBK                Number of cluster
  -f FASTAFILE          Path of the original fasta file used for the
                        computation of the distance matrix
  --interactive         Allow the user to run the script in an interactive
                        mode and change clustering parameter on the fly
                        (require -t)
  --large {memmap,h5py}
                        Used in combination with joblib for large dataset
  --noX                 Instead of showing pictures, store them in png
  -o OUTPUTDIR
  -q IN_FREQ_FILE, --infreq IN_FREQ_FILE
                        kmer frequencies input file[default:None]. If
                        provided, the clustering is performed on the kmer
                        frequency matrix instead of on the contig distance
                        matrix.

```



Extract DNA segments with homogeneous olignonucleotide composition from a genome assembly using ContaLocate
-----------------------------------------------------------------------------------------------------------

* Learns a compositional profile for the host and the contaminant, previoulsy identified with phyloligo.py / phyloselect.{py,R}.
* Generates a GFF3 map of the contaminant positions in the genome.
 
Once you have explored your assembly's oligonucleotide composition, identified and selected -potentially partial- contaminant material, use ContaLocate to target species-specific contaminant DNA according to a double parametrical threshold.

```bash
contalocate.R -i genome.fasta -r genome_host.fa -c genome_conta_1.fa 
```

The training set for the host genome can be omitted if the amount of contaminant is negligible. In this case, the profile of the host will be trained on the whole genome, including the contaminant.
```bash
contalocate.R -i genome.fasta -c genome_conta_1.fa 
```


The set up of the thresholds can be manually enforced. The user will interactively prompted to set the thresholds given the distribution of windows divergence.
```bash
contalocate.R -i genome.fasta -c genome_conta_1.fa -m
```

Parameters:
```
    -i|--genome              Multifasta of the genome assembly (required)
    -r|--host_learn          Host training set (optional)
    -c|--conta_learn         Contaminant training set (optional) if missing and sliding window parameters are given, the sliding windows composition will be compared to the whole genome composition to contrast potential HGTs (prokaryotes and simple eukaryotes only)
    -t|--win_step            Step of the sliding windows analysis to locate the contaminant (optional) default: 500bp or 100bp
    -w|--win_size            Length of the sliding window to locate the contaminant (optional) default: 5000bp 
    -W|--outputdir           path to outputdir directory
    -d|--dist                Divergence metric used to compare profiles: (KL), JSD or Eucl
    -m|--manual_threshold    You will be asked to manually set the thresholds
    -h|--help                This help

```





Install
-------

A new install from conda is available!
 check it out:

```Bash
conda install phyloligo -c itsmeludo
```


You might want to install it in an environnement:
```bash
conda create -n phyloligo
conda activate phyloligo
conda install phyloligo -c itsmeludo
```

Legacy install procedures as follow ;)

PhylOligo softwares need python 3.4 or newer and several R and python packages.

If python or R are not installed on your system:
```Bash
apt-get install python3-dev python3-setuptools r-base git emboss samtools
```

* Clone/download the git repository.

```Bash
git clone https://github.com/itsmeludo/PhylOligo.git
```
or download it from https://github.com/itsmeludo/PhylOligo

* Install python scripts and dependencies

If you have administrator rights or if you are working in a python virtual environment:

```Bash
git clone https://github.com/itsmeludo/PhylOligo.git
cd PhylOligo
pip3 install .
```

You can also install it locally using:

```Bash
git clone https://github.com/itsmeludo/PhylOligo.git
cd PhylOligo
pip3 install . --user
```

Or to install it locally in a folder of your choice:
```Bash
pip3 install . --prefix /my/local/folder
```


If locally installed, be sure to add the local directory with executable in your executable path.
On linux:

```Bash
export PATH=$HOME/.local/bin:$PATH

phyloligo.py -h
```

* Alternatively python requirements can be installed

If you want to install the dependencies separately use:
```Bash
pip3 install -r requirements.txt
```

* Install R scripts and dependencies

In R, as root or user
```R
install.packages(c("ape","getopt","gplots"))
```

Link the programs into a directory listed in your $PATH
```Bash
ln -s "`pwd`/src/phyloselect.R" <a directory in your $PATH variable>
ln -s `pwd`/src/contalocate.R <a directory in your $PATH variable>
```

* EMBOSS / Samtools 

```Bash
apt-get install emboss samtools
```

* List of Dependencies:

    * Python 3.x
        * [BioPython](biopython.org)
        * [sklearn](http://scikit-learn.org/stable/install.html)
        * [Numpy](numpy.org)
        * [matplotlib](http://matplotlib.org)
        * [hdbscan](https://pypi.python.org/pypi/hdbscan)
        * [Cython](http://cython.org)
        * [h5py](http://www.h5py.org)
    * R 3.x
        * [ape](http://ape-package.ird.fr)
        * [gplots](https://cran.r-project.org/web/packages/gplots/index.html)
        * [getopt](https://cran.r-project.org/web/packages/getopt/getopt.pdf)
    * [EMBOSS](http://emboss.sourceforge.net/download)
    * [Samtools](http://www.htslib.org/)
    * X11 (only required to run phyloselect.R)

    

Test and examples
-----------------

Pipeline example with the M. oryzae fasta file in the data folder.


First example:
* kmer = 4
* Jensen-Shannon divergence (JSD)
* parallelization with joblib using 8 cpu

```Bash
phyloligo.py -i data/M.oryzae_TH12.fasta --method joblib -p 1111 -s both -d JSD -c 8 -w data/example1/ -o data/example1/M.oryzae_TH12_ex1_JSD_joblib.mat
```

Compare the result matrix with the reference matrix
```Bash
phyloligo_comparemat.py --mat1 data/references/M.oryzae_TH12_JSD_ref.mat --format1 numpy --mat2 data/example1/M.oryzae_TH12_ex1_JSD_joblib.mat --format2 numpy
```

Other method can be used to compute the distance matrix:
* for large dataset
```Bash
phyloligo.py -i data/M.oryzae_TH12.fasta --method joblib -p 1111 -s both -d JSD -c 8 -w data/example2/ -o data/example2/M.oryzae_TH12_ex2_JSD_joblib_h5py.mat --large h5py
phyloligo_comparemat.py --mat1 data/references/M.oryzae_TH12_JSD_ref.mat --format1 numpy --mat2 data/example2/M.oryzae_TH12_ex2_JSD_joblib_h5py.mat --format2 h5py
```
or
```Bash
phyloligo.py -i data/M.oryzae_TH12.fasta --method joblib -p 1111 -s both -d JSD -c 8 -w data/example3/ -o data/example3/M.oryzae_TH12_ex3_JSD_joblib_memmap.mat --large memmap
phyloligo_comparemat.py --mat1 data/references/M.oryzae_TH12_JSD_ref.mat --format1 numpy --mat2 data/example3/M.oryzae_TH12_ex3_JSD_joblib_memmap.mat --format2 memmap
```
* using scoop 
```Bash
python -m scoop -n 8 phyloligo.py -i data/M.oryzae_TH12.fasta --method scoop -p 1111 -s both -d JSD -w data/example4/ -o data/example4/M.oryzae_TH12_ex4_JSD_scoop.mat
phyloligo_comparemat.py --mat1 data/references/M.oryzae_TH12_JSD_ref.mat --format1 numpy --mat2 data/example4/M.oryzae_TH12_ex4_JSD_scoop.mat --format2 numpy
```

Scoop can also be used to distribute the worker on a SGE cluster.
However, the cluster must be set to allow ssh connection between nodes. 
Please see with your IT administrator.






