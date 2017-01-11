# PhylOligo
Bioinformatics / Genomics 
Explore oligonucleotide composition similarity between assembly contigs or scaffolds to detect contaminant DNA.


Generate the all-by-all contig distance matrix
----------------------------------------------
```bash
phyloligo.py -d JSD -i genome.fasta -o genome.JSD.mat -u 64
```
Parameters:
```
    -i      --assembly          Multifasta of the genome assembly
    -k      --lgMot             word lenght / kmer length / k [default:4]
    -s      --strand            Strand used to compute microcomposition. {both,plus,minus} [default:both]
    -d      --distance          How to compute distance between two signatures : Eucl : Euclidean[default:Eucl], JSD : Jensen-Shannon divergence{Eucl,JSD}          
            --freq-chunk-size   The size of the chunk to use in scoop to compute frequencies
            --dist-chunk-size   The size of the chunk to use in scoop to compute distances
            --method            Which computaional optimisation to use to compute distances {scoop1,scoop2,joblib}            
            --large             Used in combination with joblib for large datasets
    -c      --cpu               How many threads to use for contigs microcomposition computation[default:4]                
    -o      --out               Output file[default:phyloligo.out]
    -w      --workdir           Working directory
    -h      --help              Show this help message and exit
```



 

Regroup contigs by compositional similarity on a tree and explore branching
---------------------------------------------------------------------------

```bash
phyloselect.R -d -m -c 0.95 -s 4000 -t BIONJ -f c -w 20  -i genome.JSD.mat -a genome.fasta -o genome_conta 

```

Parameters:
```
    -i     --matrix                      all-by-all contig distance matrix, tab separated (required)
    -a     --assembly                    multifasta file of the contigs (required)
    -f     --tree_draw_method            tree building type. [phylogram, cladogram, fan, unrooted, radial] by default cladogram.
    -t     --tree_building_method        tree drawing type [NJ, UPGMA, BIONJ, wardD, wardD2, Hsingle, Hcomplete, WPGMA, WPGMC, UPGMC] by default NJ.
    -m     --matrix_heatmap              Should a matrix heatmap should be produced
    -c     --distance_clip_percentile    Threshold to exclude very distant contigs based on the distance distribution. Use if the tree is squashed by repeats or degenerated/uninformative contigs [0.97]
    -s     --contig_min_size             Min length in bp of contigs to use in the matrix and tree. Use if the tree is squashed by repeats or degenerated/uninformative contigs [4000]
    -d     --dump_R_session              Should the R environment be saved for later exploration? the filename will be generated from the outfile parameter or its default value
    -g     --max_perc                    max edge assembly length percentage displayed (%)
    -l     --min_perc                    min edge assembly length percentage displayed (%)
    -k     --keep_perc                   ratio of out-of-range percentages to display (%)
    -o     --outfile                     outfile name, default:phyloligo.out
    -b     --branchlength                display branch length
    -w     --branchwidth                 Branch width factor [40]
    -v     --verbose                     say what the program do. Not implemented yet.
    -h     --help                        Yep, does that.
```

note: PhyloSelect uses the library Ape and its interactive clade selection function on a tree plot with the mouse. X11 is therefore required. If the program has to run on a server -typically for memory reasons- please use the -X option of ssh to allow X11 forwarding.


Regroup contigs by compositional similarity: hierarchical DBSCAN and MDS display with t-SNE
-------------------------------------------------------------------------------------------

```bash
phyloselect.py -i genome.JSD.mat -t -m hdbscan --noX -o genome_conta


```
Parameters:
```
    -i                      The input matrix file
    -t                      Perform t-SNE for visualization and pre-clustering
    -p                      Change the perplexity value for t-SNE
    -m                      Method to use to compute cluster on transformed distance matrix{hdbscan,kmedoids}
           --minclustersize Set the minimal cluster size of an HDBSCAN cluster
           --minsamples     Set the minimal sample size of an HDBSCAN cluster
    -k                      Number of clusters if kmedoids is used
    -f                      Path of the original fasta file used for the computation of the distance matrix                 
           --interactive    Allow the user to run the script in an interactive mode and change clustering parameter on the fly (require -t)
           --large          Used in combination with joblib for large dataset
           --noX            Instead of showing pictures, store them in pdf
    -o                      OUTPUTDIR
    -h      --help          Show this help message and exit
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
    -i    --genome            Multifasta of the genome assembly (required)
    -r    --host_learn        Host training set (optional)
    -c    --conta_learn       Contaminant training set (optional) if missing and sliding window parameters are given, the sliding windows composition will be compared to the whole genome composition to contrast potential HGTs (prokaryotes and simple eukaryotes only)
    -t    --win_step          Step of the sliding windows analysis to locate the contaminant (optional) default: 500bp or 100bp
    -w    --win_size          Length of the sliding window to locate the contaminant (optional) default: 5000bp 
    -d    --dist              Divergence metric used to compare profiles: (KL), JSD or Eucl
    -m    --manual_threshold  You will be asked to manually set the thresholds
    -h    --help              What it says
```


Install
-------

PhyloOligo softwares need python 3.4 or newer and several R and python packages.

If python or R are not installed on your system:
```Bash
apt-get install python3-dev python3-setuptools r-base
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

Or to install in locally in a folder of your choice:
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






