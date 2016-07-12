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
phyloselect.R -i genome.Eucl.mat -a genome.fasta
```


Install
-------

Dependencies
* Python 3.x
 * test
* R 3.x

