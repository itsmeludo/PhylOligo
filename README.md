# PhylOligo
===========
Explore oligonucleotide composition similarity between assembly contigs or scaffolds to detect contaminant DNA.



Generate the all-by-all contig distance matrix
----------------------------------------------
```
phyloligo.py -d Eucl -i genome.fasta -o genome.Eucl.mat -u 64
```

Regroup contigs by compositional similarity on a tree and explore branching
---------------------------------------------------------------------------

```
phyloselect.R -i phyloligo.out -a TH12-rn_prefix.bin
```
