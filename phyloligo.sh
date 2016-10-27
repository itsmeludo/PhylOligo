#!/bin/bash

#$ -S /bin/bash
#$ -q workq
#$ -cwd
#$ -pe parallel_smp 16
#$ -l mem=5G
#$ -l h_vmem=8G


python -m scoop -n 8 phyloligo.py  -i ../PhylOligo_Mor/data/M.oryzae_TH12.fasta -k 5 -d JSD -o ../PhylOligo_Mor/results/M.oryzae_TH12_JSD_scoop2.mat

