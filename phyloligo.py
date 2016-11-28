#!/usr/bin/env python3
""" This program builds a phylotree based on oligonucleotide profiles (microcomposition) distances for sequences given as arguments.
See the help by calling the program without any argument.

Author: Ludovic V. Mallet, PhD
2016.03.22
licence: GPLv3
Version: Alpha.1
Garanteed with misatkes. <- Including this one.
"""

#dependencies: 
  #Biopython
  #numpy
  #cython
  
from scoop import futures
  
import os, sys, re, argparse, pickle, shlex, shutil
import uuid
import time
import logging
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq

import multiprocessing, subprocess

## clustering
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import check_pairwise_arrays
from sklearn.utils.extmath import row_norms, safe_sparse_dot
from sklearn.utils import gen_even_slices
from sklearn.externals.joblib import Parallel, delayed
from sklearn.externals.joblib import dump, load

#import hdbscan
#from scipy.stats import entropy
#from numpy.linalg import norm

np.seterr(divide='ignore', invalid='ignore')


#### DISTANCES ####



    
#### frequencies computation ####


#### Sequences ####

def get_cmd():
    """ get command line argument
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--assembly", action="store", required=True, dest="genome", 
                        help="multifasta of the genome assembly")
    parser.add_argument("-k", "--lgMot", action="store", dest="k", type=int, default=4, 
                        help="word lenght / kmer length / k [default:%(default)d]")
    parser.add_argument("-s", "--strand", action="store", dest="strand", default="both", choices=["both", "plus", "minus"],
                        help="strand used to compute microcomposition. [default:%(default)s]")
    parser.add_argument("-d", "--distance", action="store", dest="dist", default="Eucl", choices=["Eucl", "JSD"], 
                        help="how to compute distance between two signatures : Eucl : Euclidean[default:%(default)s], JSD : Jensen-Shannon divergence")
    parser.add_argument("--freq-chunk-size", action="store", dest="freqchunksize", type=int, default=250,
                        help="the size of the chunk to use in scoop to compute frequencies")
    parser.add_argument("--dist-chunk-size", action="store", dest="distchunksize", type=int, default=250,
                        help="the size of the chunk to use in scoop to compute distances")
    parser.add_argument("--method", action="store", choices=["scoop1", "scoop2", "joblib"], dest="mthdrun", help="don't use scoop to compute distances use joblib")
    parser.add_argument("--large", action="store_true", dest="large", help="used in combination with joblib for large dataset", default=False)
    parser.add_argument("-c", "--cpu", action="store", dest="threads_max", type=int, default=4, 
                        help="how many threads to use for windows microcomposition computation[default:%(default)d]")
    parser.add_argument("-o", "--out", action="store", dest="out_file", default="phyloligo.out", 
                        help="output file[default:%(default)s]")
    parser.add_argument("-w", "--workdir", action="store", dest="workdir", help="working directory", required=True)
    

    params = parser.parse_args()

    # check sampling parameters
    #if not 0<=params.sampling<=100:
        #print("Error, sampling parameters must be between 0 and 100", file=sys.stderr)
        #sys.exit(1)
        
    params.workdir = os.path.abspath(params.workdir)
        
    return params

def main():
    params = get_cmd()
    
    # compute word frequency of each sequence
    print("Computing frequencies")
    frequencies, freq_name = compute_frequencies(params.mthdrun, params.large,
                                      params.genome, params.k, params.strand, 
                                      params.distchunksize, params.threads_max, params.workdir)
        
    # compute pairwise distances
    print("Computing Pairwise distances")
    res = compute_distances(params.mthdrun, params.large, frequencies, freq_name, params.out_file, 
                            params.dist, params.threads_max, params.freqchunksize, params.workdir)
        
    # save result in a numpy matrix
    if not (params.mthdrun == "joblib" and params.large):
        print("Writing")
        np.savetxt(params.out_file, res, delimiter="\t")
        
    return 0
       
if __name__ == "__main__":
    start = time.time()
    ret = main()      
    stop = time.time()
    print("Exec time: {}".format(stop-start))
    sys.exit(0)
    
      

