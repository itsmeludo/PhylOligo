#!/usr/bin/env python3
""" pre-processing of reads

- filter short reads
- (AND/OR) filter read at X percentile of the median length
- (AND/OR) perform a sub sampling of the reads
- (AND/OR) randomize read order

"""

import os, sys, argparse

from Bio import SeqIO
from Bio.Seq import Seq

import numpy as np

def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store",  required=True,dest="inputfasta")
    #parser.add_argument("-l", action="store", dest="labels")
    
    parser.add_argument("-p", action="store", dest="percentile", type=float,
                        help="remove sequences of size not in Xth percentile ")
    parser.add_argument("-m", action="store", dest="min_seqsize", type=int,
                        help="remove sequences shorter than the provided minimal size")
    parser.add_argument("-s", action="store", dest="sampling", type=float, default=0, 
                        help="percentage of read to sample")
    parser.add_argument("-r", action="store_true", dest="randorder", default=False,
                        help="the order of the sequences are randomized")
    
    parser.add_argument("-o", action="store", dest="outputfasta")
    params = parser.parse_args()
    return params

def main():
    params = get_cmd()

    data = list(SeqIO.parse(params.inputfasta, "fasta"))
    idx = np.arange(len(data))
    
    if params.min_readsize!= None:
        # remove short reads
        selected = list()
        for i in idx:
            if len(data[i].seq) > params.min_readsize:
                selected.append(i)
        idx = selected[:]
        
    if params.percentile:
        # remove reads above and below the Xth percentile from median size
        sizes = np.array([len(data[i].seq) for i in idx])
        threshold = params.percentile / 2
        selected = np.argwhere((sizes > np.percentile(sizes,threshold)) *
                               (sizes < np.percentile(sizes, 100-threshold))).flatten()

        idx = selected[:]
    
    if params.sampling != 0:
        # subsampling of read
        size = int((len(data) * params.sampling)/100.)
        selected = np.random.choice(idx, size, replace=False)
        idx = selected[:]
    
    if params.randorder:
        # randomize order
        np.random.shuffle(idx)
            
    with open(params.outputfasta, "w") as outf:
        for i in idx:
            outf.write(">{}\n{}\n".format(data[i].name, data[i].seq))
          
    return 0

if __name__ == "__main__":
    ret = main()
    sys.exit(ret)
    
    