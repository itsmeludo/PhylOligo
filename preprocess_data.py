#!/usr/bin/env python
""" pre-processing of reads
"""

import os, sys, argparse

from Bio import SeqIO
from Bio.Seq import Seq

import numpy as np

def get_cmd():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="intputfasta")
    parser.add_argument("-l", action="store", dest="labels")
    parser.add_argument("-s", action="store", dest="sampling", type=float, default=0)
    parser.add_argument("-o", action="store", dest="outputfasta")
    params = parser.parse_args()
    return params


def main():
    params = get_cmd()

    data = list(SeqIO.parse(genome, "fasta"))
    labels = list()
    with open(param.labels) as inf:
        for line in inf:
            labels.append(line.strip())
            
    idx = np.arange(len(labels))
    size = int((len(labels) * params.sampling)/100.)
    with open(params.outputfasta, "w") as outf:
        selected = np.random.choices(idx, size, replace=False)
        for i in selected:
            outf.write(">{}\n{}\n".format(data[i].name, data[i].seq))
          
    return 0

if __name__ == "__main__":
    ret = main()
    sys.exit(ret)
    
    