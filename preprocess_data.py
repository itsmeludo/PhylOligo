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

    
          
    return 0

if __name__ == "__main__":
    ret = main()
    sys.exit(ret)
    
    