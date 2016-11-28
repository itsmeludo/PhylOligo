#!/usr/bin/env python3
"""
"""
from scoop import futures
import os, sys, re, argparse, pickle, shlex
## clustering
import numpy as np

def KL(a, b):
    """ compute the KL distance
    """
    d = a * np.log(a/b)
    d[np.isnan(d)]=0 
    d[np.isinf(d)]=0
    res = np.sum(d)
    return res

def Eucl(a, b):
    """ compute Euclidean distance
    """
    d = pow(a-b,2)
    d[np.isnan(d)]=0
    d[np.isinf(d)]=0
    return np.sqrt(np.sum(d))

def JSD(a, b):
    """ Compute JSD distance
    """
    h = 0.5 * (a + b)
    d = 0.5 * (KL(a, h) + KL(b, h))
    return d

def compute_JSD_unpack(params):
    """ unpack parameters to compute distances
    """
    i, j, freqi, freqj = params
    d = JSD(freqi, freqj)
    return i, j, d

def compute_Eucl_unpack(params):
    """ unpack parameters to compute distances
    """
    i, j, freqi, freqj = params
    d = Eucl(freqi, freqj)
    return i, j, d

pathin = sys.argv[1]
pathout = sys.argv[2]
dist = sys.argv[3]

with open(pathin, "rb") as inf:
    freqchunk = pickle.load(inf)
    
if __name__ == "__main__":
    if dist == "Eucl":
        res = list(futures.map(compute_Eucl_unpack, freqchunk))
    elif dist == "JSD":
        res = list(futures.map(compute_JSD_unpack, freqchunk))
    else:
        print("Error, unknown metric parameter passed to phylo_batchdist.py: {}".format(dist), file=sys.stderr)
        sys.exit(1)
    with open(pathout, "wb") as outf:
        pickle.dump(res, outf)
    sys.exit(0)
    
    
        