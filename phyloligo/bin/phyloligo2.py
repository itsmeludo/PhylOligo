#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# START LICENCE ##############################################################
#
# Cluster read based on oligonucleotide profiles (microcomposition) distances 
# Copyright (C) 2016  Ludovic Mallet, Tristan Bitard-Feidel
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# END LICENCE ##############################################################
""" This program cluster sequences based on oligonucleotide profiles (microcomposition).
See the help by calling the program without any argument.
"""

__author__  = "Tristan Bitard-Feildel, Ludovic Mallet"
__email__   = "tristan.bitard-feildel@impmc.upmc.fr, Ludovic.mallet@inrafr"
__year__    = 2016
__licence__ = "GPLv3"
__version__ = 0.1
  
import os, sys, argparse, time
#from phyloligo import phylodist, phylofreq
from phyloligo import compute_distances, compute_frequencies
import numpy as np

np.seterr(divide='ignore', invalid='ignore')


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

