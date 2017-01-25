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
__email__   = "tristan.bitard-feildel@impmc.upmc.fr, Ludovic.mallet@inra.fr"
__year__    = 2016
__licence__ = "GPLv3"
__version__ = 0.1
  
from scoop import futures

import os, sys, argparse, re
import time, tempfile
import numpy as np
import shlex, subprocess, pickle, shutil
import h5py

import Bio.Cluster
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from itertools import product

import numpy as np
from sklearn.externals.joblib import Parallel, delayed
from sklearn.externals.joblib import dump, load
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import check_pairwise_arrays
from sklearn.utils import gen_even_slices

# TODO
# to be tested and added as option
SCALING = 1

np.seterr(divide='ignore', invalid='ignore')

#### Utilities

def remove_folder(folder):
    try:
        shutil.rmtree(folder)
    except:
        print("Failed to delete folder: {}".format(folder))

def read_seq_chunk_pos(genome, chunksize):
    """ read a first chunk of fasta sequences to be processed
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    chunksize: int
        the number of fasta entries to read
    
    Return:
    -------
    seqchunk: list
        a list of SeqRecord 
    """
    seqchunk = list()
    start = 0
    for record in SeqIO.parse(genome, "fasta"):
        seqchunk.append(str(record.seq))
        if len(seqchunk) == chunksize:
            yield start, start + chunksize, seqchunk
            start += chunksize
            seqchunk = list()
    # the last chunk 
    if seqchunk != []:
        yield start, start+len(seqchunk), seqchunk

def read_seq_chunk(genome, chunksize, pattern, strand):
    """ read a first chunk of fasta sequences to be processed
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    chunksize: int
        the number of fasta entries to read
    
    Return:
    -------
    seqchunk: list
        a list of SeqRecord 
    """
    seqchunk = list()
    pattern=str(pattern)
    for record in SeqIO.parse(genome, "fasta"):
        seqchunk.append((str(record.seq), pattern, strand))
        if len(seqchunk) == chunksize:
            yield seqchunk
            seqchunk = list()
    # the last chunk 
    if seqchunk != []:
        yield seqchunk


def select_strand (seq, strand="both"):
    """ choose which strand to work on
    
    Parameters:
    -----------
    seq: Seq
        a str object containing a nucleotide sequence, that can be converted into a Bio.Seq object
    strand: string
        select wich strand to use
    
    Return:
    -------
    seq: string
        the sequence strand
    """
    Bioseq_record = Seq(seq)
    if (strand == "both"):
        new_seq = str(str(seq)+str(Bioseq_record.reverse_complement()))
    elif (strand == "minus"):
        new_seq = str(Bioseq_record.reverse_complement())
    elif (strand == "plus"):
        new_seq = str(seq)
    else:
        print("Error, strand parameter of selectd_strand() should be choose from {'both', 'minus', 'plus'}", file=sys.stderr)
        sys.exit(1)
    return new_seq

def get_nb_records(fasta, file_format="fasta"):
    """ compute the number of sequences in a fasta file
    """
    return sum(1 for record in SeqIO.parse(fasta, file_format))


#### Compute distances 

## clustering

def posdef_check_value(d):
    d[np.isnan(d)]=0    
    d[np.isinf(d)]=0
    
## Distance functions

def KL(a, b):
    """ compute the KL distance
    """
    if a.ndim == 1 and b.ndim == 1:
        d = a * np.log(a/b)
        posdef_check_value(d)
        res = np.sum(d)
    elif a.ndim == 2 and b.ndim == 2:
        X, Y = check_pairwise_arrays(a, b)
        X = X[:,np.newaxis]
        d = X * np.log(X/Y)
        posdef_check_value(d)
        res = np.sum(d, axis=2).T
    else:
        print("Dimension erro in KL, a={}, b={}".format(a.ndim, b.ndim), file=sys.stderr)
        sys.exit(1)
    return res

def Eucl(a, b):
    """ compute Euclidean distance
    """
    d = pow(a-b,2)
    posdef_check_value(d)
    return np.sqrt(np.sum(d))

def JSD(a, b):
    """ Compute JSD distance
    """
    if a.ndim == 1 and b.ndim == 1:
        h = 0.5 * (a + b)
        d = 0.5 * (KL(a, h) + KL(b, h))
    elif a.ndim==2 and b.ndim == 1:
        h = 0.5 * (a[np.newaxis,:] + b)
        d1 = a[np.newaxis,:] * np.log(a[np.newaxis,:]/h)
        posdef_check_value(d1)
        d1 = np.sum(d1, axis=2)
        d2 = b * np.log(b/h)
        posdef_check_value(d2)
        d2 = np.sum(d2, axis=2)
        d = 0.5 * (d1 + d2)
    else:
        h = 0.5 * (a[np.newaxis,:] + b[:, np.newaxis])
        d1 = a[np.newaxis,:] * np.log(a[np.newaxis,:]/h)
        posdef_check_value(d1)
        d1 = np.sum(d1, axis=2)
        d2 = b[:, np.newaxis] * np.log(b[:, np.newaxis]/h)
        posdef_check_value(d2)
        d2 = np.sum(d2, axis=2)
        d = 0.5 * (d1 + d2)
        #d = d.T
    return d


def KT(a, b):
    """ Compute Kendall_Tau distance
    """
    return(1 - Bio.Cluster.distancematrix((a,b), dist="k")[1][0])


def weighted_rank(a, b):
    return 0

unpack_distances = {"JSD": JSD,
                    "Eucl": Eucl,
                    "KT": KT,
                    }

def compute_unpack(params):
    """ unpack parameters to compute distances
    """
    i, j, freqi, freqj, distname = params
    d = unpack_distances[distname](freqi, freqj)
    return i, j, d

#def compute_JSD_unpack(params):
    #""" unpack parameters to compute distances
    #"""
    #i, j, freqi, freqj = params
    #d = JSD(freqi, freqj)
    #return i, j, d

#def compute_KT_unpack(params):
    #""" unpack parameters to compute distances
    #"""
    #i, j, freqi, freqj = params
    #d = KT(freqi, freqj)
    #return i, j, d

#def compute_Eucl_unpack(params):
    #""" unpack parameters to compute distances
    #"""
    #i, j, freqi, freqj = params
    #d = Eucl(freqi, freqj)
    #return i, j, d


def distances_loc(output, X, s, metric):
    """ quick wrapper around function call
    """
    loc_h5py[metric](output, X, s)

def euclidean_distances_loc(output, X, s):
    distances = euclidean_distances(X, X[s])
    output[s] = distances.T
    
def JSD_loc(output, X, s):
    X, Y = check_pairwise_arrays(X, X[s])
    d = JSD(X, Y)
    output[s] = d
    
def KT_loc(output, X, s):
    X, Y = check_pairwise_arrays(X, X[s])
    d = KT(X, Y)
    output[s] = d   
    
loc_h5py = {"Eucl": euclidean_distances_loc,
             "JSD": JSD_loc,
             "KT": KT_loc}



def distances_h5py(output_dir, input, s, metric):
    """ quick wrapper around various distances
    """
    dist_h5py[metric](output_dir, input, s)

def euclidean_distances_h5py(output_dir, input, s):
    with h5py.File(input, "r") as hf:
        X = hf.get("frequencies")
        dist = euclidean_distances(X, X[s]).T
    
    output = os.path.join(output_dir, "distance_{}_{}".format(s.start, s.stop))
    with h5py.File(output, "w") as hf:
        distances = hf.create_dataset("distances", dist.shape, dtype="float32")
        distances[...] = dist[:]
    

def JSD_h5py(output_dir, input, s):
    with h5py.File(input, "r") as hf:
        X = hf.get("frequencies")
        X, Y = check_pairwise_arrays(X, X[s])
        dist = JSD(X, Y)
    
    output = os.path.join(output_dir, "distance_{}_{}".format(s.start, s.stop))
    with h5py.File(output, "w") as hf:
        distances = hf.create_dataset("distances", dist.shape, dtype="float32")
        distances[...] = dist[:]
    #hf = h5py.File(output, "r+")
    #distances = hf.get("distances")
    #distances[s] = d
    #hf.close()
    
def KT_h5py(output_dir, input, s):
    with h5py.File(input, "r") as hf:
        X = hf.get("frequencies")
        X, Y = check_pairwise_arrays(X, X[s])
        dist = KT(X, Y)
    
    output = os.path.join(output_dir, "distance_{}_{}".format(s.start, s.stop))
    with h5py.File(output, "w") as hf:
        distances = hf.create_dataset("distances", dist.shape, dtype="float32")
        distances[...] = dist[:]
    #hf = h5py.File(output, "r+")
    #distances = hf.get("distances")
    #distances[s] = d
    #hf.close()   

dist_h5py = {"Eucl": euclidean_distances_h5py,
             "JSD": JSD_h5py,
             "KT": KT_h5py}

def compute_distances_scoop(frequencies, chunksize, metric="Eucl"):
    """ compute pairwises distances
    
    Parameters
    ----------
    frequencies: np.array
        a list of frequencies, each row corresponding to a sample each column to a "word"
    metric: string
        distance method to use ('JSD', 'Eucl', 'KT')
    
    Return:
    -------
    distances: np.array
        (n_samples, n_samples) distance matrix
    """
    #scoop.logger.info("Starting distance computation")
    
    if metric not in ["JSD", "Eucl", "KT"]:
        # should already have been verified, but just in case ...
        print("Error, unknown method {}".format(metric), file=sys.stderr)
        sys.exit(1)
        
    distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
    for freqchunk in make_freqchunk(frequencies, metric, chunksize):
        res = futures.map(compute_unpack, freqchunk)
        for i, j, d in res:
            distances[i, j] = distances[j, i] = d
            
    #if metric == "Eucl":
        #distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
        #for freqchunk in make_freqchunk(frequencies, chunksize):
            #res = futures.map(compute_Eucl_unpack, freqchunk)
            #for i, j, d in res:
                #distances[i, j] = distances[j, i] = d
                
    #elif metric == "JSD":
        #distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
        #for freqchunk in make_freqchunk(frequencies, chunksize):
            #res = futures.map(compute_JSD_unpack, freqchunk)
            #for i, j, d in res:
                #distances[i, j] = distances[j, i] = d
                
    #elif metric == "KT":
        #distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
        #for freqchunk in make_freqchunk(frequencies, chunksize):
            #res = futures.map(compute_KT_unpack, freqchunk)
            #for i, j, d in res:
                #distances[i, j] = distances[j, i] = d
                
    return distances

def compute_distances_joblib(frequencies, metric="Eucl", n_jobs=1):
    """ compute pairwises distances
      
    Parameters
    ----------
    frequencies: np.array
        a list of frequencies, each row corresponding to a sample each column to a "word"
    metric: string
        distance method to use ('JSD', 'Eucl', 'KT')
    n_jobs: int
        number of parallel job to execute
    
    Return:
    -------
    distances: np.array
        (n_samples, n_samples) distance matrix
    """
    call_dist = {"Eucl": Eucl, "JSD": JSD}
    
    if not callable(metric) and metric not in call_dist:
        print("Error, unknown metric methodfor joblib: {}".format(metric), file=sys.stderr)
        sys.exit(1)
    
    if callable(metric):
        distances = pairwise_distances(frequencies, metric=metric, n_jobs=n_jobs)
    else:
        distances = pairwise_distances(frequencies, metric=call_dist[metric], n_jobs=n_jobs)
    
    return distances

def compute_distances_memmap(frequencies, freq_name, output, metric="Eucl", n_jobs=1):
    """ compute pairwises distances
    
    Parameters
    ----------
    freq_name: string
        path of the local frequency array
    metric: string
        distance method to use ('JSD', 'Eucl', 'KL')
    n_jobs: int
        number of parallel job to execute
    
    Return:
    -------
    distances: np.array
        (n_samples, n_samples) distance matrix
    """

    # Pre-allocate a writeable shared memory map as a container for the results of the parallel computation
    distances = np.memmap(output, dtype=frequencies.dtype, shape=(frequencies.shape[0], frequencies.shape[0]), mode='w+')
    
    # close and reopen it for reference
    del distances
    distances = np.memmap(output, dtype=frequencies.dtype, shape=(frequencies.shape[0], frequencies.shape[0]), mode='r+')

    if metric not in ["Eucl", "JSD", "KT"]:
        print("Error, unknown method {}".format(metric), file=sys.stderr)
        sys.exit(1)
        
    fd = delayed(distances_loc)
    Parallel(n_jobs=n_jobs, verbose=0)(fd(distances, frequencies, s, metric) for s in gen_even_slices(frequencies.shape[0], n_jobs*SCALING))
        
    folder = os.path.dirname(freq_name)
    remove_folder(folder)
    
    #if metric == "Eucl":
        ## execute parallel computation of euclidean distance
        #fd = delayed(euclidean_distances_loc)
        #Parallel(n_jobs=n_jobs, verbose=0)(fd(distances, frequencies, s) for s in gen_even_slices(frequencies.shape[0], n_jobs*SCALING))
            
        #folder = os.path.dirname(freq_name)
        #remove_folder(folder)
        ##print(freq_name)
        
    #elif metric == "JSD":
        #fd = delayed(JSD_loc)
        #Parallel(n_jobs=n_jobs, verbose=0)(fd(distances, frequencies, s) for s in gen_even_slices(frequencies.shape[0], n_jobs*SCALING))
        
        #folder = os.path.dirname(freq_name)
        #remove_folder(folder)
        
    #elif metric == "KT":
        #fd = delayed(KT_loc)
        #Parallel(n_jobs=n_jobs, verbose=0)(fd(distances, frequencies, s) for s in gen_even_slices(frequencies.shape[0], n_jobs*SCALING))
        
        #folder = os.path.dirname(freq_name)
        #remove_folder(folder)
            
    #else:
        #print("Error, unknown method {}".format(metric), file=sys.stderr)
        #sys.exit(1)

def join_distance_results(dist_folder, dist_name, size, n_jobs):
    """ join the different distances computation into a single matrix
    
    Parameters
    ----------
    dist_folder: string
        path to the distance folder storing all distances computation
    dist_name: string
        path to the final matrix
    size: int
        size of the final matrix
    n_jobs: int
        number of jobs used, the names of the files storing distances 
        depend on the number of jobs and on the size of the matrix
    """
    with h5py.File(dist_name, "w") as hf:
        distances = hf.create_dataset("distances", (size, size), dtype="float32")
        for res in os.listdir(dist_folder):
            tmp = res.split("_")
            start, stop = int(tmp[1]), int(tmp[2])
            pathin = os.path.join(dist_folder, "distance_{}_{}".format(start, stop))
            with h5py.File(pathin, "r") as inf:
                distances[start: stop] = inf.get("distances").value[:]

def compute_distances_h5py(freq_name, dist_name, metric="Eucl", n_jobs=1):
    """ compute pairwises distances
    
    Parameters
    ----------
    freq_name: string
        path of the local frequency array
    metric: string
        distance method to use ('JSD', 'Eucl', 'KL')
    n_jobs: int
        number of parallel job to execute
    
    Return:
    -------
    distances: np.array
        (n_samples, n_samples) distance matrix
    """
    #scoop.logger.info("Starting distance computation")
    #folder = tempfile.f()
    #dist_name = os.path.join(folder, output)
    
    folder = os.path.dirname(freq_name)
    
    with h5py.File(freq_name, "r") as hf:
        frequencies = hf.get("frequencies")
        size = frequencies.shape[0]
    
    dist_folder = os.path.join(folder, "distances")
    if not os.path.isdir(dist_folder):
        os.makedirs(dist_folder)
        
    if metric not in ["Eucl", "JSD", "KT"]:
        print("Error, unknown method {}".format(metric), file=sys.stderr)
        sys.exit(1)
        
    fd = delayed(distances_h5py)
    Parallel(n_jobs=n_jobs, verbose=0)(fd(dist_folder, freq_name, s, metric) for s in gen_even_slices(size, n_jobs*SCALING)) # some scaling
    
    #if metric == "Eucl":
        ## execute parallel computation of euclidean distance
        #fd = delayed(euclidean_distances_h5py)
        #Parallel(n_jobs=n_jobs, verbose=0)(fd(dist_folder, freq_name, s) for s in gen_even_slices(size, n_jobs*SCALING)) # some scaling
        
    #elif metric == "JSD":
        #fd = delayed(JSD_h5py)
        #Parallel(n_jobs=n_jobs, verbose=0)(fd(dist_folder, freq_name, s) for s in gen_even_slices(size, n_jobs*SCALING)) # some scaling
   
    #elif metric == "KT":
        #fd = delayed(KT_h5py)
        #Parallel(n_jobs=n_jobs, verbose=0)(fd(dist_folder, freq_name, s) for s in gen_even_slices(size, n_jobs*SCALING)) # some scaling
           
    # join distance results
    join_distance_results(dist_folder, dist_name, size, n_jobs)
    # remove folder storing separated distances
    remove_folder(folder)

def compute_distances(mthdrun, large, frequencies, freq_name, out_file, dist, threads_max, freqchunksize, workdir):
    """ choose which function to call to compute distances
    """
    res=None
    if mthdrun == "joblib":
        if large == "memmap":
            #res = compute_distances_memmap(frequencies, freq_name, out_file, metric=dist, n_jobs=threads_max)
            compute_distances_memmap(frequencies, freq_name, out_file, metric=dist, n_jobs=threads_max)
        if large == "h5py":
            #res = compute_distances_h5py(freq_name, out_file, metric=dist, n_jobs=threads_max)
            compute_distances_h5py(freq_name, out_file, metric=dist, n_jobs=threads_max) 
        else:
            res = compute_distances_joblib(frequencies, metric=dist, n_jobs=threads_max)
    elif mthdrun == "scoop":
        res = compute_distances_scoop(frequencies, freqchunksize, metric=dist)
    else:
        print("Error, method {} is not implemented for pairwise distances computation".format(mthdrun), file=sys.stderr)
    return res

#### Compute frequencies

def make_freqchunk(frequencies, metric, chunksize):
    """ prepare frequencies to be parallelized
    """
    chunk = list()
    for i in range(len(frequencies)):
        freqi = frequencies[i]
        for j in range(i, len(frequencies)):
            freqj = frequencies[j]
            chunk.append((i, j, freqi, freqj, metric))
            if len(chunk) == chunksize:
                yield chunk
                chunk = list()
    # the last chunk 
    if chunk != []:
        yield chunk

# TODO to be removed
#def cut_sequence_and_count(seq, ksize):
    #""" cut sequence in words of size k
    
    #Parameters:
    #-----------
    #seq: string
        #The nucleotide sequence
    #ksize: int
        #The size of the kmer
    
    #Return:
    #-------
    #count_words: dict
        #a dictionary storing the number of observations of each word
    #kword_count: int
        #the total number of words
    #"""
    #seq_words = list()
    ## re.split: excludes what is not a known characterised nucleotide 
    #for subseq in re.split('[^ACGT]+', seq): 
        #if (len(subseq) >= ksize):
            ##get all kword in sub-sequence
            #seq_words.extend(subseq[i: i+ksize] for i in range(len(subseq)-(ksize-1)))  
    #count_words = Counter(seq_words)
    #kword_count = sum(count_words.values())
    #return count_words, kword_count

def cut_sequence_and_count_pattern(seq, pattern):
    """ cut sequence in spaced-words of size k 
    
    Parameters:
    -----------
    seq: string
        The nucleotide sequence
    pattern: string
        the binary space pattern to extract spaced-words example: 1001010001 
        ksize is inferred from the number of '1' in the pattern
        
    Return:
    -------
    count_words: dict
        a dictionary storing the number of observations of each word
    kword_count: int
        the total number of words
    """
    seq_words = list()
    #ksize=pattern.count('1')
    pattern=str(pattern)
    target_index = [i  for i,j in enumerate(pattern) if j=="1"]

    # re.split: excludes what is not a known characterised nucleotide 
    for subseq in re.split('[^ACGT]+', seq): 
        if (len(subseq) >= len(pattern)):
            #get all kword in sub-sequence
            seq_words.extend("".join(list(map( lambda x: subseq[i+x], target_index))) for i in range(len(subseq)-(len(pattern)-1)))
    count_words = Counter(seq_words)
    kword_count = sum(count_words.values())
    return count_words, kword_count

def count2freq(count_words, kword_count, ksize):
    """ transform raw count into a feature vector of frequencies
    
    Parameters:
    -----------
    count_words: dict
        a dictionary storing the number of observations of each word
    kword_count: int
        the total number of words
    ksize: int
        the size of the kmer
        
    Return:
    -------
    features: np.array
        a feature vector of the frequency of each word in the read
    """
    features = list()
    if kword_count > 0:
        # iterate over letter product
        for letter_word in product(("C","G","A","T"), repeat=ksize):
            kword = "".join(letter_word)
            if kword in count_words:
                features.append(count_words[kword]/kword_count)
            else:
                features.append(0)
    else:
        features = [0 for kword in range(4**ksize)]
    return np.array(features)

def compute_frequency(seq, pattern="1111", strand="both"):
    """ compute kmer frequency, ie feature vector of each read
    
    Parameters:
    -----------
    seq: string
        The nucleotide sequence
    pattern: string
        the binary space pattern to extract spaced-words example: 1001010001 
        ksize is inferred from the number of '1' in the pattern
    strand: string
        which strand to used
        
    Return:
    -------
    features: np.array
        a feature vector of the frequency of each word in the read
    """
    seq = select_strand(seq, strand)
    # we work on upper case letters
    seq = seq.upper()
    # raw word count
    #count_words, kword_count = cut_sequence_and_count(seq, ksize)
    pattern=str(pattern)
    ksize = pattern.count("1")
    count_words, kword_count = cut_sequence_and_count_pattern(seq, pattern)
    # create feature vector
    features = count2freq(count_words, kword_count, ksize)
    return features

def compute_frequency_memmap(frequency, i, seq, pattern="1111", strand="both"):
    """ function used in joblib to write in the memmap array frequency
    
    Parameters:
    -----------
    frequency: np.memmap
        a memmap array to store computed frequencies
    i: int 
        index of the sequence
    seq: string
        The nucleotide sequence
    pattern: string
        the binary space pattern to extract spaced-words example: 1001010001 
        ksize is inferred from the number of '1' in the pattern
    strand: string
        which strand to used
        
    """
    seq = select_strand(seq, strand)
    # we work on upper case letters
    seq = seq.upper()
    # raw word count
    #count_words, kword_count = cut_sequence_and_count(seq, ksize)
    pattern=str(pattern)
    ksize = pattern.count("1")
    count_words, kword_count = cut_sequence_and_count_pattern(seq, pattern)
    # create feature vector
    frequency[i] = count2freq(count_words, kword_count, ksize)
    
def compute_frequency_h5py(freq_name_folder, i, seq, pattern="1111", strand="both"):
    """ function used in joblib to write in the h5py array frequency
    
    Parameters:
    -----------
    freq_name: string
        path to the folder file holding h5py frequencies
    i: int 
        index of the sequence
    seq: string
        The nucleotide sequence
    pattern: string
        the binary space pattern to extract spaced-words example: 1001010001 
        ksize is inferred from the number of '1' in the pattern
    strand: string
        which strand to used
        
    """
    seq = select_strand(seq, strand)
    # we work on upper case letters
    seq = seq.upper()
    # raw word count
    #count_words, kword_count = cut_sequence_and_count(seq, ksize)
    pattern=str(pattern)
    ksize = pattern.count("1")
    count_words, kword_count = cut_sequence_and_count_pattern(seq, pattern)
    # create feature vector
    freq_name = os.path.join(freq_name_folder, "frequencies_{}".format(i))
    with h5py.File(freq_name, "w") as hf:
        frequency = hf.create_dataset("frequencies", (4**ksize,), dtype="float32")
        frequency[...] = count2freq(count_words, kword_count, ksize)[:]
        #print(frequency)
    
    
def compute_frequency_h5py_chunk(freq_name_folder, seqchunk, pattern, strand, start, stop):
    """ function used in joblib to write in the h5py array frequency
    
    Parameters:
    -----------
    freq_name: string
        path to the folder file holding h5py frequencies
    i: int 
        index of the sequence
    seq: string
        The nucleotide sequence
    pattern: string
        the binary space pattern to extract spaced-words example: 1001010001 
        ksize is inferred from the number of '1' in the pattern
    strand: string
        which strand to used
        
    """
    pattern=str(pattern)
    ksize = pattern.count("1")
    size = stop-start
    freqs = np.empty((size, 4**ksize), dtype="float32")
    for i, seq in enumerate(seqchunk):
        seq = select_strand(seq, strand)
        # we work on upper case letters
        seq = seq.upper()
        # raw word count
        #count_words, kword_count = cut_sequence_and_count(seq, ksize)
        count_words, kword_count = cut_sequence_and_count_pattern(seq, pattern)
        # create feature vector
        freqs[i] = count2freq(count_words, kword_count, ksize)[:]
    freq_name = os.path.join(freq_name_folder, "frequencies_{}_{}".format(start, stop))
    #print(size, 4**ksize)
    #print(freqs.shape)
    with h5py.File(freq_name, "w") as hf:
        frequency = hf.create_dataset("frequencies", (size, 4**ksize), dtype="float32")
        frequency[...] = freqs[:]
        #print(frequency)
    
def frequency_pack(params):
    """ compute kmer frequency, ie feature vector of each read
    
    Parameters:
    -----------
    seq: string
        The nucleotide sequence
    ksize: int
        The size of the kmer
    strand: string
        which strand to used
        
    Return:
    -------
    features: np.array
        a feature vector of the frequency of each word in the read
    """
    features = compute_frequency(*params)
    return features


#### Differents colors of parallelism for compute_frequencies

def compute_frequencies_scoop(genome, pattern, strand, chunksize):
    """ compute frequencies
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    pattern: string
        the binary space pattern to extract spaced-words example: 1001010001 
        ksize is inferred from the number of '1' in the pattern
    strand: string
        select genome strand ('both', 'plus', 'minus')
    chunksize: int
        define chunk of sequences per worker
        
    Return:
    -------
    frequencies: numpy.array
        the samples x features matrix storing NT composition of fasta sequences
    """
    #scoop.logger.info("Starting frequencies computation")
    frequencies = list()
    for seqchunk in read_seq_chunk(genome, chunksize, pattern, strand):
        chunkfreq = futures.map(frequency_pack, seqchunk)
        frequencies.extend(chunkfreq)
                           
    frequencies = np.array(frequencies)
    return frequencies

def compute_frequencies_joblib(genome, pattern, strand, nbthread):
    """ compute frequencies
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    pattern: string
        the binary space pattern to extract spaced-words example: 1001010001 
        ksize is inferred from the number of '1' in the pattern
    strand: string
        select genome strand ('both', 'plus', 'minus')
    workdir: strng
        temporary working directory to store pickled chunk
        
    Return:
    -------
    frequencies: numpy.array
        the samples x features matrix storing NT composition of fasta sequences
    """
    fd = delayed(compute_frequency)
    frequencies = Parallel(n_jobs=nbthread, verbose=0)(
        fd(str(record.seq), pattern, strand) for record in SeqIO.parse(genome, "fasta"))
    
    #for ar in frequencies:
        #print(ar.shape)
    frequencies = np.vstack(frequencies)
    #print("Frequency shape: {}".format(frequencies.shape))
    if len(frequencies.shape) < 2:
        frequencies.reshape(-1, 1)
    return frequencies

def compute_frequencies_joblib_memmap(genome, pattern, strand, nbthread, workdir):
    """ compute frequencies
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    pattern: string
        the binary space pattern to extract spaced-words example: 1001010001 
        ksize is inferred from the number of '1' in the pattern
    strand: string
        select genome strand ('both', 'plus', 'minus')
    workdir: strng
        temporary working directory to store pickled chunk
        
    Return:
    -------
    frequencies: numpy.array
        the samples x features matrix storing NT composition of fasta sequences
    """
    folder = tempfile.mkdtemp(dir=workdir)
    freq_name = os.path.join(folder, 'frequencies')
    pattern=str(pattern)
    ksize = pattern.count("1")
    # Pre-allocate a writeable shared memory map as a container for the frequencies
    nb_record = get_nb_records(genome)
    frequencies = np.memmap(freq_name, dtype=np.float32, shape=(nb_record, 4**ksize), mode='w+')
    #print(freq_name)
    del frequencies
    frequencies = np.memmap(freq_name, dtype=np.float32, shape=(nb_record, 4**ksize), mode='r+')
    #dump(frequencies, freq_name)
    #frequencies = load(freq_name, mmap_mode = "r+") 
    
    fd = delayed(compute_frequency_memmap)
    Parallel(n_jobs=nbthread, verbose=0)(fd(frequencies, i, str(record.seq), pattern, strand) 
                                         for i, record in enumerate(SeqIO.parse(genome, "fasta")))
    
    return frequencies, freq_name

def join_freq_results(folder, nb_record, ksize):
    """ join results of frequency computations
    """
    freq_name = os.path.join(folder, "frequencies_results")
    freq_folder = os.path.join(folder, "frequencies")
    with h5py.File(freq_name, "w") as outhf:
        frequencies = outhf.create_dataset("frequencies", (nb_record, 4**ksize), dtype="float32")
        for res in os.listdir(freq_folder):
            tmp = res.split("_")
            start, stop = int(tmp[1]), int(tmp[2])
            with h5py.File(os.path.join(freq_folder, res), "r") as inhf:
                freq = inhf.get("frequencies")
                frequencies[start: stop] = freq.value[:]
    return freq_name

def compute_frequencies_joblib_h5py(genome, pattern, strand, nbthread, workdir):
    """ compute frequencies
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    pattern: string
        the binary space pattern to extract spaced-words example: 1001010001 
        ksize is inferred from the number of '1' in the pattern
    strand: string
        select genome strand ('both', 'plus', 'minus')
    workdir: strng
        temporary working directory to store pickled chunk
        
    Return:
    -------
    frequencies: numpy.array
        the samples x features matrix storing NT composition of fasta sequences
    """
    folder = tempfile.mkdtemp(dir=workdir)
    freq_name_folder = os.path.join(folder, 'frequencies')
    if not os.path.isdir(freq_name_folder):
        os.makedirs(freq_name_folder)
    
    # Pre-allocate a writeable shared memory map as a container for the frequencies
    nb_record = sum(1 for record in SeqIO.parse(genome, "fasta"))
    
    
    chunksize = nb_record // (nbthread*SCALING) # some scaling value
    if nb_record % (nbthread*SCALING) != 0:
        chunksize += 1
        
    fd = delayed(compute_frequency_h5py_chunk)
    Parallel(n_jobs=nbthread, verbose=0)(fd(freq_name_folder, seqchunk, pattern, strand, start, stop)
                                         for start, stop, seqchunk in read_seq_chunk_pos(genome, chunksize))
                                         #str(record.seq), pattern, strand) 
                                         #for i, record in enumerate(SeqIO.parse(genome, "fasta")))
    
    # now join the results
    pattern=str(pattern)
    ksize = pattern.count("1")
    freq_name = join_freq_results(folder, nb_record, ksize)
    #print(folder)
    return None, freq_name


def compute_frequencies(mthdrun, large, genome, pattern, strand, 
                        distchunksize, threads_max, workdir):
    """ choose which function to call to compute frequencies
    """
    freq_name = None
    if mthdrun == "scoop":
        frequencies = compute_frequencies_scoop(genome, pattern, strand, distchunksize)
    elif mthdrun == "joblib":
        if large == "memmap":
            frequencies, freq_name = compute_frequencies_joblib_memmap(genome, pattern, strand, threads_max, workdir)
        elif large == "h5py":
            frequencies, freq_name = compute_frequencies_joblib_h5py(genome, pattern, strand, threads_max, workdir)
        else:
            frequencies = compute_frequencies_joblib(genome, pattern, strand, threads_max) 
    else:
        print("Method {} is unknown".format(mthdrun), file=sys.stderr)
            
    return frequencies, freq_name
    
    
def get_cmd():
    """ get command line argument
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--assembly", action="store", required=True, dest="genome", 
                        help="multifasta of the genome assembly")
    parser.add_argument("-k", "--lgMot", action="store", dest="pattern", default=4, type=int,
                        help="word lenght / kmer length / k [default:%(default)d]")
    #parser.add_argument("-k", "--lgMot", action="store", dest="k", default="1111", 
                        #help="word lenght / kmer length / k [default:%(default)d]")
    parser.add_argument("-s", "--strand", action="store", dest="strand", default="both", choices=["both", "plus", "minus"],
                        help="strand used to compute microcomposition. [default:%(default)s]")
    parser.add_argument("-d", "--distance", action="store", dest="dist", default="Eucl", choices=["Eucl", "JSD", "KT"], 
                        help="how to compute distance between two signatures : Eucl : Euclidean[default:%(default)s], JSD : Jensen-Shannon divergence, KT: Kendall's tau")
    parser.add_argument("--freq-chunk-size", action="store", dest="freqchunksize", type=int, default=250,
                        help="the size of the chunk to use in scoop to compute frequencies")
    parser.add_argument("--dist-chunk-size", action="store", dest="distchunksize", type=int, default=250,
                        help="the size of the chunk to use in scoop to compute distances")
    parser.add_argument("--method", action="store", choices=["scoop", "joblib"], default= "joblib",dest="mthdrun", help="don't use scoop to compute distances use joblib", required=True)
    parser.add_argument("--large", action="store", dest="large", choices=["None", "memmap", "h5py"], help="used in combination with joblib for large dataset", default="None")
    parser.add_argument("-c", "--cpu", action="store", dest="threads_max", type=int, default=4, 
                        help="how many threads to use for windows microcomposition computation[default:%(default)d]")
    parser.add_argument("-o", "--out", action="store", dest="out_file", default="phyloligo.out", 
                        help="output file[default:%(default)s]")
    parser.add_argument("-q", "--outfreq", action="store", dest="out_freq_file", default="phyloligo.freq", 
                        help="kmer frequencies output file [default:infile_%(default)s]")
    parser.add_argument("-w", "--workdir", action="store", dest="workdir", default=".", help="working directory")
    parser.add_argument("-p", "--pattern", action="store", dest="pattern", default="1111", help="spaced-word pattern string, only containing 1s and 0s, i.e. '100101001', default='1111'")
    
    
    params = parser.parse_args()

    params.workdir = os.path.abspath(params.workdir)
        
    return params

def main():
    params = get_cmd()
    
    # quick hack to turn the k-mer parameter into a pattern (whithout joker)
    if type(params.pattern) == int:
        params.pattern = str("1") * params.pattern
    
    print("Using pattern {}".format(params.pattern))
    if not os.path.isdir(params.workdir):
        os.makedirs(params.workdir)
    
    # compute word frequency of each sequence
    print("Computing frequencies")
    frequencies, freq_name = compute_frequencies(params.mthdrun, params.large,
                                            params.genome, params.pattern, params.strand, 
                                            params.distchunksize, params.threads_max, params.workdir)
        
    # compute pairwise distances
    print("Computing Pairwise distances")
    res = compute_distances(params.mthdrun, params.large, frequencies, freq_name, params.out_file, 
                                      params.dist, params.threads_max, params.freqchunksize, params.workdir)
        
    # save result in a numpy matrix
    if (params.out_freq_file):
        print("Writing frequency matrix")
        np.savetxt(params.out_freq_file, frequencies, delimiter="\t")
        
    # save result in a numpy matrix
    if not (params.mthdrun == "joblib" and params.large != "None"):
        print("Writing distance matrix")
        np.savetxt(params.out_file, res, delimiter="\t")
        
    return 0
       
if __name__ == "__main__":
    #start = time.time()
    ret = main()      
    #stop = time.time()
    #print("Exec time: {}".format(stop-start))
    sys.exit(0)

