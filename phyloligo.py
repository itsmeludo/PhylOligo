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
  
import os, sys, re, math, argparse
import time
import logging

from Bio import SeqIO
from Bio.Seq import Seq

import numpy as np
import multiprocessing
from collections import Counter
from itertools import product

## clustering
from sklearn.metrics.pairwise import pairwise_distances
#import hdbscan

np.seterr(divide='ignore', invalid='ignore')

def KL(a,b):
    """ compute the KL distance
    """
    d = a * np.log(a/b)
    d[np.isnan(d)]=0 
    d[np.isinf(d)]=0
    return (np.sum(d))*10000

def Eucl(a,b):
    """ compute Euclidean distance
    """
    d = pow(a-b,2)
    d[np.isnan(d)]=0
    d[np.isinf(d)]=0
    return np.sqrt(np.sum(d))*10000

def JSD(a,b):
    """ Compute JSD distance
    """
    h = (a + b)/2
    d = (KL(a,h)/2)+(KL(b,h)/2)
    return d

def vector_to_matrix(profile):
    """ transform a vector of profile to a matrix
    """
    return list((zip(*(iter(profile),)*int(math.sqrt(len(profile))))))

dist_methods = {
    "KL": KL,
    "Eucl": Eucl,
    "JSD": JSD,
    }

def chunkitize(toseparate, nb_chunks):
    """ separate a list into sublist of the same size
    
    Parameters:
    -----------
    toseparate: list
        to be cut into sublists
    nb_chunks: int
        the number of sublists
    
    Return:
    -------
    out: list
        a list of n=nb_chunks sublist
    
    >>> chunkitize([1, 2, 3, 4, 5], 2)
    ... [[1, 2], [3, 4, 5]]
    """
    chunk_size= int((len(toseparate) / float(nb_chunks)))
    out = [toseparate[chunk_size*i:(chunk_size*i+chunk_size) if (i!=nb_chunks-1) else len(toseparate)] for i in range(0, nb_chunks)]
    return out


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

def cut_sequence(seq, ksize):
    """ cut sequence in words of size k
    
    Parameters:
    -----------
    seq: string
        The nucleotide sequence
    ksize: int
        The size of the kmer
    
    Return:
    -------
    count_words: dict
        a dictionary storing the number of observations of each word
    kword_count: int
        the total number of words
    """
    seq_words = list()
    # re.split: excludes what is not a known characterised nucleotide 
    for subseq in re.split('[^ACGTacgt]+', seq): 
        if (len(subseq) >= ksize):
            #get all kword in sub-sequence
            seq_words.extend(subseq[i: i+ksize] for i in range(len(subseq)-(ksize-1)))  
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
        letter_list = list(product(("C","G","A","T"), repeat=ksize))
        #word_universe = list(map("".join,letter_list))
        for letter_word in letter_list:
            kword = "".join(letter_word)
            if kword in count_words:
                features.append(count_words[kword]/kword_count)
            else:
                features.append(0)
    else:
        features = np.array([np.inf for kword in letter_list])
    return features

def frequency(seq, ksize=4, strand="both"):
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
    seq = select_strand(seq, strand)
    # we work on upper case letters
    seq = seq.upper()
    # raw word count
    count_words, kword_count = cut_sequence(seq, ksize)
    # create feature vector
    features = count2freq(count_words, kword_count, ksize)
    return features


#### OLD CODE BELOW #### TODO REMOVE AFTER TESTING ###

def pairwise_distances_old (args):
    """ compute pairwise distance
    """
    res=[]
    for i, j, dist, ksize, strand in args:
        res.append((i, j, dist(frequency(str(records[i].seq), ksize, strand), 
                               frequency(str(records[j].seq), ksize, strand))))
    return res


def parallel_distance_old(genome, nb_thread, dist, ksize, strand):
    """ compute distance on multithread
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    nb_thread: int
        number of parallel instance
    dist: string
        distance method to use ('JSD', 'Eucl', 'KL')
    ksize: int
        kmer size to use
    strand: string
        select genome strand ('both', 'plus', 'minus')
        
    Return:
    -------
    symmetrical_distance_matrix: numpy.array
        the matrix storing pairwize distances
    """
    args=[]
    dist_pairs=[]
    global records
    records = list(SeqIO.parse(genome, "fasta"))
    for i in range(len(records)):
        for j in range(i+1,len(records)):
            args.append((i, j, dist_methods[dist], ksize, strand))
    #print(args)
    #adds a finer granularity for jobs dispatch to allow backfill of cores after the shorter jobs among the "nb_thread" initially requested finished and the longest still run on few cores.
    if (int(len(args)/nb_thread) >= nb_thread*nb_thread*parallel_core_granularity_factor):
        parallel_args_set = chunkitize(args, int(nb_thread*(nb_thread/2)*parallel_core_granularity_factor))
    else:
        # if there is only a few fasta entries, adding granularity will cost time
        parallel_args_set = chunkitize(args, nb_thread) 
    #print(parallel_args_set)
    pool = multiprocessing.Pool(processes=nb_thread)
    res = pool.map(pairwise_distances_old, parallel_args_set)
    pool.close()
    pool.join()
    symmetrical_distance_matrix = np.zeros((len(records),len(records)))
    for i in res:
        for g in i:
            dist_pairs.append(g)
            symmetrical_distance_matrix[(g[0],g[1])]=g[2]
            symmetrical_distance_matrix[(g[1],g[0])]=g[2]
    #print(dist_pairs)
    #print(symmetrical_distance_matrix)
    
    return symmetrical_distance_matrix

#### OLD CODE ABOVE #### TODO REMOVE AFTER TESTING ###

def compute_frequencies(genome, nb_thread, ksize, strand):
    """ compute frequencies
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    nb_thread: int
        number of parallel instance
    ksize: int
        kmer size to use
    strand: string
        select genome strand ('both', 'plus', 'minus')
        
    Return:
    -------
    frequencies: numpy.array
        the samples x features matrix storing NT composition of fasta sequences
    """
    
    # compute frequencies # TODO parallelization of frequencies computation
    frequencies = list()
    for record in  SeqIO.parse(genome, "fasta"):
        frequencies.append(frequency(str(record.seq), ksize, strand))
    
    return np.array(frequencies)

def compute_distances(frequencies, metric="Eucl", n_jobs=1):
    """ compute pairwises distances
    
    Parameters
    ----------
    frequencies: np.array
        a list of frequencies, each row corresponding to a sample each column to a "word"
    metric: string
        distance method to use ('JSD', 'Eucl', 'KL')
    n_jobs: int
        number of parallel job to execute
    
    Return:
    -------
    distances: np.array
        (n_samples, n_samples) distance matrix
    """
    if metric == "Eucl":
        # use euclidean distance of sklearn
        distances = pairwise_distances(frequencies, metric="euclidean", n_jobs=n_jobs)
        distances[np.isnan(distances)]=0
        distances[np.isinf(distances)]=0
        distances *= 10000
    elif metric == "EuclLoc":
        distances = pairwise_distances(frequencies, metric=Eucl, n_jobs=n_jobs)
    elif metric == "JSD":
        distances = pairwise_distances(frequencies, metric=JSD, n_jobs=n_jobs)
    elif metric == "KL":        
        distances = pairwise_distances(frequencies, metric=KL, n_jobs=n_jobs)
    else:
        print("Error, unknown method {}".format(metric), file=sys.stderr)
        sys.exit(1)
    return distances

def get_cmd():
    """ get command line argument
    """
    #Utilisation = "%prog [-i FILE] [options]"
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--assembly", action="store", required=True, dest="genome", 
                        help="multifasta of the genome assembly")
    #parser.add_argument("-c", "--conta", action="store", dest="conta", help="multifasta of the contaminant species training set")
    #parser.add_argument("-r", "--host", action="store", dest="host", help="optional multifasta of the host species training set")
    #parser.add_argument("-n", "--n_max_freq_in_windows", action="store", type=float, dest="n_max_freq_in_windows", default=0.4, help="maximum proportion of N tolerated in a window to compute the microcomposition anyway [0~1]. Too much 'N's will reduce the number of kmer counts and will artificially coarse the resolution of the frequencies. If your assembly contains many stretches of 'N's, consider rising this parameter and shortening the windows step in order to allow computation and signal in the output, it might cost you computational time though. Windows which do not meet this criteria will be affected a distance of 'nan'")
    parser.add_argument("-k", "--lgMot", action="store", dest="k", type=int, default=4, 
                        help="word lenght / kmer length / k [default:%(default)d]")
    #parser.add_argument("-w", "--windows_size", action="store", dest="windows_size", type=int, help="Sliding windows size (bp)[default:%default]")
    #parser.add_argument("-t", "--windows_step", action="store", dest="windows_step", type=int, help="Sliding windows step size(bp)[default:%default]")
    parser.add_argument("-s", "--strand", action="store", dest="strand", default="both", choices=["both", "plus", "minus"],
                        help="strand used to compute microcomposition. [default:%(default)s]")
    parser.add_argument("-d", "--distance", action="store", dest="dist", default="Eucl", choices=["KL", "Eucl", "JSD", "EuclLoc"], 
                        help="how to compute distance between two signatures : KL: Kullback-Leibler, Eucl : Euclidean[default:%(default)s], JSD : Jensen-Shannon divergence")
    parser.add_argument("-u", "--cpu", action="store", dest="threads_max", type=int, default=4, 
                        help="how many threads to use for windows microcomposition computation[default:%(default)d]")
    #parser.add_argument("-g", "--granularity", action="store", dest="parallel_core_granularity_factor", type=float, default=1, 
                        #help="Factor to refine the granularity of parallel threads. The higher the factor, the greater the number of smaller bins.[default:%(default)d]")
    #parser.add_argument("-a", "--sampling", action="store", dest="sampling", type=float, default=0,
                        #help="performs first a read sampling on the specified percentage of read, can be used for quick analyses, "
                        #"tests, and large datasets [default:%(default)f, input value must be between 0 and 100]")
    parser.add_argument("-o", "--out", action="store", dest="out_file", default="phyloligo.out", 
                        help="output file[default:%(default)s]")

    params = parser.parse_args()

    # check sampling parameters
    #if not 0<=params.sampling<=100:
        #print("Error, sampling parameters must be between 0 and 100", file=sys.stderr)
        #sys.exit(1)
        
    return params

def main():
    params = get_cmd()
    #### OLD 
    #global parallel_core_granularity_factor
    #parallel_core_granularity_factor=1
    #res = parallel_distance_old(params.genome, params.threads_max, params.dist, params.k, params.strand)
    ####
    frequencies = compute_frequencies(params.genome, params.threads_max, params.k, params.strand)
    res = compute_distances(frequencies, metric=params.dist, n_jobs=params.threads_max)
    # save result in a numpy matrix
    np.savez_compressed(params.out_file, res, delimiter="\t")
    
    return 0
       
if __name__ == "__main__":
    start = time.time()
    ret = main()      
    stop = time.time()
    print("Exec time: {}".format(stop-start))
    sys.exit(0)
    
      

