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
  
import os, sys, re, argparse, pickle, shlex
import uuid
import time
import logging

from Bio import SeqIO
from Bio.Seq import Seq

import multiprocessing, subprocess
from collections import Counter
from itertools import product

## clustering
import numpy as np
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.externals.joblib import Parallel, delayed
from sklearn.utils import gen_even_slices

#import hdbscan
#from scipy.stats import entropy
#from numpy.linalg import norm

np.seterr(divide='ignore', invalid='ignore')

#### DISTANCES ####

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

def make_freqchunk(frequencies, chunksize):
    """ prepare frequencies to be parallelized
    """
    chunk = list()
    for i in range(len(frequencies)):
        freqi = frequencies[i]
        for j in range(i, len(frequencies)):
            freqj = frequencies[j]
            chunk.append((i, j, freqi, freqj))
            if len(chunk) == chunksize:
                yield chunk
                chunk = list()
    # the last chunk 
    if chunk != []:
        yield chunk

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
    
def compute_distances_pickle(frequencies, chunksize, metric="Eucl", n_jobs=1, workdir="/tmp/"):
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
    #scoop.logger.info("Starting distance computation")
    if metric == "Eucl":
        pathin, pathout = os.path.join(workdir, str(uuid.uuid4())), os.path.join(workdir, str(uuid.uuid4()))
        distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
        for freqchunk in make_freqchunk(frequencies, chunksize):
            with open(pathin, "wb") as outf:
                pickle.dump(freqchunk, outf)
            # subprocess the computation
            cmd = "python3 -m scoop -n {} phylo_batchdist.py {} {} {}".format(n_jobs, pathin, pathout, metric)
            cmd = shlex.split(cmd)
            try:
                res = subprocess.check_call(cmd)
            except:
                print("Error running phylo_batchdist.py on {} {} /{}".format(pathin, pathout, metric), file=sys.stderr)
                sys.exit(1)
                
            with open(pathout, "rb") as inf:
                res = pickle.load(inf)
            #res = futures.map(compute_Eucl_unpack, freqchunk)
            for i, j, d in res:
                distances[i, j] = distances[j, i] = d
        os.remove(pathin)
        os.remove(pathout)
    elif metric == "JSD":
        pathin, pathout = os.path.join(workdir, str(uuid.uuid4())), os.path.join(workdir, str(uuid.uuid4()))

        distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
        for freqchunk in make_freqchunk(frequencies, chunksize):
            with open(pathin, "wb") as outf:
                pickle.dump(freqchunk, outf)
            # subprocess the computation
            cmd = "python3 -m scoop -n {} phylo_batchdist.py {} {} {}".format(n_jobs, pathin, pathout, metric)
            cmd = shlex.split(cmd)
            try:
                res = subprocess.check_call(cmd)
            except:
                print("Error running phylo_batchdist.py on {} {} /{}".format(pathin, pathout, metric), file=sys.stderr)
                sys.exit(1)    
            with open(pathout, "rb") as inf:
                res = pickle.load(inf)
            #res = futures.map(compute_JSD_unpack, freqchunk)
            for i, j, d in res:
                distances[i, j] = distances[j, i] = d
        os.remove(pathin)
        os.remove(pathout)
    else:
        print("Error, unknown method {}".format(metric), file=sys.stderr)
        sys.exit(1)
    return distances

def compute_distances(frequencies, chunksize, metric="Eucl", localrun=False, n_jobs=1):
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
    #scoop.logger.info("Starting distance computation")
    if metric == "Eucl":
        if localrun:
            distances = pairwise_distances(frequencies, metric=Eucl, n_jobs=n_jobs)
        else:
            distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
            for freqchunk in make_freqchunk(frequencies, chunksize):
                res = futures.map(compute_Eucl_unpack, freqchunk)
                for i, j, d in res:
                    distances[i, j] = distances[j, i] = d
    elif metric == "JSD":
        if localrun:
            distances = pairwise_distances(frequencies, metric=JSD, n_jobs=n_jobs)
        else:
            distances = np.zeros((len(frequencies), len(frequencies)), dtype=float)
            for freqchunk in make_freqchunk(frequencies, chunksize):
                res = futures.map(compute_JSD_unpack, freqchunk)
                for i, j, d in res:
                    distances[i, j] = distances[j, i] = d
    
    else:
        print("Error, unknown method {}".format(metric), file=sys.stderr)
        sys.exit(1)
    return distances

#### frequencies computation ####

def cut_sequence_and_count(seq, ksize):
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
    for subseq in re.split('[^ACGT]+', seq): 
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

def compute_frequency(seq, ksize=4, strand="both"):
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
    count_words, kword_count = cut_sequence_and_count(seq, ksize)
    # create feature vector
    features = count2freq(count_words, kword_count, ksize)
    return features

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

def compute_frequencies(genome, ksize, strand, chunksize):
    """ compute frequencies
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    ksize: int
        kmer size to use
    strand: string
        select genome strand ('both', 'plus', 'minus')
        
    Return:
    -------
    frequencies: numpy.array
        the samples x features matrix storing NT composition of fasta sequences
    """
    #scoop.logger.info("Starting frequencies computation")
    # compute frequencies # TODO parallelization of frequencies computation
    frequencies = list()
    for seqchunk in read_seq_chunk(genome, chunksize, ksize, strand):
        chunkfreq = futures.map(frequency_pack, seqchunk)
        frequencies.extend(chunkfreq)
                           
    return np.array(frequencies)

def compute_frequencies_pickle(genome, ksize, strand, chunksize, nbthread, workdir):
    """ compute frequencies
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    ksize: int
        kmer size to use
    strand: string
        select genome strand ('both', 'plus', 'minus')
    workdir: strng
        temporary working directory to store pickled chunk
        
    Return:
    -------
    frequencies: numpy.array
        the samples x features matrix storing NT composition of fasta sequences
    """
    #scoop.logger.info("Starting frequencies computation")
    # compute frequencies # TODO parallelization of frequencies computation
    frequencies = list()
    pathin, pathout = os.path.join(workdir, str(uuid.uuid4())), os.path.join(workdir, str(uuid.uuid4()))
    for seqchunk in read_seq_chunk(genome, chunksize, ksize, strand):
        with open(pathin, "wb") as outf:
            pickle.dump(seqchunk, outf)
        
        cmd = "python3 -m scoop -n {} phylo_batchfreq.py {} {}".format(nbthread, pathin, pathout)
        cmd = shlex.split(cmd)
        try:
            ret = subprocess.check_call(cmd)
        except:
            print("Error running phylo_batchfreq.py on {} {}".format(pathin, pathout))
            sys.exit(1)
        #chunkfreq = futures.map(frequency_pack, seqchunk)
        with open(pathout, "rb") as inf:
            chunkfreq = pickle.load(inf)
        frequencies.extend(chunkfreq)
    os.remove(pathin)
    os.remove(pathout)
    return np.array(frequencies)

def compute_frequencies_joblib(genome, ksize, strand, chunksize, nbthread):
    """ compute frequencies
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    ksize: int
        kmer size to use
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
        fd(str(record.seq), ksize, strand) for record in SeqIO.parse(genome, "fasta"))
    
    #for ar in frequencies:
        #print(ar.shape)
        
    return np.vstack(frequencies)

#### Sequences ####

def read_seq_chunk(genome, chunksize, ksize, strand):
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
    for record in SeqIO.parse(genome, "fasta"):
        seqchunk.append((str(record.seq), ksize, strand))
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
    parser.add_argument("--method", action="store", choices=["scoop1", "scoop2", "joblib"], dest="mthdrun", help="don't use scoop to compute distances use joblib", default=False)
    
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
    if params.mthdrun == "scoop1":
        frequencies = compute_frequencies(params.genome, params.k, params.strand, params.distchunksize)
    elif params.mthdrun == "scoop2":
        frequencies = compute_frequencies_pickle(params.genome, params.k, params.strand, params.distchunksize, params.threads_max, params.workdir)
    elif params.mthdrun == "joblib":
        frequencies = compute_frequencies_joblib(params.genome, params.k, params.strand, params.distchunksize, params.threads_max)
        
    # compute pairwise distances
    if params.mthdrun == "joblib":
        res = compute_distances(frequencies, params.freqchunksize, metric=params.dist, localrun=True, n_jobs=params.threads_max)
    elif params.mthdrun == "scoop1":
        res = compute_distances(frequencies, params.freqchunksize, metric=params.dist, localrun=False, n_jobs=params.threads_max)
    elif params.mthdrun == "scoop2":
        res = compute_distances_pickle(frequencies, params.freqchunksize, metric=params.dist, n_jobs=params.threads_max, workdir=params.workdir)
    # save result in a numpy matrix
    np.savetxt(params.out_file, res, delimiter="\t")
    return 0
       
if __name__ == "__main__":
    start = time.time()
    ret = main()      
    stop = time.time()
    print("Exec time: {}".format(stop-start))
    sys.exit(0)
    
      

