#!/usr/bin/env python3
""" 
"""
  
from scoop import futures
  
import os, sys, re, pickle
import time
from Bio.Seq import Seq

from collections import Counter
from itertools import product

## clustering
import numpy as np

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


### IN GLOBAL SCOOP MEMORY ###

pathin = sys.argv[1]
pathout = sys.argv[2]

# global scoop memory
with open(pathin, "rb") as inf:
   seqchunk = pickle.load(inf)
   
   
   
if __name__ == "__main__":
    chunkfreq = list(futures.map(frequency_pack, seqchunk))
    with open(pathout, "wb") as outf:
        pickle.dump(chunkfreq, outf)
        
    sys.exit(0)
    