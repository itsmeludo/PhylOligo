#!/usr/bin/env python3

from scoop import futures

import os, sys, re
import tempfile
import shlex, subprocess, pickle
import h5py
from .phyloutils import select_strand, get_nb_records

from Bio import SeqIO

from collections import Counter
from itertools import product

# parallelism
import numpy as np
from sklearn.externals.joblib import Parallel, delayed
from sklearn.externals.joblib import dump, load

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

def compute_frequency_memmap(frequency, i, seq, ksize=4, strand="both"):
    """ function used in joblib to write in the memmap array frequency
    
    Parameters:
    -----------
    frequency: np.memmap
        a memmap array to store computed frequencies
    i: int 
        index of the sequence
    seq: string
        The nucleotide sequence
    ksize: int
        The size of the kmer
    strand: string
        which strand to used
        
    """
    seq = select_strand(seq, strand)
    # we work on upper case letters
    seq = seq.upper()
    # raw word count
    count_words, kword_count = cut_sequence_and_count(seq, ksize)
    # create feature vector
    frequency[i] = count2freq(count_words, kword_count, ksize)
    
def compute_frequency_h5py(freq_name_folder, i, seq, ksize=4, strand="both"):
    """ function used in joblib to write in the h5py array frequency
    
    Parameters:
    -----------
    freq_name: string
        path to the folder file holding h5py frequencies
    i: int 
        index of the sequence
    seq: string
        The nucleotide sequence
    ksize: int
        The size of the kmer
    strand: string
        which strand to used
        
    """
    seq = select_strand(seq, strand)
    # we work on upper case letters
    seq = seq.upper()
    # raw word count
    count_words, kword_count = cut_sequence_and_count(seq, ksize)
    # create feature vector
    freq_name = os.path.join(freq_name_folder, "frequencies_{}".format(i))
    with h5py.File(freq_name, "w") as hf:
        frequency = hf.create_dataset("frequencies", (4**ksize,), dtype="float32")
        frequency[...] = count2freq(count_words, kword_count, ksize)[:]
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

def compute_frequencies_scoop(genome, ksize, strand, chunksize):
    """ compute frequencies
    
    Parameters:
    -----------
    genome: string
        path to the genome file
    ksize: int
        kmer size to use
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
    for seqchunk in read_seq_chunk(genome, chunksize, ksize, strand):
        chunkfreq = futures.map(frequency_pack, seqchunk)
        frequencies.extend(chunkfreq)
                           
    frequencies = np.array(frequencies)
    return frequencies

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
    pathin = tempfile.mktemp(dir=workdir)
    pathout = tempfile.mktemp(dir=workdir)
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
        with open(pathout, "rb") as inf:
            chunkfreq = pickle.load(inf)
        frequencies.extend(chunkfreq)

    os.remove(pathin)
    os.remove(pathout)
    frequencies = np.array(frequencies)

    return frequencies

def compute_frequencies_joblib(genome, ksize, strand, nbthread):
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
    frequencies = np.vstack(frequencies)
    #print("Frequency shape: {}".format(frequencies.shape))
    if len(frequencies.shape) < 2:
        frequencies.reshape(-1, 1)
    return frequencies

def compute_frequencies_joblib_memmap(genome, ksize, strand, nbthread):
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
    folder = tempfile.mkdtemp()
    freq_name = os.path.join(folder, 'frequencies')
    
    
    # Pre-allocate a writeable shared memory map as a container for the frequencies
    nb_record = get_nb_records(genome)
    frequencies = np.memmap(freq_name, dtype=np.float32, shape=(nb_record, 4**ksize), mode='w+')
    #print(freq_name)
    del frequencies
    frequencies = np.memmap(freq_name, dtype=np.float32, shape=(nb_record, 4**ksize), mode='r+')
    #dump(frequencies, freq_name)
    #frequencies = load(freq_name, mmap_mode = "r+") 
    
    fd = delayed(compute_frequency_memmap)
    Parallel(n_jobs=nbthread, verbose=0)(fd(frequencies, i, str(record.seq), ksize, strand) 
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
            idx = int(res.split("_")[1])
            with h5py.File(os.path.join(freq_folder, res), "r") as inhf:
                freq = inhf.get("frequencies")
                frequencies[idx] = freq.value[:]
    return freq_name
    


def compute_frequencies_joblib_h5py(genome, ksize, strand, nbthread):
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
    folder = tempfile.mkdtemp()
    freq_name_folder = os.path.join(folder, 'frequencies')
    if not os.path.isdir(freq_name_folder):
        os.makedirs(freq_name_folder)
    
    # Pre-allocate a writeable shared memory map as a container for the frequencies
    nb_record = sum(1 for record in SeqIO.parse(genome, "fasta"))
    
    #with h5py.File(freq_name, "w") as hf:
        #frequencies = hf.create_dataset("frequencies", (nb_record, 4**ksize), dtype="float32")
    
    fd = delayed(compute_frequency_h5py)
    Parallel(n_jobs=nbthread, verbose=0)(fd(freq_name_folder, i, str(record.seq), ksize, strand) 
                                         for i, record in enumerate(SeqIO.parse(genome, "fasta")))
    #for i, record in enumerate(SeqIO.parse(genome, "fasta")):
        #compute_frequency_h5py(freq_name_folder, i, str(record.seq), ksize, strand)
    
    # now join the results
    freq_name = join_freq_results(folder, nb_record, ksize)
    #print(folder)
    return None, freq_name


def compute_frequencies(mthdrun, large, genome, k, strand, 
                        distchunksize, threads_max, workdir):
    """ choose which function to call to compute frequencies
    """
    freq_name = None
    if mthdrun == "scoop1":
        frequencies = compute_frequencies_scoop(genome, k, strand, distchunksize)
    elif mthdrun == "scoop2":
        frequencies = compute_frequencies_pickle(genome, k, strand, distchunksize, threads_max, workdir)
    elif mthdrun == "joblib":
        if large == "memmap":
            frequencies, freq_name = compute_frequencies_joblib_memmap(genome, k, strand, threads_max)
        elif large == "h5py":
            frequencies, freq_name = compute_frequencies_joblib_h5py(genome, k, strand, threads_max)
        else:
            frequencies = compute_frequencies_joblib(genome, k, strand, threads_max) 
    else:
        print("Method {} is unknown".format(mthdrun), file=sys.stderr)
            
    return frequencies, freq_name
    