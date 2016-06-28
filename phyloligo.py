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
  
from optparse import OptionParser
import argparse

from Bio import SeqIO
from Bio.Seq import Seq

import os, sys, re, math
import numpy
import multiprocessing
from collections import Counter
from itertools import product

numpy.seterr(divide='ignore', invalid='ignore')

def KL(a,b):
    """ compute the KL distance
    """
    #with numpy.errstate(invalid='ignore'):
    d = a * numpy.log(a/b)
    d[numpy.isnan(d)]=0 
    d[numpy.isinf(d)]=0
    return (numpy.sum(d))*10000

def Eucl(a,b):
    """ compute Euclidean distance
    """
    #with numpy.errstate(invalid='ignore'):
    d = pow(a-b,2)
    d[numpy.isnan(d)]=0
    d[numpy.isinf(d)]=0
    return numpy.sqrt(numpy.sum(d))*10000

def JSD(a,b):
    """ Compute JSP distance
    """
    #with numpy.errstate(invalid='ignore'):
    h = (a + b)/2
    d = (KL(a,h)/2)+(KL(b,h)/2)
    return d

def vector_to_matrix(profile):
    """ transform a vector of profile to a matrix
    """
    return list((zip(*(iter(profile),)*int(math.sqrt(len(profile))))))

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
        a Biopython seq object
    strand: string
        select wich strand to use
    
    Return:
    -------
    seq: string
        the sequence strand
    """
    Bioseq_record=Seq(seq)
    if (strand == "both"):
        return str(str(seq)+str(Bioseq_record.reverse_complement())).upper()
    elif (strand == "minus"):
        return str(Bioseq_record.reverse_complement()).upper()
    elif (strand == "plus"):
        return str(seq).upper()
    else:
        print("Error, strand parameter of selectd_strand() should be choose from {'both', 'minus', 'plus'}", file=sys.stderr)
        sys.exit(1)


def frequency (seq, ksize=4, strand="both"):
    """ compute kmer frequency
    """
    seq=select_strand(seq,strand)
    seq_words=list()
    d=dict()
    for s in re.split('[^ACGTacgt]+',seq): #excludes what is not a known characterised nucleotide 
        seq_letters=list(s)
        #print(len(s))
        if (len(seq_letters) >=ksize):
            # launch k-1 times the word generation with 1 nucleotide shift every iteration to have overlapping words.
            for i in range(ksize-1):
                #generate words from seq_letters
                words=list((zip(*(iter(seq_letters),)*ksize)))
                seq_words.extend(list(map(''.join,words))) # adds the words for this subsequence and frame to the total list of words
                seq_letters.pop(0) # shift one to compute overlapping words at the next iteration
    c = Counter(seq_words)
    word_count=sum(c.values())
    #print(word_count)
    if(word_count > 0):
        ret=list()
        word_universe=list(map("".join,(list(product(("C","G","A","T"),("C","G","A","T"),("C","G","A","T"),("C","G","A","T"))))))
        for w in word_universe:
            if(c.get(w) != None):
                ret.append(c.get(w)/word_count)
            else:
                ret.append(0)
        return ret  
    return 1


def pairwise_distance (args):
    """ compute pairwise distance
    """
    res=[]
    for i,j,dist,ksize,strand in args:
        res.append((i,j,globals()[dist](numpy.array(frequency(str(records[i].seq))),numpy.array(frequency(str(records[j].seq))))))
    return res


def parallel_distance(genome, nb_thread, dist, ksize, strand):
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
            args.append((i, j, dist, ksize, strand))
    #print(args)
    parallel_args_set = chunkitize(args, nb_thread) 
    #print(parallel_args_set)
    pool = multiprocessing.Pool(processes=nb_thread)
    res = pool.map(pairwise_distance, parallel_args_set)
    pool.close()
    pool.join()
    symmetrical_distance_matrix = numpy.zeros((len(records),len(records)))
    for i in res:
        for g in i:
            dist_pairs.append(g)
            symmetrical_distance_matrix[(g[0],g[1])]=g[2]
            symmetrical_distance_matrix[(g[1],g[0])]=g[2]
    #print(dist_pairs)
    #print(symmetrical_distance_matrix)
    
    return symmetrical_distance_matrix


def get_cmd():
    """ get command line argument
    """
    #Utilisation = "%prog [-i FILE] [options]"
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--assembly", action="store", required=True, dest="genome", help="multifasta of the genome assembly")
    #parser.add_argument("-c", "--conta", action="store", dest="conta", help="multifasta of the contaminant species training set")
    #parser.add_argument("-r", "--host", action="store", dest="host", help="optional multifasta of the host species training set")
    #parser.add_argument("-n", "--n_max_freq_in_windows", action="store", type=float, dest="n_max_freq_in_windows", default=0.4, help="maximum proportion of N tolerated in a window to compute the microcomposition anyway [0~1]. Too much 'N's will reduce the number of kmer counts and will artificially coarse the resolution of the frequencies. If your assembly contains many stretches of 'N's, consider rising this parameter and shortening the windows step in order to allow computation and signal in the output, it might cost you computational time though. Windows which do not meet this criteria will be affected a distance of 'nan'")
    parser.add_argument("-k", "--lgMot", action="store", dest="k", type="int", default=4, help="word wise / kmer lenght / k [default:%default]")
    #parser.add_argument("-w", "--windows_size", action="store", dest="windows_size", type=int, help="Sliding windows size (bp)[default:%default]")
    #parser.add_argument("-t", "--windows_step", action="store", dest="windows_step", type=int, help="Sliding windows step size(bp)[default:%default]")
    parser.add_argument("-s", "--strand", action="store", dest="strand", default="both", choices=["both", "plus", "minus"], help="strand used to compute microcomposition. leading, lagging ou both [default:%default]")
    parser.add_argument("-d", "--distance", action="store", dest="dist", default="JSD", choices=["KL", "Eucl", "JSD"], help="how to compute distance between two signatures : KL: Kullback-Leibler, Eucl : Euclidienne[default:%default], JSD : Jensen-Shannon divergence")
    parser.add_argument("-u", "--cpu", action="store", dest="threads_max", type=int, default=4, help="how many threads to use for windows microcomposition computation[default:%default]")
    parser.add_argument("-o", "--out", action="store", dest="out_file", default="phyloligo.out", help="output file[default:%default]")

    params = parser.parse_args()
    if not 0 <= params.n_max_freq_in_windows <= 1.0:
        print("Errorm parameter '-n', '-n_max_freq_in_windows' should be between 0 and 1", file=sys.stderr)
        sys.exit(1)
    return params

def main():
    params = get_cmd()

    #print("A genome was provided")
    res = parallel_distance(params.genome, params.thread_max, params.dist, params.k, params.strand)
    print(res)
    numpy.savetxt(params.out_file, res, delimiter="\t")
    sys.exit(0)
       
if __name__ == "__main__":
    main()      
      

