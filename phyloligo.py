#!/usr/bin/env python3

### Author: Ludovic V. Mallet, PhD
### 2016.03.22
### licence: GPLv3
### Version: Alpha.1
### Garanteed with misatkes. <- Including this one.


# This program builds a phylotree based on oligonucleotide profiles (microcomposition) distances for sequences given as arguments.
# See the help by calling the program without any argument.


#dependencies: 
  #Biopython
  #numpy
  #cython
  
from optparse import OptionParser

from Bio import SeqIO
from Bio.Seq import Seq

import os
import re
import math
import numpy
import multiprocessing
from collections import Counter
from itertools import product


numpy.seterr(divide='ignore', invalid='ignore')

def KL(a,b):
  #with numpy.errstate(invalid='ignore'):
  d = a * numpy.log(a/b)
  d[numpy.isnan(d)]=0 
  d[numpy.isinf(d)]=0
  return (numpy.sum(d))*10000

def Eucl(a,b):
  #with numpy.errstate(invalid='ignore'):
  d = pow(a-b,2)
  d[numpy.isnan(d)]=0
  d[numpy.isinf(d)]=0
  return numpy.sqrt(numpy.sum(d))*10000

def JSD(a,b):
  #with numpy.errstate(invalid='ignore'):
  h = (a + b)/2
  d = (KL(a,h)/2)+(KL(b,h)/2)
  return d



def vector_to_matrix(profile):
  return list((zip(*(iter(profile),)*int(math.sqrt(len(profile))))))




def chunkitize(liste, chunks):
  out=list()
  chunk_size= int((len(liste) / float(chunks)))
  for i in range(0,chunks):
    out.append(liste[chunk_size*i:(chunk_size*i+chunk_size)if(i!=chunks-1)else len(liste)])
  return out

#chunkitize([1, 2, 3, 4, 5], 2)




def select_strand (seq,strand="both"):
  Bioseq_record=Seq(seq)
  if(strand == "both"):
    return str(str(seq)+str(Bioseq_record.reverse_complement())).upper()
  elif(strand == "minus"):
    return str(Bioseq_record.reverse_complement()).upper()
  elif(strand == "plus"):
    return str(seq).upper()
  else:
    return 1


def frequency (seq,ksize=4,strand="both"):
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
  else:
    return 1


def pairwise_distance (args):
    res=[]
    for i,j,dist,ksize,strand in args:
        res.append((i,j,globals()[dist](numpy.array(frequency(str(records[i].seq))),numpy.array(frequency(str(records[j].seq))))))
    return res


def parallel_distance(genome,dist="JSD",ksize=4,strand="both"):
    args=[]
    dist_pairs=[]
    global records
    records = list(SeqIO.parse(options.genome, "fasta"))
    for i in range(len(records)):
        for j in range(i+1,len(records)):
            args.append((i,j,dist,ksize,strand))
    #print(args)
    parallel_args_set = chunkitize(args,options.threads_max) 
    #print(parallel_args_set)
    pool = multiprocessing.Pool(processes=options.threads_max)
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






if __name__ == "__main__":

    Utilisation = "%prog [-i FILE] [options]"
    parser = OptionParser(usage = Utilisation)
    parser.add_option("-i","--assembly", dest = "genome", help = "multifasta of the genome assembly")
    #parser.add_option("-c","--conta", dest = "conta", help = "multifasta of the contaminant species training set")
    #parser.add_option("-r","--host", dest = "host", help = "optional multifasta of the host species training set")
    #parser.add_option("-n","--n_max_freq_in_windows", type = "float", dest = "n_max_freq_in_windows", default = "0.4", help = "maximum proportion of N tolerated in a window to compute the microcomposition anyway [0~1]. Too much 'N's will reduce the number of kmer counts and will artificially coarse the resolution of the frequencies. If your assembly contains many stretches of 'N's, consider rising this parameter and shortening the windows step in order to allow computation and signal in the output, it might cost you computational time though. Windows which do not meet this criteria will be affected a distance of 'nan'")
    parser.add_option("-k","--lgMot", dest = "k", type = "int", default = 4, help = "word wise / kmer lenght / k [default:%default]")
    #parser.add_option("-w","--windows_size", dest ="windows_size", type = "int", help = "Sliding windows size (bp)[default:%default]")
    #parser.add_option("-t","--windows_step", dest ="windows_step", type = "int", help = "Sliding windows step size(bp)[default:%default]")
    parser.add_option("-s","--strand", dest = "strand", default = "both", help = "strand used to compute microcomposition. leading, lagging ou both [default:%default]")
    parser.add_option("-d","--distance", dest = "dist", default = "JSD", help = "méthode de distance entre 2 signatures : KL: Kullback-Leibler, Eucl : Euclidienne[default:%default], JSD : Jensen-Shannon divergence")
    parser.add_option("-u","--cpu", dest = "threads_max", type = "int", default = 4, help = "how many threads to use for windows microcomposition computation[default:%default]")
    parser.add_option("-o","--out", dest = "out_file", default = "phyloligo.out", help = "output file[default:%default]")

    (options,argument) = parser.parse_args()


    if(not options.genome):
        parser.error("An input fasta file (-i ) is mandatory")
        exit()
    else:
        #print("A genome was provided")
        res = parallel_distance(options.genome,options.dist,options.k,options.strand)
        print(res)
        numpy.savetxt(options.out_file, res, delimiter="\t")
      
       
      
      

