#!/usr/bin/env python3
"""
This program computes oligonucleotide profiles (microcomposition) for sequences given as arguments.
Depending on the set of parameters used, the program guess what is possible to do and will perform accordingly. Use the proper (and minimal) parameter set you need to achieve what you want.
Luckily, the program will tell you what it does, and how you can change its behaviour by providing him with supplementary parameters.
See the help by calling the program without any argument.


dependencies: 
  Biopython
  numpy
  cython

as root/admin if wanted installed globally
aptitude/yum install python3-dev python3-setuptools


easy_install -U setuptools
easy_install3 -U setuptools

pip install biopython
pip3 install biopython 

pip install cython
pip3 install cython

pip install numpy
pip3 install numpy
"""

__author__ = "Ludovic V. Mallet, PhD"
__copyright__ = ""
__date__ = "2016.03.22"
__licence__ = "GPLv3"
__version__ = "0.1"
__status__ = "alpha"
__email__ = ""

from scoop import futures

import os, re, math, sys, argparse, time
import tempfile
import numpy as np
import multiprocessing, shlex, subprocess, pickle, shutil

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



np.seterr(divide='ignore', invalid='ignore')
#TBF: wtf way to long ! (-:
#minimum_number_of_windows_per_fasta_entry_to_use_multiple_cpu_for_this_entry=20
#LM Fine by me, this very variable makes the implementation very much complicated and triggers the multicore or serial mode depending on the number of windows in a contig. therefore I wanted to be explicit on what it does. But I guess this mere sentence explqins it now ^^.
min_nb_w_per_fasta_for_mul_cpu=20

## Distance functions

def posdef_check_value(d):
    d[np.isnan(d)]=0    
    d[np.isinf(d)]=0

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
    return np.sqrt(np.sum(d))*1000 # Scaling
    
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
    return d*1000 # Scaling

def vector_to_matrix(profile):
    return list((zip(*(iter(profile),)*int(math.sqrt(len(profile))))))
  
def read_seq_chunk_pos(sequences, chunksize):
    """ read a first chunk of fasta sequences to be processed
    
    Parameters:
    -----------
    sequences: list
        list of nucleotide sequences
    chunksize: int
        the number of fasta entries to read
    
    Return:
    -------
    seqchunk: list
        a list of SeqRecord 
    """
    seqchunk = list()
    start = 0
    for seq in sequences:
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
    sequences: list
        list of nucleotide sequences
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

def cut_sequence_and_count_pattern(seq, pattern, strand):
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
    seq = select_strand(seq, strand)
    # we work on upper case letters
    seq = seq.upper()
    
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
    return count_words

def count2freq(count_words, ksize):
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
    kword_count = sum(count_words.values())
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

def compute_frequency(seq, n_max_freq_in_windows=1.0, pattern="1111", strand="both"):
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
    
    pattern=str(pattern)
    ksize = pattern.count("1")
    if((seq.count('N')/len(seq)) <= float(n_max_freq_in_windows)):
        count_words = cut_sequence_and_count_pattern(seq, pattern, strand)
        # create feature vector
        features = count2freq(count_words, ksize)
    else:
        features = np.array([np.nan] * ksize**4)
    return features

def compute_whole_composition(genome, pattern, strand, nb_jobs=1):
    """ for each sequence of the genome, separately compute words numbers,
    aggregate results and convert to frequency
    """
    pattern=str(pattern)
    ksize = pattern.count("1")
    fd = delayed(cut_sequence_and_count_pattern)
    counts = Parallel(n_jobs=nb_jobs, verbose=0)(fd(str(record.seq), pattern, strand) 
                                            for record in SeqIO.parse(genome, "fasta"))
    # aggregate
    count_words = Counter()
    for par_counts in counts:
        for word in par_counts:
            count_words[word] += par_counts[word]
    # create feature vector
    frequency = count2freq(count_words, ksize)
    return frequency

    
def compute_distance_joblib(mth_dist, mcp, seq, pattern, strand , n_max_freq_in_windows):
    freq = compute_frequency(seq, n_max_freq_in_windows, pattern, strand)
    if mth_dist == "JSD":
        return JSD(freq, mcp)
    else:
        return Eucl(freq, mcp)


def compute_distances(mthdrun, large, mth_dist, mcp, sequences, pattern, strand, njobs, n_max_freq_in_windows):
    
    if mthdrun == "joblib":
        if large == "None":
            fd = delayed(compute_distance_joblib)
            res = Parallel(n_jobs=njobs, verbose=0)(fd(mth_dist, mcp, sequences[i], pattern, strand , n_max_freq_in_windows) 
                                                    for i in range(len(sequences)))
            res = np.array(res)
            #print(res.shape)
    return res

def make_genome_chunk(genome, windows_size, windows_step, options, nbchunk=500):
    """ create chunk of genome sequences
    """
    chunk_info = list()
    chunk_sequences = list()
    for record in SeqIO.parse(genome, "fasta"):
        seq = str(record.seq)
        seq_id = record.id
        
        if len(seq) < windows_size: #only enough to compute one window, no sliding,
            chunk_info.append([seq_id, 0, int(len(seq))])
            chunk_sequences.append(seq)
                
            if len(chunk_info) == nbchunk:
                yield chunk_info, chunk_sequences
                chunk_info, chunk_sequences = list(), list()
                
        elif len(seq) < min_nb_w_per_fasta_for_mul_cpu*windows_step:
            #not many windows in this contig, so better launching it in serial rather than in parallel
            for s in range(0, len(seq)-windows_size, windows_step):
                if s==0:
                    # to avoid border effects being a problem in the code, we use only the simple formula 
                    # start+windows_size/2 (+/-) windows_step/2 to find the significant center part of the windows. 
                    # when several windows overlap, this centered part, long as the windows step is the most representative of the windows,
                    # not representing as much other part of this window that are overlapped by other windows. 
                    # BUT: this simple formula has border effects, so we manually correct the start of the first window and 
                    # the stop of the last window to match the contig borders.
                    displayed_start=1
                else:
                    displayed_start=int(s+windows_size/2-windows_step/2)

                if s==len(seq)-windows_size:
                    displayed_stop=len(seq)
                else:
                    displayed_stop=int(s+windows_size/2+windows_step/2)

                window=seq[s:s+windows_size]
                chunk_info.append([seq_id, displayed_start, displayed_stop])
                chunk_sequences.append(window)
                
                if len(chunk_info) == nbchunk:
                    yield chunk_info, chunk_sequences
                    chunk_info, chunk_sequences = list(), list()
                    
        else:
            for s in range(0,len(seq)-windows_size,windows_step):
                start, stop = int(s+windows_size/2-windows_step/2),int(s+windows_size/2+windows_step/2)
                if start == (windows_size/2-windows_step/2):
                    displayed_start=1
                else:
                    displayed_start=start
                if stop-windows_step/2+windows_size/2 >= len(seq)-windows_step and stop-windows_step/2+windows_size/2 <= len(seq):
                    displayed_stop = len(seq)
                else:
                    displayed_stop = stop
                    
                chunk_info.append([seq_id, displayed_start, displayed_stop])
                chunk_sequences.append(seq[s:s+windows_size])
                if len(chunk_info) == nbchunk:
                    yield chunk_info, chunk_sequences
                    chunk_info, chunk_sequences = list(), list()

    # last chunk to return
    if len(chunk_info) != 0:
        yield chunk_info, chunk_sequences

def sliding_windows_distances(genome, mcp_comparison, mth_dist="JSD", pattern="1111", windows_size=5000, windows_step=500, options=None):
                              
    """ Stuff
    
    Parameters:
    ===========
    seq: Seq object
        the processed sequence
    seq_id: string 
        the id of the processed sequence 
    mcp_comparison: np.array 
        the kmer frequencies of the query (host or genome) 
    position: bool
        ???
    mth_dist: string
        the distance method
    windows_size: int 
        the size of the sequence to process 
    windows_step: int 
        the number of nucleotide to slid the window
    pattern: string
        kmer pattern with gap
    """
    
    cnt = 0
    t1_tot = time.time()
    for chunk_info, sequences in make_genome_chunk(genome, windows_size, windows_step, options, 50000):
        res = list()
        if mth_dist == "JSD":
            vec_dist = compute_distances("joblib", "None", mth_dist, mcp_comparison, 
                                        sequences, pattern, options.strand, 
                                        options.threads_max, options.n_max_freq_in_windows)
            for i in range(vec_dist.shape[0]):
                res.append(chunk_info[i]+[vec_dist[i]])
        else:
            #freq  = compute_frequencies("joblib", "None", sequences, pattern, 
                                    #options.strand, options.threads_max, 
                                    #options.n_max_freq_in_windows, options.workdir)
            #t2 = time.time()
            vec_dist = pairwise_distances(freq, mcp_comparison, n_jobs=options.threads_max)
            #print(np.allclose(test_dist, vec_dist))
            #print(vec_dist.shape)
            for i in range(vec_dist.shape[0]):
                res.append(chunk_info[i]+[vec_dist[i]])
        yield res 
        #t3 = time.time()
        #print(cnt, t2-t1, t3-t2)
        #cnt += 1
    #t2_tot = time.time()
    #print("Done in", t2_tot - t1_tot)
    
def concate_sequence(inputfile, file_format="fasta"):
    """ agglomerate all sequences in memory in a single Seq object
    
    Parammeters:
    ============
    inputfile: string
        path to the sequence file
    file_format: string
        format of the sequence file (see BioPython SeqIO.parse for a list of available format
        
    Return:
    =======
    whole_seq: string
        a unique string of sequence separated by N, k-mer with N are not taken into account (boundaries of sequences)
    """
    seq = []
    for record in SeqIO.parse(inputfile, file_format):
        seq.append(str(record.seq))
    # N is not interpreted to compute frequencies, so by putting one between two sequences, 
    # we avoid to account for the chimeric words that would be created by juxtaposing the 2 sequences.
    whole_seq = "N".join(seq)
    return whole_seq

def get_cmd():
    """ read command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--assembly", action="store", required=True, 
                    dest="genome", help="multifasta of the genome assembly")
    parser.add_argument("-c", "--conta", action="store", dest="conta", 
                    help="multifasta of the contaminant species training set")
    parser.add_argument("-r", "--host", action="store", dest="host", 
                    help="optional host species training set in multifasta")
    parser.add_argument("-n", "--n_max_freq_in_windows", action="store", 
                    type=float, dest="n_max_freq_in_windows", default = 0.4,
                    help = "maximum proportion of N tolerated in a window to compute the microcomposition anyway [0~1]. "
                    "Too much 'N's will reduce the number of kmer counts and will artificially coarse the frequency resolutions")
                    #"If your assembly contains many stretches of 'N's, "
                    #"consider rising this parameter and shortening the windows"
                    #" step in order to allow computation and signal in the "
                    #"output, it might cost you computational time though. "
                    #"Windows which do not meet this criteria will be affected "
                    #"a distance of 'nan'")
    parser.add_argument("-k", "--lgMot", action="store", dest="k", type=int, default=4, 
                    help="word wise/ kmer lenght/ k [default:%(default)d]")
    parser.add_argument("-p", "--pattern", action="store", dest="pattern",
                    help="pattern to use for frequency computation")
    parser.add_argument("-w", "--windows_size", action="store", dest="windows_size", type=int, default=5000,
                    help="Sliding windows size (bp)")
    parser.add_argument("-t", "--windows_step", action="store", dest="windows_step", type=int, default=500, 
                    help="Sliding windows step size(bp)")
    parser.add_argument("-d", "--distance", action="store", dest="dist", choices=["JSD", "Eucl"],
                    default="JSD", help="distance method between two signatures: "
                    "Eucl : Euclidienne, JSD : Jensen-Shannon divergence [default:%(default)s]")
    parser.add_argument("-s", "--strand", action="store", default="both",  choices=["both", "plus", "minus"],
                    help="strand used to compute microcomposition. [default:%(default)s]")
    parser.add_argument("-u", "--cpu", action="store", dest="threads_max", type=int, default=4, 
                    help="how maany threads to use for windows microcomposition computation[default:%(default)d]")
    parser.add_argument("-W", "--workdir", action="store", dest="workdir", default=".", help="working directory")

    options = parser.parse_args()

    return options

def main():

    # get parameters
    options = get_cmd()
    
    print("Genome : {}".format(options.genome))
    
    base_genome = os.path.basename(options.genome)
    
    # preparing output file
    if not os.path.isdir(options.workdir):
        os.makedirs(options.workdir)
    
    # read target sequence (host or genome)
    if (not options.conta):
        #no contaminant, genome is the target
        target = options.genome
        print("Contaminant : {}".format(None))
        output = os.path.join(options.workdir, base_genome + ".mcp_windows_vs_whole_" + options.dist+".dist")
    else:
        base_conta = os.path.basename(options.conta)
        print("Contaminant : {} ".format(options.conta))
        output = base_genome+".mcp_hostwindows_vs_"
        if options.host:
            base_host = os.path.basename(options.host)
            # the host is target
            print("Host : {}".format(options.host))
            target = options.host
            output = os.path.join(options.workdir, output+"host_"+base_host+"_"+options.dist+".dist")
        else: 
            #contaminant is provided but no host, genome is the target
            print("Host : None, using whole genome".format(options.host))
            output = os.path.join(options.workdir, output+"wholegenome_"+options.dist+".dist")
            target = options.genome
        
    if type(options.pattern) == int and not options.k:
        options.pattern = "1" * options.pattern
    elif not options.pattern and options.k:
        options.pattern = "1" * options.k
        
    # one vector shape of ksize**4
    genome =  compute_whole_composition(options.genome, options.pattern, options.strand, nb_jobs=options.threads_max)
    
    if (not options.conta):
        
        if (not options.windows_size and not options.windows_step):
            print("Warning, no sliding window parameters (-w and -t )\n"
                "The signature will be computed from the whole genome\n"
                "Computing signature from the whole genome", file = sys.stderr)
            
            output = os.path.join(options.workdir, base_genome+".microcomposition.mat")
            with open(output, 'w') as outf:
                outf.write(str(vector_to_matrix(genome)))
            
            sys.exit(0)
        elif (options.windows_size or options.windows_step):
            print("Computing microcomposition signaure and distances to genome")
            

    else:
        # one vector shape of ksize**4
        conta =  compute_whole_composition(options.conta, options.pattern, options.strand, nb_jobs=options.threads_max)
          
    
    with open(output, 'w') as outf:
        for res in sliding_windows_distances(options.genome, mcp_comparison=genome, mth_dist=options.dist, pattern=options.pattern,
                                             windows_size=options.windows_size, windows_step=options.windows_step, options=options):
            for t in res:
                outf.write(str("\t".join(map(str,t)))+"\n")

    if options.conta:
        output = os.path.join(options.workdir, base_genome+".mcp_hostwindows_vs_"+"conta_"+base_conta+"_"+options.dist+".dist")
        with open(output, 'w') as outf:
            for res in sliding_windows_distances(options.genome, mcp_comparison=conta, mth_dist=options.dist, pattern=options.pattern,
                                             windows_size=options.windows_size, windows_step=options.windows_step, options=options):
                for t in res:
                    outf.write(str("\t".join(map(str,t)))+"\n")
                
    sys.exit(0)

if __name__ == "__main__":
    main()

