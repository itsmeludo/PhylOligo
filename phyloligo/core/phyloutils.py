#!/usr/bin/env python3

import sys
import shutil

from Bio import SeqIO
from Bio.Seq import Seq

def remove_folder(folder):
    try:
        shutil.rmtree(folder)
    except:
        print("Failed to delete folder: {}".format(folder))

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

def get_nb_records(fasta, file_format="fasta"):
    """ compute the number of sequences in a fasta file
    """
    return sum(1 for record in SeqIO.parse(fasta, file_format))
