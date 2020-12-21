

from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
from Bio import SeqUtils

table={
        "F":["UUU", "UUC"], 
        "L":["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
        "S":["UCU", "UCC", "UCA", "UCG"],
        "Y":["UAU", "UAC"], 
        "C":["UGU", "UGC"], 
        "W":["UGG"],
        "P":["CCU", "CCC", "CCA", "CCG"],
        "H":["CAU", "CAC"], 
        "Q":["CAA", "CAG"],
        "R":["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "I":["AUU", "AUC", "AUA"], 
        "M":["AUG"],
        "T":["ACU", "ACC", "ACA", "ACG"],
        "N":["AAU", "AAC"], 
        "K":["AAA", "AAG"],
        "S":["AGU", "AGC"],
        "V":["GUU", "GUC", "GUA", "GUG"],
        "A":["GCU", "GCC", "GCA", "GCG"],
        "D":["GAU", "GAC"], 
        "E":["GAA", "GAG"],
        "G":["GGU", "GGC", "GGA", "GGG"]
    }

identifier = ['Seq_compl']
for record in SeqIO.parse('short_protein.fasta', "fasta"):
    id = record.id
    if id in identifier:
        print(record.seq)

def count_tRNA(file, identifier):
    product = 1
    for record in SeqIO.parse(file, 'fasta'):
        id = record.id
        if id in identifier:
            for aa in record.seq:
                if isinstance(table[aa], list):          #'list' give a bouillian response True if the values of the dict are a list, False if not
                    product *= len(table[aa])
            return (product)
        else:
            print('error')

import itertools
from itertools import product

def print_seq (file, identifier):
    for record in SeqIO.parse(file, 'fasta'):
        id = record.id
        seq = ()
        if id in identifier:
            seq = list(record.seq)
            return seq

def give_combination (file, identifier):
    key_list = print_seq(file, identifier)
    res = [table[key] for key in key_list]
    if count_tRNA(file, identifier) <= 5000:
        print(list(product(*res)))
    elif count_tRNA(file, identifier) >= 5000:
        print('Too many possibilities')

identifier = ['seq_compl']

count_tRNA('short_protein.fasta', identifier)

give_combination('short_protein.fasta', identifier)
