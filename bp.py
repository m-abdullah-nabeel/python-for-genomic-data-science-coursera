#! /usr/bin/python3
from Bio.Seq import Seq
from Bio import SeqIO


for seq_record in SeqIO.parse("dna.example.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    print("...")
