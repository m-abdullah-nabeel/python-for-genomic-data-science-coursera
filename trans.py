#! /usr/bin/python3
from Bio.Seq import Seq

# my_seq = Seq("AGTACACTGGT")
# print(my_seq)
# print(my_seq.reverse_complement())

coding_dna = Seq("TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG")
print(coding_dna.translate())
