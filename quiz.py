from FindPattern import PatternCount
print PatternCount("ACTGTACGATGATGTGTGTCAAAG","TGT")

from FindPattern import FrequentWords
print FrequentWords("CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT",3)

from Replication import ReverseComplement
print ReverseComplement("CCAGATC")

from FindPattern import SkewArray
print max(SkewArray("GATACACTTCCCGAGTAGGTACTG").values())


from FindPattern import HammingDistance
print HammingDistance("CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG","ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT")


a=list(range(5))
b=a
a[2]=12
print b

text = "aaaab"
print text[0:0+5]

Motifs = [
"TCGGGGGTTTTT",
"CCGGTGACTTAC",
"ACGGGGATTTTC",
"TTGGGGACTTTT",
"AAGGGGACTTCC",
"TTGGGGACTTCC",
"TCGGGGATTCAT",
"TCGGGGATTCCT",
"TAGGGGAACTAC",
"TCGGGTATAACC"
]
from Motifs import EntropyScore
print EntropyScore(Motifs)

print 0.3*0.3*1*0.1*0.5*0.9

import random


prob = {'AA': 0.2, 'TT':0.2, 'CC':0.1, 'GG':0.1, 'AT':0.4}


prob1 = {'AA': 0.2, 'TT':0.2, 'CC':0.1, 'GG':0.1, 'AT':0.4}

print Motifs[0:0]

from Motifs import *
Dna=["AAGCCAAA",
"AATCCTGG",
"GCTACTTG",
"ATGTTTTG"]
motifs = ["CCA",
"CCT",
"CTT",
"TTG"]

print Motifs(Profile(motifs),Dna)

prob1 = {'1':0.22,'2':0.54,'3':0.58,'4':0.36,'5':0.3}
prob2 = {'1':0.15,'2':0.6,'3':0.225,'4':0.225,'5':0.3}
print Normalize(prob2)
