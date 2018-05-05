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