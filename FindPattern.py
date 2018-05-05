# Input:  Astring Text and a DNA string Pattern
# Output:  A number of Pattern occurance in Text
def PatternCount(Text,Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

import re
def PatternMatching(Pattern, Genome):
    return [m.start() for m in re.finditer('(?='+Pattern+')', Genome)]
# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        # add each key to words whose corresponding frequency value is equal to m
        if freq[key]==m:
            words.append(key)
    return words
# Copy your FrequencyMap() function here.
def FrequencyMap(Text, k):
    # your code here
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = freq.get(Pattern,0)+1
    # hint: your code goes here!
    return freq

# Input:  Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)
def SlowSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
    return array


# Input:  Strings Genome and symbol
# Output: FasterSymbolArray(Genome, symbol)
## In the improved version, we just need to excute PatternCount once
## Which significantly reduce the running time of the algorithm from
## n*n//2 which is O(n^2) to n+n//2 which is O(n).
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

# Input:  A String Genome
# Output: The skew array of Genome as a list.
def SkewArray(Genome):
    # your code here
    skew = {}
    skew[0] = 0
    for i in range(1,len(Genome)+1):
        if (Genome[i-1] == 'G'):
            skew[i] = skew[i-1] + 1
        elif (Genome[i-1] == 'C'):
            skew[i] = skew[i-1] - 1
        else:
            skew[i] = skew[i-1]
    return skew

# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
	# generate an empty list positions
    positions = [] # output variable
    # set a variable equal to SkewArray(Genome)
    skew = SkewArray(Genome)
    # find the minimum value of all values in the skew array
    ori_pos = min(skew.values())
    # range over the length of the skew array and add all positions achieving the min to positions
    for i in range(0,len(skew)):
    	if (skew[i] == ori_pos):
    		positions.append(i)
    return positions

# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    ham_dis = 0
    for i in range(len(p)):
        if (p[i] != q[i]):
            ham_dis += 1
    return ham_dis


# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)],Pattern) <= d:
            positions.append(i)
    # your code here
    return positions

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches    
def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)],Pattern)<=d:
            count = count+1
    return count