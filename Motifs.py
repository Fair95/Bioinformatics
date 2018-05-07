# Input:  A set of kmers Motifs
# Output: Count(Motifs)

def Count(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
	    count[symbol] = []
	    for j in range(k):
	        count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    for symbol in "ACGT":
	    profile[symbol] = []
	    for j in range(k):
	        profile[symbol].append(0)
    
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            profile[symbol][j] += 1.0/t
    return profile

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def Score(Motifs):
    k = len(Motifs)
    t = len(Motifs[0])
    score = 0
    consensus = Consensus(Motifs)
    for i in range(k):
        for j in range(t):
            if (Motifs[i][j] != consensus[j]):
                score+=1
    return score

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    k = len(Text)
    prob = 1
    for j in range(k):
        prob *= Profile[Text[j]][j]
    return prob

# Input:  String Text, an integer k, and profile matrix Profile
# Output: ProfileMostProbablePattern(Text, k, Profile)
def ProfileMostProbablePattern(Text, k, Profile):
    t = len(Text)
    max = -1
    pattern = ""
    for i in range(0,t-k+1):
        if (Pr(Text[i:i+k],Profile)>max):
            max = Pr(Text[i:i+k],Profile)
            pattern = Text[i:i+k]
    return pattern

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
	## Our proposed greedy motif search algorithm, GreedyMotifSearch,
	## starts by setting BestMotifs equal to the first k-mer from each string in Dna. 
	## These strings will serve as the best-scoring motifs found thus far.
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    ## It then ranges over all possible k-mers in Dna[0], trying each one as Motifs[0]
    for i in range(n-k+1):
        Motifs = []
        ## Firstly we initiate by setting Motifs[0] as the first k-mer in Dna[0]
        Motifs.append(Dna[0][i:i+k])
        ## Find the best Motifs[1] based in Motifs[0],
        ## then Motifs[2] based on Motifs[0] and Motifs[1]
        ## and continue until Motifs[j] is assigned.
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        ## Then we check if for the given first k-mer, the score is better (lower)
        ## and in next interation we check the second k-mer in Dna[0]
        ## continue until we reach the last possible k-mer in Dna[0]
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

## Entropy Function summation from 0 to i of (-pi*log(pi,2))
## Improved version of Score that address the distribution of data
def EntropyScore(Motifs):
	import math
   	k = len(Motifs)
   	t = len(Motifs[0])
   	score = 0
   	profile = Profile(Motifs)
   	for char in "ACGT":
   		for i in range(t):
   			pi = profile[char][i]
   			if (pi != 0):
   				score -= pi * math.log(pi,2)
   	return score


## Improved version of Count which instead of treating event with
## probability of zero, it becomes an event that very unlikely to happen
## by using Laplace's Rule of Succession
# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
	    count[symbol] = []
	    for j in range(k):
	        count[symbol].append(1)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    k = len(Motifs[0])
    t = len(Motifs)
    dic = {Symbol:[1/(t+4)]*k for Symbol in 'ATGC'}
    for cid in range(k):
        for row in Motifs:
            dic[row[cid]][cid] += 1/(t+4)
    return dic

    # t = len(Motifs)
    # k = len(Motifs[0])
    # profile = {} # output variable
    # for symbol in "ACGT":
	   #  profile[symbol] = []
	   #  for j in range(k):
	   #      profile[symbol].append(1.0/(t+4.0))
    
    # for i in range(t):
    #     for j in range(k):
    #         symbol = Motifs[i][j]
    #         profile[symbol][j] += 1.0/(t+4.0)
    # return profile

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    ## It then ranges over all possible k-mers in Dna[0], trying each one as Motifs[0]
    for i in range(n-k+1):
        Motifs = []
        ## Firstly we initiate by setting Motifs[0] as the first k-mer in Dna[0]
        Motifs.append(Dna[0][i:i+k])
        ## Find the best Motifs[1] based in Motifs[0],
        ## then Motifs[2] based on Motifs[0] and Motifs[1]
        ## and continue until Motifs[j] is assigned.
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        ## Then we check if for the given first k-mer, the score is better (lower)
        ## and in next interation we check the second k-mer in Dna[0]
        ## continue until we reach the last possible k-mer in Dna[0]
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

#################################################
### randomize algorithm approach starts here  ###
#################################################
# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)

## Give a profile and strings of Dna, find the most propable pattern for
## each string of Dna, i.e. for each Dna[i].
def Motifs(Profile, Dna):
    motifs = []
    for i in range(len(Dna)):
        motifs.append(ProfileMostProbablePattern(Dna[i], len(Profile['A']),Profile))
    return motifs


# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient
def RandomMotifs(Dna, k, t):
    import random
    motifs = []
    for i in range(t):
        start = random.randint(0,len(Dna[0])-k)
        end = start + k
        motifs.append(Dna[i][start:end])
    return motifs

# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: RandomizedMotifSearch(Dna, k, t)
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 

## We repeat the above function N times until the score drops to find
## the best random generated motif
def RepeatedRandomizedMotifSearch(Dna, k, t, N):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(N):
        Motifs = RandomizedMotifSearch(Dna, k, t)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs

# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities

def Normalize(Probabilities):
    # your code here
    total = sum(Probabilities.values())
    for key in Probabilities:
        Probabilities[key] /= total
    return Probabilities

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    import random
    kmer = '' # output variable
    prob = Normalize(Probabilities)
    accum = 0.0
    roll = random.uniform(0.0,1.0)
    for key in prob:
        accum+=prob[key]
        if (roll<accum):
        	kmer = key
        	break
    return kmer

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
def GibbsSampler(Dna, k, t, N):
    import random
    BestMotifs = [] # output variable
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    for j in range(N):
        i = random.randint(0,t-1)
        profile = ProfileWithPseudocounts(M[0:i]+M[i+1:len(M)])
        M[i] = ProfileGeneratedString(Dna[i],profile,k)
        if Score(M)<Score(BestMotifs):
            BestMotifs = M
    return BestMotifs

def RepeatedGibbsSampler(Dna, k, t, N, R):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(R):
        Motifs = GibbsSampler(Dna, k, t, N)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs