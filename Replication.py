# Input:  A string Pattern
# Output: The reverse of Pattern
def Reverse(Pattern):
    # your code here
    return Pattern[::-1]


# Input:  A DNA string Pattern
# Output: The complementary string of Pattern (with every nucleotide replaced by its complement).

def Complement(Pattern):
    # your code here
    di = {
    	"A":"T",
    	"T":"A",
    	"G":"C",
    	"C":"G"
	}
    result = ""
    for char in Pattern.upper():
        result+= di[char]
    return result

# Input:  A DNA string Pattern
# Output: The reverse complement of Pattern
def ReverseComplement(Pattern):   
    # your code here
    import string
    # str.maketrans() is built in (Python 3.6, anyway); it maps any character in the first string to the one in the           same position of the second
    trantab = string.maketrans("ATCG", "TAGC")
    # How to reverse a string in Python. You either know it or you don't, but once you know it you'll never forget it
    rev = Pattern[::-1]
    # Use the translation table defined above to do the complementary bit
    return rev.translate(trantab)




