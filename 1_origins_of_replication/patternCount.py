#!/usr/bin/env python

import string, sys

COMP_TRANS = string.maketrans("ACGT", "TGCA")
BASE_CODE = [ 'A','C','G','T' ]

def symbolToNumber(base):
    if base=='A': return 0
    elif base=='C': return 1
    elif base=='G': return 2
    elif base=='T': return 3
    else: return -1

# converts DNA sequence to a numerical hash code
def patternToNumber(pattern):
    BASE_CODE = { 'A':0, 'C':1, 'G':2, 'T':3 }
    pattern = pattern.upper()[::-1]
    multiplier = 1
    retVal = 0
    
    for i in range(len(pattern)):
        retVal += symbolToNumber( pattern[i] ) * multiplier
        multiplier *= 4
        
    return retVal
    
# converts numerical hash code back to DNA sequence
def numberToPattern_old(number, patternLength):
    divisor = 4 ** (patternLength - 1)
    retVal = ''
    
    for i in range(patternLength):
        retVal += BASE_CODE[ number / divisor ]
        number = number % divisor
        divisor /= 4
        
    return retVal
    
def numberToPattern(index, k):
    if k==1:
        return BASE_CODE[index]
    prefixIndex = index / 4
    remainder = index % 4
    prefixPattern = numberToPattern(prefixIndex, k-1)
    symbol = BASE_CODE[remainder]
    return prefixPattern + symbol
    
        
# compute frequency
def computingFreqs(seq, k):
    freqArray = [0]*(4**k)
    for i in range(len(seq)-k+1):
        pattern = seq[i:(i+k)]
        #print pattern
        j = patternToNumber(pattern)
        freqArray[j] += 1
        
    return freqArray
        
        
# freqArray = computingFreqs('ACAAGCTTACGCGCGGCATTCCTAGGTAAAGACTCTGCGCAGAGCATACATCACTCTATGTGGTATGGTGACCCTGGTGACCTAGAAGGCATGGTGCTTCATGTAAAAGAGAGAAAGTAACTCCCACGAATCGTACAGTCAGACACTCACGGACTGCCCCAGTCCCGGCGTTGTACCTGTACGGTTCCGTTCGAGTGTCCAAGTCTAAACCCCTCACGTATGACTGAGAGGTGCTTACATGTTACTGGAGTCTGGACCCTCTTTAGACTGGAGGCCATTATATCTGCGGCCGTGTCCGTTGTTTAATGGACGTTTAAGCAGCTGTGTGCCGGGACGGCCGCGTGCAGTGTGGAAGGCCAAATGTTGGTGGCGGGTAAAATCTGATCGGCGTTGTACATTTAGCGTATATAGCGCTCAACTAGTGACCTCCACTTCTTCGTCGAATCTGACTACGCTGATCCATAGCGAAGCTATGACTATGCAGGACGGTGCCGGCCACAGATACGAAGGACATTTGAACTCTTCTTAGAAGAATTCCCCCCGCAGCCCTCAAAGCCGGTGCGAACCGCTGCAAGCTATTGCTAAGCCCTTCTATGGATTCGGGATAATAATTACGTAGACGAGTATTAATCAAAAGAAAGTTCCGCATACGCGCGTTCTTCTCAGGTCGCTTTCACCCCCTACACTATAGATGGTCATACTTTAACTGTCCGCAGCACGCTTATAAGAGAGGTTTGAGAATTAGCCTGCCTTCACTC',7)
# print " ".join([str(x) for x in freqArray])
# sys.exit(0)

def patternCount(seq, pattern):
    count = 0
    for i in range(len(seq)-len(pattern)+1):
        #print seq[i:(i+len(pattern))]
        if seq[i:(i+len(pattern))]==pattern:
            count += 1
    
    return count
    
#print patternCount(seq, pattern)


def frequentWords(seq, k):
    freqPatterns = set()
    count = [0]*(len(seq)-k+1)
    
    for i in range(len(seq)-k+1):
        pattern = seq[i:(i+k)]
        count[i] = patternCount(seq, pattern)
    maxCount = max(count)
    print count, maxCount
    
    for i in range(len(seq)-k+1):
        if count[i]==maxCount:
            freqPatterns.add(seq[i:(i+k)])
            
    return sorted(freqPatterns)
    
# reverse complements a string
def revComp(seq):
    return seq.translate(COMP_TRANS)[::-1]
    
#print revComp('ATAACTAACAATCTACTTAACCTTACGATCATGGCCCACCCTCGGATG')


def findPatternPositions(pattern, genome):
    matchPositions = []
    n = len(pattern)
    
    for i in range(len(genome)-n):
        if genome[i:(i+n)]==pattern:
            matchPositions.append(i)
        
    return matchPositions
    
# vcFile = open('Vibrio_cholerae.txt','r')
# vc_genome = vcFile.read()
# vcFile.close()

# matchPositions = findPatternPositions('ATGATCAAG', vc_genome)
# print ' '.join([str(x) for x in matchPositions])


def findPatternClump(genome, k, L, t):
    patterns = []
    clump = [0]*(4**k)
    freqArray = computingFreqs(genome[:L], k)
    
    for i in range(len(freqArray)):
        if freqArray[i] >= t:
            clump[i] = 1
            
    for i in range(len(genome)-L):
        j = patternToNumber( genome[i:(k+i)] )
        freqArray[j] = freqArray[j] - 1
        j = patternToNumber( genome[(i+L-k+1):(i+L+1)] )
        freqArray[j] = freqArray[j] + 1
        if freqArray[j] >= t:
            clump[j] = 1
            
    for i in range(len(clump)):
        if clump[i]==1:
            pattern = numberToPattern(i, k)
            patterns.append(pattern)
            print i, pattern, freqArray[i]
            
    return patterns

#print findPatternClump('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA', 5, 50, 4)

# eColiFile = open('E-coli.txt', 'r')
# eColiGenome = eColiFile.read()
# eColiFile.close()
# eColiClumps = findPatternClump(eColiGenome, 9, 500, 3)
# print len(eColiClumps)

#print findPatternClump('GCACAAGGCCGACAATAGGACGTAGCCTTGAAGACGACGTAGCGTGGTCGCATAAGTACAGTAGATAGTACCTCCCCCGCGCATCCTATTATTAAGTTAATT', 4, 30, 3)


def skew(seq):
    skewVals = [0]*(len(seq)+1)
    
    for i in range(1, len(seq)+1):
        if seq[i-1]=='G':
            skewVals[i] = skewVals[i-1] + 1
        elif seq[i-1]=='C':
            skewVals[i] = skewVals[i-1] - 1
        else:
            skewVals[i] = skewVals[i-1]
            
    return skewVals
        
# skewVals = skew('GAGCCACCGCGATA')
# print ' '.join([ str(x) for x in skewVals ])
#print skew('GATACACTTCCCAGTAGGTACTG')

def minimumSkew(seq):
    skewVals = skew(seq)
    minVal = min(skewVals)
    retList = []
    
    for i in range(len(skewVals)):
        if skewVals[i]==minVal:
            retList.append(i)
        
    return retList
    
# print minimumSkew('CCCCGTTTCATGCAACATTTGGATCTAGCGCTCCGCACGGATGGCGAAGACAGGGGAATAGCTCGGCAGTTAGGGTCTGGGCTAAAATATAAGTGGGCGCCATAACATCAAAGAAGTTACAAGACCCTGTAAGCCCAGATCGGTTTATAGGAGAGGTTGCGTCGTACCAAGCCGGGTCGCCA')

def hammingDist(seq1, seq2):
    distance = 0
    
    for i in range(len(seq1)):
        if seq1[i]!=seq2[i]:
            distance += 1
            
    return distance
    
# print hammingDist('TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC', 'GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA')


def approximatePatternMatchPosition(pattern, text, d):
    positions = []
    
    for i in range(len(text) - len(pattern) + 1):
        if hammingDist(pattern, text[i:(i+len(pattern))]) <= d:
            positions.append(i)
            
    return positions
    

# pos = approximatePatternMatchPosition('AAAAA', 'AACAAGCTGATAAACATTTAAAGAG', 2)
# print ' '.join([ str(x) for x in pos ])

def approximatePatternCount(text, pattern, d):
    count = 0
    
    for i in range(len(text) - len(pattern) + 1):
        if hammingDist(pattern, text[i:(i+len(pattern))]) <= d:
            count += 1
            
    return count
   
# print approximatePatternCount('GGGCACTAGTACCTCAGACTGAACGCCGCACGAAGCTTTCAAGTCGAAAGTCCTGGAGGTCTGAGTCTTGTATGGTGTTAGACTCAAGTCTCATAATCCTGGCAAAAATACCTGAGACGGCAGATGCTTGTGGTATGCGTCGTCGTAGCAGGCGCACGTACACATGCCAGCAGAACCTGAGACATTCAATTATTCGATTCGTATTACGAGTTCATCGCAAAACAAAGGTCCGCAGCCCTTCTAAGTAACCAGTAGATAAGCACCACAGTACCAACTGTTGGTACTGCCCGGGACACGTATGAGTAGAGTTATCCATCCGGGCTTGGGGCCTGTCCGTCCCCGGGAGTAGATCTTAGGGGATTAGACTCTAGA', 'AGTAGAT', 3)

#print approximatePatternCount('CATGCCATTCGCATTGTCCCAGTGA','CCC',2)


def mostFreqApproxMatch(text, k, d):
    counts = [0]*(4**k)
    
    for i in range(len(counts)):
        pattern = numberToPattern(i, k)
        counts[i] = approximatePatternCount(text, pattern, d)
        if i % 10000==0:
            print '%d seqs counted' % i
    print 'done counting frequencies'
    
    maxCount = max(counts)
    maxCountPatterns = []
    for i in range(len(counts)):
        if counts[i]==maxCount:
            pattern = numberToPattern(i, k)
            maxCountPatterns.append(pattern)
            
    return maxCountPatterns
    
# print ' '.join( mostFreqApproxMatch('GCCCGCCCCTTCGCCCGCCCAAGGCGCAAGGCCGCTTCAAGGCCGCTTCCTTCGCGCAAGAAGCTTCCTTCAAGGCCCGCCGGCCGGCCGGCCGGCCCGCGCGCGCGCCCGCCCGCCCGCCGAAGGCCGGCCGGCCCGCCCGCCCGCCCGCGCAAGGCGCAAGGCCGCTTCGCGCCTTCGCCCGCCCAAGCTTCAAGCTTCAAGAAGCTTCGCCCGCCGAAGGCCGCTTCGCCCGCCCGCGCCTTCGCCCGCCGAAGCTTCGCGCAAGCTTCAAGGCCCGCGCGCCCGCCCGCCGAAGAAGCTTCGCCCGCGCGCGCGCCCCTTCGCCGCTTCGCCG', 9, 3) )

def mostFrequentApproxMatchRevComp(text, k, d):
    counts = [0]*(4**k)
    
    for i in range(len(counts)):
        pattern1 = numberToPattern(i, k)
        pattern2 = revComp(pattern1)
        j = patternToNumber(pattern2)
        idx = min(i,j)
        counts[ min(i,j) ] = approximatePatternCount(text, pattern1, d) + \
                             approximatePatternCount(text, pattern2, d)
        if i % 10000==0:
            print '%d seqs counted' % i
    print 'done counting frequencies'
    
    maxCount = max(counts)
    maxCountPatterns = []
    for i in range(len(counts)):
        if counts[i]==maxCount:
            pattern = numberToPattern(i, k)
            maxCountPatterns.append(pattern)
            maxCountPatterns.append( revComp(pattern) )
            
    return maxCountPatterns
            
# print mostFrequentApproxMatchRevComp('TTCATGTTGTCGATGTTTCGATGTGTCTTGTCGTCTGTCAGTCGTCGTCGTCTGTCATGTTTCGACATGTCGACATGTCACGATGTTGTTGTGTCCACGACGAGTCGTCCGATGTTGTGTCGTCCAGTCTGTCACAGTCGTCCATTTTCAGTCTGTCATGTGTCTGTTGTTGTTGTTTTGTCGATGTTGTTTCGATGTTTCGACACA', 8, 2)


def question9(pattern, d):
    numCloseSeqs = 0
    
    for i in range(4**len(pattern)):
        seq = numberToPattern(i, len(pattern))
        if hammingDist(pattern, seq) <= d:
            numCloseSeqs += 1
            
    return numCloseSeqs
    
print question9('TGCA', 3)
