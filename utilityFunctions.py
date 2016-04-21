#!/usr/bin/env python

import string

BASE_CODE = [ 'A','C','G','T' ]
COMP_TRANS = string.maketrans("ACGT", "TGCA")
GEN_CODE_RNA = { 'AAA':'K','AAC':'N','AAG':'K','AAU':'N','ACA':'T','ACC':'T','ACG':'T','ACU':'T','AGA':'R','AGC':'S','AGG':'R','AGU':'S','AUA':'I','AUC':'I','AUG':'M','AUU':'I','CAA':'Q','CAC':'H','CAG':'Q','CAU':'H','CCA':'P','CCC':'P','CCG':'P','CCU':'P','CGA':'R','CGC':'R','CGG':'R','CGU':'R','CUA':'L','CUC':'L','CUG':'L','CUU':'L','GAA':'E','GAC':'D','GAG':'E','GAU':'D','GCA':'A','GCC':'A','GCG':'A','GCU':'A','GGA':'G','GGC':'G','GGG':'G','GGU':'G','GUA':'V','GUC':'V','GUG':'V','GUU':'V','UAA':' ','UAC':'Y','UAG':' ','UAU':'Y','UCA':'S','UCC':'S','UCG':'S','UCU':'S','UGA':' ','UGC':'C','UGG':'W','UGU':'C','UUA':'L','UUC':'F','UUG':'L','UUU':'F' }
GEN_CODE_DNA = { 'AAA':'K','AAC':'N','AAG':'K','AAT':'N','ACA':'T','ACC':'T','ACG':'T','ACT':'T','AGA':'R','AGC':'S','AGG':'R','AGT':'S','ATA':'I','ATC':'I','ATG':'M','ATT':'I','CAA':'Q','CAC':'H','CAG':'Q','CAT':'H','CCA':'P','CCC':'P','CCG':'P','CCT':'P','CGA':'R','CGC':'R','CGG':'R','CGT':'R','CTA':'L','CTC':'L','CTG':'L','CTT':'L','GAA':'E','GAC':'D','GAG':'E','GAT':'D','GCA':'A','GCC':'A','GCG':'A','GCT':'A','GGA':'G','GGC':'G','GGG':'G','GGT':'G','GTA':'V','GTC':'V','GTG':'V','GTT':'V','TAA':' ','TAC':'Y','TAG':' ','TAT':'Y','TCA':'S','TCC':'S','TCG':'S','TCT':'S','TGA':' ','TGC':'C','TGG':'W','TGT':'C','TTA':'L','TTC':'F','TTG':'L','TTT':'F' }

# reverse complements a string
def revComp(seq):
    return seq.translate(COMP_TRANS)[::-1]

# translate DNA/RNA to protein
def translate(seq, genCode, source='DNA'):
    seq = re.sub('T','U',seq)
    retVal = ''
    
    for pos in range(0, len(seq), 3):
        if source=='DNA':
            aa = GEN_CODE_DNA[ seq[pos:(pos+3)] ]
        else:
            aa = GEN_CODE_RNA[ seq[pos:(pos+3)] ]
        if aa is None:
            break
        retVal += aa
        
    return retVal

# return Hamming distance between two sequences
def hammingDist(seq1, seq2):
    distance = 0
    
    for i in range(len(seq1)):
        if seq1[i]!=seq2[i]:
            distance += 1
            
    return distance
    
def symbolToNumber(base):
    if base=='A': return 0
    elif base=='C': return 1
    elif base=='G': return 2
    elif base=='T': return 3
    else: return -1
    
    
def numberToPattern(index, k):
    if k==1:
        return BASE_CODE[index]
    prefixIndex = index / 4
    remainder = index % 4
    prefixPattern = numberToPattern(prefixIndex, k-1)
    symbol = BASE_CODE[remainder]
    return prefixPattern + symbol