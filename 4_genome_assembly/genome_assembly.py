#!/usr/bin/env python

import math, random, re, sys
sys.path.append('..')

from utilityFunctions import *


def stringComposition(seq, k):
    retVal = []
    for pos in range(len(seq)-k+1):
        retVal.append( seq[pos:(pos+k)] )
        
    return sorted(retVal)
    
# print '\n'.join(stringComposition('CAATCCAAC', 5))


def stringByGenomePath(kMers):
    retVal = kMers[0]
    for i in range(1, len(kMers)):
        retVal += kMers[i][-1]
        
    return retVal
    
# print stringByGenomePath(['ACCGA','CCGAA','CGAAG','GAAGC','AAGCT'])


def prefix(pattern):
    return pattern[:-1]
    
def suffix(pattern):
    return pattern[1:]

def overlapGraph(patterns):
    patterns = sorted(patterns)
    for i in range(len(patterns)):
        matches = []
        
        for j in range(len(patterns)):
            if i==j: continue
            if suffix(patterns[i])==prefix(patterns[j]):
                matches.append( patterns[j] )
                
        if len(matches) > 0:
            print '%s -> %s' % (patterns[i], ' '.join(matches))
            

# overlapGraph(['0000','0001','0010','0011','0100','0101','0110','0111','1000','1001','1010','1011','1100','1101','1110','1111'])

def extendBinaryString(strings):
    retVal = []
    for string in strings:
        retVal.append(string + '0')
        retVal.append(string + '1')
    return retVal

    
def repeated_kMer(string, k):
    if len(string) <= k: return False
    patternSet = set()
    
    for pos in range(len(string)-k+1):
        if string[pos:(pos+k)] in patternSet:
            return True
        patternSet.add( string[pos:(pos+k)] )
        
    return False


def universalBinaryStrings(k):
    strings = ['']
    strLength = 2**k + k - 1
    
    while len(strings) > 0 and len(strings[0]) < strLength:
        strings = extendBinaryString(strings)
        good_idx = []
        
        for i in range(len(strings)):
            if not repeated_kMer(strings[i], k):
                good_idx.append(i)
                
        strings = [ strings[x] for x in good_idx ]
        
    return strings
        
# print len(universalBinaryStrings(3))


def deBruijnAdjancenyList(string, k):
    adjHash = dict()
    
    for pos in range(len(string)-k+1):
       node1 = string[pos:(pos+k-1)]
       node2 = string[(pos+1):(pos+k)]
       if node1 in adjHash:
           adjHash[ node1 ] = '%s,%s' % (adjHash[ node1 ], node2)
       else:
           adjHash[ node1 ] = node2
       
    for node in sorted(adjHash):
        print '%s -> %s' %(node, adjHash[ node ])
                                             
#deBruijnAdjancenyList('AAGATTCTCTAAGA', 4)


def deBruijnFromKmers(kmers):
    adjHash = dict()
    
    for kmer in kmers:
        node1 = kmer[:-1]
        node2 = kmer[1:]
        adjHash[ node1 ] = adjHash.get(node1, []) + [ node2 ]
           
    # for node in sorted(adjHash):
        # print '%s -> %s' % (node, ','.join(adjHash[ node ]))
        
    return adjHash
        
#deBruijnFromKmers(['GAGG','CAGG','GGGG','GGGA','CAGG','AGGG','GGAG'])



def addCycle(cycle, graph):
    startNode = cycle[-1]
    
    while len(graph[ startNode ]) > 0:
        idx = random.randint(0, len(graph[ startNode ])-1)
        startNode = graph[ startNode ].pop(idx)
        cycle.append(startNode)
        
    return cycle

    
def reorderCycle(cycle, graph):
    reorderedCycle = []
    
    for i in range(len(cycle)-1):
        if len(graph[ cycle[i] ]) > 0:
            reorderedCycle = cycle[i:-1] + cycle[:(i+1)]
            break
        
    return reorderedCycle

    
def graphFullyExplored(graph):
    for startNode in graph:
        if len(graph[ startNode ]) > 0:
            return False
    return True

    
def eulerianCycle(graph):
    startNode = random.choice(graph.keys())
    cycle = addCycle([startNode], graph)
    
    while not graphFullyExplored(graph):
        cycle = reorderCycle(cycle, graph)
        cycle = addCycle(cycle, graph)
        
    return cycle


def readAdjacencyFile(infile):
    adjHash = dict()
    
    adjFile = open(infile, 'r')
    for line in adjFile:
        start, stop = line.rstrip().split(' -> ')
        adjHash[ start ] = stop.split(',')
    adjFile.close()
    
    return adjHash
    
# graph = readAdjacencyFile('testEulerian3.txt')
# print '->'.join(eulerianCycle(graph))


def fixUnbalancedNodes(graph):
    countHash = dict()
    numUnbalancedNodes = 0
    
    for startNode in graph:
        countHash[ startNode ] = countHash.get(startNode, 0) + \
                                 len(graph[ startNode ])
        for endNode in graph[ startNode ]:
            countHash[ endNode ] = countHash.get(endNode, 0) - 1
    
    startNode = endNode = None
    for node in countHash:
        if countHash[ node ] != 0:
            numUnbalancedNodes += 1
            if countHash[ node ] > 0:
                startNode = node
            else:
                endNode = node
    
    if numUnbalancedNodes!=2:
        print 'Number of unbalanced nodes != 2 (%d)' % numUnbalancedNodes
        sys.exit()

    graph[ endNode ] = graph.get(endNode, []) + [ startNode ]
    
    return (startNode, endNode)
        

def findEulerianPath(graph):
    startNode, endNode = fixUnbalancedNodes(graph)
    cycle = eulerianCycle(graph)
    
    for i in range(len(cycle)-1):
        if cycle[i]==endNode and cycle[i+1]==startNode:
            cycle = cycle[(i+1):-1] + cycle[:(i+1)]
            break
            
    return cycle
    
#graph = readAdjacencyFile('testEulerianPath3.txt')
#print '->'.join(findEulerianPath(graph))


def deBruijnFromReads(reads):
    graph = deBruijnFromKmers(reads)
    path = findEulerianPath(graph)
    return stringByGenomePath(path)
    
# print deBruijnFromReads(['CTTA','ACCA','TACC','GGCT','GCTT','TTAC'])


def universalCircularStrings(k):
    strings = ['']
    
    while len(strings) > 0 and len(strings[0]) < k:
        strings = extendBinaryString(strings)
        
    graph = deBruijnFromKmers(strings)
    cycle = eulerianCycle(graph)
    return stringByGenomePath(cycle[:(1-k)])
        
#print universalCircularStrings(9)


def pairedComposition(string, k, d):
    retVal = []
    for pos in range(len(string)-2*k-d+1):
        pat1, pat2 = (string[pos:pos+k], string[pos+k+d:pos+2*k+d])
        retVal.append('(%s|%s)' % (pat1, pat2))
        
    return sorted(retVal)
        
        
#print ' '.join(pairedComposition('TAATGCCATGGGATGTT', 3, 2))


def stringSpelledByGappedPatterns(pairedStrings, d):
    prefixes = []
    suffixes = []
    
    for pairedString in pairedStrings:
        read1, read2 = pairedString.rstrip().split('|')
        prefixes.append(read1)
        suffixes.append(read2)
        
    k = len(prefixes[0])
    read1 = stringByGenomePath(prefixes)
    read2 = stringByGenomePath(suffixes)
    
    for pos in range(k+d, len(read1)):
        if read1[pos]!=read2[pos-k-d]:
            return 'Prefix/Suffix strings disagree'
            
    return read1 + read2[(len(read1)-k-d):]
    
#print stringSpelledByGappedPatterns(['GACC|GCGC','ACCG|CGCC','CCGA|GCCG','CGAG|CCGG','GAGC|CGGA'], 2)
        
        
#def stringReconFromReadPairs(readPairs, d):
    
    
    
    
    
    
    
    
# Quiz questions
print deBruijnFromReads(['AAAT','AATG','ACCC','ACGC','ATAC','ATCA','ATGC','CAAA','CACC','CATA','CATC','CCAG','CCCA','CGCT','CTCA','GCAT','GCTC','TACG','TCAC','TCAT','TGCA'])

print not repeated_kMer('1110001011', 3)
print not repeated_kMer('0101001101', 3)
print not repeated_kMer('1001101100', 3)
print not repeated_kMer('1100011011', 3)
print not repeated_kMer('1101000111', 3)
print not repeated_kMer('0100011101', 3)


quizFile = open('quiz7.txt','r')
adjHash = dict()
pairs = []

for line in quizFile:
    pairs.append(line.rstrip())
    read1, read2 = line.rstrip().split('|')
    key = read1[:2] + '|' + read2[:2]
    # val = read1[1:] + '|' + read2[1:]
    # adjHash[ key ] = adjHash.get(key, []) + [ val ]
    adjHash[ key ] = adjHash.get(key, []) + [line.rstrip()]
quizFile.close()

print adjHash
for pair in pairs:
    read1, read2 = pair.split('|')
    key = read1[1:] + '|' + read2[1:]
    print pair, adjHash.get(key, '')
    
    
    