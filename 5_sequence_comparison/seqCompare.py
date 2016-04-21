#!/usr/bin/env python

import re, sys

sys.path.append('..')
#from utilityFunctions import readAdjacencyFile

DOWN, RIGHT, DIAG, FREE_RIDE = 1, 2, 3, 4

def dp_change(money, coins):
    minNumCoins = [0]*(money+1)
    
    for m in range(1, len(minNumCoins)):
        minNumCoins[m] = sys.maxint
        
        for i in range(len(coins)):
            if m >= coins[i]:
                if minNumCoins[m-coins[i]]+1 < minNumCoins[m]:
                    minNumCoins[m] = minNumCoins[m-coins[i]]+1
    
    return minNumCoins.pop()
    
#print dp_change(40, [50,25,20,10,5,1]) #2
#print dp_change(17948, [20,14,8,5,3,1]) #898


def readManhattanFile(file):
    right, down = [], []
    
    infile = open(file, 'r')
    vals = infile.readline().rstrip().split()
    n, m = [ int(x) for x in vals ]
    
    for i in range(n):
        vals = infile.readline().rstrip().split()
        down.append([ int(x) for x in vals ])
        
    infile.readline()
    
    for i in range(n+1):
        vals = infile.readline().rstrip().split()
        right.append([ int(x) for x in vals ])
        
    return n, m, down, right
    
    
def manhattanDist(n, m, down, right):
    s = [ [0]*(m+1) for x in range(n+1) ]
    
    for i in range(n):
        s[i+1][0] = s[i][0] + down[i][0]
    for j in range(m):
        s[0][j+1] = s[0][j] + right[0][j]
    for i in range(n):
        for j in range(m):
            s[i+1][j+1] = max(s[i][j+1] + down[i][j+1],
                              s[i+1][j] + right[i+1][j])
                          
    return s[n][m]

# n, m, down, right = readManhattanFile('manhattan3.txt')
# print manhattanDist(n, m, down, right)

def LCS_backtrack(v, w):
    s = [ [0]*(len(w)+1) for x in range(len(v)+1) ]
    backtrack = [ [0]*(len(w)+1) for x in range(len(v)+1) ]
    
    for i in range(len(v)):
        for j in range(len(w)):
            s[i+1][j+1] = max(s[i][j+1], s[i+1][j],
                              (s[i][j] + 1) if v[i]==w[j] else 0)
            if s[i+1][j+1]==s[i][j+1]:
                backtrack[i+1][j+1] = DOWN
            if s[i+1][j+1]==s[i+1][j]:
                backtrack[i+1][j+1] = RIGHT
            if s[i+1][j+1]==s[i][j] + 1 and v[i]==w[j]:
                backtrack[i+1][j+1] = DIAG
                
    #print s
    return backtrack
    
#print  LCS_backtrack('ACG','TCT')

# output longest common sequence with alignments (optional)
def outputLCS(backtrack, v, w, returnAlignments=False):
    LCS = ''
    vAlign = str()
    wAlign = str()
    i, j = len(v), len(w)
    
    while i>0 or j>0:
        #print "i:%d j:%d backtrack[i][j]:%d" % (i, j, backtrack[i][j])
        if backtrack[i][j]==DOWN or j==0:
            i = i-1
            vAlign = v[i] + vAlign
            wAlign = '-' + wAlign
        elif backtrack[i][j]==RIGHT or i==0:
            j = j-1
            vAlign = '-' + vAlign
            wAlign = w[j] + wAlign
        elif backtrack[i][j]==DIAG:
            i, j = i-1, j-1
            LCS = v[i] + LCS
            vAlign = v[i] + vAlign
            wAlign = w[j] + wAlign
    
    if returnAlignments:
        return LCS, vAlign, wAlign
    return LCS

# v, w = 'AACCTTGG','ACACTGTGA'
# backtrack = LCS_backtrack(v,w)
# print backtrack
# print outputLCS(backtrack, v, w)


# return the topological ordering of the passed graph
def topologicalOrdering(graph, finalNode=None):
    topoList = []
    startNodes = set(graph.keys())
    backtrack = dict()
    
    # create dict of child nodes to parents
    # and remove nodes with parents from startNodes
    for startNode in graph:
        for endNode in graph[ startNode ]:
            backtrack[ endNode ] = backtrack.get(endNode, [])
            backtrack[ endNode ].append(startNode)
            if endNode in startNodes:
                startNodes.remove(endNode)
    
    # generate topological ordering
    while len(startNodes) > 0:
        a = startNodes.pop()
        topoList.append(a)
        
        if finalNode is not None and a==finalNode:
            break
        
        if a in graph:
            for endNode in graph[a]:
                backtrack[endNode].remove(a)
                if len(backtrack[endNode])==0:
                    startNodes.add(endNode)
                
    return topoList

# topo = topologicalOrdering(readAdjacencyFile('topo2.txt'))
# print ', '.join(topo)


# read the weighted directional acyclic graph file and trim
# all extraneous nodes that aren't directly downstream of the
# starting node
def readDagFile(infile):
    # read graph information from file and create backtrack dict
    input = open(infile, 'r')
    startNode = input.readline().rstrip()
    endNode = input.readline().rstrip()
    
    edgeWeights = dict()
    parentNodes, childNodes = [],[]
    
    # read parentNodes, childNodes, weights into lists/dicts
    for line in input:
        parent, child, weight = re.split('\D+', line.rstrip())
        edgeWeights[ parent + ":" + child ] = int(weight)
        parentNodes.append(parent)
        childNodes.append(child)
    input.close()
    
    # prune the tree of nodes not downstream of startNode
    while True:
        parentSet = set(parentNodes)
        childSet = set(childNodes)
        goodNodes = []
        for i in range(len(parentNodes)):
            if parentNodes[i]==startNode or parentNodes[i] in childNodes:
                goodNodes.append(i)
            else:
                print 'pruning %s->%s' % (parentNodes[i], childNodes[i])
        if len(goodNodes)==len(parentNodes):
            break
        parentNodes = [ parentNodes[x] for x in goodNodes ]
        childNodes = [ childNodes[x] for x in goodNodes ]
    
    # generate forward and backtrack graphs
    graph, backtrack = dict(), dict()
    for i in range(len(parentNodes)):
        graph[ parentNodes[i] ] = graph.get(parentNodes[i], [])
        graph[ parentNodes[i] ].append(childNodes[i])
        backtrack[ childNodes[i] ] = backtrack.get(childNodes[i], [])
        backtrack[ childNodes[i] ].append(parentNodes[i])
        
    return startNode, endNode, graph, backtrack, edgeWeights
    
    
    
# determine the path through a weighted directional acyclic graph
def getDagPath(infile):
    startNode, endNode, graph, backtrack, edgeWeights = readDagFile(infile)
    
    # calculate cumulative scores per node
    topoOrder = topologicalOrdering(graph, finalNode=endNode)
    nodeScores = dict()
    for node in topoOrder:
        if node in backtrack:
            scores = []
            for parentNode in backtrack[ node ]:
                weight = edgeWeights[ parentNode + ':' + node ]
                scores.append(nodeScores[parentNode] + weight)
            score = max(scores)
        else:
            score = 0
        nodeScores[ node ] = score
        
    # print topoOrder
    # print graph
    # print backtrack
    # print nodeScores
    
    # backtrack
    nodes = [endNode]
    node = endNode
    while node!=startNode:
        bestNode, maxScore = '', -sys.maxint
        for parentNode in backtrack[ node ]:
            score = nodeScores[parentNode] + edgeWeights[ parentNode + ':' + node ]
            print node, nodeScores[node], parentNode, nodeScores[parentNode], edgeWeights[ parentNode + ':' + node ], score
            if score > maxScore:
               bestNode, maxScore = parentNode, score
        nodes.append(bestNode)
        node = bestNode
        
    print nodeScores[ endNode ]
    print '->'.join(reversed(nodes))

#getDagPath('DAG3.txt')


# reads substitution matrix (e.g. BLOSUM62) from file
# and returns dictionary of substitution penalties
def readSubMatrix(filePath):
    scoreMat = dict()
    
    infile = open(filePath, 'r')
    items = infile.readline().strip().split()
    for i, line in enumerate(infile):
        tokens = line.strip().split()
        if tokens[0]!= items[i]:
            print 'items in %s at line %d do not agree: %s %s' % \
                  (filePath, i, tokens[0], items[i])
            sys.exit(1)
        for j in range(len(items)):
            scoreMat[ tokens[0] + ':' + items[j] ] = int(tokens[j+1])
    infile.close()
    
    return scoreMat
    
    
def matrix_scored_alignment(v, w, scoreMat, sigma=5, local=False):
    s = [ [0]*(len(w)+1) for x in range(len(v)+1) ]
    backtrack = [ [0]*(len(w)+1) for x in range(len(v)+1) ]
    
    for i in range(len(v)):
        s[i+1][0] = s[i][0] - sigma
        
        for j in range(len(w)):
            s[0][j+1] = s[0][j] - sigma
            
            mu = scoreMat[ v[i] + ':' + w[j] ]
            scores = [ s[i][j+1] - sigma, 
                       s[i+1][j] - sigma,
                       s[i][j] + mu ]
            if local:
                scores.append(0)
            s[i+1][j+1] = max(scores)
            
            if local and s[i+1][j+1]==0:
                backtrack[i+1][j+1] = FREE_RIDE
            elif s[i+1][j+1]==s[i][j+1] - sigma:
                backtrack[i+1][j+1] = DOWN
            elif s[i+1][j+1]==s[i+1][j] - sigma:
                backtrack[i+1][j+1] = RIGHT
            elif s[i+1][j+1]==s[i][j] + mu:
                backtrack[i+1][j+1] = DIAG
            
    return s, backtrack

# output longest common sequence and alignments (optional)
# for the local alignment case
def outputLCS_local(backtrack, s, v, w, returnAlignments=False):
    LCS = ''
    vAlign = ''
    wAlign = ''
    
    max = -sys.maxint
    i, j = len(s), len(s[0])
    for x in range(len(s)):
        for y in range(len(s[x])):
            if s[x][y] > max:
                max = s[x][y]
                i, j = x, y
    
    while i>0 or j>0:
        #print "i:%d j:%d backtrack[i][j]:%d" % (i, j, backtrack[i][j])
        if backtrack[i][j]==FREE_RIDE:
            break
        elif backtrack[i][j]==DOWN or j==0:
            i = i-1
            vAlign = v[i] + vAlign
            wAlign = '-' + wAlign
        elif backtrack[i][j]==RIGHT or i==0:
            j = j-1
            vAlign = '-' + vAlign
            wAlign = w[j] + wAlign
        elif backtrack[i][j]==DIAG:
            i, j = i-1, j-1
            LCS = v[i] + LCS
            vAlign = v[i] + vAlign
            wAlign = w[j] + wAlign
    
    if returnAlignments:
        return LCS, max, vAlign, wAlign
    return LCS

    
# blosum62 = readSubMatrix('BLOSUM62.txt')
# v, w = 'PLEASANTLY','MEANLY'
# s, backtrack = matrix_scored_alignment(v, w, blosum62)
# #print s, backtrack
# print s[len(v)][len(w)]
# LCS, vAlign, wAlign = outputLCS(backtrack, v, w, returnAlignments=True)
# print vAlign
# print wAlign


# pam250 = readSubMatrix('PAM250_1.txt')
# v, w = 'MEANLY','PENALTY'
# s, backtrack = matrix_scored_alignment(v, w, pam250, local=True)
# #print s, backtrack
# LCS, max, vAlign, wAlign = outputLCS_local(backtrack, s, v, w, returnAlignments=True)
# print '%s\n%s\n%s\n' % (max, vAlign, wAlign)


# editMat = readSubMatrix('editDistanceMatrix.txt')
# v, w = 'PLEASANTLY','MEANLY'
# s, backtrack = matrix_scored_alignment(v, w, editMat, sigma=1)
# LCS, vAlign, wAlign = outputLCS(backtrack, v, w, returnAlignments=True)
# print s[len(v)][len(w)]
# print vAlign
# print wAlign

def outputLCS_fitted(backtrack, s, v, w, returnAlignments=False):
    LCS, vAlign, wAlign = '', '', ''
    i, j = len(v), len(w)
    max = -sys.maxint
    for x in range(len(v)):
        if s[x][len(w)] > max:
            i, max = x, s[x][len(w)]
    
    while j>0:
        #print "i:%d j:%d backtrack[i][j]:%d" % (i, j, backtrack[i][j])
        if backtrack[i][j]==DOWN:
            i = i-1
            vAlign = v[i] + vAlign
            wAlign = '-' + wAlign
        elif backtrack[i][j]==RIGHT or i==0:
            j = j-1
            vAlign = '-' + vAlign
            wAlign = w[j] + wAlign
        elif backtrack[i][j]==DIAG:
            i, j = i-1, j-1
            LCS = v[i] + LCS
            vAlign = v[i] + vAlign
            wAlign = w[j] + wAlign
    
    if returnAlignments:
        return LCS, max, vAlign, wAlign
    return LCS

    
def fitted_alignment(v, w, scoreMat, sigma=5):
    s = [ [0]*(len(w)+1) for x in range(len(v)+1) ]
    backtrack = [ [0]*(len(w)+1) for x in range(len(v)+1) ]
    
    for i in range(len(v)):
        for j in range(len(w)):
            s[0][j+1] = -sys.maxint
            s[len(v)][j] = -sys.maxint
            
            mu = scoreMat[ v[i] + ':' + w[j] ]
            scores = [ s[i][j+1] - sigma, 
                       s[i+1][j] - sigma,
                       s[i][j] + mu ]
            s[i+1][j+1] = max(scores)
            
            if s[i+1][j+1]==s[i][j+1] - sigma:
                backtrack[i+1][j+1] = DOWN
            elif s[i+1][j+1]==s[i+1][j] - sigma:
                backtrack[i+1][j+1] = RIGHT
            elif s[i+1][j+1]==s[i][j] + mu:
                backtrack[i+1][j+1] = DIAG
            
    return s, backtrack
    
# dnaSubst = readSubMatrix('dnaSubstMat.txt')
# v, w = 'GTAGGCTTAAGGTTA','TAGATA'
# s, backtrack = fitted_alignment(v, w, dnaSubst, sigma=1)
# #print s, backtrack
# LCS, max, vAlign, wAlign = outputLCS_fitted(backtrack, s, v, w, returnAlignments=True)
# print '%s\n%s\n%s\n' % (max, vAlign, wAlign)


def outputLCS_overlap(backtrack, s, v, w, returnAlignments=False):
    LCS, vAlign, wAlign = '', '', ''
    i, j = len(v), len(w)
    
    max = -sys.maxint
    for col in range(len(w)):
        if s[len(v)][col] >= max:
            j, max = col, s[len(v)][col]
            
    print 'i:%d j:%d' % (i, j)
    while j>0:
        #print "i:%d j:%d backtrack[i][j]:%d" % (i, j, backtrack[i][j])
        if backtrack[i][j]==DOWN:
            i = i-1
            vAlign = v[i] + vAlign
            wAlign = '-' + wAlign
        elif backtrack[i][j]==RIGHT or i==0:
            j = j-1
            vAlign = '-' + vAlign
            wAlign = w[j] + wAlign
        elif backtrack[i][j]==DIAG:
            i, j = i-1, j-1
            LCS = v[i] + LCS
            vAlign = v[i] + vAlign
            wAlign = w[j] + wAlign

    if returnAlignments:
        return LCS, max- s[i][j], vAlign, wAlign
    return LCS

    
def overlap_alignment(v, w, sigma=2, mu=2):
    s = [ [0]*(len(w)+1) for x in range(len(v)+1) ]
    backtrack = [ [0]*(len(w)+1) for x in range(len(v)+1) ]
    
    for i in range(len(v)):
        for j in range(len(w)):
            scores = [ s[i][j+1] - sigma, 
                       s[i+1][j] - sigma,
                       s[i][j] + (1 if v[i]==w[j] else -mu) ]
            s[i+1][j+1] = max(scores)
            
            if s[i+1][j+1]==s[i][j+1] - sigma:
                backtrack[i+1][j+1] = DOWN
            elif s[i+1][j+1]==s[i+1][j] - sigma:
                backtrack[i+1][j+1] = RIGHT
            elif s[i+1][j+1]==s[i][j] + (1 if v[i]==w[j] else -mu):
                backtrack[i+1][j+1] = DIAG
            
    return s, backtrack


# v, w = 'PAWHEAE','HEAGAWGHEE'
# s, backtrack = overlap_alignment(v, w, sigma=2)
# print s, backtrack
# LCS, max, vAlign, wAlign = outputLCS_overlap(backtrack, s, v, w, returnAlignments=True)
# print '%s\n%s\n%s\n' % (max, vAlign, wAlign)


def affine_alignment(v, w, scoreMat, sigma=2, epsilong=1, mu=2):
    s = [ [0]*(len(w)+1) for x in range(len(v)+1) ]
    backtrack = [ [0]*(len(w)+1) for x in range(len(v)+1) ]
    
    for i in range(len(v)):
        for j in range(len(w)):
            scores = [ s[i][j+1] - sigma, 
                       s[i+1][j] - sigma,
                       s[i][j] + (1 if v[i]==w[j] else -mu) ]
            s[i+1][j+1] = max(scores)
            
            if s[i+1][j+1]==s[i][j+1] - sigma:
                backtrack[i+1][j+1] = DOWN
            elif s[i+1][j+1]==s[i+1][j] - sigma:
                backtrack[i+1][j+1] = RIGHT
            elif s[i+1][j+1]==s[i][j] + (1 if v[i]==w[j] else -mu):
                backtrack[i+1][j+1] = DIAG
            
    return s, backtrack
    
blosum62 = readSubMatrix('BLOSUM62.txt')
v, w = 'PRTEINS','PRTWPSEIN'
s, backtrack = affine_alignment(v, w, blosum62)
print '%s\n%s\n' % (s, backtrack)
print s[len(v)][len(w)]
LCS, vAlign, wAlign = outputLCS(backtrack, v, w, returnAlignments=True)
print vAlign
print wAlign