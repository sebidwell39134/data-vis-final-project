import igraph as ig
import numpy as np
import itertools

def generateSpanningTreeDict(graph):
    #HIGHLY inefficient in general, but good for sparse and small graphs like the Petersen Graph
    #basic theory is to randomly select all possible arrangements of the appropriate number of edges
    #and select only those ones which are connected.  Otherwise they're not a spanning tree.
    spanningTreeSet = {}
    edgeList = graph.to_tuple_list()
    vertexCount = graph.vcount()
    for edgeSubset in itertools.combinations(edgeList, vertexCount - 1):
        if(isSpanningTree(edgeSubset, vertexCount)):
            spanningTreeSet.add(edgeSubset)
    spanningTreeDict = dict()
    index = 0
    for tree in spanningTreeSet:
        index += 1
        treeGraph = ig.TupleList(tree, directed = False, weights = False)
        spanningTreeDict.append({index:[tree, treeGraph]})
    return spanningTreeDict

def treeSetToMultipleMetaGraphs(spanningTreeDict, chromaticNumber, graph):
    #generate a graph, with each tree as a node, and edges between trees if they're "adjacent"
    #that is to say, if they differ by exactly one edge
    #This is done by constructing an appropriate dict of lists and using the igraph method to make a graph from that
    metaGraphDictComplete = dict.fromkeys(spanningTreeDict.keys)
    metaGraphDictGreedy = dict.fromkeys(spanningTreeDict.keys)
    errorDict = dict.fromkeys(spanningTreeDict.keys)
    for key in errorDict.keys:
        keyError = incorrectEdges(spanningTreeDict[key], chromaticNumber, graph)
        errorDict[key].append(keyError)
    for key in spanningTreeDict.keys:
        metaGraphDictComplete[key] = []
        metaGraphDictGreedy[key] = []
        for iteratedKey in spanningTreeDict.keys:
            if len(set(spanningTreeDict[key][0]).difference(set(spanningTreeDict[iteratedKey][0]))) == 2: 
            #this equals 2 if and only if there is exactly one edge in each not in the other
                metaGraphDictComplete[key].append(iteratedKey)
                if errorDict[iteratedKey] < errorDict[key]:
                    metaGraphDictGreedy[key].append(iteratedKey)
    metaGraphComplete = ig.ListDict(metaGraphDictComplete, directed=True)
    metaGraphGreedy = ig.ListDict(metaGraphDictGreedy, directed=True)
    return metaGraphComplete, metaGraphGreedy, errorDict

def isSpanningTree(edgeList, vertexCount):
    tempGraph = ig.TupleList(edgeList)
    if tempGraph.ecount() == vertexCount - 1 and tempGraph.vcount() == vertexCount and tempGraph.isConnected():
        return True
    else:
        return False

def incorrectEdges(spanningTuple, chromaticNumber, graph):
    incorrectEdges = 0
    spanningTreeGraph = spanningTuple[1]
    vertex = graph.vs[0]  #arbitrarily defined to be the tree's root vertex, but always the same vertex of the underlying graph
    for offset in range(chromaticNumber):
        colorClass = []
        for factor in range(spanningTreeGraph.diameter()//chromaticNumber): #I want to floor the result here
            colorClass.append(spanningTreeGraph.neighborhood(vertex, order = chromaticNumber * (factor + 1) + offset, mindist = chromaticNumber * (factor + 1) + offset))
        for vertexPair in itertools.combinations(colorClass, 2):
            if vertexPair[0] in graph.neighbors(vertexPair[1]):
                incorrectEdges += 1
    return incorrectEdges