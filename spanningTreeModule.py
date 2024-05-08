import igraph as ig
import numpy as np
import itertools

def generateSpanningTreeDict(graph):
    #HIGHLY inefficient in general, but good for sparse and small graphs like the Petersen Graph
    #basic theory is to randomly select all possible arrangements of the appropriate number of edges
    #and select only those ones which are connected.  Otherwise they're not a spanning tree.
    spanningTreeSet = set()
    edgeList = graph.to_tuple_list()
    vertexCount = graph.vcount()
    for edgeSubset in itertools.combinations(edgeList, vertexCount - 1):
        if(isSpanningTree(edgeSubset, vertexCount)):
            spanningTreeSet.add(edgeSubset)
    spanningTreeDict = dict()
    index = 0
    for tree in spanningTreeSet:
        treeGraph = ig.Graph.TupleList(tree, directed = False, weights = False)
        spanningTreeDict[index] = [tree, treeGraph]
        index += 1
    return spanningTreeDict

def treeSetToMultipleMetaGraphs(spanningTreeDict, chromaticNumber, graph):
    #generate a graph, with each tree as a node, and edges between trees if they're "adjacent"
    #that is to say, if they differ by exactly one edge
    #This is done by constructing an appropriate dict of lists and using the igraph method to make a graph from that
    metaGraphDictComplete = dict.fromkeys(spanningTreeDict.keys())
    metaGraphDictGreedy = dict.fromkeys(spanningTreeDict.keys())
    errorDict = dict.fromkeys(spanningTreeDict.keys())
    for key in errorDict.keys():
        keyError = incorrectEdges(spanningTreeDict[key], chromaticNumber, graph)
        errorDict[key] = keyError
    for key in spanningTreeDict.keys():
        metaGraphDictComplete[key] = []
        metaGraphDictGreedy[key] = []
        for iteratedKey in spanningTreeDict.keys():
            if len(set(spanningTreeDict[key][0]).difference(set(spanningTreeDict[iteratedKey][0]))) == 1: 
            #this equals 2 if and only if there is exactly one edge in each not in the other
                metaGraphDictComplete[key].append(iteratedKey)
                if errorDict[iteratedKey] < errorDict[key]:
                    metaGraphDictGreedy[key].append(iteratedKey)
    metaGraphComplete = ig.Graph.ListDict(metaGraphDictComplete, directed=False)
    metaGraphGreedy = ig.Graph.ListDict(metaGraphDictGreedy, directed=True)
    return metaGraphComplete, metaGraphGreedy, errorDict

def isSpanningTree(edgeList, vertexCount):
    tempGraph = ig.Graph.TupleList(edgeList)
    return (tempGraph.ecount() == vertexCount - 1 and tempGraph.vcount() == vertexCount and tempGraph.is_connected())

def incorrectEdges(spanningTuple, chromaticNumber, graph):
    incorrectEdges = 0
    spanningTreeGraph = spanningTuple[1]
    vertex = 0  #arbitrarily defined to be the tree's root vertex, but always the same vertex of the underlying graph
    for offset in range(chromaticNumber):
        colorClass = []
        for factor in range(-(spanningTreeGraph.diameter()//-chromaticNumber)): #I want to ceiling the result here so double-negation is used
            factorColorClassIndices = spanningTreeGraph.neighborhood(vertex, order = chromaticNumber * factor + offset, mindist = chromaticNumber * factor + offset)
            factorColorClassNames = [spanningTreeGraph.vs[index]['name'] for index in factorColorClassIndices]
            colorClass = colorClass + factorColorClassNames
        for vertexPair in itertools.combinations(colorClass, 2):
            if vertexPair[0] in graph.neighbors(vertexPair[1]):
                    incorrectEdges += 1
    return incorrectEdges