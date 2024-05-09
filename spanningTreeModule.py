import igraph as ig
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

def treeSetToMetaGraph(spanningTreeDict, chromaticNumber, graph):
    #Generate a graph, with each tree as a node, and edges between trees if they're "adjacent"
    #That is to say, if they differ by exactly one edge
    #Further, the edge is only made if the next tree is "better" at coloring the original graph
    metaGraphDictGreedy = dict.fromkeys(spanningTreeDict.keys())
    errorDict = dict.fromkeys(spanningTreeDict.keys())
    for key in errorDict.keys():
        keyError = incorrectEdges(spanningTreeDict[key], chromaticNumber, graph)
        errorDict[key] = keyError
    for key in spanningTreeDict.keys():
        metaGraphDictGreedy[key] = []
        for iteratedKey in spanningTreeDict.keys():
            if len(set(spanningTreeDict[key][0]).difference(set(spanningTreeDict[iteratedKey][0]))) == 1 and errorDict[iteratedKey] < errorDict[key]: 
                metaGraphDictGreedy[key].append(iteratedKey)
    metaGraphGreedy = ig.Graph.ListDict(metaGraphDictGreedy, directed=True)
    return metaGraphGreedy, errorDict

def isSpanningTree(edgeList, vertexCount):
    #A simple check but nicer as a function
    tempGraph = ig.Graph.TupleList(edgeList)
    return (tempGraph.ecount() == vertexCount - 1 and tempGraph.vcount() == vertexCount and tempGraph.is_connected())

def incorrectEdges(spanningTuple, chromaticNumber, graph):
    #Apply the flooding algorithm and find the number of edges in the original graph between vertices of the same color
    incorrectEdges = 0
    spanningTreeGraph = spanningTuple[1]
    vertex = 0  #the tree's root vertex, which will always be the same vertex of the underlying graph
    for offset in range(chromaticNumber):
        colorClass = []
        for factor in range(spanningTreeGraph.diameter() + 1):
            distance = chromaticNumber * factor + offset
            factorColorClassIndices = spanningTreeGraph.neighborhood(vertex, order = distance, mindist = distance)
            factorColorClassNames = [spanningTreeGraph.vs[index]['name'] for index in factorColorClassIndices]
            colorClass = colorClass + factorColorClassNames
        for vertexPair in itertools.combinations(colorClass, 2):
            if vertexPair[0] in graph.neighbors(vertexPair[1]):
                incorrectEdges += 1
    return incorrectEdges