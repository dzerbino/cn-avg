# Copyright (c) 2012, Daniel Zerbino
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# (1) Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer. 
# 
# (2) Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in
# the documentation and/or other materials provided with the
# distribution.  
# 
# (3)The name of the author may not be used to
# endorse or promote products derived from this software without
# specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#!/usr/bin/env python

"""Producing an initial flow history underlying a metagenomic net flow change"""

import random
import copy
import math
from collections import Counter
from cnavg.flows.edge import Edge
from cnavg.flows.cycle import Cycle
from cnavg.flows.flows import Event
from cnavg.history.history import History

import cnavg.avg.graph as avg
import cnavg.avg.balanced as balanced
import cnavg.avg.module as module
import cnavg.history.euclidian as euclidian
import cnavg.history.constrained as constrained

MIN_FLOW=1e-10
MIN_CYCLE_FLOW=1e-2

########################################
## Closing off pseudo-telomeres:
########################################
def extractPseudoTelomereCycle(module, edges):
	firstNode = edges[0].start
	lastNode = edges[-1].finish
	if firstNode == lastNode:
		return module, edges
	else:
		nextNode = module[lastNode].twin
		edge1 = Edge(lastNode, nextNode, edges[0].value, 0)	
		module.removeEdgeFlow(edge1)
		followingNode = module[nextNode].partner
		edge2 = Edge(nextNode, followingNode, -edges[0].value, -1)	
		module.removeEdgeFlow(edge2)
		return extractPseudoTelomereCycle(module, edges + [edge1, edge2])

def closePseudoTelomere(data, index):
	module, history = data
	PT = random.choice(list(module.pseudotelomeres))
	twin = module[PT].twin
	edge1 = Edge(PT, twin, module[PT].segment[index], index)
	module.removeEdgeFlow(edge1)
	edge2 = Edge(twin, module[twin].partner, -edge1.value, -1)
	module.removeEdgeFlow(edge2)
	module, edges = extractPseudoTelomereCycle(module, [edge1, edge2])
	history.absorbEvent(Event(Cycle(edges)))
	return module, history

def closePseudoTelomeres(module, history):
	if len(module.pseudotelomeres) == 0:
		return module, history
	else:
		pseudotelomere = random.choice(list(module.pseudotelomeres))
		segmentCount = len(module[pseudotelomere].segment)
		return reduce(closePseudoTelomere, range(segmentCount), (module, history))

#############################################
## Pre-compute signed node adjacencies
#############################################

nodeEdgesTable = None

def nodePairAdjacency(node, nodeB, module):
	# Conjugate flow!
	return (node, nodeB, -module[node].edges[nodeB], -1)

def nodePairSegment(node, index, module):
	return (node, module[node].twin, module[node].segment[index], index)

def nodeEdges(node, module):
	adjacencies = [nodePairAdjacency(node, X, module) for X in module[node].edges]
	segments = [nodePairSegment(node, X, module) for X in range(len(module[node].segment))]
	return adjacencies + segments

def computeNodeEdges(module):
	return dict((X, nodeEdges(X, module)) for X in module.nodes())

#############################################
## Search for small edge
#############################################

def positiveNeighbourhood(node):
	return filter(lambda X: X[2] > MIN_FLOW, nodeEdgesTable[node])

def negativeNeighbourhood(node):
	return filter(lambda X: X[2] < -MIN_FLOW, nodeEdgesTable[node])

def phasedNeighbourhood(node, value):
	if value > 0:
		return [X[1] for X in positiveNeighbourhood(node)] 
	else:
		return [X[1] for X in negativeNeighbourhood(node)] 

def oppositeNeighbourhood(node, value):
	if value > 0:
		return [X[1] for X in negativeNeighbourhood(node)] 
	else:
		return [X[1] for X in positiveNeighbourhood(node)] 

def adjacencies():
	return sum(nodeEdgesTable.values(), [])

def minimumEdge(module):
	edges = filter(lambda X: X[0] != X[1] and abs(X[2]) > MIN_FLOW, adjacencies())
	if len(edges) == 0:
		return None
	else:
		res = min(edges, key=lambda X: abs(X[2]))
		return Edge(res[0], res[1], res[2], res[3])

#############################################
## Dijkstra
#############################################

def dijkstra(origin, value, graph, blockTwin):
    distances = dict(((node, [-1,-1]) for node in graph))
    todo = [origin] 

    distances[origin][0] = 0

    while len(todo) > 0: 
	node = todo.pop(0)
        newdist = distances[node][0] + 1 

	for node2 in oppositeNeighbourhood(node, value):
	    if blockTwin and node == origin and node2 == graph[origin].twin:
		continue
	    if distances[node2][1] > -1:
		continue
	    distances[node2][1] = newdist;

	    for node3 in phasedNeighbourhood(node2, value):
		if distances[node3][0] > -1:
		    continue
		else:
		    distances[node3][0] = newdist
		    todo.append(node3)

    return distances

#############################################
## Heuristic propagation
#############################################

def signedEdges(node, module, sign):
	return filter(lambda X: X[2] * sign > MIN_FLOW, nodeEdgesTable[node])

def nodeDistances(node, module, distances, sign, phase):
	if phase:
		return map(lambda X: distances[X[1]][1], signedEdges(node, module, sign))
	else:
		return map(lambda X: distances[X[1]][0], signedEdges(node, module, sign))

def minDist(node, module, distances, sign, phase):
	candidates = filter(lambda X: X >= 0, nodeDistances(node, module, distances, sign, phase))
	if len(candidates) == 0:
		return None
	return min(candidates)

def nextNodes(node, module, distances, sign, phase):
	dist = minDist(node, module, distances, sign, phase)
	if dist is None:
		return None
	edges = signedEdges(node, module, sign)
	if phase:
		return filter(lambda X: distances[X[1]][1] == dist, edges)
	else:
		return filter(lambda X: distances[X[1]][0] == dist, edges)

def chooseNextNode(node, module, distances, sign, phase):
	vals = nextNodes(node, module, distances, sign, phase)
	if vals is None:
		return None
	else:
		return random.choice(vals)

def signf(val):
	if val > 0:
		return 1
	else:
		return -1

def extendCycle(edgeList, module, distances, sign):
	edgeData = chooseNextNode(edgeList[-1].finish, module, distances, signf(edgeList[0].value * sign), sign > 0)
	if edgeData is None:
		return edgeList, module, False
	edge = Edge(edgeData[0], edgeData[1], edgeData[2], edgeData[3])
	edge.value = edgeList[0].value * sign
	edgeList.append(edge)

	if len(edgeList) % 2 == 0 and edge.finish == edgeList[0].start:	
		return edgeList, module, True
	else:
		return extendCycle(edgeList, module, distances, -sign)

def getAbsEdgeValue(adjacencyIndex, module):
	start, end, index = adjacencyIndex
	if index < 0:
		return abs(module[start].edges[end])
	else:
		return abs(module[start].segment[index])

def extractCycle(edge, module):
	distances = dijkstra(edge.start, edge.value, module, blockTwin=(edge.index >= 0))
	edgeList, module, success = extendCycle([edge], module, distances, -1)
	if success:
		edge_counts = Counter(E.adjacencyIndex() for E in edgeList)
		if max(edge_counts.values()) == 1:
			cycle = Cycle(edgeList)
		else:
			value = min(float(getAbsEdgeValue(E, module)) / edge_counts[E] for E in edge_counts)
			cycle = Cycle(edgeList, value=math.copysign(value, edge.value))

		return Event(cycle)
	else:
		return None


def pickOutCycle(module):
	global nodeEdgesTable
	nodeEdgesTable = computeNodeEdges(module)
	edge = minimumEdge(module) 
	if edge is None:
		# Job finished
		return None
	else:
		return extractCycle(edge, module)

def pickOutCycles(module, history):
	while True:
		event = pickOutCycle(module)
		if event is not None:
			if len(history.events) % 100 == 0:
				print 'CYCLE', len(history.events)
		        map(module.removeEdgeFlow, event.cycle)
			history.absorbEvent(event)
		else:
			return history

########################################
## Filter out low value cycles
########################################

def highFlowHistory(history, cactusHistory, net):
	res = euclidian.EuclidianHistory(copy.copy(history.module))
	res.module.reset()
	cactusHistory.update(net, res)
	events = list(sorted(history.events, key=lambda X: -X.ratio))
	for event in events:
		if event.ratio > MIN_CYCLE_FLOW:
			cactusHistory.absorbEvent(res, event)
			res.module.addEventFlow(event)
	cactusHistory.updateCNVs(net, res)
	return res

########################################
## Finding an initial net history
########################################

def seedHistory(cactusHistory, net, cnvs):
	M = module.Module(net, cactusHistory.cactus, cnvs)
	H = History(M)
	MC = copy.copy(M)
	MC, H1 = closePseudoTelomeres(MC, H)
	H2 = pickOutCycles(MC, H1)
	return highFlowHistory(H2, cactusHistory, net)

########################################
## Finding an initial cactus graph history
########################################

def propagateInitialHistory_Chain(chain, history, cnvs):
	return reduce(lambda X, Y: propagateInitialHistory_Net(Y, X, cnvs), history.cactus.chains2Nets[chain], history) 

def propagateInitialHistory_Net(net, history, cnvs):
	seedHistory(history, net, cnvs)
	return reduce(lambda X,Y: propagateInitialHistory_Chain(Y, X, history.chainCNVs[Y]), history.cactus.nets2Chains[net], history)

###############################################
## Master function
###############################################

def initialHistory(cactus):
	print "Extracting initial history from Cactus"
	return propagateInitialHistory_Net(cactus.rootNet, constrained.ConstrainedHistory(cactus), [])

###############################################
## Unit test
###############################################

def main():
	G = avg.randomNearEulerianGraph(10)
	B = balanced.BalancedAVG(G)
	C = cactus.Cactus(B)
	N = normalized.NormalizedCactus(C)
	O = oriented.OrientedCactus(N)
	H = initialHistory(O)
	print H
	H.validate()
	print H.rearrangementCost()

if __name__ == "__main__":
	import cnavg.cactus.graph as cactus
	import cnavg.cactus.oriented as oriented
	import cnavg.cactusSampling.sampling as normalized
	main()
