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

import sys
import cnavg.avg.graph as avg
import cnavg.avg.gabp.gabp as gabp
import numpy as np
import graph as cactus
import cnavg.cactusSampling.sampling as normalized

TOL = 1e2
FUDGE_FACTOR = 1
SEG_FACTOR = FUDGE_FACTOR / 10
EDGE_FACTOR = FUDGE_FACTOR * 10
MAX_NOISE = 1e3
MIN_NOISE = 1e-3

"""Balancing of the flows in a Cactus graph"""

##############################################
## Mapping of flow elements onto integer IDs 
##############################################

class Mapping(object):
	def __init__(self):
		self.chains = dict()
		self.edges = dict()
		self.size = 0

	def addChain(self, A):
		self.chains[A] = self.size
		self.size += 1

	def getChain(self, A):
		return self.chains[A]

	def addEdge(self, A, B):
		m = min([A, B])
		M = max([A, B])
		self.edges[(m,M)] = self.size
		self.size += 1

	def getEdge(self, A, B):
		m = min([A, B])
		M = max([A, B])
		return self.edges[(m, M)]

	def chainStr(self, A):
		return "\t".join([str(A), str(self.chains[A])])

	def chainsStr(self):
		return "\n".join(map(lambda X: self.chainStr(X), self.chains))

	def edgeStr(self, pair):
		return "\t".join([str(pair[0]), str(pair[1]), str(self.edges[(pair[0],pair[1])])])

	def edgesStr(self):
		return "\n".join(map(lambda X: self.edgeStr(X), self.edges))

	def __str__(self):
		return "\n".join([self.edgesStr(), self.chainsStr()])

def prepareEdgeMapping (mapping, tonode, fromnode):
	if fromnode <= tonode:
		mapping.addEdge(tonode, fromnode)
	return mapping

def prepareNodeMapping(mapping, node, graph):
	return reduce(lambda M, N: prepareEdgeMapping(M, N, node), graph[node].edges.keys(), mapping)

def prepareGraphNodesMapping(graph, mapping):
	return reduce(lambda M,X: prepareNodeMapping(M,X, graph), graph.nodes(), mapping)

def prepareChainMapping(mapping, chain):
	if chain not in mapping.chains:
		mapping.addChain(chain)
	return mapping

def prepareGraphChainsMapping(graph, mapping):
	return reduce(prepareChainMapping, graph.chains, mapping)

def prepareGraphMapping(graph):
	return prepareGraphChainsMapping(graph, prepareGraphNodesMapping(graph, Mapping()))

##############################################
## Adding in estimates
##############################################

def addEdgeEstimate(estimates, A, B, graph, mapping):
	estimates[mapping.getEdge(A, B)] = graph[A].edges[B]
	return estimates

def initialNodeEstimates(estimates, node, graph, mapping):
	return reduce(lambda E,N: addEdgeEstimate(E, node, N, graph, mapping), graph[node].edges.keys(), estimates)

def initialNodesEstimates(graph, mapping, estimates):
	return reduce(lambda E,N: initialNodeEstimates(E, N, graph, mapping), graph.nodes(), estimates)

def meanChainCoverage(graph, chain, index):
	return sum(map(lambda X: X.copynumber(graph, index) * X.length(), chain)) / chain.length()

def meanTotalChainCoverage(graph, chain):
	return sum(meanChainCoverage(graph, chain, index) for index in range(chain[0].ploidy(graph)))

def addChainEstimate(estimates, chain, graph, mapping):
	estimates[mapping.getChain(chain)] = meanTotalChainCoverage(graph, chain)
	return estimates

def initialChainsEstimates(graph, mapping, estimates):
	return reduce(lambda E,N: addChainEstimate(E, N, graph, mapping), graph.chains, estimates)

def initialEstimates(graph, mapping):
	return initialChainsEstimates(graph, mapping, initialNodesEstimates(graph, mapping, [0 for X in range(mapping.size)]))

##############################################
## Adding in precisions
##############################################

def square(X):
	return X * X

def variance(stddev):
	return square(stddev)

def precision(stddev):
	return 1 / variance(stddev)

def addEdgePrecision(precisions, A, B, graph, mapping):
	if graph[A].edges[B] >= 0:
		precisions[mapping.getEdge(A, B)] = precision(float(max([graph[A].edges[B], EDGE_FACTOR, MIN_NOISE])))
	else:
		precisions[mapping.getEdge(A, B)] = precision(float(MAX_NOISE))
	assert type(precisions[mapping.getEdge(A, B)]) is float
	return precisions

def initialNodePrecisions(precisions, node, graph, mapping):
	return reduce(lambda E,N: addEdgePrecision(E, node, N, graph, mapping), graph[node].edges.keys(), precisions)

def initialNodesPrecisions(graph, mapping, precisions):
	return reduce(lambda E,N: initialNodePrecisions(E, N, graph, mapping), graph.nodes(), precisions)

def addChainPrecision(precisions, chain, graph, mapping):
	precisions[mapping.getChain(chain)] = precision(max(min(meanTotalChainCoverage(graph, chain) / chain.length(), SEG_FACTOR), MIN_NOISE))
	return precisions

def initialChainsPrecisions(graph, mapping, precisions):
	return reduce(lambda P, C: addChainPrecision(P, C, graph, mapping), graph.chains, precisions)

def initialPrecisions(graph, mapping):
	return initialChainsPrecisions(graph, mapping, initialNodesPrecisions(graph, mapping, [0 for X in range(mapping.size)]))

##############################################
## Creating a matrix image of the graph
##############################################

class LPProblem:
	def __init__(self, columns):
		self.columns = columns
		self.matrix = None
		self.constraints = []
		self.constraintPrecisions = []
		self.estimatedValues = None
		self.estimatePrecisions = None

	def addConstraint(self, array, target, stddev):
		if self.matrix is None:
			self.matrix = np.vstack([array])
		else:
			self.matrix = np.append(self.matrix, np.vstack([array]), axis=0)

	 	prior =  target - sum((array * self.estimatedValues))
		self.constraints.append(prior)
		self.constraintPrecisions.append(variance(stddev))

	def __str__(self):
		return "\n".join(map(str, [self.matrix, self.constraints, self.constraintPrecisions, self.estimatedValues, self.estimatePrecisions]))

def edgeConstraint(array, node, tonode, mapping):
	if node != tonode:
		array[mapping.getEdge(node, tonode)] = 1
	else:
		array[mapping.getEdge(node, node)] = 2
	return array

def nodeConstraint(node, graph, mapping):
	array = reduce(lambda A,E: edgeConstraint(A, node, E, mapping), graph[node].edges.keys(),  np.zeros(mapping.size))
	array[mapping.getChain(graph.blockChain[graph.nodeBlock[node]])] = -1
	return array

def prepareNodeProblem(problem, node, graph, mapping):
	problem.addConstraint(nodeConstraint(node, graph, mapping), 0, FUDGE_FACTOR)
	return problem

def prepareGraphProblem(graph, mapping, problem):
	return reduce(lambda P, N: prepareNodeProblem(P, N, graph, mapping), graph.nodes(), problem)

##############################################
## Creating a graph image of the flow solution
##############################################

def updateEdge(graph, tonode, fromnode, values, mapping):
	graph[fromnode].edges[tonode] += values[mapping.getEdge(tonode, fromnode)]
	return graph

def updateNode(graph, node, values, mapping):
	return reduce(lambda A,B: updateEdge(A, B, node, values, mapping), graph[node].edges.keys(), graph)
	
def updateGraphNodes(graph, values, mapping):
	return reduce(lambda A,B: updateNode(A, B, values, mapping), graph.nodes(), graph)

def updateBlock(graph, block, chain, index, newValue):
	A, B = list(block.nodes)
	graph[A].segment[index] = newValue
	graph[B].segment[index] = newValue
	return graph

def updateChainPloidy(graph, chain, values, mapping, index, newValue):
	return reduce(lambda A,B: updateBlock(A, B, chain, index, newValue), chain, graph)

def updateChain(graph, chain, values, mapping):
	oldValues = [meanChainCoverage(graph, chain, index) for index in range(chain[0].ploidy(graph))]
	oldTotal = sum(oldValues)
	if oldTotal > 0:
		delta = values[mapping.getChain(chain)]
		newValues = [ X + delta * X / oldTotal for X in oldValues] 
		return reduce(lambda g,i: updateChainPloidy(g, chain, values, mapping, i, newValues[i]), range(chain[0].ploidy(graph)), graph)
	else:
		return updateChainPloidy(graph, chain, values, mapping, 0, values[mapping.getChain(chain)])

def updateGraphChains(graph, values, mapping):
	return reduce(lambda A,B: updateChain(A, B, values, mapping), graph.chains, graph)

def updateGraph(graph, values, mapping):
	return updateGraphNodes(updateGraphChains(graph, values, mapping), values, mapping)

##############################################
## Correction for overdemanded nodes
##############################################

def correctIncongruity(node, graph):
	seg_val = sum(graph[node].segment)
	edge_sum = sum(graph[node].edges.values())
	
	if edge_sum != 0 and seg_val < edge_sum:
		for X in graph[node].edges:
			graph.changeLiftedEdge(node, X, graph[node].edges[X] * (seg_val - edge_sum) / edge_sum)

def correctInsufficientEdges(node, graph):
	seg_val = sum(graph[node].segment)
	edge_sum = sum(graph[node].edges.values())
	if len(graph[node].edges) == 1 and seg_val <= sum(graph[graph[node].partner].segment) and seg_val > edge_sum:
		graph.changeLiftedEdge(node, graph[node].partner, seg_val - edge_sum)
		
def correctIncongruities(graph):
	map(lambda X: correctIncongruity(X, graph), graph.nodes())
	map(lambda X: correctInsufficientEdges(X, graph), graph.nodes())

##############################################
## Metric
##############################################
def nodeError(node, graph):
	return graph.blockChain[node.block] - sum(graph[node].edges.values())

def nodeErrors(graph):
	return map(lambda X: nodeError(X, graph), graph.nodes())

def graphError(graph):
	return sum(nodeErrors(graph))
	
##############################################
## Offset
##############################################

def weightedChainOffset(cactus, block):
	return sum(block.copynumber(cactus, index) for index in range(block.ploidy(cactus))) * block.length() 

def weightedSegmentOffsets3(cactus, chain):
	return [weightedChainOffset(cactus, block) for block in chain]

def weightedSegmentOffsets2(cactus, chain): 
	return sum(weightedSegmentOffsets3(cactus, chain))

def weightedSegmentOffsets(graph):
	return [weightedSegmentOffsets2(graph, chain) for chain in  graph.chains]

def segmentLengths(graph):
	return map(lambda X: X.length(), graph.chains)

def meanOffset(graph):
	return sum(weightedSegmentOffsets(graph)) / sum(segmentLengths(graph))

def offsetAncestral(val, mean, ploidy):
	return (val / mean *  ploidy - 1) * 2 / ploidy

def offsetNovel(val, mean, ploidy):
	return (val / mean ) * 2

def edgeOffset(nodeFlow, dest, mean):
	if dest == nodeFlow.partner:
		return offsetAncestral(nodeFlow.edges[dest], mean, 1)
	else:
		return offsetNovel(nodeFlow.edges[dest], mean, 1)

def segmentOffset(nodeFlow, index, mean, ploidy):
	if nodeFlow.node.chr is not None:
		return offsetAncestral(nodeFlow.segment[index], mean, ploidy)
	else:
		return offsetNovel(nodeFlow.segment[index], mean, ploidy)
	

def computeNodeOffsets(nodeFlow, mean):
	for index in range(len(nodeFlow.segment)):
		nodeFlow.segment[index] = segmentOffset(nodeFlow, index, mean, len(nodeFlow.segment))
	nodeFlow.edges = dict(map(lambda X: (X, edgeOffset(nodeFlow, X, mean)), nodeFlow.edges))

def computeOffsets(graph):
	mean = meanOffset(graph)
	return map(lambda X: computeNodeOffsets(X, mean), graph.values())

##############################################
## Master function
##############################################

class BalancedCactus(cactus.Cactus):
	"""Cactus graph with balanced copy number flow"""
	def __init__(self, graph):
		print 'Gaussian Belief propagation resolution of flow in Cactus'
		assert type(graph) is normalized.NormalizedCactus
		self.copy(graph)
		mapping = prepareGraphMapping(graph)
		problem = LPProblem(mapping.size)
		problem.estimatedValues = initialEstimates(graph, mapping)
		problem.estimatePrecisions = initialPrecisions(graph, mapping)
		problem = prepareGraphProblem(graph, mapping, problem)
		values, precisions = gabp.runGaBP(TOL, problem.matrix, problem.constraints, problem.constraintPrecisions, problem.estimatedValues, problem.estimatePrecisions)
		updateGraph(self, values, mapping)
		computeOffsets(self)
		correctIncongruities(self)

	##############################################
	## Metric
	##############################################
	def nodeError(self, node):
		return abs(sum(self[node].segment) - sum(self[node].edges.values()))

	def nodeErrors(self):
		return map(self.nodeError, self)

	def graphError(self):
		return sum(self.nodeErrors())

#############################################
## Unit test
#############################################
def main():
	graph = avg.randomNearEulerianGraph(10)
	C = cactus.Cactus(graph)
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print C
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	NC = normalized.NormalizedCactus(C)
	print NC
	BC = BalancedCactus(NC)
	print BC
	print BC.graphError()

if __name__=='__main__':
	main()
