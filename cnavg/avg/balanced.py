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

"""Definition of copy-number balancing across a sequence graph"""
 
import sys
import graph as avg
#import gabp.gabp as gabp
import numpy as np
import cnavg.basics.leastSquares as leastSquares

TOL = 1e2
FUDGE_FACTOR = 1
SEG_FACTOR = FUDGE_FACTOR / 10
EDGE_FACTOR = FUDGE_FACTOR * 10
MAX_NOISE = 1e3
MIN_NOISE = 1e-3

##############################################
## Mapping of flow elements onto integer IDs 
##############################################

class Mapping(object):
	"""
	An index of all the adjacencies and segments in a graph
	"""
	def __init__(self):
		self.segments = dict()
		self.edges = dict()
		self.size = 0

	def addSegment(self, A):
		self.segments[A] = self.size
		self.size += 1

	def getSegment(self, A, B):
		return self.segments[min([A,B])]

	def addEdge(self, A, B):
		m = min([A, B])
		M = max([A, B])
		self.edges[(m,M)] = self.size
		self.size += 1

	def getEdge(self, A, B):
		m = min([A, B])
		M = max([A, B])
		return self.edges[(m, M)]

	def segmentStr(self, A):
		return "\t".join([str(A), str(self.segments[A])])

	def segmentsStr(self):
		return "\n".join(map(lambda X: self.segmentStr(X), self.segments))

	def edgeStr(self, pair):
		return "\t".join([str(pair[0]), str(pair[1]), str(self.edges[(pair[0],pair[1])])])

	def edgesStr(self):
		return "\n".join(map(lambda X: self.edgeStr(X), self.edges))

	def __str__(self):
		return "\n".join(sorted([self.edgesStr(), self.segmentsStr()]))

def prepareEdgeMapping (mapping, tonode, fromnode):
	if fromnode <= tonode:
		mapping.addEdge(tonode, fromnode)
	return mapping

def prepareNodeMapping(mapping, node, graph):
	if node <= graph[node].twin:
		mapping.addSegment(node)
	return reduce(lambda M, N: prepareEdgeMapping(M, N, node), graph[node].edges.keys(), mapping)

def prepareGraphMapping(graph):
	"""
	Builds a mapping for the given sequence graph
	"""
	return reduce(lambda M,X: prepareNodeMapping(M,X, graph), graph.nodes(), Mapping())

##############################################
## Adding in estimates
##############################################

def addSegmentEstimate(estimates, node, graph, mapping):
	estimates[mapping.getSegment(node, graph[node].twin)] = float(sum(graph[node].segment))
	return estimates

def addEdgeEstimate(estimates, A, B, graph, mapping):
	if graph[A].edges[B] >= 0:
		estimates[mapping.getEdge(A, B)] = float(graph[A].edges[B])
	else: 
		estimates[mapping.getEdge(A, B)] = 0
	return estimates

def addNodeEdgesEstimates(estimates, node, graph, mapping):
	return reduce(lambda E,N: addEdgeEstimate(E, node, N, graph, mapping), graph[node].edges.keys(), estimates)

def initialNodeEstimates(estimates, node, graph, mapping):
	estimates = addSegmentEstimate(estimates, node, graph, mapping)
	return addNodeEdgesEstimates(estimates, node, graph, mapping)

def initialEstimates(graph, mapping):
	"""
	Prepares vector of observed values for the Gaussian propagation routine.
	"""
	return reduce(lambda E,N: initialNodeEstimates(E, N, graph, mapping), graph.nodes(), [0 for x in range(mapping.size)])

##############################################
## Adding in precisions
##############################################

def square(X):
	return X * X

def variance(stddev):
	return square(stddev)

def precision(stddev):
	return 1 / variance(stddev)

def addSegmentPrecision(precisions, node, graph, mapping):
	precisions[mapping.getSegment(node, graph[node].twin)] = precision(max(min(sum(graph[node].segment) / graph.segmentLength(node), SEG_FACTOR), MIN_NOISE))
	return precisions

def addEdgePrecision(precisions, A, B, graph, mapping):
	if graph[A].edges[B] >= 0:
		precisions[mapping.getEdge(A, B)] = precision(float(max([graph[A].edges[B], EDGE_FACTOR, MIN_NOISE])))
	else:
		precisions[mapping.getEdge(A, B)] = precision(float(MAX_NOISE))
	assert type(precisions[mapping.getEdge(A, B)]) is float
	return precisions

def addNodeEdgesPrecisions(precisions, node, graph, mapping):
	return reduce(lambda E,N: addEdgePrecision(E, node, N, graph, mapping), graph[node].edges.keys(), precisions)

def initialNodePrecisions(precisions, node, graph, mapping):
	precisions = addSegmentPrecision(precisions, node, graph, mapping)
	return addNodeEdgesPrecisions(precisions, node, graph, mapping)

def initialPrecisions(graph, mapping):
	"""
	Prepares vector of precisions (=1/variance) of the different observed values and constraints in the Gaussian factor graph.
	"""
	return reduce(lambda E,N: initialNodePrecisions(E, N, graph, mapping), graph.nodes(), [0 for x in range(mapping.size)])

##############################################
## Creating a matrix image of the graph
##############################################

class LPProblem:
	"""A formalization of a Gaussian propagation factor graph"""
	def __init__(self, columns):
		self.columns = columns
		self.matrix = None
		self.constraints = []
		self.constraintPrecisions = []
		self.estimatedValues = None
		self.estimatePrecisions = None

	def addConstraint(self, array, stddev):
		if self.matrix is None:
			self.matrix = np.vstack([array])
		else:
			self.matrix = np.append(self.matrix, np.vstack([array]), axis=0)

	 	prior = - sum((array * self.estimatedValues))
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
	array = reduce(lambda A,E: edgeConstraint(A, node, E, mapping), graph[node].edges.keys(), np.zeros(mapping.size))
	array[mapping.getSegment(node, graph[node].twin)] = -1
	return array

def prepareNodeProblem(problem, node, graph, mapping):
	problem.addConstraint(nodeConstraint(node, graph, mapping), FUDGE_FACTOR)
	return problem

def prepareGraphProblem(graph, mapping, problem):
	"""Formulates a graph balancing problem into a Linear Programming problem"""
	return reduce(lambda P, N: prepareNodeProblem(P, N, graph, mapping), graph.nodes(), problem)

##############################################
## Debugging option
##############################################
def copyCatResolution(graph):
	"""Direct resolution of the problem using linear algebra"""
	mapping = prepareGraphMapping(graph)
	problem = LPProblem(mapping.size)
	problem.estimatedValues = initialEstimates(graph, mapping)
	problem.estimatePrecisions = initialPrecisions(graph, mapping)
	problem = prepareGraphProblem(graph, mapping, problem)

	A = problem.matrix
	y = problem.constraints
	sigma_y = problem.constraintPrecisions
	x = [0 for X in problem.estimatedValues]
	sigma_x = problem.estimatePrecisions
	MU = np.hstack([np.diag(sigma_y), A])
	ML = np.hstack([A.T, np.diag(sigma_x)])
	M = np.vstack([MU, ML])
	res = np.linalg.solve(M, np.concatenate((y,x)))
	return res[len(y):], mapping

##############################################
## Debugging option number 2
##############################################
def copyCatResolution2(graph):
	"""Direct resolution of the problem using weighted least squares algorithm"""
	"""See http://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares"""
	mapping = prepareGraphMapping(graph)
	problem = LPProblem(mapping.size)
	problem.estimatedValues = initialEstimates(graph, mapping)
	problem.estimatePrecisions = initialPrecisions(graph, mapping)
	problem = prepareGraphProblem(graph, mapping, problem)

	x = [0 for X in problem.estimatedValues]
	res = leastSquares.solve(x, problem.estimatePrecisions, problem.matrix, problem.constraints, problem.constraintPrecisions)
	return res, mapping

##############################################
## Master function
##############################################

class BalancedAVG(avg.Graph):
	"""A sequence graph characterised by balanced flow (i.e. the Laplacian of the conjugate flow is null)"""
	def __init__(self, graph):
		self.copy(graph)
		values, mappings = copyCatResolution2(graph)
		self.updateGraph(values, mappings)
		self.correctIncongruities()
		return

		#print 'Gaussian Belief propagation resolution of flow'
		#mapping = prepareGraphMapping(self)
		#problem = LPProblem(mapping.size)
		#problem.estimatedValues = initialEstimates(self, mapping)
		#problem.estimatePrecisions = initialPrecisions(self, mapping)
		#problem = prepareGraphProblem(self, mapping, problem)

		#values, precisions = gabp.run(TOL, problem.matrix, problem.constraints, problem.constraintPrecisions, [0 for X in problem.estimatedValues], problem.estimatePrecisions)
		#self.updateGraph(values, mapping)

		#self.correctIncongruities()

	##############################################
	## Creating a graph image of the flow solution
	##############################################

	def updateEdge(self, tonode, fromnode, values, mapping):
		if self[fromnode].edges[tonode] >= 0:
			self[fromnode].edges[tonode] += values[mapping.getEdge(tonode, fromnode)]
		else:
			self[fromnode].edges[tonode] = values[mapping.getEdge(tonode, fromnode)]

	def updateNode(self, node, values, mapping):
		oldSum = sum(self[node].segment)
		if oldSum > 0:
			for index in range(len(self[node].segment)):
				self[node].segment[index] += values[mapping.getSegment(node, self[node].twin)] * self[node].segment[index] / oldSum
		else:
			self[node].segment[0] = values[mapping.getSegment(node, self[node].twin)]
		map(lambda X: self.updateEdge(X, node, values, mapping), self[node].edges)
		
	def updateGraph(self, values, mapping):
		"""Change the flow values to those produced by the LP method"""
		map(lambda X: self.updateNode(X, values, mapping), self.nodes())

	##############################################
	## Correction for overdemanded nodes
	##############################################

	def correctIncongruity(self, node):
		seg_val = sum(self[node].segment)
		edge_sum = sum(self[node].edges.values())
		
		if edge_sum != 0 and seg_val < edge_sum:
			for X in self[node].edges:
				self.changeLiftedEdge(node, X,  - self[node].edges[X] * (edge_sum - seg_val) / edge_sum)
			
	def correctIncongruities(self):
		"""Ensure that segment flow >= bond flow, under the assumption that missing breakpoints are much more common than overestimated copy number"""
		map(self.correctIncongruity, self.nodes())

	##############################################
	## Metric
	##############################################

	def nodeError(self, node):
		return abs(sum(self[node].segment) - sum(self[node].edges.values()))

	def nodeErrors(self):
		return map(self.nodeError, self)

	def graphError(self):
		"""Measure the amount the imbalance on the graph flow"""
		return sum(self.nodeErrors())

	##############################################
	## Validation
	##############################################

	def validateNode(self, node):
		if self[node].twin is not None:
			assert len(self[node].edges) > 1 or sum(self[node].segment) > sum(self[self[node].partner].segment) or sum(self[node].edges.values()) >= sum(self[node].segment) - 1e-3, "\n".join(map(str, [node, sum(self[node].segment), sum(self[self[node].partner].segment), sum(self[node].edges.values())]))
		return True

	def validate(self):
		super(BalancedGraph, self).validate()
		assert all(map(self.validateNode, self.nodes()))
	
#############################################
## Unit test
#############################################
def main():
	graph = avg.randomNearEulerianGraph(10)
	print graph
	balancedB = BalancedAVG(graph)
	print balancedB
	print balancedB.graphError()

if __name__=='__main__':
	main()
