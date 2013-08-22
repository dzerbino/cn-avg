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

"""Definition of sequence graphs"""

import sys
import copy
import random

from cnavg.avg.node import Node
from cnavg.avg.nodeFlow import NodeFlow

try:
    import cPickle 
except ImportError:
    import pickle as cPickle

###############################################
## Sequence graph
###############################################

class Graph(dict):
	"""A sequence graph"""
	def __init__(self):
		super(Graph, self).__init__()
		self.telomeres = set()

	def nodes(self):
		"""Returns list of nodes"""
		return self.keys()

	def createNode(self, chr, pos, orientation, name=None):
		"""Creates a node and adds it to the graph"""
		node = Node(len(self), chr, pos, orientation, name)
		self.addNode(node)
		return node

	def addNode(self, node):
		"""Adds a node to the graph"""
		if node not in self:
			self[node] = NodeFlow(node)

	def addLiftedEdge(self, A, B, value):
		"""Creates (if necessary) and increments the flow of a bond edge"""
		self.addNode(A)
		self.addNode(B)
		self[A].addEdge(B, value)
		self[B].addEdge(A, value)
		return self

	def setLiftedEdge(self, A, B, value):
		"""Creates (if necessary) and updates the flow of a bond edge"""
		self.addNode(A)
		self.addNode(B)
		self[A].setEdge(B, value)
		self[B].setEdge(A, value)
		return self

	def changeLiftedEdge(self, A, B, value):
		"""Increments the flow on a bond edge"""
		assert A in self and B in self
		self[A].addEdge(B, value)
		self[B].addEdge(A, value)
		return self

	def changeSegment(self, node, index, increment):
		"""Increments the flow of the segment edge incident on the given node"""
		assert node in self
		self[node].segment[index] += increment
		self[self[node].twin].segment[index] += increment

	def createSegment(self, A, B, values):
		"""Declares a set of homologous segment edges between two nodes and sets their flows"""
		self.addNode(A)
		self.addNode(B)
		self[A].twin = B
		self[B].twin = A
		if (A is not B):
			self[A].segment = list(values)
			self[B].segment = list(values)
		else:
			self[A].segment = [2 * X for X in values]

	def createAdjacency(self, A, B):
		"""Creates a bond edge and declares it ancestral"""
		assert A.orientation != B.orientation
		self.addNode(A)
		self.addNode(B)
		self[A].partner = B
		self[B].partner = A
		self.addLiftedEdge(A, B, 0)

	def validateNodeAdjacency(self, node):
		"""
		Validation function
		"""
		return self[self[node].partner].partner == node

	def validateNodeSegment(self, node):
		"""
		Validation function
		"""
		if self[node].twin is not None:
			assert self[self[node].twin].twin == node 
			assert all(X[0] == X[1] for X in zip(self[node].segment, self[self[node].twin].segment))
			#if node not in self.telomeres:
			if False:
				if node < self[node].twin:
					if not (node.orientation == True and self[node].twin.orientation == False):
						print node
						print self[node].twin
						assert False
				else:
					if not (node.orientation == False and self[node].twin.orientation == True):
						print node
						print self[node].twin
						assert False
			return True
		else:
			print node
			assert False
			return True

	def validateNodeEdge(self, node, mate):
		"""
		Validation function
		"""
		return node in self[mate].edges and self[node].edges[mate] == self[mate].edges[node]

	def validateNodeEdges(self, node):
		"""
		Validation function
		"""
		return all(map(lambda X: self.validateNodeEdge(node, X), self[node].edges))

	def validateNode(self, node):
		"""
		Validation function
		"""
		self.validateNodeSegment(node) 
		self.validateNodeEdges(node) 
		self.validateNodeAdjacency(node)
		return True

	def validate(self):
		"""
		Validation function
		"""
		return all(map(lambda X: self.validateNode(X), self))

	def coverageStats(self):
		return "\n".join(map(NodeFlow.coverageStats, filter(lambda X: X.node < X.twin, sorted(self.values()))))

	def __str__(self):
		"""
		GraphViz output.
		"""
		return "\n".join(['digraph G {'] 
				 + map(str, self.values())
				 + ['}'])

	def __copy__(self):
		res = Graph()
		res.copy(self)
		return res

	def copy(self, origin):
		for node in origin:
			self[node] = copy.copy(origin[node])
		self.telomeres = copy.copy(origin.telomeres)

	def twin(self, node):
		assert node in self
		return self[node].twin

	def partner(self, node):
		assert node in self
		return self[node].partner

	def segment(self, node, index):
		assert node in self
		return self[node].segment[index]

	def edge(self, node, node2):
		assert node in self
		return self[node].edges[node2]

	def edges(self, node):
		assert node in self
		return self[node].edges.keys()

	def balanced(self):
		return all(map(lambda X: X.balanced(), self.values()))

	def segmentLength(self, node):
		return self[node].segmentLength()

	def printCNVs(self):
		return "\n".join(filter(lambda X: X is not None, (self[X].printCNVs(self) for X in sorted(self))))

	def printMajority(self):
		return "\n".join(filter(lambda X: X is not None, (self[X].printMajority(self) for X in sorted(self))))

	def printMinority(self):
		return "\n".join(filter(lambda X: X is not None, (self[X].printMinority(self) for X in sorted(self))))

def testGraph():
	graph = Graph()
	node = map(lambda X: graph.createNode("1", X * 50, X % 2 == 0), range(0, 10, 2))
	node += map(lambda X: graph.createNode("1", (X+1) * 50 - 1, X % 2 == 0), range(1, 10, 2))
	node.sort()
	map(lambda X: graph.createSegment(node[X], node[X+1], [1]), range(0,10,2))
	map(lambda X: graph.createAdjacency(node[X], node[X+1]), range(1,9,2))
	graph.createAdjacency(node[0], node[-1])

	graph.addLiftedEdge(node[0], node[9], 1)
	graph.addLiftedEdge(node[1], node[3], 1)
	graph.addLiftedEdge(node[2], node[4], 1)
	graph.addLiftedEdge(node[5], node[8], 1)
	graph.addLiftedEdge(node[6], node[7], 1)
	graph.changeSegment(node[6], 1, -1)

	return graph

def othernodes(node, graph):
	return filter(lambda X: X != node and X != graph[node].twin, graph.nodes())

def otherchoice(node, graph):
	return random.choice(othernodes(node, graph))

def linearGraph(cardinal):
	"""Creates a circular sequence graph with a set number of segments"""
	graph = Graph()
	node = map(lambda X: graph.createNode("1", X * 50000, X % 2 == 0), range(0, cardinal, 2))
	node += map(lambda X: graph.createNode("1", (X+1) * 50000 - 1, X % 2 == 0), range(1, cardinal, 2))
	node.sort()
	map(lambda X: graph.createSegment(node[X], node[X+1], [0]), range(0,cardinal,2))
	map(lambda X: graph.createAdjacency(node[X], node[X+1]), range(1,cardinal - 1,2))
	graph.createAdjacency(node[0], node[-1])
	graph.telomeres.add(node[0])
	graph.telomeres.add(node[-1])
	return graph

def randomGraph(cardinal):
	"""Creates a random sequence graph"""
	graph = linearGraph(cardinal)

	for node in graph.nodes():
		edges = random.randrange(1,3)
		for index in range(edges):
			graph.addLiftedEdge(node, otherchoice(node, graph), 1)

	graph.validate()
	return graph

def randomEulerianGraph(cardinal):
	"""Creates a random sequence graph with balanced flow"""
	graph = linearGraph(cardinal)

	for cycle in range(random.randrange(1,cardinal/3)):
		start = random.choice(graph.nodes())
		graph.changeSegment(start, 0, 1)
		node = graph[start].twin
		constant = random.random()

		next = otherchoice(node, graph)
		graph.addLiftedEdge(node, next, 1)
		graph.changeSegment(next, 0, 1)
		node = graph[next].twin

		while random.random() < constant:
			next = otherchoice(node, graph)
			graph.addLiftedEdge(node, next, 1)
			graph.changeSegment(next, 0, 1)
			node = graph[next].twin

		if node == start:
			next = otherchoice(node, graph)
			graph.addLiftedEdge(node, next, 1)
			graph.changeSegment(next, 0, 1)
			node = graph[next].twin

		graph.addLiftedEdge(node, start, 1)

	return graph

def addNodeNoise(graph, node):
	""" Internal """
	for nodeB in graph[node].edges:
		graph.changeLiftedEdge(node, nodeB, max(-graph[node].edges[nodeB], random.gauss(0, graph[node].edges[nodeB])/5))
	index = 0
	graph.changeSegment(node, index, random.gauss(0, graph[node].segment[index] / graph[node].segmentLength()))
	return graph	

def addNoise(graph):
	""" Internal """
	return reduce(addNodeNoise, graph.nodes(), graph)

def randomNearEulerianGraph(cardinal):
	"""Creates a random sequence graph with near-balanced flow (gaussian noise function)"""
	return addNoise(randomEulerianGraph(cardinal))

#######################################################
## Test
#######################################################
def main():
	graph = testGraph()
	print graph
	print graph.coverageStats()
	print linearGraph(100)
	print randomGraph(100)
	print randomEulerianGraph(100)
	print randomNearEulerianGraph(10)

if __name__ == '__main__':
	main()
