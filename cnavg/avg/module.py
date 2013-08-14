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

"""Definition of module around nets"""

import sys
import cnavg.avg.graph as avg
from cnavg.avg.node import StubNode
import cnavg.cactus.graph as cactus
import cnavg.cactus.oriented as oriented
import cnavg.cactusSampling.sampling as normalized
import cnavg.cactus.balanced as balanced
import copy


class Module(avg.Graph):
	"""Sequence graph associated to a Cactus net"""
	def __init__(self, net=None, graph=None, cnvs=None):
		super(Module, self).__init__()

		if net is None:
			return
		self.net = net

		for group in net.groups:
			for node in group:
				self.addNetEnd(node, graph)

		self.telomeres = set(filter(lambda X: X in graph.telomeres, self.nodes()))

		if net in graph.headChain:
			self.pseudotelomeres = set(filter(lambda X: X in self.nodes(), [N for B in graph.headChain[net] for N in B.nodes]))
		else:
			self.pseudotelomeres = set()

		self.close(graph, cnvs)

		if not all(self[X].twin is not None for X in self.nodes()):
			print self
			print "\n".join(map(str,filter(lambda X: self[X].twin is None, sorted(self.nodes()))))
			assert False

		self.stub = StubNode(len(graph))
		self.addNode(self.stub)
		self[self.stub].partner = self.stub
		for node in self.nodes(): 
			if node != self.stub:
				excess = sum(self.segments[node]) - sum(self[node].edges.values())
				if excess != 0:
					self.addLiftedEdge(node, self.stub, excess)

		stubFlow = sum(self[self.stub].edges.values())
		self.createModuleSegment(self.stub, self.stub,[stubFlow/2])
		assert self.balanced()

	def copy(self, other):
		""" Copy data onto other instance """
		super(Module, self).copy(other)
		self.pseudotelomeres = copy.copy(other.pseudotelomeres)
		self.stub = other.stub
		self.segments = dict()
		for X in other.segments:
			self.segments[X] = copy.copy(other.segments[X])

	def __copy__(self):
		res = Module()
		res.copy(self)
		return res

	def segmentNodeString(self, node, index):
		""" GraphViz output """ 
		return "\t%s -> %s [style=dashed, color=green, label=%f]" % (node.ID, self[node].twin.ID, -self.segments[node][index])

	def segmentNodeStrings(self, node):
		""" GraphViz output """ 
		if self[node].twin is not None and node <= self[node].twin:
			return map(lambda X: self.segmentNodeString(node, X), range(len(self.segments[node]))) 
		else:
			return []

	def segmentStrings(self):
		""" GraphViz output """ 
		return sum(map(lambda X: self.segmentNodeStrings(X), self), [])

	def __str__(self):
		""" GraphViz output """ 
		return "\n".join(['digraph G {'] 
				 + ['node [style=filled,color=purple]']
				 + ["\t" + str(X.ID) for X in self.pseudotelomeres]
				 + ['node [style=filled,color=red]']
				 + ["\t" + str(X.ID) for X in self.telomeres]
				 + ['node [style=line,color=black]']
				 + map(str, self.values())
				 + self.segmentStrings()
				 + ['}'])

	def __hash__(self):
		return id(self)

	def addModuleSegment(self, nodeA, values):
		""" Define segment flows incident onto nodeA """
		if nodeA not in self.segments:
			self.segments[nodeA] = values
		else:
			self.segments[nodeA].extend(values)

	def createModuleSegment(self, nodeA, nodeB, values):
		self.addModuleSegment(nodeA, values)
		if nodeB != nodeA:
			self.addModuleSegment(nodeB, values)
		if self[nodeA].twin is None:
			self.createSegment(nodeA, nodeB, values)
		else:
			self.changeSegment(nodeA, values)

	def changeModuleSegment(self, nodeA, nodeB, value, index):
		self.segments[nodeA][index] += value
		self.segments[nodeB][index] += value

	def removeEdgeFlow(self, edge):
		if edge.index == -1:
			self.changeLiftedEdge(edge.start, edge.finish, -edge.value)
		else:
			self.changeModuleSegment(edge.start, edge.finish, edge.value, edge.index)

	def addNetEnd(self, node, graph):
		self.addNode(node)

		if graph[node].partner in self.nodes():
			self.createAdjacency(node, graph[node].partner)
		for dest in graph[node].edges:
			if dest in self:
				self.addLiftedEdge(node, dest, graph[node].edges[dest])

	def close(self, graph, cnvs):
		""" Connect pseudotelomeres as virtual twins """
		self.segments = dict()
		nodes = filter(lambda X: X.orientation == True, sorted(self.nodes()))
		if len(nodes) < 1:
			return
		else:
			PT = list(self.pseudotelomeres)
			assert len(PT) == 0 or len(PT) == 2
			if len(PT) == 2:
				head = PT[0]
				next = PT[1]
				self.createSegment(head, next, [0])
				self.segments[head] = [-X[0].value for X in cnvs]
				self.segments[next] = [-X[0].value for X in cnvs]

			next = nodes.pop(0) 
			if next in self.telomeres:
				chromosomeEnd = self[next].partner
			else:
				chromosomeEnd = None
			for node in nodes:
				partner = self[node].partner
				if next.chr == node.chr:
					self.createModuleSegment(next, partner, graph[next].segment)
				else:
					self.createModuleSegment(next, chromosomeEnd, graph[next].segment)
					chromosomeEnd = partner
				next = node
			if chromosomeEnd is not None:
				self.createModuleSegment(next, chromosomeEnd, graph[next].segment)

	def nodeFlow(self, node):
		""" Flow imbalance around a given node """
		if self[node].twin != node:
			return sum(self[node].edges.values()) - sum(self.segments[node])
		else:
			return sum(self[node].edges.values()) - 2 * sum(self.segments[node])

	def nodeBalanced(self, node):
		""" Check whether flow is balanced (approximately) around a given node """
		if self.nodeFlow(node) >= 1e-3:
			print node
			print sum(self[node].edges.values())
			print sum(self.segments[node])
		assert abs(self.nodeFlow(node)) < 1e-3
		return True

	def balanced(self):
		""" Check whether flow is balanced """
		assert all(self.nodeBalanced(X) for X in self.nodes())
		return True

def netModulePairs(graph):
	""" Enumerate Nets paired with their module in a Cactus graph """
	return map(lambda X: (X, Module(X, graph, [])), graph.nets)

###############################################
## Master function
###############################################

def extractGraphModules(graph):
	"""Returns a dictionary which maps each Net of a Cactus graph to its own module"""
	return dict(netModulePairs(graph))

###############################################
## Test 
###############################################
def main():
	G = avg.randomNearEulerianGraph(10)
	C = cactus.Cactus(G)
	N = normalized.NormalizedCactus(C)
	B = balanced.BalancedCactus(N)
	O = oriented.OrientedCactus(B)
	print C
	print '>>>>>>>>>>>>>>>>>>>>>>'
	print "\n".join(map(str, extractGraphModules(O).values()))

if __name__ == "__main__":
	main()
