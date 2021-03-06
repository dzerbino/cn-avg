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

"""
Definition of NodeFlow, a static node's dynamic context within a graph
"""

import sys
import copy
import random

from cnavg.avg.node import Node

class NodeFlow(object):
	"""
	Flow data of the edges incident on a given node
	"""
	def __init__(self, node):
		self.node = node
		self.edges = dict()
		self.twin = None
		self.segment = []
		self.partner = None
		self.selfLoops = False

	def __copy__(self):
		new = NodeFlow(self.node) 
		new.twin = self.twin
		new.segment = copy.copy(self.segment)
		new.partner = self.partner
		for dest in self.edges:
			new.edges[dest] = self.edges[dest]
		if hasattr(self, 'selfLoops'):
			new.selfLoops = self.selfLoops
		else:
			new.selfLoops = (self.node in self.edges)
		return new

	def __cmp__(self, other):
		"""
		Ordering by node position.
		"""
		return cmp(self.node, other.node)

	def adjacencyString(self):
		"""
		GraphViz output.
		"""
		return "\t%s -> %s [color=black]" % (self.node.ID, self.partner.ID)

	def segmentString(self):
		"""
		GraphViz output.
		"""
		return "\n".join("\t%s -> %s [color=red, label=%f]" % (self.node.ID, self.twin.ID, val) for val in self.segment)

	def liftedEdgeString(self, dest):
		"""
		GraphViz output.
		"""
		return "\t%s -> %s [color=green, label=%f]" % (self.node.ID, dest.ID, self.edges[dest])

	def __str__(self):
		"""
		String representation (debugging purposes).
		"""
		if self.twin is not None and self.node <= self.twin:
			segmentStr = [self.segmentString()]
		else:
			segmentStr = []

		if self.node <= self.partner:
			adjacencyStr = [self.adjacencyString()]
		else:
			adjacencyStr = []

		liftedEdgeStr = [self.liftedEdgeString(dest) for dest in self.edges if self.node <= dest] 
		return "\n".join([str(self.node)] + segmentStr + adjacencyStr + liftedEdgeStr)

	def segmentValue(self, index):
		"""
		Returns the flow of an incident segment edge.
		"""
		if index >= len(self.segment):
			return "NA"
		else:
			return str(self.segment[index])

	def coverageStats(self):
		""" 
		Returns string with coverage stats (debugging purposes)
		"""
		return "%s\t%i\t%i\t%f\t%f\t" % (self.node.chr, self.node.pos, self.twin.pos) + "\t".join(map(self.segmentValue, range(2)))

	def addEdge(self, node, flow):
		"""
		Adds a bond between two nodes (one way only).
		"""
		if node in self.edges:
			self.edges[node] += flow
		else:
			self.edges[node] = flow
		if node is self.node:
			self.selfLoops = True

	def setEdge(self, node, flow):
		"""
		Sets the bond flow between two nodes (one-way only).
		"""
		self.edges[node] = flow

	def flow(self):
		"""
		Returns the sum of incident bonds flows (i.e. the sum of incident segment flows).
		"""
		return sum(self.edges.values())

	def balanced(self):
		"""
		Tests whether the incident flows are (approximately) balanced.
		"""
		if abs(self.flow() - sum(self.segment)) >= 1e-5:
			print self
			print self.flow()
			print self.segment
		assert abs(self.flow() - sum(self.segment)) < 1e-5
		return True

	def segmentLength(self):
		"""
		Length of the incident segment edges
		"""
		if self.node > self.twin:
			return self.node - self.twin
		else:
			return self.twin - self.node

	def inRegion(self, chr, start, finish):
		"""
		Tests whether the incident segment edges overlap with a given chromosomal region
		"""
		segChr = self.node.chr
		segStart = min(self.node, self.twin)
		segFinish = max(self.node, self.twin)
		return segChr == chr and segStart.pos <= finish and segFinish.pos >= start

	def ploidy(self):
		""" 
		Returns the number of incident segment edges.
		"""
		return len(self.segment)

	def printCNVs(self, graph):
		""" 
		Returns the segment segment in BED format.
		"""
		if self.node.chr == "None" or self.node in graph.telomeres or self.twin in graph.telomeres:
			return None
		if self.node.pos < self.twin.pos and len(self.segment) > 0:
			return "%s\t%i\t%i\t%f" % (self.node.chr, self.node.pos, self.twin.pos, sum(self.segment))
		else:
			return None

	def printMajority(self, graph):
		""" 
		Returns the flow of the majority segment in BED format.
		"""
		if self.node.chr == "None" or self.node in graph.telomeres or self.twin in graph.telomeres:
			return None
		if self.node.pos < self.twin.pos and len(self.segment) > 1:
			return "%s\t%i\t%i\t%f" % (self.node.chr, self.node.pos, self.twin.pos, max(self.segment))
		else:
			return None

	def printMinority(self, graph):
		""" 
		Returns the flow of the minority segment in BED format.
		"""
		if self.node.chr == "None" or self.node in graph.telomeres or self.twin in graph.telomeres:
			return None
		if self.node.pos < self.twin.pos and len(self.segment) > 1:
			return "%s\t%i\t%i\t%f" % (self.node.chr, self.node.pos, self.twin.pos, min(self.segment))
		else:
			return None
