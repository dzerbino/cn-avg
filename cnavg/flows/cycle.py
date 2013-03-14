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

"""Cycle traversal underlying an alternating flow"""

import copy

from cnavg.flows.edge import Edge

###############################################
## Cycle
###############################################
class Cycle(list):
	###########################################
	## Basic
	###########################################
	def __init__(self, iter, value=None, conserve=False):
		if not conserve:
			for edge in iter:
				self.append(copy.copy(edge))
		else:
			for edge in iter:
				self.append(edge)

		if value is not None:
			self.setRatio(value)
		elif len(self) > 0:
			self.value = self[0].value
		else:
			self.value = None
		if self.value == None:
			print self
			assert False

	def __copy__(self):
		return Cycle(copy.copy(X) for X in self)

	def __hash__(self):
		return id(self)

	def nodes(self):
		return [X.start for X in self]

	def reverse(self):
		return Cycle([X.reverse() for X in reversed(self)])

	def startAt(self, index, conserve=False):
		return Cycle(self[index:] + self[:index], conserve=conserve)

	def startAtEdge(self, edge):
		nodes = set(edge.nodes())
		index = edge.index
		hits = filter(lambda X: set(X[1].nodes()) == nodes and X[1].index == index, enumerate(self))
		hit = hits[0][0]
		return self.startAt(hit, conserve=True)

	def __add__(self, other):
		return Cycle(list(self) + list(other), self.value)

	def setRatio(self, ratio):
		for index in range(len(self)):
			if index % 2 == 0:
				self[index].value = ratio
			else:
				self[index].value = -ratio
		self.value = ratio

	###########################################
	## Stats
	###########################################

	def getCNVs(self, hash, event, graph):
		return reduce(lambda X, Y: Y.getCNVs(X, event, graph), self, hash)

	###########################################
	## Validation
	###########################################
	def isClosed(self):
		if len(self) < 1:
			print len(self)
			print self
		assert self[-1].finish == self[0].start
		return True

	def alternates(self):
		assert all(self[X].value * self[X + 1].value < 0 for X in range(len(self) - 1))
		return True

	def isConnected(self):
		assert all(self[X].finish == self[X + 1].start for X in range(len(self) - 1))
		return True

	def validate(self):
		assert self.isClosed() 
		assert len(self) % 2 == 0 
		assert self.isConnected() 
		assert self.alternates() 
		return True

	###########################################
	## Display
	###########################################

	def __str__(self):
		return "\n".join(["CYCLE " + str(self.value) + " " + str(len(self))] + [X.shortString() for X in self])

	def dot(self):
		return "\n".join(X.dot() for X in self)

	def braneyText(self, historyID, netID, cycleID, order, complexity, ptr):
		return "\n".join(self[X].braneyText(historyID, netID, cycleID, X, order, complexity, ptr) for X in range(len(self)))

	def simplifyStubsAndTrivials(self, cactus):
		nodes = [X.start for X in self]
		indices = [X.index for X in self]
		values = [X.value for X in self]

		while len(nodes) > 2 and ((nodes[0].chr == "None" and nodes[-1].chr == "None")
		       or (nodes[0] in cactus and nodes[0] not in cactus.telomeres and cactus[nodes[0]].partner == nodes[-1]) and indices[0] > -1 and indices[-2] > -1):
			nodes = nodes[1:-1]
			indices = indices[1:-1]
			values = values[1:-1]

		silent = False
		nodesOut = []
		valuesOut = []
		indicesOut = []
		for i in range(len(nodes) -1):
			if silent:
				silent = False
			elif (nodes[i].chr == "None" and nodes[i+1].chr == "None") or (nodes[i] in cactus and nodes[i] not in cactus.telomeres and cactus[nodes[i]].partner == nodes[i+1] and indices[i+1] > -1 and indices[i-1] > -1):
				silent = True
			else:
				nodesOut.append(nodes[i])
				valuesOut.append(values[i])
				indicesOut.append(indices[i])
				silent = False

		if len(nodes) == 0:
			return None
		if not silent:
			nodesOut.append(nodes[-1])
			valuesOut.append(values[-1])
			indicesOut.append(indices[-1])

		if len(nodesOut) == 0:
			return None

		return Cycle([Edge(nodesOut[i-1], nodesOut[i], valuesOut[i-1], index=indicesOut[i-1]) for i in range(len(nodesOut))])

	def doesNotContainNodes(self, nodes):
		return all(X.doesNotContainNodes(nodes) for X in self)
