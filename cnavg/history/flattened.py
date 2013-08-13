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
Flattening history so that they are no longer split by net
"""

import copy

from constrained import ConstrainedHistory

from cnavg.flows.edge import Edge
from cnavg.flows.cycle import Cycle
from cnavg.flows.flows import Event
from cnavg.history.history import History

import cnavg.avg.graph as avg
import cnavg.cactus.graph as cactus
import cnavg.cactus.balanced as balanced
import cnavg.cactusSampling.sampling as normalized
import cnavg.avg.balanced as balancedAVG
import cnavg.cactus.oriented as oriented
import cnavg.historySampling.cycleCover as cycleCover

def flattenGraph(cactusHistory):
	""" Returns a flattened copy of the CactusHistory """
	new = FlattenedHistory(cactusHistory.cactus)
	new.copy(cactusHistory)
	oldHistories = new.netHistories.values()
	new.netHistories = new._unifyCactusHistory()
	for history in oldHistories:
		for event in history.events:
			if _pseudotelomeric(event, history.module) or _telomeric(event, history.module):
				new.popEvent(event)
	return new

def _telomeric(event, module):
	""" Checks whether any of the nodes traversed by an event flow are telomeres """
	return any(X.start in module.telomeres and X.finish in module.telomeres and X.index > -1 for X in event.cycle)

def _pseudotelomeric(event, module):
	""" Checks whether any of the nodes traversed by an event flow are pseudotelomeres """
	return any(X.start in module.pseudotelomeres and X.finish in module.pseudotelomeres and X.index > -1 for X in event.cycle)


class FlattenedHistory(ConstrainedHistory):
	""" A history where events are not partitioned within nets, instead unified so as to span the whole sequence graph """

	##########################################
	## Basic
	##########################################
	def __init__(self, cactus):
		super(FlattenedHistory, self).__init__(cactus)

	def __copy__(self):
		new = FlattenedHistory(self.cactus)
		new.copy(self)
		return new

	##########################################
	## Unifying cycles
	##########################################
	def _otherTelomere(self, node):
		net = self.cactus.nodeNet(node)		
		module = self.netHistories[net].module
		return module[node].twin

	def _enumerateChain_Block(self, list, finish):
		last = list[-1]
		intermediary = self._otherTelomere(last)
		block = self.cactus.nodeBlock[intermediary]
		newNode = block.otherNode(intermediary)
		if newNode == finish:
			return list
		else:
			return self._enumerateChain_Block(list + [newNode], finish)

	def _enumerateChain(self, edge):
		""" Enumerate blocks in chain finishing at edge """
		block = self.cactus.nodeBlock[edge.start]
		start = block.otherNode(edge.start)
		finish = edge.finish
		return self._enumerateChain_Block([start], finish)

	def _blockEdge(self, node, value):
		""" Replacement edge to represent block """
		return Edge(node, self.cactus[node].twin, value, index=0)

	def _extendEdgeIntoNet(self, node, edge):
		""" Expand edge as path through net """
		net = self.cactus.nodeNet(node)		
		if net == self.cactus.rootNet:
			assert node in self.netHistories[net].module
			assert False
		chain = self.cactus.headChain[net]
		chainCNVs = self.chainCNVs[chain]
		module = self.netHistories[net].module
		twin = module[node].twin
		value = edge.value

		cnvs = filter(lambda X: X[1][0] == edge, enumerate(chainCNVs))
		if len(cnvs) != 1:
			print "NET",id(net)
			print "CHAIN",id(chain)
			print "CHAINCNVS",id(chainCNVs)
			print len(cnvs)
			print edge, id(edge)
			print len(chainCNVs)
			print edge.start
			headNet = self.cactus.headNet[chain]
			history = self.netHistories[headNet]
			module2 = history.module
			print module2[edge.start].twin
			print module2[edge.finish].twin
			for cnv in chainCNVs:
				print '>>>>>>>>>>>>>>>>>>>>>'
				print cnv
				print cnv[0], id(cnv[0])
		assert len(cnvs) == 1
		cnv = cnvs[0]
		index = cnv[0]

		aIndex = (min(node, twin), max(node, twin), index)
		table = self.netHistories[net].overlapTable

		if aIndex not in table or len(table[aIndex]) == 0:
			return [Edge(node, twin, value, index)] + [self._blockEdge(module[node].twin, edge.value)]

		chosenEvents = table[aIndex].keys()
		if len(chosenEvents) != 1:
			print chainCNVs
			print edge.value
			print list(enumerate(chainCNVs))
			print cnvs
			print cnv
			print index
			print module
			print len(chosenEvents)
			for event in chosenEvents:
				print event.dot()
		assert len(chosenEvents) == 1
		chosenEvent = chosenEvents[0]
		chosenCycle = chosenEvent.cycle
		childCycle = chosenCycle.startAtEdge(Edge(node, twin, value, index=index))
		recursedCycle = self._unifyCycle(childCycle[1:])
		if childCycle[0].start == twin:
			return recursedCycle + [self._blockEdge(twin, value)]
		else:
			return recursedCycle.reverse() + [self._blockEdge(twin, value)]

	def _extendEdgeIntoChain(self, edge):
		""" Expand edge into a path through chain """
		return sum(map(lambda X: self._extendEdgeIntoNet(X, edge), self._enumerateChain(edge)), [])

	def _extendEdge(self, edge):
		""" Return edge or equivalent path through subordinate nets """
		if edge.index == -1:
			# Simple adjacency
			return [edge]
		elif edge.start in self.netHistories[self.cactus.nodeNet(edge.start)].module.pseudotelomeres:
			# Going to an upper cycle (counter recursion)
			assert False
		elif self.cactus.chains2Nets[self.cactus.nodeChain(edge.start)] is None or len(self.cactus.chains2Nets[self.cactus.nodeChain(edge.start)]) == 0:
			# No sub-chain, end of recursion
			return [edge]
		else:
			# Recursion
			return [self._blockEdge(edge.start, edge.value)] + self._extendEdgeIntoChain(edge)

	def _unifyCycle(self, cycle):
		""" Returned a flattened version of the cycle """
		return Cycle(sum(map(lambda X: self._extendEdge(X), cycle), []))

	def _unifyEvent(self, event):
		""" Return a flattened version of the event """
		return Event(self._unifyCycle(event.cycle))

	def _unifyNet(self, net):
		""" Return flattened versions of the net's events """
		history = self.netHistories[net]
		module = history.module
		new = History(module)
		for event in history.events:
			if not _pseudotelomeric(event, module):
				newEvent = self._unifyEvent(event)
				newEvent.setRatio(event.cycle[0].value)
				self.slideIn_Event(event, newEvent)
				self.popEvent(event)
				new.absorbEvent(newEvent)
		return new


	def _unifyCactusHistory(self):
		""" Return flattened versions of the history's events """
		return dict((net, self._unifyNet(net)) for net in self.cactus.nets)

def main():
	graph = avg.randomNearEulerianGraph(10)
	C = cactus.Cactus(graph)
	NC = normalized.NormalizedCactus(C)
	BC = balanced.BalancedCactus(NC)
	OC = oriented.OrientedCactus(BC)
	H = cycleCover.initialHistory(OC)
	H.validate()
	print H
	c = H.rearrangementCost()
	FH = flattenGraph(H)
	print FH
	print FH.cactus.telomeres
	FH.validate()

if __name__ == "__main__":
	main()
