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

"""History with look up table to speed up the search for overlapping flow changes"""

import copy
import random

import debug

from cnavg.flows.edge import Edge
import cnavg.history.history as history

""" Definition of overlap graphs """

NON_STUB_COEFF = 1e-2
""" Probability mass of an edge not-connecting to the stub wrt the weight of stub edges """

#########################################
## Overlaps
#########################################

class Overlap:
	""" Overlap between two events """
	def __init__(self, localEvent, localCut, remoteEvent, remoteCut):
		self.localEvent = localEvent
		self.localCut = localCut
		self.remoteEvent = remoteEvent
		self.remoteCut = remoteCut

	def __str__(self):
		return "%i:%i -> %i:%i" % (self.localEvent, self.localCut, self.remoteEvent, self.remoteCut)

def outAdjacencyOverlaps(eventA, adjacencyOverlapsA, eventB, adjacencyOverlapsB):
	i = random.randrange(len(adjacencyOverlapsA))
	j = random.randrange(len(adjacencyOverlapsB))
	edgeA = eventA.cycle[adjacencyOverlapsA[i]]
	edgeB = eventB.cycle[adjacencyOverlapsB[j]]

	return Overlap(eventA, adjacencyOverlapsA[i], eventB, adjacencyOverlapsB[j])

def selfAdjacencyOverlaps(event, adjacencyOverlaps, stub):
	if stub:
		return [(1, Overlap(event, adjacencyOverlaps[i], event, adjacencyOverlaps[j]))
			for i in range(len(adjacencyOverlaps) - 1)
			for j in range(i+1, len(adjacencyOverlaps))
			if event.cycle[adjacencyOverlaps[i]].start == event.cycle[adjacencyOverlaps[j]].start]
	else:
		return [(NON_STUB_COEFF, Overlap(event, adjacencyOverlaps[i], event, adjacencyOverlaps[j]))
			for i in range(len(adjacencyOverlaps) - 1)
			for j in range(i+1, len(adjacencyOverlaps))
			if event.cycle[adjacencyOverlaps[i]].start == event.cycle[adjacencyOverlaps[j]].start]

def square(X):
	return X * X

##########################################
## Overlap History
##########################################

def adjacencyIndex(edge):
	return (min(edge.start, edge.finish), max(edge.start, edge.finish), edge.index)

def addEdgeToOverlapTable(overlapTable, edgeIndex, event):
	edge = event.cycle[edgeIndex]
	index = adjacencyIndex(edge)
	if index not in overlapTable:
	    overlapTable[index] = dict()
	if event not in overlapTable[index]:
	    overlapTable[index][event] = []

	overlapTable[index][event].append(edgeIndex)
	return overlapTable

def addEventToOverlapTable(overlapTable, event):
	return reduce(lambda X, Y: addEdgeToOverlapTable(X, Y, event), range(len(event.cycle)), overlapTable)

def removeEdgeFromOverlapTable(overlapTable, edge, event):
	index = adjacencyIndex(edge)
	if event in overlapTable[index]:
	    del overlapTable[index][event]
	return overlapTable

def removeEventCycleFromOverlapTable(overlapTable, event):
	return reduce(lambda X, Y: removeEdgeFromOverlapTable(X, Y, event), event.cycle, overlapTable)

def copyAdjacencyOverlaps(adjacencyOverlaps, eventsMappings):
	return dict((eventsMappings[X], copy.copy(adjacencyOverlaps[X])) for X in adjacencyOverlaps)

def _findLotteryWinner(target, weightedValues):
	index = 0
	weight, value = weightedValues[index]
	while index < len(weightedValues) and weight <= target:
		weight, value = weightedValues[index]
		target -= weight
		index += 1
	return value

def _weightedRandomSelection(weightedList):
	total = sum(X[0] for X in weightedList)
	target = random.uniform(0, total)
	return _findLotteryWinner(target, weightedList)

class OverlapHistory(history.History):
	###########################################
	## Basics
	###########################################
	def __init__(self, module):
		super(OverlapHistory, self).__init__(module)
		self.overlapTable = dict()

	def __copy__(self):
		new = OverlapHistory(self.module)
		new.copy(self)
		return new

	def copy(self, other):
		super(OverlapHistory, self).copy(other)
		eventsMapping = dict(zip(other.events, self.events))
		self.overlapTable = dict((X, copyAdjacencyOverlaps(other.overlapTable[X], eventsMapping)) for X in other.overlapTable)

	def embalm(self):
		self.overlapTable = None

	###########################################
	## Output
	###########################################
	def _overlapString(self, edge):
		return str(edge) + "\n" + "\n".join(str(event.cycle[index]) for event in self.overlapTable[edge] for index in self.overlapTable[edge][event])


	def _overlapStrings(self):
		return '\n>>>>>>>>>>>>>>>>>\n'.join(map(self._overlapString, self.overlapTable))

	def __str__(self):
		return "\n".join([super(OverlapHistory, self).__str__(),
				  'OVERLAPTABLE',
				  self._overlapStrings()])

	###########################################
	## Inherited API
	###########################################
	def pop(self, event):
		super(OverlapHistory, self).pop(event)
		self.overlapTable = removeEventCycleFromOverlapTable(self.overlapTable, event)
		return event

	def absorbEvent(self, event):
		super(OverlapHistory, self).absorbEvent(event)
		self.overlapTable = addEventToOverlapTable(self.overlapTable, event)
		return self

	###########################################
	## Core functionality
	###########################################
	# NOTE: this function is the heart of the heurisitic. It decides what gets merged and what gets split
	def adjacencyOverlaps(self, adjacency, indirect):
		# Do not take into account overlaps in the stub-stub edge as it represents no material evidence
		if adjacency[0].chr == "None" and adjacency[1].chr == "None":
			return []
		adjacencyTable = self.overlapTable[adjacency]
		listLengths = [len(adjacencyTable[X]) for X in adjacencyTable]
		totalLength = sum(listLengths)
		totalOverlaps = square(totalLength) - len(listLengths)
		internalOverlaps = sum(square(X) for X in listLengths) - len(listLengths)
		externalOverlaps = totalOverlaps - internalOverlaps

		stub = (adjacency[0].chr == "None" or adjacency[1].chr == "None")

		if indirect:
			selfOverlaps = []
		else:
			selfOverlaps = sum(map(lambda X: selfAdjacencyOverlaps(X, adjacencyTable[X], stub), adjacencyTable), [])

		events = adjacencyTable.keys()
		if len(events) > 1:
			if debug.DEBUG:
				candidateBs = filter(lambda X: X not in self.untouchables, events)
			else:
				candidateBs = filter(lambda X: X not in self.untouchables and X.ratio > debug.RATIO_CUTOFF, events)
			if len(candidateBs) > 0:
				B = random.choice(candidateBs)
				A = random.choice(filter(lambda X: X != B, events))
				if stub:
					overlaps = [(externalOverlaps, outAdjacencyOverlaps(A, adjacencyTable[A], B, adjacencyTable[B]))]
				else:
					overlaps = [(NON_STUB_COEFF * externalOverlaps, outAdjacencyOverlaps(A, adjacencyTable[A], B, adjacencyTable[B]))]
			else:
				overlaps = []
		else:
			overlaps = []

		return selfOverlaps + overlaps

	def overlaps(self, indirect):
		return sum([self.adjacencyOverlaps(X, indirect) for X in self.overlapTable],[])

	def overlap(self, indirect):
		overlaps = self.overlaps(indirect)
		if len(overlaps) == 0:
			return None
		else:
			return _weightedRandomSelection(overlaps)

	###########################################
	## Stats
	###########################################
	def density(self):
		if len(self.overlapTable) > 0:
			if debug.DEBUG:
				return sum(max(len(filter(lambda X: X not in self.untouchables, self.overlapTable[X].keys())) - 1, 0) for X in self.overlapTable if X[0].chr != "None" or X[1].chr != "None" ) / float(len(self.overlapTable))
			else:
				return sum(max(len(filter(lambda X: X not in self.untouchables and X.ratio > debug.RATIO_CUTOFF, self.overlapTable[X].keys())) - 1, 0) for X in self.overlapTable if X[0].chr != "None" or X[1].chr != "None" ) / float(len(self.overlapTable))
		else:
			return 0

	###########################################
	## Validation
	###########################################
	def _untouchablesUntouched_index(self, node, index):
		twin = self.module[node].twin
		edgeIndex = (min(twin, node), max(twin, node), index)
		if edgeIndex not in self.overlapTable:
			return True
		assert len(self.overlapTable[edgeIndex].keys()) <= 1
		if len(self.overlapTable[edgeIndex].keys()) == 1:
			assert self.overlapTable[edgeIndex].keys()[0] in self.untouchables
		return True

	def _untouchablesUntouched(self):
		if len(self.module.pseudotelomeres) > 0:
			assert all(self._untouchablesUntouched_index(PT, X) for PT in self.module.pseudotelomeres for X in range(len(self.module.segments[PT])))
		return True

	def _validateSegments(self):
		for aIndex in self.overlapTable:
			edgeTable = self.overlapTable[aIndex]
			node, twin, index = aIndex
			if index == -1:
				continue
			total = sum(event.cycle[edgeIndex].value for event in edgeTable for edgeIndex in edgeTable[event])

			if total == 0 and abs(self.module.segments[node][index]) > 1e-1:
				#print self.module
				#print self
				print node, twin, index
				print total
				print self.module.segments[node][index]
			assert not (total == 0 and abs(self.module.segments[node][index]) > 1e-1)
		return True

	def _validateInheritedCycles(self, node):
		module = self.module
		twin = module[node].twin
		for segment in range(len(self.module.segments[node])):
			if self.module.segments[node][segment] < 1e-2:
				continue
			aIndex = (min(node, twin), max(node, twin), segment)
			if aIndex not in self.overlapTable or len(self.overlapTable[aIndex].keys()) != 1:
				print self.module
				print self.module.segments[node][segment]
				print self
				print len(self.overlapTable[aIndex].keys())
			assert len(self.overlapTable[aIndex].keys()) == 1

	def validate(self):
		super(OverlapHistory, self).validate()
		assert self._validateSegments()
		assert self._untouchablesUntouched()
		for x in self.module.pseudotelomeres:
			self._validateInheritedCycles(x)
		return True
