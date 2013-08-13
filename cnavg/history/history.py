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

"""Flow history"""
import copy

from cnavg.flows.flows import Event

#############################################
## History
#############################################

def median(list):
	return sorted(list)[int(len(list)/2)]

class History(object):
	#############################################
	## Basic
	#############################################
	def __init__(self, module):
		self.module = module
		self.events = list() 
		self.eventIndex = dict()
		self.untouchables = dict()

	def copy(self, other):
		self.events = map(copy.copy, other.events)
		self.eventIndex = dict((X[1], X[0]) for X in enumerate(self.events))
		self.untouchables = dict((self.events[other.eventIndex[X]], other.untouchables[X]) for X in other.untouchables)

	def __copy__(self):
		new = History(self.module)
		new.copy(self)
		return new

	def absorbEvent(self, event):
		""" Add event to history """
		self.eventIndex[event] = len(self.events)
		self.events.append(event)
		index = self.isUntouchable(event)
		if index >= 0:
			self.untouchables[event] = index

	def pop(self, event):
		""" Remove event from history """
		self.events.remove(event)
		if event in self.untouchables:
			del self.untouchables[event]
		index = self.eventIndex[event]
		for sibling in self.events:
			if self.eventIndex[sibling] > index:
				self.eventIndex[sibling] -= 1
		del self.eventIndex[event]

	#############################################
	## The untouchables 
	#############################################

	def _isPseudoTelomeric(self, edge):
		if edge.index >= 0 and edge.start in self.module.pseudotelomeres and edge.finish in self.module.pseudotelomeres:
			return edge.index
		else:
			return -1	

	def _hasAPseudoTelomericEdge(self, cycle):
		return max(self._isPseudoTelomeric(X) for X in cycle)
	
	def isUntouchable(self, event):
		""" Tests whether the event is actually the extension of an event on a higher level net onto the current net """
		return self._hasAPseudoTelomericEdge(event.cycle)

	#############################################
	## Validation
	#############################################

	def _numberingIsCorrect(self):
		assert all(event == self.events[self.eventIndex[event]] for event in self.events)
		return True

	def validate(self):
		""" Validation function """
		assert self._numberingIsCorrect()
		assert all(map(Event.validate, self.events))
		return True
		
	#############################################
	## Stats
	#############################################

	def getCNVs(self, graph):
		return reduce(lambda X, Y: Y.getCNVs(X, graph), self.events, dict())

	def eventLengths(self):
		return map(Event.length, self.events)

	def medianLength(self):
		return median(self.eventLengths())

	def maxLength(self):
		if len(self.events) > 0:
			return max(self.eventLengths())
		else:
			return 0

	def length(self):
		""" Number of events """
		return len(self.events)

	#############################################
	## Display
	#############################################
	def __str__(self):
		return "\n".join(["HISTORY"] + map(str, self.events))

	def dot(self):
		""" GraphViz representation """
		return "\n".join(X.dot() for X in self.events)

	def braneyText(self, historyID, netID, ordering, complexity):
		""" Braney representation """
		return "\n".join(X[1].braneyText(historyID, netID, str(X[0]), ordering, complexity) for X in enumerate(set(self.events) - set(self.untouchables)))

	def simplifyStubsAndTrivials(self, cactusHistory):
		new = History(self.module)
		for event in self.events:
			newEvent = event.simplifyStubsAndTrivials(cactusHistory.cactus)
			if newEvent is not None:
				new.absorbEvent(newEvent)
				cactusHistory.slideIn_Event(event, newEvent)
				cactusHistory.popEvent(event)
		return new

	def removeLowRatioEvents(self, ratio, cactusHistory):
		""" Filters out events with low ratios """
		new = History(self.module)
		for event in self.events:
			if abs(event.ratio) > ratio:
				new.absorbEvent(event)
			else:
				cactusHistory.popEvent(event)
		return new

	def selectForRegion(self, nodes):
		""" Selects event which do not contain nodes """
		return filter(lambda X: X.doesNotContainNodes(nodes), self.events)

#############################################
## Cactus History
#############################################
class CactusHistory(object):
	""" CN-AVG history across an entire Cactus graph """
	#############################################
	## Basic
	#############################################
	def __init__(self, cactus):
		self.cactus = cactus
		self.netHistories = dict()
		self.chainCNVs = dict()
		self.complexity = None

	def __copy__(self):
		result = CactusHistory(self.cactus)
		result.copy(self)
		return result

	def copy(self, origin):
		self.netHistories = copy.copy(origin.netHistories)
		self.chainCNVs = copy.copy(origin.chainCNVs)
		self.complexity = None

	#############################################
	## Main Operations
	#############################################
	def update(self, net, history):
		""" Assign local history to a net """
		if net in self.netHistories:
			popped = self.netHistories[net]
			if self.netHistories[net] == history:
				return
			else:
				self.popHistory(self.netHistories[net])
		self.netHistories[net] = history

	def updateCNVs(self, net, history):
		""" Update CNVs assigned to a chain given a new history on a given net """
		newChainCNVs = history.getCNVs(self.cactus)
		for x in self.cactus.nets2Chains[net]:
			if x in newChainCNVs:
				self.chainCNVs[x] = sorted(newChainCNVs[x], key=lambda X: X[0].value)
			else:
				self.chainCNVs[x] = []

	def absorbEvent(self, history, event):
		""" Add event to local history """
		history.absorbEvent(self, event)

	def pop(self, history, event):
		""" Remove event from local history """
		history.pop(event)

	def popHistory(self, history):
		pass

	def embalmDeadHistories(self, other):
		""" Embalm histories which do not need to store their data """
		old = set(self.netHistories.values())
		new = set(other.netHistories.values())
		dead = old - new
		for X in dead:
			X.embalm()

	#############################################
	## Stats
	#############################################
	def netHistoryCosts(self):
		return map(lambda X: X.rearrangementCost(), self.netHistories.values())

	def rearrangementCost(self):
		if self.complexity is None:
			self.complexity = sum(self.netHistoryCosts())
		print self
		print self.complexity
		return self.complexity

	def _validateSegmentsLong(self, module, node):
		twin = module[node].twin
		pair = set([node, twin])

		# Ignore stubs
		if node not in self.cactus:
			return True

		if self.cactus.chains2Nets[self.cactus.nodeChain(node)] is not None and len(self.cactus.chains2Nets[self.cactus.nodeChain(node)]) > 0:
			return True

		if node in module.telomeres or node in module.pseudotelomeres:
			return True

		for haplotype in range(len(module.segments[node])):
			edges = [edge for history in self.netHistories.values() for event in history.events for edge in event.cycle if set([edge.start, edge.finish]) == pair and edge.index == haplotype]
			total = sum(edge.value for edge in edges) 

			if abs(total + module.segments[node][haplotype]) > 1e-1:
				#print self.module
				#print self
				print node, twin, haplotype
				print total
				print len(edges)
				print [edge.value for edge in edges]
				print [id(edge) for edge in edges]
				print module.segments[node][haplotype]
				assert False
		return True

	def validate(self):
		assert all(self._validateSegmentsLong(history.module, node) for history in self.netHistories.values() for node in history.module) 
		assert all(X.validate() for X in self.netHistories.values())
		return True

	def length(self):
		return sum(X.length() for X in self.netHistories.values())

	def maxLength(self):
		if len(self.netHistories) > 0:
			return max(X.maxLength() for X in self.netHistories.values())
		else:
			return 0

	#############################################
	## Display
	#############################################
	def __str__(self):
		return '\n'.join(['CACTUSHISTORY'] + map(str, self.netHistories.values()))

	def medianLength(self):
		return median(sum((X.eventLengths() for X in self.netHistories.values()),[]))

	def stats(self):
		return "\t".join(map(str, [self.rearrangementCost(), self.errorCost(), self.length(), self.medianLength(), self.maxLength()]))

	def dot(self):
		""" GraphViz representation """
		return "\n".join(X.dot() for X in self.netHistories.values())

	def removeLowRatioEvents(self, cutoff):
		""" Filter out events which do not pass a ratio cutoff """
		new = copy.copy(self)
		new.netHistories = dict((X, self.netHistories[X].removeLowRatioEvents(cutoff, new)) for X in self.netHistories)
		return new

	def simplifyStubsAndTrivials(self):
		new = copy.copy(self)
		new.netHistories = dict((X, self.netHistories[X].simplifyStubsAndTrivials(new)) for X in self.netHistories)
		return new

	def selectForRegion(self, chr, start, end):
		""" Select events which affect a given region """
		new = copy.copy(self)
		nodeFlows = filter(lambda X: X.inRegion(chr, start, end), new.cactus.values())
		nodes = [X.node for X in nodeFlows]
		for net in new.netHistories:
			for event in new.netHistories[net].selectForRegion(nodes):
				new.pop(self.netHistories[net], event)
		return new
