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

"""Keeping track of the tree of events in a history"""
import sys
import copy

import history

class ScheduledHistory(history.CactusHistory):
	#################################
	## Basics
	#################################
	def __init__(self, module):
		super(ScheduledHistory, self).__init__(module)
		self.initialize()

	def initialize(self):
		self.parent = dict()
		self.children = dict()
		self.ancestors = dict()
		self.descendants = dict()
		self.roots = set()

	def __copy__(self):
		new = ScheduledHistory(self.cactus)
		new.copy(self)
		return new

	def copy(self, other):
		super(ScheduledHistory, self).copy(other)
		self.parent = copy.copy(other.parent) 
		self.children = dict((X, copy.copy(other.children[X])) for X in other.children) 
		self.ancestors = dict((X, copy.copy(other.ancestors[X])) for X in other.ancestors) 
		self.descendants = dict((X, copy.copy(other.descendants[X])) for X in other.descendants) 
		self.roots = copy.copy(other.roots) 

	#################################
	## Ensuring tree-ness
	#################################
	def addAncestry(self, parent, event, singleParent=False, singleChild=False):
		descendants = set([event])
		if not singleChild:
			descendants |= self.descendants[event]
		ancestors = set([parent])
		if not singleParent:
			ancestors |= self.ancestors[parent]
		for descendant in descendants:
			assert descendant in self.ancestors, descendant == parent
			self.ancestors[descendant] |= ancestors
		for ancestor in ancestors:
			assert ancestor in self.descendants, ancestor == parent
			self.descendants[ancestor] |= descendants

	def addDescent(self, parent, child, singleParent=False, singleChild=False):
		if child is None:
			return
		elif parent is None:
			self.rootify(child)
			return
		if child in self.ancestors[parent]:
			print child.ratio
			print parent.ratio
			print self.descentTreeString()
			assert False
		self.parent[child] = parent
		self.children[parent].add(child)
		if parent not in self.ancestors[child]:
			self.addAncestry(parent, child, singleParent=singleParent, singleChild=singleChild)
		if child in self.roots:
			self.roots.remove(child)

	#################################
	## Moves
	#################################
	def swapIn(self, event, parent, child):
		self.addDescent(parent, event, singleChild=True)
		self.addDescent(event, child, singleParent=True)
		if parent is not None and child is not None:
			self.children[parent].remove(child)

	def swapUnder(self, event, parent):
		if parent is not None:
			todo = [X for X in self.children[parent]]
			if event not in self.children[parent]:
				self.addDescent(parent, event, singleChild=True)
			for child in todo:
				if child != event:
					self.addDescent(event, child, singleParent=True)
					self.children[parent].remove(child)
		else:
			todo = [X for X in self.roots if X != event]
			for child in todo:
				self.addDescent(event, child, singleParent=True)
			self.rootify(event)

	def rootify(self, event):
		self.parent[event] = None
		self.ancestors[event] = set()
		self.roots.add(event)

	def spacify(self, event):
		self.rootify(event)
		self.descendants[event] = set()
		self.children[event] = set()

	def absorbEvent(self, history, event):
		super(ScheduledHistory, self).absorbEvent(history, event)
		if event not in history.untouchables and event in history.events: # Event can be bumped out at Euclidian stage
			self.spacify(event)

	def swapOut(self, event):
		for descendant in self.descendants[event]:
			self.ancestors[descendant].remove(event)
		for ancestor in self.ancestors[event]:
			self.descendants[ancestor].remove(event)
		parent = self.parent[event]
		if parent is not None:
			self.children[parent].remove(event)
		for child in self.children[event]:
			self.addDescent(parent, child)
		self.spacify(event)

	def popEvent(self, event):
		# Untouchable events are not tracked
		if event in self.parent:
			self.swapOut(event)
			del self.parent[event]
			del self.children[event]
			del self.ancestors[event]
			del self.descendants[event]
			if event in self.roots:
				self.roots.remove(event)

	def popHistory(self, history):
		for event in history.events:
			self.popEvent(event)

	def pop(self, history, event):
		super(ScheduledHistory, self).pop(history, event)
		self.popEvent(event)

	def getTopEvent(self, history, event):
		if event in history.untouchables:
			index = history.untouchables[event]
			topChain = self.cactus.headChain[history.module.net] 
			topNet = self.cactus.headNet[topChain]
			topHistory = self.netHistories[topNet]

			topEdge, topEvent = self.chainCNVs[topChain][index]
			return self.getTopEvent(topHistory, topEvent)
		else:
			return event

	def getAncestors(self, history, event):
		return self.ancestors[self.getTopEvent(history, event)]

	def slideIn_Event(self, oldEvent, newEvent):
		if oldEvent in self.parent:
			self.spacify(newEvent)
			self.swapIn(newEvent, self.parent[oldEvent], oldEvent)

	def slideIn(self, oldHistory, newHistory):
		for index in range(len(oldHistory.events)):
			if oldHistory.events[index] not in self.parent:
				assert oldHistory.events[index] in oldHistory.untouchables
				assert newHistory.events[index] in newHistory.untouchables
			self.slideIn_Event(oldHistory.events[index], newHistory.events[index])

	def forceOrder(self, iter):
		""" Emergency re-ordering function """
		self.initialize()
		last = None
		for head in iter:
			self.spacify(head)
			if last is not None:
				self.addDescent(last, head)
			last = head

	#################################
	## Statistics
	#################################
	def forks(self):
		return filter(lambda X: len(self.children[X]) > 1, self.children)

	def countForks(self):
		return len(self.forks())

	def stats(self):
		return super(ScheduledHistory, self).stats() + "\t".join(map(str, [self.countForks()]))

	#################################
	## Display
	#################################
	def descentTreeString_Branch(self, parent, child):
		return "%i -> %i" % (id(parent), id(child))

	def descentTreeString_Node(self, event):
		return "%i [label=%f]" % (id(event), event.ratio)

	def descentTreeString_Event(self, event):
		return "\n".join([self.descentTreeString_Branch(event, X) for X in self.children[event]] + [self.descentTreeString_Node(event)])
	
	def descentTreeString(self):
		return "\n".join([self.descentTreeString_Event(X) for X in self.children if len(self.children[X]) > 0] + [self.descentTreeString_Node(X) for X in self.children if len(self.children[X]) == 0])

	def __str__(self):
		return "\n".join([super(ScheduledHistory, self).__str__(), self.descentTreeString()])

	#################################
	## Newick
	#################################
	def newick_remainderLeaf(self, event):
		return "%i:%f" % (id(event), event.ratio)

	def newick_node(self, event):
		if len(self.children[event]) > 0:
			fork = "(" + ",".join([self.newick_node(child) for child in self.children[event]] + [self.newick_remainderLeaf(event)]) + ")"
		else:
			fork = "(" + ",".join([self.newick_remainderLeaf(event), "n" + self.newick_remainderLeaf(event)]) + ")"

		if self.parent[event] is None:
			length = str(1 - event.ratio)
		else:
			length = str(self.parent[event].ratio - event.ratio)

		return ":".join([fork, length])

	def newick(self):
		if len(self.roots) == 1:
			return self.newick_node(list(self.roots)[0]) + ";"
		elif len(self.roots) > 1:
			return "(" + ",".join([self.newick_node(event) for event in self.roots]) + ");"

		
	#################################
	## Validate
	#################################
	def validateEvent(self, history, event):
		if event is None:
			return True
		if event not in self.parent:
			assert event in history.untouchables, id(event)
			return True

		return True

	def validateHistory(self, history):
		assert all(self.validateEvent(history, X) for X in history.events)
		return True

	def validateTree(self, event):
		assert event in self.parent
		assert event in self.children
		assert event in self.ancestors
		assert event in self.descendants
		if self.parent[event] is not None:
			assert event in self.children[self.parent[event]]
		else:
			assert event in self.roots
		assert event not in self.ancestors[event]
		parent = self.parent[event]
		while parent is not None:
			assert parent in self.ancestors[event]
			parent = self.parent[parent]

		for child in self.children[event]:
			assert child in self.parent and self.parent[child] == event

		for descendant in self.descendants[event]:
			assert descendant in self.ancestors and event in self.ancestors[descendant]

		for ancestor in self.ancestors[event]:
			assert ancestor in self.descendants and event in self.descendants[ancestor]
		return True

	def validate(self):
		assert all(self.validateTree(event) for event in self.parent)
		assert super(ScheduledHistory, self).validate()
		assert all(self.validateHistory(history) for history in self.netHistories.values())
		events  = set(E for H in self.netHistories.values() for E in H.events)
		for X in self.parent:
			if X not in events:
				print id(X)
		assert all(X in events for X in self.parent), str(self.parent)
		return True

#########################################
## Unit test
#########################################
def main():
	pass

if __name__ == '__main__':
	main()
