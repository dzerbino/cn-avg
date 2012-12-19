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

import math
import random
import scheduled
import numpy as np

import debug

""" Computing and detecting rearrangment costs plus any discrepancy in a history """

def _originalComponents_Node(data, node, module):
	if node.orientation:
		return data
	else:
		partner = module[node].partner
		if partner is node:
			return data

		componentLoops, nodeComponents = data
		newComponent = frozenset([node, partner])
		componentLoops[newComponent] = 0
		nodeComponents[node] = newComponent
		nodeComponents[partner] = newComponent
		return componentLoops, nodeComponents

def _originalComponents(module):
	return reduce(lambda X,Y: _originalComponents_Node(X, Y, module), module, (dict(), dict()))

def undoComponent_Node(data, node):
	componentLoops, nodeComponents = data
	if node in nodeComponents:
		component = nodeComponents[node]
		print 'DELETING', len(component), "\t".join(str(X.ID) for X in component)
		for node2 in component:
			del nodeComponents[node2]
		if component in componentLoops:
			del componentLoops[component]
	return componentLoops, nodeComponents

def undoComponents_Edge(data, edge):
	return reduce(undoComponent_Node, edge.nodes(), data)

def expandComponent(data, node, netHistory, activeBonds, createdBonds):
	component, L = data
	component.add(node)
	nodeFlow = netHistory.module[node]

	# Check for self loop
	if nodeFlow.selfLoops and createdBonds[netHistory.mappings.getBond(node, node)]:
		loops = 1
	else:
		loops = 0

	# Propagate component
	neighbours = [X for X in nodeFlow.edges if X not in component and activeBonds[netHistory.mappings.getBond(node, X)]]
	return reduce(lambda X, Y: expandComponent(X, Y, netHistory, activeBonds, createdBonds), neighbours, (component, L + loops))

def redoComponents_Node(data, node, netHistory, activeBonds, createdBonds):
	componentLoops, nodeComponents = data
	if node not in nodeComponents:
		component, loops = expandComponent((set(), 0), node, netHistory, activeBonds, createdBonds)
		if len(component) > 1 or loops > 0:
			print 'NEW', len(component)
			component = frozenset(component)
			componentLoops[component] = int(math.ceil(loops*0.5))
			for node2 in component:
				nodeComponents[node2] = component
	return componentLoops, nodeComponents

def redoComponents_Edge(data, edge, netHistory, activeBonds, createdBonds):
	return reduce(lambda X, Y: redoComponents_Node(X, Y, netHistory, activeBonds, createdBonds), edge.nodes(), data)

def updateComponents(netHistory, activeBonds, createdBonds, cycle, data):
	data = reduce(undoComponents_Edge, cycle, data)
	return reduce(lambda X, Y: redoComponents_Edge(X, Y, netHistory, activeBonds, createdBonds), cycle, data)

class ConstrainedHistory(scheduled.ScheduledHistory):
	######################################
	## Basic
	######################################
	def __copy__(self):
		new = ConstrainedHistory(self.cactus)
		new.copy(self)
		return new
		
	######################################
	## Find possible insertion points
	######################################
	def childInsertionPoints(self, event, new):
		ratio = new.ratio
		if ratio > event.ratio or event == new:
			return []
		else:
			ratios = [X.ratio for X in self.children[event]]
			if len(ratios) == 0:
				return [(event, None)]
			else:
				total = sum(ratios)

				insertionPoints = []
				if total + ratio <= event.ratio:
					insertionPoints.append((event, None))
				if ratio >= max(ratios):
					insertionPoints.append((event, new))
				inside = map(lambda X: (event, X), filter(lambda X: X.ratio <= ratio and total - X.ratio + ratio <= event.ratio, self.children[event]))
				return inside + insertionPoints

	def rootInsertionPoints(self, new):
		ratio = new.ratio
		total = sum(X.ratio for X in self.roots)
		
		insertionPoints = []
		if total <= 1:
			insertionPoints.append((None, None))
		if total - ratio <= ratio:
			insertionPoints.append((None, new))
		inside = map(lambda X: (None, X), filter(lambda X: X != new and X.ratio <= ratio and total - X.ratio <= 1, self.roots))
		return inside + insertionPoints

	def insertionPoints(self, event):
		return self.rootInsertionPoints(event) + sum(map(lambda X: self.childInsertionPoints(X, event), self.children), [])

	######################################
	## Find desperate insertion points
	######################################
	def childDesperateInsertionPoints(self, event, new, parent):
		ratio = new.ratio
		if event == parent:
			return []
		elif ratio > event.ratio or event == new:
			return []
		else:
			insideNodes = filter(lambda X: X.ratio <= ratio, self.children[event])
			if len(insideNodes) > 1:
				insideNode = max(insideNodes, key=lambda X: X.ratio)
				return [(event, insideNode)]
			else:
				return []

	def rootDesperateInsertionPoints(self, new):
		ratio = new.ratio
		insideNodes = filter(lambda X: X != new and X.ratio <= ratio, self.roots)
		if len(insideNodes) > 1:
			insideNode = max(insideNodes, key=lambda X: X.ratio)
			return [(None, insideNode)]
		else:
			return []

	def desperateInsertionPoints(self, event, parent):
		return self.rootDesperateInsertionPoints(event) + sum(map(lambda X: self.childDesperateInsertionPoints(X, event, parent), self.children), [])

	######################################
	## Schedule event
	######################################

	def scheduleEvent(self, event, parent=None):
		if len(self.parent) == 1:
			return
		candidateInsertions = self.insertionPoints(event)
		if len(candidateInsertions) != 0:
			insertion = random.choice(candidateInsertions)
			if insertion[1] != event:
				self.swapIn(event, insertion[0], insertion[1])
			else:
				self.swapUnder(event, insertion[0])
		else: 
			candidateInsertions = self.desperateInsertionPoints(event, parent)
			insertion = random.choice(candidateInsertions)
			self.swapIn(event, insertion[0], insertion[1])

	def absorbEvent(self, history, event):
		self.correctSchedulingErrors()
		super(ConstrainedHistory, self).absorbEvent(history, event)
		if event in self.parent:
			self.scheduleEvent(event)
			self.correctSchedulingErrors()

	######################################
	## Correct errors in scheduling
	######################################
	def childErrors(self, event):
		# NOTE: 0.1 is criteria for keeping an event in the default pipeline
		badBoys = [X for X in self.children[event] if X.ratio > event.ratio]
		if len(badBoys) > 0:
			return badBoys
		if (not debug.DEBUG) and sum(X.ratio for X in self.children[event] if X.ratio >= 0.1) > event.ratio:
			return [min(filter(lambda X: X.ratio >= 0.1, self.children[event]), key=lambda X: X.ratio)]
		elif sum(X.ratio for X in self.children[event]) > event.ratio:
			return [min(self.children[event], key=lambda X: X.ratio)]
		else:
			return []

	def rootErrors(self):
		if sum(X.ratio for X in self.roots) > 1:
			return list(self.roots)
		else:
			return []

	def schedulingErrors(self):
		# NOTE: 0.1 is criteria for keeping an event in the default pipeline
		if debug.DEBUG:
			return self.rootErrors() + sum(map(self.childErrors, (X for X in self.children.keys())), []) 
		else:
			return self.rootErrors() + sum(map(self.childErrors, (X for X in self.children.keys() if X.ratio > 0.1)), []) 

	def correctSchedulingError(self, event):
		parent = self.parent[event]
		self.swapOut(event)
		self.scheduleEvent(event, parent)

	def correctSchedulingErrors(self):
		errors = self.schedulingErrors()
		counter = 0
		while len(errors) != 0:
			error = random.choice(errors)
			self.correctSchedulingError(error)
			errors = self.schedulingErrors()
			counter += 1
			if counter > 1000:
				# Getting desperate, just start over with linear tree
				self.forceOrder(sorted(self.parent.keys(), key=lambda X: -X.ratio))

	######################################
	## Computing rearrangement cost
	######################################

	def expandComponent(self, data, node, netHistory, previousBonds, createdBonds):
		visited, edges, ancestralEdges, overlaps = data
		if node in visited:
			return data
		else:
			visited.add(node.ID)
			nodeFlow = netHistory.module[node]

			# Check for self loop
			if nodeFlow.selfLoops and createdBonds[netHistory.mappings.getBond(node, node)]:
				selfEdge = 1
			else:
				selfEdge = 0

			# Propagate component
			oldNeighbours = set(X for X in nodeFlow.edges if X.ID not in visited and previousBonds[netHistory.mappings.getBond(node, X)])
			a = len(oldNeighbours)
			o = max(a - 1, 0)
			newNeighbours = set(X for X in nodeFlow.edges if X.ID not in visited and createdBonds[netHistory.mappings.getBond(node, X)])
			neighbours = list(oldNeighbours | newNeighbours)
			newEdges = len(neighbours)

			return reduce(lambda X, Y: self.expandComponent(X, Y, netHistory, previousBonds, createdBonds), neighbours, (visited, edges + selfEdge + newEdges, ancestralEdges + a, overlaps + o))
		
	def countComponents_Node(self, data, node, netHistory, previousBonds, createdBonds, cycleLength):
		visited, cost = data
		if node.ID in visited:
			return data
		else:
			componentNodes, edges, ancestralEdges, overlaps = self.expandComponent((set(), 0, 0, 0), node, netHistory, previousBonds, createdBonds)
			if edges == len(componentNodes) and edges == cycleLength:
				componentLoops = 0
			else:
				componentLoops = int(math.ceil((edges - len(componentNodes) + 1)*0.5))

			#print 'COMPLEXITY', ancestralEdges, componentLoops, overlaps, max(ancestralEdges - 1 + componentLoops - overlaps, 0)
			return visited | componentNodes, cost + max(ancestralEdges - 1 + componentLoops - overlaps, 0)

	def countComponents(self, cycle, netHistory, previousBonds, createdBonds):
		return reduce(lambda X,Y: self.countComponents_Node(X, Y, netHistory, previousBonds, createdBonds, len(cycle)), (X.start for X in cycle), (set(), 0))[1]

	def rearrangementCost_Net(self, net):
		netHistory = self.netHistories[net]
		originalVector = netHistory.originalVector()
		bond = netHistory.bondIndices()
		segment = (bond == False)
		queue = [(X, originalVector) for X in self.roots]
		total = 0

		# I wish I could do this in proper recursion FP style, but Python is not very good with deep recursions
		while len(queue) > 0:
			event, previousVector = queue.pop(0)

			if (not debug.DEBUG) and event.ratio < 0.1:
				continue

			if event in netHistory.eventIndex:
				if event.cycle.value > 0:
					eventVector = -netHistory.eventVector(event) 
				else:
					eventVector = netHistory.eventVector(event) 
				newVector = previousVector + eventVector
				previousBonds = bond & (previousVector < 0)
				createdBonds = bond & (newVector < 0) & (previousVector >= 0)
				cost = self.countComponents(event.cycle, netHistory, previousBonds, createdBonds)
				segmentErrors = (eventVector != 0) & segment & (previousVector <= 0) 
				bondDestructions = (eventVector != 0) & bond & (newVector > 0)
				total += cost + 2 * np.sum(segmentErrors | bondDestructions)
				#print 'COMPLEXITY', n, C, L, O, np.sum(segmentErrors | bondDestructions), max(n - C + L - O + 2 * np.sum(segmentErrors | bondDestructions), 0) 
			else:
				newVector = previousVector
			queue.extend((X, newVector) for X in self.children[event])
		return total

	def rearrangementCost(self):
		# Yeah I know, ugly side effect, but avoids re-computing stuff.
		# You only compute complexity on the final version of a history. 
		# If you need to modify it after that, create a copy then carry on.
		if self.complexity is None:
			self.complexity = sum(map(self.rearrangementCost_Net, self.netHistories))
		return self.complexity
