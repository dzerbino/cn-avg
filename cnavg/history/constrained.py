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

""" Computing and detecting rearrangment costs plus any discrepancy in a history """

import math
import random
import scheduled
import numpy as np
import copy

import debug

def sumTrios(list):
	return reduce(lambda S, X: (S[0] + X[0], S[1] + X[1], S[2] + X[2]), list, (0, 0, 0))

class ConstrainedHistory(scheduled.ScheduledHistory):
        ######################################
        ## Basic
        ######################################
	def __init__(self, module):
		super(ConstrainedHistory, self).__init__(module)
		self.eventCosts = None

        def __copy__(self):
                new = ConstrainedHistory(self.cactus)
                new.copy(self)
                return new

	def copy(self, other):
		super(ConstrainedHistory, self).copy(other)
		self.eventCosts = copy.copy(other.eventCosts)
                
        ######################################
        ## Find possible insertion points
        ######################################
        def childInsertionPoints(self, event, new):
                ratio = new.ratio

		if debug.INTEGER_HISTORY:
			if event == new:
				return []
			elif len(self.children[event]) == 0:
				return [(event, None)]
			else:
				return [(event, new)]

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
                
		if debug.INTEGER_HISTORY:
			if len(self.roots) == 0:
				return [(None, None)]
			else:
				return [(None, new)]

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
                if (not debug.DEBUG) and sum(X.ratio for X in self.children[event] if X.ratio >= debug.RATIO_CUTOFF) > event.ratio:
                        return [min(filter(lambda X: X.ratio >= debug.RATIO_CUTOFF, self.children[event]), key=lambda X: X.ratio)]
                elif debug.DEBUG and sum(X.ratio for X in self.children[event]) > event.ratio:
                        return [min(self.children[event], key=lambda X: X.ratio)]
                else:
                        return []

	def parentErrors(self, event):
		if self.parent[event] is not None and event.ratio > self.parent[event].ratio:
			return [event]
		else:
			return []

        def rootErrors(self):
                if sum(X.ratio for X in self.roots) > 1:
                        return list(self.roots)
                else:
                        return []

        def schedulingErrors(self):
                if debug.DEBUG:
			validNodes = self.children.keys()
                else:
			validNodes = filter(lambda X: X.ratio > debug.RATIO_CUTOFF, self.children)
		return self.rootErrors() + sum(map(self.childErrors, validNodes), []) + sum(map(self.parentErrors, validNodes), [])

        def correctSchedulingError(self, event):
                parent = self.parent[event]
                self.swapOut(event)
                self.scheduleEvent(event, parent)

        def correctSchedulingErrors(self):
		if debug.INTEGER_HISTORY:
			return
                errors = self.schedulingErrors()
                counter = 0
                while len(errors) > 0:
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
                visited, ancestralEdges, overlaps = data
		if node.ID in visited:
			return data
		nodeFlow = netHistory.module[node]

		# Propagate component
		newNeighbours = []
		ancestralNeighbours = 0
		newAncestralEdges = 0
		for partner in nodeFlow.edges:
			if previousBonds[netHistory.mappings.getBond(node, partner)]:
				if partner != node:
					ancestralNeighbours += 1
				else:
					ancestralNeighbours += 2
					
				if partner.ID not in visited:
					newAncestralEdges += 1
					newNeighbours.append(partner)
			elif createdBonds[netHistory.mappings.getBond(node, partner)] and partner.ID not in visited:
				newNeighbours.append(partner)

		o = max(ancestralNeighbours - 1, 0)

		visited.add(node.ID)
		return reduce(lambda X, Y: self.expandComponent(X, Y, netHistory, previousBonds, createdBonds), newNeighbours, (visited, ancestralEdges + newAncestralEdges, overlaps + o))
                
        def computeLower_Node(self, data, node, netHistory, previousBonds, createdBonds, cycleLength):
                visited, cost = data
                if node.ID in visited:
                        return data
                else:
                        componentNodes, ancestralEdges, overlaps = self.expandComponent((set(), 0, 0), node, netHistory, previousBonds, createdBonds)
                        #return visited | componentNodes, cost + max(math.ceil(len(componentNodes)/2) - 1 - math.ceil(overlaps/2), 0)
                        return visited | componentNodes, cost + max(len(componentNodes) - ancestralEdges - 1, 0)

        def computeLower(self, cycle, netHistory, previousBonds, createdBonds):
                return reduce(lambda X,Y: self.computeLower_Node(X, Y, netHistory, previousBonds, createdBonds, len(cycle)), (X.start for X in cycle), (set(), 0))[1]

	def sumIncidents(self, netHistory, node, vector):
		return sum(lambda X: vector[netHistory.mappings.getBond(node, X)], netHistory.module[node].edges) 

	def computeImbalance_Segment(self, netHistory, node, twin, duplication, denovo):
		Dnode = self.sumIncidents(netHistory, node, duplication)
		Dtwin = self.sumIncidents(netHistory, node, duplication)
		Cnode = self.sumIncidents(netHistory, node, denovo)
		Ctwin = self.sumIncidents(netHistory, node, denovo)
		return max([Dnode - Dtwin - Ctwin, Dtwin - Dnode - Cnode, 0])

	def computeImbalance_Edge(self, data, edge, netHistory, duplication, denovo):
		visited, cost = data
		node = edge.start
		twin = netHistory.module[node].twin
		if node.ID in visited or twin in visited:
			return data
		else:
			return visited | set([node]), cost + self.computeImbalance_Segment(netHistory, node, twin, duplication, denovo)

	def computeImbalances(self, cycle, netHistory, duplication, denovo):
		return reduce(lambda X, Y: self.computeImbalance_Edge(X, Y, netHistory, duplication, denovo), cycle, (set(), 0))[1]

        def rearrangementCost_Net(self, net):
                netHistory = self.netHistories[net]
                originalVector = netHistory.originalVector()
                bond = netHistory.bondIndices()
                segment = np.logical_not(bond)
		real = np.logical_not(netHistory.stubIndices())
		negativeBonds = netHistory.negativeBondsIndices()
		legalBonds = np.logical_not(negativeBonds)
                queue = [(X, originalVector) for X in self.roots]
                totalLower = 0
		totalUpper = 0
		totalError = 0

                localEvents = dict((self.getTopEvent(netHistory, event), event) for event in netHistory.events)

                # I wish I could do this in proper recursion FP style, but Python is not very good with deep recursions
                while len(queue) > 0:
                        event, previousVector = queue.pop(0)

                        if (not debug.DEBUG) and event.ratio < debug.RATIO_CUTOFF:
                                continue

                        if event in localEvents:
                                localEvent = localEvents[event]
                                eventVector = netHistory.eventVector(localEvent) 
                                newVector = previousVector + eventVector

                                previousBonds = bond & (previousVector < 0)
                                deletedBonds = bond & (newVector > 0) & (previousVector >= 0)
				pseudoPreviousBonds = previousBonds | deletedBonds
                                createdBonds = bond & (newVector < 0) & (previousVector >= 0)
                                lower = self.computeLower(localEvent.cycle, netHistory, pseudoPreviousBonds, createdBonds)

				duplication = np.minimum(-previousVector, -eventVector)
				denovo = np.where(previousVector >= 0, -eventVector, 0) 
				supraduplications = np.sum(np.where(bond, np.maximum(2 * previousVector - newVector, 0), 0))
				imbalances = self.computeImbalances(localEvent.cycle, netHistory, duplication, denovo)
				upper = supraduplications + imbalances
				# If all the edges of the event cycle are bonds, then the flow is bond cycle 
				if np.sum((eventVector != 0) & bond) == len(localEvent.cycle) and np.sum(createdBonds) == len(localEvent.cycle) / 2:
					upper -= 1

				assert lower <= upper, "%i > %i = %i - 1\n%i bonds\n%i created\n%i ancestral\n%s" % (lower, upper, upper + 1, np.sum((eventVector != 0) & bond), np.sum(createdBonds), np.sum((eventVector != 0) & previousBonds), localEvent)
				self.eventCosts[event][0] += upper
				self.eventCosts[event][1] += lower
				totalLower += lower
				totalUpper += upper

                                segmentErrors = (eventVector != 0) & segment & ((previousVector == 0) | ((previousVector >= 0) & (newVector < 0))) & real
                                effectiveBondDestructions = (eventVector != 0) & bond & ((previousVector <= 0) & (newVector > 0)) 
				bondDestructions = effectiveBondDestructions & legalBonds
				totalError += 2 * np.sum(segmentErrors | bondDestructions)
                        else:
                                newVector = previousVector
                        queue.extend((X, newVector) for X in self.children[event])
                return totalUpper, totalLower, totalError
	
	def computeCost(self):
		self.eventCosts = dict((X, [0,0]) for X in self.parent)
		self.upper, self.lower, self.error = sumTrios(map(self.rearrangementCost_Net, self.netHistories))

        def errorCost(self):
                if self.error is None:
			self.computeCost()
                return self.error

        def rearrangementCost(self):
                # Yeah I know, ugly side effect, but avoids re-computing stuff.
                # You only compute complexity on the final version of a history. 
                # If you need to modify it after that, create a copy then carry on.
                if self.error is None:
			self.computeCost()
                return self.upper + self.error, self.lower + self.error

	def validate(self):
                assert super(ConstrainedHistory, self).validate()
		assert len(self.schedulingErrors()) == 0
