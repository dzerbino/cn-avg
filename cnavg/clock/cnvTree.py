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

import rpy2
from rpy2.robjects import IntVector
from rpy2.robjects import FloatVector
from numpy import array
from numpy import ones

PARAMETER = float(50)
dpois = rpy2.robjects.r.dpois

###################################################
## Copy number change
###################################################
def computeCopyNumberChange(event, overlaps):
	if event in overlaps:
		if event.cycle[overlaps[event][0]].value > 0:
			return -len(overlaps[event]) * math.ceil(abs(event.cycle[0].value))
		else:
			return len(overlaps[event]) * math.ceil(abs(event.cycle[0].value))
	else:
		return 0

###################################################
## Enumerating CNV SliceTrees
###################################################

class CNVSliceTree(list):
	def __init__(self, iter, count, event, prior=None):
		super(CNVSliceTree, self).__init__(iter)
		self.count = count
		self.event = event
		self.prior = prior

	# Because list not hashable
	def __hash__(self):
		return id(self)

	def __or__(self, other):
		return CNVSliceTree(super(CNVSliceTree, self).__add__(other), self.count * other.count, self.event)

	def nodes(self):
		# Careful, nodes can be carried over across slices, so set() removes doubles
		return set(sum(map(list, self), []))

	def originators(self):
		return filter(lambda X: X.parent is None, self.nodes())

def organizeByEvents(eventSliceTrees, sliceTrees):
	if sliceTrees[0].event in eventSliceTrees:
		eventSliceTrees[sliceTrees[0].event].extend(sliceTrees)
	else:
		eventSliceTrees[sliceTrees[0].event] = sliceTrees
	return eventSliceTrees

def crossProduct(list, newList):
	if len(list) == 0:
		return newList
	else:
		return [X | Y for X in list for Y in newList]

class CNVSlice(set):
	def __init__(self, iter, parent, count, event):
		super(CNVSlice, self).__init__(iter)
		self.count = count
		self.weight = float(sum(X.copynumber for X in self))
		if parent is not None:
			parent.children.append(self)
		self.children = []
		self.event = event

	# Because set not hashable
	def __hash__(self):
		return id(self)

	def __eq__(self, other):
		return id(self) == id(other)

	def extractSliceTrees(self):
		if len(self.children) == 0:
			return [CNVSliceTree([self], self.count, self.event)]
		elif len(self.children) == 1:
			return [CNVSliceTree(list(X) + [self], X.count * self.count, self.event) for X in self.children[0].extractSliceTrees()]
		else:
			# Recursion 
			childSliceTrees = map(lambda X: X.extractSliceTrees(), self.children)
			# Organize by creator event to avoid self-products
			childEventSliceTrees = reduce(organizeByEvents, childSliceTrees, dict())
			# Cardinal product of child SliceTrees
			sliceTrees = reduce(crossProduct, childEventSliceTrees.values(), [])
			# Iterative step
			return [CNVSliceTree(list(X) + [self], X.count * self.count, self.event) for X in sliceTrees]

	def cnvNodeWeight(self, cnvNode):
		if self.weight == 0:
			return 1
		else:
			return cnvNode.copynumber / self.weight

	def add(self, element):
		super(CNVSlice, self).add(element)
		self.weight += element.copynumber

	def remove(self, element):
		super(CNVSlice, self).remove(element)
		self.weight -= element.copynumber

class CNVNode(object):
	def __init__(self, copynumber, creator, parent):
		self.copynumber = copynumber
		self.creator = creator
		self.parent = parent

	def shadow(self):
		return CNVNode(self.copynumber - 1, None, self)

	def __str__(self):
		return "\t".join(map(str, [self.copynumber, self.creator, self.parent.event]))

def newHistorySlice_Node(previousSlice, event, change, cnvNode):
	if cnvNode.copynumber == 0:
		return []
	else:
		res = CNVSlice(previousSlice, previousSlice, cnvNode.copynumber, event)
		res.remove(cnvNode)

		if cnvNode.copynumber - 1 > 0:
			res.add(cnvNode.shadow())

		if change > 0:
			res.add(CNVNode(change + 1, event, cnvNode))
		else:
			res.add(CNVNode(0, event, cnvNode))

		return [res]

def newHistorySlice(previousSlice, event, change):
	return sum((newHistorySlice_Node(previousSlice, event, change, cnvNode) for cnvNode in previousSlice), [])

def virginBirth(event, change):
	if change < 0:
		return CNVNode(0, event, None)
	elif change > 0:
		return CNVNode(change + 1, event, None)
	else:
		assert False

def newHistorySlices(previousSlices, event, change):
	if change == 0:
		return [CNVSlice(X, X, 1, event) for X in previousSlices]
	else:
		res = sum((newHistorySlice(slice, event, change) for slice in previousSlices), [])
		if len(res) > 0:
			return res
		else:
			# If all the nodes have copy number 0, then the result is an empty list the history is messed up but we 
			# try to pick up the pieces as we can	
			return [CNVSlice(list(slice) + [virginBirth(event, change)], slice, 1, event) for slice in previousSlices]

def getChildren(children, node):
	if node.parent is not None:
		children[node.parent].append(node)
	return children

def updateLikelihoods(likelihoods, sliceTree, data):
	if data[0] not in likelihoods:
		likelihoods[data[0]] = dict()
	if sliceTree not in likelihoods[data[0]]:
		likelihoods[data[0]][sliceTree] = dict()
	likelihoods[data[0]][sliceTree][event] = data[1]
	return likelihoods

###################################################
## CNVTree
###################################################

class CNVTree(object):
	###################################################
	## Overall construction
	###################################################

	def __init__(self, cactusHistory, overlaps, originalCopyNumber):
		# Recoding atomic info
		events = cactusHistory.parent.keys()
		self.cactusHistory = cactusHistory
		if overlaps is None:
			self.copyNumberChanges = None
		else:
			self.copyNumberChanges = dict((event, computeCopyNumberChange(event, overlaps)) for event in events)
			self.impactChanges = dict((event, event.ratio * self.copyNumberChanges[event]) for event in events)

			# Reconstructing possible molecular histories
			self.rootNode = CNVNode(originalCopyNumber, None, None)
			self.rootSlice = CNVSlice([self.rootNode], None, 1, None)
			map(lambda X: self.computeHistorySlices(X,[self.rootSlice]), cactusHistory.roots)
			preSliceTrees = self.rootSlice.extractSliceTrees()
			totalSliceTreeCount = float(sum(X.count for X in preSliceTrees))
			self.sliceTrees = [CNVSliceTree(list(X), X.count, X.event, X.count/totalSliceTreeCount) for X in preSliceTrees]

			# Computing expected impacts
			self.children = dict((X, self.getSliceTreeChildren(X)) for X in self.sliceTrees)
			self.cumulativeImpactChanges = dict((sliceTree, self.cumulativeImpactChanges_SliceTree(sliceTree)) for sliceTree in self.sliceTrees)
			self.totalImpacts = dict((sliceTree, self.totalImpact(sliceTree)) for sliceTree in self.sliceTrees)

	###################################################
	## Building the damn thing
	###################################################

	def computeHistorySlices(self, event, previousSlices):
		eventSlices = newHistorySlices(previousSlices, event, self.copyNumberChanges[event])
		map(lambda X: self.computeHistorySlices(X, eventSlices), self.cactusHistory.children[event])

	def getSliceTreeChildren(self, sliceTree):
		return reduce(getChildren, sliceTree.nodes(), dict((node, []) for node in sliceTree.nodes()))

	def cumulativeImpactChanges_CNVNode(self, cumulativeImpactChanges, cnvNode, sliceTree):
		cumulativeImpactChanges = reduce(lambda X, Y: self.cumulativeImpactChanges_CNVNode(X,Y,sliceTree), self.children[sliceTree][cnvNode], cumulativeImpactChanges)
		cumulativeImpactChanges[cnvNode] = sum(cumulativeImpactChanges[X] for X in self.children[sliceTree][cnvNode]) + sum(C.creator.ratio * self.copyNumberChanges[C.creator] for C in self.children[sliceTree][cnvNode] if C.creator is not None)
		return cumulativeImpactChanges

	def cumulativeImpactChanges_SliceTree(self, sliceTree):
		return reduce(lambda X, Y: self.cumulativeImpactChanges_CNVNode(X,Y,sliceTree), sliceTree.originators(), dict())

	def totalImpact_Node(self, sliceTree, cnvNode):
		if cnvNode.creator is not None:
			return max(cnvNode.creator.ratio + self.cumulativeImpactChanges[sliceTree][cnvNode], 0)
		else:
			return max(1 + self.cumulativeImpactChanges[sliceTree][self.rootNode], 0)

	def totalImpact(self, sliceTree):
		return max(sum(self.totalImpact_Node(sliceTree, cnvNode) for cnvNode in sliceTree.originators()), 0.00001)

	###################################################
	## Likelihoods
	###################################################

	def expectedImpact(self, slice, sliceTree, cnvNode):
		if slice.event is not None:
			return max(slice.event.ratio + self.cumulativeImpactChanges[sliceTree][cnvNode], 0)
		else:
			return max(1 + self.cumulativeImpactChanges[sliceTree][self.rootNode], 0)

	def expectedRatio(self, slice, sliceTree, cnvNode):
		return self.expectedImpact(slice, sliceTree, cnvNode) / float(self.totalImpacts[sliceTree])

	def expectedCoverages(self, cnvNode, slice, sliceTree, totals):
		return FloatVector(self.expectedRatio(slice, sliceTree, cnvNode) * totals)

	# Poisson likelihoods (R function then numpy array)
	def conditionedLikelihoods(self, cnvNode, slice, sliceTree, counts, totals):
		return array(dpois(counts, self.expectedCoverages(cnvNode, slice, sliceTree, totals)))
	
	def segregatingLikelihood(self, vals):
		# TODO The likelihood function of a segregating SNP could be refined
		if vals[1] == 0:
			return 0
		else:
			return PARAMETER * math.exp(-(vals[0]/float(vals[1])) * PARAMETER)

	def segregatingLikelihoodArray(self, counts, totals):
		return map(self.segregatingLikelihood, zip(counts, totals))

	# Weighted sum
	def likelihoodArray(self, slice, sliceTree, counts, totals):
		if slice is self.rootSlice:
			return self.segregatingLikelihoodArray(counts,totals)
		else:
			return sum(slice.cnvNodeWeight(cnvNode) * self.conditionedLikelihoods(cnvNode, slice, sliceTree, counts, totals) for cnvNode in slice)

	# Going through events
	# Returns dictionary: snv -> likelihood
	def sliceLikelihoods(self, slice, sliceTree, snvs, counts, totals):
		return dict(zip(snvs, self.likelihoodArray(slice, sliceTree, counts, totals)))

	# Turns event -> snv -> likelihood
	# Into array -> likelihood
	def snvLikelihoods(self, snv, timings, eventLikelihoods):
		return array([eventLikelihoods[timing][snv] for timing in timings])

	# Going through slices
	# Returns dictionary: snv -> array -> likelihood
	def sliceTreeLikelihoods(self, sliceTree, snvs, timings, counts, totals):
		# Computing sliceTree stats
		eventLikelihoods = dict((slice.event, self.sliceLikelihoods(slice, sliceTree, snvs, counts, totals)) for slice in sliceTree)
		return dict((snv, self.snvLikelihoods(snv, timings, eventLikelihoods)) for snv in snvs)

	def defaultLikelihoods_SNV(self, timing, counts, totals):
		if timing is None:
			return array(dpois(counts, FloatVector(totals)))
		else:
			return array(dpois(counts, FloatVector(timing.ratio * totals)))

	# Returns dictionary: snv -> array -> likelihood
	def defaultLikelihoods_Event(self, snvs, timing, counts, totals):
		return dict(zip(snvs, self.defaultLikelihoods_SNV(timing, counts, totals)))

	# Returns dictionary: snv -> array -> likelihood
	def defaultLikelihoods(self, snvs, timings, counts, totals):
		eventLikelihoods = dict((timing, self.defaultLikelihoods_Event(snvs, timing, counts, totals)) for timing in timings)
		return dict((snv, self.snvLikelihoods(snv, timings, eventLikelihoods)) for snv in snvs)

	# Returns dictionary: sliceTree -> snv -> array -> likelihood
	def likelihoods(self, counts_list, totals_list, snvs, events, phasedRatio):
		timings = events + [None]
		# Creating objects to avoid multiple redundant creations
		counts = IntVector(counts_list)
		totals = array(totals_list) * phasedRatio
		if self.copyNumberChanges is None:
			return dict([(CNVSliceTree([], 1, None, 1), self.defaultLikelihoods(snvs, timings, counts, totals))])
		else:
			return dict((sliceTree, self.sliceTreeLikelihoods(sliceTree, snvs, timings, counts, totals)) for sliceTree in self.sliceTrees)

	###################################################
	## Convenience
	###################################################

	def eventString(self, event):
		string = '%i [label="%i,%f,%f,%f,%f"]' % (id(event), self.copyNumberChanges[event], event.ratio, self.impactChanges[event], self.cumulativeImpactChanges[event])
		string += "\n" + "\n".join("%i -> %i" % (id(event), id(X)) for X in self.cactusHistory.children[event]) 
		return string

	def __str__(self):
		string = "digraph G {\n"
		string += "\n".join(map(self.eventString, self.cactusHistory.parent.keys())) + "\n"
		string += "}\n"
		return string

##############################################
## Unit testing
##############################################

def main():
	sys.setrecursionlimit(5000)
	if len(sys.argv) == 1:
		G = avg.randomNearEulerianGraph(10)
		C = cactus.Cactus(G)
		N = normalizedCactus.NormalizedCactus(C)
		O = oriented.OrientedCactus(N)
		H = cycleCover.initialHistory(O)
	else:
		file = open(sys.argv[1])
		H = None
		while True:
			try:
				tmp = pickle.load(file)
				if tmp is not None:
					H = tmp
			except:
				break
		file.close()
		
	FH = flattened.flattenGraph(H)
	S = FH.simplifyStubsAndTrivials()
	F = S.removeLowRatioEvents(0.1)

	if len(sys.argv) == 1:
		print F

	overlapTable = dict()
	for e in F.parent:
		overlapTable = overlap.addEventToOverlapTable(overlapTable, e)
	counts = []
	totalCounts = []
	for edgeIndex in overlapTable.keys():
		if edgeIndex[2] > -1:
			t = CNVTree(F, overlapTable[edgeIndex], 1)
			counts.append(len(t.sliceTrees))
			totalCounts.append(sum(X.count for X in t.sliceTrees))
			print [len(X.originators()) for X in t.sliceTrees]
			print [len(X.nodes()) for X in t.sliceTrees]
	print "\n".join(["\t".join(X) for X in zip(map(str, counts), map(str, totalCounts))])

if __name__ == "__main__":
	import sys 
	import cPickle as pickle
	import cnavg.avg.graph as avg
	import cnavg.cactus.graph as cactus
	import cnavg.cactus.oriented as oriented
	import cnavg.cycleSampling.cycleCover as cycleCover
	import cnavg.cactusSampling.sampling as normalizedCactus
	import cnavg.history.flattened as flattened
	import cnavg.history.overlap as overlap
	main()
