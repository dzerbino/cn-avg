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

"""Produces random evolutionary histories"""

import sys
import random
import cnavg.avg.graph
import cnavg.basics.partialOrderSet

BRANCHPROB = 0
MEAN_INDEL_LENGTH = 10
MEAN_TANDEMS = 1

""" Probability that a branch point occurs at a given node """

class HistoryBranch(object):
        """Branch in a history"""
        #########################################
        ## Basics
        #########################################
        def __init__(self):
                self.children = []

        def __str__(self):
                return "\n".join([self._label()] + map(str, self.children))

        #########################################
        ## Stats
        #########################################
        def _enumerate(self):
                return [self] + sum((X._enumerate() for X in self.children), [])

        def _cost(self):
                return self._operationCost() + sum(X._cost() for X in self.children)

        def _testPosition(self, position):
                if len(self.children):
                        return any(child._testPosition(position) for child in self.children)
                else:
                        return True

        #########################################
        ## GraphViz representation
        #########################################
        def _dot(self, weights):
                """ GraphViz output """
                return "\n".join(["digraph G {" , "node [shape=rectangle]", self._dot2(weights), "}"])

        def _dot2(self, weights):
                """ GraphViz output """
                return "\n".join([self._dotString(weights)] + [X._dot2(weights) for X in self.children])

        def _dotString(self, weights):
                """ GraphViz output """
                return "\n".join([self._dotLabel(weights)] + [self._dotEdge(X) for X in self.children])

        def _dotLabel(self, weights):
                label = str(self.genome)
                if len(self.children) == 0:
                        label += " (%f)" % weights[self]
                return '%i [label="%s"]' % (id(self), label)

        def _dotEdge(self, child):
                return '%i -> %i [label="%s"]' % (id(self), id(child), child._dotBlurb())

        #########################################
        ## Braney representation 
        #########################################

	def _braneyText(self, history, avg, cost, poset):
		for child in self.children:
			poset.addElement(child)
			poset.addConstraint(self, child)

		return "\n".join([self._braneyBlurb(history, avg, cost, poset)] + [X._braneyText(history, avg, cost, poset) for X in self.children])

#########################################
## Top Branch
#########################################

class InitialBranch(HistoryBranch):
        """Root branch of an evolutionary history"""

        def __init__(self, length):
                super(InitialBranch, self).__init__()
                self.genome = range(1,length+ 1)

        def _label(self):
                return "INIT\n" + str(self.genome)

        def _operationCost(self):
                return 0

	def _braneyBlurb(self, history, avg, cost, poset):
		return ""

#########################################
## Transformation branches
#########################################

class Operation(HistoryBranch):
        """Non-root branch in an evolutionary history"""

        def __init__(self, parent):
                super(Operation, self).__init__()
                self.parent = parent
                parent.children.append(self)
                self.genome = self._product(parent.genome)

        def _operationCost(self):
                if self._persists():
                        return 1
                else:
                        return 0

	def followsAncestralSequence(self, history, index):
		return (self.parent.genome[index % len(self.parent.genome)] % len(history.root.genome) == (self.parent.genome[(index - 1) % len(self.parent.genome)] + 1) % len(history.root.genome)
		       and self.parent.genome[index % len(self.parent.genome)] * self.parent.genome[(index - 1) % len(self.parent.genome)] > 0)

class Identity(Operation):
        """Nothing happens"""
        def __init__(self, parent):
                super(Identity, self).__init__(parent)

        def _product(self, genome):
                return genome

        def _label(self):
                return "ID"

        def _dotBlurb(self):
                return "ID"

	def _braneyBlurb(self, history, avg, cost, poset):
		return ""

        def _operationCost(self):
                return 0

        def _persists(self):
                return True

def signChar(value, revComp = False):
	if revComp:
		if value > 0:
			return '-'
		else:
			return '+'
	else:
		if value > 0:
			return '+'
		else:
			return '-'

def getNode(nodes, index, _3prime):
	if index > 0:
		if _3prime:
			return nodes[2 * index - 1]
		else:
			return nodes[2 * index - 2]
	if index < 0:
		if _3prime:
			return nodes[-2 * index - 2]
		else:
			return nodes[-2 * index - 1]

def orientString(node):
	if node.orientation:
		return '+'
	else:
		return '-'

class Inversion(Operation):
        """Inversion branch"""
        def __init__(self, parent, start, length):
                self.start = start
                self.length = length
                super(Inversion, self).__init__(parent)

        def _product(self, genome):
                if self.start + self.length < len(genome):
                        return genome[:self.start] + [-X for X in reversed(genome[self.start:self.start+self.length])] + genome[self.start + self.length:]
                else:
                        return genome[:self.start + self.length - len(genome)] + [-X for X in reversed(genome[self.start + self.length - len(genome):self.start])] + genome[self.start:]

        def _label(self):
                return "INV\t%i\t%i\n%s" % (self.start, self.length, str(self.genome))

        def _dotBlurb(self):
                str = "INV %i,%i" % (self.start, self.length)
                if self._persists():
                        return str
                else:
                        return "*" + str

        def _testPosition(self, position):
                position = position % len(self.genome)
                if len(self.children):
                        if position >= self.start and position < self.start + self.length:
                                newPosition = 2*self.start + self.length - position - 1
                        else:
                                newPosition = position
                        return any(child._testPosition(newPosition) for child in self.children)
                else:
                        return True

        def _persists(self):
                return (self._testPosition(self.start) and self._testPosition(self.start + self.length)) or (self._testPosition(self.start - 1) and self._testPosition(self.start + self.length - 1))

	def _braneyBlurb(self, history, avg, cost, poset):
		nodes = sorted(avg.nodes())
		rank = poset.depth[self]

		startA = getNode(nodes, self.parent.genome[(self.start - 1) % len(self.parent.genome)], True)
		startB = getNode(nodes, self.parent.genome[(self.start) % len(self.parent.genome)], False) 
		finishB = getNode(nodes, self.parent.genome[(self.start + self.length - 1) % len(self.parent.genome)], True)
		finishA = getNode(nodes, self.parent.genome[(self.start + self.length) % len(self.parent.genome)], False)
		prevalence = history.weights[self]

		return "\n".join([
		       "\t".join(map(str, ['A', startA.chr, startA.pos, orientString(startA), startB.chr, startB.pos, orientString(startB),-prevalence,prevalence, 0,0,rank,0,rank,cost,cost,1,1,id(self)])),
		       "\t".join(map(str, ['A', startB.chr, startB.pos, orientString(startB), finishA.chr, finishA.pos, orientString(finishA),prevalence,prevalence, 0,0,rank,1,rank,cost,cost,1,1,id(self)])),
		       "\t".join(map(str, ['A', finishA.chr, finishA.pos, orientString(finishA), finishB.chr, finishB.pos, orientString(finishB),-prevalence,prevalence, 0,0,rank,2,rank,cost,cost,1,1,id(self)])),
		       "\t".join(map(str, ['A', finishB.chr, finishB.pos, orientString(finishB), startA.chr, startA.pos, orientString(startA),prevalence,prevalence, 0,0,rank,3,rank,cost,cost,1,1,id(self)])),
		       ])

class Duplication(Operation):
        """Duplication branch"""
        def __init__(self, parent, start, length, count):
                self.start = start
                self.length = length
		self.count = count
                super(Duplication, self).__init__(parent)

        def _product(self, genome):
                if self.start + self.length < len(genome):
                        return genome[:self.start] + sum((genome[self.start:self.start+self.length] for X in range(self.count)), [])  + genome[self.start:]
                else:
                        return genome + sum((genome[:self.start+self.length-len(genome)] + genome[self.start:] for X in range(self.count)), [])

        def _label(self):
                return "DUP\t%i\t%i\t%i\n%s" % (self.start, self.length, self.count, str(self.genome))

        def _dotBlurb(self):
                str = "DUP %i,%i,%i" % (self.start, self.length, self.count)
                if self._persists():
                        return str
                else:
                        return "*" + str

        def _testPosition(self, position):
                position = position % len(self.genome)
                if len(self.children):
                        if position >= self.start and position < self.start + self.length:
                                return any(child._testPosition(position) for child in self.children) or any(child._testPosition(position + self.length) for child in self.children)
                        if position >= self.start + self.length:
                                return any(child._testPosition(position + self.length) for child in self.children)
                        else:
                                return any(child._testPosition(position) for child in self.children)
                else:
                        return True

        def _persists(self):
                return any(self._testPosition(X) for X in range(self.start, self.start + self.length))

	def _braneyBlurb(self, history, avg, cost, poset):
		nodes = sorted(avg.nodes())
		prevalence = history.weights[self] * self.count
		lines = []
		edgeIndex = 0
		rank = poset.depth[self]

		previous = getNode(nodes, self.parent.genome[(self.start + self.length - 1) % len(self.parent.genome)], True)
		i = 0
		while i < self.length:
			start = getNode(nodes, self.parent.genome[(self.start + i) % len(self.parent.genome)], False)
			shift = 1;
			while i + shift < self.length and self.followsAncestralSequence(history, self.start + i + shift):
				shift += 1
			finish = getNode(nodes, self.parent.genome[(self.start + i + shift - 1) % len(self.parent.genome)], True)
			lines.append("\t".join(map(str, ['A', previous.chr, previous.pos, orientString(previous), start.chr, start.pos, orientString(start),prevalence,prevalence, 0,0,rank,edgeIndex,rank,cost,cost,1,1,id(self)])))
			lines.append("\t".join(map(str, [start.chr, start.pos, finish.pos, -prevalence,prevalence, 0,0,rank,edgeIndex+1,rank,cost,cost,1,1,id(self)])))
			previous = finish
			edgeIndex += 2
			i += shift

		return "\n".join(lines)

class Deletion(Operation):
        """Deletion branch"""
        def __init__(self, parent, start, length):
                self.start = start
                self.length = length
                super(Deletion, self).__init__(parent)

        def _product(self, genome):
                if self.start + self.length < len(genome):
                        return genome[:self.start] + genome[self.start+self.length:]
                else:
                        return genome[(self.start + self.length) - len(genome):self.start]

        def _label(self):
                return "DEL\t%i\t%i\n%s" % (self.start, self.length, str(self.genome))

        def _dotBlurb(self):
                str = "DEL %i,%i" % (self.start, self.length)
                if self._persists():
                        return str
                else:
                        return "*" + str

        def _testPosition(self, position):
                if len(self.genome) == 0:
                        return False
                position = position % len(self.genome)
                if position < self.start or position >= self.start + self.length:
                        if len(self.children):
                                if position >= self.start + self.length:
                                        return any(child._testPosition(position - self.length) for child in self.children)
                                else:
                                        return any(child._testPosition(position) for child in self.children)
                        else:
                                return True
                else:
                        return False

        def _persists(self):
                return self._testPosition(self.start - 1) or self._testPosition(self.start + self.length)

	def _braneyBlurb(self, history, avg, cost, poset):
		nodes = sorted(avg.nodes())
		prevalence = history.weights[self]
		edgeIndex = 0
		lines = []
		rank = poset.depth[self]

		leftFlank = getNode(nodes, self.parent.genome[(self.start - 1) % len(self.parent.genome)], True)
		previous = leftFlank
		i = 0
		while i < self.length:
			start = getNode(nodes, self.parent.genome[(self.start + i) % len(self.parent.genome)], False)
			shift = 1;
			while i + shift < self.length and self.followsAncestralSequence(history, self.start + i + shift):
				shift += 1

			finish = getNode(nodes, self.parent.genome[(self.start + i + shift - 1) % len(self.parent.genome)], True)
			lines.append("\t".join(map(str, ['A', previous.chr, previous.pos, orientString(previous), start.chr, start.pos, orientString(start),-prevalence,prevalence, 0,0,rank,edgeIndex,rank,cost,cost,1,1,id(self)])))
			lines.append("\t".join(map(str, [start.chr, start.pos, finish.pos, prevalence,prevalence, 0,0,rank,edgeIndex + 1,rank,cost,cost,1,1,id(self)])))
			previous = finish
			edgeIndex += 2
			i += shift

		rightFlank = getNode(nodes, self.parent.genome[(self.start + self.length) % len(self.parent.genome)], False)
		lines.append("\t".join(map(str, ['A', previous.chr, previous.pos, orientString(previous), rightFlank.chr, rightFlank.pos, orientString(rightFlank),-prevalence,prevalence, 0,0,rank,edgeIndex,rank,cost,cost,1,1,id(self)])))
		lines.append("\t".join(map(str, ['A', rightFlank.chr, rightFlank.pos, orientString(rightFlank), leftFlank.chr, leftFlank.pos, orientString(leftFlank),prevalence,prevalence, 0,0,rank,edgeIndex+1,rank,cost,cost,1,1,id(self)])))

		return "\n".join(lines)

#########################################
## Evolutionary History
#########################################
class History(object):
        """Evolutionary history"""

        def __init__(self, root):
                self.root = root

        def cost(self):
                return self.root._cost()

        def enumerate(self):
                return self.root._enumerate()

        def __str__(self):
                return str(self.root)

        def dot(self):
                """ GraphViz output """
                return self.root._dot(self.weights)

	def braneyText(self, avg):
		return self.root._braneyText(self, avg, self.cost(), cnavg.basics.partialOrderSet.PartialOrderSet([self.root]))

#########################################
## Weighted Evolutionary History
#########################################

def _removeInitialPath(avg):
        nodes = sorted(avg.nodes())
        previousNode = nodes[-1]
                
        for index in range(len(nodes)/2):
                nextNode = nodes[2 * index]
                avg.setLiftedEdge(previousNode, nextNode, -1)
                avg.changeSegment(nextNode, 0, -1)
                previousNode = avg[nextNode].twin
        return avg

class WeightedHistory(History):
        """Evolutionary history with weighted branches"""
        def __init__(self, root, weights):
                super(WeightedHistory, self).__init__(root)
                self.weights = weights

        def _normalizeBranches(self, weights):
                total = float(sum(weights.values()))
                return dict((X, weights[X]/total) for X in weights)

        def _weightBranches(self):
                return self._normalizeBranches(self._weightBranchesAtRandom())

	def _weightAllBranches(self, branch):
		if branch in self.weights:
			return self.weights[branch]
		else:
			return sum(map(self._weightAllBranches, branch.children))

        def _avg_Branch(self, avg, branch):
                if len(branch.children) == 0 and len(branch.genome) > 0:
                        weight = self.weights[branch]
                        nodes = sorted(avg.nodes())
                        if branch.genome[-1] > 0:
                                previousNode = nodes[2 * branch.genome[-1] - 1]
                        else:
                                previousNode = nodes[2 * -branch.genome[-1] - 2]
                                
                        for index in range(len(branch.genome)):
                                if branch.genome[index] > 0:
                                        nextNode = nodes[2 * branch.genome[index] - 2]
                                else:
                                        nextNode = nodes[2 * -branch.genome[index] - 1]
                                
                                avg.addLiftedEdge(previousNode, nextNode, weight)
                                avg.changeSegment(nextNode, 0, weight)
                                previousNode = avg[nextNode].twin
                return avg

        def avg(self):
                """Produce resultant sequence graph with flow values"""
                avg = _removeInitialPath(cnavg.avg.graph.linearGraph(len(self.root.genome) * 2))
                return reduce(self._avg_Branch, self.enumerate(), avg)

#########################################
## Random Evolutionary History
#########################################
def _addChildBranch(branch, inv_prob = .7, dup_prob = .1, del_prob = .1, id_prob = .1):
	if len(branch.genome) > 1:
		choice = random.random()
	else:
		choice = 1
        start = random.randrange(len(branch.genome))

        if choice < inv_prob:
                length = random.randrange(1, len(branch.genome))
                Inversion(branch, start, length)
        elif choice < inv_prob + dup_prob:
                # Separate length prob for different operations (e.g. long distance duplications followed by deletions make things moot)
                length = int(random.expovariate(1/float(MEAN_INDEL_LENGTH)))
                if length >= len(branch.genome):
                        length = len(branch.genome) - 1
                if length == 0:
                        length = 1
		count = max(1, int(random.expovariate(1/float(MEAN_TANDEMS))))
		Duplication(branch, start, length, count)
        elif choice < inv_prob + dup_prob + del_prob:
                length = int(random.expovariate(1/float(MEAN_INDEL_LENGTH)))
                if length >= len(branch.genome):
                        length = len(branch.genome) - 1
                if length == 0:
                        length = 1
                Deletion(branch, start, length)
        else:
                Identity(branch)

def _birthDeathModel():
        choice = random.random()
        if choice < BRANCHPROB:
                return 2
        else:
                return 1

def _extendHistory_Branch(branch):
        if len(branch.genome) == 0:
                return
        for i in range(_birthDeathModel()):
                _addChildBranch(branch)

def _extendHistory(branch, counter):
        newBranches = _extendHistory_Branch(branch)
        if counter > 1:
                map(lambda X: _extendHistory(X, counter - 1), branch.children)

class RandomHistory(WeightedHistory):
        def __init__(self, length, maxDepth):
                root = InitialBranch(length)
                _extendHistory(root, maxDepth)
                super(RandomHistory, self).__init__(root, None)

#########################################
## Random Weighted Evolutionary History
#########################################

class RandomIntegerHistory(RandomHistory):
        def _weightBranchesAtRandom(self):
                return dict((X,1) for X in self.enumerate() if len(X.children) == 0)

        def __init__(self, length, maxDepth, indel_length=MEAN_INDEL_LENGTH):
		global MEAN_INDEL_LENGTH
		MEAN_INDEL_LENGTH = indel_length
		global BRANCHPROB
		BRANCHPROB = 0
                super(RandomIntegerHistory, self).__init__(length, maxDepth)
                self.weights = self._weightBranches()
		self.weights = dict((X, self._weightAllBranches(X)) for X in self.enumerate())

#########################################
## Random Cancer History
#########################################

class RandomCancerHistory(RandomHistory):
	def __init__(self, length, maxDepth, indel_length=MEAN_INDEL_LENGTH):
		global MEAN_INDEL_LENGTH
		MEAN_INDEL_LENGTH = indel_length
		global BRANCHPROB
		BRANCHPROB = 0
                super(RandomCancerHistory, self).__init__(length, maxDepth)
		germline = Identity(self.root)
		somatic = filter(lambda X:len(X.children) == 0, self.enumerate())[0]
                _extendHistory(somatic, maxDepth)
		somatic2 = filter(lambda X:len(X.children) == 0, self.enumerate())[0]
                _extendHistory(somatic2, maxDepth)
		somatic3 = filter(lambda X:len(X.children) == 0, self.enumerate())[0]

		self.weights = dict([(somatic3, .60), (Identity(somatic2), .10), (Identity(somatic), .10), (germline, .20)])
		self.weights = dict((X, self._weightAllBranches(X)) for X in self.enumerate())

#########################################
## Unit test
#########################################
def main():
        debug.DEBUG = True
        debug.PLOIDY = 1
        #history = RandomCancerHistory(3, 4, 1)
        history = RandomIntegerHistory(10, 1)
        G = history.avg()
        C = cactus.Cactus(G)
        O = oriented.OrientedCactus(C)
        H = cycleCover.initialHistory(O)
        print C
        print history
        print history.dot()
        print history.cost()
	print history.braneyText(G)
        print H
        print H.rearrangementCost()
	return
        for NH in H.netHistories.values():
                print NH.dot()

if __name__ == '__main__':
	import cnavg.cactus.graph as cactus
	import cnavg.cactus.oriented as oriented 
	import cnavg.historySampling.cycleCover as cycleCover
	import cnavg.history.euclidian as euclidian
	import cnavg.history.constrained as constrained
	import cnavg.history.debug as debug
        main()
