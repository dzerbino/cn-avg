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

"""Handling a collection of breakends"""

import sys
import vcf
import cnavg.avg.graph as avg

def mergeOverlappingBreakends(list, breakend):
	if len(list) == 0 or list[-1] < breakend:
		list.append(breakend)
	elif all(X is not breakend.partner for X in list):
		list[-1].absorb(breakend)
	elif all(X is not list[-1].partner for X in list):
		reject = list.pop()
		breakend.absorb(reject)  
		list.append(breakend)
	else:
		list.remove(breakend.partner)
		list[-1].absorb(breakend)

	return list

#########################################################
## Graph Structure
#########################################################

class BreakendGraph(list):
	""" A collection of breakends organized into a graph """
	def validate(self):
		assert all(map(lambda X: X.validate(), self))
		assert all(X.partner is None or X.partner in self for X in self)
		assert all(mate is None or mate in self for X in self for mate in X.mates)

	def searchBreakend(self, region):
		""" Search for breakend which covers a given region """ 
		return filter(lambda X: X == region, self)

	def _smallerNodes(self, region):
		return filter(lambda X: X < region, self)

	def _insertionPoint(self, region):
		return len(self._smallerNodes(region))

	def addBreakend(self, breakend):
		""" Insert breakend """
		self.insert(self._insertionPoint(breakend), breakend)

	#########################################################
	## Consolidating missing labels
	#########################################################

	def _consolidateEmptyBreakend(self, breakend):
		breakend.consolidateEmptyBreakend(self)

		if breakend.partner is None:
			breakend.createPartner(self)

		if breakend.mates is None:
			breakend.createMates(self)

	def consolidate(self):
		""" Ensuring all breakends in the graph have a mate and a partner """
		print "\tConsolidating breakends"
		map(lambda X: self._consolidateEmptyBreakend(X), self)

	#########################################################
	## Merge Overlapping breakends 
	#########################################################

	def mergeOverlapping(self):
		new = BreakendGraph(reduce(mergeOverlappingBreakends, sorted(self), []))
		for X in list(new):
			if X.partner not in new:
				assert len(X.mates) == 0, id(X)
				new.remove(X)
		return new

	#########################################################
	## Merging CNV info 
	#########################################################
	def assignCap(self, cnv, nodeindex):
		if nodeindex > 0:
			nodeindex -= 1

		while nodeindex < len(self) and cnv > self[nodeindex]:
			nodeindex += 1

		while nodeindex < len(self) and cnv == self[nodeindex]:
			assert cnv.val >= 0
			self[nodeindex].segment = cnv.val
			if self[nodeindex].start <= cnv.softStart and self[nodeindex].finish >= cnv.start:
				if cnv.startCap is None and self[nodeindex].orientation:
					cnv.startCap = self[nodeindex]
				else:
					pass
			elif self[nodeindex].start <= cnv.finish and self[nodeindex].finish >= cnv.softFinish:
				if not self[nodeindex].orientation:
					cnv.finishCap = self[nodeindex]
					break
			nodeindex += 1

		# Have to assign aritificial breakends
		if cnv.startCap is None:
			cnv.addStartCap(self)
		if cnv.finishCap is None:
			cnv.addFinishCap(self)

		return nodeindex

	def imposeSegmentFlow(self, cnv, nodeindex):
		while nodeindex < len(self) and self[nodeindex] != cnv.startCap:
			nodeindex += 1

		while nodeindex < len(self) and self[nodeindex] != cnv.finishCap:
			assert cnv.val >= 0
			self[nodeindex].segment = cnv.val
			nodeindex += 1
		# Last pass for finish cap:
		assert cnv.val >= 0
		cnv.finishCap.segment = cnv.val
		nodeindex += 1

		return nodeindex

	def assignCaps(self, cnvs):
		reduce(lambda L, E: self.assignCap(E, L), cnvs, 0)

	def imposeSegmentFlows(self, cnvs):
		reduce(lambda L, E: self.imposeSegmentFlow(E, L), cnvs, 0)

	def incorporateCNVs(self, cnvs):
		print "Mergeing CNVs to breakends"
		self.assignCaps(cnvs)
		self.imposeSegmentFlows(cnvs)

	#########################################################
	## Convert to graph
	#########################################################

	def _closeChromosome(self, chr, graph, chromosomeLength):
		chrBreakends = sorted(filter(lambda x: x.chr == chr and x.orientation == False, self), key = lambda X: X.node)
		if len(chrBreakends) == 0:
			return

		telomere1 = graph.createNode(chr, -1, True, name= chr + ".5prime")
		telomere2 = graph.createNode(chr, chromosomeLength + 1, False, name = chr + ".3prime")

		graph.telomeres.add(telomere1)
		graph.telomeres.add(telomere2)
		graph.createAdjacency(telomere1, telomere2)
		graph.addLiftedEdge(telomere1, telomere2, -1)

		assert all(sum(x.segment) >= 0 for x in self)
		assert all(sum(x.segment) >= 0 for x in chrBreakends)

		A = chrBreakends.pop(0)
		graph.createSegment(telomere1, A.node, A.segment)
		B = A.partner

		while len(chrBreakends) > 0:
			A = chrBreakends.pop(0)
			assert sum(A.segment) >= 0
			graph.createSegment(A.node, B.node, A.segment)
			B = A.partner

		graph.createSegment(telomere2, B.node, B.segment)

	def avg(self):
		print "Converting to AVG" 

		graph = avg.Graph()

		map(lambda X: X.attachNode(graph), self)
		map(lambda X: X.connectNode(graph), self)
		map(lambda X: self._closeChromosome(X, graph, self.lengths[X]), self.lengths)

		return graph

	def __str__(self):
		return "\n".join(map(str, self))

########################################################
## Test function
########################################################
def main():
	breakends = vcf.parse('../data/test.vcf', '../data/toto.lengths')
	print breakends
	print breakends.avg()
	breakends.validate()

if __name__ == "__main__":
	main()
