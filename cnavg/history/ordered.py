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

"""History with arbitrary linear ordering, for display purposes"""

from cnavg.history.scheduled import ScheduledHistory
from cnavg.basics.partialOrderSet import PartialOrderSet
import cnavg.history.flattened
import debug


class OrderedHistory(ScheduledHistory):
	""" History with arbitrary linear ordering compatible with tree structure """

	def __init__(self, cactusHistory):
		super(OrderedHistory, self).__init__(cactusHistory.cactus)
		self.copy(cactusHistory)
		self.ordering = PartialOrderSet(X for X in self.parent)
		for event in self.parent:
			if self.parent[event] is not None and (debug.DEBUG or event.ratio > debug.RATIO_CUTOFF): 
				assert self.ordering.addConstraint(self.parent[event], event)
	
	def braneyText(self, ID, cost=None):
		if cost is None:
			cost = self.rearrangementCost()
		return "\n".join(filter(lambda X: len(X) > 0, (X[1].braneyText(ID, X[0], self.ordering, cost) for X in enumerate(self.netHistories.values()))))

def prettify(H, i=0):
	""" Transform CactusHistory into BraneyText """
	c = H.rearrangementCost()
	FH = cnavg.history.flattened.flattenGraph(H)
	S = FH.simplifyStubsAndTrivials()
	F = S.removeLowRatioEvents(debug.RATIO_CUTOFF)
	O = OrderedHistory(F)
	return O.braneyText(i,c)

def main():
	if len(sys.argv) > 1:
		HH = pickle.load(open(sys.argv[1]))
		for i in range(len(HH)):
			print prettify(HH[i], i)
	else:
		G = avg.randomEulerianGraph(10)
		C = cactus.Cactus(G)
		N = normalizedCactus.NormalizedCactus(C)
		O = oriented.OrientedCactus(N)
		H = cnavg.cycleSampling.cycleCover.initialHistory(O)
		H.validate()
		print prettify(H, 0)

if __name__ == '__main__':
	import sys
	import cnavg.avg.graph as avg
	import cnavg.cactus.graph as cactus
	import cnavg.cactusSampling.sampling as normalizedCactus
	import cnavg.cactus.oriented as oriented
	import cnavg.cycleSampling.cycleCover
	main()
