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

import sys
import cPickle as pickle
import random

import cnavg.cactus.oriented
import cnavg.cactusSampling.sampling
import cnavg.cactus.graph
import cnavg.cycleSampling.cycleCover as cycleCover
import cnavg.history.flattened as flattened
import cnavg.cycleSampling.sampleGraphCycles

def filterHistory(history):
	#return history.removeLowRatioEvents(0.1).simplifyStubs().selectForRegion("7",55054218,55242525)
	H = history.removeLowRatioEvents(0.1).simplifyStubs()
	events = sum((H.netHistories[N].events for N in H.netHistories), [])

	nodeFlowsA = filter(lambda X: X.inRegion("12",67488237,67525479), H.cactus.values())
	nodesA = [X.node for X in nodeFlowsA]
	eventsA = filter(lambda X: not X.doesNotContainNodes(nodesA), events)

	nodeFlowsB = filter(lambda X: X.inRegion("12",56427776,56432497), H.cactus.values())
	nodesB = [X.node for X in nodeFlowsB]
	eventsB = filter(lambda X: not X.doesNotContainNodes(nodesB), eventsA)

	return len(eventsB)

def filterHistories(histories):
	return map(filterHistory, histories)

def main():
	for file in sys.argv[1:]:
		try:
			H = pickle.load(open(file))
		except:
			continue
			
		# Preparing file for progressive write
		print H.printMajority()
		return
		pickle_file = open('complex.pickle', "wb")
		pickle.dump(H, pickle_file)
		braney_file = open("complex.braney", "w")
		braney_file.write("%s\n" % H.braneyText(0))
		SH = cnavg.cycleSampling.sampleGraphCycles.sample(H, 50000, file=pickle_file, braney=braney_file)
		return
		file = open("original", "w")
		file.write(H.braneyText(0, 999))
		file.close()

		file2 = open("flattened", "w")
		F = flattened.flattenGraph(H)
		file2.write(F.braneyText(0, 999))
		file2.close()

		file4 = open("simplified", "w")
		S = F.simplifyStubsAndTrivials()
		file4.write(S.braneyText(0, 999))
		file4.close()

		file3 = open("filtered", "w")
		FLT = S.removeLowRatioEvents(0.1)
		file3.write(FLT.braneyText(0, 999))
		file3.close()

if __name__ == "__main__":
	main()
