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

"""Sampling histories at the local level"""

import sys
import random
import math
import copy

import cnavg
import cnavg.avg.graph as avg
import cnavg.cactus.graph as cactus
from cnavg.flows.flows import Event
from cnavg.flows.cycle import Cycle
from cnavg.flows.edge import Edge
from cnavg.history.history import History
from cnavg.history.history import CactusHistory
from cnavg.history.overlap import OverlapHistory
from cnavg.history.overlap import Overlap
from cnavg.history.euclidian import EuclidianHistory
import cnavg.cactus.oriented as orientedCactus
import cnavg.cactusSampling.sampling as normalized
import cycleCover
import cnavg.flows.unitFlows as unitFlows

ROUNDING_ERROR = 1e-2

#########################################
## Wrapper for simplification
#########################################

def addEvents(cactusHistory, history, events):
	for event in unitFlows.simplifyEventCycles(events):
		cactusHistory.absorbEvent(history, event)

#########################################
## Getting boundaries
#########################################

def getStartOfRepeat(cycle1, start1, pos1, cycle2, start2, pos2):
	next1 = (pos1 - 1) % len(cycle1)
	next2 = (pos2 - 1) % len(cycle2)
	if next1 == start1 or next2 == start2:
		return pos1, pos2
	if cycle1[next1].start != cycle2[next2].start:
		return pos1, pos2
	elif cycle1[next1].index != cycle2[next2].index:
		return pos1, pos2
	else:
		return getStartOfRepeat(cycle1, start1, next1, cycle2, start2, next2)

def getEndOfRepeat(cycle1, start1, pos1, cycle2, start2, pos2):
	next1 = (pos1 + 1) % len(cycle1)
	next2 = (pos2 + 1) % len(cycle2)
	if next1 == start1 or next2 == start2:
		return pos1, pos2
	elif cycle1[next1].finish != cycle2[next2].finish:
		return pos1, pos2
	elif cycle1[next1].index != cycle2[next2].index:
		return pos1, pos2
	else:
		return getEndOfRepeat(cycle1, start1, next1, cycle2, start2, next2)

#########################################
## Split
#########################################

def splitCycles(cycleA, indexA, cycleB, indexB, stub):
	b1, b2 = getEndOfRepeat(cycleA, indexA, indexA, cycleB, indexB, indexB) 
	cycleA = cycleA.startAt(b1+1)
	cycleB = cycleB.startAt(b2+1)
	a1, a2 = getStartOfRepeat(cycleA, 0, 0, cycleB, 0, 0) 

	shortCut = Cycle([Edge(cycleA[a1].start, stub, 0)])
	if a1 % 2 == 1:
		shortCut.append(Edge(stub, stub, 0))
	shortCut.append(Edge(stub, cycleA[0].start, 0))

	cycleAPrime = Event(Cycle(cycleA[:a1] + shortCut, cycleA.value))
	cycleBPrime = Event(Cycle(cycleB[:a2] + shortCut, cycleB.value))
	newCycles = [cycleAPrime, cycleBPrime]

	# Note: a1 and a2 have the same parity because a1 + len(shortCut) and a2 + len(shortCut) both even numbers. 
	if a1 % 2 == 0: 
		residue = cycleB.value + cycleA.value
	else:
		residue = -(cycleB.value + cycleA.value)
	if abs(residue) > ROUNDING_ERROR:
		newCycles.append(Event(Cycle(cycleA[a1:] + shortCut.reverse(), residue)))
	return newCycles

#########################################
## Merge
#########################################

def mergeCycles(cycleA, indexA, cycleB, indexB):
	b1, b2 = getEndOfRepeat(cycleA, indexA, indexA, cycleB, indexB, indexB) 
	cycleA = cycleA.startAt(b1+1)
	cycleB = cycleB.startAt(b2+1)
	a1, a2 = getStartOfRepeat(cycleA, 0, 0, cycleB, 0, 0) 

	merged = Event(Cycle(cycleA[:a1]) + Cycle(cycleB[:a2]).reverse())
	merged.setRatio(cycleA.value)

	residue = cycleB.value + cycleA.value
	if abs(residue) > ROUNDING_ERROR:
		return [merged, Event(Cycle(cycleB, residue))]
	else:
		return [merged]

#########################################
## Ergodic moves
#########################################

def rescheduleRandomEvent(cactusHistory, history):
	if len(history.events) > 0:
		event = random.choice(history.events)
		cactusHistory.pop(history, event)
		addEvents(cactusHistory, history, [event])
	return history
	
def changeEventsInHistory(cactusHistory, history, overlap):
	cactusHistory.pop(history, overlap.localEvent)
	cactusHistory.pop(history, overlap.remoteEvent)

	edgeA = overlap.localEvent.cycle[overlap.localCut]
	edgeB = overlap.remoteEvent.cycle[overlap.remoteCut]

	if edgeA.nodes() == edgeB.nodes(): 
		cycleB = overlap.remoteEvent.cycle
	else:
		cycleB = overlap.remoteEvent.cycle.reverse()
		overlap.remoteCut = len(cycleB) - 1 - overlap.remoteCut

	if edgeA.start.chr == "None" or edgeA.finish.chr == "None":
		print 'MERGE'
		newEvents = mergeCycles(overlap.localEvent.cycle, overlap.localCut, cycleB, overlap.remoteCut)
	else:
		print 'SPLIT'
		newEvents = splitCycles(overlap.localEvent.cycle, overlap.localCut, cycleB, overlap.remoteCut, history.module.stub)
	addEvents(cactusHistory, history, newEvents)
	return history
	
def createNewHistory(cactusHistory, history, indirect):
	history2 = copy.copy(history)
	overlap = history2.overlap(indirect)
	cactusHistory.slideIn(history, history2)
	if overlap is None:
		print 'RESCHEDULE'
		return rescheduleRandomEvent(cactusHistory, history2)
	else:
		return changeEventsInHistory(cactusHistory, history2, overlap)

#############################################
## Unit test
#############################################
def main():
	G = avg.randomEulerianGraph(10)
	C = cactus.Cactus(G)
	N = normalized.NormalizedCactus(C)
	O = orientedCactus.OrientedCactus(N)
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print O
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	CH = cycleCover.initialHistory(O)
	print CH
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	N = random.choice(O.nets) 
	NH = createNewHistory(CH, CH.netHistories[N])
	NH.validate()

if __name__ == "__main__":
	main()
