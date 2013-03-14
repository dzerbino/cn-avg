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

"""Sampling through histories at the global level"""

import cnavg.avg.graph as avg
import cnavg.cactus.graph as cactus
import cnavg.cactusSampling.sampling as normalizedCactus
import cnavg.cactus.oriented as oriented
import cnavg.history.history as history
import sampleModuleCycles
import cycleCover
import random
import math
import copy
import time
import cPickle as pickle

import cnavg.history.ordered

TEMPERATURE = 1
DIRECT_CUT = 0

TIMER_END = None
TIMER_LENGTH = 36000

#############################################
## Select optimal solutions
#############################################

def rearrangementCosts(histories):
	return map(lambda X: X.rearrangementCost(), histories)

def minimumRearrangementCost(histories):
	return min(rearrangementCosts(histories))

def optimalSolutions(histories):
	costs = map(lambda X: X.rearrangementCost(), histories)
	target = min(costs)
	optimals = filter(lambda X: X[0] == target, zip(costs, histories))
	return map(lambda X: X[1], optimals)

########################################
## Propagation through the graph 
########################################

def changedCNVs(new, old):
	if len(new) != len(old):
		return True

	return any(abs(X[0][0].value - X[1][0].value) > 1e-6 for X in zip(new, old))

def modifyCactusHistory_Net(oldhistory, history, net):
	for chain in history.cactus.nets2Chains[net]:
		if changedCNVs(history.chainCNVs[chain], oldhistory.chainCNVs[chain]):
			modifyCactusHistory_Chain(oldhistory, history, chain)
	return history

def modifyCactusHistory_Chain(oldhistory, history, chain):
	for net in history.cactus.chains2Nets[chain]:
		cycleCover.seedHistory(history, net, history.chainCNVs[chain])
		modifyCactusHistory_Net(oldhistory, history, net)
	return history

def modifyCactusHistory_Initiate(oldhistory, history, net, indirect):
	newLocalHistory = sampleModuleCycles.createNewHistory(history, history.netHistories[net], indirect)
	history.update(net, newLocalHistory)
	history.updateCNVs(net, newLocalHistory)
	for chain in history.cactus.nets2Chains[net]:
		if changedCNVs(history.chainCNVs[chain], oldhistory.chainCNVs[chain]):
			modifyCactusHistory_Chain(oldhistory, history, chain)
	return history

########################################
## Choose a net to modify
########################################

def enumerateNets(history):
	return [(X[1].density(), X[0]) for X in history.netHistories.iteritems()]

def weightedChoice(weightedOptions):
	total = sum(X[0] for X in weightedOptions)
	target = random.uniform(0, total)
	for option in weightedOptions:
		if option[0] >= target:
			return option[1]
		target -= option[0]
	assert False

def chooseNet(history):
	return weightedChoice(enumerateNets(history))

########################################
## MC exploration 
########################################

def modifyCactusHistory(history, indirect):
	if len(history.netHistories) > 0:
		return modifyCactusHistory_Initiate(history, copy.copy(history), chooseNet(history), indirect)
	else:
		return history

def MCTest(newScore, oldScore, temperature):
	print 'TEST',oldScore,newScore
	return newScore <= oldScore or random.random() < math.exp((oldScore - newScore) / temperature)

def chooseNewHistory(history, temperature, depth=0, index=0):
	newHistory = modifyCactusHistory(history, index > DIRECT_CUT)
	if newHistory is None:
		return None
	elif MCTest(newHistory.rearrangementCost(), history.rearrangementCost(), temperature):
		print 'HISTORY ', index, depth, len(history.parent), newHistory.rearrangementCost(), newHistory.errorCost(), time.asctime()
		# Memory optimization trick -> Frees unused matrices and overlap tables
		history.embalmDeadHistories(newHistory)
		return newHistory
	else:
		print 'REFUSE'
		return chooseNewHistory(history, temperature * 1.1, depth=depth+1, index=index)

MINHISTORY = None

def addNewHistory(histories, index, file=None, stats=None, braney=None, tree=None):
	# Just do not process if run time > TIMER_LENGTH
	global TIMER_END
	if time.time() > TIMER_END:
		return histories
	
	# Generate new history from the last one
	global TEMPERATURE
	newHistory = chooseNewHistory(histories[-1], TEMPERATURE, depth=1, index=index)
	if stats is not None:
		stats.write("%s\n" % newHistory.stats())
	if braney is not None:
		braney.write("%s\n" % cnavg.history.ordered.prettify(newHistory, index))
	if tree is not None:
		tree.write("%s\n" % newHistory.newick())
	if file is None: 
		histories.append(newHistory)
		return histories 
	else:
		# Really, who need to pickle files when there's a perfectly good zipped braney file next to it?
		#pickle.dump(newHistory, file)
		# Cheap trick: only store minimal history
		global MINHISTORY
		if MINHISTORY is None or newHistory.rearrangementCost() < MINHISTORY.rearrangementCost():
			MINHISTORY = newHistory
		return [newHistory]

########################################
## Master function
########################################

def sample(cactusHistory, size, file=None, stats=None, braney=None, tree=None):
	print 'Sampling history space of Cactus graph'
	global TIMER_END
	TIMER_END = time.time() + TIMER_LENGTH
	res = reduce(lambda X, Y: addNewHistory(X, Y, file, stats, braney, tree), range(size), [cactusHistory])
	if file is not None:
		pickle.dump(res[0], file)
	return res

########################################
## Unit test
########################################

def main():
	G = avg.randomEulerianGraph(10)
	C = cactus.Cactus(G)
	N = normalizedCactus.NormalizedCactus(C)
	orientedCactus = oriented.OrientedCactus(N)
	print orientedCactus
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	res = sample(cycleCover.initialHistory(orientedCactus), 20)[-1]
	print res
	res.validate()

if __name__ == "__main__":
	main()
