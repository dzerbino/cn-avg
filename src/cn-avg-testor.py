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

import cnavg.simulator.simulator as simulator
import cnavg.avg.balanced as balancedAVG
import cnavg.cactus.graph as cactus
import cnavg.cactusSampling.sampling as normalized
import cnavg.cactus.oriented as oriented
import cnavg.cactus.balanced as balancedCactus
import cnavg.historySampling.cycleCover as cycleCover
import cnavg.historySampling.sampleGraphCycles as sampleGraphCycles


import cnavg.history.debug
cnavg.history.debug.DEBUG = True
cnavg.history.debug.PLOIDY = 1


def _sampleAVG(avg, size):
	""" Simplified pipeline """
	C = cactus.Cactus(avg)
	NC = normalized.NormalizedCactus(C)
	OC = oriented.OrientedCactus(NC)
	H = cycleCover.initialHistory(OC)
	return sampleGraphCycles.sample(H, size)

def _compareHistories(realHistory, sampledHistories):
	""" Reporting """
	X = map(lambda X: X.rearrangementCost(), sampledHistories)
	m = min(X)
	return realHistory.cost(), m, len(filter(lambda x : x == m, X)), X

def testHistory(length, maxDepth, iterations):
	"""Tests the efficiency of the CN-AVG pipeline on random evolutionary histories"""
	realHistory = simulator.RandomWeightedHistory(length, maxDepth)
	sampledHistories = _sampleAVG(realHistory.avg(), iterations)
	return "\t".join(map(str, _compareHistories(realHistory, sampledHistories)))

def main():
	print "\n".join([testHistory(100, 10, 100) for X in range(100)])

if __name__ == "__main__":
	main()
