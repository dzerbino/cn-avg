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

"""
Expectation Maximisation timing of the SNVs wrt CNVs
"""

import sys 
import random
import math
import numpy as np
import rpy2
from rpy2.robjects import FloatVector
from rpy2.robjects import IntVector
dchisq = rpy2.robjects.r.dchisq

from operator import mul

MAXITERATIONS = 100
TOLERANCE = 1e-1
MIN_PROB = 1e-10

########################################
## Generating configurations
########################################

def expandSliceTreeOptions_Solution(solution, phase, sliceTree):
	# Returns list -> tuple -> phase x sliceTree
	return solution + [(phase, sliceTree)]

def expandSliceTreeOptions_SliceTree(solutions, phase, sliceTree):
	# Return list -> list -> tuple -> phase x sliceTree
	return [expandSliceTreeOptions_Solution(solution, phase, sliceTree) for solution in solutions]

def expandSliceTreeOptions_Phase(solutions, phase, phaseLikelihoods):
	# Receives dictionary: sliceTree -> event -> snv -> likelihood
	# Return list -> list -> tuple -> phase x sliceTree
	return sum((expandSliceTreeOptions_SliceTree(solutions, phase, sliceTree) for sliceTree in phaseLikelihoods), [])

def generateSliceTreeOptions_Block2(blockLikelihoods):
	# Return list -> dict -> phase x sliceTree
	return map(dict, reduce(lambda X, Y: expandSliceTreeOptions_Phase(X, Y, blockLikelihoods[Y]), blockLikelihoods, [[]]))

def generateSliceTreeOptions_Block(likelihoods, block):
	# Return block x list -> dict -> phase x sliceTree
	return block, generateSliceTreeOptions_Block2(likelihoods[block]) 

def generateSliceTreeOptions(likelihoods):
	# Return block -> list -> dict -> phase x sliceTree
	return dict(generateSliceTreeOptions_Block(likelihoods, block) for block in likelihoods)
 
########################################
## EM convergence criterion
########################################

def chiSquare(X, Y):
	return sum((X - Y) * (X - Y) / np.maximum(Y, MIN_PROB))

def chiSquareTest(X, Y, snvCount):
	return dchisq(FloatVector([snvCount * chiSquare(X, Y)]), IntVector([len(X)-1]))[0]

def mean(V):
	return sum(V) / len(V)

def stop_EM(new, old):
	res = mean(np.abs(new - old) / np.maximum(old, MIN_PROB))
	print res
	return res < TOLERANCE

########################################
## Expectation Maximization of event parameters
########################################

def EM_eventParameters_Phase(likelihoods, eventParameters):
	# Returns array
	return likelihoods * eventParameters

def normalized(V):
	# Returns array
	return V / sum(V)

def EM_eventParameters_SNV(snv, likelihoods, chosenSliceTrees, eventParameters):
	# Returns array
	return normalized(sum(EM_eventParameters_Phase(likelihoods[phase][chosenSliceTrees[phase]][snv], eventParameters) for phase in chosenSliceTrees))

def EM_eventParameters_Block(snvs, likelihoods, chosenSliceTrees, eventParameters):
	# Returns array
	return sum(EM_eventParameters_SNV(snv, likelihoods, chosenSliceTrees, eventParameters) for snv in snvs)

def EM_eventParameters(blockSNVs, likelihoods, chosenSliceTrees, eventParameters, snvCount):
	# Returns array
	return sum(EM_eventParameters_Block(blockSNVs[block], likelihoods[block], chosenSliceTrees[block], eventParameters) for block in blockSNVs) / snvCount

########################################
## Pre-calculations for the choice of sliceTree
########################################

#### Computing Joint probs
def computeJointProbs_SNV(likelihoods, eventParameters):
	return likelihoods * eventParameters

def computeJointProbs_Phase(likelihoods, eventParameters): 
	return dict((snv, computeJointProbs_SNV(likelihoods[snv], eventParameters)) for snv in likelihoods)

def computeJointProbs(likelihoods, chosenSliceTrees, eventParameters):
	return dict((phase, computeJointProbs_Phase(likelihoods[phase][chosenSliceTrees[phase]], eventParameters)) for phase in likelihoods)

#### Computing sum of joint probs for each locus
def computeNormalizationFactors_SNV(factors, snv, jointProbs):
	if snv not in factors:
		factors[snv] = 0
	factors[snv] += sum(jointProbs[snv])
	return factors

def computeNormalizationFactors_Phase(factors, jointProbs):
	return reduce(lambda X, Y: computeNormalizationFactors_SNV(X, Y, jointProbs), jointProbs, factors)

def computeNormalizationFactors(jointProbs):
	return reduce(computeNormalizationFactors_Phase, jointProbs.values(), dict())

#### Computing posterior probs by normalization
def normalizeJointProbs_SNV(jointProbs, normalizationFactor):
	return jointProbs / normalizationFactor
	
def normalizeJointProbs_Phase(jointProbs, normalizationFactors):
	return dict((snv, normalizeJointProbs_SNV(jointProbs[snv], normalizationFactors[snv])) for snv in jointProbs)

def normalizeJointProbs(jointProbs, normalizationFactors):
	return dict((phase, normalizeJointProbs_Phase(jointProbs[phase], normalizationFactors)) for phase in jointProbs)

def computePosteriors(likelihoods, chosenSliceTrees, eventParameters):
	jointProbs = computeJointProbs(likelihoods, chosenSliceTrees, eventParameters)
	normalizationFactors = computeNormalizationFactors(jointProbs)
	return normalizeJointProbs(jointProbs, normalizationFactors)

########################################
## Choice of sliceTree
########################################

def sliceTreePrior(chosenSliceTrees):
	return reduce(mul, (X.prior for X in chosenSliceTrees.values()), 1)
	
def EM_sliceTree_prior(chosenSliceTrees):
	return math.log(np.maximum(sliceTreePrior(chosenSliceTrees), MIN_PROB))

def EM_sliceTree_SNV(likelihoods, posteriors):
	return sum(posteriors * np.log(np.maximum(likelihoods, MIN_PROB)))

def EM_sliceTree_Phase(likelihoods, posteriors):
	return sum(EM_sliceTree_SNV(likelihoods[snv], posteriors[snv]) for snv in likelihoods)

def EM_sliceTree_likelihood(likelihoods, posteriors, sliceTreeOption):
	return sum(EM_sliceTree_Phase(likelihoods[phase][sliceTreeOption[phase]], posteriors[phase]) for phase in likelihoods)

def EM_sliceTree_SliceTreeOption(likelihoods, posteriors, sliceTreeOption):
	# Returns sliceTreeOption x float
	return sliceTreeOption, EM_sliceTree_prior(sliceTreeOption) + EM_sliceTree_likelihood(likelihoods, posteriors, sliceTreeOption)

def argmax(V):
	M = max([X[1] for X in V])
	return random.choice(filter(lambda X: X[1] == M, V))[0]

def EM_sliceTree_Block(likelihoods, sliceTreeOptions, chosenSliceTrees, eventParameters):
	# Returns sliceTrees
	posteriors = computePosteriors(likelihoods, chosenSliceTrees, eventParameters)
	return argmax([EM_sliceTree_SliceTreeOption(likelihoods, posteriors, sliceTreeOption) for sliceTreeOption in sliceTreeOptions])

def EM_sliceTree(likelihoods, sliceTreeOptions, chosenSliceTrees, eventParameters):
	# Returns block -> sliceTreeOption
	return dict((block, EM_sliceTree_Block(likelihoods[block], sliceTreeOptions[block], chosenSliceTrees[block], eventParameters)) for block in likelihoods)
		
#######################################
## Master function
#######################################

def EMStep(blockSNVs, likelihoods, sliceTreeOptions, oldParameters, snvCount):
	eventParameters, chosenSliceTrees = oldParameters
	newEventParameters = EM_eventParameters(blockSNVs, likelihoods, chosenSliceTrees, eventParameters, snvCount)
	newSliceTrees = EM_sliceTree(likelihoods, sliceTreeOptions, chosenSliceTrees, eventParameters)
	return newEventParameters, newSliceTrees

def EMIterations(data, index, likelihoods, blockSNVs, snvCount, sliceTreeOptions):
	oldParameters, stop = data
	if stop:
		return data
	else:
		parameters = EMStep(blockSNVs, likelihoods, sliceTreeOptions, oldParameters, snvCount)
		return parameters, stop_EM(parameters[0], oldParameters[0])

def EM2(likelihoods, events, blockSNVs, sliceTreeOptions, snvCount):
	startingEventParameters = normalized(np.random.random(len(events) + 1))
	startingSliceTrees = dict((block, random.choice(sliceTreeOptions[block])) for block in blockSNVs)
	startingParameters = startingEventParameters, startingSliceTrees
	return reduce(lambda X, Y: EMIterations(X, Y, likelihoods, blockSNVs, snvCount, sliceTreeOptions), range(MAXITERATIONS), (startingParameters, False))[0]

def EM(likelihoods, events, blockSNVs):
	print 'EM sampling'
	snvCount = sum(len(snvs) for snvs in blockSNVs.values())
	print "Counted %i snvs" % snvCount
	sliceTreeOptions = generateSliceTreeOptions(likelihoods)
	res = [EM2(likelihoods, events, blockSNVs, sliceTreeOptions, snvCount) for X in range(10)]
	print '>>>>>>>>>>>'
	for X in res:
		for Y in res:
			if Y is not X:
				#print chiSquareTest(X[0], Y[0], snvCount)
				#print "\n".join(map(lambda X: "\t".join(map(str, X)), zip(X[0], Y[0])))
				print "COR", rpy2.robjects.r.cor(FloatVector(np.log(X[0])), FloatVector(np.log(Y[0])))[0]
				print "LOG", rpy2.robjects.r.cor(FloatVector(X[0]), FloatVector(Y[0]))[0]
	assert False

########################################
## Test
########################################
def main():
	assert False
	
if __name__ == "__main__":
	main()
