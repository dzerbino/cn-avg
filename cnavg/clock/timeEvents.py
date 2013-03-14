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

""" Merging the SNV and CNV data"""

import vcf
import random
import cnavg.history.debug

from cnavg.basics.coords import Position
from cnavg.clock.EM import EM
from cnavg.clock.cnvTree import CNVTree

########################################
## Computing likelihoods
########################################

def computeLikelihoods_Block(block, blockTrees, blockSNVs, events, cactus):
	# Returns dictionary: phase -> sliceTree -> snv -> array -> likelihood
	# Note: + 1 is for original haploid copy, since copynumber only records a change 
	snvs = blockSNVs[block]
	counts = [snv.mutant for snv in snvs]
	totals = [snv.total for snv in snvs]
	totalCopyNumber = sum(block.copynumber(cactus, phase) + 1 for phase in range(block.ploidy(cactus)))
	return dict((phase, blockTrees[block][phase].likelihoods(counts, totals, snvs, events, (block.copynumber(cactus, phase) + 1)/ totalCopyNumber)) for phase in range(block.ploidy(cactus)))

def computeLikelihoods(blockTrees, blockSNVs, events, cactus):
	print 'Computing likelihoods'
	# Returns dictionary: block -> phase -> sliceTree -> snv -> array -> likelihood
	return dict((block, computeLikelihoods_Block(block, blockTrees, blockSNVs, events, cactus)) for block in blockSNVs) 

########################################
## Block molecular SNV trees
########################################

def cactusHistoryOverlapTable(cactusHistory):
	overlapTable = dict()
	for e in cactusHistory.parent:
		overlapTable = overlap.addEventToOverlapTable(overlapTable, e)
	return overlapTable

def phasedBlockIndex(block, phase):
	return (min(block.nodes), max(block.nodes), phase)

def phasedBlockEvents(block, phase, overlapTable):
	index = phasedBlockIndex(block, phase)
	if index in overlapTable:
		return overlapTable[(min(block.nodes), max(block.nodes), phase)]
	else:
		return None

def createBlockTrees_Block(block, overlapTable, cactusHistory):
	# 1 is assuming normal ploidy of the germline genome
	return dict((phase, CNVTree(cactusHistory, phasedBlockEvents(block, phase, overlapTable), 1)) for phase in range(block.ploidy(cactusHistory.cactus)))

def createBlockTrees(cactusHistory, blocks):
	print 'Creating block trees'
	overlapTable = cactusHistoryOverlapTable(cactusHistory)
	return dict((block, createBlockTrees_Block(block, overlapTable, cactusHistory)) for block in blocks)

########################################
## Parsing SNVs
########################################

class SNV(Position):
	def __init__(self, chr, pos, mutant, total, record):
		super(SNV, self).__init__(chr, pos)
		self.mutant = int(mutant)
		self.total = int(total)
		self.record = record

	def __str__(self):
		return "\t".join([super(SNV, self).__str__(), str(self.mutant), str(self.total)])

def filteredSNV(record):
	return 'VT' in record.INFO and record.INFO['VT'] == 'SNP' and 'DB' not in record.INFO and record.FILTER is None and 'SS' in record.INFO and record.INFO['SS'] == 'Somatic'

def readVCFRecord(snvs, record):
	if filteredSNV(record):
		call = record.genotype('PRIMARY')
		snvs.append(SNV(record.CHROM, record.POS, int(call.data['DP'] * call.data['FA']), call.data['DP'], record))
	return snvs

def readVCFFile(filename):
	return reduce(readVCFRecord, vcf.Reader(open(filename), prepend_chr=True), [])
	
########################################
## Assigning SNVs to blocks
########################################

def assignSNV(data, snv):
	blockSNVs, blockList = data

	# Skip blocks
	while len(blockList) > 0 and snv > max(blockList[0].nodes):
		blockList.pop(0)

	# Assign if block found
	if len(blockList) > 0 and snv >= min(blockList[0].nodes):
		block = blockList[0]
		if block not in blockSNVs:
			blockSNVs[block] = []
		blockSNVs[block].append(snv)

	return data

def assignSNVs(blocks, snvs):
	return reduce(assignSNV, sorted(snvs), (dict(), sorted(blocks)))[0]

########################################
## Newick Tree
########################################

def newickLeaf(event, priors):
	return str(id(event)) + ":" + str(priors[None])

def newickFork(cactusHistory, event, priors):
	return ",".join([newick(cactusHistory, child, priors) for child in cactusHistory.children[event]] + [newickLeaf(event, priors)])

def newick(cactusHistory, event, priors):
	return "(" + newickFork(cactusHistory, event, priors) + "):" + str(priors[event])

def newicks(cactusHistory, priors):
	return "\n".join(newick(cactusHistory, event, priors) + ";" for event in cactusHistory.roots)

########################################
## Master function
########################################

def timeEvents(cactusHistory, snvs):
	blocks = [B for C in cactusHistory.cactus.chains for B in C if B.ploidy(cactusHistory.cactus)]
	blockTrees = createBlockTrees(cactusHistory, blocks)
	blockSNVs = assignSNVs(blocks, snvs)
	events = cactusHistory.parent.keys()
	likelihoods = computeLikelihoods(blockTrees, blockSNVs, events, cactusHistory.cactus)		
	return EM(likelihoods, events, blockSNVs)

##############################################
## Unit testing
##############################################

def main():
	sys.setrecursionlimit(5000)
	print 'Looking for pickled file'
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
	assert H is not None

	print 'Cleaning up graph'
	FH = flattened.flattenGraph(H)
	S = FH.simplifyStubsAndTrivials()
	F = S.removeLowRatioEvents(debug.RATIO_CUTOFF)
	print 'Reading VCF'
	snvs = readVCFFile(sys.argv[2])
	print 'Timing events'
	priors, assignments = timeEvents(F, snvs)
	#print newicks(F, priors)
	O = OrderedHistory(F)
	for event in O.parent:
		print event.ratio, sum(priors[parent] for parent in F.ancestors[event]), O.ordering.depth[event]

if __name__ == "__main__":
	import sys 
	import cPickle as pickle
	import cnavg.history.overlap as overlap
	import cnavg.history.flattened as flattened
	from cnavg.postprocess.annotate import annotate
	from cnavg.history.ordered import OrderedHistory
	main()
