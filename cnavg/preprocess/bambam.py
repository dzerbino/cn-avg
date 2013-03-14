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

"""Parsing BamBam output"""

import sys
import re
import os

from cnavg.basics.coords import Position
from breakend import Breakend
from cnv import CNV
from breakendGraph import BreakendGraph
from segmentCNVs import segmentCNVs
import cbs.cbs

import rpy2.robjects
from rpy2.robjects import IntVector

READLENGTH = 50
BUFFER_LENGTH = 1000

#############################################
## Chromosomal lengths 
#############################################

def parseChromLengths(file):
	print "\tParsing chromosome lengths"

	res = dict()
	for line in file:
		items = line.strip().split()
		res[items[0]] = int(items[1])
	return res 
	
#############################################
## Breakends stuff 
#############################################

def parseBreaksLine(line, breakends):
	items = line.strip().split()
	if float(items[11]) < 0:
		return
	if items[0][:3] != "chr":
		items[0] = "chr" + items[0]
	breakend = Breakend(items[0], items[1], items[3] == '-', "BND" + str(len(breakends)) + ":" + items[0] + ":" + items[1] + items[3] + ">>" + items[4] + ":" + items[5] + items[7], finish=int(items[2]))
	breakend.mates = [None]

	if items[4][:3] != "chr":
		items[4] = "chr" + items[4]
	breakend.remoteChr = [items[4]]
	breakend.remoteStart = [int(items[5])]
	breakend.remoteFinish = [int(items[6])]
	breakend.remoteOrientation = [items[7] == '-']

	breakend.germline = (items[8] == 'GERMLINE')
	breakend.adjacency_cov = [-1, float(items[9])]

	breakends.append(breakend)

def parseBreaksLineNew(line, breakends):
	items = line.strip().split()
	if float(items[8]) < 0:
		return
	coords = items[0].split(':')
	chr = coords[0]
	if chr[:3] != "chr":
		chr = "chr" + chr
	pos = coords[1].split('-')
	breakend = Breakend(chr, pos[0], items[2] == '-', "BND" + str(len(breakends)) + ":" + items[0] + ":" + items[1], finish=int(pos[1]))
	breakend.mates = [None]

	coords2 = items[1].split(':')
	chr2 = coords2[0]
	if chr2[:3] != "chr":
		chr2 = "chr" + chr2
	pos2 = coords2[1].split('-')
	breakend.remoteChr = [chr2]
	breakend.remoteStart = [int(pos2[0])]
	breakend.remoteFinish = [int(pos2[1])]

	breakend.remoteOrientation = [items[3] == '-']
	breakend.adjacency_cov = [-1, float(items[5])]
	breakend.germline = (items[4] == 'germline')
	breakends.append(breakend)

def overlappingGermlineSomatic(A, B):
	return (A is not None and B is not None and A == B and A.remoteChr[0] == B.remoteChr[0] and A.remoteStart[0] < B.remoteFinish[0] and A.remoteFinish[0] > B.remoteStart[0] and A.remoteOrientation[0] == B.remoteOrientation[0] and A.germline != B.germline)

def filterOutGermlineBreakends(data, breakend):
	somatics, germline = data
	if breakend.germline: 
		if len(somatics) > 0 and overlappingGermlineSomatic(somatics[-1], breakend):
			somatics.pop(-1)
		return somatics, breakend
	else:
		if not overlappingGermlineSomatic(germline, breakend):
			somatics.append(breakend)
		return somatics, germline

def removeGermlineBreakends(breakends):
	return reduce(filterOutGermlineBreakends, breakends, ([], None))[0]

def removeNonChromosomalBreakpoints(graph, lengths):
	return BreakendGraph(filter(lambda X: X.chr in lengths and X.remoteChr[0] in lengths, graph))

def parseBreaksFile(file):
	print "\tParsing breakends"
	breakends = BreakendGraph()
	for line in file:
		parseBreaksLine(line, breakends)
	breakends.sort()
	filtered = removeGermlineBreakends(breakends)
	return BreakendGraph(filtered)
	
#############################################
## CNV stuff
#############################################

def parseCNVLine(items, cnvs):
	if items[1][:3] != "chr":
		items[1] = "chr" + items[1]
	cnvs.append(CNV(items[1], items[2], items[3], [float(items[4])], ""))
	cnvs[-1].base = int(items[5])
	return cnvs

def parseBBLine(cnvs, line):
	items = line.strip().split()
	if items[0] == 'CNV':
		return parseCNVLine(items, cnvs)
	else:
		return cnvs

def parseBBFile(cnvs, file):
	return reduce(parseBBLine, open(file), cnvs)

def parseCNVData(files):
	print "\tParsing CNV data"
	assert len(files) > 1
	return sorted(reduce(parseBBFile, files, []))

#############################################
## Adding buffers to replace missing colums
#############################################

def addMissingCNV(cnvs, cnv):
	if len(cnvs) == 0:
		return [cnv]
	while cnv.chr == cnvs[-1].chr and cnv.start - cnvs[-1].finish > BUFFER_LENGTH:
		cnvs.append(CNV(cnvs[-1].chr, cnvs[-1].finish, cnvs[-1].finish + BUFFER_LENGTH, [0], "EMPTY"))
	cnvs.append(cnv)
	return cnvs

def addMissingCNVs(cnvs):
	return reduce(addMissingCNV, cnvs, [])

#############################################
## Phasing stuff
#############################################

def getAlleleFrequencies_Line(results, line, cnvsByChrom, offsets):
	items = line.strip().split('\t')
	if items[0] != 'HAP':
		return results

	# Find covering CNV
	chr = items[1]
	if len(chr) < 3 or chr[:3] != "chr":
		chr = "chr" + chr
	pos = int(items[2])
	snp = Position(chr, pos)
	if chr not in offsets:
		print offsets.keys()
		print chr
		assert False
		return results
	cnv = cnvsByChrom[chr][offsets[chr]] 
	while cnv < snp:
		if offsets[chr] >= len(cnvsByChrom[chr]) - 1:
			return results
		cnv = cnvsByChrom[chr][offsets[chr] + 1]
		if cnv.chr != chr:
			return results
		offsets[chr] += 1
	if cnv > snp:
		return results

	# Store results
	germline_haplotype = re.split('[ /:]', items[6])
	germline_calls = germline_haplotype[0::3]
	somatic_haplotype = re.split('[ /:]', items[7])
	somatic_calls = dict(zip(somatic_haplotype[0::3], somatic_haplotype[1::3]))
	selected = []
	for call in germline_calls:
		if call not in somatic_calls:
			selected.append(0)
		else:
			selected.append(int(somatic_calls[call]))
	results[cnv].append(sorted(selected))
	return results

def getAlleleFrequencies_File(results, file, cnvsByChrom, offsets):
	return reduce(lambda r,l: getAlleleFrequencies_Line(r, l, cnvsByChrom, offsets), open(file), results)

rpy2.robjects.r('''
	F <- function(X, Y, variance) { 
		A = c(X[1:(length(X)/2)], Y[(length(X)/2 + 1):length(X)])
		B = c(Y[1:(length(X)/2)], X[(length(X)/2 + 1):length(X)])
		C = A + B
		var(A[C > 0]/C[C > 0])
	}
	''')
	
def tTest(X, Y, variance):
	if len(X) < 2:
		return False
	# Index 0 is to pull out the value from its list wrapper
	try:
		return rpy2.robjects.r.F(IntVector(X),IntVector(Y))[0] > 1.5 * variance
	except rpy2.rinterface.RRuntimeError:
		return False

def splitPhases(cnv, results, variance):
	values = results[cnv]
	majority = map(max, values)
	minority = map(min, values)
	if tTest(majority, minority, variance):
		ratios = [float(X[0])/(X[0] + X[1]) for X in zip(majority, minority) if (X[0] + X[1]) > 0]
		meanRatio = sum(ratios) / len(ratios)
		value = cnv.val[0]
		majVal = value * meanRatio
		minVal = value * (1 - meanRatio)
		newValues = [majVal, minVal]
		return CNV(cnv.chr, cnv.start, cnv.finish, newValues, cnv.name, softStart = cnv.softStart, softFinish = cnv.softFinish)
	else:
		return cnv

def splitByChroms(cnvsByChrom, cnv):
	if cnv.chr not in cnvsByChrom:
		cnvsByChrom[cnv.chr] = []
	cnvsByChrom[cnv.chr].append(cnv)
	return cnvsByChrom

def computeAllelicVariance(cnvMarkers):
	markers = sum(cnvMarkers.values(), [])
	valid = filter(lambda X: sum(X) > 0, markers)
	shuffled = [float(X[0])/sum(X) for X in valid[0:int(len(valid)/2)]] + [float(X[1])/sum(X) for X in valid[int(len(valid)/2):len(valid)]]
	expected = sum(shuffled) / len(shuffled)
	return sum(X*X for X in shuffled)/len(shuffled) - expected * expected

def phaseCNVs(files, cnvs, chroms):
	print "Phasing CNVs"
	offsets = dict((X, 0) for X in chroms)
	cnvsByChrom = reduce(splitByChroms, cnvs, dict())
	print '\tReading Heterozygous counts'
	cnvMarkers = reduce(lambda r,f: getAlleleFrequencies_File(r, f, cnvsByChrom, offsets), files, dict((X,[]) for X in cnvs))
	print '\tCalculating allelic ratio variance'
	variance = computeAllelicVariance(cnvMarkers)
	print 'Variance = ', variance
	print '\tT-tests across %i blocks' % len(cnvs)
	return map(lambda X: splitPhases(X, cnvMarkers, variance), cnvs)

#############################################
## Master function
#############################################

def parse(bbfiles, breaksfile, lengthsfile, snpsfiles=None):
	print "Parsing BamBam data"

	# Reading files
	lengths = parseChromLengths(lengthsfile)
	breakends = parseBreaksFile(breaksfile)
	breakends = removeNonChromosomalBreakpoints(breakends, lengths)
	breakends.consolidate()
	breakends.lengths = lengths
	print "Found %i breakends" % len(breakends) 

	# CNV data
	cnvs = segmentCNVs(parseCNVData(bbfiles))
	print "Found %i cnv regions" % len(cnvs) 

	# Quick normalization
	meanCov = sum(sum(X.val) * X.length() for X in cnvs) * READLENGTH / float(sum(X.length() for X in cnvs))
	for breakend in breakends:
		breakend.adjacency_cov = [X / meanCov for X in breakend.adjacency_cov]

	# Putting it all together
	if snpsfiles is None:
		cnvs = phaseCNVs(bbfiles, cnvs, lengths.keys())
	else:
		cnvs = phaseCNVs(snpsfiles, cnvs, lengths.keys())
	breakends.incorporateCNVs(cnvs)
	return breakends
		
#############################################
## Unit Test
#############################################

def main():
	print "\n".join(map(str, addMissingCNVs([CNV("1", 10, 20, [1], ""), CNV("2",10,20,[1],""), CNV("2", 2000, 3000, [1], ""), CNV("2", 6000, 7000, [1], "")])))
	return
	print parse("../data/test.bb", "../data/test.breaks", "../data/toto.lengths")

if __name__ == "__main__":
	main()
