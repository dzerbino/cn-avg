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
import re
import vcf

import cnavg.avg.graph as avg
import breakendGraph

###############################################
## Parsing the chromosome lengths
###############################################

def parseChromLengths(file):
	res = dict()
	for line in open(file):
		items = line.strip().split()
		res[items[0]] = int(items[1])
	return res 


###############################################
## Parsing the VCF file
###############################################

def parseVCFBreakend(breakends, record):
	if record.FILTER != 'PASS':
		return breakends

	if "SVTYPE" not in record.INFO or record.INFO["SVTYPE"] != "BND":
		return breakends

	assert len(record.ALT) > 0
	assert all(X.reconnects for X in record.ALT)
	breakend = breakendGraph.Breakend(record.CHROM, record.POS, record.ALT[0].orientation, record.ID)
	breakend.remoteChr = [X.chr for X in record.ALT]
	breakend.remotePos = [X.pos for X in record.ALT]
	breakend.remoteOrientation = [X.remoteOrientation for X in record.ALT]
	breakend.adjacency_cov = [ 0 for x in breakend.remoteChr]
	breakend.adjacency_cov.append(0)

	if "MATEID" in record.INFO:
		breakend.mates = record.INFO["MATEID"]
	if "PARID" in record.INFO:
		breakend.partner = record.INFO['PARID']

	assert len(record.samples) > 1

	depths = sampleData[format['BDP']]
	for allele in genotypes:
		breakend.adjacency_cov[genotypes[allele]] += float(depths[allele])

	breakend.segment_cov = int(sampleData[format['DP']])
	
	breakends[breakend.ID] = breakend
	return breakends

def parseVCFFile(vcffile):
	return reduce(parseVCFBreakend, vcf.Reader(open(vcffile)), dict())

########################################################
## Replace name by pointers
########################################################

def ConsolidateBreakend(breakend, breakends):
	if breakend.partner is not None:
		breakend.partner = breakends[breakend.partner] 

	if breakend.mates is not None:
		breakend.mates = map(lambda X: breakends[X], breakend.mates) 

def AlignPointers(breakends):
	map(lambda X: ConsolidateBreakend(breakends[X], breakends), breakends)

########################################################
## Master function
########################################################

def parse(vcffile, lengthsFile):
	breakends = dict()
	AlignPointers(breakends)
	graph = breakendGraph.BreakendGraph(sorted(breakends.values()))
	graph.consolidate()
	graph.validate()
	graph.lengths = parseChromLengths(lengthsFile)
	return graph

########################################################
## Test function
########################################################
def main():
	breakends = parse('../data/test.vcf', '../data/toto.lengths')
	print breakends

if __name__ == "__main__":
	main()
