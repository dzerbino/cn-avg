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
import cbs.cbs as cbs
from cnv import CNV

def mergeCNVs_Forward(A, B):
	for i in range(A.ploidy()):
		A.val[i] = (A.val[i] * A.length() + B.val[i] * B.length() / 2) / (A.length() + B.length() / 2)
	assert A.val > 0
	if B.softStart is not None:
		A.finish = (B.softStart + B.finish) / 2 
	else:
		A.finish = B.finish
	if B.softFinish is not None:
		A.softFinish = (B.start + B.softFinish) / 2
	else:
		A.softFinish = B.start
	A.numMarks += B.numMarks / 2

def mergeCNVs_Backward(A, B):
	for i in range(A.ploidy()):
		B.val[i] = (A.val[i] * A.length() / 2 + B.val[i] * B.length()) / (A.length() / 2 + B.length())
	assert B.val > 0
	if A.softFinish is not None:
		B.start = (A.start + A.softFinish) / 2
	else:
		B.start = A.start

	if A.softStart is not None:
		B.softFinish = (A.softStart + A.finish) / 2
	else:
		B.softFinish = A.finish
	B.numMarks += A.numMarks / 2

def filterCNVs(cnvs):
	filtered = []
	for index in range(len(cnvs)):
		cnv = cnvs[index]
		if cnv.numMarks < 5:
			if len(filtered) > 0:
				valsA = sum(filtered[-1].val)
			else:
				valsA = 0
			valsB = sum(cnv.val)
			if len(filtered) > 0 and filtered[-1].chr == cnv.chr and filtered[-1].ploidy() == cnv.ploidy() and valsB > 0 and valsA / float(valsB) > 0.5 and valsA / float(valsB) < 2:
				mergeCNVs_Forward(filtered[-1], cnv)
			else:
				filtered += [cnv]

			if index + 1 < len(cnvs) and cnvs[index+1].chr == cnv.chr and cnv.ploidy() == cnvs[index+1].ploidy():
				mergeCNVs_Backward(cnv, cnvs[index+1])
		else:
			filtered += [cnv]
	return filtered

def glueCNV(cnvs, index):
	prev = cnvs[index-1]
	cnv = cnvs[index]
	if prev.chr == cnv.chr and prev.finish < cnv.start and prev.finish > cnv.start - 50:
		prev.finish = cnv.start

def glueCNVs(cnvs):
	map(lambda X: glueCNV(cnvs, X), range(1, len(cnvs))) 

def computeBorderMargins(segments, sortedCNVs):
	print "Computing CNV border margins"
	cnvindex = 0
	for segment in segments:
		while len(sortedCNVs) > 0 and sortedCNVs[cnvindex] < segment:
			cnvindex += 1
			assert cnvindex < len(sortedCNVs)
		
		if cnvindex > 0 and sortedCNVs[cnvindex - 1].chr == segment.chr:
			segment.start = sortedCNVs[cnvindex - 1].start
		else:
			segment.start = sortedCNVs[cnvindex].start
		segment.softStart = sortedCNVs[cnvindex].finish

		while cnvindex < len(sortedCNVs) and sortedCNVs[cnvindex] == segment:
			cnvindex += 1
	
		# Backing up to last valid index
		assert cnvindex > 0
		cnvindex -= 1
		segment.softFinish = sortedCNVs[cnvindex].start
		if segment.softFinish == segment.start:
			segment.softStart = sortedCNVs[cnvindex].start + 1
			segment.softFinish = sortedCNVs[cnvindex].finish - 1
		if segment.softFinish <= segment.softStart:
			segment.softFinish = segment.softStart + 1
		if cnvindex < len(sortedCNVs) - 1 and sortedCNVs[cnvindex + 1].chr == segment.chr:
			segment.finish = sortedCNVs[cnvindex + 1].finish
		else:
			segment.finish = sortedCNVs[cnvindex].finish
		assert segment.finish > segment.start
		if segment.softFinish >= segment.finish:
			print cnvindex
			print len(sortedCNVs) - 1
			print sortedCNVs[cnvindex + 1].chr
			print segment.chr
			print sortedCNVs[cnvindex].start
			print segment.softStart
			print segment.softFinish
			if cnvindex < len(sortedCNVs) - 1:
				print sortedCNVs[cnvindex + 1].finish
			print sortedCNVs[cnvindex].finish
		segment.validate()

	return segments

def uniq(list, elem):
	if len(list) == 0 or list[-1] != elem:
		list.append(elem)
	return list

###########################################
## Master function
###########################################

def median(X):
	return (X.start + X.finish)/2

def segmentCNVs(cnvs):
	print "Segmentation of CNV data"
	sortedCNVs = sorted(cnvs)	
	glueCNVs(sortedCNVs)

	chrom = [X.chr for X in sortedCNVs]
	pos = map(median, sortedCNVs)
	vals = [X.val[0] for X in sortedCNVs]

	cbsCNVS = cbs.run(chrom, pos, vals)
	sortedCBSCNVs = sorted(cbsCNVS)
	filtered = filterCNVs(sortedCBSCNVs)

	return computeBorderMargins(filtered, sortedCNVs)
	
###########################################
## Unit test
###########################################
def main():
	region1 = CNV("chr1", 1000, 2000, [1.0], "A")
	region2 = CNV("chr1", 2000, 3000, [1.0], "B")
	region3 = CNV("chr1", 3000, 4000, [5.0], "C")
	region4 = CNV("chr1", 4000, 5000, [5.0], "D")
	print "\n".join(map(str, segmentCNVs([region1, region2, region3, region4])))

if __name__ == "__main__":
	main()
