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
import subprocess
import tempfile
import os
import os.path
from cnavg.preprocess.cnv import CNV

chromosomeNames = map(lambda X: "chr" + X, map(str, range(1, 23)) + ['X', 'Y'])

################################################
## Conversion functions for input
################################################

def writeCBSLine(index, chrom, pos, val, file):
	# 'Cause CBS has hard coded chromosome names... 
	if chrom[:3] == "chr":
		chrom = chrom[3:]
	file.write("%i\t%s\t%i\t%f\n" % (index, chrom, pos, val))

def writeCBSInput(chrom, pos, vals, output):
	file = open(output, "w")
	file.write("SNP\tChromosome\tPhysicalPosition\tsample1\n")
	for index in range(len(chrom)):
		writeCBSLine(index, chrom[index], pos[index], vals[index], file)
	file.close()

################################################
## Conversion functions for output
################################################
def parseCBSLine(line, chromosomeNames):
	items = line.strip().split()
	cnv = CNV(chromosomeNames[int(items[1]) - 1], items[2], items[3], [float(items[5])], "CNV:" + chromosomeNames[int(items[1]) - 1] + ":" + str((int(items[2]) + int(items[3])) / 2))
	cnv.numMarks = int(items[4])
	return cnv

def parseCBSFile(input, chromosomeNames):
	regions = []
	file2 = open(input)	
	file2.readline()
	for line in file2:
		regions.append(parseCBSLine(line, chromosomeNames))	
	file2.close()
	return regions

def uniq(list, elem):
	if len(list) == 0 or list[-1] != elem:
		list.append(elem)
	return list

################################################
## Master function 
################################################
def run(chrom, pos, vals):
	if not os.path.exists('CBS_OUT'):
		print 'Running CBS...'
		file1, input = tempfile.mkstemp(dir='.')

		writeCBSInput(chrom, pos, vals, input)
		if subprocess.Popen(['cbs.R', input, 'CBS_OUT'], stdout=sys.stdout, stderr=subprocess.STDOUT).wait() != 0:
		    sys.exit("CBS did not complete")
		os.remove(input)
	
	print 'Reading CBS output...'
	return parseCBSFile('CBS_OUT', chromosomeNames)

################################################
## Unit test
################################################

def main():
	chrom = ['1', '1', '1', '1']
	pos = [1000, 2000, 3000, 4000]
	vals = [1, 1, 1, 2]
	print "\n".join(map(str, run(chrom, pos, vals)))

if __name__=="__main__":
	main()
