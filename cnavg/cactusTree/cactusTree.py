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
Integrated jobTree pipeline.
"""

import sys
import os
import glob
import subprocess
import shutil
import tempfile
import cPickle as pickle
import gzip
import tempfile

from optparse import OptionParser

from jobTree.src.bioio import logger
from jobTree.src.bioio import getBasicOptionParser
from jobTree.src.bioio import getTempFile
from jobTree.src.bioio import getLogLevelString
from jobTree.src.bioio import system

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

IDEAL_JOB_RUNTIME = 500
MEM4G=4290000000
MEM2G=3000000000

class MyParser(OptionParser):
	""" Command line options """
	def __init__(self):
		OptionParser.__init__(self)
		self.add_option('--vcf', '-v', dest='vcffile', help='A VCF (ver >= 4.1) file')
		self.add_option('--bambam', '-b', dest='bambam', help='Common BamBam file string')
		self.add_option('--snpsfile', '-p', dest='snpsfiles', help='Common BamBam file string')
		self.add_option('--breaks', '-k', dest='breaks', help='A BamBam breaks file')
		self.add_option('--bam', '-B', dest='bam', help='Bam file')
		self.add_option('--lengths', '-l', dest='chromLengths', help='Chromosome lengths')
		self.add_option('--dir', '-d', dest='dir', help='Working directory')
		self.add_option('--size', '-s', dest='size', type=int, help='Histories per job', default=300)
		self.add_option('--number', '-n', dest='number', type=int, help='Number of jobs', default=500)
		self.add_option('--debug', '-g', dest='debug', action='store_true', help='Debug switch for whatever')
		self.add_option('--url', '-u', dest='url', help='URL prefix of home directory', default="https://tcga1.kilokluster.ucsc.edu/~dzerbino/")
		Stack.addJobTreeOptions(self)

def do(command, outf=None, inf=None, append=False):
	""" Convience wrapper around subprocess.call """
	if outf is not None:
		if not append:
			output = open(outf,"w+")
		else:
			output = open(outf,"a+")
	else:
		output = tempfile.TemporaryFile()
	error = tempfile.TemporaryFile()
	if inf is not None:
		input = open(inf)
	else:
		input = None
    	ret = subprocess.call(command, stdout=output, stderr=error, stdin=input, shell=True)
	if ret != 0:
		print inf
		print command
		print outf
		output.seek(0)
		for line in output:
			print line,
		error.seek(0)
		for line in error:
			print line,
		error.close()
		output.close()
		sys.exit(ret)
	error.close()
	output.close()
	if input is not None:
		input.close()

def buildGraph(options):
	""" Builds initial Cactus graph """
    	command = "cn-avg.py --dir " + options.dir
	command += " --bambam " + " ".join(options.bambam)
	if options.snpsfiles is not None:
		command += " --snps" + options.snpsfiles
	command += " --breaks " + options.breaks
	command += " --lengths " + options.chromLengths
	do(command)

class SetupPhase(Target):
    """ Beginning of pipeline, builds cactus, launches parallel samplers """
    def __init__(self, options):
        Target.__init__(self, time=10, memory = MEM2G)
	self.options = options
        
    def run(self):
        buildGraph(self.options)
	for index in range(self.options.number):
	    self.addChildTarget(SampleCycles(self.options, index))
        self.setFollowOnTarget(Reduce(self.options))

class SampleCycles(Target):
    """ Sampler: normalizes the cactus graph randomly, then performs MCMC search through possible histories. """
    def __init__(self, options, index):
    	Target.__init__(self, time=600, memory=MEM4G)
	self.options = options
	self.index = str(index)

    def run(self):
    	do("cn-avg.py -d " + self.options.dir + " -i " + self.index + " -s " + str(self.options.size))

class Reduce(Target):
    """ Collects sampled histories and selects optimal ones for a separate file """
    def __init__(self, options):
    	Target.__init__(self, time=600, memory=MEM2G)
	self.options = options

    def run(self):
        if self.options.size > 0 and self.options.number > 0:
		temp = tempfile.mkstemp()[1]
		do("gzip -dc %s/HISTORIES_*.braney | awk '/^A/ {if ($14 < min || ! min) min = $14} END {print min}'" % (self.options.dir), outf=temp)
		score = open(temp).readline().strip()
		files = glob.glob("%s/HISTORIES_*.braney" % (self.options.dir))
		out = open("%s/OPTIMAL_HISTORIES.braney" % (self.options.dir), "w")
		counter = -1
		lastID = -1
		lastFile = "" 
		for file in files:
			input = gzip.open(file)
			for line in input:
				items = line.strip().split()
				if items[0] == "A" and items[13] == score:
					if items[8] != lastID or file != lastFile:
						counter += 1
					lastID = items[8]
					lastFile = file
					items[8] = str(counter)
					out.write("\t".join(items) + "\n")
				elif items[0] != "A" and items[9] == score:
					if items[4] != lastID or file != lastFile:
						counter += 1
					lastID = items[4]
					lastFile = file
					items[4] = str(counter)
					out.write("\t".join(items) + "\n")
			input.close()
		out.close()

        self.setFollowOnTarget(Process(self.options))

class Process(Target):
    """ Prepares different tracks for visual inspection on the browser """
    def __init__(self, options):
    	Target.__init__(self, time=10, memory=MEM2G)
	self.options = options

    def run(self):
    	self.options.tracks = os.path.join(self.options.dir, "tracks")
        if not os.path.exists(self.options.tracks):
		os.mkdir(self.options.tracks)
        self.addChildTarget(PrepareBAMBAMGraph(self.options))
        self.addChildTarget(PrepareCBSGraph(self.options))
	self.addChildTarget(PrepareBreakendGraph(self.options))
        self.addChildTarget(PrepareCycleGraph(self.options))
	self.addChildTarget(PrepareCactusCoverageGraphs(self.options))

	os.chdir(self.options.tracks)
	logger.debug(os.path.basename(self.options.bam))
	if os.path.exists(os.path.basename(self.options.bam) + ".bai"):
		os.remove(os.path.basename(self.options.bam) + ".bai")
	if os.path.exists(os.path.basename(self.options.bam)):
		os.remove(os.path.basename(self.options.bam))
	do("ln -s %s" % self.options.bam)
	do("ln -s %s.bai" % self.options.bam)
        self.setFollowOnTarget(PrepareIndex(self.options))

chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

class PrepareBAMBAMGraph(Target):
    """ Prepares wiggle file of Bambam coverage bins """
    def __init__(self, options):
    	Target.__init__(self, time=1800, memory=MEM2G)
	self.options = options

    def run(self):
    	os.chdir(self.options.tracks)
	out = open("bambam.cov.bg", "w")
	lastChr = None
	lastPos = None
	for file in map(open, self.options.bambam):
		for line in file:
			items = line.strip().split()
			if items[0] == "CNV":
				if items[1][:3] != "chr":
					items[1] = "chr" + items[1]
				if lastChr is not None and lastChr == items[1] and lastPos >= int(items[2]):
					items[2] = str(lastPos + 1)
				out.write("\t".join(map(str,items[1:5])) + "\n")
				lastChr = items[1]
				lastPos = int(items[3])
		file.close()
	out.close()
	do("bedGraphToBigWig bambam.cov.bg " + self.options.chromLengths + " bambam.cov.bw")

class PrepareCBSGraph(Target):
    """ Prepares wiggle file of CBS coverage segments """
    def __init__(self, options):
    	Target.__init__(self, time=1800, memory=MEM2G)
	self.options = options

    def run(self):
    	os.chdir(self.options.tracks)
    	file = open(os.path.join(self.options.dir, "CBS_OUT"))
	out = open("cbs.cov.bg", "w")
    	for line in file:
		items = line.strip().split()
		if items[2] != "loc.start":
			if items[1] == "23":
				items[1] = "chrX"
			elif items[1] == "24":
				items[1] = "chrY"
			else:
				items[1] = "chr" + items[1]
			out.write("\t".join(map(str, items[1:4] + [items[5]])) + "\n")
	file.close()
	out.close()
	do("bedGraphToBigWig cbs.cov.bg " + self.options.chromLengths + " cbs.cov.bw")

class PrepareCactusCoverageGraphs(Target):
    """ Prepares wiggle file of normalized and balanced flow change on cactus """
    def __init__(self, options):
    	Target.__init__(self, time=1800, memory=MEM2G)
	self.options = options

    def run(self):
    	os.chdir(self.options.tracks)
    	file = open(os.path.join(self.options.dir, "CACTUS_0"))
	C = pickle.load(file)
	file = open("cactus.cnv.bg", "w")
	file.write(C.printCNVs())
	file.close()
	do("bedGraphToBigWig cactus.cnv.bg " + self.options.chromLengths + " cactus.cnv.bw")
	file = open("majority.cnv.bg", "w")
	file.write(C.printMajority())
	file.close()
	do("bedGraphToBigWig majority.cnv.bg " + self.options.chromLengths + " majority.cnv.bw")
	file = open("minority.cnv.bg", "w")
	file.write(C.printMinority())
	file.close()
	do("bedGraphToBigWig minority.cnv.bg " + self.options.chromLengths + " minority.cnv.bw")

class PrepareBreakendGraph(Target):
    """ Prepares adjacency file of Bambam breakend calls """ 
    def __init__(self, options):
    	Target.__init__(self, time=100, memory=MEM2G)
	self.options = options

    def run(self):
    	os.chdir(self.options.tracks)
    	file = open(self.options.breaks)
	out = open("breaks.braney", "w")
    	for line in file:
		items = line.strip().split()

		if items[3] == "-":
			sgnA = "+"
		else:
			sgnA = "-"

		if items[7] == "-":
			sgnB = "+"
		else:
			sgnB = "-"

		out.write("\t".join(map(str,["A", items[0], items[1], sgnA, items[4], items[5], sgnB, 1000, '0', '0','0','0','0','0'])) + "\n")
	file.close()
	out.close()
    	do("braneyConversion.py %s 0" % self.options.chromLengths, inf="breaks.braney", outf="breaks.adj")

class PrepareCycleGraph(Target):
    """ Prepares adjacency file of an optimal history """
    def __init__(self, options):
    	Target.__init__(self, time=100, memory=MEM2G)
	self.options = options

    def run(self):
    	os.chdir(self.options.tracks)
    	do("braneyConversion.py %s %i" % (self.options.chromLengths, 0), inf=os.path.join(self.options.dir, "OPTIMAL_HISTORIES.braney"), outf="cnavg.adj")

class PrepareIndex(Target):
    """ Prepares text file with UCSC browser track options """
    def __init__(self, options):
    	Target.__init__(self, time=10, memory=MEM2G)
	self.options = options

    def run(self):
    	home = self.options.url + "/" + os.path.basename(self.options.dir) + "/"
    	os.chdir(self.options.tracks)
	file = open("index.txt", "w")
	file.write("track name=cnavg type=adjacency visibility=hide\n")
	file.write(home + "cnavg.adj\n")
	file.write("track name=bambamadj type=adjacency visibility=hide\n")
	file.write(home + "breaks.adj\n")
	file.write("track name=minority.cnv type=bigWig alwaysZero=on visibility=hide bigDataUrl=" + home + "minority.cnv.bw\n")
	file.write("track name=majority.cnv type=bigWig alwaysZero=on visibility=hide bigDataUrl=" + home + "majority.cnv.bw\n")
	file.write("track name=cactus.cnv type=bigWig alwaysZero=on visibility=hide bigDataUrl=" + home + "cactus.cnv.bw\n")
	file.write("track name=cbs.cov type=bigWig alwaysZero=on yLineMark=1 yLineOnOff=on visibility=hide bigDataUrl=" + home + "cbs.cov.bw\n")
	file.write("track name=bambam.cov type=bigWig alwaysZero=on yLineMark=1 yLineOnOff=on visibility=hide bigDataUrl=" + home + "bambam.cov.bw\n")
	file.write("track name=bam type=bam visibility=hide pairEndsByName=true bigDataUrl=" + home + os.path.basename(self.options.bam) + "\n")
	file.close()

############################################################
# Main script
############################################################
def main():
    """
    Runs the integrated pipeline.
    """
    options, args = MyParser().parse_args()

    if options.bam is None:
	    assert False, "No Bam file!"
    else:
	    options.bam = os.path.abspath(options.bam)

    if options.dir is None:
	    assert False, "No working directory!"
    else:
	    options.dir = os.path.abspath(options.dir)

    if options.jobTree is None:
	    options.jobTree = options.dir + "jobTree"

    if options.chromLengths is not None:
	    options.chromLengths = os.path.abspath(options.chromLengths)

    if options.bambam is not None:
	    options.bambam = map(os.path.abspath, glob.glob(options.bambam))
    
    if options.breaks is not None:
	    options.breaks = os.path.abspath(options.breaks)
    
    ##########################################
    # Off we go...
    ##########################################
    i = Stack(SetupPhase(options)).startJobTree(options)
    if i:
	raise RuntimeError("The jobTree contained %i failed jobs" % i)

if __name__ == '__main__':
    from cnavg.cactusTree.cactusTree import *
    main()
