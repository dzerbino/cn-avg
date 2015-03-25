#!/usr/bin/env python
# Copyright (c) 2012, Daniel Zerbino
# All rights reserved.
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

import sys
import os
import os.path
import argparse
import cPickle as pickle
import random
import glob
import gzip 

import cnavg.preprocess.vcf as vcf
import cnavg.preprocess.bambam as bambam

import cnavg.avg.balanced as balancedAVG
import cnavg.cactus.graph as cactus
import cnavg.cactusSampling.sampling as normalized
import cnavg.cactus.oriented as oriented
import cnavg.cactus.balanced as balancedCactus
import cnavg.historySampling.cycleCover as cycleCover
import cnavg.historySampling.sampleGraphCycles as sampleGraphCycles
import cnavg.history.flattened as flattened
import cnavg.history.ordered as ordered
import cnavg.history.debug as debug
from cnavg.history.ordered import prettify

def _parseOptions():
	print "Parsing options"

	parser = argparse.ArgumentParser(description="Process a VCF files to sample possible historical explanations")
	parser.add_argument('--vcf', '-v', dest='vcffile', type=file, help='A VCF (ver >= 4.1) file')
	parser.add_argument('--bambam', '-b', dest='bambam', nargs='*', help='BamBam files')
	parser.add_argument('--snps', '-p', dest='snpsfiles', nargs='*', help='SNPs files (optional)')
	parser.add_argument('--index', '-i', dest='index', type=int, help='ID of sampling run')
	parser.add_argument('--breaks', '-k', dest='breaks', type=file, help='A BamBam breaks file')
	parser.add_argument('--lengths', '-l', dest='chromLengths', type=file, help='Chromosome lengths')
	parser.add_argument('--dir', '-d', dest='dir', help='Working directory')
	parser.add_argument('--debug', '-g', dest='debug', action='store_true', help='Debug switch for whatever')
	parser.add_argument('--continue', '-c', dest='cont', action='store_true', help='Continue sampling for 24 hours')
	parser.add_argument('--integer', '-n', dest='integer', action='store_true', help='Integer switch for idealized integer histories')
	parser.add_argument('--size', '-s', dest='size', type=int, default=100, help='Number of sampled histories')
	parser.add_argument('--temp', '-t', dest='temp', type=float, default=1, help='Starting temperature of MCMC sampling')
	parser.add_argument('--simulation', dest='simulation', action='store_true', help='Simuated histories')
	return parser.parse_args()

def _parseGraph(options):
	print "Parsing input files"

	if options.bambam is not None and options.breaks is not None and options.chromLengths is not None:
		options.bambam = sum(map(glob.glob, options.bambam), [])
		assert len(options.bambam) > 0, options.bambam
		breakends = bambam.parse(options.bambam, options.breaks, options.chromLengths, options.snpsfiles)
	elif options.vcffile is not None and options.chromLengths is not None:
		breakends = vcf.parse(options.vcffile, options.chromLengths)
	else:
		if options.vcffile is None:
			print "No VCF"
		else:
			print "VCF: %s" % options.vcffile

		if options.chromLengths is None:
			print "No chromosome lengths"
		else:
			print "Chromosome lengths: %s" % options.chromLengths

		if options.bambam is None:
			print "No BamBam files"
		else:
			print "BamBam files: %s" % options.bambam

		if options.breaks is None:
			print "No BamBam break file"
		else:
			print "Breaks lengths: %s" % options.breaks

		sys.exit("Not enough files")

	breakends.validate()
	return breakends.avg()

def main():
	options = _parseOptions()

	sampleGraphCycles.TEMPERATURE = options.temp

	if options.dir is not None:
		if not os.path.exists(options.dir):
			os.mkdir(options.dir)
		os.chdir(options.dir)

	if options.simulation:
		debug.RATIO_TO_OFFSET = False

	if options.index is None:
		## Initial graph construction
		G = _parseGraph(options)
		B = balancedAVG.BalancedAVG(G)
		C = cactus.Cactus(B)
		pickle.dump(C, open('CACTUS', "wb"))
	else:
		H = None

		if options.debug:
			## Picking up from where we started
			OC = pickle.load(open('CACTUS_%i' % options.index))
			random.setstate(pickle.load(open("STATE_%i" % options.index)))
		elif options.cont:
			## Picking up from where we stopped
			OC = pickle.load(open('CACTUS_%i' % options.index))
			## Going through all the histories to the last in file 
			file= open('HISTORIES_%i' % (options.index))
			while True:
				try:
					H = pickle.load(file)
				except:
					break
			file.close()
		else:
			## Just moving from there 
			pickle.dump(random.getstate(), open("STATE_%i" % options.index, "wb"))
			C = pickle.load(open('CACTUS'))

			## Sampling possible cactus
			NC = normalized.NormalizedCactus(C)
			if debug.RATIO_TO_OFFSET:
				BC = balancedCactus.BalancedCactus(NC)
			else:
				BC = NC
			OC = oriented.OrientedCactus(BC)

			## Saving sampled cactus
			pickle.dump(OC, open('CACTUS_%i' % options.index, "wb"))
		

		# Moving into historical space
		if options.integer:
			debug.INTEGER_HISTORY = True
		if H is None:
			H = cycleCover.initialHistory(OC)
		FH = flattened.flattenGraph(H)
		S = FH.simplifyStubsAndTrivials()
		F = S.removeLowRatioEvents(debug.RATIO_CUTOFF)
		O = ordered.OrderedHistory(F)

		# Preparing file for progressive write
		if options.cont:
			stats_file = open("HISTORY_STATS_%li" % options.index, "a")
			pickle_file = open('HISTORIES_%i' % options.index, "ab")
			braney_file = gzip.open("HISTORIES_%i.braney" % options.index, "a")
		else:
			stats_file = open("HISTORY_STATS_%li" % options.index, "w")
			pickle_file = open('HISTORIES_%i' % options.index, "wb")
			braney_file = gzip.open("HISTORIES_%i.braney" % options.index, "w")
		stats_file.write("%s\n" % H.stats())
		#pickle.dump(H, pickle_file)
		braney_file.write("%s\n" % O.braneyText(0, H.rearrangementCost()))
		#tree_file = open("HISTORY_TREES_%li" % options.index, "w")
		#tree_file.write("%s\n" % O.newick())
		tree_file = None

		# Sampling
		SH = sampleGraphCycles.sample(H, options.size, pickle_file, stats_file, braney_file, tree_file)

		# Cleaning up
		stats_file.close()
		pickle_file.close()
		braney_file.close()
		#tree_file.close()

		## Removing temp file
		if os.path.exists("STATE_%i" % options.index):
			os.remove("STATE_%i" % options.index)

if __name__ == "__main__":
	main()
