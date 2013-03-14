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
	parser.add_argument('--size', '-s', dest='size', type=int, default=100, help='Number of sampled histories')
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
	if options.dir is not None:
		if not os.path.exists(options.dir):
			os.mkdir(options.dir)
		os.chdir(options.dir)

	if options.index is None:
		## Initial graph construction
		G = _parseGraph(options)
		B = balancedAVG.BalancedAVG(G)
		C = cactus.Cactus(B)
		pickle.dump(C, open('CACTUS', "wb"))
	else:
		if not options.debug:
			## Just moving from there 
			pickle.dump(random.getstate(), open("STATE_%i" % options.index, "wb"))
			C = pickle.load(open('CACTUS'))

			## Sampling possible cactus
			NC = normalized.NormalizedCactus(C)
			BC = balancedCactus.BalancedCactus(NC)
			OC = oriented.OrientedCactus(BC)

			## Saving sampled cactus
			pickle.dump(OC, open('CACTUS_%i' % options.index, "wb"))
		else:
			## Picking up from where we started
			print 'DEBUG'
			OC = pickle.load(open('CACTUS_%i' % options.index))
			random.setstate(pickle.load(open("STATE_%i" % options.index)))

		# Moving into historical space
		H = cycleCover.initialHistory(OC)

		# Preparing file for progressive write
		stats_file = open("HISTORY_STATS_%li" % options.index, "w")
		stats_file.write("%s\n" % H.stats())
		pickle_file = open('HISTORIES_%i' % options.index, "wb")
		#pickle.dump(H, pickle_file)
		braney_file = gzip.open("HISTORIES_%i.braney" % options.index, "w")
		braney_file.write("%s\n" % prettify(H))
		tree_file = open("HISTORY_TREES_%li" % options.index, "w")
		tree_file.write("%s\n" % H.newick())

		# Sampling
		SH = sampleGraphCycles.sample(H, options.size, pickle_file, stats_file, braney_file, tree_file)

		# Cleaning up
		stats_file.close()
		pickle_file.close()
		braney_file.close()
		tree_file.close()

		## Removing temp file
		if os.path.exists("STATE_%i" % options.index):
			os.remove("STATE_%i" % options.index)

if __name__ == "__main__":
	main()
