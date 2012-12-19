#!/usr/bin/env python

import sys
import os
import vcf
import cnavg.avg.graph as avg
import argparse
import bambam
import segmentCNVs
import cbs.cbs

##########################################
## Options parser
##########################################
def arguments():
	parser = argparse.ArgumentParser(description='Read BamBam files and detect inconsistencies')
	parser.add_argument('--breaks', dest='breaks', type=file, help='a BB breaks file')
	parser.add_argument('--cnvs', dest='cnvs', type=str, nargs='+', help='BB calls')
	parser.add_argument('--lengths', dest='lengths', type=file, help='chromosome lengths')
	parser.add_argument('--output', '-o', dest='output', type=str, help='Output filename', default='discrepancies.txt')
	parser.add_argument('--dir', '-d', dest='dir', type=str, help='Output filename', default='.')
	return parser.parse_args()

##########################################
## Core function
##########################################
def discrepancies_CNV(cnv, nodeindex, breakends, file):
	if nodeindex > 0:
		nodeindex -= 1

	while nodeindex < len(breakends) and cnv > breakends[nodeindex]:
		nodeindex += 1

	while nodeindex < len(breakends) and cnv == breakends[nodeindex]:
		assert cnv.val >= 0
		breakends[nodeindex].segment = cnv.val
		if breakends[nodeindex].start <= cnv.softStart and breakends[nodeindex].finish >= cnv.start:
			if cnv.startCap is None and breakends[nodeindex].orientation:
				cnv.startCap = breakends[nodeindex]
			else:
				pass
		elif breakends[nodeindex].start <= cnv.finish and breakends[nodeindex].finish >= cnv.softFinish:
			if not breakends[nodeindex].orientation:
				cnv.finishCap = breakends[nodeindex]
				break
		nodeindex += 1

	# Have to assign aritificial breakends
	if cnv.startCap is None:
		file.write("\t".join(map(str, [cnv.chr, cnv.start, cnv.softStart])) + "\n")
	if cnv.finishCap is None:
		file.write("\t".join(map(str, [cnv.chr, cnv.softFinish, cnv.finish])) + "\n")

	return nodeindex

def discrepancies(lengths, breakends, cnvs, file):
	reduce(lambda L, E: discrepancies_CNV(E, L, breakends, file), cnvs, 0)

##########################################
## Wrapper
##########################################
def main():
	args = arguments()
	os.chdir(args.dir)

	# Lengths
	lengths = bambam.parseChromLengths(args.lengths)

	# Breakends
	breakends = bambam.parseBreaksFile(args.breaks)
	breakends = bambam.removeNonChromosomalBreakpoints(breakends, lengths)

	# CNV data
	cnvs = segmentCNVs.segmentCNVs(bambam.parseCNVData(args.cnvs))

	# Business logic
	discrepancies(lengths, breakends, cnvs, open(args.output, 'w'))
	if not os.path.exists('tracks'):
		os.mkdir('tracks')
	os.rename(args.output, os.path.join('tracks', args.output))

if __name__ == "__main__":
	main()
