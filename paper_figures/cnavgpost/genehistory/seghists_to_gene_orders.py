#!/usr/bin/env python 
import sys
import argparse
import cnavgpost.genehistory.bedFileModule as bedmod
import cnavgpost.genehistory.segments_history_module as sgh 

SEGINDEX=0

def get_overlapping_seghists(bed, mysegs): 
	goodsegs=[]
	global SEGINDEX
	i=SEGINDEX
	s=mysegs[i]
	# find my starting point to look for overlaps 
	while (not s.comes_before(bed)) and i>0: 
		i=i-1
		s=mysegs[i]
	for s in mysegs[i:]: 
#		sys.stderr.write("bed\t%sseg\t%s\t%d\t%d\nsegindex: %d, i: %d\n" % (bed, s.chr, s.start, s.end, SEGINDEX, i))
		s=mysegs[i] 
		if s.overlaps(bed): 
			goodsegs.append(s)
			i=i+1
		elif s.comes_before(bed):
			i=i+1
		elif s.comes_after(bed):
			SEGINDEX=i
			return goodsegs
	SEGINDEX=i
	return goodsegs
		
def get_best_overlapping_seghist(bed, mysegs): 
	seghists=get_overlapping_seghists(bed, mysegs) 
#	sys.stderr.write("There are %d overlapping\n" % len(seghists))
	maxlsh=0
	maxlbp=0
	bestsh=0
	bestbp=0
	for sh in seghists: 
		isbreak=False
		if sh.start == sh.end: isbreak=True
		ltot=get_total_likelihood(sh.CNprofiles)
		if isbreak and ltot > maxlbp: 
			bestbp=sh
			maxlbp=ltot
		elif not isbreak and ltot> maxlsh: 
			bestsh=sh
			maxlsh=ltot
	if bestsh != 0: 
		return bestsh
	else: 
		return bestbp
	
def get_total_likelihood(CNprofs): 
	mytot=0
	for c in CNprofs: 
		mytot+=c.likelihood
	return mytot

def main(segfn, bedfn, outfn, pval=0.5): 
	mysegs=sgh.read_in_segs(segfn)
	# get rid of chr if necessary, 
	for s in mysegs: 
		s.chr = s.chr.replace("chr", "")
	mysegs=sorted(mysegs, key=lambda x: (x.chr, x.start, x.end))
	outfh=open(outfn, 'w')
	outfh.write("chr\tstart\tend\tname\tLtotal\tLmax\tpreval\tprevalsd\n")
	for l in open(bedfn, 'r'): 
		bed=bedmod.BedEntry(l)
		seghist=get_best_overlapping_seghist(bed, mysegs)
		if seghist != 0:
			CNprofiles=seghist.CNprofiles
			ltot=get_total_likelihood(CNprofiles)
		else: ltot=0
		if ltot > pval: 
			bestCN=sorted(CNprofiles, key=lambda x: x.likelihood, reverse=True)[0] 
			mypreval=bestCN.pvals[0]
			myprevalsd=bestCN.pvalsd[0]
			mylscore=bestCN.likelihood
		else: 
			(mypreval, myprevalsd, mylscore)=(0,0, 0)
		outfh.write("\t".join(map(str, [
			bed.chr, bed.start, bed.end, bed.name,
			ltot, mylscore, mypreval, myprevalsd]))+"\n")


if __name__== "__main__": 
	parser=argparse.ArgumentParser(description='give seghists and a bedfile of genes, it will score each gene by the time that it gets first hit.')
	parser.add_argument('seghists', help='a file of seghists')
	parser.add_argument('bedfile', help='a bedfile of genes to score.')
	parser.add_argument('outfn', help='The output file')
	args=parser.parse_args()
	main(args.seghists, args.bedfile, args.outfn)
