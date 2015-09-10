#!/usr/bin/env python 

import sys
import argparse
import cnavgpost.genehistory.bedFileModule as bedmod
import cnavgpost.genehistory.segments_history_module as sgh 


BEDINDEX=0

def get_overlapping_beds(s, beds): 
	goodbeds=[]
	global BEDINDEX
	i=BEDINDEX
	bed=beds[i]
	while (not s.comes_after(bed)) and i>0:
		i=i-1
		bed=beds[i]
	for bed in beds[i:]:
#		sys.stderr.write("bed\t%sseg\t%s\t%d\t%d\nsegindex: %d, i: %d\n" % (bed, s.chr, s.start, s.end, BEDINDEX, i))
		if s.overlaps(bed):
			goodbeds.append(bed)
			i=i+1
		elif s.comes_after(bed):
			i=i+1
		elif s.comes_before(bed):
			BEDINDEX=min(i, len(beds)-1)
			return goodbeds
	BEDINDEX=min(i, len(beds)-1)
	return goodbeds

def label_segment(seg, beds): 
	mybeds=get_overlapping_beds(seg, beds)
	labels=[]
	for b in mybeds: 
		labels.append(b.name)
	return ",".join(sorted(set(labels)))

def read_in_bedfile(bedfn): 
	mybeds=[]
	for l in open(bedfn, 'r'): 
		mybeds.append(bedmod.BedEntry(l))
	return sorted(mybeds, key=lambda x: (x.chr, x.start, x.end))

def main(segfn, bedfn, useall, outfn): 
	mysegs=sgh.read_in_segs(segfn)
	for s in mysegs:
		s.chr = s.chr.replace("chr", "")
	mysegs=sorted(mysegs, key=lambda x: (x.chr, x.start, x.end))
	mybeds=read_in_bedfile(bedfn)
	bestsegs=[]
	outfh=open(outfn, 'w')
	for s in mysegs: 
		label=label_segment(s, mybeds)
		if not useall:  
			CNprofs=s.CNprofiles
			bestcnp=sorted(CNprofs, key=lambda x: x.likelihood, reverse=True)[0]
			s.CNprofiles=[bestcnp]
		outstr=""
		for line in str(s).split("\n"):
			if line !="": 
				outstr+=line.strip() + "\t"+label+"\n"
		outfh.write(outstr)
	outfh.close()
	

if __name__=="__main__": 
	parser=argparse.ArgumentParser(description='filter and label segment histories')
	parser.add_argument('--useall', default=False, action='store_true', help='if this option is set, all CNprofiles for seghists will be labeled and printed, not just the best ones.')
	parser.add_argument('seghists', help='a file of seghists')
	parser.add_argument('bedfile', help='a bedfile of genes to put as labels.')
	parser.add_argument('outfn', help='The output file (seghists.labeled).')
	args=parser.parse_args()
	main(args.seghists, args.bedfile, args.useall, args.outfn)

