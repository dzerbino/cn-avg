#!/usr/bin/env python
import sys, os
import cnavgpost.mergehistories.event_cycles_module as histseg
import pysam 
import cPickle as pickle
import re

def merge_segments_into_regions(event): 
	if len(event.segs) ==0: 
		event.make_segs_from_str()
	mylocs=[]
	for seg in event.segs: 
		if seg.adj: 
			mylocs.append((seg.chr, seg.start, seg.start+1))
			mylocs.append((seg.chr2, seg.end, seg.end+1))
		else: 
	#	elif event.cnval != 1: # only count a region as affected if it has a copy number change, not if it is simply inverted.  
			mylocs.append((seg.chr, seg.start, seg.end))
	mycoords=[]
	while len(mylocs) >0:
		(chr, start, end) = mylocs.pop(0)
		merger=True
		while merger:
			merger=False
			newlocs=[]
			for j in xrange(len(mylocs)): 
				(chr2, start2, end2)= mylocs[j]
				if (chr == chr2 and start <=end2 and end >= start2): 
					start=min(start, start2)	
					end=max(end, end2)	
					merger=True
				else:
					newlocs.append([chr2, start2, end2])
			mylocs=newlocs
		mycoords.append([chr, start, end])
	return mycoords

def consolidate_genes(genelist): 
	mynames=set([])
	for gene in genelist:
		(chr, start, end, name) = gene.split("\t")
		if name not in mynames: 
			mynames.add(name)
	return sorted(list(mynames))

def annotate_events(allevents, mytabix): 
	myannotations=[] # this will be a list of tuples: (eventid, comma-sep list of gene names)
	for myevent in allevents: 
		coords= merge_segments_into_regions(myevent)	
		mygenes=[]
		for coord in coords: 
			if coord[1] == -1: coord[1]=0
			if coord[0] != "None": 
				region="%s:%d-%d" % (coord[0], coord[1], coord[2])
				try: 
					mygenes+=mytabix.fetch(region)  # make sure the chromosome format is the same as what is used in the tabix file 
				except ValueError:
					sys.stderr.write("Error in tabix fetch for %s, with %s, %s\n" % (args.tabixfile, region, str(ValueError))) 
		genesetstr="None"
		if len(mygenes) >0: 
			geneset=consolidate_genes(mygenes)
			genesetstr=",".join(map(str, geneset))
		myannotations.append((myevent.id, genesetstr))
	return myannotations

def main(allevents, mytabix, outputfh):
	annotations=annotate_events(allevents, mytabix)
	for myann in annotations: 
		(id, genesetstr)= myann
		outputfh.write(id + "\t" + genesetstr + "\n")

if __name__ == '__main__': 
	import argparse
	parser=argparse.ArgumentParser(description='Given a tabixfile with gene annotations and an .evnt file with rearrangment cycles, it will intersect the genes with the cycles to create groups of genes that are linked in rearrangements.')
	parser.add_argument('pevnt', help='the pickled event file (.pevnt)') 
	parser.add_argument('tabixfile', help='a tabix file of genes')
	args=parser.parse_args()
	allevents=pickle.load(open(args.pevnt, 'rb'))
	mytabix=pysam.Tabixfile(args.tabixfile, 'r')
	main(allevents, mytabix, sys.stdout)

