#!/usr/bin/env python 
#This is used to look at individual edges (segments or adjacencies) across histories as opposed to the cyclic events.   

import sys, os
import copy
import cPickle as pickle 
import cnavgpost.mergehistories.event_cycles_module as histseg
import argparse
import numpy as np
import re
bindir="/inside/home/tballing/cnavg-study/cnavgmbin"

def unique_loc_edges(alledges): 
	for e in alledges: 
		#e.segstr=histseg.remove_signs_from_segstr(e.segstr)[0]
		e.cnval=0 #change this to 0 because it will be uninformative in this context
	sortededges=sorted(alledges, key=lambda x: (x.segstr)) 
	unique_edges=[]
	prevedge=sortededges[0]
	for edge in sortededges[1:]:
		if prevedge.segstr==edge.segstr: 
			prevedge.merge_Event_data(edge) #modifies prevedge
		else: 
			unique_edges.append(prevedge)
			prevedge=edge
	unique_edges.append(prevedge)
	return unique_edges
	
def score_edges_within_pevents(allevents, historyScores, totalp, prev_error=0.05, ignore_cn=True): 
	prevalence_error=prev_error
	sys.stderr.write("number of events: %d\n" % (len(allevents)))
	sys.stderr.write("ignore_cn: %s\n" % ignore_cn) 
	alledges=[]
	for event in allevents:
		event.unpack() 
		for seg in event.segs: 
			edge=copy.deepcopy(event)
			edge.segs=[seg]
			edge.make_segstr()
			alledges.append(edge)
	if ignore_cn:
		unique_edges=unique_loc_edges(alledges)
	else: 
		unique_edges=histseg.unique_c_events_list(alledges)
	sys.stderr.write("totalp: %s\n" % (str(totalp)))
	for edge in unique_edges: 
		edge.update(historyScores)
		edge.likelihood=histseg.compute_likelihood_histories(edge.histories, historyScores, totalp)
		edge.trim()
	return unique_edges 

def add_score_edges_options(parser): 
	parser.add_argument('pevnts', help='a .pevnts file.')
	parser.add_argument('historystats', help='The file with historystats')
	parser.add_argument('--outpickle', help='pickle the edges and write them to this file.')
	parser.add_argument('--edges', help='write edges to this file in text (not pickled).')
	parser.add_argument('--prevalence_error', help='the difference in prevalences to be considered the same.', type=float, default=0.05)
	parser.add_argument('--ignore_cn', help='merge together edges with different CN values.', default=False, action='store_true')
	parser.add_argument('--totalp', help='total probability of the histories', type=float)
	parser.add_argument('--binwidth', help='the multiplier between history ids of independent runs', default=histseg.Global_BINWIDTH, type=int)	

if __name__ == "__main__": 
	parser = argparse.ArgumentParser(description='Given an .pevnts file, it will split events into edges, combine equivalent edges (segments or adjacencies), and score them by likelihood.')
	add_score_edges_options(parser)
	args=parser.parse_args()
	histseg.Global_BINWIDTH=args.binwidth
	allevents=pickle.load(open(args.pevnts, 'rb'))
	historyScores=np.loadtxt(args.historystats, dtype=int)
	totalp=0
	if args.totalp: 
		totalp=args.totalp
	else: 
		totalp = histseg.compute_likelihood_histories(historyScores[:,0], historyScores) 
	alledges = score_edges_within_pevents(allevents, historyScores, totalp, args.prevalence_error, args.ignore_cn)
	if args.edges:
		outfile=open(args.edges, 'w') 
		for edge in alledges: 
			outfile.write(str(edge))
	if args.outpickle: 
		pickle.dump(alledges, open(args.outpickle, 'wb'), pickle.HIGHEST_PROTOCOL)


