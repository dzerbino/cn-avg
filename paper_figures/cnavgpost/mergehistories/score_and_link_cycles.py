#!/usr/bin/env python

import sys, os
import cnavgpost.mergehistories.event_cycles_module as histseg 
from cnavgpost.mergehistories.history_node_links_module import * 
import cPickle as pickle

def add_event_link_options(parser): 
	parser.add_argument('--cnavg', help='the CN-AVG output directory for a sample')
	parser.add_argument('--events', help='the file to write events to')  
	parser.add_argument('--links', help='the file to write links to')
	parser.add_argument('--inpickle', help='read the evnts from a pickled file.') 
	parser.add_argument('--outpickle', help='pickle the evnts and write them to this file.') 
	parser.add_argument('--totalp', help='total probability of histories', type=float)
	parser.add_argument('--binwidth', help='the multiplier for each history id to distinguish independent simulations.', type=float)
	parser.add_argument('--historystats', help='The file to write the history stats to')

def score_and_link_cycles(args):
	if args.binwidth: 
		histseg.Global_BINWIDTH=args.binwidth
	if not args.historystats and args.cnavg: 	
		historyScores=histseg.combine_history_statsfiles(args.cnavg)
	elif not os.path.isfile(args.historystats) and args.cnavg: 
		historyScores=histseg.combine_history_statsfiles(args.cnavg)
		np.savetxt(args.historystats, historyScores, fmt="%d", delimiter='\t')
	elif os.path.isfile(args.historystats):
		historyScores=np.loadtxt(args.historystats, dtype=int)
	if not historyScores: 
		sys.exit("Need to use --historystats or --cnavg option")
	totalp=0
	if args.totalp: 
		totalp=args.totalp
	else: 
		totalp = histseg.compute_likelihood_histories(historyScores[:,0], historyScores)
	allevents=[]
	if args.cnavg: 
		sys.stderr.write("using cnavg dir: %s\n" % (args.cnavg))
		allevents=histseg.get_events_from_cnavgdir(args.cnavg, historyScores, totalp)
	elif args.inpickle and os.path.isfile(args.inpickle): 
		sys.stderr.write("using pickle file\n")
		allevents=pickle.load(open(args.inpickle, 'rb'))
	sys.stderr.write("there are %d events\n" % (len(allevents)))
	if args.outpickle: 
		for event in allevents: 
			event.trim()
		eventfile= open(args.outpickle, 'wb')
		pickle.dump(allevents, eventfile, pickle.HIGHEST_PROTOCOL)
	if args.events: 
		eventfile=open(args.events, 'w')
		for evnt in allevents: 
			eventfile.write("%s" % (str(evnt)))
	# link the events...
	if args.links: 
		if not allevents: 
			sys.exit("Need events to link!  use --inpickle or --cnavg or --events")
		if not totalp: 
			sys.exit("Need a --totalp or --cnavg or --historystats options")	
		eventlinks = link_events_by_order_within_histories(allevents)
		linkfile=open(args.links, 'w')
		for link in eventlinks: 
			link.likelihood=histseg.compute_likelihood_histories(link.histories, historyScores)
			linkfile.write("%s" % (str(link)))


if __name__ == '__main__': 
	import argparse
	parser = argparse.ArgumentParser(description='Will essentially reformat braney files to have one event (rearrangment cycle) per line, and will assign each event a likelihood score based on the history complexity values.  Will also create links between cycles that co-occur in the same histories.')
	add_event_link_options(parser) 
	args=parser.parse_args()
	score_and_link_cycles(args)

