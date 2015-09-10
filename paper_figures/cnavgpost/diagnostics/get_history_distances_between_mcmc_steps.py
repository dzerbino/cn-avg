#!/usr/bin/env python 

import sys, os 
import pickle 
import cnavgpost.mergehistories.event_cycles_module as histseg

def get_historyid_vs_not_eventcounts(refhistoryid, events, outfh): 
	refevents=0
	inboth=0
	otherevents=0
	for evnt in events: 
		if refhistoryid in evnt.histories: 
			refevents+=1
			if len(evnt.histories) >1: 
				inboth+=1
		else: 
			otherevents+=1
	outfh.write("%d\tNA\t%d\t%d\t%d\n" % (refhistoryid, refevents, otherevents, inboth)) 
	
def get_eventcounts_for_history(refhistoryid, events): 
	refevents=0
	for evnt in events: 
		if refhistoryid in evnt.histories: 
			refevents+=1
	return refevents

def get_history_distances_between_mcmc_steps(events, refhistoryid, refhistoryid2, numsteps, stepsize1, stepsize2, outfh): 
	outfh.write("historyA\thistoryB\tnum_eventsA\tnum_eventsB\tnum_both\n")
	if numsteps==0:
		get_historyid_vs_not_eventcounts(refhistoryid, events, outfh) 
	elif numsteps>0:
		for i in xrange(numsteps):
			otherevents=0
			inboth=0
			if refhistoryid2: 
				otherid = refhistoryid2 + i*(stepsize2)
				historyid = refhistoryid + i*(stepsize1)
				refevents=0
			else: 
				otherid=refhistoryid + i*(stepsize2)
				historyid=refhistoryid
				refevents=get_eventcounts_for_history(historyid, events)
			for evnt in events: 
				if otherid in evnt.histories: 
					otherevents+=1
					if historyid in evnt.histories: 
						inboth+=1
				if refhistoryid2: 
					if historyid in evnt.histories: 
						refevents+=1
			outfh.write("%d\t%d\t%d\t%d\t%d\n" % (historyid, otherid, refevents, otherevents, inboth)) 
			
if __name__ == '__main__': 
	import argparse
	parser=argparse.ArgumentParser(description='compute the number of events that are different across a set of histories')

	parser.add_argument('--braneyfile', help='a HISTORIES*.braney file')
	parser.add_argument('--pevntsfile', help='a *.pevnts file')
	parser.add_argument('--refhistoryid', help='the historyid that you want to take as step 0.', default=2500, type=int)
	parser.add_argument('--refhistoryid2', help='a second historyid that you want to take as step 0.', type=int)
	parser.add_argument('--numsteps', help='the number of steps to take.', default=0, type=int)
	parser.add_argument('--stepsize', help='the step size from the reference id(s) in the mcmc that you want to go.', default=1, type=int)
	parser.add_argument('--stepsize2', help='the step size from the second history id(s) in the mcmc that you want to go.', default=1, type=int)
	args=parser.parse_args()
	if args.braneyfile and not args.pevntsfile: 
		events=histseg.make_events_from_braneyfn(args.braneyfile)
	elif args.pevntsfile: 
		events=pickle.load(open(args.pevntsfile, 'rb'))
		#for e in events: 
		#	e.histories=histseg.listout_ranges(e.histRanges)
	get_history_distances_between_mcmc_steps(events, args.refhistoryid, args.refhistoryid2, args.numsteps, args.stepsize, args.stepsize2, sys.stdout) 
