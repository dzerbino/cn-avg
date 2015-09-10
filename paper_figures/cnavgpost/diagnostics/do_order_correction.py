#!/usr/bin/env python 

import argparse 
import cPickle as pickle
import sys, os 
import cnavgpost.mergehistories.event_cycles_module as ecycles 
import numpy as np
import re


def count_earlylate_with_correction(events, historyScores, outfn1, outfn2): 
	numhists=historyScores.shape[0]
	mymax=len(events)+1
	myranks=np.ones((numhists, len(events))) * mymax 
	myTPranks=np.ones((numhists, len(events))) * mymax 
#	sys.stderr.write("myranks: %s, myTPranks %s\n" % (str(myranks.shape), str(myTPranks.shape)))
	simhist=0
	for j in xrange(len(events)): 
		e=events[j]
		#e.histories=ecycles.listout_ranges(e.histRanges)
		hindices = ecycles.historyids_to_indices(e.histories, historyScores)
		for h in xrange(len(e.histories)): 
			i=hindices[h]
			myord=float(e.orders[h])
			myranks[i,j]=myord
			if simhist in e.histories: 
				myTPranks[i,j]=myord 
	# change the orders into ranks.  The ranks will be different if we just look at the TP events in a history vs if we include all of the events. 
	hlengths=np.sum(myranks<mymax, axis=1).astype('float')
	trueonlylen=np.sum(myTPranks<mymax, axis=1).astype('float')
	np.savetxt(outfn2, np.vstack((hlengths, trueonlylen)).T, fmt="%d", delimiter='\t', header="length\tlen_onlytrue")
	#process the ranks of the events including the whole history
	xord=myranks[hlengths>0,:].argsort()
	xranks=xord.argsort()
	cranks=xranks.astype('float')
	hlengths=hlengths[hlengths>0]
	for i in xrange(hlengths.shape[0]): 
		cranks[i,:]=cranks[i,:]/(hlengths[i]-1)
	#process the ranks for the histories including only the TP events
	xord=myTPranks[trueonlylen>0,:].argsort()
	xranks=xord.argsort()
	cTPranks=xranks.astype('float')	
	trueonlylen=trueonlylen[trueonlylen>0]
	for i in xrange(trueonlylen.shape[0]): 
		cTPranks[i,:]=cTPranks[i,:]/(trueonlylen[i]-1)
	# skip history 0 because that's the simulated history
	truth=cranks[0,:]
	cranks=cranks[1:,:]
	cTPranks=cTPranks[1:,:]
	earlycnts=np.sum(cranks<0.5, axis=0)
	latecnts=np.sum(np.logical_and(cranks>0.5, cranks<=1), axis=0)
	tpearlycnts=np.sum(cTPranks<0.5, axis=0)
	tplatecnts=np.sum(np.logical_and(cTPranks>0.5, cTPranks<=1), axis=0)
	totcnts=np.sum(cranks<=1, axis=0)
	outfh1=open(outfn, 'w')
	outfh1.write("EventID\tEvent_type\tearly\tlate\tearlyTP\tlateTP\tTotal\tTruth\n")
	for j in xrange(len(events)): 
		e=events[j]
		outfh1.write("%s\t%s\t%s\n" % (e.id, e.determineEventType(), "\t".join(map(str, (earlycnts[j], latecnts[j], tpearlycnts[j], tplatecnts[j], totcnts[j], truth[j])))))	
	outfh1.close()


def count_early_vs_late(event, historylengths, simulation):
	#event.histories=ecycles.listout_ranges(event.histRanges)
	hindices = ecycles.historyids_to_indices(event.histories, historylengths) 
	histlens= historylengths[hindices,1]
	early=0
	late=0
	mid=0
	truth=-1
	for i in xrange(len(event.histories)): 
		h=event.histories[i]
		myord=event.orders[i]
		hlen=histlens[i]
		mytime=myord/hlen
		if h==0 and simulation:
			if mytime >0.5: 
				truth=1
			else: 
				truth=0
		else: 
			if mytime <=0.5: 
				early=early+1
			else:
				late=late+1
	return(early, late, truth)	

def get_history_lengths(events, histScores, simhist=None, usepreval=False): 
	numhists=histScores.shape[0]
	mymax=len(events)+500
	myranks=np.ones((numhists, len(events))) * mymax 
	for j in xrange(len(events)):
		e=events[j]
		hindices = ecycles.historyids_to_indices(e.histories, histScores)
		if usepreval: 
			if simhist is not None and simhist in e.histories: 
				myranks[hindices,j]=np.array(e.prevals, dtype=float)
			else: 
				myranks[hindices,j]=np.array(e.prevals, dtype=float)
		else: 
			if simhist is not None and simhist in e.histories: 
				myranks[hindices,j]=np.array(e.orders, dtype=float)
			else: 
				myranks[hindices,j]=np.array(e.orders, dtype=float)
    #change the orders into ranks. 
	hlengths=np.sum(myranks<mymax, axis=1).astype('float')
	hlengths[hlengths==1]=1.1
	xord=myranks.argsort(axis=1)
	xranks=xord.argsort(axis=1)
	cranks=xranks.astype('float') / (hlengths[:,None] -1)
	return (hlengths.astype('int'), cranks)

def main(events, pvalcutoff, histScores, statsout=None, datout=None, usepreval=False): 
	(hlengths, oranks)=get_history_lengths(events, histScores)
	if usepreval: 
		(hlengths, pranks)=get_history_lengths(events, histScores, usepreval=True)
	histScores=np.hstack((histScores, np.atleast_2d(hlengths).T))
	for j in xrange(len(events)):
		e=events[j]
		if e.likelihood > pvalcutoff:
			e.orders=oranks[np.logical_and(oranks[:,j]<=1, hlengths>0),j].tolist()
			if usepreval: 
				e.prevals=pranks[np.logical_and(pranks[:,j]<=1, hlengths>0),j].tolist()
			e.compute_timing_wmeansd(histScores)
	if statsout is not None: 
		np.savetxt(statsout, histScores, fmt='%d', delimiter='\t')
	if datout is not None: 
		outfh=open(datout, 'w')		
		outfh.write(ecycles.write_dat_header())
		for e in events:
			if e.likelihood > pvalcutoff:
				outfh.write(ecycles.write_event_dat(e))
		

if __name__ == '__main__': 
	parser=argparse.ArgumentParser(description='given a .pevnts file, will do pairwise comparison of all events, counting the number of times they come in a certain order') 
	parser.add_argument('pevntsfile', help='a .pevnts file')
	parser.add_argument('historystats', help='historystats.txt file for the sample')
	parser.add_argument('--statsout', help='file to write updated historystats to.')
	parser.add_argument('--datout', help='file to write event data to.')
	parser.add_argument('--cutoff', help='only look at events with a likelihood above this cutoff', default=0, type=float)
	args=parser.parse_args()
	histScores=np.loadtxt(args.historystats, dtype='int') 
	events=pickle.load(open(args.pevntsfile, 'rb'))
	#for e in events:
	#	e.histories=ecycles.listout_ranges(e.histRanges)
	main(events, args.cutoff, histScores, args.statsout, args.datout)
