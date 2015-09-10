#!/usr/bin/env python 

import argparse 
import cPickle as pickle
import sys, os 
import cnavgpost.mergehistories.event_cycles_module as ecycles 
import numpy as np
import re

def get_events_order_counts(myevents, outfn, simulation, histScores=None, shufstart=0):
	outfh=open(outfn, 'w')
	outfh.write("eventAid\teventBid\tAtype\tBtype\ttruth\tabcost1\tbacost1\ttot1\tabcnt1\tbacnt1\tabcost2\tbacost2\ttot2\tabcnt2\tbacnt2\n")
	for a in xrange(len(myevents)):
		eventA = myevents[a]
		for b in xrange(a+1, len(myevents)): 
			eventB = myevents[b]
#			(ab, ba, tot, truth) = get_order_counts(eventA, eventB, simulation) 
			results= get_hcost_difference(eventA, eventB, simulation, histScores, shufstart) 
			#outfh.write("%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n" % (eventA.id, eventB.id, ab, ba, tot, truth, eventA.determineEventType(), eventB.determineEventType()))
			outfh.write("\t".join(map(str, [eventA.id, eventB.id, eventA.determineEventType(), eventB.determineEventType(), results]))+"\n")
	outfh.close()

def get_hcost_difference(eventA, eventB, simulation, historyScores, shufstart): 
	shuffledat={'abcost': 0, 'bacost': 0, 'tot':0, 'abcnt': 0, 'bacnt':0}
	noshuffdat={'abcost': 0, 'bacost': 0, 'tot':0, 'abcnt': 0, 'bacnt':0}
	truth=-1
	ai=np.in1d(np.array(eventA.histories), np.array(eventB.histories))
	if sum(ai)>0: 
		hids=np.array(eventA.histories)[ai]
		# hids now has the history ids where both eventA and eventB exist together 
		if simulation and 0 in hids:
			hids=hids[hids!=0] 
			atruth=eventA.orders[eventA.histories.index(0)]	
			btruth=eventB.orders[eventB.histories.index(0)]
			if atruth<btruth: 
				truth=1
			elif btruth < atruth: 
				truth=2
			else: truth=0
		# split hids into those before shuffling and after shuffling. 
		(hids1, hids2) = split_hids_before_after_iteration(hids, shufstart)
		for (hids, dat) in [(hids1, noshuffdat), (hids2, shuffledat)]:
			ai=np.in1d(np.array(eventA.histories), hids)
			bi=np.in1d(np.array(eventB.histories), hids)
			aorders=np.array(eventA.orders)[ai]
			borders=np.array(eventB.orders)[bi]
			abi=np.array(eventA.histories)[ai][aorders<borders]
			dat['abcnt']=len(abi)
			if len(abi)>0: 
				hi=ecycles.historyids_to_indices(abi, historyScores)
				dat['abcost']=np.mean(historyScores[hi,1:3])
			bai=np.array(eventA.histories)[ai][borders<aorders]
			dat['bacnt']=len(bai)
			if len(bai)>0: 
				hi=ecycles.historyids_to_indices(bai, historyScores)
				dat['bacost']=np.mean(historyScores[hi,1:3])
			eqhids=np.array(eventA.histories)[ai][borders==aorders]
			if len(eqhids)>0: 
				hi=ecycles.historyids_to_indices(eqhids, historyScores)
				dat['eqcost']=np.mean(historyScores[hi,1:3])
			dat['tot']=sum(ai)
	keys=['abcost', 'bacost', 'tot', 'abcnt', 'bacnt']
	mystr=str(truth)
	for dat in [shuffledat, noshuffdat]: 
		for k in keys: 
			mystr+="\t"+str(dat[k])
	return(mystr) 

def split_hids_before_after_iteration(hids, tsplit): 
	itr=np.fmod(hids, ecycles.Global_BINWIDTH)
	hids1=hids[itr<tsplit]
	hids2=hids[itr>=tsplit]
	return(hids1, hids2)


def get_order_counts(eventA, eventB, simulation): 
	ab=0
	ba=0
	truth=-1
	#eventA.histories=ecycles.listout_ranges(eventA.histRanges)	
	#eventB.histories=ecycles.listout_ranges(eventB.histRanges)
	aindices=[]
	bindices=[]
	for ai in xrange(len(eventA.histories)): 
		h=eventA.histories[ai]
		if h in eventB.histories: 
			aindices.append(ai)
			bindices.append(eventB.histories.index(h))
	aorders=np.array(eventA.orders)[aindices]
	borders=np.array(eventB.orders)[bindices]
	ab=sum(aorders<borders)
	ba=sum(aorders>borders)
	tot=len(aorders)
	if simulation and (0 in eventA.histories) and (0 in eventB.histories): 
		atruth=eventA.orders[eventA.histories.index(0)]	
		btruth=eventB.orders[eventB.histories.index(0)]
		tot=tot-1 # don't count the simulated history
		if atruth<btruth: 
			truth=1
			ab=ab-1 
		elif btruth<atruth: 
			truth=2
			ba=ba-1
		else: 
			truth=0
	return (ab, ba, tot, truth)	

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

def get_shuffle_history_scores(histscores, shuftime, simulation): 
	if simulation: 
		hids=histscores[histscores[:,0]>0,0]
	itr=np.fmod(hids, ecycles.Global_BINWIDTH)
	(hids1, hids2)=(hids[itr<shuftime], hids[itr>=shuftime])
	results=[]
	for hids in (hids1, hids2): 
		i=ecycles.historyids_to_indices(hids, histscores)
		a=histscores[i,1:3]
		b=np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
		x, idx = np.unique(b, return_index=True)
		results.append(a[idx])
	return results
	

def main(pevntsfile, outdir, simulation, pvalcutoff, histstatsfn, shufstart): 
	sys.stderr.write("pvalcutoff is %f\n" % pvalcutoff)
	useEdges=re.search(".pedgs", pevntsfile)
	historyScores=np.loadtxt(histstatsfn, dtype='int') #You really only need this file in order to tell how many histories there are.
	results=get_shuffle_history_scores(historyScores, shufstart, simulation) 
	np.savetxt("cnranges.txt", results[0], fmt='%d', delimiter='\t')
	np.savetxt("shranges.txt", results[1], fmt='%d', delimiter='\t')
	sys.exit()
	events=pickle.load(open(pevntsfile, 'rb'))
	myevents=[]
	for e in events:
		#e.histories=ecycles.listout_ranges(e.histRanges) 
		if e.likelihood > pvalcutoff:
			myevents.append(e)
	else:
		myevents=events
	sys.stderr.write("There are %d events with pvalue gt %f\n" % (len(myevents), pvalcutoff))
	for e in myevents: 
		sys.stderr.write(e.id+"\n")
	if True: 
		if useEdges:  
			outfn=os.path.join(outdir, "edges_ordcnts.dat")
		else: 
			outfn=os.path.join(outdir, "event_ordcnts.dat")
		get_events_order_counts(myevents, outfn, simulation, historyScores, shufstart)
	if False: 
		if useEdges:  
			outfn1=os.path.join(outdir, "edges_earlycnts.dat")
		else:
			outfn1=os.path.join(outdir, "event_earlycnts.dat")
		outfn2=os.path.join(outdir, "histlengths.dat")
		count_earlylate_with_correction(myevents, historyScores, outfn1, outfn2) 
		

if __name__ == '__main__': 
	parser=argparse.ArgumentParser(description='given a .pevnts file, will do pairwise comparison of all events, counting the number of times they come in a certain order') 
	parser.add_argument('pevntsfile', help='a .pevnts file')
	parser.add_argument('historystats', help='historystats.txt file for the sample')
	parser.add_argument('outdir', help='The directory to write the results to. A file called event_ordcnts.dat will be made and, if historylengths is specified, one called event_earlycnts.dat')
	parser.add_argument('--simulation', help='whether the history is a simulation', action='store_true') 
	parser.add_argument('--cutoff', help='only look at events with a likelihood above this cutoff', default=0, type=float)
	parser.add_argument('--shufstart', default=0, type=int, help='time when order shuffling begins and cn-avg sampleing ends.')
	args=parser.parse_args()
	main(args.pevntsfile, args.outdir, args.simulation, args.cutoff, args.historystats, args.shufstart)
