#!/usr/bin/env python 

# For use in determining the number of TP and FP, etc, of events across a CNAVG sampling of histories given that one of the histories is the truth.  This history is what all other will be compared to.  
import os, sys
import cPickle as pickle
import argparse
import numpy as np
import re
import cnavg.linear_decomp as linear_decomp
import cnavgpost.mergehistories.event_cycles_module as histseg
import cnavgpost.diagnostics.do_order_correction as do_order_correction 


class EdgeSimulationData:
	def __init__(self, event, histScores, totalp, refhistoryid=0): 
		self.event=event
		#(segstr, self.sign)=histseg.remove_signs_from_segstr(event.segstr)
		self.sign=1
		if event.cnval <0: 
			self.sign=-1
		segstr=event.segstr
		self.cnval=event.cnval*self.sign		
		self.isTrue=0  # this will be 0 if edge is FP, 1 if TP, 2 if TN, -1 if FN, and 3 if it's a linear combination of true events.  
		self.refindex=-1
		self.refpreval=-1
		self.reforder=-1
		self.avecost=-1
		if refhistoryid in event.histories: 
			self.update_for_true_event(event, refhistoryid, histScores, totalp)
			if len(event.histories)>1:
				self.isTrue=1
			else: 
				self.isTrue=-1
		else: 
			self.isTrue=0
			event.likelihood=histseg.compute_likelihood_histories(event.histories, histScores, totalp)
			

	def __str__(self):
		event=self.event 
		prevals=",".join(map(str, [self.refpreval, event.prevalmean, event.prevalsd]))
		orders=",".join(map(str, [self.reforder, event.ordermean, event.ordersd]))
		mytype=event.CharacterizeEvent()[0]
		length=event.get_Event_length()
		mystr="\t".join(map(str, [
			event.id, 
			mytype,
			self.avecost, 
			event.likelihood, 
			event.cnval, 
			self.isTrue, 
			length, 
			prevals, 
			orders, 
			event.numhists])) + "\n"	
		return mystr

	def update_for_true_event(self, event, refhistoryid, histScores, totalp):  
		refindex=event.histories.index(refhistoryid)
		# need to pop off the simulated history values and recompute the likelihood for the event
		if len(event.histories)>1: 
			event.histories.pop(refindex)
			refpreval=event.prevals.pop(refindex)
			reforder=event.orders.pop(refindex)
			event.likelihood = histseg.compute_likelihood_histories(event.histories, histScores, totalp)	
			event.numhists=len(event.histories)
			event.compute_timing_wmeansd(histScores)
			event.histories.insert(refindex, refhistoryid)	# reinsert these 
			event.prevals.insert(refindex, refpreval)	
			event.orders.insert(refindex, reforder)	
			upperc=event.uppercosts.pop(refindex)
			lowerc=event.lowercosts.pop(refindex)
			avecost=np.mean(np.array(event.uppercosts+event.lowercosts))
			event.uppercosts.insert(refindex, upperc)	
			event.lowercosts.insert(refindex, lowerc)
		else: 
			refpreval=event.prevals[refindex]
			reforder = event.orders[refindex]
			avecost=np.mean(np.array(event.uppercosts+event.lowercosts))
			event.likelihood=0
		(self.refpreval, self.reforder, self.avecost) = (refpreval, reforder, avecost)

def analyze_simulation(events, refhistoryid, histScores, datout_fh, stats_fh, breaks_fh, outdir):
#	do_order_correction.main(events, 0, histScores, statsout=os.path.join(outdir, "historystats.txt"))
	#do_order_correction.main(events, 0, historyScores, usepreval=True)
	do_order_correction.main(events, 0, histScores)
	#make the cost of the refhistoryid 0 so that is doesn't get included in the likelihood calculation 
	histScores[np.where(histScores[:,0] == refhistoryid),:]=0	 
	totalp=histseg.compute_likelihood_histories(histScores[:,0], histScores)
	sys.stderr.write("totalp is %f\n" % totalp)
	types=histseg.Global_EVENTTYPES 
	TP=np.zeros(len(types), dtype=int)
	FP=np.zeros(len(types), dtype=int)
	TN=np.zeros(len(types), dtype=int)
	FN=np.zeros(len(types), dtype=int)
	explained=np.zeros(len(types), dtype=int)
	FNesims=[]
	FPesims=[]
	myEdgeSimData=[] 
	for event in events:
		myedgesim=EdgeSimulationData(event, histScores, totalp, refhistoryid)
		etype=event.determineEventType()
		if myedgesim.isTrue==1: 
			TP[0]+=1
			TP[etype]+=1
		elif myedgesim.isTrue==-1: 
			FN[0]+=1
			FN[etype]+=1
			FNesims.append(myedgesim)
		elif myedgesim.isTrue==0:  
			FP[0]+=1
			FP[etype]+=1
			FPesims.append(myedgesim)
		myEdgeSimData.append(myedgesim)
	check_for_linear_decomp(FNesims, FPesims, explained)
	for i in xrange(len(TN)): 
		FN[i]=FN[i]-TN[i]	
	
	if datout_fh: 
		header="event_id\tevent_type\tavecost\tLscore\tCNval\ttrue\tlength\tprevals\torders\tnumhists\n"
		datout_fh.write(header)
		for edgesim in myEdgeSimData:
			datout_fh.write(str(edgesim))
	
	if stats_fh: 
		stats_fh.write("type\ttotal\tAmp\tDel\tAdj\n")
		stats_fh.write("TP\t%s\nFP\t%s\nFN\t%s\nTN\t%s\nEX\t%s\n" % 
			("\t".join(map(str, TP)), 
			"\t".join(map(str, FP)), 
			"\t".join(map(str, FN)), 
			"\t".join(map(str, TN)),
			"\t".join(map(str, explained)) ))
		tp=TP[0]+explained[0]
		fn=FN[0]-explained[0]
		f1score = float(2*tp)/float(2*tp+fn+FP[0])
		stats_fh.write("F1Score:\t%s\n" % (str(f1score)))
	
	if breaks_fh: 
		breakpoints=histseg.get_breakpoints(edges, refhistoryid)
		for loc in breakpoints.keys(): 
			(n, t) = breakpoints[loc]
			breaks_fh.write("%s\t%d\t%d\n" % (loc, n, t))
		breaks_fh.write("Breakpoints: %d\n" % len(breakpoints))


#check if the FP event is a linear combination of True events and, 
# likewise if the FN events can be explained by the FP events.
def check_for_linear_decomp(FNesims, FPesims, explained):
	fpexp=[] # a list of the FP events that can be explained. 
	if len(FNesims) >0: 
		truelist=[] 
		for esim in FNesims: 
			truelist.append(create_edge_tuplelist(esim.event))
		RV=linear_decomp.ReferenceVectors(truelist)
		for fpesim in FPesims: 
			fplist=create_edge_tuplelist(fpesim.event)
			if RV.canExplain(fplist): 
				explained[0]+=1
				fpevent=fpesim.event
				explained[fpevent.determineEventType()]+=1
				fpesim.isTrue=3
		for esim in FNesims: 
			if create_edge_tuplelist(esim.event) in RV.used_elements: 
				esim.isTrue=-3
			
		#test_FN_events_for_linear_decomp(FNesims, fpexp)

def test_FN_events_for_linear_decomp(FNesims, fpexp): 	
# now correct the TN events that can be explained by FP explained events
	if len(fpexp) >0: 
	#find the combination of FN events that explains each FP event
	# do this by a leaving one out until you get to the reduced set that can explain the event
		for i in xrange(len(FNesims)): 
			FNsubset=FNesims[:i]+FNesims[(i+1):]
			truelist=[]
			for esim in FNsubset:
				truelist.append(create_edge_tuplelist(esim.event))
			RV=linear_decomp.ReferenceVectors(truelist)
			for esim in fpexp: 
				fnlist=create_edge_tuplelist(esim.event)
				if RV.canExplain(fnlist):
					test_FN_events_for_linear_decomp(FNsubset, fpexp)
				else: 
					FNesims[i].isTrue=-3


def create_edge_tuplelist(event):
	segstr=event.segstr
	locs=segstr.split(',')
	newlocs=[]
	sign=1
	if event.cnval <0: sign =-1
	for (i, loc) in enumerate(locs):
		seg=False
		m=re.search('(\w+):(-?\d+)\(([+|-])\)-(\w+):(-?\d+)\(([+|-])\)', loc)
		if m:
			(chr1, s, st1, chr2, e, st2)=m.group(1,2,3,4,5,6)
			# order the ends
		else: 
			m=re.search('(\w+):(-?\d+)-(\-?\d+)', loc)
			if m: 
				(chr1, s, e) = m.group(1,2,3)
				(chr2, st1, st2)=(chr1, "+", "+")
				seg=True
		if (chr1 > chr2) or (chr1==chr2 and s>e): 
			(x,y,z)=(chr1, s, st1)
			(chr1, s, st1)=(chr2, e, st2)
			(chr2,e,st2)=(x, y,z)
		mysign =sign 
		if (i % 2) ==1: mysign = sign * -1
		if seg: newstr="%s:%s-%s" % (chr1, s, e)
		else: newstr="%s:%s(%s)-%s:%s(%s)" % (chr1, s, st1, chr2, e, st2)
		newlocs.append((newstr, mysign))
	return(newlocs)


def checkForCancellingEdges(edgesims):
	sortededgesims=sorted(edgesims, key=lambda x: (x.segstr, x.edge.cnval))  
	preseg=sortededgesims[0]
	numpairs=0
	TN=[0,0,0,0]
	for edgesim in sortededgesims[1:]: 
		if edgesim.segstr == preseg.segstr and (edgesim.cnval == -1 * preseg.cnval):
			numpairs+=1
			edgesim.isTrue=2
			preseg.isTrue=2
			TN[0]+=2
			TN[preseg.type]+=1
			TN[edgesim.type]+=1		
		preseg=edgesim
	return TN

def add_options(parser):
	parser.add_argument('pedgs', help='a .pedgs file')
	parser.add_argument('trueID', help='the ID of the true history', type=int)
	parser.add_argument('historystats', help='The file of history stats')
	parser.add_argument('--datout', help='the file to write data to', type=argparse.FileType('w')) #, default=sys.stdout) 
	parser.add_argument('--stats', help='the file to write stats to.', type=argparse.FileType('w')) #, default=sys.stderr) 
	parser.add_argument('--breaks', help='the file to write breaklocations to.', type=argparse.FileType('w')) #, default=sys.stderr) 
	parser.add_argument('--binwidth', dest='binwidth', help='the multiplier between history ids of independent runs', default=histseg.Global_BINWIDTH, type=int)	

if __name__== '__main__': 
	parser = argparse.ArgumentParser(description='does analysis for simulation tests')
	add_options(parser)
	args=parser.parse_args()
	histseg.Global_K=0  #this used to be 1 initially.  
	histseg.Global_BINWIDTH=args.binwidth
	allevents=pickle.load(open(args.pedgs, 'rb'))
	historyScores=np.loadtxt(args.historystats, dtype=int)
	analyze_simulation(allevents, args.trueID, historyScores, args.datout, args.stats, args.breaks, "./" )

