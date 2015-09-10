# Module for merging histories created from the CN-AVG pipeline together
# Events are equivalent if they have identical coordinates across all segments and adjacencies and have the same CN value.  Basically, it has to be the same graph structure with the same Copy number change (The CN will be an integer).

import sys, os
import re, glob, subprocess
import itertools 
import math 
import copy
import cPickle as pickle 
from collections import Counter
import numpy as np
import pysam, gzip
from cnavgpost.mergehistories.braney_lines_module import * 

Global_BINWIDTH=10000
Global_MAXCOST=300 
Global_K=0
Global_SPLITOFFS=True # include duplicate events if an event occurs twice in the same history
Global_SPLITCYCLES=False#do we want to split all figure 8 types into smaller cycles? 
Global_EVENTTYPES = ['any', 'amp', 'del', 'oth', 'amdl']

# an Event is made up of multiple Braneysegs and some extra info
class Event:
	def __init__(self, bseglist):
		if isinstance(bseglist, str): 
			data=bseglist.strip().split('\t')
			self.id=data[0]  # The id is usually the Run id and the first iteration of the run that the event is found in.  It has the form [int].[int]
			self.likelihood=data[1]  # This is a float. 
			(self.numhists, self.numsims) = map(int, data[2].split(','))
			(self.prevalmean, self.prevalsd)=map(float, data[3].split(','))
			(self.ordermean, self.ordersd)=map(float, data[4].split(','))
			self.cnval=float(data[5])
			self.segstr=data[6]
			self.indyRunCounts=data[7]
			self.numsegs=int(data[8])
			self.numadjs=0
			self.histories=[] 
			self.prevals=[] 
			self.orders=[] 
			self.segs=[]
			self.uppercosts=[] 
			self.lowercosts=[]
			self.dupsremoved=True
			
		elif isinstance(bseglist, list):
			self.id="%d.%d" % (bseglist[0].historyid, bseglist[0].order)
			self.likelihood=0
			(self.numhists, self.numsims)=(1,1)
			(self.prevalmean, self.prevalsd)=(bseglist[0].preval, "NA")
			(self.ordermean, self.ordersd)=(bseglist[0].order, "NA")
			self.cnval=bseglist[0].cnval
			self.segstr=""
			(self.numsegs, self.numadjs) = (0 , 0)	
			self.histories=[bseglist[0].historyid]
			self.indyRunCounts=""
			self.prevals=[bseglist[0].preval]
			self.orders=[bseglist[0].order]
			self.segs=bseglist
			self.uppercosts=[bseglist[0].upperEventCost]
			self.lowercosts=[bseglist[0].lowerEventCost]
			self.dupsremoved=False
	
	def __str__(self):
		fstr="%s\t%s\t%d,%d\t%s,%s\t%s,%s\t%f\t%s\t%s\t%d\t%d\n" % (self.id, str(self.likelihood), self.numhists, self.numsims, self.prevalmean, str(self.prevalsd), str(self.ordermean), str(self.ordersd), self.cnval, self.segstr, self.indyRunCounts, self.numsegs, self.numadjs)
		return fstr
	
	def __eq__(self, other): 
		if self.segstr == "" : self.make_segstr()
		if other.segstr == "" : other.make_segstr()
		if self.segstr == other.segstr and self.cnval == other.cnval: 
			return True
		else:
			return False

	# update is done after equivalent events across multiple histories have been merged together.  Then the likelihood score and other stats for this event can be calculated. 	
	def update(self, historyScores, totalp=1):
		if self.segstr=="": 
			self.make_segstr()
		if self.numsegs==0: 
			self.count_segs()
		self.sort_values_by_history()
		self.id = "%d.%d" % (self.histories[0], self.orders[0]) 
		self.numhists=len(self.histories)
		(self.numsims, self.indyRunCounts) = getIndRunCounts(self.histories)
		if historyScores is not None: 
			self.compute_timing_wmeansd(historyScores)
			self.likelihood=compute_likelihood_histories(self.histories, historyScores, totalp)

	def addseg(self, bseg):
		self.segs.append(bseg)
	
	#Trimming of the event is done to make it more compact before it is written out to a pickled file. 
	def trim(self): 
		if self.segstr=="":
			self.make_segstr()
		self.segs=[]

	# Unpacking of an event is done after it's been read in from a pickled file.  This basically undoes the trim (above). 
	def unpack(self):
		if not self.segs:
			self.make_segs_from_str()
	
	# This makes a string from the genomic coordinates of the event (event being a flow in the CN-AVG).  
	# The CN for the event will have the sign of whatever the first edge is.  The sign of the edges alternate after that. 	
	def make_segstr(self):
		self.remove_dup_adj()
		self.merge_adjacent_segs()
		self.ordersegs()
		seg0=self.segs[0]
		self.cnval = seg0.cnval 
		mysegs=[]
		for seg in self.segs: 
			if seg.adj: 
				s="%s:%d(%s)-%s:%d(%s)" % (seg.chr, seg.start, seg.st1, seg.chr2, seg.end, seg.st2)	
			else: 
				s="%s:%d-%d" % (seg.chr, seg.start, seg.end)
			mysegs.append(s)	
		self.segstr=",".join(mysegs) 
	
	# This orders the segments of the event, putting the lowest genomic coordinate first and then going around the rearrangement cycle. 
	def ordersegs(self):
		#first order segs so the start is lower than the end 
		for seg in self.segs: #If there's only one segment, we are dealing with .edges and want them to have the lower coordinate first for consistency  
			if (seg.chr ==seg.chr2 and seg.start > seg.end) or (seg.chr>seg.chr2): 
				seg.flip_ends()	
		#If there's only one segment, we don't have the worry about ordering it.  
		if len(self.segs)>1:  
			firstseg=sorted(self.segs, key=lambda x: (x.chr, x.start, x.chr2, x.end, x.st1, x.st2))[0]
			sortedsegs=sorted(self.segs, key=lambda x: x.cycleorder)
			ifirst=sortedsegs.index(firstseg)
			second = sortedsegs[ifirst-len(sortedsegs)+1]
			last = sortedsegs[ifirst-1]
			direction=get_direction(firstseg, second, last)
			if direction==1: 
				self.segs=sortedsegs[ifirst:]+sortedsegs[:ifirst]
			elif direction==-1: 
				self.segs=sortedsegs[:(ifirst+1)][::-1]+sortedsegs[(ifirst+1):][::-1]
			else: #direction==0
				self.segs=sortedsegs[ifirst:]+sortedsegs[:ifirst]
				sys.stderr.write("Don't know which way to go for this event!\n%s" % str(self))
				for s in self.segs: 
					sys.stderr.write("%d\t%d\t%d\t%d\n" % (s.seg, s.start, s.end, s.cycleorder))
			# want to flip the ends of segments so they read in order. 
			firstseg=self.segs[0]
			for seg in self.segs[1:]: 
				if (firstseg.chr2== seg.chr2 and firstseg.end==seg.end): 
					seg.flip_ends()
				firstseg=seg
	
	def merge_adjacent_segs(self):
		sortedsegs=sorted(self.segs, key=lambda x: x.cycleorder)
		l=len(sortedsegs)
		i=0
		while i < l:  
			seg = sortedsegs[i]
			if (seg.adj and (seg.chr == seg.chr2) and 
				(seg.start== seg.end+1 or seg.start==seg.end-1)):
				pseg=sortedsegs[i-1]
				if (i+1)==l: 
					nseg=sortedsegs[0]
				else: 
					nseg=sortedsegs[i+1]
				if nseg.seg and pseg.seg:
					minloc=min(nseg.start, nseg.end, pseg.start, pseg.end)
					maxloc=max(nseg.start, nseg.end, pseg.start, pseg.end)
					pseg.start=minloc
					pseg.end=maxloc
					sortedsegs.pop(i) #get rid of the adjacency
					if (i+1)==l:
						sortedsegs.pop(0) #get rid of the next segment
					else: 
						sortedsegs.pop(i)
					l=len(sortedsegs)
				else: 
					i+=1
			else: 
				i+=1
		self.segs=sortedsegs

	def make_segs_from_str(self):
		mysegs=self.segstr.split(",")
		dummysegline="%s\t%s\t%s\t%f\t%f\t%d\t0\t0\t%d\t%d\t0\t0\t0\t0\t0\n"
		dummyadjline="A\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%d\t0\t0\t%d\t%d\t0\t0\t0\t0\t0\n"
		self.segs=[]
		cycleorder=0
		minhist=0
		for (i, loc) in enumerate(mysegs): 
			m=re.search('(\w+):(-?\d+)-(-?\d+)', loc)			
			if m: 
				coords=m.group(1,2,3)
				bseg=Braney_seg(dummysegline % (coords[0], coords[1], coords[2], self.cnval*self.prevalmean, self.prevalmean, minhist, cycleorder, self.ordermean))
				#switch the sign of the CN for odd segs in the cycle
				if (i % 2) == 1: bseg.cnval = bseg.cnval * -1  
			else: 
				m=re.search('(\w+):(-?\d+)\(([+|-])\)-(\w+):(-?\d+)\(([+|-])\)', loc)
				if m: 	
					coords=m.group(1,2,3,4,5,6)
					bseg=Braney_seg(dummyadjline % (coords[0], coords[1], coords[2], coords[3], coords[4], coords[5], self.cnval * self.prevalmean, self.prevalmean, minhist, cycleorder, self.ordermean))
					if (i % 2) == 1: bseg.cnval = bseg.cnval * -1  
				else: 
					m=re.search('([+|-])/(\w+):(-?\d+)-(-?\d+)', loc)			
					if m: 
						coords=m.group(2,3,4)
						bseg=Braney_seg(dummysegline % (coords[0], coords[1], coords[2], self.cnval, self.prevalmean, minhist, cycleorder, self.ordermean))
						sign=m.group(1)
						if sign=="+": bseg.cnval = bseg.cnval
					else: 
						m=re.search('([+|-])/(\w+):(-?\d+)\(([+|-])\)-(\w+):(-?\d+)\(([+|-])\)', loc)	
						if m: 
							coords=m.group(2,3,4,5,6,7)
							bseg=Braney_seg(dummyadjline % (coords[0], coords[1], coords[2], coords[3], coords[4], coords[5], self.cnval, self.prevalmean, minhist, cycleorder, self.ordermean))
							sign=m.group(1)
							if sign=="-": bseg.cnval = bseg.cnval * -1
			self.segs.append(bseg)
			cycleorder+=1
	
	def add_Event_data(self, other):
		self.histories+= other.histories
		self.uppercosts+= other.uppercosts
		self.lowercosts += other.lowercosts
		self.prevals += other.prevals 
		self.orders += other.orders 
		
	def merge_Event_data(self, other): 
		iToAdd=get_index_of_non_intersecting_items(self.histories, other.histories)
		self.histories+=list(np.array(other.histories)[iToAdd])
		self.uppercosts+=list(np.array(other.uppercosts)[iToAdd])
		self.lowercosts+=list(np.array(other.lowercosts)[iToAdd])
		self.prevals+=list(np.array(other.prevals)[iToAdd])
		self.orders+=list(np.array(other.orders)[iToAdd])
	
	def sort_values_by_history(self):
		histarray=np.array(self.histories) 
		myorderi = list(np.argsort(np.array(self.histories)))
		self.histories=[self.histories[i] for i in myorderi]
		self.uppercosts=[self.uppercosts[i] for i in myorderi]
		self.lowercosts=[self.lowercosts[i] for i in myorderi]
		self.prevals=[self.prevals[i] for i in myorderi]
		self.orders=[self.orders[i] for i in myorderi]

	def remove_dup_adj(self): 
		if not self.dupsremoved:
			nodup_segs=[]
			for seg in self.segs: 
				addseg=True
				if seg.adj: # if the braney_seg is an adjacency, see that it's not a duplicate of one that's been added already (there are two lines for every adjacency in .braney files)
					for adj in nodup_segs: 	
						addseg = (addseg and seg.cycleorder != adj.cycleorder) # seg.is_dup(adj))
					if addseg: 
						nodup_segs.append(seg)
				else: nodup_segs.append(seg)
			self.segs=nodup_segs
		self.dupsremoved=True

#Global_EVENTTYPES = ['any', 'amp', 'del', 'oth', 'amdl']
	def determineEventType(self): 
		mytype=0
		if len(self.segs)==0: 
			self.make_segs_from_str()
		for seg in self.segs: 
			if seg.seg: 
				if seg.cnval <0 and mytype != 2: 
					mytype=1
				elif seg.cnval >0 and mytype !=1: 
					mytype=2
				else: 
					mytype=4
		if mytype==0:
			mytype=3
		return mytype
	
	def CharacterizeEvent(self): 
		mytype='any'
		mychr=None
		intra=True
		numfuse=0
		numfiss=0
		numsegs=0
		if len(self.segs) ==0: 
			self.make_segs_from_str()
		if len(self.segs)>4: 
			type='oth'
		for seg in self.segs:
			if mychr is None and seg.chr != "None": 
				mychr=seg.chr
			elif (seg.chr != "None") and (mychr != seg.chr): 
				intra=False
			if seg.seg:
				numsegs+=1
				if seg.cnval<0 and mytype !='del':#the sign is opposite for segments
					mytype='amp'
				elif seg.cnval>0 and mytype !='amp': 
					mytype='del'
				else: 
					mytype='amdl'
			elif seg.isbreakpoint(): 
				if seg.cnval<0: 
					numfuse+=1
				else: 
					numfiss+=1
		if numsegs==0 and len(self.segs)==4: 
			mytype=determine_4adj_cycle_type(self.segs, intra)
		return (mytype, intra, numfuse, numfiss, numsegs, len(self.segs)) 
	
	def get_Event_length(self):
		seglen=0
		if len(self.segs)==0: 
			self.make_segs_from_segstr()
		for seg in self.segs: 
			if seg.seg: 
#				seglen=max(seglen, seg.end-seg.start+1)
				seglen+=(abs(seg.end-seg.start)+1)
#			else:
#				if seg.chr==seg.chr2:  
#					adjlen=max(adjlen, seg.end-seg.start)
		if seglen>0:
			return seglen
		else: 
			return 0	

# this just counts the number of segments vs adjacencies for an event. 
	def count_segs(self): 
		if not self.dupsremoved: 
			self.remove_dup_adj() 
		myadjs=[]
		(numsegs, numadjs)= (0,0)
		if len(self.segs) == 0: 
			self.make_segs_from_str()
		for seg in self.segs: 
			if seg.adj: 
				numadjs+=1
			else: 
				numsegs+=1
		self.numsegs=numsegs
		self.numadjs=numadjs
	
	def get_chr_list(self): 
		mychr=[]
		for seg in self.segs:
			if seg.chr != "None":  
				mychrs.append(seg.chr)
			if seg.chr2 != "None": 
				mychrs.append(seg.chr)
		return "\t".join(map(str, set(mychrs)))

	def check_overlap(self, chr, start, end):
		for seg in self.segs:
			if seg.seg: 
				seg.order_ends()
				if (seg.chr == chr and seg.start < end and seg.end>start): 
					return True
			elif seg.adj: 
				if ((seg.chr == chr and (seg.start <= end and seg.start >= start)) or 
					(seg.chr2 == chr and (seg.end <= end and seg.end >= start))):
					return True 
		return False
	
	def multiline_str(self):
		mystr=""
		for seg in self.segs: 
			mystr+="%s\t%s\t%d\t%s\t%s\t%f\n" % (str(seg), self.id, len(self.histories), self.hranges, str(self.likelihood))
		return mystr

	def compute_likelihood(self, historyScores, totalp=1):
		self.likelihood=compute_likelihood_histories(self.histories, historyScores, totalp)
	
	def compute_timing_wmeansd(self, historyScores):
		if len(self.histories)>0: 
			hindices = historyids_to_indices(self.histories, historyScores)
			costsarray=historyScores[hindices, 1]
			w=np.exp(-1*Global_K*np.array(costsarray))
			pvals=np.array(self.prevals, dtype=float)
			waverage=np.average(pvals, weights=w)
			var=np.average((pvals - waverage)**2, weights=w)
			self.prevalmean=waverage
			self.prevalsd=math.sqrt(var)
			# do the same calculation for the ordering of the event within the history
			pvals=np.array(self.orders, dtype=float)
			waverage=np.average(pvals, weights=w)
			var=np.average((pvals - waverage)**2, weights=w)
			self.ordermean=waverage
			self.ordersd=math.sqrt(var)
		else: 
			self.prevalmean=-1
			self.prevalsd=-1
			self.ordermean=-1
			self.ordersd=-1
			
def sort_segs_in_cycle(seglist):
	segs=sorted(seglist, key=lambda x: x.cycleorder)
	currentseg=segs.pop(0)
	myorderedsegs=[currentseg]
	i=0
	maxiter=len(segs)**2
	tmpsegs=[]
	while (len(segs)>0 and maxiter >=0):
		seg=segs.pop(0)
		maxiter= maxiter -1 
		if (currentseg.chr2 == seg.chr2 and currentseg.end==seg.end): 
			seg.flip_ends()
			myorderedsegs.append(seg)
			currentseg=seg
		elif (currentseg.chr2 == seg.chr and currentseg.end==seg.start): 
			myorderedsegs.append(seg)
			currentseg=seg
		else: 
			segs.append(seg)
	if maxiter <0: 
		sys.stderr.write("len of segs: %d, and segstr: %s" % (len(segs), currentseg))
		sys.exit(-1)
	return myorderedsegs

def remove_signs_from_segstr(segstr):
	locs=segstr.split(',')
	newlocs=[]
	sign=""
	for loc in locs:
		m=re.search('([+|-])/(\w+):(-?\d+)\(([+|-])\)-(\w+):(-?\d+)\(([+|-])\)', loc)
		if m:	#(chr1, s, st1, chr2, e, st2)=m.group(2,3,4,5,6,7)
			newlocs.append("%s:%s(%s)-%s:%s(%s)" % m.group(2,3,4,5,6,7))
			if sign=="": sign=m.group(1) #take the sign of the first segment
		else:
			m=re.search('([+|-])/(\w+):(-?\d+)-(\-?\d+)', loc)
			if m: #(chr1, s, e) = m.group(2,3,4)
				newlocs.append("%s:%s-%s" % m.group(2,3,4))
				if sign=="": sign=m.group(1)
	mysign=1
	if sign=="-": mysign=-1
	return (",".join(newlocs), mysign)

# this will modify e1 and e2. 
def cancel_Event_data(e1, e2): 
	keepi2=get_index_of_non_intersecting_items(e1.histories, e2.histories)
	keepi1=get_index_of_non_intersecting_items(e2.histories, e1.histories)
	for (e, keepi) in [(e1, keepi1), (e2, keepi2)]: 
		e.histories=list(np.array(e.histories)[keepi])
		e.uppercosts=list(np.array(e.uppercosts)[keepi])
		e.lowercosts=list(np.array(e.lowercosts)[keepi])
		e.prevals=list(np.array(e.prevals)[keepi])
		e.orders=list(np.array(e.orders)[keepi])
	

# If the given event (aka cycle) has figure 8 structures, it will split the event into the constitutative smaller cycles. This is used to break figure 8s up in the simulations so that the events are comparable to the true history in which figure 8s aren't allowed.  
def split_cycle(event):
	if len(event.segs) == 0:
		event.make_segs_from_str()
	for ia in xrange(len(event.segs)):
		sega=event.segs[ia]
		for ib in xrange(ia): 
			segb=event.segs[ib]
	#		sys.stderr.write("ia %d, sega.start: %d, ib %d, segb.start: %d\n" % (ia, sega.start, ib, segb.start))
			if ((((ia - ib) % 2) == 0) and sega.start == segb.start and sega.chr == segb.chr): 
				newevent=copy.deepcopy(event)
				newevent.segs=event.segs[ib:ia]
				newevent.make_segstr()
				newevent.count_segs()
				event.segs=event.segs[:ib]+event.segs[ia:]
				event.make_segstr()
				event.count_segs()
				sys.stderr.write("splitting at %d, %d\nnewevent: %d oldevent: %d\n" % (ib, ia, len(newevent.segs), len(event.segs)))
				return[newevent, event]
	return [event]

def split_cycles(events): 
	to_split=events
	split=[]
	while len(to_split)>0: 
		cycle= to_split.pop()
		decomposition=split_cycle(cycle)
		if len(decomposition)==1: 
			#cycle couldn't be broken up
			split.extend(decomposition)
		else: 
			# cycle was broken up, and we need to check the components now
			to_split.extend(decomposition)
	return split 

def listout_ranges(ranges): 
	myilist=[]
	for (s,e) in ranges: 
		for i in range(s,e+1):
			myilist.append(i)
	return myilist

def ranges(ilist): 
	G=(list(x) for _,x in itertools.groupby(sorted(ilist), lambda x, c=itertools.count(): next(c)-x))
	return ",".join("-".join(map(str, (g[0], g[-1])[:len(g)])) for g in G)

def get_index_of_non_intersecting_items(list1, list2): 
	myis=[]
	for i in xrange(len(list2)): 
		if list2[i] not in list1: 
			myis.append(i)
	return myis

def get_index_of_intersecting_items(list1, list2): 
	myis=[]
	for i in xrange(len(list2)): 
		if list2[i] in list1: 
			myis.append(i)
	return myis

def getRanges(vals):
	myranges=[] 
	vals.sort()
	rangestart=vals[0]
	rangeend=vals[0]
	for i in xrange(len(vals)-1): 
		if vals[i+1] > vals[i]+1: 
			rangeend=vals[i]
			myranges.append((rangestart, rangeend))
			rangestart=vals[i+1]
	rangeend=vals[-1]
	myranges.append((rangestart, rangeend))
	return myranges 
	
def getIndRunCounts(histories):
		binwidth=Global_BINWIDTH
		numbins=int(max(histories)/binwidth)+1
		hranges=[0]*numbins
		for h in histories: 
			bval=int(h/binwidth)
			hranges[bval]+=1
		vals=[]
		numsims=0
		for i in xrange(numbins): 
			vals.append("%d:%d" % (i, hranges[i]))
			if hranges[i]>0: numsims+=1
		hranges=",".join(vals)
		return (numsims, hranges)	

def get_cnvalue_from_seglist(bseglist):
	cnval=None
	for bseg in bseglist: 
		if bseg.seg: 
			cnval=round(bseg.cnval/bseg.preval, 2)
	if cnval==None: 
		cnval=round(abs(bseg.cnval/bseg.preval), 2)
	return cnval 

 
#### get_overlapping_events #######################
# Input: a chromosome, start and end, and the filename of a .evnts file 
# Output: a list of Events, where either an adjancency or a segment overlaps the region 
def get_overlapping_events(chr, start, end, evntsfn ):
	overlapping_events=[]
	for line in open(evntsfn, 'r'):
		myevent=Event(line)
		if (myevent.check_overlap(chr, start, end)):
			overlapping_events.append(myevent)

#### get_overlapping_events_tabix #######################
# Input: a chromosome, start and end, and the filename of a .braney file in tabix form
# Output: a list of Events, where either an adjancency or a segment overlaps the region 
def get_overlapping_events_tabix(chr, start, end, tabixfn, tabixfn2):
	eventid=""
	bsegs=[]
	eventlines=[]
	myevents=[] # a list of Events that contain a segment overlapping our genomic region
	tabixfile=pysam.Tabixfile(tabixfn)
	bseglines=[]
	try: 
		bseglines=tabixfile.fetch(reference=chr, start=start, end=end)
	except ValueError: 
		sys.stderr.write("Error in tabix fetch for %s, with %s:%d-%d, %s\n" % (tabixfn, chr, start, end, str(ValueError)))
	for bline in bseglines: 
		bsegs.append(Braney_seg(bline))
	badjlines=[]
	if (tabixfn2): 
		tabixfile2=pysam.Tabixfile(tabixfn2)
		try: 
			badjlines=tabixfile2.fetch(reference=chr, start=start, end=end)
		except ValueError: 
			sys.stderr.write("Error in tabix fetch for %s, %s\n" % (tabixfn2, str(ValueError)))
		for bline in badjlines: 
			bsegs.append(Braney_seg(bline))
	if len(bsegs)>0:
		bsegs.sort(key=lambda x: x.ptrid)
		for bseg in bsegs:
			if eventid == "": 
				eventid=bseg.ptrid
				eventlines.append(bseg) 
			elif bseg.ptrid == eventid: 
				eventlines.append(bseg)
			elif bseg.ptrid != eventid: 
				myevents.append(Event(eventlines))
				eventlines=[bseg]
				eventid=bseg.ptrid
		myevents.append(Event(eventlines))
	return myevents

# Input: a file handle to an eventsfile sorted by genomic description
# Output: a list of events where duplicate ones are merged together. 
def unique_c_events_sorted_file(evntsfile):
	unique_events=[]
	line1=evntsfile.readline()
	eventA=Event(line1)
	for line in evntsfile:
		eventB=Event(line)
		if eventA == eventB:
			eventA.add_Event_data(eventB) 
		else:
			eventA.update(None)
			unique_events.append(eventA)
			eventA=eventB
	unique_events.append(eventA)
	return unique_events

def get_split_offs(finalevents):
	splitoffs=[]
	for evnt in finalevents:
		tmpevents=split_off_duplicate(evnt)
		n=len(tmpevents)
		while n>0: 
			xlist=split_off_duplicate(evnt)
			tmpevents+=xlist
			n=len(xlist)
		splitoffs+=tmpevents
	return splitoffs 

def split_off_duplicate(event):
	hcounts=Counter(event.histories).items()
	dups=[h for h, k in Counter(event.histories).items() if k>1]
	if len(dups)>=1:
		newevent=copy.deepcopy(event)
		newevent.histories=[]
		newevent.prevals=[]
		newevent.orders=[]
		newevent.uppercosts=[]
		newevent.lowercosts=[]
		myprevals=[]
		for h in dups:
			histilist=[i for i, j in enumerate(event.histories) if j==h]
			pvals=[(event.prevals[i], i) for i in histilist]
			pvals.sort(key=lambda x: x[0])
			(p, i) = pvals[0]
			newevent.histories.append(event.histories.pop(i))
			newevent.prevals.append(event.prevals.pop(i))
			newevent.orders.append(event.orders.pop(i))
			newevent.uppercosts.append(event.uppercosts.pop(i))
			newevent.lowercosts.append(event.lowercosts.pop(i))
			if len(event.id.split('.'))>1:
				(idhist, order) = map(int, event.id.split('.'))
				newevent.id = "%d.%d" % (idhist, min(newevent.orders))
		return [newevent]
	else: 
		return []

def unique_c_events_list(eventslist):
	sys.stderr.write("There are %d events before splitting cycles\n" % (len(eventslist))) 
	if Global_SPLITCYCLES: 
		eventslist=split_cycles(eventslist)
	eventslist=sorted(eventslist, key=lambda x: (x.segstr, x.cnval, x.prevals[0]))
	unique_events=[]
	eventA=eventslist[0]
	for eventB in eventslist[1:]:
		if eventA == eventB:
			eventA.add_Event_data(eventB)
		else: 
			unique_events.append(eventA)
			eventA=eventB
	unique_events.append(eventA)
	sys.stderr.write("There are %d events after merging\n" % len(unique_events))
	# check that the event doesn't need to be split again because it happens twice in some histories
	finalevents=unique_events
	splitoffs=get_split_offs(unique_events)
	sys.stderr.write("There are %d splitoffs\n" % len(splitoffs))
	if Global_SPLITOFFS: 
		finalevents+=splitoffs
	return finalevents 
	
def make_events_from_braneyfn(braneyfn): 
	sys.stderr.write("working on %s\n" % braneyfn)
	braneyf=gzip.open(braneyfn, 'rb')
	eventorder=-1
	histid=0
	myevents=[]
	eventlines=[]
	batchbreak=3000 # after x histories, stop and merge them together.  
	breakcounter=0	# use to keep track of the number of histories we've done
	histcount=0  #count the number of histories we've done, and what is used as the history id for an event. 
	eventcount=0
	peventcount=0
	for braneyline in braneyf: 
		if braneyline.strip() != '':
			braneyseg=Braney_seg(braneyline)
			if eventorder==-1: 
				eventorder=braneyseg.order #the order is unique within but not between histories
				histid=braneyseg.historyid
				eventlines.append(braneyseg)
			elif braneyseg.order == eventorder and braneyseg.historyid == histid:
				eventlines.append(braneyseg)
			elif braneyseg.order != eventorder or braneyseg.historyid != histid: 
				myevent=Event(eventlines)
				myevent.histories=[histcount]  #use histcount instead of the given history id because sometimes when cn-avg.py is run in -c mode histories are added to the end of .braney files, but the history ids start from 0 again. 
				myevent.make_segstr()
				myevents.append(myevent)
				eventcount+=1
				if histid != braneyseg.historyid: 
					histcount +=1
					breakcounter+=1
					peventcount=eventcount
					eventcount=0
				if breakcounter >= batchbreak:
					uniqueevents=unique_c_events_list(myevents)
					sys.stderr.write("%d\t%d\t%d\t%d\n" % (histid, peventcount, len(myevents), len(uniqueevents)))
					myevents=uniqueevents
					breakcounter=0
				eventlines=[braneyseg]
				eventorder=braneyseg.order
				histid=braneyseg.historyid
	# append the last event 
	myevent=Event(eventlines)
	myevent.histories=[histcount]  
	myevent.make_segstr()
	myevents.append(myevent)
	uniqueevents=unique_c_events_list(myevents)
	sys.stderr.write("num_events: %d, filtered: %d\n" % (len(myevents), len(uniqueevents)))
	return (uniqueevents)

def get_events_from_cnavgdir(cnavgdir, historyScores, totalp=0):
	allevents=[]
	braneyfiles=glob.glob(cnavgdir+"/"+"HISTORIES_*.braney")
	sys.stderr.write("braneyfiles: %s\n" % (str(braneyfiles)))
	for braneyfn in braneyfiles:
		sim=int(re.match(".*HISTORIES_(\d+)\.braney", braneyfn).group(1))
		events=make_events_from_braneyfn(braneyfn)
		pickle.dump(events, open("sim.%d.pevnts" % sim, 'w'), pickle.HIGHEST_PROTOCOL) 
		for evnt in events:
			for i in xrange(len(evnt.histories)):
				evnt.histories[i] = evnt.histories[i] + sim*Global_BINWIDTH
			(idhist, order) = map(int, evnt.id.split('.'))
			evnt.id = "%d.%d" % (idhist+(sim*Global_BINWIDTH), order)
		allevents=unique_c_events_list(allevents+events) 
	finalevents=allevents
	if totalp==0:
		# get the total likelihood of all events for computing marginal likelihoods
		totalp=compute_totalp(historyScores)
	for evnt in finalevents:
		evnt.update(historyScores)
		evnt.likelihood=compute_likelihood_histories(evnt.histories, historyScores, totalp)
		evnt.trim()
	return finalevents

def merge_pevnts_files(pevntsfiles, outputfile, historyScores, totalp): 
	allevents=pickle.load(open(pevntsfiles[0], 'rb'))
	for pevntsfile in pevntsfiles[1:]:
		addinevents=pickle.load(open(pevntsfile, 'rb'))
		allevents=unique_c_events_list(allevents+addinevents)
		sys.stderr.write("addinevents %d, allevents: %d\n" % (len(addinevents), len(allevents)))
	for evnt in allevents:
		evnt.update(historyScores)
		evnt.likelihood=compute_likelihood_histories(evnt.histories, historyScores, totalp)
		evnt.trim()
	#return(allevents)
	pickle.dump(allevents, open(outputfile, 'wb'), pickle.HIGHEST_PROTOCOL)

def combine_history_statsfiles(cnavgdir): 
	statsfiles=glob.glob(cnavgdir+"/"+"HISTORY_STATS*")
	sys.stderr.write("statsfiles: %s\n" % (str(statsfiles)))
	mystats=np.array([])
	mysims=[]
	runlens=[]
	for statsfile in statsfiles:
		sim=int(re.match(".*HISTORY_STATS_(\d+)", statsfile).group(1))
		mysims.append(sim)
		print "sim is %d" % sim
		historystats=np.loadtxt(statsfile, dtype=int) 
		runlens.append(historystats.shape[0])
	runlen=max(runlens)
	for sim in mysims:
		statsfile=os.path.join(cnavgdir, "HISTORY_STATS_%d" % sim)
		historystats=np.loadtxt(statsfile, dtype=int, ndmin=2) 
		sys.stderr.write("dim of historystats is %s\n" % str(historystats.shape))
		if mystats.size==0: 
			mystats=np.zeros(((max(mysims)+1)*runlen, historystats.shape[1]+1), dtype=int)
		hids=np.array(range(historystats.shape[0]))+ sim*Global_BINWIDTH
		i = sim*runlen
		mystats[i:(i+historystats.shape[0]),:] = np.hstack((np.atleast_2d(hids).T, historystats))
	return mystats

def get_historyScores(statsfile): 
	historystats=np.fromregex(statsfile, statsFileRegex, dtype=int)
	hids=np.array(range(historystats.shape[0]))
	historyScores=np.hstack((np.atleast_2d(hids).T, historystats))
	return historyScores

def compute_likelihood_histories(historyids, historyScores, denom=1):
	if len(historyids)>0: 
		historyids=np.unique(historyids)
		hindices = historyids_to_indices(historyids, historyScores)
		costsarray=np.mean(historyScores[hindices,1:3], axis=1)
		maskedcosts=np.ma.masked_where(costsarray==0, costsarray)
		x=np.sum(np.exp(-1*Global_K*maskedcosts))
		if denom==1: 
			denom=compute_totalp(historyScores)
		return x/denom
	else: 
		return 0

def compute_totalp(historyScores): 
	allh=historyScores[:,1]>0
	historyids = historyScores[allh,0]
	hindices = historyids_to_indices(historyids, historyScores)
	costsarray=np.mean(historyScores[hindices,1:3], axis=1)
	maskedcosts=np.ma.masked_where(costsarray==0, costsarray)
	x=np.sum(np.exp(-1*Global_K*maskedcosts))
	return x 

def historyids_to_indices(historyids, historyScores): 
	hids=np.array(historyids, dtype=int)
	itr=np.fmod(hids, Global_BINWIDTH)
	sim=np.round(hids/Global_BINWIDTH)
	runlen=max(np.fmod(historyScores[:,0], Global_BINWIDTH))+1
	newi=itr+sim*runlen
	return newi.astype('int')

def get_runlen(historyScores): 
	runlen=max(np.fmod(historyScores[:,0], Global_BINWIDTH))+1
	nruns=historyScores.shape[0]/runlen
	return(runlen, nruns)

def get_breakpoints(events, historyid): 
	breaklocs={}
	for event in events:
		event.unpack() 
		istrue=0
		ispred=1
		if historyid != "": 
			if historyid in event.histories: 
				istrue=1
				if len(event.histories)==1:
					ispred=0
		for seg in event.segs: 
			locs=("%s\t%d" % (seg.chr, seg.start), "%s\t%d" % (seg.chr2, seg.end))
			for loc in locs: 
				if loc not in breaklocs.keys(): 
					breaklocs[loc]=[ispred, istrue]
				else: 
					x=breaklocs[loc]
					x[0] = x[0]+ ispred
					x[1]= x[1]+istrue
	return breaklocs

def get_event_costs_over_time(events, historyScores, run):
	newevent.id = "%d.%d" % (newevent.histories[0], newevent.orders[0]) 
	return newevent

def merge_events_by_type(events, historyScores=None): 
	if Global_SPLITCYCLES: 
		events=split_cycles(events)
	 #order events by segstr
	sevents=sorted(events, key=lambda x: (x.segstr))
	unique_events=[]
	eventA=sevents[0]
	for eventB in sevents[1:]: 
		if eventB.segstr == eventA.segstr and ((eventB.cnval * eventA.cnval) >0) :  
			eventA.add_Event_data(eventB)
		else: 
			unique_events.append(eventA)
			eventA=eventB
	unique_events.append(eventA)
	finalevents=unique_events
	splitoffs=get_split_offs(finalevents)
	sys.stderr.write("There are %d splitoffs\n" % len(splitoffs))
	if Global_SPLITOFFS: 
		finalevents+=splitoffs
	for e in finalevents: 
		e.update(historyScores)
	return(finalevents)

#write an event in a different format
def write_dat_header(): 
	return "event_id	event_type	avecost	Lscore	CNval	length	numsegs	prevmean	prevsd	ordmean	ordsd	numhists\n"	

def write_event_dat(event, dat=False): 
	mystr="\t".join(map(str, [
	event.id,
	event.determineEventType(),
	np.mean(np.array(event.uppercosts + event.lowercosts)),
	event.likelihood,
	event.cnval,
	event.get_Event_length(),
	event.numsegs,
	event.prevalmean, event.prevalsd,
	event.ordermean, event.ordersd,
	event.numhists]
	)) + "\n"
	return mystr

