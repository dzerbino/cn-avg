# Module for reading .braney files, which are the text output from the cn-avg.py 
# The .braney files contain multiple histories sampled over many iterations, typically.  

import sys, os
import re
import subprocess
import math

class Braney_seg:
	def __init__(self, myline):
		data=myline.strip().split('\t')
		if myline[0] != 'A': # check that it is not an adjacency line
			self.adj=False
			self.seg=True
			self.chr=data[0]
			self.start=int(data[1])
			self.st1=None
			self.chr2=self.chr
			self.end=int(data[2])
			self.st2=None
			self.preval=float(data[4])
			self.cnval=round(float(data[3])/self.preval)
			self.historyid=int(data[5])
			self.cycleorder=int(data[8])
			self.order=int(data[9])
			self.upperHistCost=int(data[10]) 
			self.lowerHistCost=int(data[11])
			self.upperEventCost=int(data[12]) 
			self.lowerEventCost=int(data[13]) 
			self.ptrid=data[14]
		else: # an adjacency line
			self.adj=True
			self.seg=False
			self.chr=data[1]
			self.start=int(data[2])
			self.st1=data[3]
			self.chr2=data[4]
			self.end=int(data[5])
			self.st2=data[6]
			self.preval=float(data[8])
			self.cnval=round(float(data[7])/self.preval)
			self.historyid=int(data[9])
			self.cycleorder=int(data[12])
			self.order=int(data[13])		
			self.upperHistCost=int(data[14]) 
			self.lowerHistCost=int(data[15]) 
			self.upperEventCost=int(data[16]) 
			self.lowerEventCost=int(data[17]) 
			self.ptrid=data[18]
			self.order=float(data[13])		
	#	self.order_ends()

	def __str__(self):
		if self.seg:
			return "%s\t%d\t%d\t%f\t%f\t%d\t0\t0\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" % (self.chr, self.start, self.end, self.cnval*self.preval, self.preval, self.historyid, self.cycleorder, self.order, self.upperHistCost, self.lowerHistCost, self.upperEventCost, self.lowerEventCost, self.ptrid)
		else: 
			return "A\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%f\t%d\t0\t0\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" % (self.chr, self.start, self.st1, self.chr2, self.end, self.st2, self.cnval * self.preval, self.preval, self.historyid, self.cycleorder, self.order, self.upperHistCost, self.lowerHistCost, self.upperEventCost, self.lowerEventCost, self.ptrid)

	def __eq__(self, other):
		vals_to_check=['adj', 'chr','chr2', 'start', 'end', 'st1', 'st2', 'cnval']
		self_dict=self.__dict__
		other_dict=other.__dict__
		equal=True
		for kval in vals_to_check: 
			if self_dict[kval]!= other_dict[kval]:
				equal=False
		return equal

	def same_coords(self, other): 
		vals_to_check=['adj', 'chr','chr2', 'start', 'end', 'st1', 'st2']
		self_dict=self.__dict__
		other_dict=other.__dict__
		equal=True
		for kval in vals_to_check: 
			if self_dict[kval]!= other_dict[kval]:
				equal=False
		return equal

	def is_dup(self, other): #here we are assuming we know that the segments are in the same event and history. id self.historyid == other.historyid and self.order=other.order 
		return (self.cycleorder == other.cycleorder)  

	def is_dup_OLD(self, other): 
		if (
			((self.adj and other.adj) and 
			 (self.cnval == other.cnval) and 
			 (self.preval==other.preval) and 
			 (self.cycleorder == other.cycleorder))
			and 
			(((self.chr == other.chr2) and 
			  (self.start == other.end) and 
			  (self.end == other.start) and
			  (self.st1 == other.st2) and 
			  (self.st2== other.st1)) 
				or 
			 	((self.chr == other.chr) and 
			 	(self.chr2==other.chr2) and 
				(self.end==other.end) and 
				(self.start == other.start) and 
				(self.st1 == other.st1) and 
				(self.st2 == other.st2)
				)
				)
			): 
			return True
		else: 
			return False
	
	def order_ends(self):
		if self.chr > self.chr2: 
			(x, y, z) = (self.chr, self.start, self.st1)
			(self.chr, self.start, self.st1) =(self.chr2, self.end, self.st2)
			(self.chr2, self.end, self.st2) =(x, y, z)
		elif (self.chr == self.chr2 and self.start > self.end): 
			(y,z)=(self.start, self.st1)
			(self.start, self.st1) =(self.end, self.st2)
			(self.end, self.st2)=(y,z)
	
	def flip_ends(self): 
		(x, y, z) = (self.chr, self.start, self.st1)
		(self.chr, self.start, self.st1) =(self.chr2, self.end, self.st2)
		(self.chr2, self.end, self.st2) =(x, y, z)
		

	def adjacency_cross(self, other):
		self.order_ends()
		other.order_ends() 
		if ((self.chr == self.chr2 == other.chr == other.chr2) and 
			(((self.start > other.start and self.start < other.end and self.end > other.end)) or 
			((self.start < other.start and self.end < other.end and self.end > other.start)))): 
			return True
		else: return False

	def isbreakpoint(self): 
		self.order_ends()
		return  (self.adj and (self.start == (self.end-1)) and (self.chr == self.chr2))
	
	def connected_to(self, other): 
		if ((self.chr == other.chr and self.start == other.start) or  
			(self.chr == other.chr2 and self.start == other.end) or 
			(self.chr2 == other.chr and self.end == other.start) or 
			(self.chr2 == other.chr2 and self.end == other.end)): 
			return True
		else: 
			return False

def compute_likelihood(costs, totp): 
	mysum=0
	for c in costs: 
		mysum+=math.exp(-c)
	likelihood=mysum/totp
	return likelihood

def make_tabix_from_braney(braneyfn, outdir):
	newfn="%s/%s.gz" % (outdir, os.path.basename(os.path.normpath(braneyfn)))
	if not os.path.isfile(newfn):
		subprocess.call(["gunzip -dc %s | grep -v ^A | sort -k1,1 -k2,2n | awk '$3>$2' | bgzip > %s" % (braneyfn, newfn)], shell=True)
	if not os.path.isfile(newfn+".tbi"):
		subprocess.call(["tabix -0 -s 1 -b 2 -e 3 %s" % (newfn)], shell=True)
	newfna="%s/%s-adj.gz" % (outdir, os.path.basename(os.path.normpath(braneyfn)))
	if not os.path.isfile(newfna):
		subprocess.call(["gunzip -dc %s | grep ^A | sort -k2,2 -k3,3n | bgzip > %s" % (braneyfn, newfna)], shell=True)
	if not os.path.isfile(newfna+".tbi"):
		subprocess.call(["tabix -0 -s 2 -b 3 -e 3 %s" % (newfna)], shell=True)
	return (newfn, newfna)

def get_total_history_prob(braneyfn):
	costs = subprocess.check_output(["gunzip -dc %s | grep -v ^A | grep -v '^$' | cut -f6,11 | uniq | cut -f2" % (braneyfn)], shell=True)
	mysum=0
	for cost in costs.strip().split("\n"):
		c=float(cost)
		mysum+= math.exp(-c)
	return mysum

def same_endpoints(seg1, seg2): 
	return ((seg1.chr == seg2.chr and seg1.start ==seg2.start and 
		seg1.chr2==seg2.chr2 and seg1.end ==seg2.end) or 
		(seg1.chr == seg2.chr2 and seg1.start ==seg2.end and 
		seg1.chr2==seg2.chr and seg1.end ==seg2.start)) 
	
def get_direction(firstseg, second, last): 
	end_matches_second=False
	end_matches_last=False
	start_matches_second=False
	start_matches_last=False
	direction=0
	if ((firstseg.chr2 == second.chr and firstseg.end == second.start) 
		or (firstseg.chr2 == second.chr2 and firstseg.end == second.end)):
		end_matches_second=True
	if ((firstseg.chr2 == last.chr and firstseg.end == last.start)
		or (firstseg.chr2 == last.chr2 and firstseg.end == last.end)):
		end_matches_last=True
	if ((firstseg.chr == second.chr and firstseg.start == second.start)
		or (firstseg.chr == second.chr2 and firstseg.start == second.end)):
		start_matches_second=True
	if ((firstseg.chr == last.chr and firstseg.start == last.start)
		or (firstseg.chr == last.chr2 and firstseg.start == last.end)):
		start_matches_last=True
	if end_matches_second and start_matches_last: 
		direction=1
	elif end_matches_last and start_matches_second: 
		direction=-1
	return direction

def determine_4adj_cycle_type(segs, intra): # should be a list of 4 adjacencies
	for s in segs: 
		s.order_ends()
	(adj1, adj2, adj3, adj4) = segs[:4]
	if intra: # see if it's an inversion, reverse inversion, fusion or fission
		if ((adj1.adjacency_cross(adj3) and adj1.cnval >0) or 
			(adj2.adjacency_cross(adj4) and adj2.cnval>0)):
			type='inv'
		elif ((adj1.adjacency_cross(adj3) and adj1.cnval <0) or 
			(adj2.adjacency_cross(adj4) and adj2.cnval<0)):
			type='rinv'
		elif ((adj1.isbreakpoint() and adj3.isbreakpoint() and adj1.cnval<0) or 
			(adj2.isbreakpoint() and adj4.isbreakpoint() and adj2.cnval<0)): 
			type='fiss'
		elif ((adj1.isbreakpoint() and adj3.isbreakpoint() and adj1.cnval>0) or 
			(adj2.isbreakpoint() and adj4.isbreakpoint() and adj2.cnval>0)): 
			type='fuse'
		else: 
			type='oth'
	else: # it's an interchromosomal set of segments
		if ((adj1.isbreakpoint() and adj3.isbreakpoint() and adj1.cnval <0) or 			
			(adj2.isbreakpoint() and adj4.isbreakpoint() and adj2.cnval <0)): 
			type='tran'
		elif ((adj1.isbreakpoint() and adj3.isbreakpoint() and adj1.cnval >0) or 			
			(adj2.isbreakpoint() and adj4.isbreakpoint() and adj2.cnval >0)): 
			type='rtran'
		else: 
			type='oth'
	return type

