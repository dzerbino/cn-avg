# Module for linking events or segments, etc. according to order and co-occurence within histories 

import sys, os 
from cnavgpost.mergehistories.braney_lines_module import compute_likelihood

class Event_link: 
	def __init__(self, hist_objects): 
		if isinstance(hist_objects, str):
			data=myline.strip().split('\t')
			self.node1=data[0]
			self.node2=data[1]
			self.label1=data[4]
			self.label2=data[5]
			self.numhists=int(data[2])
			self.likelihood=float(data[3])
			#self.costs=[]
			self.histories=[]
		elif isinstance(hist_objects, tuple): # NOTE: the intersecthi are the indices of the second object that are in the first (so the indices should be with respect to object b)
			node1=hist_objects[0]
			node2=hist_objects[1]
			self.node1=node1.id
			self.node2=node2.id
			self.label1=str(node1.prevalmean)
			self.label2=str(node2.prevalmean)
			intersecthi=hist_objects[2]
			self.numhists=len(intersecthi) 
			if self.numhists >0:
				#self.costs=[node2.costs[i] for i in intersecthi]
				self.histories=[node2.histories[i] for i in intersecthi]
			else: 
				#self.costs=[]
				self.histories=[]
			self.likelihood=0

	def __str__(self): 
		mystr="%s\t%s\t%d\t%s\n" % (str(self.node1), str(self.node2), self.numhists, str(self.likelihood))
		return mystr
	
	def switch_nodes(self):
		x=self.node1
		self.node1=self.node2
		self.node2=x
		x=self.label1
		self.label1=self.label2
		self.label2=x

def link_events_by_order_within_histories(events): 
	#sort the events by prevalence first
	eventlist=sorted(events, key=lambda x: (x.prevalmean))	
	mylinks=[]
	ai=0
	while ai < len(eventlist): 
		eventa=eventlist[ai]
		ahistoriesleft=eventa.histories
		bi = ai+1
		while (bi < len(eventlist)) and (len(ahistoriesleft) >0): 
			eventb=eventlist[bi]
			intersecthi=get_index_of_intersecting_items(eventb.histories, list(ahistoriesleft))
			if len(intersecthi) >0: 
				mylink=Event_link((eventa, eventb, intersecthi))
				ahistoriesleft=set(ahistoriesleft) - set(mylink.histories)
				mylinks.append(mylink)
			bi=bi+1	
		ai=ai+1
	return mylinks

def get_index_of_intersecting_items(list1, list2): 
	myis=[]
	for i in xrange(len(list1)): 
		x=list1[i]
		if x in list2: 
			myis.append(i)
	return myis

	
