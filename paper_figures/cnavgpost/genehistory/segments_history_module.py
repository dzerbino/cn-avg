# Module for merging different edges together into a set of CN changes for a segment of the genome. 
import numpy as np
import cnavgpost.mergehistories.event_cycles_module as ecycles 
import copy, sys

class EventSegment: 
	def __init__(self, event, chrm, start, end): 
		self.chr=chrm
		self.start=start
		self.end=end
		self.histories=event.histories
		self.orders=event.orders
		self.prevals=event.prevals
		self.cnval=event.cnval
	
	def __str__(self): 
		return "%s\t%d\t%d\t%d\t%d\n" % (self.chr, self.start, self.end, self.cnval, len(self.histories))
		

class SegmentHistory: 
	def __init__(self, eseg): 
		if isinstance(eseg, EventSegment):
			self.chr=eseg.chr
			self.start=eseg.start 
			self.end=eseg.end 
			self.esegs=[eseg]
			self.CNprofiles=[]
		elif isinstance(eseg, str): 
			dat=eseg.strip().split("\t")
			self.chr=dat[0]
			self.start=int(dat[1])
			self.end=int(dat[2])
			self.esegs=[]
			self.CNprofiles=[CNprofile(eseg)]


	def __str__(self):
		fstr=""
		if len(self.CNprofiles)>0:  
			for p in self.CNprofiles:  
				mystr="\t".join(map(str, [
            		self.chr, self.start, self.end,
            		p.likelihood,
            		p.numhists,
            		",".join(map(str, p.cnvals)),
            		",".join(map(str, p.pvals)),
            		",".join(map(str, p.pvalsd))
					])) + "\n"
				fstr+=mystr
		else: 
			fstr= "\t".join(map(str, [self.chr, self.start, self.end, len(self.esegs)])) + "\n" 
		return fstr

	def samelocus(self, other):
		return ((self.chr == other.chr) and (self.start == other.start) and (self.end==other.end))

	def overlaps(self, other): 
		return (self.chr == other.chr and self.start <= other.end and self.end >= other.start)

	def comes_before(self, other):
		return ((self.chr < other.chr) or (self.chr == other.chr and self.end < other.start))

	def comes_after(self, other):
		return ((self.chr > other.chr) or (self.chr == other.chr and self.start > other.end))

	def appendvals(self, eseg): 
		self.esegs.append(eseg)

class CNprofile: 
	def __init__(self, line=""):
		if line=="": 
			self.likelihood=0
			self.numhists=0
			self.pvals=[]
			self.pvalsd=[]
			self.cnvals=[]
		else: 
			dat=line.strip().split("\t")
			self.likelihood=float(dat[3])
			self.numhists=int(dat[4])
			self.cnvals=map(float, dat[5].split(','))	
			self.pvals=map(float, dat[6].split(','))	
			self.pvalsd=map(float, dat[7].split(','))	
			

def create_CNprofiles_from_Edges(self, histScores, totalp=0):
	edgelist=self.esegs 
	mycnprofs=[]
	if totalp==0: totalp=ecycles.compute_totalp(historyScores)
	(hprofiles, pprofiles) = create_profile_matrices(edgelist, histScores)
	(ridx, profiles) = get_unique_rows(hprofiles)
	goodi=np.where(np.sum(profiles != 0, axis=1) > 0)[0]
	for i in goodi:
		mycnp=CNprofile() 
		mycnp.likelihood=ecycles.compute_likelihood_histories(histScores[ridx==i,0], histScores, totalp)
		mycnp.numhists=np.sum(ridx==i) 
		cnvals=profiles[i,:]
		pvals = pprofiles[ridx==i,:]
		pvalsd=np.std(pvals, axis=0)  
		pvalsm=np.mean(pvals, axis=0)
		mycnp.pvals=pvalsm[cnvals!=0].tolist()
		mycnp.pvalsd=pvalsd[cnvals!=0].tolist()
		mycnp.cnvals=cnvals[cnvals!=0].tolist()
		mycnprofs.append(mycnp)
	self.CNprofiles=mycnprofs
#	self.edges=[]

def create_profile_matrices(edgelist, histScores): 
	hprofiles=np.zeros((histScores.shape[0], len(edgelist)))
	pprofiles=np.zeros((histScores.shape[0], len(edgelist)))
	for i,e in enumerate(edgelist): 
		hi=ecycles.historyids_to_indices(e.histories, histScores)
		hprofiles[hi,i]=e.cnval
		pprofiles[hi,i]=np.array(e.prevals)	
	# get equivalent profiles after removing 0s
	order_by_prevals(hprofiles, pprofiles)
	pprofiles=collapse_zeros(pprofiles)
	hprofiles=collapse_zeros(hprofiles)
	return(hprofiles, pprofiles)


def order_by_prevals(hprofiles, pprofiles): 
	sorti=np.argsort(pprofiles, axis=1)
	for i in xrange(hprofiles.shape[0]): 
		myi=sorti[i,:][::-1]
		pprofiles[i,:]=pprofiles[i,myi]
		hprofiles[i,:]=hprofiles[i,myi]
#	return(hprofiles, pprofiles) 

def collapse_zeros(hprofiles): 
	maxrows=max(np.sum(hprofiles!=0, axis=1))
	newm=np.zeros((maxrows, hprofiles.shape[1]))
	(ridx, uprofiles) = get_unique_rows(hprofiles != 0)
	for p in xrange(uprofiles.shape[0]):
		subarray=hprofiles[ridx==p,:]
		profile=uprofiles[p,:]
		for col in xrange(subarray.shape[1]):
			if profile[col] == False and (col+1) < profile.shape[0]:
				subarray[:,col:(subarray.shape[1]-1)]=subarray[:,(col+1):]
				subarray[:,-1]=0
		hprofiles[ridx==p,:]=subarray
	newm=hprofiles[:,:maxrows]
	return newm

def get_unique_rows(hprofiles):
	a=hprofiles
	b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
	u, uidx, ridx = np.unique(b, return_index=True, return_inverse=True)
	uprofiles=hprofiles[uidx]
	# get rid of the nothing profiles where the CN doesn't change at all. 
	return(ridx, uprofiles)

# given an edge and a segment history, it will merge the edge into the seghist, splitting the seghist if necessary, but will not return a segment outside the boundaries of the original segment. 
def merge_in_edge(edge, seghist):
	new_segments=[]
	if seghist.overlaps(edge):
		if (seghist.start >= edge.start and seghist.end <= edge.end): 
			seghist.appendvals(edge)
			return [seghist]
		elif (seghist.start >= edge.start) and (seghist.end > edge.end): 
			newseg=copy.deepcopy(seghist)
			newseg.end=edge.end
			newseg.appendvals(edge)
			seghist.start=edge.end+1
			new_segments.append(newseg)
			new_segments.append(seghist)
		elif (seghist.start < edge.start) and (seghist.end <= edge.end): 
			newseg=copy.deepcopy(seghist)
			newseg.end=edge.start-1
			new_segments.append(newseg)
			seghist.appendvals(edge)
			seghist.start=edge.start
			new_segments.append(seghist)
		elif (edge.start > seghist.start) and (seghist.end > edge.end): 
			newseg=copy.deepcopy(seghist)
			newseg.end=edge.start-1
			new_segments.append(newseg)
			newseg=copy.deepcopy(seghist)
			newseg.start=edge.end+1
			seghist.appendvals(edge)
			seghist.start=edge.start
			seghist.end=edge.end
			new_segments.append(seghist)
			new_segments.append(newseg)
	return new_segments


# This will create segments for the largest region, including the edge even if the edge expands outside the boundary of the segment history. 
def incorporate_edge(edge, seghist): 
	new_segments=[]
	if seghist.overlaps(edge):
		if (seghist.start == edge.start and seghist.end == edge.end): 
			seghist.appendvals(edge)
			return [seghist]
		elif seghist.start > edge.start:
			newseg=SegmentHistory(edge)
			newseg.end=seghist.start-1
			new_segments.append(newseg)
			if seghist.end <= edge.end: 
				seghist.appendvals(edge)
				new_segments.append(seghist)
				if seghist.end < edge.end: 
					newseg=SegmentHistory(edge)
					newseg.start=seghist.end+1
					new_segments.append(newseg)
			else: # seghist.end > edge.end: 
				newseg=copy.deepcopy(seghist)
				newseg.end=edge.end
				new_segments.append(newseg)
				seghist.appendvals(edge)
				seghist.start=edge.end+1
				new_segments.append(seghist)	
				
		elif seghist.start < edge.start:
			newseg=copy.deepcopy(seghist)
			newseg.end=edge.start-1
			new_segments.append(newseg)
			seghist.start=edge.start
			if seghist.end <= edge.end:
				seghist.appendvals(edge)
				new_segments.append(seghist)
				if seghist.end < edge.end: 
					newseg=SegmentHistory(edge)
					newseg.start=seghist.end+1
					new_segments.append(newseg)
			else: # seghist.end > edge.end: 
				newseg=copy.deepcopy(seghist)
				newseg.end=edge.end
				new_segments.append(newseg)
				seghist.appendvals(edge)
				seghist.start=edge.end+1
				new_segments.append(seghist)
		elif seghist.start == edge.start:
			if seghist.end < edge.end: 
				seghist.appendvals(edge)
				new_segments.append(seghist)
				newseg=SegmentHistory(edge)
				newseg.start=seghist.end+1
				new_segments.append(newseg)
			else: # seghist.end > edge.end: 
				newseg=copy.deepcopy(seghist)
				newseg.end=edge.end
				new_segments.append(newseg)
				seghist.appendvals(edge)
				seghist.start=edge.end+1
	return(new_segments)


# The Segment file should be ordered by genome coordinate so that the same coordinates are sequential.  Also, they won't overlap, although that doesn't matter for this. 
def read_in_segs(segfn): 
	mysegs=[]
	currseg=None
	for l in open(segfn): 
		segh=SegmentHistory(l)
		if currseg is not None:
			if segh.samelocus(currseg):
				currseg.CNprofiles += segh.CNprofiles
			else: 
				mysegs.append(currseg)
				currseg=segh
		else: 
			currseg=segh
	mysegs.append(currseg)
	return mysegs



	
