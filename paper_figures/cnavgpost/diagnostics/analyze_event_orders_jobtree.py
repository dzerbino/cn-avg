#!/usr/bin/env python 

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger
from sonLib.bioio import system
import os, sys, glob, re
import subprocess
import cPickle as pickle
import cnavgpost.mergehistories.event_cycles_module as histseg
import cnavgpost.diagnostics.get_event_order_counts as ordcnts
import cnavgpost.simulations.correct_simulation_events as correctevents 
import cnavgpost.cnavg_post_analysis_jobtree as cpaj
import numpy as np

class Setup(Target):
	def __init__(self, options):
		Target.__init__(self)
		self.options=options

	def run(self):
		opts=self.options
		self.logToMaster("setting up...")
		mylines=open(opts.samplelist, 'r').readlines()
		for line in mylines:
			(cnavgpost, sampleid) = line.strip().split('\t')[:2]
			self.logToMaster("working on %s in %s" % (sampleid, cnavgpost))
			statsfn=os.path.join(cnavgpost, "historystats.txt")
			outputdir=cnavgpost
			pevntsfile=os.path.join(outputdir, sampleid + ".pevnts")	
			self.addChildTarget(RunEventOrderAnalysis(sampleid, cnavgpost, statsfn, opts)) 

class RunGetShuffleHistoryRanges(Target):
	def __init__(self, sampleid, cnavgpost, statsfn, options): 
		Target.__init__(self)
		self.sampleid = sampleid
		self.cnavgpost=cnavgpost
		self.statsfn=statsfn
		self.opts=options		

	def run(self): 
		opts=self.opts
		histScores=np.loadtxt(self.statsfn, dtype='int')
		(cndat, shdat)=ordcnts.get_shuffle_history_scores(histScores, opts.shufstart, opts.simulation)
		np.savetxt(os.path.join(self.cnavgpost, "cnranges.txt"), cndat, fmt='%d', delimiter='\t')	
		np.savetxt(os.path.join(self.cnavgpost, "shranges.txt"), shdat, fmt='%d', delimiter='\t')	

class RunEventOrderAnalysis(Target): 
	def __init__(self, sampleid, cnavgpost, statsfn, options): 
		Target.__init__(self)
		self.sampleid=sampleid
		self.cnavgpost=cnavgpost
		self.statsfn=statsfn
		self.opts=options
	
	def run(self):
		opts=self.opts
		pvalcutoff=opts.cutoff 
		outdir=self.cnavgpost
		useEdges=opts.edges
		if useEdges: 
			pevntsfile=os.path.join(outdir, "%s.pedgs" % self.sampleid)
			outfn=os.path.join(outdir, "edges_ordcnts.dat") 
		else: 
			pevntsfile=os.path.join(outdir, "%s.pevnts" % self.sampleid)
			outfn=os.path.join(outdir, "events_ordcnts.dat")
		histScores=np.loadtxt(self.statsfn, dtype='int')
		# filter the events
		events=pickle.load(open(pevntsfile, 'rb'))
		myevents=[]
		for e in events:
			#e.histories=histseg.listout_ranges(e.histRanges)
			if e.likelihood > pvalcutoff:
					myevents.append(e)
		else:
			myevents=events
		self.logToMaster("There are %d events with pvalue gt %f\n" % (len(myevents), pvalcutoff))
		ordcnts.get_events_order_counts(myevents, outfn, opts.simulation, histScores, opts.shufstart)
		if False: #I'm still debugging this so it's not included
			if useEdges:
				outfn1=os.path.join(outdir, "edges_earlycnts.dat")
			else:
				outfn1=os.path.join(outdir, "event_earlycnts.dat")
			outfn2=os.path.join(outdir, "histlengths.dat")
			self.addChildTarget(ordcnts.count_earlylate_with_correction(myevents, historyScores, outfn1, outfn2))


class PrintEventCountsPerHistory(Target): 
	def __init__(self, sampleid, cnavgpost, opts): 
		Target.__init__(self)
		self.opts=opts
		self.sampleid=sampleid
		self.cnavgpost=cnavgpost

	def run(self): 
		opts=self.opts
		cnavgout="/bigdrive/intsims/cnavgout"
		braneyfiles=glob.glob(os.path.join(cnavgout, self.sampleid, "HISTORIES_*.braney"))
		self.logToMaster("%s, braneyfiles: %s\n" % (cnavgout, braneyfiles))
		sys.stderr.write("%s, braneyfiles: %s\n" % (cnavgout, braneyfiles))
		runlens=[]
		myarrays=[]
		mysims=[]
		for bfile in braneyfiles: 
			sim=int(re.match(".*/HISTORIES_(\d+).braney", bfile).group(1))
			mysims.append(sim)
			histlengths=subprocess.check_output("zcat %s | grep ^A | cut -f10,14 | sort -u | cut -f1 | sort -n | uniq -c | awk '{print $2\"\\n\"$1}'" % (bfile), shell=True)
			myarray=np.fromstring(histlengths, dtype=int, sep='\n')
			runlen=myarray.size/2
			myarray=np.reshape(myarray, (runlen, 2))
			myarray[:,0] = myarray[:,0] + sim*histseg.Global_BINWIDTH
			myarrays.append(myarray)
			runlens.append(runlen)
			self.logToMaster("runlens is %d, and runlen is %d\n" % (len(runlens), runlen))
		runlen=max(runlens)
		mystats=np.array([], dtype=int)
		for s in xrange(len(mysims)):
			sim=mysims[s]  
			histstats=myarrays[s]
			i=sim*runlen
			if mystats.size==0: 
				mystats=np.zeros(((max(mysims)+1)*runlen, histstats.shape[1]), dtype=int)
			mystats[i:(i+histstats.shape[0]),:]=histstats
		np.savetxt(os.path.join(self.cnavgpost, "historylengths.txt"), mystats, fmt='%d', delimiter="\t")
	

def main():
	parser = OptionParser(usage = "run get_event_order_counts on the list of samples in samplelist.") 
	parser.add_option("--samplelist", dest="samplelist", help="The list of CNAVG outputs and sample ids. Should have the form <directory><ID>.")
	parser.add_option("--edges", dest="edges", action='store_true', help="Do the analysis using the .pedges instead of .pevnts files.")
	parser.add_option('--simulation', dest="simulation", action='store_true', help='whether the history is a simulation')
	parser.add_option('--cutoff', dest="cutoff", help='only look at events with a likelihood above this cutoff', default=0, type=float)
	parser.add_option('--shufstart', dest="shufstart", default=0, type=int, help='time when order shuffling begins and cn-avg sampleing ends.')
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	i = Stack(Setup(options)).startJobTree(options)
	if i:
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)


if __name__=="__main__":
	from analyze_event_orders_jobtree import *
	main()

