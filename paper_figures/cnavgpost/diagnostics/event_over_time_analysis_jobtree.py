#!/usr/bin/env python 

''' 
Analysis on CN-AVG output to test how well histories are being sampled across runs 
'''
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger
import os, sys, glob
import cPickle as pickle 
import numpy as np
import cnavgpost.mergehistories.event_cycles_module as histseg
import cnavgpost.diagnostics.likelihood_over_time_analysis as lota

class Setup(Target): 
	def __init__(self, options):
		Target.__init__(self)
		self.options=options

	def run(self):
		opts=self.options
		self.logToMaster("setting up...")
		histseg.Global_K=0
		mylines=open(opts.samplelist, 'r').readlines()
		for line in mylines:
			(cnavgout, sampleid) = line.strip().split('\t')[:2]
			self.logToMaster("working on %s in %s" % (sampleid, cnavgout))
			if opts.edges:
				subdir="lota_edge" 
				outputdir=os.path.join(cnavgout, "lota_edge")
			else: 	
				subdir="lotadir" 
				outputdir=os.path.join(cnavgout, "lotadir")
			if not os.path.exists(outputdir):
				os.makedirs(outputdir)
				self.addChildTarget(SetupLota(cnavgout, sampleid, outputdir, opts))
		outputdir=opts.outputdir
		self.setFollowOnTarget(CombineEventCounts(outputdir, mylines, subdir))

class SetupLota(Target):
	def __init__(self,cnavgout, sampleid,  outputdir, options): 
		Target.__init__(self) 
		self.options=options
		self.cnavgout=cnavgout
		self.sampleid=sampleid
		self.outputdir=outputdir

	def run(self):
		opts=self.options
		cnavgout=self.cnavgout
		sampleid=self.sampleid
		outputdir=self.outputdir 
		if opts.edges: 
			pevntsfile=os.path.join(cnavgout, "%s.pedgs" % sampleid)
		else: 
			pevntsfile=os.path.join(cnavgout, "%s.pevnts" % sampleid)
		events=pickle.load(open(pevntsfile, 'rb'))
		#filter out the subset of events that you want to look at - the True ones in this case
		TPevents=[]
		FPevents=[]
		Pcount=0
		for event in events: 
			if opts.simulation:
				if 0 in event.histories:
					Pcount+=1
					if len(event.histories)>1:
						TPevents.append(event)
				else: 
					FPevents.append(event)
			else: 
				Pcount+=1
				TPevents.append(event)
#		pickle.dump(truevents, os.path.join(outputdir, "true.pevnts"), pickle.HIGHEST_PROTOCOL) 
		historystats=os.path.join(cnavgout, "historystats.txt")
		historyScores=np.loadtxt(historystats, dtype=int)
		lota.likelihood_over_time_analysis(TPevents, historyScores, outputdir, "TP", 0.5, opts.simulation) 
		if opts.simulation: 
			lota.likelihood_over_time_analysis(FPevents, historyScores, outputdir, "FP", 0.5, opts.simulation) 
		lota.likelihood_over_time_analysis(events, historyScores, outputdir, "Tot", 0, opts.simulation) 
		myfh=open(os.path.join(outputdir, "true_count.txt"), 'w')
		myfh.write("%d\n" % Pcount)
		myfh.close()

class CombineEventCounts(Target):
	def __init__(self, outputdir, samples, subdir): 
		Target.__init__(self)
		self.outputdir=outputdir
		self.samples=samples
		self.subdir=subdir
	
	def run(self):
		mydata={}
		#mydata={"event_counts_ks": np.array([]), "nonzero_counts_ks": np.array([]), "zero_counts_ks": np.array([]), "event_counts_ns": np.array([]), "nonzero_counts_ns": np.array([]), "zero_counts_ns": np.array([])}
		totPcount=0
		mylines=self.samples
		self.logToMaster("running CombineEventCounts")
		for line in mylines:
			(cnavgout, sampleid) = line.strip().split('\t')[:2]
			datfiles=glob.glob(os.path.join(cnavgout, self.subdir, "*s.dat"))
			for datfile in datfiles:
				dat=np.loadtxt(datfile, dtype=int)
				datkey=os.path.basename(datfile).split(".dat")[0]
				if datkey in mydata.keys():
					currlen=mydata[datkey].shape[0]
					if (currlen > dat.shape[0]):
						nanarray=np.empty((currlen-dat.shape[0], dat.shape[1]))
						nanarray[:]=np.NAN
						dat=np.vstack((dat, nanarray))
					elif (currlen < dat.shape[0]):
						nanarray=np.empty((dat.shape[0]-currlen, dat.shape[1]))
						nanarray[:]=np.NAN
						mydata[datkey]=np.vstack((mydata[datkey], nanarray))
					mydata[datkey]=mydata[datkey]+dat
				else: 
					mydata[datkey]=dat
			Pcountfh=open(os.path.join(cnavgout, self.subdir, "true_count.txt"))
			Pcount=int(Pcountfh.readline().strip())
			totPcount+=Pcount
		myfh=open(os.path.join(self.outputdir, "Pcount.sum.txt"), 'w')
		myfh.write("%d\n" % totPcount)
		myfh.close()
		for datkey in mydata.keys():
			outputf=os.path.join(self.outputdir, "%s.sum.dat" % datkey)
			self.logToMaster("outputf is %s" % outputf)
			np.savetxt(outputf, mydata[datkey].astype(int), delimiter="\t", fmt="%d")
		# combine the lscores_*.txt files by cating them together. 
		if (True): 
			mydata={}
			for line in mylines:
				(cnavgout, sampleid) = line.strip().split('\t')[:2]
				lscorefiles=glob.glob(os.path.join(cnavgout, self.subdir, "lscores_*.txt.gz"))
				for datfile in lscorefiles: 
					dat=np.loadtxt(datfile)
					datkey=os.path.basename(datfile).split(".txt")[0]
					if datkey in mydata.keys(): 
						mydata[datkey]=np.hstack((mydata[datkey], dat))
					else: 
						mydata[datkey]=dat
			for datkey in mydata.keys(): 
				outputf=os.path.join(self.outputdir, "%s.sum.txt.gz" % datkey)
				self.logToMaster("outputf is %s" % outputf)
				np.savetxt(outputf, mydata[datkey], delimiter="\t")
			


def main(): 
	parser = OptionParser(usage = "event_over_time_analysis_jobtree.py ...")
	parser.add_option("--samplelist", dest="samplelist", help="The list of CNAVG outputs and sample ids. Should have the form <directory><ID>...")
	parser.add_option("--outputdir", dest="outputdir", help='Where you want the summary files to be put.  This directory should exist.')
	parser.add_option("--binwidth", dest="binwidth", help='the multiplier between history ids of independent runs', default=histseg.Global_BINWIDTH, type="int")
	parser.add_option('--simulation', dest="simulation", default=False, action="store_true", help="do simulation analysis.")
	parser.add_option('--edges', dest="edges", default=False, action="store_true", help="do analysis on *.pedgs files instead of *.pevnts files.")
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	histseg.Global_BINWIDTH=options.binwidth
	i = Stack(Setup(options)).startJobTree(options)
	if i: 
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == "__main__":
	from event_over_time_analysis_jobtree import *
	main()

