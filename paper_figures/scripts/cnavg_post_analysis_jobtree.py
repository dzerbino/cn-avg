#!/usr/bin/env python

# This uses jobTree (https://github.com/benedictpaten/jobTree) to run analysis on the CNAVG output.  

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger
from sonLib.bioio import system

import subprocess, os
import argparse
import numpy as np
import pysam 
import re, glob
import cPickle as pickle 

from cnavgpost.mergehistories.score_edges_within_pevents import *
from cnavgpost.mergehistories.score_and_link_cycles import *
import cnavgpost.mergehistories.event_cycles_module as histseg
import cnavgpost.mergehistories.create_pevnts_file_jobtree as pevntsjobtree
import cnavgpost.genehistory.pick_and_label_best_seghists as pick_and_label_best_seghists 
import cnavgpost.genehistory.make_seghists_from_edges as make_seghists_from_edges
import cnavgpost.genehistory.seghists_to_gene_orders as seghists_to_gene_orders 
import cnavgpost.diagnostics.mcmc_mixing_analysis_jobtree as mcmcjobtree
import cnavgpost.simulations.analyze_simulation as analyze_simulation
import cnavgpost.genehistory.annotate_events as annotate_events

#======= Main Setup ===================
class Setup(Target):
	def __init__(self, options): 
		Target.__init__(self)
		self.options=options
		self.events=[]
		self.totalp=0
		self.edges=[]
		self.historyScores=[]	
	
	def run(self):
		self.logToMaster("Setting up...") 
		opts=self.options
		histseg.Global_BINWIDTH=opts.binwidth
		sampleid=opts.sampleid
		outputdir=opts.outputdir
		subprocess.call("mkdir -p %s" % outputdir, shell=True)

		historystatsfile=os.path.join(outputdir, "historystats.txt")
		if not os.path.exists(historystatsfile): 
			self.logToMaster("Creating historystats.txt...%s" % historystatsfile) 
			logger.info("historystatsfile: %s" % historystatsfile)
			pevntsjobtree.CombineHistoryStatsfiles(opts, historystatsfile).run()
		self.historyScores=np.loadtxt(historystatsfile, dtype=int)
		self.totalp=histseg.compute_likelihood_histories(self.historyScores[:,0], self.historyScores)
		logger.info("Global_BINWIDTH: %d" % histseg.Global_BINWIDTH)
		logger.info("totalp is %s" % str(self.totalp))	

		pevntsfile=os.path.join(outputdir, opts.sampleid + ".pevnts")
		if opts.pevnts or not os.path.exists(pevntsfile): 
			self.logToMaster("Creating pevntsfile...%s" % pevntsfile) 
			logger.info("pevntsfile: %s" % pevntsfile)
			self.addChildTarget(pevntsjobtree.CreatePevntsFile(pevntsfile, self.historyScores, self.totalp, opts))
		self.setFollowOnTarget(DoAnalysisOfMergedEvents(opts))


class DoAnalysisOfMergedEvents(Target):
	def __init__(self, options): 
		Target.__init__(self)
		self.options=options
		self.events=[]
		self.totalp=0
		self.edges=[]
		self.historyScores=[]	
	
	def run(self):
		self.logToMaster("Setting up...") 
		opts=self.options
		histseg.Global_BINWIDTH=opts.binwidth
		sampleid=opts.sampleid
		outputdir=opts.outputdir
	
		historystatsfile=os.path.join(outputdir, "historystats.txt")
		if not os.path.exists(historystatsfile): 
			self.logToMaster("Creating historystats.txt...%s" % historystatsfile) 
			logger.info("historystatsfile: %s" % historystatsfile)
			pevntsjobtree.CombineHistoryStatsfiles(opts, historystatsfile).run()
		self.historyScores=np.loadtxt(historystatsfile, dtype=int)
		self.totalp=histseg.compute_likelihood_histories(self.historyScores[:,0], self.historyScores)
		#check that the *.pevnts file exists. 
		pevntsfile=os.path.join(outputdir, opts.sampleid + ".pevnts")
		if not os.path.exists(pevntsfile): 
			sys.exit("The required %s file does not exist." % pevntsfile)
	
		pedgesfile=os.path.join(outputdir, sampleid + ".pedgs")
		if opts.pedges or not os.path.exists(pedgesfile):
			self.logToMaster("Creating pedgesfile...%s" % pedgesfile) 
			logger.info("pedgesfile: %s" % pedgesfile)
			CreatePedgesFile(pickle.load(open(pevntsfile, 'rb')), pedgesfile, self.historyScores, self.totalp, False).run()
		
		seghistfile=os.path.join(outputdir, "seghists.txt")
		if opts.sgh or not os.path.exists(seghistfile): 
			self.logToMaster("Creating seghists file ... %s" % seghistfile)
			make_seghists_from_edges.main(pickle.load(open(pedgesfile, 'rb')), self.historyScores, seghistfile)
			
		# label the seghists if an annotation file is given
		if opts.bedfile: 
			labeledfn=os.path.join(outputdir, "seghists.labeled")
			if not os.path.exists(labeledfn):
				pick_and_label_best_seghists.main(seghistfile, opts.bedfile, True, labeledfn)
			geneordfn=os.path.join(outputdir, "geneords.txt")
			if opts.geneords or not os.path.exists(geneordfn):
				seghists_to_gene_orders.main(seghistfile, opts.bedfile, geneordfn)

		mrgpeventsfile=os.path.join(outputdir, sampleid + ".pmevnts")
		if not os.path.exists(mrgpeventsfile):
			self.logToMaster("Creating mpevnts...%s" % mrgpeventsfile) 
			logger.info("mrgpeventsfile: %s" % mrgpeventsfile)
			CreateMergedEventsFile(pickle.load(open(pevntsfile, 'rb')), mrgpeventsfile, self.historyScores).run()
		
		mrgpedgesfile=os.path.join(outputdir, sampleid + ".pmedgs")
		if not os.path.exists(mrgpedgesfile):
			self.logToMaster("Creating mrgpegesfile...%s" % mrgpedgesfile) 
			logger.info("mrgpedgesfile: %s" % mrgpedgesfile)
			CreateMergedEventsFile(pickle.load(open(pedgesfile, 'rb')), mrgpedgesfile, self.historyScores).run()
	
		breaksfile=os.path.join(outputdir, "breakpoints.txt")
		if not os.path.exists(breaksfile): 
			self.logToMaster("Creating breaksfile...%s" % breaksfile) 
			breaklocs=histseg.get_breakpoints(pickle.load(open(pedgesfile, 'rb')), opts.trueID)
			breaklocs2=histseg.get_breakpoints(pickle.load(open(mrgpedgesfile, 'rb')), opts.trueID)
			breaksfh=open(breaksfile, 'w')
			for loc in sorted(breaklocs.keys()):
				(n, t) = breaklocs[loc]
				(n2, t2) = breaklocs2[loc]
				breaksfh.write("%s\t%d\t%d\t%d\t%d\n" % (loc, n, t, n2, t2))	
		
		# Creating links is no longer an option.  
		#linksfile =os.path.join(outputdir, sampleid +".links")
		#if opts.links and not os.path.exists(linksfile): 
		#	self.logToMaster("Creating linksfile...%s" % linksfile) 
		#	logger.info("linksfile: %s" % linksfile)
		#	self.addChildTarget(CreateLinksFile(pevntsfile, linksfile, self.totalp))		

		#Annotating Events is no longer an option. Seghists are annotated instead.
		#annotationfile=os.path.join(outputdir, "evnts.ann")
		#only create the annotations file here if we aren't doing gene ranking.  Otherwise the gene rank option will create the annotation file for itself.  
		#if opts.ann and not opts.generank and not os.path.exists(annotationsfile): 
		#	logger.info("annotationfile: %s" % annotationfile)
		#	if not self.events: 
		#		self.events=pickle.load(open(pevntsfile, 'rb'))
		#	self.addChildTarget(CreateAnnotationFile(self.events, opts.tabixfile, annotationfile))
	
		# generank isn't an option - geneords is done instead using seghists. 		
		#generankfile=os.path.join(outputdir, "generanks.txt")
		# annotation file comes before generankfile (gene ranking depends on annotations.) 
		#if opts.generank and not os.path.exists(generankfile): 
		#	self.logToMaster("Creating generankfile: %s" % generankfile)
		#	logger.info("generankfile: %s" % generankfile)
		#	if not self.events: 
		#		self.events=pickle.load(open(pevntsfile, 'rb'))
		#		self.historyScores=np.loadtxt(historystatsfile, dtype=int)
		
		if opts.mcmcmix: 	
			self.logToMaster("Setting up MCMC analysis")
			mcmcdir=os.path.join(outputdir, "mcmcdata")
			mcmcdat=os.path.join(mcmcdir, "edge_counts.dat")
			mcmcdir=os.path.join(outputdir, "mcmcdata")
			mcmcdat=os.path.join(mcmcdir, "edge_counts.dat")
			if not os.path.exists(mcmcdir) or not os.path.exists(mcmcdat):
				subprocess.call("mkdir -p %s" % mcmcdir, shell=True)	
				opts.pevnts=pevntsfile 
				opts.pedges=pedgesfile
				self.addChildTarget(mcmcjobtree.SetupMCMC(opts, mcmcdir))

		if opts.simulation: 
			self.logToMaster("Setting up Simulation analysis")
			simoutput=os.path.join(outputdir, "events.stats")
			if ((not os.path.exists(simoutput)) or (os.path.getsize(simoutput) == 0)): 
				self.addChildTarget(SimAnalysisJob(pevntsfile, opts.trueID, self.historyScores, "events", outputdir, opts.binwidth))
			simoutput2=os.path.join(outputdir, "edges.stats")
			if ((not os.path.exists(simoutput2)) or (os.path.getsize(simoutput2) == 0)): 
				self.addChildTarget(SimAnalysisJob(pedgesfile, opts.trueID, self.historyScores, "edges", outputdir, opts.binwidth))
			simoutput3=os.path.join(outputdir, "mrgedges.stats")
			if ((not os.path.exists(simoutput3)) or (os.path.getsize(simoutput3) == 0)): 
				self.addChildTarget(SimAnalysisJob(mrgpedgesfile, opts.trueID, self.historyScores, "mrgedges", outputdir, opts.binwidth))


class SimAnalysisJob(Target): 
	def __init__(self, peventsfile, trueID, historyScores, outname, outputdir, binwidth): 
		Target.__init__(self)
		self.events=pickle.load(open(peventsfile, 'rb'))
		self.outname=outname
		self.trueID = trueID
		self.historyScores=historyScores
		self.outputdir=outputdir
		self.binwidth=binwidth
	
	def run(self): 
		datout=open(os.path.join(self.outputdir, "%s.dat" % self.outname), 'w')
		statout=open(os.path.join(self.outputdir, "%s.stats" % self.outname), 'w')
		#breaksfile=open(os.path.join(self.outputdir, "breakpoints.txt"), 'w')
		breaksfile=""
		histseg.Global_BINWIDTH=self.binwidth
		analyze_simulation.analyze_simulation(self.events, self.trueID, self.historyScores, datout, statout, breaksfile, self.outputdir)

class CreateLinksFile(Target): 
	def __init__(self, pevntsfile, linksfile, totalp): 
		Target.__init__(self)
		self.pevntsfile=pevntsfile
		self.linksfile=linksfile
		self.totalp=totalp
	
	def run(self): 
		self.logToMaster("CreateLinksFile\n") 
		opts=self.options
		myargs=['--inpickle', self.pevntsfile, '--links', self.linksfile, '--totalp', str(self.totalp)]
		parser = argparse.ArgumentParser()
		add_event_link_options(parser)
		args=parser.parse_args(myargs)
		score_and_link_cycles(args)

class CreatePedgesFile(Target): 
	def __init__(self, events, pedgesfile, historyScores, totalp, ignore_cn): 
		Target.__init__(self)
		self.events=events
		self.pedgesfile=pedgesfile
		self.historyScores=historyScores
		self.totalp=totalp
		self.ignore_cn=ignore_cn	
	
	def run(self):
		self.logToMaster("CreatePedgesFile, ignore_cn: %s\n" % self.ignore_cn) 
		edges=score_edges_within_pevents(self.events, self.historyScores, self.totalp, ignore_cn=self.ignore_cn)
		pickle.dump(edges, open(self.pedgesfile, 'wb'), pickle.HIGHEST_PROTOCOL) 
	
class CreateMergedEventsFile(Target): 
	def __init__(self, events, outpfn, historyScores):
		Target.__init__(self)
		self.events=events
		self.outpfn=outpfn
		self.historyScores=historyScores
	
	def run(self):
		self.logToMaster("CreateMergedEventsFile\n")
		events=self.events
		for e in events: 
			e.unpack() 
		mrged=histseg.merge_events_by_type(events, self.historyScores)
		for e in mrged: 
			e.trim()
		pickle.dump(mrged, open(self.outpfn, 'wb'), pickle.HIGHEST_PROTOCOL)
 
class CreateAnnotationFile(Target):
	def __init__(self, evnts, tabixfile, annotationfile): 
		Target.__init__(self)
		self.evnts=evnts
		self.tabixfile=tabixfile
		self.annfile=annotationfile
			
	def run(self): 
		self.logToMaster("CreateAnnotationFile\n") 
		outputfh=open(self.annfile, 'w')
		allevents=self.evnts
		mytabix=pysam.Tabixfile(self.tabixfile, 'r')
		annotate_events.main(allevents, mytabix, outputfh)	

def add_analysis_options(parser): 
	group=OptionGroup(parser, "CNAVG Post Analysis Options")
	group.add_option("--sampleid", dest="sampleid", help="the sampleid")
	group.add_option("--cnavgout", dest="cnavgout", help="The CN-AVG output directory for the sample")
	group.add_option("--outputdir", dest="outputdir", help="The output directory for the analysis")
	group.add_option('--binwidth', dest='binwidth', help='the multiplier between history ids of independent runs', default=histseg.Global_BINWIDTH, type="int")
	group.add_option('--pevnts', dest="pevnts", default=False, action="store_true", help="Rewrite .pevnts file.  This file by default is created if it doesn't exist.")
	group.add_option('--pedges', dest="pedges", default=False, action="store_true", help="Rewrite .pedges file.  This file is created by default if it doesn't exist.")
	group.add_option('--sgh', dest="sgh", default=False, action="store_true", help="Rewrite .seghists.txt file.  This file is created by default if it doesn't exist.")
	group.add_option('--bedfile', dest="bedfile", help="The bed file containing gene annotations you are interested in.  If given, the seghists will be annotated using these genes.")
	group.add_option('--geneords', dest="geneords", default=False, action="store_true", help="Rewrite geneords.txt file.  This is created by default if given a bedfile.")
	#group.add_option('--tabixfile', dest="tabixfile", help="The tabix file containing gene annotations you are interested in.  If given, the seghists will be annotated using these genes.")
	#group.add_option('--generank', dest="generank", default=False, action="store_true", help="Create .gnrnk file if it doesn't exist.")
	#group.add_option('--links', dest="links", default=False, action="store_true", help="create or rewrite .links file")
#	group.add_option('--ann', dest="ann", default=False, action="store_true", help="create .annotation file if it doesn't exist.")
	group.add_option('--mcmcmix', dest="mcmcmix", default=False, action="store_true", help="do analysis to look at the mixing across and within runs.")
	group.add_option('--simulation', dest="simulation", default=False, action="store_true", help="do simulation analysis.")
	group.add_option('--trueID', dest="trueID", default=0, help="The history id of the true simulated history.", type="int")
	parser.add_option_group(group)

def main(): 
	parser = OptionParser(usage = "cn-avg_post_analysis_jobtree.py --cnavgout CN-AVG_outputdir --sampleid Sampleid --jobTree jobtreedir")
	
	add_analysis_options(parser)
	mcmcjobtree.add_mcmc_options(parser)
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	i = Stack(Setup(options)).startJobTree(options)
	if i: 
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == "__main__": 
	from cnavg_post_analysis_jobtree import *
	main() 
