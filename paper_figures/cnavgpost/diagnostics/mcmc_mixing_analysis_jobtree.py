#!/usr/bin/env python 

''' 
Analysis on CN-AVG output to test how well histories are being sampled across runs 
''' 
import os, sys
import re
import subprocess, glob
import cPickle as pickle
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack 
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger
import numpy as np
import cnavgpost.diagnostics.get_history_distances_between_mcmc_steps as mcmcdist 
import cnavgpost.mergehistories.event_cycles_module as histseg 


#======== MAIN COMMAND ============
class SetupMCMC(Target): 
	def __init__(self, options, outputdir): 
		Target.__init__(self)
		self.options=options
		self.outputdir=outputdir
		
	def run(self): 
		self.logToMaster("Setting up MCMC...")
		args=self.options
		#global_dir=self.getGlobalTempDir()
		global_dir=self.outputdir
		cmd = "mkdir -p %s" % global_dir 
		logger.info(cmd)
		subprocess.call(cmd, shell=True)
		edgedatadir=os.path.join(global_dir, "edgedata")
		evntdatadir=os.path.join(global_dir, "eventdata")
		if args.pevnts: 
			subprocess.call("mkdir -p %s" % evntdatadir, shell=True)
			logger.info("loading events: %s\n" % args.pevnts)
			events=pickle.load(open(args.pevnts, 'rb'))
			for e in events: 
				e.unpack()
		if args.pedges:
			subprocess.call("mkdir -p %s" % edgedatadir, shell=True)
			logger.info("loading edges: %s\n" % args.pedges)
			edges=pickle.load(open(args.pedges, 'rb'))
			for e in edges: 
				e.unpack()
		logger.info("numruns %d" % (args.numruns))
		if args.simulation: 
			myis = range(1, args.numruns)
		else: 
			myis=range(args.numruns)
		for i in myis:
			logger.info("Adding Child %d" % i)
			i2 = i+1
			if i== myis[len(myis)-1]:
				i2=myis[0]
			if args.pevnts:  
				self.addChildTarget(DoMCMCforSingleRun(i, args, events, evntdatadir))
				self.addChildTarget(DoMCMCforIndepRuns(i, i2, args, events, evntdatadir))
			if args.pedges: 
				self.addChildTarget(DoMCMCforSingleRun(i, args, edges, edgedatadir))
				self.addChildTarget(DoMCMCforIndepRuns(i, i2, args, edges, edgedatadir))
			if args.simulation: 
				if args.pedges: 
					self.addChildTarget(DoMCMCforSimulatedRuns(i, args, edges, edgedatadir))
				if args.pevnts: 
					self.addChildTarget(DoMCMCforSimulatedRuns(i, args, events, evntdatadir))
		self.setFollowOnTarget(CombineMCMCdataRuns(global_dir, evntdatadir, edgedatadir))

class DoMCMCforSimulatedRuns(Target): 
	def __init__(self, i, options, events, datadir): 
		Target.__init__(self)
		self.options=options
		self.i=i
		self.outdir=datadir
		self.events=events

	def run(self): 
		args=self.options
		refhistoryid=args.trueID
		refhistoryid2=histseg.Global_BINWIDTH* self.i + args.refhistoryid 
		outfile=os.path.join(self.outdir, "runsim_%d.dat" % (self.i))
		if not os.path.exists(outfile):
			outfh=open(outfile, 'w')
			logger.info("running %d: get_history_distances_between_mcmc_steps(events, %s, %s, %s, %s, %s) > %s" % (self.i, refhistoryid, refhistoryid2, args.numsteps, 0, args.stepsize, outfile))
			mcmcdist.get_history_distances_between_mcmc_steps(self.events, refhistoryid, refhistoryid2, args.numsteps, 0, args.stepsize, outfh) 

class DoMCMCforIndepRuns(Target): 
	def __init__(self, i, i2, options, events, datadir): 
		Target.__init__(self)
		self.options=options
		self.i=i
		self.i2=i2
		self.outdir=datadir
		self.events=events

	def run(self): 
		args=self.options
		refhistoryid=args.refhistoryid+ histseg.Global_BINWIDTH*self.i 
		i2 = self.i +1
		if self.i== args.numruns-1:
			i2=0 
		refhistoryid2=args.refhistoryid + histseg.Global_BINWIDTH* self.i2 
		outfile=os.path.join(self.outdir, "runs_%d.%d.dat" % (self.i, self.i2))
		if not os.path.exists(outfile):
			outfh=open(outfile, 'w')
			logger.info("running %d: get_history_distances_between_mcmc_steps(events, %s, %s, %s, %s, %s) > %s" % (self.i, args.refhistoryid, refhistoryid2, args.numsteps, args.stepsize, args.stepsize, outfile))
			mcmcdist.get_history_distances_between_mcmc_steps(self.events, refhistoryid, refhistoryid2, args.numsteps, args.stepsize, args.stepsize, outfh) 

class DoMCMCforSingleRun(Target): 
	def __init__(self, i, options, events, datadir): 
		Target.__init__(self)
		self.options=options
		self.i=i
		self.outdir=datadir
		self.events=events

	def run(self): 
		args=self.options
		refhistoryid=args.refhistoryid + histseg.Global_BINWIDTH*self.i 
		outfile=os.path.join(self.outdir, "run_%d.dat" % self.i)
		if not os.path.exists(outfile):
			outfh=open(outfile, 'w')
			logger.info("running %d: get_history_distances_between_mcmc_steps(events, %s, %s, %s, %s, %s) > %s" % (self.i, refhistoryid, "" , args.numsteps, args.stepsize, args.stepsize, outfile))
			mcmcdist.get_history_distances_between_mcmc_steps(self.events, refhistoryid, "", args.numsteps, args.stepsize, args.stepsize, outfh) 


class CombineMCMCdataRuns(Target): 
	def __init__(self, global_dir, evntdatadir, edgedatadir):
		Target.__init__(self)
		self.global_dir=global_dir
		self.evntdatadir=evntdatadir
		self.edgedatadir=edgedatadir

	def run(self):
		global_dir=self.global_dir 
		datafiles=glob.glob(self.evntdatadir+"/"+"run_*.dat")
		outfile=os.path.join(global_dir, "evntmixing_dep.dat")
		combine_datafiles(datafiles, outfile)
		datafiles=glob.glob(self.evntdatadir+"/"+"runs_*.*.dat")
		outfile=os.path.join(global_dir, "evntmixing_ind.dat")
		combine_datafiles(datafiles, outfile)
		outfile=os.path.join(global_dir, "evnt_counts.dat")
		combine_data_counts(datafiles, outfile)
		datafiles=glob.glob(self.evntdatadir+"/"+"runsim_*.dat")
		outfile=os.path.join(global_dir, "evntmixing_sim.dat")
		combine_datafiles(datafiles, outfile)
		datafiles=glob.glob(self.edgedatadir+"/"+"run_*.dat")
		outfile=os.path.join(global_dir, "edgemixing_dep.dat")
		combine_datafiles(datafiles, outfile)
		datafiles=glob.glob(self.edgedatadir+"/"+"runs_*.*.dat")
		outfile=os.path.join(global_dir, "edgemixing_ind.dat")
		combine_datafiles(datafiles, outfile)
		outfile=os.path.join(global_dir, "edge_counts.dat")
		combine_data_counts(datafiles, outfile)
		datafiles=glob.glob(self.edgedatadir+"/"+"runsim_*.dat")
		outfile=os.path.join(global_dir, "edgemixing_sim.dat")
		combine_datafiles(datafiles, outfile)
	
	
def combine_datafiles(datafiles, outfile): 
	mydata=np.array([], dtype=float)
	myfmts=[]
	myheader=[]
	for i in xrange(len(datafiles)): 
		datafile=datafiles[i] 
		label=re.search('(run.*)\.dat', datafile).group(1)
		logger.info("datafile is " + datafile)
		data=np.loadtxt(datafile, skiprows=1)
	#	logger.info("data is %s" % str(data.shape))
		if i ==0:
			mydata=np.mod(np.int_(data[:,0:2]), histseg.Global_BINWIDTH)
			myfmts.append("%d")	
			myfmts.append("%d")
			myheader.append("historyID_1")	
			myheader.append("historyID_2")	
		x=data[:,4]/(data[:,2]+data[:,3]-data[:,4])
		x=x.reshape(data.shape[0], 1)
		myfmts.append("%f")
		myheader.append(label)	
		logger.info("header is %s, format is %d" % ("\t".join(myheader), len(myfmts)))
		mydata=np.append(mydata, x, axis=1)
		logger.info("x is %s, mydata is %s" % (str(x.shape), str(mydata.shape)))
	if len(datafiles) >0: 
		np.savetxt(outfile, mydata, fmt="\t".join(myfmts), header="\t".join(myheader), delimiter="\t")	

def combine_data_counts(datafiles, outfile):
	mydata=np.array([], dtype=int)
	for i in xrange(len(datafiles)):
		datafile=datafiles[i]
		data=np.loadtxt(datafile, skiprows=1)
		if i==0:
			mydata=np.int_(data[:,2])
		else: 
			mydata=np.vstack((mydata, data[:,2]))
	np.savetxt(outfile, np.transpose(mydata), fmt="%d", delimiter="\t") 
			

def add_mcmc_options(parser): 
	group = OptionGroup(parser, "MCMC analysis options")
	group.add_option('--refhistoryid', dest='refhistoryid', help='the maximum number of steps in the sampling chain', default=2500, type="int")
	group.add_option('--numsteps', dest="numsteps", help='the number of steps to take.', default=0, type="int")
	group.add_option('--stepsize', dest="stepsize", help='the step size from the reference id(s) in the mcmc that you want to go.', default=-1, type="int")
	group.add_option('--numruns', dest="numruns", help='the number of independent histories.' , default=10, type="int")
	parser.add_option_group(group)

def main(): 
	parser = OptionParser(usage = "mcmc_mixing_analysis_jobtree.py --pevnts sample.pevnts --refhistoryid=2500 --numsteps=1000 --stepsize=1 --logInfo --jobTree=/outputdir --batchSystem=singleMachine")
	parser.add_option("--pevnts", dest="pevnts", help="a .pevnts file", type="string")
	parser.add_option("--pedges", dest="pedges", help="a .pedges file", type="string")
	parser.add_option("--outputdir", dest="outputdir", help="where you want the final output to go", type="string")
	parser.add_option("--simulation", dest="simulation", default=False, action="store_true", help="flag to indicate that it's a simulated history")
	parser.add_option('--trueID', dest="trueID", help='the id of the true history for simulations.', default=0, type="int")
	parser.add_option('--binwidth', dest='binwidth', help='the multiplier between history ids of independent runs', default=histseg.Global_BINWIDTH, type="int")
	add_mcmc_options(parser)
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	histseg.Global_BINWIDTH=options.binwidth
	i = Stack(SetupMCMC(options, options.outputdir)).startJobTree(options)
	if i: 
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == "__main__": 
	from mcmc_mixing_analysis_jobtree import *
	main()
	
