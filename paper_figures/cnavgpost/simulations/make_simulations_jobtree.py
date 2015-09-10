#!/usr/bin/env python

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger

import subprocess, os, sys 

import cnavg.simulator.simulator
import cnavg.cactus.oriented
import cnavg.cactus.graph
import cPickle as pickle

CNAVG_HOMEEXE="/home/tballing/bin/"

class SetupSim(Target):
	def __init__(self, options):
		Target.__init__(self)
		self.options=options
	
	def run(self):
		opts=self.options
		specfile=opts.specfile
		outputdir=opts.outputdir
		id=opts.id 
		for specs in open(specfile, 'r'): 
			(blocks, events) = map(int, specs.strip().split())
			for i in xrange(opts.reps): 
				self.addChildTarget(RunCnavgSimulation(blocks, events, i, outputdir, opts))			

class Chdir: 
	def __init__(self, newPath):
		self.savedPath=os.getcwd()
		os.chdir(newPath)
	
	def __del__(self):
		os.chdir(self.savedPath)

class RunCnavgSimulation(Target): 
	def __init__(self, blocks, events, i, outputdir, options): 
		Target.__init__(self)
		self.options=options
		self.blocks=blocks
		self.events=events
		self.outputdir=outputdir
		self.i = i
	
	def run(self):
		opts=self.options 
		simoutdir=os.path.join(os.path.abspath(self.outputdir), "sim%d.%d.%d.%d" % (self.blocks, self.events, opts.id, self.i))
		subprocess.call("mkdir -p %s" % simoutdir, shell=True)
		cwd=os.getcwd()
		os.chdir(simoutdir)
		if opts.integer: 
			H = cnavg.simulator.simulator.RandomIntegerHistory(self.blocks, self.events, indel_length=10)
		else: 
			H = cnavg.simulator.simulator.RandomCancerHistory(self.blocks, self.events, indel_length=10)
		A = H.avg()
		C = cnavg.cactus.graph.Cactus(A)
		G = cnavg.cactus.oriented.OrientedCactus(C)
		
		file = open('true.braney', 'w')
		file.write(H.braneyText(G))
		file.close()

		file = open('CACTUS', 'w')
		pickle.dump(G, file)
		file.close()
		os.chdir(cwd)
		
		if opts.timeshuffle: 
			for j in xrange(1,(opts.runs+1)): 
				self.addChildTarget(RunShuffleForSim(simoutdir, j, opts))
		else: 			
			for j in xrange(1,(opts.runs+1)): 
				self.addChildTarget(RunCnavgForSim(simoutdir, j, opts))

class RunCnavgForSim(Target):
	def __init__(self, simoutdir, j, options): 
		Target.__init__(self)
		self.simoutdir=simoutdir
		self.j=j
		self.options=options

	def run(self):
		cwd=os.getcwd()
		os.chdir(self.simoutdir)
		opts=self.options
		if self.options.integer: 
			self.logToMaster("running:\n%s/cn-avg.py -n -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps))
			if subprocess.call("%s/cn-avg.py -n -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps), shell=True) !=0: 
				sys.exit("cn-avg.py did not complete.\n")
		else: 
			self.logToMaster("running:\n%s/cn-avg.py -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps))
			if subprocess.call("%s/cn-avg.py -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps), shell=True) != 0: 
				sys.exit("cn-avg.py did not complete.\n")
		os.chdir(cwd)

class RunShuffleForSim(Target):
	def __init__(self, simoutdir, j, options): 
		Target.__init__(self)
		self.simoutdir=simoutdir
		self.j=j
		self.options=options

	def run(self):
		cwd=os.getcwd()
		os.chdir(self.simoutdir)
		opts=self.options
		if self.options.integer: 
			self.logToMaster("running:\n%s/cn-avg.py -n -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps))
			if subprocess.call("%s/cn-avg.py -n -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps), shell=True) !=0: 
				sys.exit("cn-avg.py did not complete.\n")
			#subprocess.call("rm HISTORIES_%d.braney" % (self.j), shell=True)
			#subprocess.call("rm HISTORY_STATS_%d" % (self.j), shell=True)
			self.logToMaster("running:\n%s/cn-avg-timeshuffle.py -c -n -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps))
			if subprocess.call("%s/cn-avg-timeshuffle.py -c -n -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps), shell=True) !=0: 
				sys.exit("cn-avg-timeshuffle.py did not complete.\n")
		else: 
			self.logToMaster("running:\n%s/cn-avg.py -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps))
			if subprocess.call("%s/cn-avg.py -n -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps), shell=True) !=0: 
				sys.exit("cn-avg.py did not complete.\n")
			#subprocess.call("rm HISTORIES_%d.braney" % (self.j), shell=True)
			#subprocess.call("rm HISTORY_STATS_%d" % (self.j), shell=True)
			self.logToMaster("running:\n%s/cn-avg-timeshuffle.py -c -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps))
			if subprocess.call("%s/cn-avg-timeshuffle.py -c -d . -i %d -s %d" % (CNAVG_HOMEEXE, self.j, opts.steps), shell=True) != 0: 
				sys.exit("cn-avg-timeshuffle.py did not complete.\n")
		os.chdir(cwd)


def add_simulation_options(parser): 
	group = OptionGroup(parser, "Make Simulation Options")
	group.add_option("--specfile", dest='specfile', help='A file specifying the number of blocks and events to put in each simulation', type='string')
	group.add_option("--reps", dest='reps', type="int", default=5, help="The number of reps to do for each simulation")
	group.add_option("--outputdir", dest='outputdir', help="The directory to write to")
	group.add_option("--integer", dest='integer', default=False, action="store_true", help="Use integer simulations rather than metagenomic.")
	group.add_option("--id", dest='id', type="int", default=1, help="an id for the simulation runs.")
	group.add_option("--steps", dest='steps', type="int", default=1000, help="The number of steps per iteration to run.")
	group.add_option("--runs", dest='runs', type="int", default=10, help="The number of runs per sample to do.")
	group.add_option("--timeshuffle", dest='timeshuffle', default=False, action="store_true", help="run cn-avg-timeshuffle.py instead of cn-avg.py")
	parser.add_option_group(group)	

def main(): 
	parser= OptionParser(usage = "make_simulations_jobtree.py ....")
	add_simulation_options(parser)
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	i = Stack(SetupSim(options)).startJobTree(options)
	if i: 
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)


if __name__== "__main__": 
	from make_simulations_jobtree import * 
	main()

 
