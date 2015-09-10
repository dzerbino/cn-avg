#!/usr/bin/env python 

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger
from sonLib.bioio import system
import os, sys

import cnavgmbin.diagnostics.mcmc_mixing_analysis_jobtree as mcmcjobtree
import cnavg_post_analysis_jobtree as postjobtree

class Setup(Target): 
	def __init__(self, options): 
		Target.__init__(self)
		self.options=options
		self.homedir=options.outputdir
	
	def run(self): 
		opts=self.options
		self.logToMaster("setting up...")
		mylines=open(opts.samplelist, 'r').readlines()
		for line in mylines: 
			(cnavgout, sampleid) = line.strip().split('\t')
			self.logToMaster("working on %s in %s" % (sampleid, cnavgout))
			opts.sampleid=sampleid
			opts.cnavgout=cnavgout
			opts.outputdir=os.path.join(self.homedir, opts.sampleid)
			self.addChildTarget(postjobtree.Setup(opts))
	

def main(): 
	parser = OptionParser(usage = "note: This currently doesn't work.  run_cn-avg_post_analysis_jobtree.py --samplelist samplelist.txt")
	parser.add_option("--samplelist", dest="samplelist", help="The list of CNAVG outputs and sample ids. Should have the form <directory><ID>")
	postjobtree.add_analysis_options(parser)
	mcmcjobtree.add_mcmc_options(parser)
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	i = Stack(Setup(options)).startJobTree(options)
	if i:
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)


if __name__=="__main__": 
	from run_cnavg_post_analysis_jobtree import *
	main()
