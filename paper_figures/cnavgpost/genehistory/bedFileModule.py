#!/usr/bin/env python 

import sys, os
import cnavgpost.mergehistories.event_cycles_module as histseg
import argparse
import cPickle as pickle 
import subprocess, re
import numpy as np

class BedEntry: 
	def __init__(self, bedline):
		dat=bedline.strip().split("\t")
		(chr, start, end) = dat[:3]
		self.chr=re.sub("chr", "", chr)
		self.start=int(start)
		self.end=int(end)
		self.name=""
		self.score=""
		if len(dat)>3: 
			self.name=dat[3]
		if len(dat)>4:
			self.strand=dat[4]
		if len(dat)>5:
			self.score=float(dat[5])
	
	def __str__(self): 
		#return "%s:%d-%d" % (self.chr, self.start, self.end)
		return "\t".join(map(str, (self.chr, self.start, self.end, self.name, self.score)))

	def comes_before(self, other): 
		return ((self.chr < other.chr) or (self.chr == other.chr and self.start < other.start))		
		
	def comes_after(self, other): 
		return ((self.chr > other.chr) or (self.chr==other.chr and self.start>other.start))

	def overlaps(self, region): 
		return ((self.chr == region.chr) and (self.start <= region.end) and (self.end >= region.start))

	def overlap(self, region): 
		overlap=0
		if (self.chr == region.chr): 
			overlap=min(region.end, self.end) - max(region.start, self.start)
		return overlap 
	
