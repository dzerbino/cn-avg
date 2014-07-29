# Copyright (c) 2012, Daniel Zerbino
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# (1) Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer. 
# 
# (2) Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in
# the documentation and/or other materials provided with the
# distribution.  
# 
# (3)The name of the author may not be used to
# endorse or promote products derived from this software without
# specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#!/usr/bin/env python

"""Net flow change between two genomes"""

import copy
import math

from cnavg.flows.cycle import Cycle

#############################################
## Event 
#############################################
class Event(object):
	#############################################
	## Basic
	#############################################
	def __init__(self, cycle, ratio=None):
		self.cycle = cycle 
		if ratio is None:
			self.ratio = abs(self.cycle[0].value) / math.ceil(abs(self.cycle[0].value))
		else:
			self.setRatio(ratio)

	def __copy__(self):
		return Event(copy.copy(self.cycle))

	#############################################
	## Stats
	#############################################
	def length(self):
		return len(self.cycle)

	def validate(self):
		if len(self.cycle) == 0:
			print self
		assert self.cycle.validate()
		return True

	def getCNVs(self, hash, graph):
		return self.cycle.getCNVs(hash, self, graph)
	
	def setRatio(self, value):
		self.cycle.setRatio(value)
		self.ratio = abs(value) / math.ceil(abs(value))

	#############################################
	## Display
	#############################################
	def __str__(self):
		return "\n".join(["EVENT %f %i" % (self.ratio, id(self)), str(self.cycle)])

	def dot(self):
		return "\n".join(["digraph G {","\trankdir=LR",self.cycle.dot(),"}"])

	def braneyText(self, historyID, netID, cycleID, ordering, complexity, eventCosts=None):
		if eventCosts is not None and self in eventCosts:
			return self.cycle.braneyText(historyID, netID, cycleID, ordering.depth[self], complexity, eventCosts[self][0], eventCosts[self][1], id(self), self.ratio)
		else:
			return self.cycle.braneyText(historyID, netID, cycleID, ordering.depth[self], complexity, '.', '.', id(self), self.ratio)

	def simplifyStubsAndTrivials(self, cactus):
		cycle = self.cycle.simplifyStubsAndTrivials(cactus)
		if cycle is not None:
			return Event(cycle)
		else:
			return None

	def doesNotContainNodes(self, nodes):
		return self.cycle.doesNotContainNodes(nodes)
