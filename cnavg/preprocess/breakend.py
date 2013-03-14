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

"""VCF style representation of breakend data"""

import sys
import cnavg.basics.coords as coords

###############################################
## Rearrangement lingo
###############################################

class Breakend(coords.OrientedRegion):
	def __init__(self, chr, start, orientation, ID, finish=None):
		if finish is None:
			super(Breakend, self).__init__(chr, start, start + 1, orientation)
		else:
			super(Breakend, self).__init__(chr, start, finish, orientation)
		self.ID = ID

		self.partner = None
		self.mates = []
		self.remoteChr = []
		self.remoteStart = []
		self.remoteFinish = []
		self.remoteOrientation = []
		self.adjacency_cov = []
		self.segment = [0]

	def __str__(self):
		items = ["[%s] %s" % (self.ID, super(Breakend, self).__str__())] 

		if self.partner is not None and self.partner > self:
			items += ["[%s] %s <-> [%s] %s | %f" % (self.ID, super(Breakend, self).__str__(), self.partner.ID, super(Breakend, self.partner).__str__(), self.adjacency_cov[0])]
		for mateIndex in range(len(self.mates)):
			mate = self.mates[mateIndex]
			if mate is not None and mate > self:
				items += ["[%s] %s <-> [%s] %s | %f" % (self.ID, super(Breakend, self).__str__(), mate.ID, super(Breakend, mate).__str__(), self.adjacency_cov[mateIndex + 1])]

		return "\n".join(items)

	#########################################################
	## Consolidating missing labels
	#########################################################

	def createPartner(self, graph):
		partnerID = str(self.ID) + "_partner"
		if self.orientation:
			partner = Breakend(self.chr, self.start - 1, False, partnerID, finish=self.finish - 1)
		else:
			partner = Breakend(self.chr, self.start + 1, True, partnerID, finish=self.finish + 1)
		hits = filter(lambda X: X.partner is None, graph.searchBreakend(partner))
		if len(hits) == 1:
			if len(hits) >= 2:
				print 'COLLISION'
				print self
				print '>>>>>>>>>>>>>>>>>'
				print '\n'.join(map(str, hits))
			assert len(hits) < 2
			self.partner = hits[0] 
			self.partner.partner = self
			if len(self.adjacency_cov) == 0:
				self.adjacency_cov.append(-1)
			partner.adjacency_cov = [self.adjacency_cov[0]]
			if not max(self.start, self.partner.start) <= min(self.finish, self.partner.finish):
				print self
				print self.partner
				assert False
			cutpoint = (max(self.start, self.partner.start) + min(self.finish, self.partner.finish)) / 2
			if self.orientation:
				self.start = cutpoint + 1
				self.partner.finish = cutpoint
				assert cutpoint + 1 < self.finish
				assert cutpoint > self.partner.start
			else:
				self.partner.start = cutpoint + 1
				self.finish = cutpoint
				assert cutpoint + 1 < self.partner.finish
				assert cutpoint > self.start
			self.partner.start = max(self.partner.start, partner.start)

			if (not self.orientation and self.start >= self.partner.start) or (self.orientation and self.start < self.partner.start):
				print self
				print self.partner
				sys.exit(1)
			return

		partner.partner = self
		if len(self.adjacency_cov) == 0:
			self.adjacency_cov.append(-1)
		partner.adjacency_cov = [self.adjacency_cov[0]]

		self.partner = partner
		if self.orientation:
			assert self > self.partner
		else:
			assert self < self.partner
		graph.addBreakend(partner)

	def createMate(self, index, graph):
		mate = Breakend(self.remoteChr[index], self.remoteStart[index], self.remoteOrientation[index], str(self.ID) + "_mate" + str(index), finish=self.remoteFinish[index])
		
		hits = graph.searchBreakend(mate)
		if len(hits) > 0:
			if len(hits) >= 2:
				hits = filter(lambda X: id(X) != id(self), hits)
				assert len(hits) > 0
			if len(hits) >= 2:
				hitDistances = sorted((abs(X.position() - self.position()), X) for X in hits)
				hits = [hitDistances[0][1]]
			assert len(hits) < 2
			mate = hits[0]
			self.mate = mate 
			self.mates[index] = mate

			remoteIndex = filter(lambda X: mate.remoteChr[X] == self.chr and mate.remoteOrientation[X] == self.orientation and mate.remoteStart[X] < self.finish and mate.remoteFinish[X] > self.start, range(len(mate.mates)))
			assert len(remoteIndex) < 2
			if len(remoteIndex) == 1:
				mate.mates[remoteIndex[0]] = self
				mate.adjacency_cov[remoteIndex[0]] = self.adjacency_cov[index]
			else:
				mate.mates.append(self)
				mate.remoteChr.append(self.chr)
				mate.remoteStart.append(self.start)
				mate.remoteFinish.append(self.finish)
				mate.remoteOrientation.append(self.orientation)
				mate.adjacency_cov.append(self.adjacency_cov[index])
		else:
			self.mates[index] = mate
			mate.mates = [self]
			mate.remoteChr = [self.chr]
			mate.remoteStart = [self.start]
			mate.remoteFinish = [self.finish]
			mate.remoteOrientation = [self.orientation]
			mate.adjacency_cov = [-1, self.adjacency_cov[index]]
			mate.segment = [0]

			mate.createPartner(graph)
			graph.addBreakend(mate)

	def consolidateEmptyBreakend(self, breakendGraph):
		if self.partner is None:
			self.createPartner(breakendGraph)

		for counter in range(len(self.mates)):
			if self.mates[counter] is None:
				self.createMate(counter, breakendGraph)

	########################################################
	## Conversion to AVG
	########################################################

	def attachNode(self, graph):
		self.node = graph.createNode(self.chr, self.start, self.orientation, name=self.ID)

	def connectNode(self, graph):
		graph.createAdjacency(self.node, self.partner.node)
		if self.node < self.partner.node:
			graph.addLiftedEdge(self.node, self.partner.node, self.adjacency_cov[0])

		for index in range(len(self.mates)):
			if self.node < self.mates[index].node:
				graph.addLiftedEdge(self.node, self.mates[index].node, self.adjacency_cov[index + 1])

	def validate(self):
		if self.partner is not None:
			assert self.partner.partner == self
			if self < self.partner:
				assert self.orientation == False and self.partner.orientation == True
			else:
				assert self.orientation == True and self.partner.orientation == False
		if self.mates is not None:
			assert all(self in mate.mates for mate in self.mates if mate is not None)
		assert self.finish >= self.start
		assert sum(self.segment) >= 0
		return True
