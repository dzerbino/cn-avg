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

"""Edge traversals"""

import copy

###############################################
## Edge
###############################################
class Edge(object):
	def __init__(self, start, finish, value, index=-1):
		self.start = start
		self.finish = finish
		self.value = value
		self.index = index

	def nodes(self):
		return (self.start, self.finish)

	def __str__(self):
		return "%s -> %s: %f (%i)" % (self.start.ID , self.finish.ID, self.value, self.index)

	def shortString(self):
		if self.value > 0:
			sign = '+'
		else:
			sign = '-'

		if self.index > -1:
			return sign + " %s -> %s (%i)" % (self.start.shortString(), self.finish.shortString(), self.index)
		else:
			return sign + " %s -> %s (A)" % (self.start.shortString(), self.finish.shortString())

	def reverse(self):
		return Edge(self.finish, self.start, self.value, self.index)

	def getCNVs(self, hash, event, graph):
		if self.index != -1 and graph.nodeBlock[self.start] is not None:
			chain = graph.blockChain[graph.nodeBlock[self.start]]
			if chain in hash:
				hash[chain].append((self, event))
			else:
				hash[chain] = [(self, event)]
		return hash

	def adjacencyIndex(self):
		return (min(self.start, self.finish), max(self.start, self.finish), self.index)

	def dot(self):
		if self.index == -1:
			if self.start < self.finish:
				if self.value > 0:
					return str(self.finish) + '\n\t%s -> %s [color=crimson, label="+%f"]' % (self.start.ID, self.finish.ID, self.value)
				else:
					return str(self.finish) + '\n\t%s -> %s [color=blueviolet, label="%f"]' % (self.start.ID, self.finish.ID, self.value)
			else:
				if self.value > 0:
					return str(self.finish) + '\n\t%s -> %s [color=crimson, label="+%f"]' % (self.finish.ID, self.start.ID, self.value)
				else:
					return str(self.finish) + '\n\t%s -> %s [color=blueviolet, label="%f"]' % (self.finish.ID, self.start.ID, self.value)
		if self.start < self.finish:
			if self.value > 0:
				return str(self.finish) + '\n\t%s -> %s [color=blue, label="DEL:%f"]' % (self.start.ID, self.finish.ID, self.value)
			elif self.value < 0:
				return str(self.finish) + '\n\t%s -> %s [color=red, label="DUP:%f"]' % (self.start.ID, self.finish.ID, -self.value)
		else:
			if self.value > 0:
				return str(self.finish) + '\n\t%s -> %s [color=blue, label="DEL:%f"]' % (self.finish.ID, self.start.ID, self.value)
			elif self.value < 0:
				return str(self.finish) + '\n\t%s -> %s [color=red, label="DUP:%f"]' % (self.finish.ID, self.start.ID, -self.value)


	def braneyText(self, historyID, netID, cycleID, edgeIndex, order, complexity, ptr):
		if self.start.orientation:
			startOString = '+'
		else:
			startOString = '-'

		if self.finish.orientation:
			finishOString = '+'
		else:
			finishOString = '-'

		if self.index == -1:
			str1 = "\t".join(['A', self.start.chr, str(self.start.pos), startOString, self.finish.chr, str(self.finish.pos), finishOString, str(self.value), str(historyID), str(netID), str(cycleID), str(edgeIndex), str(order), str(complexity), str(ptr)])
			str2 = "\t".join(['A', self.finish.chr, str(self.finish.pos), finishOString, self.start.chr, str(self.start.pos), startOString, str(self.value), str(historyID), str(netID), str(cycleID), str(edgeIndex), str(order), str(complexity), str(ptr)])
			return "\n".join([str1, str2])
		else:
			start = min(self.start.pos, self.finish.pos)
			finish = max(self.start.pos, self.finish.pos)
			return "\t".join([self.start.chr, str(start), str(finish), str(self.value), str(historyID), str(netID), str(cycleID), str(edgeIndex), str(order), str(complexity), str(ptr)])

	def doesNotContainNodes(self, nodes):
		return self.start not in nodes and self.finish not in nodes
