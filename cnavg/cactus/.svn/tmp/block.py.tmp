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

import cnavg.avg.graph as avg
import copy
import threeWay.components

from cnavg.cactus.group import Group
 
"""Definition of Block"""

###########################################
## Block structure 
###########################################
class Block(object):
	"""Cactus graph block, i.e. a segment edge"""
	def __init__(self, nodeA, nodeB, graph):
		self.nodes = frozenset([nodeA, nodeB])

	def __str__(self):
		return "[" + ";".join(map(str, self.nodes)) + "]"

	def __cmp__(self, other):
		# Note: actually hacky
		return cmp(min(self.nodes), min(other.nodes))

	def length(self):
		nodes = list(self.nodes)
		if nodes[0] > nodes[1]:
			return nodes[0] - nodes[1]
		else:
			return nodes[1] - nodes[0]

	def copynumber(self, cactus, index):
		return cactus[list(self.nodes)[0]].segment[index]

	def ploidy(self, cactus):
		return len(cactus[list(self.nodes)[0]].segment)

	def otherNode(self, node):
		assert node in self.nodes
		nodeA, nodeB = self.nodes
		if nodeA == node:
			return nodeB
		else:
			return nodeA

	def validate(self, graph):
		nodeA, nodeB = self.nodes
		assert graph.nodeBlock[nodeA] == self
		assert graph.nodeBlock[nodeB] == self
		assert nodeA.orientation != nodeB.orientation
		return True
