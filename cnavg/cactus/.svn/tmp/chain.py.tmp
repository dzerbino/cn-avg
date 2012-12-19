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

"""Definition of Chain"""

###########################################
## Chain structure 
###########################################
class Chain(list):
	"""A chain in a Cactus graph, i.e. a cycle of blocks which are always found in a consistent order and orientation"""
	def update(self, netA, netB, newNet):
		if self[0] == netA or self[0] == netB:
			self.pop(0)
			self.insert(0, newNet)
		elif self[-1] == netA or self[-1] == netB:
			self.pop(-1)
			self.append(newNet) 
		return self

	def __str__(self):
		return "\t".join(map(str, self))

	def nodes(self):
		return sum((list(X.nodes) for X in self), [])

	def span(self):
		nodes = sorted(self.nodes())
		chr = nodes[0].chr
		assert all(X.chr == chr for X in nodes), str(self)
		return "\t".join(map(str, [chr, nodes[0].pos, nodes[-1].pos]))

	def __hash__(self):
		return id(self)

	def lengths(self):
		return map(lambda X: X.length(), self)

	def length(self):
		return sum(self.lengths())

	def validate(self, graph):
		assert len(self) > 0
		assert len(self) > 1 or len(set(graph.groupNet[graph.nodeGroup[x]] for x in self[0].nodes)) == 1
		assert all((graph.blockChain[X] == self for X in self))
		assert all((X.validate(graph) for X in self))
		return True

	def nextNode_block(self, node):
		return filter(lambda X: id(X) != id(node) and graph.nodeNet(node) == graph.nodeNet(X) for X in self.nodes())
