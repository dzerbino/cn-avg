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

import copy

from cnavg.cactus.chain import Chain

"""Definition of Net"""

###########################################
## Net structure 
###########################################
class Net(object):
	"""A net in a cactus graph, i.e. a set of 3-edge connected blocks"""
	def __init__(self, iter):
		self.groups = frozenset(iter)

	def __str__(self):
		return str(min(self.groups))

	def dot(self, id):
		""" GraphViz output """
		return "\n".join(['\n\tsubgraph cluster%i {' % id]
		                + map(lambda X: X[1].dot(id, X[0]), enumerate(self.groups))
				+ ['\t}'])

	def nodes(self):
		""" List of nodes included in net's groups """
		return [N for G in self.groups for N in G.nodes]

	def span(self):
		""" Bed string representation of the segments spanned by the net """
		nodes = sorted(self.nodes())
		return "\t".join(map(str, [nodes[0].chr, nodes[0].pos, nodes[-1].pos]))

	"""
	def computeBiConnectedComponents(self, blocks, visited, count, dfn, low, stack, father):
		visited.add(self)
		dfn[self] = count + 1
		low[self] = count + 1
		chains = []

		for net2, block in blocks[self]:
			if net2 == self:
				visited.add(block)
				chain = Chain([block])
				chains.append(chain)
				continue
				
			if block not in visited:
				visited.add(block)
				stack.append(block)

			if net2 not in visited:
				father[net2] = self 
				chains += net2.computeBiConnectedComponents(blocks, visited, count + 1, dfn, low, stack, father)
				if low[net2] >= dfn[self]:
					chain = Chain()
					while len(stack) > 0 and (len(chain) == 0 or chain[-1] != block):
						chain.append(stack.pop(-1))
					chains.append(chain)
				low[self] = min(low[self], low[net2])

			elif net2 != father[self]:
				low[self] = min([low[self], dfn[net2]])

		return chains
	"""

	def __add__(self, other):
		return Net(self.groups | other.groups)

	def validate(self, graph):
		"""
		Validation function
		"""
		assert len(self.groups) > 0
		assert all(id(graph.groupNet[X]) == id(self) for X in self.groups)
		assert all(X.validate(graph) for X in self.groups)
		return True

	def _groupCount(self):
		return len(self.groups)

	def _chainCount(self, graph):
		return len(graph.nets2Chains[self])
	
	def _nodeCount(self):
		return sum(X.nodeCount() for X in self.groups)

	def stats(self, graph):
		"""Returns a string of basic stats on the graph"""
		return "\t".join(map(str, [self._groupCount(), self._nodeCount(), self._chainCount(graph)]))
