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

"""Definition of Cactus graphs"""

import cnavg.avg.graph as avg
import copy
import threeWay.components

from cnavg.cactus.group import Group
from cnavg.cactus.block import Block
from cnavg.cactus.chain import Chain
from cnavg.cactus.net import Net

		
def _updateMapping(mapping, pair):
	A, B, block = pair
	if A not in mapping:
		mapping[A] = []
	if B not in mapping:
		mapping[B] = []

	mapping[A].append((B, block))
	if B is not A:
		mapping[B].append((A, block))
	return mapping

def _finalDest(prev, curr, bridges):
	assert len(bridges[curr]) > 0
	if len(bridges[curr]) != 2:
		return curr
	else:
		next = filter(lambda X: X[0] != prev, bridges[curr])[0][0]
		return _finalDest(curr, next, bridges)

def _collapseBridge(nets, net, bridges):
	if len(bridges[net]) == 0 or len(bridges[net]) == 2:
		return nets
	else:
		mergedNets = [net] + map(lambda X: _finalDest(net, X[0], bridges), bridges[net])
		return filter(lambda X: X not in mergedNets, nets) + [sum(mergedNets, Net([]))]

class Cactus(avg.Graph):
	"""A Cactus graph"""
	##########################################
	## Basics
	##########################################
	def __init__(self, graph):
		super(Cactus, self).__init__()
		super(Cactus, self).copy(graph)

		# Compute basic components
		groups = self.computeGroups()
		self.nodeGroup = dict((N,G) for G in groups for N in G.nodes)
		blocks = self.computeBlocks(groups)
		self.nodeBlock = dict((N, B) for B in set(P[1] for X in blocks.values() for P in X) for N in B.nodes)

		# Collapse 3 way connected components
		self.nets = self.computeBasicNets(groups)
		self.groupNet = dict((G, N) for N in self.nets for G in N.groups)
		self.nets = self.computeNets(blocks)
		self.groupNet = dict((G, N) for N in self.nets for G in N.groups)
		self.chains = self.computeChains(blocks)
		self.blockChain = dict((B, C) for C in self.chains for B in C)

		# Collapse bridges 
		self.nets = self._collapseBridges(blocks)
		self.groupNet = dict((G, N) for N in self.nets for G in N.groups)
		self.chains = self.computeChains(blocks)
		self.blockChain = dict((B, C) for C in self.chains for B in C)

	def copy(self, origin):
		super(Cactus, self).copy(origin)
		self.nets = origin.nets
		self.chains = origin.chains
		self.nodeBlock = origin.nodeBlock
		self.blockChain = origin.blockChain
		self.nodeGroup = origin.nodeGroup
		self.groupNet = origin.groupNet

	##########################################
	## Output
	##########################################
	def dot(self):
		""" GraphViz output """
		return "\n".join(['digraph G {'] 
				 + map(str, sorted(self.values()))
				 + map(lambda X: X[1].dot(X[0]), enumerate(self.nets))
				 + ['}'])
	def __str__(self):
		return self.dot()

	def netStats(self):
		""" String of stats on nets """
		return "\n".join(map(lambda X: X.stats(self), self.nets))

	
	def chainSpans(self):
		""" Spans of chains """
		return "\n".join(X.span() for X in self.chains)

	##########################################
	## Shorthands
	##########################################
	def nodeNet(self, node):
		""" Returns the net which contains node """
		return self.groupNet[self.nodeGroup[node]]

	def nodeChain(self, node):
		""" Returns the chain which contains the node """
		return self.blockChain[self.nodeBlock[node]]

	##########################################
	## Computing groups
	##########################################
	def expandGroup(self, list, node, visited):
		if node in visited:
			return list
		visited[node] = True
		list.append(node)
		
		list = self.expandGroup(list, self[node].partner, visited)
		return reduce(lambda X,Y: self.expandGroup(X, Y, visited),  self[node].edges.keys(), list)

	def computeGroups(self):
		groups = []
		visited = dict()
		for node in self.nodes():
			if node not in visited:
				groups.append(Group(self.expandGroup([], node, visited)))
		return groups

	##########################################
	## Computing Blocks
	##########################################
	def groupPairs3(self, group):
		return [(group, self.nodeGroup[self[node].twin], Block(node, self[node].twin, self)) for node in group if node < self[node].twin]

	def groupPairs2(self, groups):
		return map(lambda X: self.groupPairs3(X), groups)

	def groupPairs(self, groups):
		return sum(self.groupPairs2(groups), [])

	def computeBlocks(self, groups):
		return reduce(_updateMapping, self.groupPairs(groups), dict())

	##########################################
	## Computing Initial nets
	##########################################
	def telomericGroups(self):
		return set(self.nodeGroup[X] for X in self.telomeres)

	def computeBasicNets(self, groups):
		telomericGroups = self.telomericGroups()
		return [Net([group]) for group in groups if group not in telomericGroups] + [Net(telomericGroups)]

	##########################################
	## Computing Net Blocks
	##########################################
	def updateNetGroupGroupMapping(self, netBlocks, net2, block, net):
		if net not in netBlocks:
			netBlocks[net] = []
		netBlocks[net].append((net2, block))
		return netBlocks

	def updateNetGroupMapping(self, netBlocks, group, net, blocks):
		if group in blocks:
			return reduce(lambda A,B: self.updateNetGroupGroupMapping(A, self.groupNet[B[0]], B[1], net), blocks[group], netBlocks)
		else:
			return netBlocks

	def updateNetMapping(self, netBlocks, net, blocks):
		return reduce(lambda A,B: self.updateNetGroupMapping(A, B, net, blocks), net.groups, netBlocks)

	def computeNetBlocks(self, blocks):
		return reduce(lambda A,B: self.updateNetMapping(A, B, blocks), self.nets, dict())

	##########################################
	## Computing 3 way connected components
	##########################################
	def netConnections(self, net, netIDs, netBlocks):
		return [netIDs[X[0]] for X in netBlocks[net]]

	def connections(self, netIDs, blocks):
		netBlocks = self.computeNetBlocks(blocks)
		return map(lambda X: self.netConnections(X, netIDs, netBlocks), self.nets)

	def convertThreeWayResultLine(self, threeWayResultLine, netIDs):
		return sum((self.nets[X] for X in threeWayResultLine), Net([]))

	def convertThreeWayResults(self, netIDs, threeWayResults):
		return map(lambda X: self.convertThreeWayResultLine(X, netIDs), threeWayResults)

	def computeNets(self, blocks):
		netIDs = dict((X[1], X[0]) for X in enumerate(self.nets))
		threeWayresult = threeWay.components.compute(self.connections(netIDs, blocks))
		return self.convertThreeWayResults(netIDs, threeWayresult)

	##########################################
	## Computing chains
	##########################################

	def computeChains(self, blocks):
		# Initialization
		start = self.nets[0]
		stack = []
		dfn = dict()
		father = dict()
		fatherBlock = dict()
		chains = []
		netBlocks = self.computeNetBlocks(blocks)

		# Arbitrary depth first numbering
		stack = [(start, None, None, 0)]
		while len(stack) > 0:
			net, parent, block, depth = stack.pop()
			if net not in dfn:
				dfn[net] = depth
				father[net] = parent
				fatherBlock[net] = block
				stack.extend((X[0], net, X[1], depth + 1) for X in netBlocks[net])

		# Look for loop closures
		for net in self.nets:
			for net2, block in netBlocks[net]:
				if block is not fatherBlock[net] and dfn[net2] <= dfn[net]:
					# Aha: loop closure: follow fathers to self:
					chain = Chain([block])
					net3 = net
					while net3 is not net2:
						chain.append(fatherBlock[net3])
						net3 = father[net3]
					chains.append(chain)

		return chains
			
	##########################################
	## Collapsing bridges
	##########################################

	def removeChainedBlocks(self, chainedBlocks, blocks, net, netBlocks):
		blocks[net] = filter(lambda X: X[1] not in chainedBlocks, netBlocks[net])
		return blocks

	def bridges(self, netBlocks):
		chainedBlocks = reduce(lambda X, Y: X | set(Y), self.chains, set())
		return reduce(lambda X,Y: self.removeChainedBlocks(chainedBlocks, X, Y, netBlocks), self.nets, dict())

	def _collapseBridges(self, blocks):
		bridges = self.bridges(self.computeNetBlocks(blocks))
		return reduce(lambda X,Y: _collapseBridge(X, Y, bridges), bridges.keys(), self.nets)

	##########################################
	## Find root net
	##########################################
	def findRootNet(self):
		telomericNets = set(self.groupNet[self.nodeGroup[X]] for X in self.telomeres)
		return list(telomericNets)[0]

	##########################################
	## Validation
	##########################################
	def _neighbours(self, net, parent, netBlocks):
		return (X[0] for X in filter(lambda X: X[0] != net and X[0] != parent, netBlocks[net]))

	def _propagateBetweenNets(self, net, parent, netB, netBlocks, visited):
		if net not in visited:
			visited[net] = 1
			if net != netB:
				return reduce(lambda V,N: self._propagateBetweenNets(N, net, netB, netBlocks, V), self._neighbours(net, parent, netBlocks), visited)
			else:
				return visited
		else:
			visited[net] += 1
			return visited

	def _nothreeWayComponentBetweenNets(self, netA, netB, netBlocks):
		if netA != netB:
			visited = self._propagateBetweenNets(netA, netA, netB, netBlocks, dict())
			assert netB not in visited or visited[netB] == 2
		return True

	def _nothreeWayComponent(self):
		groups = [G for N in self.nets for G in N.groups]
		blocks = self.computeBlocks(groups) 
		netBlocks = self.computeNetBlocks(blocks)
		return all(map(lambda X: self._nothreeWayComponentBetweenNets(self.nets[0], X, netBlocks), self.nets))

	def _countMapping(self, mapping, unit):
		if unit not in mapping:
			mapping[unit] = 1
		else:
			mapping[unit] += 1
		return mapping

	def _group2NodeMappings2(self, mapping, group):
		return reduce(lambda X, Y: self._countMapping(X,Y), group, mapping)

	def _group2NodeMappings(self):
		return reduce(lambda X, Y: self._group2NodeMappings2(X,Y), [X for N in self.nets for X in N.groups], dict())

	def _allNodesInOneGroup(self):
		assert all(map(lambda X: X == 1, self._group2NodeMappings().values()))
		return True

	def _net2GroupMappings2(self, mapping, net):
		return reduce(lambda X, Y: self._countMapping(X,Y), net.groups, mapping)

	def _net2GroupMappings(self):
		return reduce(lambda X, Y: self._net2GroupMappings2(X,Y), self.nets, dict())

	def _allGroupsInOneNet(self):
		assert all(map(lambda X: X == 1, self._net2GroupMappings().values()))
		return True

	def _validateNets(self):
		assert all(map(lambda X: X.validate(self), self.nets))
		return True

	def _validateChains(self):
		assert all(map(lambda X: X.validate(self), self.chains))
		return True

	def _allNodeBlocksInSelf(self):
		assert all(X in self.blockChain for X in self.nodeBlock.values())
		return True

	def validate(self):
		""" Validation function """
		assert self._validateNets() 
		assert self._validateChains() 
		assert super(Cactus, self).validate() 
		assert self._nothreeWayComponent() 
		assert self._allNodesInOneGroup() 
		assert self._allGroupsInOneNet() 
		assert self._allNodeBlocksInSelf()
		return True

###########################################
## Unit test
###########################################

def main():
	G = avg.randomEulerianGraph(10)
	print G
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	cactus = Cactus(G)
	assert all(X in cactus.blockChain for X in cactus.nodeBlock.values())
	print cactus
	assert all(X in cactus.blockChain for X in cactus.nodeBlock.values())
	cactus.validate()
	for chain in cactus.chains:
		print chain

if __name__ == "__main__":
	main()
