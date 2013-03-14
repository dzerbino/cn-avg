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

"""Definition of OrientedCactus, where nets and chains are organized into a bi-partite tree"""

import cnavg.avg.graph as avg
import graph as cactus

###########################################
## Oriented Cactus 
###########################################

def reverseHash(hash):
	return dict((v,k) for k in hash for v in hash[k]) 

class OrientedCactus(cactus.Cactus):
	"""An oriented cactus, i.e. a cactus with an overlayed rooted hierarchy of blocks and nets"""
	###########################################
	## Basics
	###########################################
	def __init__(self, cactus):
		super(OrientedCactus, self).copy(cactus)
		self.validate()
		self.rootNet = self.findRootNet()
		self.nets2Chains, self.chains2Nets = self.orient()
		self.headNet = reverseHash(self.nets2Chains)
		self.headChain = reverseHash(self.chains2Nets)

	def copy(self, other):
		super(OrientedCactus, self).copy(other)
		self.rootNet = other.rootNet
		self.nets2Chains = other.nets2Chains
		self.chains2Nets = other.chains2Nets
		self.headNet = other.headNet
		self.headChain = other.headChain

	###########################################
	## Output
	###########################################
	def chainTreeString(self, chain, offset):
		margin = "     ".join(["" for X in range(offset)])
		childSequences ="\n".join(map(lambda X: self.netTreeString(X, offset+1), self.chains2Nets[chain]))
		return margin + "[CHAIN " + str(chain) + ":\n" + childSequences + "\n" + margin + "]"

	def netTreeString(self, net, offset):
		margin = "     ".join(["" for X in range(offset)])
		if net in self.nets2Chains:
			childSequences ="\n".join(map(lambda X: self.chainTreeString(X, offset+1), self.nets2Chains[net]))
		else:
			childSequences = ""
		return margin + "[NET " + str(net) + ":\n" + childSequences + "\n" + margin + "]"

	def treeString(self):
		return self.netTreeString(self.rootNet, 1)

	def netSpans(self):
		return "\n".join(X.span() for X in self.nets if X != self.headNet)

	def __str__(self):
		return "\n".join([super(OrientedCactus, self).__str__(), self.treeString()])

	###########################################
	## Stats
	###########################################
	def maxDepth_Chain(self, chain):
		list = self.chains2Nets[chain]
		if list is None or len(list) == 0:
			return 1
		else:
			return max(map(self.maxDepth_Net, list)) + 1

	def maxDepth_Net(self, net):
		list = self.nets2Chains[net]
		if list is None  or len(list) == 0:
			return 1
		else:
			return max(map(self.maxDepth_Chain, list)) + 1

	def maxDepth(self):
		return self.maxDepth_Net(self.rootNet)

	def maxBranches(self):
		return max(map(len, self.nets2Chains.values()))

	def depth_chain(self, chain):
		return self.depth_net(self.headNet[chain]) + 1

	def depth_net(self, net):
		if net == self.rootNet:
			return 0
		else:
			return self.depth_chain(self.headChain[net]) + 1

	def stats(self):
		return "\t".join(map(str, [self.maxDepth(), self.maxBranches()]))

	def netStats_net(self, net):
		return "\t".join(map(str, [net.stats(self), self.depth_net(net)]))

	def netStats(self):
		return "\n".join(map(self.netStats_net, self.nets))

	###########################################
	## Net 2 Chain mapping
	###########################################
	def updateNetChains_Node(self, netChains, node, chain):
		net = self.groupNet[self.nodeGroup[node]]
		if net not in netChains:
			netChains[net] = set()
		if chain not in netChains[net]:
			netChains[net].add(chain)
		return netChains

	def updateNetChains_Block(self, netChains, block, chain):
		return reduce(lambda X, Y: self.updateNetChains_Node(X, Y, chain), block.nodes, netChains)

	def updateNetChains_Chain(self, netChains, chain):
		return reduce(lambda X, Y: self.updateNetChains_Block(X, Y, chain), chain, netChains)

	def computeNetChains(self):
		return reduce(lambda X, Y: self.updateNetChains_Chain(X, Y), self.chains, dict())

	###########################################
	## Tree Building 
	###########################################
	def propagateNetOrientation(self, net, parent, netChains, trees):
		nets2Chains, chains2Nets = trees
		if net not in netChains:
			return trees
		else:
			nets2Chains[net] = set(filter(lambda X: X != parent, netChains[net]))
			return reduce(lambda X, Y: self.propagateChainOrientation(Y, net, netChains, X), nets2Chains[net], (nets2Chains, chains2Nets))

	def propagateChainOrientation(self, chain, parent, netChains, trees):
		nets2Chains, chains2Nets = trees
		chains2Nets[chain] = set(filter(lambda X: X != parent, [self.groupNet[self.nodeGroup[N]] for B in chain for N in B.nodes]))
		return reduce(lambda X, Y: self.propagateNetOrientation(Y, chain, netChains, X), chains2Nets[chain], (nets2Chains, chains2Nets))

	def orient(self):
		return self.propagateNetOrientation(self.rootNet, self.rootNet, self.computeNetChains(), (dict(), dict()))

###########################################
## Unit test
###########################################

def main():
	G = avg.randomEulerianGraph(10)
	C = cactus.Cactus(G)
	print C
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	for chain in C.chains:
		print chain
	OC = OrientedCactus(C)
	print OC
	print OC.netStats()
	print OC.stats()

if __name__ == "__main__":
	main()
