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

"""Definition of Group, a connected component of bonds"""

import cnavg.avg.graph as avg
import copy

###########################################
## Group structure 
###########################################

class Group(object):
	"""A group in a cactus graph, i.e. a connected component of bonds"""
	def __init__(self, iter):
		self.nodes = frozenset(iter)

	def __iter__(self):
		return iter(self.nodes)

	def __str__(self):
		return str(min(self))

	def dot(self, netID, groupID):
		""" GraphViz representation """
		return "\n".join(['\n\t\tsubgraph cluster%i_%i {' % (netID, groupID)]
		                + map(lambda X: "\t\t" + str(X.ID), sorted(self.nodes))
				+ ['\t\t}'])

	def validate(self, graph):
		""" Validation function """
		assert all((graph.nodeGroup[X] == self for X in self.nodes))
		return True

	def nodeCount(self):
		""" Count nodes in group """
		return len(self.nodes)
