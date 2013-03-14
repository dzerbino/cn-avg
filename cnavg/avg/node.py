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

"""
Definition of the Node object
"""

import sys
import copy
import random

import cnavg.basics.coords as coords

try:
    import cPickle 
except ImportError:
    import pickle as cPickle

class Node(coords.OrientedPosition):
	"""
	Node in a sequence graph, i.e. the end of a segment edge, a breakend.
	"""

	def __init__(self, ID, chr, pos, orientation, name):
		super(Node, self).__init__(chr, pos, orientation)
		self.ID = int(ID) # Mainly useful for GraphViz, not much else
		if name is None:
			self.name = str(self.ID)
		else:
			self.name = str(name)

	def __cmp__(self, other):
		"""
		Compare chromosome coordinates, possibly using ID if equal.
		"""
		assert self is not None
		assert other is not None
		if cmp(self.chr, other.chr):
			return cmp(self.chr, other.chr)
		elif cmp(self.pos, other.pos):
			return cmp(self.pos, other.pos)
		else:
			return cmp(self.ID, other.ID)

	# Had to redefine hash because the comparison function is not a strict ordering function.
	def __hash__(self):
		return id(self)

	def __str__(self):
		"""
		GraphViz representation (debugging purposes).
		"""
		if self.orientation:
			return '\t%s [label="[%s] %s:%i+"]' % (str(self.ID), self.name, self.chr, self.pos)
		else:
			return '\t%s [label="[%s] %s:%i-"]' % (str(self.ID), self.name, self.chr, self.pos)

	def shortString(self):
		"""
		Returns short identification string (mainly for debugging)
		"""
		if self.orientation:
			return '[%s]\t%s\t%i\t+' % (self.ID, self.chr, self.pos)
		else:
			return '[%s]\t%s\t%i\t-' % (self.ID, self.chr, self.pos)

	def __sub__(self, other):
		""" 
		A - B = length of chromsome from A to B (can be negative).
		"""
		assert self.chr == other.chr
		return self.pos - other.pos + 1

class StubNode(Node):
	"""
	A node which represents an unidentified position in the genome
	"""
	def __init__(self, ID):
		super(StubNode, self).__init__(ID, "None", 0, True, "Stub")
