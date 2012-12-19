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

import overlap

import math
import random
import numpy as np
import copy
import debug

import scipy.sparse

""" Linear algebra representation of a history """

ROUNDING_ERROR=1e-10
""" Amount by which a coefficient is considered to be equal to 0 after linear simplification """ 

##############################################
## Mapping of edges onto integer IDs 
##############################################

class Mapping(dict):
	""" Maps edges of a module to integer indices """
	##############################################
	## Basic functionalities
	##############################################
	def _addMatrixEdge(self, IDA, IDB, value):
		""" The use of the matrix allows for accelerated bond lookups """
		M = max(IDA, IDB) + 1
		if M > self.maxNodeID:
			if self.maxNodeID == 0:
				self.matrix = np.zeros((M,M))
			else:
				self.matrix = np.concatenate((self.matrix, np.zeros((M - self.maxNodeID, self.maxNodeID))), axis=0)
				self.matrix = np.concatenate((self.matrix, np.zeros((M, M - self.maxNodeID))), axis=1)
			self.maxNodeID = M
		self.matrix[IDA, IDB] = value
		self.matrix[IDB, IDA] = value

	def addEdge(self, A, B, index):
		""" Add new edge to mapping """
		self[(A,B,index)] = self.length
		self[(B,A,index)] = self.length
		if index == -1:
			self._addMatrixEdge(A.ID, B.ID, self.length)
		self.length += 1

	def getBond(self, A, B):
		""" Return index assigned to bond edge """
		return self.matrix[A.ID,B.ID]

	def getEdge(self, A, B, index):
		""" Return index assigned to edge (-1 => bond, >=0 => indexed segment edge) """
		try:
			return self[(A,B,index)]
		except KeyError:
			self.addEdge(A, B, index)
			return self.length - 1

	def nullVector(self):
		""" Returns a vector of 0 of the same length as the dimension of the module """
		return np.zeros(self.length)

	def falseVector(self):
		""" Returns a vector of False of the same length as the dimension of the module """
		return self.nullVector() == 1

	##############################################
	## Constructing a mapping from a module
	##############################################

	def _prepareEdgeMapping(self, tonode, fromnode):
		if fromnode <= tonode:
			self.addEdge(tonode, fromnode, -1)

	def _prepareSegmentEdgeMapping(self, index, tonode, fromnode):
		if fromnode <= tonode:
			self.addEdge(tonode, fromnode, index)

	def _prepareNodeMapping(self, node, module):
		map(lambda N: self._prepareEdgeMapping(N, node), module[node].edges.keys())
		map(lambda N: self._prepareSegmentEdgeMapping(N, node, module[node].twin), range(len(module.segments[node])))

	def _bondVector_Edge(self, vector, edge):
		if edge[2] == -1:
			vector[self[edge]] = True
		return vector

	def _bondVector(self):
		""" Returns a vector which determines which indices correspond to bonds """
		return reduce(lambda V,E: self._bondVector_Edge(V, E), self, self.falseVector())

	def __init__(self, module):
		super(Mapping, self).__init__()
		self.length = 0
		self.maxNodeID = 0
		self.matrix = None
		map(lambda X: self._prepareNodeMapping(X, module), module.nodes())

	##############################################
	## Stats
	##############################################

	def _edgeStr(self, pair):
		return "\t".join([str(pair[0]), str(pair[1]), str(pair[2]), str(self.getEdge(pair[0], pair[1], pair[2]))])

	def __str__(self):
		""" Returns string represenation of the the mappings of edges to indices """
		return "\n".join(map(lambda X: self._edgeStr(X), self))

	##############################################
	## Producing the vector associated to a cycle
	##############################################

	def _updateVector(self, vector, edge, ratio):
		if edge.value * ratio > 0:
			vector[self.getEdge(edge.start, edge.finish, edge.index)] += 1
		else:
			vector[self.getEdge(edge.start, edge.finish, edge.index)] -= 1
		return vector

	def vector(self, cycle):
		""" Returns a vector which corresponds to an Event Cycle """
		# Record all the edges once in case one of them is new, thus affecting the length of nullVector
		map(lambda X: self.getEdge(X.start, X.finish, X.index), cycle)
		return reduce(lambda V,E: self._updateVector(V, E, cycle.value), cycle, self.nullVector())

	##############################################
	## Producing the vector associated to the original genome flow
	##############################################

	def _updateOriginalVector(self, vector, key, module):
		if key[0] >= key[1]:
			return vector
		if key[2] >= 0:
			if module[key[0]].twin == key[1]:
				vector[self[key]] += debug.PLOIDY / module[key[0]].ploidy()
		else:
			if module[key[0]].partner == key[1] and key[0] is not key[1]:
				# This assumes homozygocity of the breakends
				vector[self[key]] -= debug.PLOIDY
		return vector

	def _originalVector(self, module):
		""" Returns a vector with the original genome flow """
		return reduce(lambda V,E: self._updateOriginalVector(V, E, module), self, self.nullVector())

	##############################################
	## Matrix of bonds assigned to each node
	##############################################
	def _overlapMatrix_Node(self, matrix, indexedNode, module):
		row, node = indexedNode
		for other in module[node].edges:
			if other is not node:
				column = self.getEdge(node, other, -1)
				matrix[row, column] = 1
		return matrix

	def _overlapMatrix(self, module):
		""" Returns a matrix with a row for each bond, a column for each edge, and 1 if the edge is an incident bond edge """
		return reduce(lambda M, N: self._overlapMatrix_Node(M, N, module), enumerate(module.keys()), scipy.sparse.lil_matrix((len(module), self.length))).tocsr()

##############################################
## Linear Algebra shorthands
##############################################

def _argmax(vector):
	# Had to reimplement because the numpy argmax returns the first occurence of the max value....
	return map(lambda X: X[0], filter(lambda X: X[1] == max(vector), enumerate(vector)))

def _normalized(vector):
	return vector / np.linalg.norm(vector)

#########################################
## Euclidian histories
#########################################

class EuclidianHistory(overlap.OverlapHistory):
	#########################################
	## Basics
	#########################################
	def __init__(self, module):
		super(EuclidianHistory, self).__init__(module)
		self.mappings = Mapping(module) 
		self.vectors = None
		self.q_base = None
		self.r_triangle = None

	def __copy__(self):
		new = EuclidianHistory(self.module)
		new.copy(self)
		return new

	def copy(self, other):
		super(EuclidianHistory, self).copy(other)
		self.mappings = other.mappings	
		if other.vectors is not None:
			self.vectors = np.array(other.vectors)
			self.q_base = np.array(other.q_base) 
			self.r_triangle = np.array(other.r_triangle)
			eventsMapping = dict(zip(other.events, self.events))
		else:
			self.vectors = None
			self.q_base = None
			self.r_triangle = None

	def embalm(self):
		""" Discards links to matrices which can then be garbage collected (horrible hack, I know) """
		super(EuclidianHistory, self).embalm
		self.mappings = None
		self.vectors = None
		self.q_base = None
		self.r_triangle = None

	#########################################
	## Access to vectors
	#########################################

	def originalVector(self):
		""" Returns the genome flow of the original genome in vector notation """
		return self.mappings._originalVector(self.module)

	def overlapMatrix(self):
		""" Returns the incidence matrix of nodes x bond edges """
		return self.mappings._overlapMatrix(self.module)

	def eventVector(self, event):
		""" Returns the vector assigned to an event Cycle """
		## If new dimensions were added to the mapping
		if self.mappings.length > np.size(self.vectors, 0):
			newDims = self.mappings.length - np.size(self.vectors, 0)
			self.vectors = np.concatenate((self.vectors, np.zeros((newDims, np.size(self.vectors, axis=1)))), axis=0)
			self.q_base = np.concatenate((self.q_base, np.zeros((newDims, np.size(self.vectors, axis=1)))), axis=0)
		return self.vectors[:,self.eventIndex[event]]

	def bondIndices(self):
		""" Returns vector which indicates which indices correspond to bond edges """
		return self.mappings._bondVector()

	#########################################
	## Linear algebra
	#########################################
	def pop(self, event):
		""" Remove event from history """
		i = self.eventIndex[event]
		self.vectors = np.delete(self.vectors, i, 1)
		if np.size(self.vectors) == 0:
			self.vectors = None
			self.q_base = None
			self.r_triangle = None
		elif i != np.shape(self.q_base)[1] - 1:
			self.q_base, self.r_triangle = np.linalg.qr(self.vectors, mode='full')
		else:
			self.q_base = np.delete(self.q_base, i, 1)
			self.r_triangle = np.delete(self.r_triangle, i, 1)
			self.r_triangle = np.delete(self.r_triangle, i, 0)
		super(EuclidianHistory, self).pop(event)

	def _ratios(self):
		return np.array([X.cycle.value for X in self.events])

        def _imposeNewRatios(self, oldRatios, ratios, cactusHistory):
        	# Correct values
                map(lambda X: X[0].cycle.setRatio(X[1]), zip(self.events, ratios))

                # Filter values
                badeggs = filter(lambda X: abs(X.cycle.value) <= ROUNDING_ERROR, self.events)
		indexes = [self.eventIndex[X] for X in badeggs]

                for X in badeggs:
                	cactusHistory.pop(self, X)

               	cactusHistory.correctSchedulingErrors() 

	def absorbEvent(self, cactusHistory, event):
		""" Add event to history """
		## Represent as Euclidian vector
		vector = self.mappings.vector(event.cycle)

		## If empty matrix
		if self.vectors is None:
			self.vectors = np.array([vector]).T
			self.q_base = np.array([_normalized(vector)]).T
			self.r_triangle = np.array([[np.linalg.norm(vector)]])
			super(EuclidianHistory, self).absorbEvent(event)
			return

		## If new dimensions were added to the mapping
		if np.size(vector) > np.size(self.vectors, 0):
			newDims = np.size(vector) - np.size(self.vectors, 0)
			self.vectors = np.concatenate((self.vectors, np.zeros((newDims, np.size(self.vectors, axis=1)))), axis=0)
			self.q_base = np.concatenate((self.q_base, np.zeros((newDims, np.size(self.vectors, axis=1)))), axis=0)

		## Project onto pre-existing vectors base
		projections = self.q_base.T.dot(vector)
		residual = vector - self.q_base.dot(projections)

		if np.linalg.norm(residual) < ROUNDING_ERROR * np.linalg.norm(vector): 
			## Linear combination! 
			## One or more vectors need to be removed...
			weights = np.append(np.linalg.solve(self.r_triangle, projections), [-1])

			# Add the vector to the mix
			self.vectors = np.append(self.vectors, np.array([vector]).T, axis=1)
			self.q_base, self.r_triangle = np.linalg.qr(self.vectors, mode='full')
			super(EuclidianHistory, self).absorbEvent(event)

			## Choose which vector(s) should be eliminated
			ratios = self._ratios()
			rhos = weights / ratios
			lameDuckIndex = random.choice(_argmax(np.absolute(rhos)))

			## Operate the shift
			maxRho = -rhos[lameDuckIndex]
			corrections = weights / maxRho
			newRatios = ratios + corrections
			self._imposeNewRatios(ratios, newRatios, cactusHistory)
		else:
			## Independent vectors, congratulations!
			normal = _normalized(residual)
			self.vectors = np.append(self.vectors, np.array([vector]).T, axis=1)
			self.q_base = np.append(self.q_base, np.array([normal]).T, axis=1)
			self.r_triangle = np.append(self.r_triangle, np.zeros([1,np.shape(self.r_triangle)[1]]), axis=0)
			projections = np.append(projections, normal.dot(vector))
			self.r_triangle = np.append(self.r_triangle, np.array([projections]).T, axis=1)
			super(EuclidianHistory, self).absorbEvent(event)

	###############################################
	## Validation
	###############################################
	def _isOrthonormal(self, matrix):
		n = np.shape(matrix)[1]
		if n > 0:
			assert np.sum(matrix.T.dot(matrix) - np.identity(n)) / (n*n) < 1e-5
		return True

	def _isTriangular(self, matrix):
		assert np.sum(np.tril(matrix, k=-1)) < 1e-5
		return True	

	def _isLinearlyIndependent(self):
		n = len(self.events)
		if n > 0:
			assert self._isOrthonormal(self.q_base)
			assert self._isTriangular(self.r_triangle)
			if not all(x != 0 for x in np.diag(self.r_triangle)):
				print self.r_triangle
			assert all(x != 0 for x in np.diag(self.r_triangle))
			if np.sum(self.q_base.dot(self.r_triangle) - self.vectors) / (n*n) >= 1e-5:
				print self.vectors
				print '>>>>>>>>>>>>'
				print self.q_base.dot(self.r_triangle)
				print '>>>>>>>>>>>>'
				print self.q_base.dot(self.r_triangle) - self.vectors
			assert np.sum(self.q_base.dot(self.r_triangle) - self.vectors) / (n*n) < 1e-5
			assert np.linalg.det(self.r_triangle) != 0
		return True

	def validate(self):
		""" Validation function """
		assert super(EuclidianHistory, self).validate() 
		assert self._isLinearlyIndependent()
		if self.vectors is not None:
			assert len(self.eventIndex) == np.shape(self.vectors)[1]
		return True

###############################################
## Unit test
###############################################

def main():
	pass

if __name__ == "__main__":
	main()
