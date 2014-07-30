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

""" Linear algebra representation of a history """

import math
import random
import numpy as np
import copy

import debug
import overlap

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
	def length(self):
		return len(self.keys)

	def _addMatrixEdge(self, IDA, IDB):
		""" The use of the matrix allows for accelerated bond lookups """
		self.matrix[(IDA, IDB)] = self.length()
		self.matrix[(IDB, IDA)] = self.length()

	def addEdge(self, A, B, index):
		""" Add new edge to mapping """
		self[(A,B,index)] = self.length()
		self[(B,A,index)] = self.length()
		if index == -1:
			self._addMatrixEdge(A.ID, B.ID)
			if B not in self.module[A].edges:
				self.module.addLiftedEdge(A, B, 0)
		self.keys.append((A,B,index))

	def getBond(self, A, B):
		""" Return index assigned to bond edge """
		return self.matrix[(A.ID,B.ID)]

	def getEdge(self, A, B, index):
		""" Return index assigned to edge (-1 => bond, >=0 => indexed segment edge) """
		try:
			return self[(A,B,index)]
		except KeyError:
			self.addEdge(A, B, index)
			return self.length() - 1

	def reverseLookup(self, value):
		return self.keys[value]

	def nullVector(self):
		""" Returns a vector of 0 of the same length as the dimension of the module """
		return np.zeros(self.length())

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

	def _prepareNodeMapping(self, node):
		map(lambda N: self._prepareEdgeMapping(N, node), self.module[node].edges)
		map(lambda i: self._prepareSegmentEdgeMapping(i, node, self.module[node].twin), range(len(self.module[node].segment)))

	def __init__(self, module):
		super(Mapping, self).__init__()
		self.maxNodeID = 0
		self.matrix = dict()
		self.keys = []
		self.module = module
		map(self._prepareNodeMapping, module)

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
		vector[self.getEdge(edge.start, edge.finish, edge.index)] += int(round(edge.value / ratio))
		return vector

	def vector(self, event):
		""" Returns a vector which corresponds to an Event Cycle """
		# Record all the edges once in case one of them is new, thus affecting the length of nullVector
		map(lambda X: self.getEdge(X.start, X.finish, X.index), event.cycle)
		return reduce(lambda V,E: self._updateVector(V, E, event.ratio), event.cycle, self.nullVector())

	##############################################
	## Producing the unitary vector associated to a cycle
	##############################################

	def _updateUnitaryVector(self, vector, edge):
		if edge.value > 0:
			vector[self.getEdge(edge.start, edge.finish, edge.index)] += 1
		else:
			vector[self.getEdge(edge.start, edge.finish, edge.index)] -= 1
		return vector

	def unitaryVector(self, event):
		""" Returns a vector which corresponds to an Event Cycle """
		# Record all the edges once in case one of them is new, thus affecting the length of nullVector
		map(lambda X: self.getEdge(X.start, X.finish, X.index), event.cycle)
		return reduce(lambda V,E: self._updateUnitaryVector(V, E), event.cycle, self.nullVector())

	##############################################
	## Producing the vector associated to the original genome flow
	##############################################

	def _updateOriginalVector(self, vector, key, module):
		if key[0] >= key[1]:
			return vector
		if key[2] >= 0:
			if key[0] in module.pseudotelomeres:
				vector[self[key]] += int(math.ceil(abs(module[key[0]].segment[key[2]])))
			else:
				vector[self[key]] += debug.PLOIDY / module[key[0]].ploidy()
		else:
			# Note conjuguate transformation
			if module[key[0]].partner == key[1]:
				# This assumes homozygocity of the breakends
				vector[self[key]] -= debug.PLOIDY
		return vector

	def _originalVector(self, module):
		""" Returns a vector with the original genome flow """
		return reduce(lambda V,E: self._updateOriginalVector(V, E, module), self, self.nullVector())

	##############################################
	## Producing the vector associated to bond edges
	##############################################

	def _bondVector_Edge(self, vector, edge):
		if edge[2] == -1:
			vector[self[edge]] = True
		return vector

	def _bondVector(self):
		""" Returns a vector which determines which indices correspond to bonds """
		return reduce(lambda V,E: self._bondVector_Edge(V, E), self, self.falseVector())

	##############################################
	## Producing the vector associated to the stub segment
	##############################################

	def _stubVector_Edge(self, vector, edge):
		if edge[0].chr == 'None' and edge[1].chr == 'None':
			vector[self[edge]] = True
		return vector

	def _stubVector(self):
		""" Returns a vector which determines which indices correspond to stubs """
		return reduce(lambda V,E: self._stubVector_Edge(V, E), self, self.falseVector())

	##############################################
	## Producing the vector associated to novel bonds with negative edges
	##############################################

	def _negativeBondsVector_Edge(self, vector, edge, module):
		assert edge[0] in module
		assert edge[1] in module
		if edge[2] == -1:
			assert edge[1] in module[edge[0]].edges
			assert edge[0] in module[edge[1]].edges

		if edge[2] == -1 and module[edge[0]].partner != edge[1] and module[edge[0]].edges[edge[1]] < 0:
			vector[self[edge]] = True
		return vector

	def _negativeBondsVector(self, module):
		""" Returns a vector which determines which indices correspond to negativeBondss """
		return reduce(lambda V,E: self._negativeBondsVector_Edge(V, E, module), self, self.falseVector())


##############################################
## Linear Algebra shorthands
##############################################

def _argmax(vector):
	max = np.amax(vector)
	return np.where(vector == max)[0]

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

	def eventVector(self, event):
		""" Returns the vector assigned to an event Cycle """
		## If new dimensions were added to the mapping
		if self.mappings.length() > np.size(self.vectors, 0):
			newDims = self.mappings.length() - np.size(self.vectors, 0)
			self.vectors = np.concatenate((self.vectors, np.zeros((newDims, np.size(self.vectors, axis=1)))), axis=0)
			self.q_base = np.concatenate((self.q_base, np.zeros((newDims, np.size(self.vectors, axis=1)))), axis=0)
		normalizedVector = self.vectors[:,self.eventIndex[event]]
		values = np.abs(normalizedVector)
		correction = 1.0 / np.min(values[values > 0])
		return np.round(normalizedVector * correction)

	def bondIndices(self):
		""" Returns vector which indicates which indices correspond to bond edges """
		return self.mappings._bondVector()

	def segmentIndices(self):
		""" Returns vector which indicates which indices correspond to segment edges """
		return self.mappings._segmentVector()

	def stubIndices(self):
		""" Returns vector which indicates which indices correspond to stub-stub edges """
		return self.mappings._stubVector()

	def negativeBondsIndices(self):
		""" Returns vector which indicates which indices correspond to negative bond edges """
		return self.mappings._negativeBondsVector(self.module)

	#########################################
	## Linear algebra
	#########################################
	def pop(self, event):
		""" Remove event from history """
		i = self.eventIndex[event]
		self.vectors = np.delete(self.vectors, i, 1)
		if len(np.shape(self.vectors)) == 0:
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
		return np.array([X.cycle[0].value for X in self.events])

        def _imposeNewRatios(self, ratios, cactusHistory, oldratios, rhos, weights, maxRho, corrections):
		badeggs = []
		for event, ratio in zip(self.events, ratios):
			if abs(ratio) <= ROUNDING_ERROR:
				print 'POPPING OUT'
				print event
				print 'MAX RHO', maxRho
				print 'OLD RATIO', oldratios[self.eventIndex[event]]
				print 'RHO', rhos[self.eventIndex[event]]
				print 'WEIGHTS', weights[self.eventIndex[event]]
				print 'CORRECTION', corrections[self.eventIndex[event]]
				badeggs.append(event)
			else:
				if abs(ratio - event.cycle[0].value) > 1e-2:
					print 'REWEIGHTING %f' % ratio
					print event
					print 'MAX RHO', maxRho
					print 'OLD RATIO', oldratios[self.eventIndex[event]]
					print 'RHO', rhos[self.eventIndex[event]]
					print 'WEIGHTS', weights[self.eventIndex[event]]
					print 'CORRECTION', corrections[self.eventIndex[event]]
				event.setRatio(ratio)

                for X in badeggs:
                	cactusHistory.pop(self, X)

               	cactusHistory.correctSchedulingErrors() 

	def absorbEvent(self, cactusHistory, event):
		""" Add event to history """
		## Represent as Euclidian vector
		if self.mappings is None:
			self.mappings = Mapping(self.module)
		vector = self.mappings.unitaryVector(event)

		## If empty matrix
		if self.vectors is None:
			self.vectors = np.array([vector]).T
			self.q_base = np.array([_normalized(vector)]).T
			self.r_triangle = np.array([[np.linalg.norm(vector)]])
			super(EuclidianHistory, self).absorbEvent(event)
			return

		#assert np.size(self.vectors, axis=1) == np.size(self.q_base, axis=1)

		## If new dimensions were added to the mapping
		if np.size(vector) > np.size(self.vectors, 0):
			newDims = np.size(vector) - np.size(self.vectors, 0)
			self.vectors = np.concatenate((self.vectors, np.zeros((newDims, np.size(self.vectors, axis=1)))), axis=0)
			self.q_base = np.concatenate((self.q_base, np.zeros((newDims, np.size(self.q_base, axis=1)))), axis=0)

		## Project onto pre-existing vectors base you can do this on the orthonormal matrix
		projections = self.q_base.T.dot(vector)
		residual = vector - self.q_base.dot(projections)
		if np.linalg.norm(residual) < ROUNDING_ERROR * np.linalg.norm(vector):
			## Note: solving Mx = y is equivalent to solving Rx = Q.Ty
			if not debug.INTEGER_HISTORY:
				solution = np.linalg.solve(self.r_triangle, projections)
			else:
				# If we're working exclusively with integer coefficients
				# then non-independent vectors will be allowed, i.e. 
				# the equation system is not singular. 
				solution = np.linalg.lstsq(self.r_triangle, projections)[0]
			weights = np.append(solution, [-1])

		if np.linalg.norm(residual) < ROUNDING_ERROR * np.linalg.norm(vector) and (not debug.INTEGER_HISTORY or all(abs(X - round(X)) < ROUNDING_ERROR for X in weights)): 
			## Linear combination! 
			## One or more vectors need to be removed...
			if debug.INTEGER_HISTORY:
				weights = map(round, weights)

			# Add the vector to the mix
			self.vectors = np.append(self.vectors, np.array([vector]).T, axis=1)
			self.q_base, self.r_triangle = np.linalg.qr(self.vectors, mode='full')
			super(EuclidianHistory, self).absorbEvent(event)

			## Choose which vector(s) should be eliminated
			ratios = self._ratios()
			rhos = weights / np.abs(ratios)
			lameDuckIndex = random.choice(_argmax(np.absolute(rhos)))

			## Operate the shift
			maxRho = -rhos[lameDuckIndex]
			corrections = weights / maxRho
			# So this is really weird an deserves to be better done, but essentially the ratio 
			# encodes for two things: norm and marker of the sign of the first edge.
			# therefore sometimes we work with only norm, sometimes we heed both.
			newRatios = np.where(ratios > 0, ratios + corrections, ratios - corrections)
			self._imposeNewRatios(newRatios, cactusHistory, ratios, rhos, weights, maxRho, corrections)
			#debug.TRIGGER_TEST = True
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
	G = avg.randomNearEulerianGraph(10)
	B = balanced.BalancedAVG(G)
	C = cactus.Cactus(B)
	N = normalized.NormalizedCactus(C)
	O = oriented.OrientedCactus(N)
	H = cycleCover.initialHistory(O)
	h = random.choice(H.netHistories.values())
	e = random.choice(h.events)
	e1 = copy.copy(e)
	e1.setRatio(e.cycle[0].value / 2)
	e2 = copy.copy(e)
	e2.setRatio(-e.cycle[0].value / 2)
	e3 = copy.copy(e)
	e3.setRatio(e.cycle[0].value / 2)
	e4 = copy.copy(e)
	e4.setRatio(e.cycle[0].value / 2)
	e5 = copy.copy(e)
	e5.setRatio(e.cycle[0].value * 7)
	e6 = copy.copy(e)
	e6.setRatio(-e.cycle[0].value * 7)

	h.validate()
	h.pop(e)
	h.absorbEvent(H, e2)
	h.absorbEvent(H, e1)
	h.absorbEvent(H, e)
	h.validate()
	h.pop(e)
	h.absorbEvent(H, e3)
	h.absorbEvent(H, e4)
	h.validate()
	h.absorbEvent(H, e5)
	h.absorbEvent(H, e6)
	h.validate()

if __name__ == "__main__":
	import cnavg.cactus.graph as cactus
	import cnavg.cactus.oriented as oriented
	import cnavg.cactusSampling.sampling as normalized
	import cnavg.avg.graph as avg
	import cnavg.avg.balanced as balanced
	import cnavg.historySampling.cycleCover as cycleCover
	debug.DEBUG = True
	main()
