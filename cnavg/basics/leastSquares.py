# Copyright (c) 2013, Daniel Zerbino
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

"""Definition of copy-number balancing across a sequence graph"""
 
import numpy as np

def solve(x, precision_x, matrix, y, precision_y):
	"""Direct resolution of the problem using weighted least squares algorithm"""
	"""See http://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares"""
	# Likelihoods: X = x +/- sigma_x = 0 +/- sigma_x
	# Constraints: A.X = y +/- sigma_y
	# precision = 1 / variance = 1 / sigma^2

	# Weighted least squares, we merge everything into one equation
	M = np.vstack([np.identity(len(x)), matrix])
	Y = np.concatenate((x, y))
	W = np.concatenate((precision_x, precision_y))
	
	# We now have: M^T.diag(W).M.X = M^T.W.Y, we factorize the scaling matrix W
	# Note: W should be a diagonal matrix, but I only convert it after 
	# Computing the square root for performance reasons 
	w = np.diag(np.sqrt(W))
	M_prime = np.dot(w, M)
	Y_prime = np.dot(w, Y)

	# We now have: M'^T.M.X = M'^T.Y', we proceed to solving the damn thing
	return np.linalg.solve(np.dot(M_prime.T, M_prime), np.dot(M_prime.T, Y_prime))
