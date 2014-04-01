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
Python wrapper for GaBP application in GraphLab
By Daniel Zerbino, based on Matlab code by Danny Bickson
"""

import sys
import os
import tempfile
import struct
import numpy as np
import subprocess

#############################################
# Convenience binary writer/reader functions:
#############################################

def writeDouble(F, x):
     F.write(struct.pack('<d', x))

def writeVector(F, vector):
     for X in vector: 
	writeDouble(F, X)

def writeInt(F, x):
     F.write(struct.pack('<i', x))

def readInt(F):
     return struct.unpack('<i', F.read(struct.calcsize('i')))[0]

def readDouble(F):
     return struct.unpack('<d', F.read(struct.calcsize('d')))[0]

def readVector(F, n):
     return [readDouble(F) for i in range(n)]

#############################################
# Convenience MatLab like functions
#############################################

def enumerate(M):
     for i in range(np.shape(M)[0]):
	for j in range(np.shape(M)[1]):
	   # Beware of Matlab numbering!!
	   yield i+1, j+1, M[i,j]

def find(M):
     return filter(lambda X: X[2] != 0, enumerate(M)) 


def save_c_gl(fn, A, y, x=None, sigma_y=None, sigma_x=None, square=False):
    """
    Script for exporting system of linear equations of the type Ax = y to graphlab format.
    Input: fn - output file name
    A - A mxn matrix (if A is square than m=n)
    y - mx1 observation vector
    x - nx1 known solution (optional, if not given will write a vector of zeros)
    sigma - a vector m+nx1 of noise levels (optional, for non-square matrices only)
    """
    m, n = np.shape(A)
    if len(y) != m:
       sys.exit('y vector should be of len as the number of A rows (%i and %i resp.)' % (len(y), m))

    if x is not None and len(x) != n:
       sys.exit('x vector length should be as the matrix A columns')

    if sigma_y is None and sigma_x is not None:
       sys.exit("sigma_y should be provided when entering sigma_x")

    if sigma_x is None and sigma_y is not None:
       sys.exit("sigma_x should be provided when entering sigma_y")

    if sigma_x is not None:
       if square:
          sys.exit('sigma noise level input is allowed only for non-square matrices')
       else:
          if len(sigma_x) != n:
             sys.exit('sigma_x length should be number of cols of A')
          if len(sigma_y) != m:
             sys.exit('sigma_y length should be number of rows of A')

    if square:
        # matrix is square, edges are non digonal entries
        vals = find((A - np.diag(np.diag(A))))
	print "Saving a square matrix A"
    else:
        # matrix is not square, edges are non zero values
        vals = find(A)
	print "Saving a non-square matrix A"

    F = open(fn, 'wb')

    #write matrix size
    writeInt(F, m)
    writeInt(F, n)

   # write y (the observation), x (the solution, if known), diag(A) the
   # variance (if known, else default variance of 1)
    writeVector(F, y)
    if x is not None:
	writeVector(F, x)
    else:
	writeVector(F, (0 for x in range(n)))

    if square:
	writeVector(F, diag(A))
    else:
        if sigma_y is not None:
	    writeVector(F, sigma_y)
	    writeVector(F, sigma_x)
        else:
            writeVector(F, (1 for x in range(m+n)))

    #write number of edges
    assert len(vals) > 0
    writeInt(F, len(vals))
    # pad with zeros for 64 bit offset
    writeInt(F, 0)

    if not square:
        offset = m
    else:
        offset = 0

    #write all edges
    for val in vals:
       writeInt(F, val[0])
       writeInt(F, val[1] + offset)
       writeDouble(F, val[2])

    F.close()

    #verify written file header
    F = open(fn,'rb')
    x = readInt(F)
    assert x == m
    F.close()

    print 'Wrote succesfully into file: %s' % fn


def load_c_gl(filename, columns):
    """
    Script for reading the output of the GaBP GraphLab program into matlab
    Returns (x, diag) tuple
    where x = inv(A)*b as computed by GaBP
    and diag = diag(inv(A)) - an approximation to the main diagonal of the inverse matrix of A.
    """
    F = open(filename, 'rb')

    x = readVector(F, columns)
    diag = readVector(F, columns)

    F.close()
    os.remove(filename)
    return x, diag

########################################################
## Wrapper Utility to be used from outside
########################################################

def runGaBP(convergence, A, y, sigma_y=None, x=None, sigma_x=None, square=False):
    """Wrapper function for the GaBP module which resolves a Gaussian propagation problem using GraphLab"""
    file, input = tempfile.mkstemp(dir='.')

    save_c_gl(input, A, y, x=x, sigma_y=sigma_y, sigma_x=sigma_x, square=square)

    args = ['gabp', '--data', input, '--threshold', str(convergence), '--algorithm', '0', '--scheduler=round_robin', '--square']
    if not square:
        args.append('false')
    else:
	args.append('true')
    print "Running " + " ".join(args)
    ret = subprocess.Popen(args, stdout=sys.stdout, stderr=subprocess.STDOUT).wait()

    if ret != 0:
	sys.exit("GaBP did not complete")

    os.remove(input)
    x2, diag = load_c_gl(input + ".out", len(x))
    return x2, diag

#########################################################
## Unit test
#########################################################
def main():
	A = np.array([[0.2785, 0.9649],[0.5469, 0.1576],[0.9575, 0.9706]])
        y = np.array([1.2434, 0.7045, 1.9281])
	sigma_y= np.array([1e-10, 1e-10, 1e-10]) 
        x = np.array([0, 0])
	sigma_x = np.array([1, 1])
        convergence = 1e-10
        x2, diag = runGaBP(convergence, A, y, sigma_y=sigma_y, x=x, sigma_x=sigma_x)

	print 'A'
	print A
	print 'y'
	print y
	print 'Initial X'
	print x
	print 'Initial Error'
	print A.dot(x) - y
	print 'Final X'
	print x2
	print 'Final Error'
	print A.dot(x2) - y
	print 'diag'	
	print diag

if __name__=='__main__':
        main()
