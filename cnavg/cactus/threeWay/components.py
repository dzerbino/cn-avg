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

import sys
import subprocess
import tempfile
import os

"""Wrapper for the threeWay executable"""

###########################################
## Input 
###########################################

def connectionNames(connection):
	return map(str, connection)

def groupString(connection):
	return "\t".join(connectionNames(connection)) + "\n"

def write3WayGroupInput(connection, file):
	file.write(groupString(connection))

def write3WayInput(connections, filename):
	"""Writes file to be read by threeWay"""
	file = open(filename, "w")
	for connection in connections:
		write3WayGroupInput(connection, file)
	file.close()
			
###########################################
## Output 
###########################################

def parse3WayOutputLine(line):
	return map(int, line.strip().split())

def parse3WayOutput(filename):
	"""Reads the output of threeWay"""
	file = open(filename)
	result = [parse3WayOutputLine(line) for line in file]
	file.close()
	os.remove(filename)
	return result

###########################################
## Master function
###########################################

def compute(connections):
	"""Wrapper to the threeWay executable, which detects three wat connected components in linear time"""
	file1, input = tempfile.mkstemp(dir='.')
	file2, output = tempfile.mkstemp(dir='.')

	write3WayInput(connections, input)
	print "Running " + " ".join(['3way', input, output])
	if subprocess.Popen(['3way', input, output], stdout=sys.stdout, stderr=subprocess.STDOUT).wait() != 0:
	    sys.exit("3way did not complete")
	print "Done"
	os.remove(input)
	
	return parse3WayOutput(output)

###########################################
## Unit test
###########################################

def main():
	pass

if __name__ == "__main__":
	main()
