#!/usr/bin/env python
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

import sys

lengths = dict()
file = open(sys.argv[1])
for line in file:
	items = line.strip().split()
	lengths[items[0]] = int(items[1])
file.close()

index = sys.argv[2]

for line in sys.stdin:
    items = line.strip().split()
    if len(items) == 0:
	    continue
    if items[0] == "A":
	    if items[9] != index:
		    continue

	    val = str(abs(float(items[7])))[:5]
	    items[7] = float(items[7])
	    if abs(items[7]) > 1:
		    items[7] = items[7] / abs(items[7])
	    items[7] = int(items[7] * 500 + 500)

	    if items[1] == "None":
		    items[1] = "chr1"
            else:
		    items[1] = items[1]
	    if items[1][:3] != "chr":
		    items[1] = "chr" + items[1]

	    if items[4] == "None":
		    items[4] = "chr1"
	    else:
		    items[4] = items[4]
	    if items[4][:3] != "chr":
		    items[4] = "chr" + items[4]

	    items[2] = int(items[2])
	    items[5] = int(items[5])
	    if items[2] <= 0:
		    items[2] += 1
		    items[5] += 1

	    if items[1] not in lengths:
		    continue

            if items[2] >= lengths[items[1]]:
		    items[2] -= items[2] - lengths[items[1]] + 1
		    items[5] -= items[2] - lengths[items[1]] + 1
            if items[5] >= lengths[items[1]]:
		    items[2] -= items[5] - lengths[items[1]] + 1
		    items[5] -= items[5] - lengths[items[1]] + 1
	    if items[2] <= 0:
		    items[2] = 1

	    name = "_".join([val] + items[9:12])
	    
	    if items[1] == items[4] and items[2] < items[5]:
		    print "\t".join(map(str, [items[1],items[2],items[5],items[4],items[7],items[3],items[6],name,items[13]]))
	    elif items[1] != items[4]:
		    print "\t".join(map(str, [items[1],items[2],items[2],items[4],items[7],items[3],items[6],name,items[13]]))

    elif len(items) > 0: 
    	    if items[5] != index:
		    continue

	    val = str(abs(float(items[3])))[:5]
	    items[3] = float(items[3])
	    if abs(items[3]) > 1:
		    items[3] = items[3] / abs(items[3])
	    items[3] = int(-items[3] * 500 + 500)

	    if items[0][:3] != "chr":
		    items[0] = "chr" + items[0]

	    items[1] = int(items[1])
	    items[2] = int(items[2])
	    if items[1] <= 0:
		    items[1] += 1
		    items[2] += 1

            if items[2] >= lengths[items[0]]:
		    items[2] -= items[2] - lengths[items[0]] + 1
		    items[1] -= items[2] - lengths[items[0]] + 1
            if items[1] >= lengths[items[0]]:
		    items[2] -= items[1] - lengths[items[0]] + 1
		    items[1] -= items[1] - lengths[items[0]] + 1
	    if items[1] <= 0:
		    items[1] = 1

	    name = "_".join([val] + items[5:8])

	    print "\t".join(map(str, [items[0],items[1],items[2],"SEGMENT",items[3],"+","+",name,items[9]]))
