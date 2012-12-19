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
#! /inside/home/bjrice/python-2.7.2/python

################################################################################
# Author: Brandon Rice (bjrice@ucsc.edu)
################################################################################

import argparse
import collections

def main():
    args = parse_args()

    # Parse analysis files & collect stats
    genes = parse_analyses( args.analysis_files )

    # Merge gene values and print to stdout
    print "\n".join(map(print_gene, genes))

def print_gene( gene ):
    out_list = [gene]

    for values in genes[gene][:6]:
        out_list.append( str( sum( values ) / len( values ) ) )

    out_list.extend( genes[gene][6] )

    return '\t'.join(out_list)

def parse_analyses( analysis_files ):
    genes = {}
    for file in analysis_files:
        for line in open( file, 'r' ):
            # Skip header
            if line.startswith('#'): continue

            line_list = line.rstrip().split('\t')

            # Initialize this gene's entry if it doens't already exist
            if not line_list[0] in genes:
                genes[line_list[0]] = [[] for i in range(6)]
                genes[line_list[0]].append( set() )

            # Add values to list
            for i in range(6):
                genes[ line_list[0] ][i].append( float(line_list[i+1]) )

            # Add cycle-sharing genes to set
            genes[ line_list[0] ][6].update( line_list[7:] )

    return genes

def parse_args():
    ''' Parses & returns command-line arguments '''
    parser = argparse.ArgumentParser()

    parser.add_argument( 'analysis_files', nargs = '+' )

    return parser.parse_args()

main()
