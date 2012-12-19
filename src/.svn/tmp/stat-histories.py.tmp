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
import math
import operator
import gzip

def main():
    args = parse_args()

    # Parse regions file
    regions = parse_regions( args.regions_file )
    
    # Analyze regions in histories
    region_results = parse_histories( regions, args.history_files, args.buffer )

    # Summarize results & print
    summarize_results( region_results )

def parse_regions( regions_file ):
    regions = collections.defaultdict( list )

    for line in open( regions_file, 'r' ):
        # Ignore comment lines
        if line.startswith('#'): continue

        region = Region( line )
        regions[region.ref].append( region )

    return regions

def parse_histories( regions, history_files, buffer ):
    regs_w_events = {}

    for history_file in history_files:
        for line in gzip.open( history_file, 'r' ):
            # Ignore adjacencies for now
            if line.startswith('A'): continue

            # Parse segment
            segment = Segment( line )

            # Store segment for regions with which it overlaps
            for region in regions[segment.ref]:
                if region.stop >= (segment.start - buffer) and \
                   region.start <= (segment.stop + buffer):
                    
                    # Store this event with this region
                    if not region.key in regs_w_events:
                        regs_w_events[region.key] = region

                    regs_w_events[region.key].segments.append( segment )

    for region in regs_w_events.values():
        region.cycle_ids = region.get_cycle_ids()

    return regs_w_events

def summarize_results( region_results ):
    # Print header
    print '# region\tdup_per_history\tmean_dup_pos\tmean_dup_value\t' + \
          'del_per_history\tmean_del_pos\tmean_del_value\t' + \
          '[gene:prop_shared_cycles:distance]'
    for region in sorted( region_results.itervalues(), \
                          key = lambda r: r.calc_dup_stats()[1] ):

        dup_stats = region.calc_dup_stats()
        del_stats = region.calc_del_stats()
        cycles_summary = region.get_cycles_sum( region_results )

        out_str = ['%s(%s)' % ( region.label, region.key )]

        # Mean dup/history, mean pos, mean value
        # Build dup stats string
        out_str.append( '%s\t%s\t%s' % tuple( map( str, dup_stats ) ) )

        # Build del stats string
        out_str.append( '%s\t%s\t%s' % tuple( map( str, del_stats ) ) )

        # Print all genes w/ > 0 shared cycles
        for gene_cycle in cycles_summary[:3]:
            if gene_cycle[1] > 0:
                out_str.append( '%s:%.4f:%s' % gene_cycle )

        print '\t'.join( out_str )

def parse_args():
    ''' Parses & returns command-line arguments '''
    parser = argparse.ArgumentParser()

    parser.add_argument( '-b', '--buffer', default = 0, type = int )
    parser.add_argument( '-r', '--regions-file' )
    parser.add_argument( '-a', '--adjacencies-file' )
    parser.add_argument( 'history_files', nargs = '+' )

    return parser.parse_args()

class Region:
    def __init__( self, region_line ):
        line_list = region_line.rstrip().split('\t')

        self.label = line_list[0]
        self.ref = line_list[1]
        self.orientation = line_list[2]
        self.start = int( line_list[3] )
        self.stop = int( line_list[4] )

        # Unique key for referring to regions in dictionaries
        self.key = '%s:%d-%d' % ( self.ref, self.start, self.stop )
        
        # Segment list is container for any segments that affect this region
        self.segments = []

        self.dup_stats = None
        self.del_stats = None

	self.cycle_ids = None

    def calc_dup_stats( self ):
        # Cache stats
        if self.dup_stats == None:
            # Calculate stats per history in order to weight by complexity
            histories = {}
            for segment in self.segments:
                # Dups have negative value
                if segment.value < 0:
                    # Initialize history if it's not already in dict
                    if not segment.history_id in histories:
                        histories[segment.history_id] = \
                            History( segment.history_id, \
                                     segment.history_complexity )

                    history = histories[segment.history_id]

                    history.dup_count += 1
                    history.dup_temp_positions.append( segment.temporal_pos )
                    history.dup_values.append( segment.value )

            # Generate weighted sums for mean calculations
            dup_sum = 0.0
            temp_sum = 0.0
            value_sum = 0.0
            scaler_sum = 0.0

            for history in histories.itervalues():
                scaler = math.exp( -history.complexity )
                scaler_sum += scaler

                dup_sum += history.dup_count * scaler
                temp_sum += min( history.dup_temp_positions ) * scaler
                value_sum += sum( history.dup_values ) * scaler

            # Calculate means
            history_count = len( histories )

            if history_count > 0:
                dup_mean = dup_sum / (scaler_sum)
                temp_mean = temp_sum / (scaler_sum)
                value_mean = -value_sum / (scaler_sum)
            else:
                dup_mean, temp_mean, value_mean = [0,1000,0]

            self.dup_stats = ( dup_mean, temp_mean, value_mean )

        return self.dup_stats

    def calc_del_stats( self ):
        # Cache stats
        if self.del_stats == None:
            # Calculate stats per history in order to weight by complexity
            histories = {}
            for segment in self.segments:
                # Dels have positive value
                if segment.value > 0:
                    # Initialize history if it's not already in dict
                    if not segment.history_id in histories:
                        histories[segment.history_id] = \
                            History( segment.history_id, \
                                     segment.history_complexity )

                    history = histories[segment.history_id]

                    history.del_count += 1
                    history.del_temp_positions.append( segment.temporal_pos )
                    history.del_values.append( segment.value )

            # Generate weighted sums for mean calculations
            del_sum = 0.0
            temp_sum = 0.0
            value_sum = 0.0
            scaler_sum = 0.0

            for history in histories.itervalues():
                scaler = math.exp( -history.complexity )
                scaler_sum += scaler

                del_sum += history.del_count * scaler
                temp_sum += min( history.del_temp_positions ) * scaler
                value_sum += sum( history.del_values ) * scaler

            # Calculate means
            history_count = len( histories )

            if history_count > 0:
                del_mean = del_sum / (scaler_sum)
                temp_mean = temp_sum / (scaler_sum)
                value_mean = value_sum / (scaler_sum)
            else:
                del_mean, temp_mean, value_mean = [0,1000,0]

            self.del_stats = ( del_mean, temp_mean, value_mean )

        return self.del_stats

    def get_cycles_sum( self, region_results ):
        # Calculate shared cycles
        cycle_summary = set()
        these_cycles = self.cycle_ids
        for other_region in region_results.itervalues():
            # Ignore other regions with this one's label
            if other_region.label == self.label: continue 
            those_cycles = other_region.cycle_ids
            shared_cycles = these_cycles & those_cycles
            proportion = len(shared_cycles)/ float(len( these_cycles ))

            # Calculate distance of these two regions
            if self.ref == other_region.ref:
                distance = str( self.start - other_region.start )
            else:
                distance = 'N/A'

            cycle_summary.add( ( other_region.label, proportion, distance ) )

        # Sort cycle summary & convert to list
        cycle_summary = sorted( list(cycle_summary), \
                        key = lambda s: s[1], \
                        reverse = True )

        return cycle_summary

    def get_cycle_ids( self ):
        return set((segment.history_id, segment.net_id, segment.cycle_id) for segment in self.segments)

class Segment:
    def __init__( self, segment_line ):
        line_list = segment_line.rstrip().split('\t')

        self.ref = line_list[0]
        self.start = int( line_list[1] )
        self.stop = int( line_list[2] )
        self.value = float( line_list[3] )
        self.history_id = int( line_list[4] )
        self.net_id = int( line_list[5] )
        self.cycle_id = int( line_list[6] )
        self.edge_id = int( line_list[7] )
        self.temporal_pos = int( line_list[8] )
        self.history_complexity = int( line_list[9] )

class Adjacency:
    def __init__( self, adjacency_line ):
        line_list = adjacency_line.rstrip().split('\t')

        self.ref1 = line_list[1]
        self.pos1 = int( line_list[2] )
        self.orientation1 = line_list[3]
        self.ref2 = line_list[4]
        self.pos2 = int( line_list[5] )
        self.orientation2 = line_list[6]
        self.value = float( line_list[7] )
        self.history_id = int( line_list[8] )
        self.net_id = int( line_list[9] )
        self.cycle_id = int( line_list[10] )
        self.edge_id = int( line_list[11] )
        self.temporal_pos = int( line_list[12] )
        self.history_complexity = int( line_list[13] )

class History:
    def __init__( self, id, complexity ):
        self.id = id
        self.complexity = complexity

        self.dup_count = 0
        self.dup_values = []
        self.dup_temp_positions = []
        
        self.del_count = 0
        self.del_values = []
        self.del_temp_positions = []

main()
