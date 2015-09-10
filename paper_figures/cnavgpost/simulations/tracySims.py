#!/usr/bin/env python

import cnavg.simulator.simulator
import cnavg.cactus.oriented
import cnavg.cactus.graph
import cPickle as pickle
import argparse 

parser=argparse.ArgumentParser()
parser.add_argument('blocks', help='the number of blocks in the genome', type=int, default=100)
parser.add_argument('events', help='the number of rearrangments to put in the history', type=int, default=10)
args=parser.parse_args()

H = cnavg.simulator.simulator.RandomCancerHistory(args.blocks, args.events, indel_length=10) #20, 15)
A = H.avg()
C = cnavg.cactus.graph.Cactus(A)
G = cnavg.cactus.oriented.OrientedCactus(C)

file = open('true.braney', 'w')
file.write(H.braneyText(G))
file.close()

file = open('CACTUS', 'w')
pickle.dump(G, file)
file.close()
