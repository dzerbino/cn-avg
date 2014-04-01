#!/usr/bin/env python

import cnavg.simulator.simulator
import cnavg.cactus.oriented
import cnavg.cactus.graph
import cPickle as pickle

H = cnavg.simulator.simulator.RandomIntegerHistory(20, 20, indel_length = 10)
A = H.avg()
C = cnavg.cactus.graph.Cactus(A)
G = cnavg.cactus.oriented.OrientedCactus(C)

file = open('true.braney', 'w')
file.write(H.braneyText(G))
file.close()

file = open('CACTUS', 'w')
pickle.dump(G, file)
file.close()
