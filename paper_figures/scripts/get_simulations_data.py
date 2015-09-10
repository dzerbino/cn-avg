#!/usr/bin/env python 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse 
import os, sys
import cnavgpost.mergehistories.event_cycles_module as histseg
#from IPython.parallel.util import interactive 
#from ../event_cycles_module import Global_EVENTTYPES
Global_types=['any', 'amplification', 'deletion', 'adjacency']
	
def make_SimData_dict(simdir, edgeorevent, pvalcutoff=0, cutcost=False):
	sys.stderr.write("loading %s\n" % simdir)
	mysimdata={}
	datfile=os.path.join(simdir, "%s.dat" % edgeorevent) 
	breaks=os.path.join(simdir, "breakpoints.txt")
	histstats=os.path.join(simdir, "historystats.txt")
	historyScores=np.loadtxt(histstats, dtype=int)
	(runlen, nruns)=histseg.get_runlen(historyScores)
	mysimdata['truescores']=np.mean(historyScores[0,1:3])
	mysimdata['avescores']=np.mean(historyScores[runlen:,1:3])
	mysimdata['minscores']=np.min(np.mean(historyScores[runlen:,1:3], axis=1))
	mysimdata['maxscores']=np.max(historyScores[runlen:,2])
	nodes=np.loadtxt(breaks, usecols=(2,3,4,5), dtype=int)
	mysimdata['truenodes']=sum(nodes[:,1]>0)
	mysimdata['prednodes']=sum(nodes[:,0]>0)
	mysimdata['totalnodes']=nodes.shape[0]
	# get TP, etc. for calculating accuracies, etc. 
	(TP, FP, FN, TN) = (0,0,0,0)
	mydat=np.loadtxt(datfile, skiprows=1, usecols=(2,3,5), ndmin=2)
	(avecost, lscore, t) = range(mydat.shape[1])
	if cutcost:
		Ps=(mydat[:,lscore]>pvalcutoff) & (mydat[:,avecost]>0)
		Ns=(mydat[:,lscore]<=pvalcutoff) | (mydat[:,avecost]==0)
	else: 
		Ps=mydat[:,lscore]>pvalcutoff
		Ns=mydat[:,lscore]<=pvalcutoff
	if pvalcutoff>0: 
		TP = ((mydat[:,t] ==1) | (mydat[:,t]==3)) & Ps
		FP = (mydat[:,t]==0) & Ps
		TN = (Ns & (mydat[:,t] == 0)) | (mydat[:,t]==2)
		FN = (mydat[:,t] ==-1) | ((mydat[:,t] == 1) & Ns)
	else: 
		TP = ((mydat[:,t] ==1) | (mydat[:,t]==3)) 
		FP = mydat[:,t] ==0
		TN = mydat[:,t] == 2 #mydat[:,t] ==0
		FN = mydat[:,t] ==-1
	mysimdata['TP']=sum(TP)
	mysimdata['FP']=sum(FP)
	mysimdata['TN']=sum(TN)
	mysimdata['FN']=sum(FN)
	return mysimdata
	
def get_accuracy_values(simulations, datafn, edgeorevent, pvalcutoff, cutcost):
	mysimulations=open(simulations, 'r').readlines()
	myfmts={'totalnodes':'%d', 'truenodes':'%d', 'truescores':'%f', 'avescores':'%f', 'minscores':'%f', 'TP':'%d', 'FP':'%d', 'TN':'%d', 'FN':'%d', 'simids':'%s', 'blocks':'%s', 'prednodes':'%d', 'maxscores':'%d'}
	mykeys=['simids', 'blocks', 'totalnodes', 'truenodes', 'prednodes', 'truescores', 'avescores', 'minscores', 'maxscores', 'TP', 'FP', 'TN', 'FN']
	fout=open(datafn, 'w')
	fout.write("#" + "\t".join(mykeys) + "\n")
	formats=[]
	for key in mykeys: 
		formats.append(myfmts[key])
	alldata=[]
	for i in xrange(len(mysimulations)): 
		line=mysimulations[i]
		(dir, simid, blocks, events)=line.strip().split('\t')
		mysimdata=make_SimData_dict(dir, edgeorevent, pvalcutoff, cutcost)
		mysimdata['blocks']= blocks
		mysimdata['simids']=simid
		tmpdata=[]	
		for key in mykeys: 
			tmpdata.append(mysimdata[key])
		fout.write("\t".join(map(str, tmpdata)) + "\n")
	fout.close()
	
#@interactive 
def main(simulations, datafn, edgeorevent, pvalcutoff, cutcost):
	get_accuracy_values(simulations, datafn, edgeorevent, pvalcutoff, cutcost)


if __name__ == "__main__": 
	parser=argparse.ArgumentParser(description="makes plots")
	parser.add_argument('simulations', help='A list of simulation directories.  Should have the form: <dir><ID><#blocks><# events>.')
	parser.add_argument('data', help='File to print data to.')
	parser.add_argument('--datatype', help='whether to analyze data for edges, merged edges, or events', choices=['edges', 'mrgedges', 'events'], default='events')
	parser.add_argument('--pvalcutoff', help='the p-value cutoff for what counts as true', default=0, type=float)
	parser.add_argument('--cutcost', help='exclude events with cost 0.', action='store_true')
	args=parser.parse_args()
	main(args.simulations, args.data, args.datatype, args.pvalcutoff, args.cutcost)


