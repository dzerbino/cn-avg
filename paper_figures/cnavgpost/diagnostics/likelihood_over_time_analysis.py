#!/usr/bin/env python 

import argparse
import numpy as np
import os, sys 
import cPickle as pickle 
import cnavgpost.mergehistories.event_cycles_module as histseg

def likelihood_over_time_analysis(events, historyScores, outputdir, filename, pvalcutoff=0.5, simulation=False, testks=True, testns=True): 
	runlen=max(np.fmod(historyScores[:,0], histseg.Global_BINWIDTH))+1
	numruns=historyScores.shape[0]/runlen
	historycosts=np.mean(historyScores[:,1:3], axis=1).reshape(runlen, numruns, order='F')
	#costsarray will have dimensions (runlen)x(numruns)x(# events) 
	costsarray=np.ones((runlen, numruns, len(events)), dtype=int)*-1
	for i in xrange(len(events)):	
		event=events[i]
		ecosts=make_event_costs_array(event, historyScores, numruns, runlen)
		costsarray[:,:,i]=ecosts
	
	if simulation: 
		historycosts=historycosts[:,1:]
		costsarray=costsarray[:,1:,:]	
	# test out different ks
	if testks: 
		myks=(0, 0.5, 1.0)
		get_data_for_diff_ks(historycosts, costsarray, myks, filename, pvalcutoff, outputdir)
	if testns: 
		myk=histseg.Global_K
		get_data_for_diff_ns(historycosts, costsarray, myk, filename, pvalcutoff, outputdir)

def get_data_for_diff_ks(historycosts, costsarray, myks, filename, pvalcutoff, outputdir): 
	kcounts=np.array([])
	for k in myks: 
		#likelihood_array with have dimensions (runlen)x(events)
		likelihood_array=get_likelihood_over_time(historycosts, costsarray, k)	
		np.savetxt(os.path.join(outputdir, "lscores_%f.txt.gz" % k), likelihood_array[likelihood_array.shape[0]-1,:])
		# eventcounts will have dim (runlen)x 3.  One column for all events, one for events with cost 0 and one for events with cost>0
		eventcounts=get_event_counts_over_time(likelihood_array, costsarray, pvalcutoff)
		if kcounts.size==0: 
			kcounts=eventcounts
		else: 
			kcounts=np.dstack((kcounts, eventcounts))
	headline="\t".join(map(str, myks))
	np.savetxt(os.path.join(outputdir, "counts_%s_ks.dat" % filename), kcounts[0,:,:], header=headline, fmt='%d')
	np.savetxt(os.path.join(outputdir, "nonzero_%s_ks.dat" % filename), kcounts[1,:,:], header=headline, fmt='%d')
	np.savetxt(os.path.join(outputdir, "zero_%s_ks.dat" % filename), kcounts[2,:,:], header=headline, fmt='%d')


def get_data_for_diff_ns(historycosts, costsarray, k, filename, pvalcutoff, outputdir):
	#test out different numbers of runs 
	ncounts=np.array([])
	numruns=costsarray.shape[1]
	for n in xrange(1,(numruns+1)):
		likelihood_array=get_likelihood_over_time(historycosts[:,:n], costsarray[:,:n,:], k)
		eventcounts=get_event_counts_over_time(likelihood_array, costsarray[:,:n,:], pvalcutoff)
		if ncounts.size==0: 
			ncounts=eventcounts
		else: 
			ncounts=np.dstack((ncounts, eventcounts))
	headline="\t".join(map(str, range(numruns)))
	np.savetxt(os.path.join(outputdir, "counts_%s_ns.dat" % filename), ncounts[0,:,:], header=headline, fmt='%d')
	np.savetxt(os.path.join(outputdir, "nonzero_%s_ns.dat" % filename), ncounts[1,:,:], header=headline, fmt='%d')
	np.savetxt(os.path.join(outputdir, "zero_%s_ns.dat" % filename), ncounts[2,:,:], header=headline, fmt='%d')

def get_event_counts_over_time(likelihood_array, costsarray, pvalcutoff):
	lcutoff=pvalcutoff 
	numtotal=np.sum(likelihood_array>lcutoff, axis=1)
	ctmp=np.amax(costsarray, axis=1)
	ctmp2=np.maximum.accumulate(ctmp, axis=0)
	nonzero=np.sum((likelihood_array>lcutoff) & (ctmp2>0), axis=1)
	zeros=np.sum((likelihood_array>lcutoff) & (ctmp2==0), axis=1)
	return(np.vstack((numtotal, nonzero, zeros)))

def get_likelihood_over_time(historycosts, costsarray, k):
	hrunscores=np.sum(np.exp(-1*k*historycosts), axis=1)
	historyscores=np.cumsum(hrunscores, axis=0)
	eventscores=np.zeros((historyscores.shape[0],costsarray.shape[2]))
	for n in xrange(costsarray.shape[2]): 
		maskedcosts=np.ma.masked_array(historycosts[:,:costsarray.shape[1]], mask=(costsarray[:,:,n]==-1))
		erunscores=np.sum(np.exp(-1*k*maskedcosts), axis=1).filled(0)
		escores=np.cumsum(erunscores, axis=0)
		eventscores[:,n]=escores/historyscores
	return eventscores
		
	
def make_event_costs_array(event, historyScores, numruns, runlen): 
#	if event.histories==[]: 
#		event.histories=histseg.listout_ranges(event.histRanges)
	hids=np.array(event.histories, dtype=int)
	his=histseg.historyids_to_indices(hids, historyScores)
	eventcosts=np.ones((runlen*numruns))* -1
	eventcosts[his]=np.mean(np.vstack((event.uppercosts, event.lowercosts)), axis=0)
	return eventcosts.reshape((runlen, numruns),order='F' )

def main(pevntsfile, historystats, outputdir, filename, pvalcutoff, simulation): 
	events=pickle.load(open(pevntsfile, 'rb'))
	historyScores=np.loadtxt(historystats, dtype=int)
	if not os.path.exists(outputdir): 
		os.makedirs(outputdir)
	likelihood_over_time_analysis(events, historyScores, outputdir, filename, pvalcutoff, simulation)	


if __name__ == '__main__': 
	parser=argparse.ArgumentParser(description='given a .pevnts file, will do analysis of how many of those events exist and what their likelihood is over the iterations') 
	parser.add_argument('--pevntsfile', help=' a *.pevnts file')
	parser.add_argument('--historystats', help='a historystats.txt file')
	parser.add_argument('--outputdir', help='where to put the output .dat files', default="./")
	parser.add_argument('--filename', help='the basename for the output files', default="evnts")
	parser.add_argument('--pvalcutoff', help='The likelihood score for what counts as true or not', default=0.5, type=float)
	parser.add_argument('--simulation', help='whether it is simulated data or not (if yes, it will ignore run 0.', action='store_true')
	args=parser.parse_args()
	main(args.pevntsfile, args.historystats, args.outputdir, args.filename, args.pvalcutoff, args.simulation)	
