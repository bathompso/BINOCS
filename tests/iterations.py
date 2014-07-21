# BINOCS number-of-iterations test
# Described in section 5.1 of Thompson et al. (in prep)
from __future__ import print_function, division
import numpy as np
import sys, binocs

# Read in data from files
options = binocs.readopt((sys.argv)[1])
info, mag = binocs.readdata(options)
oiso = binocs.readiso(options)

# Interpolate isochrone to new mass grid
singles = binocs.minterp(oiso, options['dm'])

# Adjust isochrone to empirical ridgeline, if necessary
singles = binocs.fidiso(singles, options)

# Create binary array
binary = binocs.makebin(singles, options)

# Loop through different number of iterations
nruns = [3, 10, 30, 90, 150, 200, 300, 400, 500, 600, 700, 1200]
sigmas = np.zeros([len(nruns), 4])
for n in range(len(nruns)):
	print("\nWorking on %dx..." % nruns[n]) 
	options['nruns'] = nruns[n]
	iteration_summary = np.zeros([mag.shape[0], 2, 5])
	# Complete all runs for the iteration amount
	for r in range(iteration_summary.shape[2]):
		results = binocs.sedfit(singles, binary, mag, options)
		summary = binocs.summarize(results, binary, singles)
		iteration_summary[:,0,r] = summary[:,0]
		iteration_summary[:,1,r] = summary[:,1]
	# Compute Sigma
	primary_sigma, secondary_sigma = np.zeros(mag.shape[0]), np.zeros(mag.shape[0])
	for m in range(mag.shape[0]):
		valid = [x for x in range(iteration_summary.shape[2]) if iteration_summary[m,0,x] > 0]
		if len(valid) == 0:
			primary_sigma[m], secondary_sigma[m] = -1, -1
			continue
		elif len(valid) == 1:
			primary_sigma[m], secondary_sigma[m] = 0, 0
			continue
		primary_sigma[m] = np.std(iteration_summary[m,0,:]) / np.mean([iteration_summary[m,0,x] for x in valid])
		if np.mean([iteration_summary[m,1,x] for x in valid]) > 0:
			secondary_sigma[m] = np.std(iteration_summary[m,1,:]) / np.mean([iteration_summary[m,1,x] for x in valid])
		
	# Select only valid results
	primary_sigma, secondary_sigma = primary_sigma[primary_sigma > -0.5], secondary_sigma[secondary_sigma > -0.5]
	primary_sigma, secondary_sigma = primary_sigma[primary_sigma < 1.1], secondary_sigma[secondary_sigma < 1.1]
	
	# Sort sigmas
	primary_sigma, secondary_sigma = np.sort(primary_sigma), np.sort(secondary_sigma)

	sigmas[n,0], sigmas[n,1] = np.median(primary_sigma), np.median(secondary_sigma)
	sigmas[n,2] = primary_sigma[int(0.90 * len(primary_sigma))]
	sigmas[n,3] = secondary_sigma[int(0.90 * len(secondary_sigma))]
	
for n in range(len(nruns)):
	print("%3d  %.3f %.3f %.3f %.3f" % (nruns[n], sigmas[n,0], sigmas[n,1], sigmas[n,2], sigmas[n,3]))