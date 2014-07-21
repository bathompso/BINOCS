# BINOCS threshold test
# Described in section 5.3 of Thompson et al. (in prep)
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

# Make synthetic dataset
single_synth = binocs.makesynth(mag, binary, options)
synth = np.zeros([3*single_synth.shape[0], single_synth.shape[1]])
synthlen = single_synth.shape[0]
synth[0:synthlen,:] = single_synth[:,:]
synth[synthlen:2*synthlen,:] = single_synth[:,:]
synth[2*synthlen:3*synthlen,:] = single_synth[:,:]
for f in range(17):
	synth[synthlen:2*synthlen,2*f] += 0.8
	synth[2*synthlen:3*synthlen,2*f] -= 0.8
	
print("\nCreated synthetic set with %d stars." % (synth.shape[0]))

# Loop through different threshold values
thresholds = range(5,20)
confusion_matrix = np.zeros([len(thresholds), 2])
for t in range(len(thresholds)):
	print("\nWorking on threshold %.1f..." % (thresholds[t]))

	# Run SED fitting on entire synthetic set
	results = binocs.sedfit(singles, binary, synth, options, chicut=thresholds[t])
	summary = binocs.summarize(results, binary, singles)
	
	# Missed percentage
	confusion_matrix[t,0] = len([x for x in range(synthlen) if summary[x,0] == 0]) / synthlen
	
	# Contamination percentage
	confusion_matrix[t,1] = len([x for x in range(2*synthlen) if summary[x+synthlen,0] > 0]) / (2*synthlen)
	
	# Print results to file
	of = open("binocs_thresh.%02d.dat" % (t), 'w')
	for l in range(synth.shape[0]):
		if l < synthlen: member = 1
		else: member = 0
		if summary[l,0] > 0: result = 1
		else: result = 0
		print("%8.3f %8.3f %8.3f %8.3f   %d %d" % (synth[l,2], synth[l,4], synth[l,20], synth[l,24], member, result), file=of)
	of.close()

print("\n")	
for t in range(len(thresholds)):
	print("%5.1f  %5.1f %5.1f " % (thresholds[t], confusion_matrix[t,0]*100, confusion_matrix[t,1]*100))
