# Import modules
from __future__ import print_function, division
import numpy as np
import sys, binocs2
import matplotlib.pyplot as plt
	
	
print("\n\n")
print("===========================================")
print("|                BINARYFIT                |")
print("|             by Ben Thompson             |")
print("===========================================")


# Read in data from files
options = binocs2.readopt((sys.argv)[1])
info, mag = binocs2.readdata(options)
oiso = binocs2.readiso(options)

# Interpolate isochrone to new mass grid
singles = binocs2.minterp(oiso, options['dm'])

# Adjust isochrone to empirical ridgeline, if necessary
singles = binocs2.fidiso(singles, options)

# Create binary array
binary = binocs2.makebin(singles, options)

#### INITIAL BINARY FITTING
# Run SED fitting on all stars in dataset
#print("\nComputing Initial Results.")
#results = binocs2.sedfit(singles, binary, mag, options)

# Compute Initial Results
#summary = binocs2.summarize(results, binary, singles)
#print("\nWriting intial results to '%s--binary.txt'" % (options['data']))
#out = open(options['data']+"--binary.txt", "w")
#for s in range(mag.shape[0]):
#	# Print out star to file
#	outstr = "%16s %9.5f %9.5f " % (info[s][0], info[s][1], info[s][2])
#	for i in range(17): outstr += "%6.3f " % mag[s,2*i]
#	outstr += "   %2d %2d %6.3f %6.4f %6.3f %6.4f %5.2f   %6.3f %6.4f %5.2f" % (summary[s,8], info[s][3], summary[s,0], summary[s,1], summary[s,2], summary[s,3], summary[s,4], summary[s,5], summary[s,6], summary[s,7])
#	print(outstr, file=out)
#out.close()
#print("\n\n")

#### SYNTHETIC LIBRARY FITTING
print("\nComputing Synthetic Results.")

# Create synthetic library
synth = binocs2.makesynth(mag, binary, options)

# Run SED fitting on synthetic library
synth_results = binocs2.sedfit(singles, binary, synth, options)
	
# Compute Synthetic Results
synth_summary = binocs2.summarize(synth_results, binary, singles)
print("\nWriting synthetic results to '%s--synth.txt'" % (options['data']))
out = open(options['data']+"--synth.txt", "w")
for s in range(synth.shape[0]):
	# Print out star to file
	outstr = "%9.5f %9.5f " % (binary[s,0], binary[s,1])
	for i in range(17): outstr += "%6.3f " % synth[s,2*i]
	outstr += "   %6.3f %6.4f %6.3f %6.4f %5.2f   %6.3f %6.4f %5.2f" % (synth_summary[s,0], synth_summary[s,1], synth_summary[s,2], synth_summary[s,3], synth_summary[s,4], synth_summary[s,5], synth_summary[s,6], synth_summary[s,7])
	print(outstr, file=out)
out.close()

#### SYNTHETIC ANALYSIS
print("\nAnalyzing Synthetic Results...")

# % Completion
nfit = len(synth_summary[synth_summary[:,0] > 0,0])
print("    Pct Detected: %.1f" % (100*nfit/synth_summary.shape[0]))
nsin = len(binary[binary[:,1] == 0,0])
nsinfit = len([synth_summary[x,0] for x in range(synth_summary.shape[0]) if binary[x,1] == 0 and synth_summary[x,0] > 0])
print("    Pct Singles Detected: %.1f" % (100*nsinfit/nsin))

# Minimum Mass Ratio Determination
# Synthetic Fit Threshold
minq_synth = np.zeros(max(binary[:,0])//0.1)
for m in range(len(minq_synth)):
	binqs = [synth_summary[x,2] / synth_summary[x,0] for x in range(synth_summary.shape[0]) if synth_summary[x,0] > 0 and binary[x,0]//0.1 == m and binary[x,1] == 0]
	if len(binqs) == 0: continue
	minq_synth[m] = np.mean(binqs)

# Minimum Model Threshold
minq_mod = np.zeros(max(binary[:,0])//0.1)
minmass = min(binary[:,0])
for m in range(len(minq_mod)):
	minq_mod[m] = minmass/(m/10+0.05)
	
print("\nWriting minimum mass ratio results to '%s--minq.txt'" % (options['data']))
out = open(options['data']+"--minq.txt", "w")
for m in range(len(minq_synth)):
	print("%4.1f  %5.3f  %5.3f" % (m/10, minq_synth[m], minq_mod[m]))
out.close()



