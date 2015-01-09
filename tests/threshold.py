from __future__ import print_function, division
import numpy as np
import sys, binocs, subprocess

'''
	PROGRAM:		THRESHOLD
	DESCRIPTION: Computes confusion percentage (false positives and true negatives) for BINOCS routine, for various threshold values.
	             Described in section 5.3 of Thompson et al. (in prep).
	INPUT: [From command line] BINOCS option file, specifying data file and isochrone to be used.
	OUTPUT: None
	FILE OUTPUT: "binocs_thresh.%02d.dat" -- Output summary file for each threshold value used.
	                 Columns: [Primary Mass] [Secondary Mass] [B Magnitude] [V Magnitude] [J Magnitude] [K Magnitude] [Actual Membership Flag] [Derived Membership Flag]
	                 Membership Flags: 0 = Non-Member; 1 = Member
'''

thresholds = range(5,20)

# Output files exist, read in this data and print to terminal
try:
	# Check to see whether files already exist
	out_files = subprocess.check_output("ls binocs_thresh*", shell=True).splitlines()
	
	confusion_matrix = np.zeros([len(out_files), 3])
	for f in range(len(out_files)):
		this_threshold = float(out_files[f][14:16])
		confusion_matrix[f,0] = this_threshold
	
		data = np.loadtxt(out_files[f])
		data = data[data[:,0] > 0.5,:]
		true_members, true_non = data[data[:,6] == 1,:], data[data[:,6] == 0,:]
		confusion_matrix[f,1] = true_members[true_members[:,7] == 0,:].shape[0] / true_members.shape[0] * 100
		confusion_matrix[f,2] = true_non[true_non[:,7] == 1,:].shape[0] / true_non.shape[0] * 100
		print("%4.1f  %5.1f %5.1f %5.1f" % (confusion_matrix[f,0], confusion_matrix[f,1], confusion_matrix[f,2], confusion_matrix[f,1] + confusion_matrix[f,2]))
	
except Exception as e:
	options = binocs.readopt((sys.argv)[1])
	info, mag = binocs.readdata(options)
	oiso = binocs.readiso(options)
	singles = binocs.minterp(oiso, options['dm'])
	singles = binocs.fidiso(singles, options, file_output=False)
	binary = binocs.makebin(singles, options, file_output=False)
	single_synth = binocs.makesynth(mag, binary, options)
	
	# Only allow certain filters in final dataset (otherwise it's not showing reality)
	f = [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0]
	# Check to see whether UBVRI or ugriz data exists
	if len(mag[mag[:,12] < 80,12]) == 0:
		# No ugriz data, swap columns with UBVRI
		f[0:5] = f[5:10]
		f[5:10] = 0
	# Remove unnecessary magnitudes
	for i in range(17):
		if f[i] == 0: single_synth[:,2*i], single_synth[:,2*i+1] = 99.999, 9.999
	
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
	confusion_matrix = np.zeros([len(thresholds), 2])
	for t in range(len(thresholds)):
		print("\nWorking on threshold %.1f..." % (thresholds[t]))

		# Run SED fitting on entire synthetic set
		results = binocs.sedfit(singles, binary, synth, options, chicut=thresholds[t])
		summary = binocs.summarize(results, binary, singles)
	
		# Print results to file
		of = open("binocs_thresh.%02d.dat" % (thresholds[t]), 'w')
		for l in range(synth.shape[0]):
			bidx = l % single_synth.shape[0]
			if l < synthlen: member = 1
			else: member = 0
			if summary[l,0] > 0: result = 1
			else: result = 0
			print("%9.4f %9.4f  %8.3f %8.3f %8.3f %8.3f   %d %d" % (binary[bidx,0], binary[bidx,1], synth[l,2], synth[l,4], synth[l,20], synth[l,24], member, result), file=of)
		of.close()
