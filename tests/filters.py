from __future__ import print_function, division
import numpy as np
import sys, binocs
from copy import copy

'''
	PROGRAM:		FILTERS
	DESCRIPTION: Computes 1 sigma uncertainty in mass determinations of synthetic stars.
	             Described in section 5.2 of Thompson et al. (in prep).
	INPUT: [From command line] BINOCS option file, specifying data file and isochrone to be used.
	OUTPUT: None
	FILE OUTPUT: "binocs_filter.%03d.dat" -- Output summary file for each filter combination used.
	                 Columns: [B Magnitude] [V Magnitude] [J Magnitude] [K Magnitude] [Actual Membership Flag] [Derived Membership Flag]
	                 Membership Flags: 0 = Non-Member; 1 = Member
'''

filter_combos =	[ np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]),	# 101: g[3.6]
				  np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0]),	# 111: gJ[3.6]
				  np.array([0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0]),	# 202: gr[3.6][4.5]
				  np.array([0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0]),	# 211: grJ[3.6]
				  np.array([0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0]),	# 222: grJK[3.6][4.5]
				  np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0]),	# 322: griJK[3.6][4.5]
				  np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0]),	# 332: griJHK[3.6][4.5]
				  np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0])]	# 532: ugrizJHK[3.6][4.5]

# Output files exist, read in this data and print to terminal
try:
	# Check to see whether files already exist
	out_files = subprocess.check_output("ls binocs_filter*", shell=True).splitlines()



# This has not be run before. Compute tests
except:
	options = binocs.readopt((sys.argv)[1])
	info, mag = binocs.readdata(options)
	oiso = binocs.readiso(options)
	singles = binocs.minterp(oiso, options['dm'])
	singles = binocs.fidiso(singles, options)
	binary = binocs.makebin(singles, options)
	synth = binocs.makesynth(mag, binary, options)
	
	# Check to see whether UBVRI or ugriz data exists
	if len(mag[mag[:,12] < 80,12]) == 0:
		# No ugriz data, swap columns with UBVRI
		for f in filter_combos:
			f[0:5] = f[5:10]
			f[5:10] = 0
			
	# Loop through filter combinations
	for f in filter_combos:
		# Print status
		print("\nWorking on %d%d%d..." % (np.sum(f[0:10]), np.sum(f[10:13]), np.sum(f[13:17])))
	
		# Remove filter magnitudes not being used for this run
		data = copy(synth)
		for i in range(17):
			# Remove this magnitude and error from the array
			if f[i] == 0: data[:,2*i], data[:,2*i+1] = 99.999, 9.999
				
		best_mass = np.zeros([synth.shape[0], 5, 2])
		for r in range(best_mass.shape[1]):
			results = binocs.sedfit(singles, binary, data, options, nvis=np.sum(f[0:10]), nnir=np.sum(f[10:13]), nmir=np.sum(f[13:17]))
			summary = binocs.summarize(results, binary, singles)
			
			# Compute mass uncertainties
			best_mass[:,r,0], best_mass[:,r,1] = summary[:,0], summary[:,2]
			
		# Print out results
		of = open("binocs_filter.%d%d%d.dat" % (np.sum(f[0:10]), np.sum(f[10:13]), np.sum(f[13:17])), 'w')
		for i in range(binary.shape[0]):
			outstr = "%7.4f  %s    %7.4f  %s" % (binary[i,0], ' '.join(["%7.4f" % x for x in best_mass[i,:,0]]), binary[i,1], ' '.join(["%7.4f" % x for x in best_mass[i,:,1]]))
			print(outstr, file=of)
		of.close()
		
		
		