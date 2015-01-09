from __future__ import print_function, division
import numpy as np
import sys, binocs, subprocess
from copy import copy

'''
	PROGRAM:		FILTERS
	DESCRIPTION: Computes 1 sigma uncertainty in mass determinations of synthetic stars as a function of number of required filters.
	             Described in section 5.2 of Thompson et al. (in prep).
	INPUT: [From command line] BINOCS option file, specifying data file and isochrone to be used.
	OUTPUT: None
	FILE OUTPUT: "binocs_filter_synth.%03d.dat" -- Output synthetic test summary file for each filter combination used.
	                 Columns: [Actual Primary Mass] [5x Primary Mass Estimates] [Actual Secondary Mass] [5x Secondary Mass Estimates]
'''

filternames = ['U','B','V','R','I','u','g','r','i','z','J','H','K_S','[3.6]','[4.5]','[5.8]','[8.0]']
filter_combos =	[ np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]),	# 101: g[3.6]
				  np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0]),	# 111: gJ[3.6]
				  np.array([0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0]),	# 202: gr[3.6][4.5]
				  np.array([0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0]),	# 211: grJ[3.6]
				  np.array([0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0]),	# 222: grJK[3.6][4.5]
				  np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0]),	# 322: griJK[3.6][4.5]
				  np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0]),	# 332: griJHK[3.6][4.5]
				  np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0])]	# 532: ugrizJHK[3.6][4.5]

# Read in photometry file (necessary for both modes)
options = binocs.readopt((sys.argv)[1])
info, mag = binocs.readdata(options)

# Determine nvis, nnir, nmir for this photometry file
phot = np.zeros([mag.shape[0], 3])
for i in range(mag.shape[0]):
	phot[i,0] = len(max([[x for x in range(0,10,2) if mag[i,x] < 80], [x for x in range(10,20,2) if mag[i,x] < 80]]))
	phot[i,1] = len([x for x in range(20,26,2) if mag[i,x] < 80])
	phot[i,2] = len([x for x in range(26,34,2) if mag[i,x] < 80])

# Output files exist, read in this data and print to terminal
try:
	print("\n\n")
	# Check to see whether files already exist
	synth_files = subprocess.check_output("ls binocs_filter*.dat", shell=True).splitlines()
	
	for f in range(len(synth_files)):
		# Determine number of stars possible to use in this dataset
		minvis, minnir, minmir = int(synth_files[f][14]), int(synth_files[f][15]), int(synth_files[f][16])
		npossible = len(phot[(phot[:,0] >= minvis) & (phot[:,1] >= minnir) & (phot[:,2] >= minmir)])

		# Read in data from results files for this filter combination
		synth_data = np.loadtxt(synth_files[f])

		# Determine synthetic mass uncertainties
		synth_errors = np.zeros([synth_data.shape[0]*5, 2])
		for l in range(synth_data.shape[0]):
			for i in range(5):
				# Primary Mass % Uncertainty
				if synth_data[l,i+1] > 0:
					synth_errors[l*5+i,0] = np.abs(synth_data[l,i+1]-synth_data[l,0]) / synth_data[l,0] * 100
				else: synth_errors[l*5+i,0] = -1

				# Seconday Mass % Uncertainty
				if synth_data[l,i+1] > 0 and synth_data[l,6] > 0:
					synth_errors[l*5+i,1] = np.abs(synth_data[l,i+7]-synth_data[l,6]) / synth_data[l,6] * 100
				else: synth_errors[l*5+i,1] = -1
		pri_err, sec_err = synth_errors[synth_errors[:,0] >= 0,0], synth_errors[synth_errors[:,1] >= 0,1]

		# Print summary for this file
		filtdisplay = np.array(copy(filternames))
		filtdisplay[filter_combos[f] == 0] = '.'
		print("%3s: $%10s$%-12s & %5d & %5.1f & %5.1f \\\\" % (synth_files[f][14:17], ''.join(filtdisplay[5:13]), ''.join(filtdisplay[13:]), npossible, np.mean(pri_err), np.mean(sec_err)))

except Exception as e:
	print(e)

	oiso = binocs.readiso(options)
	singles = binocs.minterp(oiso, options['dm'])
	singles = binocs.fidiso(singles, options, file_output=False)
	binary = binocs.makebin(singles, options, file_output=False)
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
			if f[i] == 0: data[:,2*i], data[:,2*i+1] = 99.999, 9.999
				
		best_mass = np.zeros([synth.shape[0], 5, 2])
		for r in range(best_mass.shape[1]):
			results = binocs.sedfit(singles, binary, data, options, nvis=np.sum(f[0:10]), nnir=np.sum(f[10:13]), nmir=np.sum(f[13:17]))
			summary = binocs.summarize(results, binary, singles)
			best_mass[:,r,0], best_mass[:,r,1] = summary[:,0], summary[:,2]
			
		# Print out synthetic results
		of = open("binocs_filter.%d%d%d.dat" % (np.sum(f[0:10]), np.sum(f[10:13]), np.sum(f[13:17])), 'w')
		for i in range(binary.shape[0]):
			outstr = "%7.4f  %s    %7.4f  %s" % (binary[i,0], ' '.join(["%7.4f" % x for x in best_mass[i,:,0]]), binary[i,1], ' '.join(["%7.4f" % x for x in best_mass[i,:,1]]))
			print(outstr, file=of)
		of.close()






