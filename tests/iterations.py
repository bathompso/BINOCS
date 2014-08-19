from __future__ import print_function, division
import numpy as np
import sys, binocs, subprocess

'''
	PROGRAM:		ITERATIONS
	DESCRIPTION: Computes per-run variation in derived stellar masses, for various numbers of iterations. 
	             Described in section 5.1 of Thompson et al. (in prep)
	INPUT: [From command line] BINOCS option file, specifying data file and isochrone to be used.
	OUTPUT: None
	FILE OUTPUT: "binocs_iter.%04d.dat" -- Output summary file for each number of iterations used.
	                 Columns: [Star Line Number] [5x Individual Run Primary Mass] [Primary Sigma] [5x Individual Run Secondary Mass] [Secondary Sigma]
	                 Sigma values are explained in section 5.1 of Thompson et al. (in prep)
'''
nruns = [3, 10, 30, 90, 150, 200, 300, 400, 500, 600, 700, 1200]

# Output files exist, read in this data and print to terminal
try:
	# Check to see whether files already exist
	out_files = subprocess.check_output("ls binocs_iter*", shell=True).splitlines()

	for f in out_files:
		this_nruns = int(f.split('.')[1])
	
		df = open(f, 'r')
		lines = df.read().splitlines()
		df.close()
		primary_sigma, secondary_sigma = np.zeros(len(lines)), np.zeros(len(lines))
		
		# Compute new sigmas for every line
		for l in range(len(lines)):
			tmp = np.array([float(x) for x in lines[l].split()])
			primary_summary = tmp[1:6]
			secondary_summary = tmp[7:12]
			
			primary_valid = primary_summary[primary_summary > 0]
			secondary_valid = secondary_summary[primary_summary > 0]
			
			if len(primary_valid) == 0:
				primary_sigma[l], secondary_sigma[l] = -1, -1
				continue
			
			primary_sigma[l] = np.std(primary_valid) / np.mean(primary_valid)
			if np.mean(secondary_valid) > 0: secondary_sigma[l] = np.std(secondary_valid) / np.mean(secondary_valid)
			
		# Fix overly bad secondary Sigmas
		secondary_sigma[secondary_sigma > 1] = 1
			
		# Compute final statistics
		primary_sigma, secondary_sigma = np.sort(primary_sigma), np.sort(secondary_sigma)
		pmed, smed = np.median(primary_sigma), np.median(secondary_sigma)
		p95, s95 = primary_sigma[int(0.97*len(primary_sigma))], secondary_sigma[int(0.97*len(secondary_sigma))]
		
		print("%4d  %.3f %.3f %.3f %.3f" % (this_nruns, pmed, p95, smed, s95))
			
		

# This has not be run before. Compute tests
except:
	options = binocs.readopt((sys.argv)[1])
	info, mag = binocs.readdata(options)
	oiso = binocs.readiso(options)
	singles = binocs.minterp(oiso, options['dm'])
	singles = binocs.fidiso(singles, options, file_output=False)
	binary = binocs.makebin(singles, options, file_output=False)

	# Loop through different number of iterations
	sigmas = np.zeros([len(nruns), 4])
	for n in range(len(nruns)):
		print("\nWorking on %dx..." % nruns[n]) 
		options['nruns'] = nruns[n]
		iteration_summary = np.zeros([mag.shape[0], 2, 5])
	
		# Open output file
		of = open('binocs_iter.%04d.dat' % nruns[n], 'w')
	
		# Complete all runs for the iteration amount
		for r in range(iteration_summary.shape[2]):
			results = binocs.sedfit(singles, binary, mag, options)
			summary = binocs.summarize(results, binary, singles)
			iteration_summary[:,0,r] = summary[:,0]
			iteration_summary[:,1,r] = summary[:,2]
			
		# Compute Sigma
		primary_sigma, secondary_sigma = np.zeros(mag.shape[0]), np.zeros(mag.shape[0])
		for m in range(mag.shape[0]):
			valid = [x for x in range(iteration_summary.shape[2]) if iteration_summary[m,0,x] > 0]
			if len(valid) == 0:
				primary_sigma[m], secondary_sigma[m] = -1, -1
				continue
			primary_sigma[m] = np.std(iteration_summary[m,0,:]) / np.mean([iteration_summary[m,0,x] for x in valid])
			if np.mean([iteration_summary[m,1,x] for x in valid]) > 0:
				secondary_sigma[m] = np.std(iteration_summary[m,1,:]) / np.mean([iteration_summary[m,1,x] for x in valid])
			
			# Print out summary statistics
			pri_out, sec_out = '', ''
			for i in range(iteration_summary.shape[2]):
				pri_out += "%6.3f " % iteration_summary[m,0,i]
				sec_out += "%6.3f " % iteration_summary[m,1,i]
			print("%7d  %s %7.4f   %s %7.4f" % (m, pri_out, primary_sigma[m], sec_out, secondary_sigma[m]), file=of)
		of.close()