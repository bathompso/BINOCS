from __future__ import print_function, division
import numpy as np
import sys, binocs, subprocess
from copy import copy

'''
	PROGRAM:		FILTERS
	DESCRIPTION: Computes 1 sigma uncertainty in mass determinations of synthetic stars.
	             Described in section 5.2 of Thompson et al. (in prep).
	INPUT: [From command line] BINOCS option file, specifying data file and isochrone to be used.
	OUTPUT: None
	FILE OUTPUT: "binocs_filter.%03d.dat" -- Output summary file for each filter combination used.
	                 Columns: [Actual Primary Mass] [5x Primary Mass Estimates] [Actual Secondary Mass] [5x Secondary Mass Estimates]
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
	out_files = subprocess.check_output("ls binocs_filter*.dat", shell=True).splitlines()
	pct_error = np.zeros([len(out_files), 11, 2])
	
	# Loop through files and compute percent uncertainties
	for f in range(len(out_files)):
		data = np.loadtxt(out_files[f])
		data_errors = np.zeros([data.shape[0], 4])
		data_errors[:,0], data_errors[:,2] = data[:,0], data[:,6]
	
		for l in range(data.shape[0]):
			# Primary percent uncertainty
			pri_err = [np.abs(data[l,i]-data[l,0]) / data[l,0] * 100 for i in range(1,6) if data[l,i] > 0]
			if len(pri_err) > 0: data_errors[l,1] = np.mean(pri_err) + np.std(pri_err)
			else: data_errors[l,1] = -1
			# Secondary percent uncertainty
			if data[l,6] == 0: continue
			sec_err = [np.abs(data[l,i]-data[l,6]) / data[l,6] * 100 for i in range(7,12) if data[l,i-6] > 0]
			if len(sec_err) > 0: data_errors[l,3] = np.mean(sec_err) + np.std(sec_err)
			else: data_errors[l,3] = -1
			
		# Collapse percent errors into grid
		for q in range(11):
			bin_pri = [data_errors[x,1] for x in range(data_errors.shape[0]) if data_errors[x,2] / data_errors[x,0] >= q/10 and data_errors[x,2] / data_errors[x,0] < (q+1)/10 and data_errors[x,1] >= 0]
			if len(bin_pri) > 0: pct_error[f,q,0] = np.mean(bin_pri)
			
			bin_sec = [data_errors[x,3] for x in range(data_errors.shape[0]) if data_errors[x,2] / data_errors[x,0] >= q/10 and data_errors[x,2] / data_errors[x,0] < (q+1)/10 and data_errors[x,3] >= 0]
			if len(bin_sec) > 0: pct_error[f,q,1] = np.mean(bin_sec)
			
	# Print results to screen
	for m in range(2):
		print("\n")
		for f in range(len(out_files)):
			print("%3s  %s" % (out_files[f][14:17], ' '.join(["%5.1f" % x for x in pct_error[f,:,m]])))
			
	# Print results to LaTeX file
	filternames = ['U','B','V','R','I','u','g','r','i','z','J','H','K_S','[3.6]','[4.5]','[5.8]','[8.0]']
	of = open('binocs_filter.tex', 'w')
	# Print Header
	print('\\begin{table*} \centering \small', file=of)
	print('\\begin{tabular}{|l|ccccccccccc|c|}', file=of)
	print('\multicolumn{1}{c}{} & \multicolumn{11}{c}{Mass Ratio} & \multicolumn{1}{c}{} \\\\', file=of)
	print('\multicolumn{1}{c}{Filters} & 0.0 & 0.1 & 0.2 & 0.3 & 0.4 & 0.5 & 0.6 & 0.7 & 0.8 & 0.9 & \multicolumn{1}{c}{1.0} & \n\t\multicolumn{1}{c}{\multirow{8}{*}{\\vspace{-0.7cm}\\begin{turn}{-90}1$\sigma$ \% Error in $M_{\\text{pri}}$ \end{turn}}} \\\\ \hline \hline', file=of)
	for m in range(2):
		for f in range(len(out_files)):
			filtdisplay = np.array(copy(filternames))
			filtdisplay[filter_combos[f] == 0] = '.'
			if m == 0: outstr = "%3s: $%s$%s & %s" % (out_files[f][14:17], ''.join(filtdisplay[5:13]), ''.join(filtdisplay[13:]), ' & '.join(["%5.1f" % x for x in pct_error[f,:,m]]))
			else: outstr = "%3s: $%s$%s & ... & %s" % (out_files[f][14:17], ''.join(filtdisplay[5:13]), ''.join(filtdisplay[13:]), ' & '.join(["%5.1f" % x for x in pct_error[f,1:,m]]))
			# First line of table has special ending
			if f == len(out_files)-1: 
				if m == 0: outstr += ' & \n\t\multirow{8}{*}{\\vspace{-0.7cm}\\begin{turn}{-90}1$\sigma$ \% Error in $M_{\\text{sec}}$ \end{turn}} \\\\ \hline \hline'
				else: outstr += ' & \\\\ \hline'
			else: outstr += ' & \\\\'
			print(outstr, file=of)
	# Print footer
	print('\end{tabular}', file=of)
	print('\caption{1$\sigma$ \% errors in mass estimates for various combinations of filters. \label{tab:filters_test}}', file=of)
	print('\end{table*}', file=of)
	of.close()



# This has not be run before. Compute tests
except:
	options = binocs.readopt((sys.argv)[1])
	info, mag = binocs.readdata(options)
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
		
		
		