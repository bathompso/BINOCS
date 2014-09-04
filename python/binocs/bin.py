# BINOCS synthetic dataset creation routines
from __future__ import print_function, division
import numpy as np
import os, sys

def makebin(iso, options, file_output=True):
	'''
	SUBROUTINE:			MAKEBIN
	DESCRIPTION: Flux-combines single isochrone stars into model binaries
	INPUT:       iso -- isochrone data
	             options -- parameters dictionary from READOPT
	             file_output -- boolean flag to determine whether file with model magnitudes should be output
	OUTPUT:      bin -- Binary star data
	                  0-1: Primary / Secondary Mass
	                  2-5: Zeros
	                  6-23: Magnitudes
	FILE OUTPUT: iso.m[dm].a[age].bin -- File containing data on binaries generated for this run.
	                  0-1: Primary / Secondary Mass
	                  2-18: UBVRIugrizJHK[3][4][5][8] Magnitudes
	'''
	
	# Create initial matrix to hold all binary models
	bmod = np.zeros([iso.shape[0]**2+iso.shape[0]+1, 23])
	
	# Loop through primary stars
	print("\nCreating synthetic binary models... ", end='')
	sys.stdout.flush()
	for p in range(iso.shape[0]):
		# Add single star to resulting array
		bmod[iso.shape[0]*p,0] = iso[p,0]
		bmod[iso.shape[0]*p,6:] = iso[p,6:]
		
		# Loop through secondary stars
		lastq = -1
		for s in range(p+1):
			# Check to make sure we're generating a different-enough mass ratio
			if iso[s,0] / iso[p,0] - 0.02 < lastq: continue
			lastq = iso[s,0] / iso[p,0]
			
			# Calculate new magnitudes
			newmags = -2.5 * np.log10( 10.0**(-1.0*iso[p,6:]/2.5) + 10.0**(-1.0*iso[s,6:]/2.5) )
			
			# Skip adding this binary if it is not different enough from other binaries
			magdiff = len(newmags[iso[p,6:] - newmags < 0.001])
			if magdiff > 3: continue
			
			# Add this binary to the dataset
			bmod[iso.shape[0]*p+s+p+1,0] = iso[p,0]
			bmod[iso.shape[0]*p+s+p+1,1] = iso[s,0]
			bmod[iso.shape[0]*p+s+p+1,6:] = newmags
			
	# Copy out only models that have magnitudes
	bin = bmod[bmod[:,0] > 0,:]
	print(" Done.")
	print("    Created %d binary models for comparison." % bin.shape[0])
	
	# Print created binaries to file
	if file_output:
		if len(options['data'].split('/')) == 1: outdir = ''
		else: outdir = '/'.join(options['data'].split('/')[0:len(options['data'].split('/'))-1]) + '/'
	
		if 'fid' not in options.keys(): basename = 'iso'
		else: basename = options['fid'].split('/')[-1].split('.')[0]
	
		binoutname = "%s%s.m%03d.a%05d.bin" % (outdir, basename, options['dm']*100, options['age']*1000)
		bo = open(binoutname, 'w')
		for b in range(bin.shape[0]):
			outstr = "%7.4f %7.4f " % (bin[b,0], bin[b,1])
			for i in range(6,23): outstr += "%6.3f " % bin[b,i]
			print(outstr, file=bo)
		bo.close()
		print("    Binary models output to '%s'" % (binoutname))
	
	return bin




def makesynth(mag, binary, options):
	'''
	SUBROUTINE:			MAKESYNTH
	DESCRIPTION: Generates synthetic star dataset for testing routines.
	INPUT:       mag -- matrix of observed UBVRIugrizJHK[3][4][5][8] magnitudes + uncertainties
	             binary -- matrix of binary star data from MAKEBIN
	             options -- parameters dictionary from READOPT
	OUTPUT:      synth -- matrix of UBVRIugrizJHK[3][4][5][8] magnitudes + uncertainties. Copy of mag structure from READDATA for the synthetic dataset.
	'''
	# Find range of magnitudes for V or g
	err_mag = 12
	if len(mag[mag[:,err_mag] < 80, err_mag]) == 0: err_mag = 4
	good_mag = mag[mag[:,err_mag] < 80, err_mag]
	mag_bins = np.arange(int(max((good_mag) - min(good_mag)) / 0.5)) * 0.5 + min(good_mag)
	
	# Find average errors for each filter
	avg_err = np.zeros([int((max(good_mag) - min(good_mag)) / 0.5), 17])
	filt_used = np.zeros(17)
	for f in range(17):
		if len(mag[mag[:,2*f] < 80, 2*f]) == 0: continue
		filt_used[f] = 1
		
		# Compute average uncertainty
		for m in range(len(mag_bins)):
			binerrs = [mag[x,2*f+1] for x in range(mag.shape[0]) if mag[x,err_mag] >= mag_bins[m] and mag[x,err_mag] < mag_bins[m]+0.5 and mag[x,2*f] < 80]
			if len(binerrs) == 0: avg_err[m,f] = -1
			else: avg_err[m,f] = np.mean(binerrs)
		
		# Fill in missing slots (at top)
		for m in range(1,len(mag_bins)):
			if avg_err[m,f] < 0 and avg_err[m-1,f] > 0: avg_err[m,f] = avg_err[m-1,f]
		
		# Fill in missing slots (at bottom)
		for m in range(len(mag_bins)-2,-1,-1):
			if avg_err[m,f] < 0 and avg_err[m+1,f] > 0: avg_err[m,f] = avg_err[m+1,f]
			
	print("    Filters to be used: %s" % (' '.join([options['filternames'][x] for x in range(len(filt_used)) if filt_used[x] == 1]))) 
	
	# Create new input array with all synthetic stars
	print("    Creating Synthetic Grid... ", end='')
	sys.stdout.flush()
	synth = np.zeros([binary.shape[0],34])
	for f in range(17): 
		if filt_used[f] == 1: synth[:,2*f] = binary[:,f+6] + options['m-M'] + options['ebv'] * 3.08642 * options['ak'][f]
		else: synth[:,2*f] = 99.999
	
	# Adjust errors for correct bin
	for m in range(len(mag_bins)):
		bin_idx = np.array([synth[x,err_mag] >= mag_bins[m] and synth[x,err_mag] < mag_bins[m]+0.5 for x in range(synth.shape[0])])
		for f in range(17): 
			if filt_used[f] == 1: synth[bin_idx,2*f+1] = avg_err[m,f]
			else: synth[bin_idx,2*f+1] = 9.999
	
	# Give errors to stars outside range
	for f in range(17): synth[synth[:,2*f+1] == 0,2*f+1] = avg_err[-1,f]
	
	# Multiply uncertainties by 2
	for f in range(17): synth[:,2*f+1] *= 2
	
	# Randomize magnitudes
	for f in range(17):
		if filt_used[f] == 0: continue
		rand1, rand2 = np.random.rand(synth.shape[0]), np.random.rand(synth.shape[0])
		synth[:,2*f] = synth[:,2*f] + np.sqrt(-2 * np.log(rand1)) * np.cos(2 * np.pi * rand2) * synth[:,2*f+1]
	print("Done.")
	
	return synth
	