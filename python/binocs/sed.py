# BINOCS SED fitting subroutine
from __future__ import print_function, division
import numpy as np
import pyopencl as cl
from time import time
import sys, os

def sedfit(singles, binary, mag, options, chicut=10.0, nvis=3, nnir=3, nmir=2, chi=False):
	'''
	SUBROUTINE:			SEDFIT
	DESCRIPTION: Nearest-neighbor comparison between star data and synthetic models
	INPUT:       singles -- isochrone data from READISO, MINTERP or FIDISO
	             binary -- synthetic binary model data from MAKEBIN
	             mag -- star magnitude array from READDATA
	             options -- parameter dictionary from READOPT
	             chicut (optional, default = 10) -- minimum value a sum element must have to be added to total
	             nvis (optional, default = 2) -- minimum number of 'good' optical filters necessary to compute masses
	             nnir (optional, default = 2) -- minimum number of 'good' near-IR filters necessary to compute masses
	             nmir (optional, default = 2) -- minimum number of 'good' mid-IR filters necessary to compute masses
	OUTPUT:      4D matrix of mass determination information. Axes are:
	                 0: Star index. Aligns with mag
	                 1: 0 = fit chi value
	                    1 = best-fit binary model index. Aligns with binary
	                 2: Iteration index
	                 3: 0 = fitting results when compared to all binaries
	                    1 = fitting results when compared to only singles
	'''
	
	# Read in specified kernel
	pwd = os.path.dirname(os.path.realpath(__file__))
	if chi: df = open("%s/kernel/binocs_chi.c" % pwd, 'r')
	else: df = open("%s/kernel/binocs.c" % pwd, 'r')
	kernelstr = df.read().replace('GOODOPTICAL', '%d' % nvis).replace('GOODNEARIR', '%d' % nnir).replace('GOODMIDIR', '%d' % nmir)
	df.close()

	# Prepare OpenCL routine
	context = cl.create_some_context()
	queue = cl.CommandQueue(context)
	program = cl.Program(context, kernelstr).build()
	binsub = program.binsub
	binsub.set_scalar_arg_dtypes([None, None, None, None, None, np.float32, np.int32])
	
	# Copy model data to device
	d_single = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.ravel(singles).astype(np.float32))
	d_binary = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.ravel(binary).astype(np.float32))
	
	# Separate star magnitudes and uncertainties
	data, err = np.zeros([mag.shape[0], mag.shape[1]//2]), np.zeros([mag.shape[0], mag.shape[1]//2])
	for i in range(mag.shape[1]):
		# If this is a magnitude, convert it to absolute
		if i%2 == 0:
			data[:,i//2] = mag[:,i] 
			data[data[:,i//2] < 80, i//2] -= options['m-M'] + options['ebv'] * 3.08642 * options['ak'][i//2]
		# Save magnitude errors
		else: err[:,(i-1)//2] = mag[:,i]
	
	# Flatten data and err for use in OpenCL kernel
	data, err = np.ravel(data), np.ravel(err)
	err[data > 80] = 0
	
	# Copy star uncertainties to device
	d_err = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=err.astype(np.float32))
	
	# Choose ETA printing frequency, based on total number of runs
	if options['nruns'] < 200: p = 10
	elif options['nruns'] < 500: p = 30
	else: p = 50
	
	results = np.zeros([mag.shape[0], 2, options['nruns'], 2])
	
	# Begin loop over runs
	start_time = time()
	for r in range(options['nruns']):
		# Print progress
		if r < 1: print("    Run %3d of %3d " % (r, options['nruns']), end='')
		elif (r % p) == 0 and r > 0 and r < options['nruns'] - 4:
			time_perloop = (time() - start_time) / r
			time_left = ((options['nruns'] - r) * time_perloop)
			if time_left < 99: print(" ETA: %3d sec.\n    Run %3d of %3d " % (round(time_left), r, options['nruns']), end='')
			elif time_left < 5900: print(" ETA: %3.1f min.\n    Run %3d of %3d " % (time_left/60, r, options['nruns']), end='')
			else: print(" ETA: %3.1f hrs.\n    Run %3d of %3d " % ( time_left/360, r, options['nruns']), end='')
		sys.stdout.flush()
			
		# Randomize magnitudes
		rand1, rand2 = np.random.rand(len(data)), np.random.rand(len(data))
		rundata = data + np.sqrt(-2.0 * np.log(rand1)) * np.cos(2.0 * np.pi * rand2) * err / 2
		
		# Copy star magnitudes to device
		d_data = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=rundata.astype(np.float32))
		
		# Create output arrays
		bestchi, bestfit = np.zeros(len(rundata)//17).astype(np.float32), np.zeros(len(rundata)//17).astype(np.int32)
		d_chi, d_fit = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, bestchi.nbytes), cl.Buffer(context, cl.mem_flags.WRITE_ONLY, bestfit.nbytes)
		
		# Compare stars to binary models
		binsub(queue, bestchi.shape, None, d_binary, d_data, d_err, d_chi, d_fit, chicut, binary.shape[0])
		queue.finish()
		cl.enqueue_copy(queue, bestchi, d_chi)
		cl.enqueue_copy(queue, bestfit, d_fit)
		
		# Save results
		for s in range(len(bestchi)):
			if bestchi[s] > 0 and bestchi[s] < 999:
				results[s, 0, r, 0] = bestfit[s]
				results[s, 1, r, 0] = bestchi[s]
			else:
				results[s, 0, r, 0] = -1
				results[s, 1, r, 0] = -1.0
		
		# Compare stars to only single models
		binsub(queue, bestchi.shape, None, d_single, d_data, d_err, d_chi, d_fit, 1.0, singles.shape[0])
		queue.finish()
		cl.enqueue_copy(queue, bestchi, d_chi)
		cl.enqueue_copy(queue, bestfit, d_fit)
		
		# Save results
		for s in range(len(bestchi)):
			if bestchi[s] > 0 and bestchi[s] < 999:
				results[s, 0, r, 1] = bestfit[s]
				results[s, 1, r, 1] = bestchi[s]
			else:
				results[s, 0, r, 1] = -1
				results[s, 1, r, 1] = -1.0
		print('.', end='')
		sys.stdout.flush()
	
	# Print out completion message
	total_time = time() - start_time
	if total_time < 100: print("\n    %3d Runs Complete in %4.1f seconds." % (options['nruns'], total_time))
	elif total_time < 6000: print("\n    %3d Runs Complete in %4.1f minutes." % (options['nruns'], total_time/60))
	else: print("\n    %3d Runs Complete in %5.1f hours.\n" % (options['nruns'], total_time/3600))
	
	return results
	
	
	

def summarize(results, binary, singles):
	'''
	SUBROUTINE:			SUMMARIZE
	DESCRIPTION: Summarizes SED results into best-fit masses and uncertainties
	INPUT:       results -- full SED fitting results. Output from SEDFIT
	             binary -- synthetic binary model data from MAKEBIN
	             singles -- isochrone data from READISO, MINTERP or FIDISO
	OUTPUT:      summary -- Summarized BINOCS results
	                 0: best-fit primary mass
	                 1: uncertainty on primary mass
	                 2: best-fit secondary mass
	                 3: uncertainty on secondary mass
	                 4: average \Sigma per filter
	                 5: best-fit single star model mass
	                 6: uncertainty in single star model mass
	                 7: average \Sigma per filter for single star model matches
	                 8: binary flag
	                     0: unknown
	                     1: single
	                     2: binary
	                    -1: non-member
	'''
	summary = np.zeros([results.shape[0], 9])
	for s in range(results.shape[0]):
		starchi, staridx, singlechi, singleidx = results[s, 1, :, 0], results[s, 0, :, 0], results[s, 1, :, 1], results[s, 0, :, 1]
	
		# Find best-fit single star
		smchi = np.median(singlechi)
		if smchi > 0:
			mass = [singles[singleidx[l],0] for l in range(len(singlechi)) if singlechi[l] > 0]
			smass = np.median(mass)
			if len(mass) > 1: umass = np.std(mass)
			else: umass = 0.0
		else:
			smass = 0.0
			umass = 0.0
	
		# Find median chi value (this will determine whether the star is considered a member or not).
		medchi = np.median(starchi)
	
		# Star is not a cluster member
		if medchi < 0:
			bflag = -1
			mpri, upri, msec, usec, medchi = 0, 0, 0, 0, 0
	
		# Star is a cluster member
		else:
			# Find best-fit primary mass
			pri = [binary[staridx[l],0] for l in range(len(starchi)) if starchi[l] > 0]
			mpri = np.median(pri)
			if len(pri) > 1: upri = np.std(pri)
			else: upri = 0.0
	
			# Find best-fit secondary mass
			sec = [binary[staridx[l],1] for l in range(len(starchi)) if starchi[l] > 0]
			msec = np.median(sec)
			if len(sec) > 1: usec = np.std(sec)
			else: usec = 0.0
	
			# Determine binarity flag
			if msec / mpri > 0.3: bflag = 2
			else: bflag = 1
	
		summary[s,:] = [mpri, upri, msec, usec, medchi, smass, umass, smchi, bflag]
		
	return summary	
	