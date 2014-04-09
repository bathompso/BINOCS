##### BINOCS SUBROUTINE PYTHON MODULE #######################################################
# This module holds the main routines used in the BINOCS project.
# c. 08/2013: Ben Thompson <bathompso@gmail.com>


# Import Necessary Modules
from __future__ import print_function, division
import numpy as np
import pyopencl as cl
from time import time
from scipy import interpolate


##### OPTION FILE READIN SUBROUTINE ########################################################
# Description:	Reads in option file to dictionary
# Input:		optname -- name of input option file
# Output:		Dictionary containing all important parameters for binaryfit routine
def readopt(optname):
	filter_names = ['U', 'B', 'V', 'R', 'I', 'SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4']
	ak = [1.531, 1.324, 1.000, 0.748, 0.482, 1.593, 1.199, 0.858, 0.639, 0.459, 0.282, 0.175, 0.112, 0.0627, 0.0482, 0.0482, 0.0482]
	
	# Get path to current option file
	optdir = '/'.join((optname.split('/'))[0:-1]) + '/'
	
	# Declare empty dictionary
	options = dict()
	options['ak'] = ak
	options['filternames'] = filter_names

	# Read in options from file
	of = open(optname, "r")
	optlines = of.read().splitlines()
	for l in optlines:
		if l.find('#') >= 0: continue
		tmp = [t.strip(' \t\n\r') for t in l.split("=")]
		if tmp[0] == "data": options['data'] = optdir+tmp[1]
		if tmp[0] == "iso":  options['iso'] = optdir+tmp[1]
		if tmp[0] == "fid":  options['fid'] = optdir+tmp[1]
		if tmp[0] == "dm":   options['dm'] = float(tmp[1])
		if tmp[0] == "age":  options['age'] = float(tmp[1])
		if tmp[0] == "m-M":  options['m-M'] = float(tmp[1])
		if tmp[0] == "ebv":  options['ebv'] = float(tmp[1])
		if tmp[0] == "nruns": options['nruns'] = int(tmp[1])
		if tmp[0] == "dr": options['dr'] = float(tmp[1])
		
	# Print out imported parameters
	print("\nParameters:")
	print("    data:", options['data'])
	print("    iso:", options['iso'])
	print("    dm =", options['dm'])
	print("    age =", options['age'])
	print("    m-M =", options['m-M'])
	print("    E(B-V) =", options['ebv'])
	
	return options



##### DATAFILE READIN SUBROUTINE ###########################################################
# Description:	Reads in star data from a magnitude file created by PAYST
# Input:		options -- parameter dictionary from readopt
# Output:		List of star data. Elements in each row are:
#				0-1: RA / Dec Coordinates
#				2-35: UBVRIugrizJHK[3.6][4.5][5.8][8.0] Magnitudes + Uncertainties
#				36: 2MASS ID, or a similarly generated ID from PAYST
#				37: RV Variability index -- (2) = Binary, (1) = Single, (0) = Unknown, 
#											(-1) = Non-Member
def readdata(options):
	# Read in star data from file
	datafile = open(options['data'], 'r')
	datalines = datafile.read().splitlines()
	datafile.close()
	data = []
	for dline in datalines:
		tmp = dline.split()
		new_tmp = [float(i) for i in tmp[1:37]]
		new_tmp.append(tmp[0])
		if tmp[48] == 'U': new_tmp.append(0)
		elif tmp[48] == 'SM': new_tmp.append(1)
		elif tmp[48] == 'BM' or tmp[48] == 'BLM': new_tmp.append(2)
		else: new_tmp.append(-1)
		data.append(new_tmp)
	print("\nRead in %d stars from file." % (len(data)))
	return data
	
	
	
##### ISOCHRONE READIN SUBROUTINE ##########################################################
# Description:	Read in isochrone data from a file created by makeiso
# Input:		options -- parameter dictionary from readopt
# Output:		List of isochrone data. Values are the same as in each line of makeiso file.
def readiso(options):
	# Read in isochrone from file
	isofile = open(options['iso'], 'r')
	isolines = isofile.read().splitlines()
	isofile.close()
	oiso = []
	for iline in isolines:
		tmp = iline.split()
		new_tmp = [float(i) for i in tmp]
		if abs(new_tmp[0] - options['age']) <= 0.001:
			oiso.append(new_tmp)
	return oiso



##### ISOCHRONE MASS INTERPOLATION SUBROUTINE ##############################################
# Description:	Interpolates original isochrone into a more fine mass grid
# Input:		original -- original isochrone data from readiso
#				dm -- Mass increment between resulting isochrone points
# Output:		1D NumPy array of isochrone data for stars in new finely-gridded isochrone.
#					Each star has 23 entries in array, same as values in makeiso line.
def minterp(original, dm):
	# Find turnoff in B-V
	toindex = -1
	minbv = 99 
	bindex = -1
	bv = [grid[8] - grid[9] for grid in original]
	for c in range(len(bv)-1):
		if toindex >= 0:
			continue
			
		# See if star is bluest yet
		if bv[c] < minbv and bv[c] < 4:
			minbv = bv[c]
			bindex = c
		
		# See if star is 0.5 mag redder than the bluest star
		if bv[c] - minbv >= 0.5:
			toindex = bindex
			
		# Emergency break (mass > 8 M_sun)
		if original[c][1] >= 8.0:
			toindex = c-4
			
		# Emergency break (mass degeneracy)
		if abs(original[c][1] - original[c+1][1]) < 0.0005:
			toindex = c-4
			
	# End turn-off loop
	if toindex < 0:
		toindex = len(bv)-1
		
	print("\nTurn-Off Found:")
	print("    V = ", original[toindex][9])
	print("    B-V = ", bv[toindex])
	
	
	# Create mass grid
	cmass = float(int(original[2][1] * 100.0)) / 100.0 + 0.01
	print("    New Grid:",cmass,"-->",original[toindex][1])
	mgrid = [grid[1] for grid in original[0:toindex+4]]
	newmass = []
	while cmass < original[toindex][1]:
		newmass.append(cmass)
		# Increment the mass counter for next pass
		if cmass < 1.0: cmass += dm / 2.0
		elif cmass < 1.5: cmass += dm
		else: cmass += 2.0 * dm
		
	# Interpolate parameters and magnitudes to new grid
	interpolated = np.zeros((len(newmass)+len(original)-toindex, len(original[0])-1))
	for i in range(len(newmass)):
		interpolated[i,0] = newmass[i]
	for i in range(2,len(original[0])-1):
		thisparam = [grid[i+1] for grid in original[0:toindex+4]]
		fit = interpolate.interp1d(mgrid, thisparam, kind=3)
		interpy = fit(newmass)
		for j in range(len(newmass)):
			interpolated[j,i] = interpy[j]
			
	# Add original isochrone points beyond the turnoff to the isochrone
	for i in range(len(original)-toindex):
		interpolated[i+len(newmass),0] = original[i+toindex][1]
		for j in range(3,len(original[0])):
			interpolated[i+len(newmass),j-1] = original[i+toindex][j]

	singles = []
	for i in range(len(interpolated)):
		for j in range(len(interpolated[i])):
			singles.append(interpolated[i,j])

	npsingles = np.zeros(len(singles))
	for i in range(len(singles)): npsingles[i] = singles[i]
	
	print("\nIsochrone contains",len(npsingles)/23,"single stars.")
	
	return npsingles



##### FIDUCIAL ISOCHRONE ADJUSTMENT SUBROUTINE ##############################################
# Description:	Adjusts original isochrone data to empirical ridgelines
# Input:		singles -- 1D Numpy array from minterp
#				options -- parameter dictionary from readopt
# Output:		1D NumPy array of isochrone data for stars with empirically-adjusted
#					magnitudes.
def fidiso(singles, options):
	# Read in parameters from dictionary
	ak = options['ak']
	
	# Check to see if this operation is necessary
	if options['fid'] == '': return singles

	# Create new list to hold adjusted isochrone
	adj = []
	for r in range(len(singles)//23):
		for c in range(6): adj.append(singles[23*r+c])
		for c in range(6, 23): adj.append(99.999)
		
	# Read in fiducial file
	ff = open(options['fid'], "r")
	fidlines = ff.read().splitlines()
	fiddata = [fidlines[x].replace('\t',' ').split() for x in range(1,len(fidlines))]
	
	# Read in fiducial magnitude and colors
	tmp = fidlines[0].replace('\t',' ').split()
	fidmag = [x for x in range(len(options['filternames'])) if (options['filternames'])[x] == tmp[0]]
	fidcol = []
	for f in range(1,len(tmp)):
		tmpcol = tmp[f].split('-')
		c = [x for x in range(len(options['filternames'])) if (options['filternames'])[x] == tmpcol[0]]
		m = [x for x in range(len(options['filternames'])) if (options['filternames'])[x] == tmpcol[1]]
		fidcol.append([c[0], m[0]])
	
	# Adjust fiducial data to absolute scale
	for l in range(len(fiddata)-1):
		# Adjust for distance and extinction
		fiddata[l][0] = float(fiddata[l][0]) - options['m-M'] - ak[fidmag[0]] / (ak[1] - ak[2]) * options['ebv']
		# Adjust for reddening
		for c in range(len(fidcol)): fiddata[l][c+1] = float(fiddata[l][c+1]) - (ak[fidcol[c][0]] - ak[fidcol[c][1]]) / (ak[1] - ak[2]) * options['ebv']
	
	# Loop through isochrone and save fiducial magnitude
	for r in range(len(singles)//23): adj[23*r+6+fidmag[0]] = singles[23*r+6+fidmag[0]]
	
	colcomplete = np.zeros(len(fidcol))
	# Loop through colors multiple times to adjust all necessary filters
	for l in range(3):
		# Loop through all colors specified in fiducial file
		for c in range(len(fidcol)):
			if colcomplete[c] == 1: continue
			# Check to see if one of the magnitudes necessary has already been solved for.
			goodmag = [i for i in range(len(adj)//23) if adj[23*i+6+fidcol[c][0]] < 80]
			goodcol = [i for i in range(len(adj)//23) if adj[23*i+6+fidcol[c][1]] < 80]
			
			# Neither magnitude has data, skip it
			if len(goodmag) == 0 and len(goodcol) == 0: continue
			
			# Both magnitudes have been completely solved for
			if len(goodmag) > 0 and len(goodcol) > 0:
				colcomplete[c] = 1
				continue
				
			# Compute interpolation for colors
			datmag = [float(f[0]) for f in fiddata if float(f[c+1]) > -1000]
			datcol = [float(f[c+1]) for f in fiddata if float(f[c+1]) > -1000]
			fit = interpolate.interp1d(datmag, datcol, kind=3)
			
			# Magnitude filter solved for, but not color
			if len(goodmag) > 0:
				for s in range(len(adj)//23):
					if adj[23*s+6+fidmag[0]] < min(datmag) or adj[23*s+6+fidmag[0]] > max(datmag): continue
					adj[23*s+6+fidcol[c][1]] = adj[23*s+6+fidcol[c][0]] - fit(adj[23*s+6+fidmag[0]])
				
			# Color filter solved for, but not magnitude
			if len(goodcol) > 0:
				for s in range(len(adj)//23):
					if adj[23*s+6+fidmag[0]] < min(datmag) or adj[23*s+6+fidmag[0]] > max(datmag): continue
					adj[23*s+6+fidcol[c][0]] = adj[23*s+6+fidcol[c][1]] + fit(adj[23*s+6+fidmag[0]])
					
	# Loop through filters and complete any missing entries
	for f in range(6, 23):
		# Find all values where this magnitude is already solved for
		goodmag = [i for i in range(len(adj)//23) if adj[23*i+f] < 80]
		
		# No fiducial for this filter
		if len(goodmag) == 0:
			for i in range(len(adj)//23): adj[23*i+f] = singles[23*i+f]
			
		# There was a fiducial for this filter
		else:
			# From the last index on, fill values
			lastidx = max(goodmag)
			for i in range(lastidx+1, len(adj)//23):
				prevdiff = singles[23*i+f] - singles[23*(i-1)+f]
				adj[23*i+f] = prevdiff + adj[23*(i-1)+f]
			# From the first index and below, fill values
			firstidx = min(goodmag)
			for i in range(firstidx-1, -1, -1):
				prevdiff = singles[23*i+f] - singles[23*(i+1)+f]
				adj[23*i+f] = prevdiff + adj[23*(i+1)+f]
	
	# Return new isochrone
	npadj = np.zeros(len(adj))
	for i in range(len(adj)): npadj[i] = adj[i]
	
	ndirsplit = options['data'].split('/')
	fidoutname = "%s/iso_%s.fid.dat" % ('/'.join(ndirsplit[0:len(ndirsplit)-1]), ndirsplit[len(ndirsplit)-2])
	fio = open(fidoutname, "w")
	for s in range(len(singles)//23):
		print("%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (options['age'], npadj[23*s], npadj[23*s+1], npadj[23*s+2], npadj[23*s+3], npadj[23*s+4], npadj[23*s+5], npadj[23*s+6], npadj[23*s+7], npadj[23*s+8], npadj[23*s+9], npadj[23*s+10], npadj[23*s+11], npadj[23*s+12], npadj[23*s+13], npadj[23*s+14], npadj[23*s+15], npadj[23*s+16], npadj[23*s+17], npadj[23*s+18], npadj[23*s+19], npadj[23*s+20], npadj[23*s+21], npadj[23*s+22]), file=fio)
	fio.close()
	print("\nAdjusted isochrone to fiducial sequence.")
	
	return npadj
	
	
	
##### SYNTHETIC BINARY CREATION SUBROUTINE #################################################
# Description:	Flux-combines single isochrone stars into model binaries
# Input:		singles -- list of isochrone data from readiso or fidiso
#				options -- parameter dictionary from readopt
# Output:		List of binary star data. Elements in each row are:
#					0: Primary Mass
#					1: Secondary Mass
#					2-18: UBVRIugrizJHK[3.6][4.5][5.8][8.0] magnitudes
def makebin(singles, options):
	binary = []
	# Loop through primary stars
	for p in range(len(singles)//23):
		# Add the single star to the array
		for f in range(23): binary.append(singles[23*p+f])
	
		lastq = -1
		# Loop through secondary stars
		for s in range(p+1):
			# Check to make sure that we are generating a different mass ratio
			if (singles[23*s] / singles[23*p] - 0.02 < lastq): continue
			
			# Calculate new magnitudes
			mags = np.zeros(17)
			diff = 0
			for m in range(len(mags)):
				mags[m] = -2.5 * np.log10( 10.0**(-1.0*singles[23*p+m+6]/2.5) + 10.0**(-1.0*singles[23*s+m+6]/2.5) )
				if (abs(mags[m] - singles[23*p+m+6]) < 0.001): diff += 1
			# Skip adding this binary if it is not different enough from other binaries	
			if (diff > 3): continue
			
			# Add this binary to dataset
			lastq = singles[23*s] / singles[23*p]
			binary.append(singles[23*p])
			binary.append(singles[23*s])
			for j in range(4): binary.append(0)
			for m in range(len(mags)):
				binary.append(mags[m])
	
	npbinary = np.zeros(len(binary))
	for i in range(len(binary)): npbinary[i] = binary[i]
	
	print("\nCreated",len(npbinary)//23,"binary models for comparison.")
	# Print created binaries to file
	ndirsplit = options['data'].split('/')
	dirsplit = options['iso'].split('/')
	namesplit = dirsplit[len(dirsplit)-1].split('.')
	if options['fid'] == '': isooutname = "%s/%s.m%03d.a%05d.bin" % ('/'.join(ndirsplit[0:len(ndirsplit)-1]), '.'.join(namesplit[0:len(namesplit)-1]), options['dm']*100, options['age']*1000)
	else: isooutname = "%s/iso_%s.m%03d.a%05d.bin" % ('/'.join(ndirsplit[0:len(ndirsplit)-1]), ndirsplit[len(ndirsplit)-2], options['dm']*100, options['age']*1000)
	isoo = open(isooutname, "w")
	for b in range(len(binary)//23):
		print("%7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (binary[23*b], binary[23*b+1], binary[23*b+6], binary[23*b+7], binary[23*b+8], binary[23*b+9], binary[23*b+10], binary[23*b+11], binary[23*b+12], binary[23*b+13], binary[23*b+14], binary[23*b+15], binary[23*b+16], binary[23*b+17], binary[23*b+18], binary[23*b+19], binary[23*b+20], binary[23*b+21], binary[23*b+22]), file=isoo)
	isoo.close()
	print("    Binary models output to '%s'\n" % (isooutname))
	
	return npbinary
	
	
	
##### SED BINARY FITTING MAIN LOOP ##########################################################
# Description:	Compares star data to synthetic binary models to determine masses.
# Input:		singles -- list of isochrone data from readiso or fidiso
#				binary -- synthetic binary models from makebin
#				data -- star data array from readdata
#				options -- parameter dictionary from readopt
# Output:		4-D array of mass determination information. 4 indexes are:
#					0: Star index, aligns with element of data array
#					1: 0 = fit chi value
#					   1 = best-fit binary model index, aligns with element of binary
#					2: Iteration index
#					3: 0 = fitting to all binary models
#					   1 = fitting to only single star models
def sedfit(singles, binary, data, options):
	# Read data out of option dictionary
	ak = options['ak']
	nruns = options['nruns']

	# Isochrone comparison kernel
	kernelstr = """
	__kernel void binsub( __global float* iso, __global float* data, __global float* chi, const float chithresh ) {
		int m = get_global_id(0);
		chi[m] = -1.0;
		float tmpchi = 0.0, thischi = 0.0, totfilt = 0.0;
		int gubv = 0, gsds = 0, gvis = 0, gnir = 0, gmir = 0;
	
		// Loop through filters and compare star to the model
		for (int f = 0; f < 17; f++){
			thischi = 1.0 / (fabs(data[f] - iso[23*m+f+6]) + 0.01);
			if (thischi > chithresh){
				if (f < 5){ gubv++; }
				else if (f < 10){ gsds++; }
				else if (f < 13){ gnir++; }
				else{ gmir++; }
				totfilt++;
				tmpchi += thischi;
			} 
			// If star is more than 2 magnitudes away on this filter it *will not* fit the star. Abort.
			else if (thischi < 0.5 && data[f] < 80) { break; }
		}
		// See which visual filter set has more matches
		if (gubv > gsds){ gvis = gubv; }
		else {gvis = gsds; }
		// See if this comparison has enough filters to be used
		if (gvis >= %d && gnir >= %d && gmir >= %d){
			tmpchi /= totfilt;
			chi[m] = tmpchi;
		}
	}
	""" % (3,3,2)
	
	# Prepare OpenCL routine
	context = cl.create_some_context()
	queue = cl.CommandQueue(context)
	program = cl.Program(context, kernelstr).build()
	binsub = program.binsub
	binsub.set_scalar_arg_dtypes([None, None, None, np.float32])
	
	# Copy necessary arrays to device for processing
	d_single = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=singles.astype(np.float32))
	d_binary = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=binary.astype(np.float32))

	# Choose ETA printing frequency, based on total number of runs
	if nruns < 200: p = 10
	elif nruns < 500: p = 30
	elif nruns < 1000: p = 50
	else: p = 100
	
	### Begin loop over fitting runs
	start_time = time()
	results = np.zeros((len(data), 2, nruns, 2))
	for r in range(nruns):
		if r < 1: print("Run %3d of %3d ..." % (r+1, nruns))
		elif (r % p) == (p-1):
			time_perloop = (time() - start_time) / r
			time_left = ((nruns - r) * time_perloop)
			if time_left < 99:
				print("Run %3d of %3d ...  ETA: %2d sec." % (r+1, nruns, round(time_left)))
			elif time_left < 5900:
				print("Run %3d of %3d ...  ETA: %2d min." % (r+1, nruns, round(time_left/60.0)))
			else:
				print("Run %3d of %3d ...  ETA: %4.1f hrs." % (r+1, nruns, time_left/3600.0))

	
		# Randomize magnitudes
		rundata = np.zeros(len(data)*17)
		for e in range(len(data)):
			rand1 = np.random.rand(17)
			rand2 = np.random.rand(17)
			for f in range(17):
				if data[e][2*f+2] > 80:
					rundata[17*e+f] = 99.999
				else:
					rundata[17*e+f] = data[e][2*f+2] - options['m-M'] - options['ebv'] * 3.08642 * ak[f]
					rundata[17*e+f] += np.sqrt(-2.0 * np.log(rand1[f])) * np.cos(2.0 * np.pi * rand2[f]) * data[e][2*f+3]
				
		### Begin loop over stars
		for s in range(len(data)):
			# Copy out this star from rundata
			thisdata = rundata[17*s:17*s+17]
			d_data = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=thisdata.astype(np.float32))
		
			### COMPARE STARS TO BINARY MODELS
			# Create output array
			thischi = np.zeros(len(binary)//23).astype(np.float32)
			d_chi = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, thischi.nbytes)
		
			# Run kernel
			binsub(queue, thischi.shape, None, d_binary, d_data, d_chi, 10.0)
			queue.finish()
			cl.enqueue_copy(queue, thischi, d_chi)
		
			# Save results
			sindex = thischi.argsort()
			if thischi[sindex[len(sindex)-1]] > 0:
				results[s, 0, r, 0] = sindex[len(sindex)-1]
				results[s, 1, r, 0] = thischi[sindex[len(sindex)-1]]
			else:
				results[s, 0, r, 0] = -1
				results[s, 1, r, 0] = -1.0
			
		
			### COMPARE STARS TO SINGLE MODELS
			# Create output array
			thischi = np.zeros(len(singles)//23).astype(np.float32)
			d_chi = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, thischi.nbytes)
		
			# Run kernel
			binsub(queue, thischi.shape, None, d_single, d_data, d_chi, 1.0)
			queue.finish()
			cl.enqueue_copy(queue, thischi, d_chi)
		
			# Save results
			sindex = thischi.argsort()
			if thischi[sindex[len(sindex)-1]] > 0:
				results[s, 0, r, 1] = sindex[len(sindex)-1]
				results[s, 1, r, 1] = thischi[sindex[len(sindex)-1]]
			else:
				results[s, 0, r, 1] = -1
				results[s, 1, r, 1] = -1.0
	
	# Print out completion message
	total_time = time() - start_time
	if total_time < 100: print("\n%3d Runs Complete in %4.1f seconds." % (nruns, total_time))
	elif total_time < 6000: print("\n%3d Runs Complete in %4.1f minutes." % (nruns, total_time/60.0))
	else: print("\n%3d Runs Complete in %5.1f hours.\n" % (nruns, total_time/3600.0))
	
	return results
	
	
	

