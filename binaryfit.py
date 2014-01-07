# Import modules
import pyopencl as cl
import sys, math
import matplotlib.pyplot as plt
import numpy as np
from time import time
from scipy import interpolate

# Isochrone comparison kernel
kernelstr = """
__kernel void binsub( __global float* iso, __global float* data, __global float* err, __global float* chi, const float chithresh ) {
	int m = get_global_id(0);
	chi[m] = -1.0;
	float tmpchi = 0.0, thischi = 0.0, totfilt = 0.0;
	int gubv = 0, gsds = 0, gvis = 0, gnir = 0, gmir = 0;
	
	// Loop through filters and compare star to the model
	for (int f = 0; f < 17; f++){
		thischi = (data[f] - iso[23*m+f+6])*(data[f] - iso[23*m+f+6]) / (err[f]*err[f]);
		if (thischi < chithresh){
			if (f < 5){ gubv++; }
			else if (f < 10){ gsds++; }
			else if (f < 13){ gnir++; }
			else{ gmir++; }
			tmpchi += thischi;
		}
	}
	// See which visual filter set has more matches
	if (gubv > gsds){ gvis = gubv; }
	else { gvis = gsds; }
	// See if this comparison has enough filters to be used
	if (gvis >= 3 && gnir >= 3 && gmir >= 2){
		totfilt = gvis + gnir + gmir;
		chi[m] = tmpchi / (totfilt - 2 - 1);
	}
}
"""





##### ISOCHRONE MASS INTERPOLATION SUBROUTINE ##############################################
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
		if bv[c] < minbv:
			minbv = bv[c]
			bindex = c
		
		# See if star is 0.5 mag redder than the bluest star
		if bv[c] - minbv >= 0.5:
			toindex = bindex
			
		# Emergency break (mass > 8 M_sun)
		if original[c][1] >= 8.0:
			toindex = c-4
			
		# Emergency break (mass degeneracy)
		if abs(original[c][1] - original[c+1][1]) < 0.001:
			toindex = c-4
			
	# End turn-off loop
	if toindex < 0:
		toindex = len(bv)-1
		
	print "\nTurn-Off Found:"
	print "    V = ", original[toindex][9]
	print "    B-V = ", bv[toindex]
	
	
	# Create mass grid
	cmass = float(int(original[2][1] * 100.0)) / 100.0 + 0.01
	print "    New Grid:",cmass,"-->",original[toindex][1]
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
	return npsingles
	
	
	
##### SYNTHETIC BINARY CREATION SUBROUTINE #################################################
def makebin(singles):
	binary = []
	# Loop through primary stars
	for p in range(len(singles)/23):
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
				mags[m] = -2.5 * math.log10( 10.0**(-1.0*singles[23*p+m+6]/2.5) + 10.0**(-1.0*singles[23*s+m+6]/2.5) )
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
	return npbinary
	
	
	
##### FIDUCIAL ISOCHRONE ADJUSTMENT SUBROUTINE #################################################
def fidiso(singles, fidname, fidmag, fidcol, ak, d, ebv):
	# Create new list to hold adjusted isochrone
	adj = []
	for r in range(len(singles)/23):
		for c in range(6): adj.append(singles[23*r+c])
		for c in range(6, 23): adj.append(99.999)
		
	# Read in fiducial file
	ff = open(fidname, "r")
	fidlines = ff.read().splitlines()
	fiddata = [f.split('\t') for f in fidlines]
	
	# Adjust fiducial data to absolute scale
	for l in range(len(fiddata)):
		# Adjust for distance and extinction
		fiddata[l][0] = float(fiddata[l][0]) - d - ak[fidmag[0]] / (ak[1] - ak[2]) * ebv
		# Adjust for reddening
		for c in range(len(fidcol)): fiddata[l][c+1] = float(fiddata[l][c+1]) - (ak[fidcol[c][0]] - ak[fidcol[c][1]]) / (ak[1] - ak[2]) * ebv
	
	# Loop through isochrone and save fiducial magnitude
	for r in range(len(singles)/23): adj[23*r+6+fidmag[0]] = singles[23*r+6+fidmag[0]]
	
	colcomplete = np.zeros(len(fidcol))
	# Loop through colors multiple times to adjust all necessary filters
	for l in range(3):
		# Loop through all colors specified in fiducial file
		for c in range(len(fidcol)):
			if colcomplete[c] == 1: continue
			# Check to see if one of the magnitudes necessary has already been solved for.
			goodmag = [i for i in range(len(adj)/23) if adj[23*i+6+fidcol[c][0]] < 80]
			goodcol = [i for i in range(len(adj)/23) if adj[23*i+6+fidcol[c][1]] < 80]
			
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
				for s in range(len(adj)/23):
					if adj[23*s+6+fidmag[0]] < min(datmag) or adj[23*s+6+fidmag[0]] > max(datmag): continue
					adj[23*s+6+fidcol[c][1]] = adj[23*s+6+fidcol[c][0]] - fit(adj[23*s+6+fidmag[0]])
				
			# Color filter solved for, but not magnitude
			if len(goodcol) > 0:
				for s in range(len(adj)/23):
					if adj[23*s+6+fidmag[0]] < min(datmag) or adj[23*s+6+fidmag[0]] > max(datmag): continue
					adj[23*s+6+fidcol[c][0]] = adj[23*s+6+fidcol[c][1]] + fit(adj[23*s+6+fidmag[0]])
					
	# Loop through filters and complete any missing entries
	for f in range(6, 23):
		# Find all values where this magnitude is already solved for
		goodmag = [i for i in range(len(adj)/23) if adj[23*i+f] < 80]
		
		# No fiducial for this filter
		if len(goodmag) == 0:
			for i in range(len(adj)/23): adj[23*i+f] = singles[23*i+f]
			
		# There was a fiducial for this filter
		else:
			# From the last index on, fill values
			lastidx = max(goodmag)
			for i in range(lastidx+1, len(adj)/23):
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
	return npadj
	
	
	
	
			
			
# Prepare OpenCL routine
context = cl.create_some_context()
queue = cl.CommandQueue(context)
program = cl.Program(context, kernelstr).build()
binsub = program.binsub
binsub.set_scalar_arg_dtypes([None, None, None, None, np.float32])
print "\n\n"
print "==========================================="
print "|                BINARYFIT                |"
print "|                  v2.0                   |"
print "|    GPU-Enabled (Requires OpenCL 1.2+)   |"
print "|                                         |"
print "|             by Ben Thompson             |"
print "==========================================="




##### MAIN PROGRAM #########################################################################
# Read in command line arguments
filter_names = ['U', 'B', 'V', 'R', 'I', 'SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4']
ak = [1.531, 1.324, 1.000, 0.748, 0.482, 1.593, 1.199, 0.858, 0.639, 0.459, 0.282, 0.175, 0.112, 0.0627, 0.0482, 0.0482, 0.0482]
clargs = sys.argv
optname = clargs[1]
fidcol = []

# Read in option file
of = open(optname, "r")
optlines = of.read().splitlines()
for l in optlines:
	if l.find('#') >= 0: continue
	tmp = [t.strip(' \t\n\r') for t in l.split("=")]
	if tmp[0] == "data": dataname = tmp[1]
	if tmp[0] == "iso":  isoname = tmp[1]
	if tmp[0] == "fid":  fidname = tmp[1]
	if tmp[0] == "fidmag":  fidmag = [f for f in range(len(filter_names)) if filter_names[f] == tmp[1]]
	if tmp[0] == "fidcol":
		fidcolname = [(t.strip(' \t\n\r')).split('-') for t in tmp[1].split(',')]
		for c,m in fidcolname:
			cidx = [f for f in range(len(filter_names)) if filter_names[f] == c]
			midx = [f for f in range(len(filter_names)) if filter_names[f] == m]
			fidcol.append([cidx[0], midx[0]])
	if tmp[0] == "dm":   dm = float(tmp[1])
	if tmp[0] == "age":  age = float(tmp[1])
	if tmp[0] == "m-M":  d = float(tmp[1])
	if tmp[0] == "ebv":  ebv = float(tmp[1])
	if tmp[0] == "nruns": nruns = int(tmp[1])
	if tmp[0] == "plot": plotstatus = int(tmp[1])


# Print out imported parameters
print "\nParameters:"
print "    data:",dataname
print "    iso:",isoname
print "    dm =",dm
print "    age =",age
print "    m-M =",d
print "    E(B-V) =",ebv

# Print run parameters to file
pf = open(dataname+"--notes.txt", "w")
print >>pf, "ISO: %s" % (isoname)
if len(fidcol) > 0: print >>pf, "FID: %s" % (fidname)
print >>pf, "AGE: %6.3f" % (age)
print >>pf, "DIS: %5.2f" % (d)
print >>pf, "EBV: %4.2f" % (ebv)
pf.close()


# Read in star data from file
datafile = open(dataname, 'r')
datalines = datafile.read().splitlines()
datafile.close()
data, id, rvflag = [], [], []
for dline in datalines:
	tmp = dline.split()
	new_tmp = [float(i) for i in tmp[1:47]]
	data.append(new_tmp)
	id.append(tmp[0])
	if tmp[48] == 'U': rvflag.append(0)
	elif tmp[48] == 'SM': rvflag.append(1)
	elif tmp[48] == 'BM' or tmp[48] == 'BLM': rvflag.append(2)
	else: rvflag.append(-1)
print "\nRead in %d stars from file." % (len(datalines))
	
	
# Read in isochrone from file
isofile = open(isoname, 'r')
isolines = isofile.read().splitlines()
isofile.close()
oiso = []
for iline in isolines:
	tmp = iline.split()
	new_tmp = [float(i) for i in tmp]
	if abs(new_tmp[0] - age) <= 0.001:
		oiso.append(new_tmp)
		
### Plot star data and original isochrone line
if plotstatus == 1:
	# Plot star data
	dmag, dcol = [], []
	for de in data:
		dmag.append(de[14])
		dcol.append(de[14] - de[16])
	imag, icol = [], []
	for ie in oiso:
		imag.append(ie[13] + d + ebv*ak[6]/0.324)
		icol.append(ie[13] - ie[14] + ebv*(ak[6]-ak[7])/0.324)
	plt.plot(dcol, dmag, "ko", icol, imag, "b-", markersize=1)
	plt.axis([-0.5, 2, 24, 12])
	plt.show()

# Interpolate isochrone to new mass grid
singles = minterp(oiso, dm)
print "\nIsochrone contains",len(singles)/23,"single stars."

# Plot newly interpolated isochrone
if plotstatus == 1:
	smag = [singles[23*m+6+6] + d + ebv*ak[6]/0.324 for m in range(len(singles)/23)]
	scol = [singles[23*m+6+6] - singles[23*m+7+6] + ebv*(ak[6]-ak[7])/0.324 for m in range(len(singles)/23)]
	plt.plot(dcol, dmag, "ko", icol, imag, "g-", scol, smag, "b-", markersize=1)
	plt.axis([-0.5, 2, 24, 12])
	plt.show()

# Adjust isochrone if necessary
if len(fidcol) > 0:
	fidsingles = fidiso(singles, fidname, fidmag, fidcol, ak, d, ebv)
	# Plot adjusted isochrones
	if plotstatus == 1:
		fmag = [fidsingles[23*m+6+6] + d + ebv*ak[6]/0.324 for m in range(len(fidsingles)/23)]
		fcol = [fidsingles[23*m+6+6] - fidsingles[23*m+7+6] + ebv*(ak[6]-ak[7])/0.324 for m in range(len(fidsingles)/23)]
		plt.plot(dcol, dmag, "ko", scol, smag, "g-", fcol, fmag, "b-", markersize=1)
		plt.axis([-0.5, 2, 24, 12])
		plt.show()
	singles = fidsingles
	print "\nAdjusted isochrone to fiducial sequence."
	ndirsplit = dataname.split('/')
	fidoutname = "%s/iso_%s.fid.dat" % ('/'.join(ndirsplit[0:len(ndirsplit)-1]), ndirsplit[len(ndirsplit)-2])
	fio = open(fidoutname, "w")
	for s in range(len(singles)/23):
		print >>fio, "%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (age, singles[23*s], singles[23*s+1], singles[23*s+2], singles[23*s+3], singles[23*s+4], singles[23*s+5], singles[23*s+6], singles[23*s+7], singles[23*s+8], singles[23*s+9], singles[23*s+10], singles[23*s+11], singles[23*s+12], singles[23*s+13], singles[23*s+14], singles[23*s+15], singles[23*s+16], singles[23*s+17], singles[23*s+18], singles[23*s+19], singles[23*s+20], singles[23*s+21], singles[23*s+22])
	fio.close()

# Create binary array
binary = makebin(singles)
print "\nCreated",len(binary)/23,"binary models for comparison."

# Print created binaries to file
ndirsplit = dataname.split('/')
dirsplit = isoname.split('/')
namesplit = dirsplit[len(dirsplit)-1].split('.')
if len(fidcol) == 0: isooutname = "%s/%s.m%03d.a%05d.bin" % ('/'.join(ndirsplit[0:len(ndirsplit)-1]), '.'.join(namesplit[0:len(namesplit)-1]), dm*100, age*1000)
else: isooutname = "%s/iso_%s.m%03d.a%05d.bin" % ('/'.join(ndirsplit[0:len(ndirsplit)-1]), ndirsplit[len(ndirsplit)-2], dm*100, age*1000)
isoo = open(isooutname, "w")
for b in range(len(binary)/23):
	print >>isoo, "%7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (binary[23*b], binary[23*b+1], binary[23*b+6], binary[23*b+7], binary[23*b+8], binary[23*b+9], binary[23*b+10], binary[23*b+11], binary[23*b+12], binary[23*b+13], binary[23*b+14], binary[23*b+15], binary[23*b+16], binary[23*b+17], binary[23*b+18], binary[23*b+19], binary[23*b+20], binary[23*b+21], binary[23*b+22])
isoo.close()
print "    Binary models output to '%s'\n" % (isooutname)

# Plot newly-created binaries
if plotstatus == 1:
	bmag = [binary[23*m+6+6] + d + ebv*ak[6]/0.324 for m in range(len(binary)/23)]
	bcol = [binary[23*m+6+6] - binary[23*m+7+6] + ebv*(ak[6]-ak[7])/0.324 for m in range(len(binary)/23)]
	plt.plot(dcol, dmag, "ko", bcol, bmag, "c+", markersize=3)
	plt.axis([-0.5, 2, 24, 12])
	plt.show()

# Copy necessary arrays to device for processing
d_single = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=singles.astype(np.float32))
d_binary = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=binary.astype(np.float32))

### Begin loop over fitting runs
start_time = time()
results = np.zeros((len(data), 2, nruns))
sresults = np.zeros((len(data), 2, nruns))
for r in range(nruns):
	if r < 1: print "Run %3d of %3d ..." % (r+1, nruns)
	elif (r % 10) == 9:
		time_perloop = (time() - start_time) / r
		time_left = ((nruns - r) * time_perloop)
		if time_left < 100:
			print "Run %3d of %3d ...  ETA: %2d sec." % (r+1, nruns, round(time_left))
		elif time_left < 6000:
			print "Run %3d of %3d ...  ETA: %2d min." % (r+1, nruns, round(time_left/60.0))
		else:
			print "Run %3d of %3d ...  ETA: %4.1f hrs." % (r+1, nruns, time_left/3600.0)

	
	# Randomize magnitudes
	rundata = np.zeros(len(data)*17)
	runerr = np.zeros(len(data)*17)
	for e in range(len(data)):
		rand1 = np.random.rand(17)
		rand2 = np.random.rand(17)
		for f in range(17):
			if data[e][2*f+2] > 80:
				rundata[17*e+f] = 99.999
				runerr[17*e+f] = 9.999
			else:
				rundata[17*e+f] = data[e][2*f+2] - d - ebv*3.08642*ak[f]
				rundata[17*e+f] += math.sqrt(-2.0 * math.log(rand1[f])) * math.cos(2.0 * math.pi * rand2[f]) * data[e][2*f+3]
				runerr[17*e+f] = data[e][2*f+3]
				
	# Plot running fit
	if plotstatus == 1:
		bmag = [binary[23*m+6+6] for m in range(len(binary)/23)]
		bcol = [binary[23*m+6+6] - binary[23*m+7+6] for m in range(len(binary)/23)]
		rmag = [rundata[17*m+6] for m in range(len(rundata)/17)]
		rcol = [rundata[17*m+6] - rundata[17*m+7] for m in range(len(rundata)/17)]
		plt.plot(rcol, rmag, "ko", bcol, bmag, "c+", markersize=3)
		plt.axis([-0.5, 2, 12, -2])
		plt.show()
				
	### Begin loop over stars
	for s in range(len(data)):
		# Copy out this star from rundata
		thisdata = rundata[17*s:17*s+17]
		d_data = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=thisdata.astype(np.float32))
		thiserr = runerr[17*s:17*s+17]
		d_err = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=thiserr.astype(np.float32))
		
		### COMPARE STARS TO BINARY MODELS
		# Create output array
		thischi = np.zeros(len(binary)/23).astype(np.float32)
		d_chi = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, thischi.nbytes)
		
		# Run kernel
		binsub(queue, thischi.shape, None, d_binary, d_data, d_err, d_chi, 7.0)
		queue.finish()
		cl.enqueue_copy(queue, thischi, d_chi)
		
		# Save results
		sidx = thischi.argsort()
		if thischi[sidx[len(sidx)-1]] >= 0:
			bidx = np.searchsorted(thischi[sidx], [-0.01,], side='right')[0]
			results[s, 0, r] = sidx[bidx]
			results[s, 1, r] = thischi[sidx[bidx]]
		else:
			results[s, 0, r] = -1
			results[s, 1, r] = -1.0
			
		
		### COMPARE STARS TO SINGLE MODELS
		# Create output array
		thischi = np.zeros(len(singles)/23).astype(np.float32)
		d_chi = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, thischi.nbytes)
		
		# Run kernel
		binsub(queue, thischi.shape, None, d_single, d_data, d_err, d_chi, 14.0)
		queue.finish()
		cl.enqueue_copy(queue, thischi, d_chi)
		
		# Save results
		sidx = thischi.argsort()
		if thischi[sidx[len(sidx)-1]] >= 0:
			bidx = np.searchsorted(thischi[sidx], [-0.01,], side='right')[0]
			sresults[s, 0, r] = sidx[bidx]
			sresults[s, 1, r] = thischi[sidx[bidx]]
		else:
			sresults[s, 0, r] = -1
			sresults[s, 1, r] = -1.0
			
			
# Print out completion message
total_time = time() - start_time
if total_time < 100: print "\n%3d Runs Complete in %4.1f seconds." % (nruns, total_time)
elif total_time < 6000: print "\n%3d Runs Complete in %4.1f minutes." % (nruns, total_time/60.0)
else: print "\n%3d Runs Complete in %5.1f hours.\n" % (nruns, total_time/3600.0)

# Now that we have completed all runs, we compute final results
# Open file for writing
print "\nWriting results to '%s--binary.txt'." % (dataname)
out = open(dataname+"--binary.txt", "w")

# Create grid for mass distribution results
#	mass -> 0 to 5 in 0.1 increments
#	q -> 0 to 1 in 0.1 increments
gmass = [y*0.1 for y in range(51)]
gq = [y*0.1 for y in range(11)]
grid = np.zeros((51, 11))

for s in range(len(data)):
	starchi = results[s, 1, :]
	staridx = results[s, 0, :]
	singlechi = sresults[s, 1, :]
	singleidx = sresults[s, 0, :]
	
	# Find best-fit single star
	smchi = np.median(singlechi)
	if smchi > 0:
		mass = [singles[23*singleidx[l]] for l in range(nruns) if singlechi[l] > 0]
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
		pri = [binary[23*staridx[l]] for l in range(nruns) if starchi[l] > 0]
		mpri = np.median(pri)
		if len(pri) > 1: upri = np.std(pri)
		else: upri = 0.0
		
		# Find best-fit secondary mass
		sec = [binary[23*staridx[l]+1] for l in range(nruns) if starchi[l] > 0]
		msec = np.median(sec)
		if len(sec) > 1: usec = np.std(sec)
		else: usec = 0.0
		
		# Determine binarity flag
		if msec / mpri > 0.4: bflag = 2
		else: bflag = 1
		
	# Print out star to file
	print >>out, "%16s %9.5f %9.5f   %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f   %2d %2d %6.3f %6.4f %6.3f %6.4f %5.2f   %6.3f %6.4f %5.2f" % (id[s], data[s][0], data[s][1], data[s][2], data[s][4], data[s][6], data[s][8], data[s][10], data[s][12], data[s][14], data[s][16], data[s][18], data[s][20], data[s][22], data[s][24], data[s][26], data[s][28], data[s][30], data[s][32], data[s][34], bflag, rvflag[s], mpri, upri, msec, usec, medchi, smass, umass, smchi)
	
	# Add star to mass grid
	if medchi > 0:
		gmidx = round(mpri / 0.1)
		if gmidx > 51: continue
		gqidx = round((msec / mpri) / 0.1)
		grid[gmidx, gqidx] += 1
		
out.close()

# Plot final mass distribution
gplotx = []
gploty = []
gsize = []
df = open(dataname+"--qdist.txt", "w")
print "\nWriting q distribution data to '%s--qdist.txt'." % (dataname)
for x in range(len(gmass)):
	totstars = sum(grid[x,:])
	if totstars == 0: continue
	for y in range(len(gq)):
		if grid[x,y] == 0: continue
		gplotx.append(gmass[x])
		gploty.append(gq[y])
		gsize.append((float(grid[x,y])/totstars)*500.0)
		print >>df, "%4.1f %3.1f %6.4f %6d" % (gmass[x], gq[y], float(grid[x,y])/float(totstars), totstars)
plt.scatter(gplotx, gploty, s=gsize, c=gsize)
plt.axis([0, 2.2, -0.02, 1.02])
plt.xlabel("Stellar Mass (M_o)")
plt.ylabel("Mass Ratio   q")
#plt.show()
df.close()


print "\n\n"

















