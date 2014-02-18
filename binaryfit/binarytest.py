# Import modules
import numpy as np
import sys, os, binocs

print "\n\n"
print "==========================================="
print "|                BINARYFIT                |"
print "|              MASS ACCURACY              |"
print "|             by Ben Thompson             |"
print "==========================================="

# Read in options
options = binocs.readopt((sys.argv)[1])

# Check to see whether mass uncertainty file has already been created
if not os.path.exists(options['data']+"--mass.txt"):

	### MODE 1: GENERATE MASS UNCERTAINTY FILE
	data = binocs.readdata(options)
	oiso = binocs.readiso(options)
	
	# Interpolate isochrone to new mass grid
	singles = binocs.minterp(oiso, options['dm'])

	# Adjust isochrone to empirical ridgeline, if necessary
	singles = binocs.fidiso(singles, options)

	# Create binary array
	binary = binocs.makebin(singles, options)

	# Find average errors for magnitude bins in V or g
	good, gidx = [x for x in range(len(data)) if data[x][14] < 80], 14
	if len(good) == 0: good, gidx = [x for x in range(len(data)) if data[x][6] < 80], 6
	gmag, gstep = [data[g][gidx] for g in good], 0.3
	gbins = [min(gmag) + x * gstep for x in range(int((max(gmag) - min(gmag))/gstep)+1)]

	# Loop through bins and find average error for all filters
	errs = np.zeros([len(gbins), 17])
	for i in range(len(gbins)):
		binidx = [x for x in range(len(data)) if data[x][gidx] >= gbins[i] and data[x][gidx] < gbins[i] + gstep]
		for f in range(17): 
			gerrs = [data[x][2*f+3] for x in binidx if data[x][2*f+3] < 8]
			if len(gerrs) > 0 and len(binidx) > 0: errs[i,f] = np.mean(gerrs)
			else: errs[i,f] = 9.999
	# Loop through errors again, and fill in empty values with nearest neighbor
	for f in range(17):
		bad = [x for x in range(len(gbins)) if errs[x,f] > 8]
		good = [x for x in range(len(gbins)) if errs[x,f] < 8]
		if len(bad) == len(gbins): continue
		for i in range(len(bad)):
			nei = np.argsort([abs(gbins[x] - gbins[bad[i]]) for x in good])[0]
			errs[bad[i],f] = errs[good[nei],f]	
		
	# Loop through data and determine what filters are available in input data
	printf = []
	for f in range(17):
		good = [x for x in range(len(data)) if data[x][f] < 80]
		if len(good) == 0: printf.append(0)
		else: printf.append(1)
	
	# OPTIONAL: Loop through data and trim possible filters to specified number
	# UBVRI OR ugriz
	if gidx == 14: 
		for f in range(5): printf[f] = 0
	else: 
		for f in range(5,10): printf[f] = 0

	# Loop through isochrone binary stars and create a new data array
	bindata, ak = [], options['ak']
	for l in range(len(binary)/23):
		binline = [0.00, -1.00]
		# Find which error bin this most closely matches
		bindiff = [abs(binary[23*l+6+(gidx/2-1)] - gbins[x]) for x in range(len(gbins))]
		bin = (np.argsort(bindiff))[0] 
		# Add magnitudes and errors to array
		rand1 = np.random.rand(17)
		rand2 = np.random.rand(17)
		for f in range(17):
			if printf[f] == 1:
				binline.append(binary[23*l+6+f] + options['m-M'] + options['ebv'] * 3.08642 * ak[f] + np.sqrt(-2.0 * np.log(rand1[f])) * np.cos(2.0 * np.pi * rand2[f]) * errs[bin,f] * 2.0)
				binline.append(errs[bin,f] * 2.0)
			else:
				binline.append(99.999)
				binline.append(9.999)
		bindata.append(binline)
	
	# Fit these single isochrone stars in the SED routine
	massresults = binocs.sedfit(singles, binary, bindata, options)

	# Print out results of isochrone single star fitting
	of = open(options['data']+"--mass.txt", 'w')
	fitmass = np.zeros([len(bindata), 2])
	for i in range(len(bindata)):
		opri = binary[23*i]
		osec = binary[23*i+1]
	
		# Fit to all models
		bfitchi = massresults[i,1,:,0]
		bfitidx = massresults[i,0,:,0]
		if np.median(bfitchi) < 0: continue
		bfitmpri = np.median([binary[23*bfitidx[x]] for x in range(len(bfitchi)) if bfitchi[x] > 0])
		bfitmsec = np.median([binary[23*bfitidx[x]+1] for x in range(len(bfitchi)) if bfitchi[x] > 0])
		bruns = len([x for x in range(len(bfitchi)) if bfitchi[x] > 0])
	
		fitmass[i,0], fitmass[i,1] = bfitmpri, bfitmsec
	
		if osec > 0: secpcterr = (osec-bfitmsec)/osec*100
		elif bfitmsec > 0: secpcterr = 100
		else: secpcterr = 0
	
		print >>of, "%6.3f %6.3f  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f  %6.3f %6.3f %6.1f %6.1f  %3d" % (opri, osec, bindata[i][2], bindata[i][4], bindata[i][6], bindata[i][8], bindata[i][10], bindata[i][12], bindata[i][14], bindata[i][16], bindata[i][18], bindata[i][20], bindata[i][22], bindata[i][24], bindata[i][26], bindata[i][28], bindata[i][30], bindata[i][32], bindata[i][34], bfitmpri, bfitmsec, (opri-bfitmpri)/opri*100, secpcterr, bruns)
	of.close()	


	# MINIMUM MASS RATIO
	# Find all single models
	sidx = [x for x in range(len(bindata)) if binary[23*x+1] == 0]
	smass = [binary[23*x] for x in sidx]
	# Compute bins of single stars in steps of 0.1 M_sun
	ds = 0.05
	sbins = [min(smass) + ds * x for x in range(int((max(smass) - min(smass)) / ds) + 1)]
	minq = np.zeros(len(sbins))
	# Loop through bins and compute maximum returned mass ratio
	print " MASS   MIN Q"
	for s in range(len(sbins)):
		minfitq = max([fitmass[x,1] / fitmass[x,0] for x in range(len(bindata)) if binary[23*x+1] == 0 and fitmass[x,0] > 0 and binary[23*x] >= sbins[s] and binary[23*x] < sbins[s] + ds])
		minmodq = min([binary[23*x+1] / binary[23*x] for x in range(len(bindata)) if binary[23*x] >= sbins[s] and binary[23*x] < sbins[s] + ds and binary[23*x+1] > 0])
		minq[s] = max([minfitq, minmodq])
		print "%5.2f   %5.3f" % (sbins[s], minq[s])

else:
	#### MODE 2: ANALYZE BINARY FILE AND PRINT RESULTS
	# Read in mass uncertainty file
	df = open(options['data']+'--mass.txt', 'r')
	lines = df.read().splitlines()
	df.close()
	massinfo = np.zeros([len(lines), 4])
	for l in range(len(lines)):
		tmp = lines[l].split()
		massinfo[l,0], massinfo[l,1], massinfo[l,2], massinfo[l,3] = float(tmp[0]), float(tmp[1]), float(tmp[19]), float(tmp[20])
	
	# Determine minimum mass ratios
	print "\n\nMinimum Mass Ratios:"
	ds = 0.05
	sbins = [min(massinfo[:,0]) + ds * x for x in range(int((max(massinfo[:,0]) - min(massinfo[:,0])) / ds) + 1)]
	minq = np.zeros(len(sbins))
	print " MASS   MIN Q"
	for s in range(len(sbins)):
		minfitq = max([massinfo[x,3] / massinfo[x,2] for x in range(len(lines)) if massinfo[x,1] == 0 and massinfo[x,0] >= sbins[s] and massinfo[x,0] < sbins[s] + ds])
		minmodq = min([massinfo[x,1] / massinfo[x,0] for x in range(len(lines)) if massinfo[x,1] > 0 and massinfo[x,0] >= sbins[s] and massinfo[x,0] < sbins[s] + ds])
		minq[s] = max([minfitq, minmodq])
		print "%5.2f   %5.3f" % (sbins[s], minq[s])
		
	# Read in binary results file
	df = open(options['data']+'--binary.txt', 'r')
	lines = df.read().splitlines()
	df.close()
	
	# Loop through binary data and adjust binarity + mass info
	of = open(options['data']+'--results.txt', 'w')
	for l in range(len(lines)):
		tmp = lines[l].split()
		
		# No mass info if this is an unknown or non-member star
		if int(tmp[20]) < 1:
			binflag = int(tmp[20])
			mpri, epri, msec, esec, chi = 0.0, 0.0, 0.0, 0.0, 0.0
		
		# New determination of binary flag
		else:
			opri, osec = float(tmp[22]), float(tmp[24])
			thissbin = (np.argsort([abs(opri - sbins[x]) for x in range(len(sbins))]))[0]
			# Star is below minimum q threshold: declared single
			if osec / opri <= minq[thissbin]:
				binflag = 1
				mpri, epri, msec, esec, chi = float(tmp[27]), float(tmp[28]), 0.0, 0.0, float(tmp[29])
			# Star is above threshold: declared binary
			else:
				binflag = 2
				mpri, epri, msec, esec, chi = float(tmp[22]), float(tmp[23]), float(tmp[24]), float(tmp[25]), float(tmp[26])
				
		# Print out new results
		outstr = "%16s %9.5f %9.5f   " % (tmp[0], float(tmp[1]), float(tmp[2]))
		for i in range(3,20): outstr += "%6.3f " % (float(tmp[i]))
		outstr += "  %2d %2d  %6.3f %6.4f %6.3f %6.4f %5.2f" % (binflag, int(tmp[21]), mpri, epri, msec, esec, chi)
		print >>of, outstr
		
	of.close()









