# Import modules
import numpy as np
import sys, binocs
	
	
print "\n\n"
print "==========================================="
print "|                BINARYFIT                |"
print "|             by Ben Thompson             |"
print "==========================================="


# Read in data from files
options = binocs.readopt((sys.argv)[1])
data = binocs.readdata(options)
oiso = binocs.readiso(options)

# Interpolate isochrone to new mass grid
singles = binocs.minterp(oiso, options['dm'])

# Adjust isochrone to empirical ridgeline, if necessary
singles = binocs.fidiso(singles, options)

# Create binary array
binary = binocs.makebin(singles, options)

# Run SED fitting on all stars in dataset
results = binocs.sedfit(singles, binary, data, options)

# Compute Initial Results
# Open file for writing
print "\nWriting results to '%s--binary.txt'." % (options['data'])
out = open(options['data']+"--binary.txt", "w")
for s in range(len(data)):
	starchi = results[s, 1, :, 0]
	staridx = results[s, 0, :, 0]
	singlechi = results[s, 1, :, 1]
	singleidx = results[s, 0, :, 1]

	# Find best-fit single star
	smchi = np.median(singlechi)
	if smchi > 0:
		mass = [singles[23*singleidx[l]] for l in range(len(singlechi)) if singlechi[l] > 0]
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
		pri = [binary[23*staridx[l]] for l in range(len(starchi)) if starchi[l] > 0]
		mpri = np.median(pri)
		if len(pri) > 1: upri = np.std(pri)
		else: upri = 0.0

		# Find best-fit secondary mass
		sec = [binary[23*staridx[l]+1] for l in range(len(starchi)) if starchi[l] > 0]
		msec = np.median(sec)
		if len(sec) > 1: usec = np.std(sec)
		else: usec = 0.0

		# Determine binarity flag
		if msec / mpri > 0.4: bflag = 2
		else: bflag = 1

	# Print out star to file
	print >>out, "%16s %9.5f %9.5f   %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f   %2d %2d %6.3f %6.4f %6.3f %6.4f %5.2f   %6.3f %6.4f %5.2f" % (data[s][36], data[s][0], data[s][1], data[s][2], data[s][4], data[s][6], data[s][8], data[s][10], data[s][12], data[s][14], data[s][16], data[s][18], data[s][20], data[s][22], data[s][24], data[s][26], data[s][28], data[s][30], data[s][32], data[s][34], bflag, data[s][37], mpri, upri, msec, usec, medchi, smass, umass, smchi)
out.close()
print "\n\n"
