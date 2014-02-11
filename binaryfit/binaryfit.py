# Import modules
import numpy as np
import sys
import binocs
	
	
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

# Read in data from files
data = binocs.readdata(dataname)
oiso = binocs.readiso(isoname, age)

# Interpolate isochrone to new mass grid
singles = binocs.minterp(oiso, dm)

# Adjust isochrone if necessary
if len(fidcol) > 0: singles = binocs.fidiso(singles, fidname, fidmag, fidcol, ak, d, ebv, age, dataname)

# Create binary array
binary = binocs.makebin(singles, dataname, isoname, dm, age, fidcol)

# Run SED fitting on all stars in dataset
results = binocs.sedfit(singles, binary, data, nruns, d, ebv, ak)

# Now that we have completed all runs, we compute final results
# Open file for writing
print "\nWriting results to '%s--binary.txt'." % (dataname)
out = open(dataname+"--binary.txt", "w")

for s in range(len(data)):
	starchi = results[s, 1, :, 0]
	staridx = results[s, 0, :, 0]
	singlechi = results[s, 1, :, 1]
	singleidx = results[s, 0, :, 1]

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
	print >>out, "%16s %9.5f %9.5f   %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f   %2d %2d %6.3f %6.4f %6.3f %6.4f %5.2f   %6.3f %6.4f %5.2f" % (data[s][36], data[s][0], data[s][1], data[s][2], data[s][4], data[s][6], data[s][8], data[s][10], data[s][12], data[s][14], data[s][16], data[s][18], data[s][20], data[s][22], data[s][24], data[s][26], data[s][28], data[s][30], data[s][32], data[s][34], bflag, data[s][37], mpri, upri, msec, usec, medchi, smass, umass, smchi)

out.close()

print "\n\n"

















