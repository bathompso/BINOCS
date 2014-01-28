# Import modules
import pyopencl as cl
import matplotlib.pyplot as plt
import sys, math, subprocess
import numpy as np
from time import time


# Coordinate matching kernel
cookernel = """__kernel void coomatch(
	__global float* ra,
	__global float* dec,
	__global float* thesera,
	__global float* thesedec,
	__global int* indexes,
	__global int* matched,
	const float radius,
	const int length)
{
	float dra, ddec, diff;
	int match = -1;
	float mindiff = 99.0;
	int j = get_global_id(0);
	for (int i = 0; i < length; i++){
		if (i == j){ continue; }
		dra = (ra[i] - thesera[j]) * cos(thesedec[j] * 3.14159265 / 180.0);
		ddec = dec[i] - thesedec[j];
		diff = sqrt( dra*dra + ddec*ddec );
		if (diff < radius && diff < mindiff){
			mindiff = diff;
			match = i;
		}
	}
	indexes[j] = match;
}
"""


# Define function to automatically extend every array we have to deal with
def extendarr( id2mass, mag, mempct, memchar, match ):
	id2mass.append('00000000+0000000')
	mempct.append(-1)
	memchar.append('U')
	match.append(-1)
	

# Define function to check for magnitude replacement
def savemag( dataarr, mag, mempct, memchar, index, filters, maxerr ):
	for f in range(len(filters)):
		# Check to see if error on magnitude is worth saving
		if float(dataarr[2*f+1]) > maxerr and filters[0] < 38: continue
		# Check to see whether errors on the new value are less than what is currently there
		if float(dataarr[2*f+1]) < mag[index, filters[f]+1]:
			mag[index, filters[f]] = float(dataarr[2*f])
			mag[index, filters[f]+1] = float(dataarr[2*f+1])
	# If this is a PM file, we need to add member information
	if filters[0] == 40:
		mempct[index] = float(dataarr[4])
		memchar[index] = dataarr[5]
	# If this is a RV file, we need to add member information
	if filters[0] == 38:
		mempct[index] = float(dataarr[2])
		memchar[index] = dataarr[3]

	
# Define variables
maxlines = 0
nfiles = 0
ctr = 0
oldctr = 0
maxerr = 0.1
masterfilters = ['U', 'B', 'V', 'R', 'I', 'SU', 'SG', 'SR', 'SI', 'SZ', 'J', 'H', 'K', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6']
ak = [1.531, 1.324, 1.000, 0.748, 0.482, 1.593, 1.199, 0.858, 0.639, 0.459, 0.282, 0.175, 0.112, 0.0627, 0.0482, 0.0482, 0.0482]

print "\n\n"
print "==========================================="
print "|                  PAYST                  |"
print "|                  v2.0                   |"
print "|    GPU-Enabled (Requires OpenCL 1.2+)   |"
print "|                                         |"
print "|             by Ben Thompson             |"
print "==========================================="


### Create OpenCL queue
context = cl.create_some_context()							# Asks user to specify which OpenCL device should be used.
queue = cl.CommandQueue(context)							# Creates a command queue for specified device
program = cl.Program(context, cookernel).build()			# Compiles kernel program on specified device
coomatch = program.coomatch
coomatch.set_scalar_arg_dtypes([None, None, None, None, None, None, np.float32, np.uint32])


# Read options from command line
comargs = sys.argv


### Catalog matching mode
if len(comargs) > 2:
	optname = comargs[1]
	minradius = float(comargs[2]) / 3600.0
	print "\n!!!  Parameters:   Option File =", optname, "    radius =", minradius*3600.0, "arcsec"
	optsplit = optname.split('/')
	if len(optsplit) > 1: prefix = "/".join(optsplit[0:len(optsplit)-1]) + "/"
	else: prefix = ""
	if len(comargs) > 3: maxerr = float(comargs[3])

	# Read parameters from option file
	optfile = open(optname, 'r')
	optlines = optfile.read().splitlines()
	optfile.close()
	for line in optlines:
		if '#' in line:
			continue
		line_split = line.split()
		if len(line_split) == 0: continue
		thisfile = open(prefix+line_split[0], 'r')
		thislines = thisfile.read().splitlines()
		thisfile.close()
		maxlines += len(thislines)
		nfiles += 1
	print "!!!     Matching", maxlines, "stars in", nfiles, "files."
	print ""

	# Generate all arrays & generic entries
	id2mass = []
	ra = np.zeros(maxlines).astype(np.float32)
	dec = np.zeros(maxlines).astype(np.float32)
	mag = np.zeros([maxlines, 44])
	for l in range(maxlines):
		for e in range(44):
			if e % 2 == 0: mag[l,e] = 99.999
			else: mag[l,e] = 9.999
	mempct = []
	memchar = []
	match = []


	# Loop through option file and match everything
	start = time()
	for o in range(len(optlines)):
		starttime = time()

		if '#' in optlines[o]:
			continue
		tmp = optlines[o].split()
		if len(tmp) == 0: continue
		f = open(prefix+tmp[0], 'r')
		datalines = f.read().splitlines()
		f.close()
	
		thisjnk = int(tmp[1])
		thisfiltchars = tmp[2:len(tmp)]
		print "Opening", tmp[0], "..."
	
		if tmp[0].find('2MASS') >= 0: is2mass = 1
		else: is2mass = 0
	
		RVPM = 0
		if len(thisfiltchars) > 0:
			print "     ", len(datalines), "stars in filters:", ' '.join(thisfiltchars)
			filters = []
			for f in thisfiltchars:    # Loop through and find filter indices
				fin = [j for j in range(len(masterfilters)) if f == masterfilters[j]]
				if len(fin) == 0:
					continue
				filters.append(2*fin[0])
		elif 'RV' in tmp[0]:
			RVPM = 1
			print "     ", len(datalines), "stars with RV data."
			filters = [38]
		elif 'PM' in tmp[0]:
			RVPM = 2
			print "     ", len(datalines), "stars with PM data."
			filters = [40, 42]
		else:
			print "Unknown file, please rename."
			exit()
		
		#### SPECIAL LOGIC FOR FIRST RUN
		if len(id2mass) < 1:
			for l in datalines:
				tmp = l.split()
				thisra = float(tmp[thisjnk])
				thisdec = float(tmp[thisjnk+1])
				# Add member info to array if PM data file, if not, just need mags
				if RVPM == 2:
					dataarr = tmp[thisjnk+2:2*len(filters)+thisjnk+4]
					if len([x for x in range(len(dataarr)/2-1) if float(dataarr[2*x+1]) < maxerr]) == 0: continue
				else:
					dataarr = tmp[thisjnk+2:2*len(filters)+thisjnk+2]
					if len([x for x in range(len(dataarr)/2) if float(dataarr[2*x+1]) < maxerr]) == 0: continue
				extendarr(id2mass, mag, mempct, memchar, match)
				ra[len(id2mass)-1] = thisra
				dec[len(id2mass)-1] = thisdec
				if is2mass == 1: id2mass[len(id2mass)-1] = tmp[0]
				savemag(dataarr, mag, mempct, memchar, len(id2mass)-1, filters, maxerr)
				match[len(id2mass)-1] = o
			print "      Added", len(id2mass), "stars."
			print "      Elapsed: %d seconds." % (time() - start)
			print ""
			continue
		
		#### LOGIC FOR NON-FIRST RUNS		
		# Extract coordinates for this file
		thisra = np.empty(len(datalines)).astype(np.float32)
		thisdec = np.empty(len(datalines)).astype(np.float32)
		for l in range(len(datalines)):
			tmp = datalines[l].split()
			thisra[l] = float(tmp[thisjnk])
			thisdec[l] = float(tmp[thisjnk+1])
		indexes = np.empty(len(thisra)).astype(np.int32)
		
		# Copy data to device for matching
		d_ra = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=ra)
		d_dec = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=dec)
		d_tra = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=thisra)
		d_tdec = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=thisdec)
		d_matched = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=indexes)
		d_index = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, indexes.nbytes)
	
		coomatch(queue, thisra.shape, None, d_ra, d_dec, d_tra, d_tdec, d_index, d_matched, minradius, len(id2mass))
	
		queue.finish()
		cl.enqueue_copy(queue, indexes, d_index)
		matched = ((indexes >= 0) & (indexes <= len(ra)))
	
		print "      Matched", len(indexes[matched]), "stars."
	
		# Save magnitudes
		for i in range(len(indexes)):
			tmp = datalines[i].split()
	
			# Add member info to array if RV/PM data file, if not, just need mags
			if RVPM > 0:
				dataarr = tmp[thisjnk+2:2*len(filters)+thisjnk+4]
			else:
				dataarr = tmp[thisjnk+2:2*len(filters)+thisjnk+2]
				if len([x for x in range(len(dataarr)/2) if float(dataarr[2*x+1]) < maxerr]) == 0: continue
	
			# If we matched another star
			if (indexes[i] >= 0): savemag(dataarr, mag, mempct, memchar, indexes[i], filters, maxerr)
			else:
				extendarr(id2mass, mag, mempct, memchar, match)
				ra[len(id2mass)-1] = thisra[i]
				dec[len(id2mass)-1] = thisdec[i]
				savemag(dataarr, mag, mempct, memchar, len(id2mass)-1, filters, maxerr)
	
		print "      Elapsed: %d seconds." % (time() - start)
		print ""
	
	# End loop through option file

	# Print out results
	print "!!! WRITING RESULTS TO OUTPUT FILE"
	namesplit = optsplit[len(optsplit)-1].split('.')
	out = open(prefix+namesplit[0]+".Merged.txt", 'w')
	for i in range(len(id2mass)):
		outstring = ""
		# If 2MASS ID does not exist, make a dummy one
		if id2mass[i] == '00000000+0000000':
			rahrs = ra[i] / 15.0
			rah = int(rahrs)
			ram = int((rahrs - rah) * 60.0)
			ras = int((((rahrs - rah) * 60.0) - ram) * 6000.0)
			decd = int(dec[i])
			decm = int((dec[i] - decd) * 60.0)
			decs = int((((dec[i] - decd) * 60.0) - decm) * 600.0)
			if dec[i] > 0: id2mass[i] = 'BT%02d%02d%04d+%02d%02d%03d' % (rah, ram, ras, decd, decm, decs)
			else: id2mass[i] = 'BT%02d%02d%04d-%02d%02d%03d' % (rah, ram, ras, abs(decd), decm, decs)
		else: id2mass[i] = '2M%s' % (id2mass[i])
		outstring += "%16s " % (id2mass[i])					# Add 2MASS ID to output
		outstring += "%9.5f %9.5f " % (ra[i], dec[i])		# Add coordinates to output
		for f in range(19):									# Add filter magnitudes to output
			outstring += "%6.3f %5.3f " % (mag[i, 2*f], mag[i, 2*f+1])
		for v in range(6):									# Add kinematic information to output
			outstring += "%8.3f " % (mag[i, v+38])
		outstring += "%5d " % (mempct[i])					# Add member information to output
		outstring += "%5s " % (memchar[i])
		print >>out, outstring
		
	out.close()	
		
		
### Catalog trimming mode
else:
	catalog = comargs[1]
	print "\n!!!  Parameters:   Catalog =",catalog
	catf = open(catalog, "r")
	lines = catf.read().splitlines()
	catf.close()
	print "!!!     Catalog contains", len(lines), "stars."
	
	# Create arrays
	id2mass = []
	ra, dec = np.zeros(len(lines)), np.zeros(len(lines))
	mags = np.zeros([len(lines), 44])
	mempct = np.zeros(len(lines))
	memchar = []
	member = np.zeros(len(lines))
	
	# Loop through file and read in data
	for l in range(len(lines)):
		tmp = lines[l].split()
		id2mass.append(tmp[0])
		ra[l], dec[l] = float(tmp[1]), float(tmp[2])
		for i in range(44): mags[l,i] = float(tmp[i+3])
		mempct[l] = float(tmp[47])
		memchar.append(tmp[48])
		member[l] = 1
		
	# Begin Menu Loop
	choice, isoplot = 100, 0
	pmag, pcol = -1, [0, 0]
	while choice != 0:
		# Check to see whether CMD can be plotted
		if pmag < 0:
			# Ask for new CMD choice
			pmagstr = raw_input('\nWhat is the CMD magnitude? ')
			pcolstr = raw_input('What is the CMD color? ')
			pmag = ([x for x in range(len(masterfilters)) if pmagstr == masterfilters[x]])[0]
			pcol[0] = ([x for x in range(len(masterfilters)) if (pcolstr.split('-'))[0] == masterfilters[x]])[0]
			pcol[1] = ([x for x in range(len(masterfilters)) if (pcolstr.split('-'))[1] == masterfilters[x]])[0]
			
		# Plot CMD
		cmdmag = [mags[x, 2*pmag] for x in range(len(lines)) if mags[x, 2*pmag] < 80 and mags[x, 2*pcol[0]] < 80 and mags[x, 2*pcol[1]] < 80 and member[x] == 1]
		cmdcol = [mags[x, 2*pcol[0]] - mags[x, 2*pcol[1]] for x in range(len(lines)) if mags[x, 2*pmag] < 80 and mags[x, 2*pcol[0]] < 80 and mags[x, 2*pcol[1]] < 80  and member[x] == 1]
		plt.clf()
		plt.plot(cmdcol, cmdmag, "ko", markersize=1)
		if isoplot == 1:
			isomag = [isodata[x,pmag+1] + isod + ak[pmag]/0.324*isoebv for x in range(len(isolines)) if isodata[x,0] == isoage]
			isocol = [isodata[x,pcol[0]+1] - isodata[x,pcol[1]+1] + (ak[pcol[0]] - ak[pcol[1]])/0.324*isoebv for x in range(len(isolines)) if isodata[x,0] == isoage]
			plt.plot(isocol, isomag, "b-")
		plt.axis([min(cmdcol)-0.1, max(cmdcol)+0.1, max(cmdmag)+0.5, min(cmdmag)-0.5])
		plt.ylabel(pmagstr)
		plt.xlabel(pcolstr)
		plt.show(block=False)
		plt.draw()
			
		# Print Menu
		print "\n\n 1) RA / Dec Cut"
		print " 2) Membership Cut"
		print " 3) A_K Cut"
		print " 4) Full Photometry Cut"
		print ""
		print " 8) Change CMD Options"
		print " 9) Reset Trimming"
		print " 0) Exit and Print Out Results"
		print "-1) Abort"
		choice = input(': ')
		print "\n"
		
		# Emergency Break
		if choice == -1: sys.exit(0)
		
		# Spatial RA/Dec Trim
		if choice == 1:
			# Pre-determined centers for all popular clusters
			clustnames = ['Berkeley39', 'Collinder261', 'IC4651', 'Melotte66', 'NGC1039', 'NGC1817', 'NGC188', 'NGC1912', 'NGC1960', 'NGC2099', 'NGC2158', 'NGC2168', 'NGC2420', 'NGC2477', 'NGC2682', 'NGC3680', 'NGC6791', 'NGC6819']
			clustras = [116.675, 189.488, 261.204, 111.596, 40.521, 78.063, 11.867, 82.167, 84.075, 88.075, 91.854, 92.250, 114.596, 118.042, 132.825, 171.408, 290.221, 295.325]
			clustdecs = [-4.600, -68.367, -49.933, -47.667, 42.762, 16.690, 85.255, 35.848, 34.140, 32.553, 24.097, 24.350, 21.573, -38.530, 11.800, -43.243, 37.772, 40.187]
			# Figure out whether cluster has a known location
			clidx = [x for x in range(len(clustnames)) if catalog.find(clustnames[x]) >= 0]
			if len(clidx) == 0:
				print "Unknown Cluster."
				clra = input('Enter Cluster RA: ')
				cldec = input('Enter Cluster Dec: ')
			else:
				clra = clustras[clidx[0]]
				cldec = clustdecs[clidx[0]]
				print "Cluster %s identified at RA = %.1f, DEC = %.1f" % (clustnames[clidx[0]], clra, cldec)
				
			# Generate Spatial Plot for Cluster
			cutset = 0
			while cutset != 1:
				gra = [ra[x] for x in range(len(lines)) if member[x] == 1]
				gdec = [dec[x] for x in range(len(lines)) if member[x] == 1]
				plt.clf()
				plt.plot(gra, gdec, "ko", markersize=1)
				plt.plot(clra, cldec, "bo", markersize=4)
				if max(gdec) - min(gdec) > 2: plt.axis([clra-1, clra+1, cldec-1, cldec+1])
				else: plt.axis([min(gra)-0.1, max(gra)+0.1, min(gdec)-0.1, max(gdec)+0.1])
				plt.show(block=False)
				plt.draw()
				
				# Ask for trimming radius
				clrad = input('Enter cluster radius (arcmin): ') / 60.0
				
				# Show trimming radius on plot
				radra, raddec = [], []
				for e in range(100):
					angle = 2.0*math.pi * (float(e)/100.0)
					radra.append(clra + np.cos(angle) * clrad)
					raddec.append(cldec + np.sin(angle) * clrad)
				plt.clf()
				plt.plot(gra, gdec, "ko", markersize=1)
				plt.plot(clra, cldec, "bo", markersize=4)
				plt.plot(radra, raddec, "b-")
				if max(gdec) - min(gdec) > 2: plt.axis([clra-1, clra+1, cldec-1, cldec+1])
				else: plt.axis([min(gra)-0.1, max(gra)+0.1, min(gdec)-0.1, max(gdec)+0.1])
				plt.show(block=False)
				plt.draw()
				
				cutchoice = raw_input('Radius ok? (y|n)  ')
				if cutchoice == 'y': cutset = 1
				
			# Loop through and deselect stars outside radius
			trimmed = 0
			for s in range(len(ra)):
				if member[s] == 0: continue
				dist = np.sqrt((ra[s] - clra)**2 + (dec[s] - cldec)**2)
				if dist > clrad:
					trimmed += 1
					member[s] = 0
			print "REMOVED", trimmed, "stars."
			
		# Membership Selection
		if choice == 2:
			print "Qualitative Selection:"
			print "  1) Select only Members"
			print "  2) Select only Single Members"
			print "  3) Deselect only Non-Members"
			print "  4) No Selection"
			charchoice = input(': ')
			print "\nQuantitative Selection:"
			print "  1) Select only Stars > %"
			print "  2) Deselect only Stars < %"
			print "  3) No Selection"
			pctchoice = input(': ')
			if pctchoice < 3: pctcut = input('Enter % cutoff: ')
			
			# Loop through stars and select necessary stars
			trimmed = 0
			for s in range(len(member)):
				if member[s] == 0: continue
				ccut, pcut = 0, 0
				# Go through qualitative selection criteria
				if charchoice == 1:
					if memchar[s].find('M') < 0: ccut = 1
				elif charchoice == 2:
					if memchar[s].find('M') < 0 or memchar[s].find('S') < 0: ccut = 1
				elif charchoice == 3:
					if memchar[s].find('N') >= 0: ccut = 1
				else:
					ccut = -1
				# Go through quantitative selection criteria
				if pctchoice == 1:
					if mempct[s] < pctcut: pcut = 1
				elif pctchoice == 2:
					if mempct[s] < pctcut and mempct[s] >= 0: pcut = 1
				else:
					pcut = -1
				# Make choice to trim or not
				if ccut > 0:
					member[s] = 0
					trimmed += 1
				elif pcut > 0:
					member[s] = 0
					trimmed += 1
					
			print "REMOVED",trimmed,"stars."
			
		# A_K Cut
		#if choice == 3:
					
		# Photometry Cut
		if choice == 4:
			nvis, nnir, nmir = np.zeros(len(ra)), np.zeros(len(ra)), np.zeros(len(ra))
			nfilt = np.zeros(16)
			for s in range(len(ra)):
				# Calculate number of good visual filters
				gjc = [x for x in range(5) if mags[s, 2*x] < 80]
				gtg = [x for x in range(5,10) if mags[s, 2*x] < 80]
				# N_VIS is maximum of either visual set
				nvis[s] = max([len(gjc), len(gtg)])
				
				# Calculate number of good 2MASS filters
				g2m = [x for x in range(10,13) if mags[s, 2*x] < 80]
				nnir[s] = len(g2m)
				
				# Calculate number of good IRAC filters
				gir = [x for x in range(13, 17) if mags[s, 2*x] < 80]
				nmir[s] = len(gir)
				
				# Adjust count for filters
				nfilt[nvis[s]+nnir[s]+nmir[s]] += 1
				# Adjust for SED filter combinations
				if nvis[s] >= 3 and nnir[s] >= 3 and nmir[s] >= 3: nfilt[13] += 1
				elif nvis[s] >= 3 and nnir[s] >= 3 and nmir[s] >= 2: nfilt[14] += 1
				elif nvis[s] >= 3 and nnir[s] >= 2 and nmir[s] >= 2: nfilt[15] += 1
				
			# Print out options for filter counts
			print "# Filters    # Stars"
			for f in range(13):
				if nfilt[f] == 0: continue
				print "   %3d       %7d" % (f, nfilt[f])
			print "\nSED FITTING"
			print "   %3d       %7d" % (333, nfilt[13])
			print "   %3d       %7d" % (332, nfilt[14])
			print "   %3d       %7d" % (322, nfilt[15])
			minfilt = input('\nMinimum number of filters: ')
			
			# Trim stars
			trimmed = 0
			for s in range(len(ra)):
				if minfilt < 100:
					if nvis[s]+nnir[s]+nmir[s] < minfilt:
						member[s] = 0
						trimmed += 1
				elif minfilt == 333:
					if nvis[s] < 3 or nnir[s] < 3 or nmir[s] < 3:
						member[s] = 0
						trimmed += 1
				elif minfilt == 332:
					if nvis[s] < 3 or nnir[s] < 3 or nmir[s] < 2:
						member[s] = 0
						trimmed += 1
				elif minfilt == 322:
					if nvis[s] < 3 or nnir[s] < 2 or nmir[s] < 2:
						member[s] = 0
						trimmed += 1
			print "REMOVED",trimmed,"stars."
			
		# Isochrone overplotting
		if choice == 5:
			if isoplot == 0:
				isonames = subprocess.check_output('ls ~/projects/isochrones/interp/*.dat', shell=True).splitlines()
				isoname = raw_input('Enter isochrone name to overplot: ')
				isoidx = [x for x in range(len(isonames)) if isonames[x].find(isoname) >= 0]
				isof = open(isonames[isoidx[0]], "r")
				isolines = isof.read().splitlines()
				isof.close()
				isodata = np.zeros([len(isolines), 18])
				for l in range(len(isolines)):
					tmp = [float(x) for x in isolines[l].split()]
					isodata[l,0] = tmp[0]
					for j in range(1,18): isodata[l,j] = tmp[j+6]
				isoplot = 1
			else:
				removechoice = raw_input('Adjust parameters? (y|n): ')
				if removechoice == 'n':
					isoplot = 0
					print "Removing isochrone ridgeline."
					continue
			isoage = input('Enter isochrone age: ')
			isod = input('Enter isochrone m-M: ')
			isoebv = input('Enter isochrone E(B-V): ')
			
		# Change CMD options
		if choice == 8:
			pmag = -1
			
		# Reset trimming options
		if choice == 9:
			for s in range(len(member)): member[s] = 1
			print "RESTORED", len(member), "stars."
			
	# Now that trimming is done, print out everything
	catsplit = catalog.split('/')
	if len(catsplit) > 1: prefix = "/".join(catsplit[0:len(catsplit)-1]) + "/"
	else: prefix = ""
	
	print "!!! WRITING RESULTS TO OUTPUT FILE"
	namesplit = catsplit[len(catsplit)-1].split('.')
	out = open(prefix+namesplit[0]+".Trimmed.txt", 'w')
	for i in range(len(id2mass)):
		if member[i] == 0: continue
		outstring = ""
		outstring += "%16s " % (id2mass[i])					# Add 2MASS ID to output
		outstring += "%9.5f %9.5f " % (ra[i], dec[i])		# Add coordinates to output
		for f in range(19):									# Add filter magnitudes to output
			outstring += "%6.3f %5.3f " % (mags[i, 2*f], mags[i, 2*f+1])
		for v in range(6):									# Add kinematic information to output
			outstring += "%8.3f " % (mags[i, v+38])
		outstring += "%5d " % (mempct[i])					# Add member information to output
		outstring += "%5s " % (memchar[i])
		print >>out, outstring
	out.close()	
		
		