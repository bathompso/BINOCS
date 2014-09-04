# BINOCS isochrone adjustment subroutines
from __future__ import print_function, division
import numpy as np
from scipy import interpolate
import os, sys

def minterp(oiso, dm):
	'''
	SUBROUTINE:			MINTERP
	DESCRIPTION: Interpolates isochrone onto a more fine mass grid
	INPUT:       oiso -- original isochrone
	             dm -- Mass increment between resulting isochrone points
	OUTPUT:      iiso -- interpolated isochrone matrix
	                0-1: Primary / Secondary Mass
	                2-5: Parameters: LogL, LogT, LogG, Mbol
	                6-22: Magnitudes
	'''
	
	print("\nInterpolating isochrone to %.3f M_sun..." % dm)
	
	# Find turnoff in B-V
	toindex, minbv, bindex = -1, 99, -1
	bv = oiso[:,7] - oiso[:,8]
	
	for c in range(len(bv)-1):
		if toindex >= 0: continue
			
		if bv[c] < minbv and bv[c] < 4:								# See if star is bluest yet
			minbv = bv[c]
			bindex = c
			
		if bv[c] - minbv >= 0.5: toindex = bindex			# Turn-off found: star is 0.5 mag redder than bluest
		elif oiso[c,0] >= 8.0: toindex = c-4						# Emergency break (mass > 8 M_sun)
		elif abs(oiso[c,0] - oiso[c+1,0]) < 0.0005: toindex = c-4	# Emergency break (mass degeneracy)
	# Turn-off is last datapoint if none was found
	if toindex < 0: toindex = len(bv)-1
		
	print("    Turn-Off Found: V = %.3f, B-V = %.4f" % (oiso[toindex,8], bv[toindex]))
	
	# Determine new mass points to interpolate to
	cmass = float(int(oiso[2,0] * 100.0)) / 100.0 + 0.01
	print("    New Grid: %.3f --> %.3f" % (cmass, oiso[toindex,1]))
	
	print("    Interpolating...", end='')
	sys.stdout.flush()
	newmass = []
	while cmass < oiso[toindex,1]:
		newmass.append(cmass)
		# Increment the mass counter for next pass
		if cmass < 1.0: cmass += dm / 2.0
		elif cmass < 1.5: cmass += dm
		else: cmass += 2.0 * dm
		
	# Create matrix for new grid
	iiso = np.zeros([len(newmass)+len(oiso)-toindex, oiso.shape[1]])
		
	# Interpolate parameters and magnitudes to new grid
	iiso[0:len(newmass),0] = newmass
	for i in range(2,oiso.shape[1]):
		fit = interpolate.interp1d(oiso[0:toindex+4,1], oiso[0:toindex+4,i], kind=3)
		interp_par = fit(newmass)
		iiso[0:len(newmass),i] = interp_par
			
	# Add original isochrone points beyond the turnoff to the isochrone
	iiso[len(newmass):,0] = oiso[toindex:,0]
	iiso[len(newmass):,2:] = oiso[toindex:,2:]
	print(" Done.")
	
	print("    Interpolated isochrone contains %d single stars." % iiso.shape[0])
	return iiso
	
	
	

def fidiso(iso, options, file_output=True):
	'''
	SUBROUTINE:			FIDISO
	DESCRIPTION: Adjusts isochrone data to empirical ridgelines
	INPUT:       iso -- isochrone data
	             options -- parameter dictionary from READOPT
	             file_output -- boolean flag to determine whether file with adjusted isochrone magnitude should be output
	OUTPUT:      fiso -- empirically-adjusted isochrone data
	                0-1: Primary / Secondary Mass
	                2-5: Parameters: LogL, LogT, LogG, Mbol
	                6-22: Magnitudes
	FILE OUTPUT: 'iso_[cluster].fid.dat' -- File containing the adjusted isochrone magnitudes. Same format as original isochrone file.
	'''
	
	# Check to see if this operation is necessary
	if 'fid' not in options.keys(): return iso

	# Create new list to hold adjusted isochrone
	fiso = np.zeros(iso.shape)
	fiso[:,0:6] = iso[:,0:6]
	fiso[:,6:] = 99.999
		
	# Read in fiducial file
	ff = open(options['fid'], "r")
	fidlines = ff.read().splitlines()
	ff.close()
	fiddata = [fidlines[x].replace('\t',' ').split() for x in range(1,len(fidlines))]
	
	# Read in fiducial magnitude and colors
	tmp = fidlines[0].replace('\t',' ').split()
	fidmag = [x for x in range(len(options['filternames'])) if options['filternames'][x] == tmp[0]][0]
	fidcol = []
	for f in range(1,len(tmp)):
		tmpcol = tmp[f].split('-')
		c = [x for x in range(len(options['filternames'])) if options['filternames'][x] == tmpcol[0]]
		m = [x for x in range(len(options['filternames'])) if options['filternames'][x] == tmpcol[1]]
		fidcol.append([c[0], m[0]])
	
	# Adjust fiducial data to absolute scale
	for l in range(len(fiddata)-1):
		# Adjust for distance and extinction
		fiddata[l][0] = float(fiddata[l][0]) - options['m-M'] - options['ak'][fidmag] / (options['ak'][1] - options['ak'][2]) * options['ebv']
		# Adjust for reddening
		for c in range(len(fidcol)): fiddata[l][c+1] = float(fiddata[l][c+1]) - (options['ak'][fidcol[c][0]] - options['ak'][fidcol[c][1]]) / (options['ak'][1] - options['ak'][2]) * options['ebv']
	
	# Fiducial magnitude values are assumed to be right
	fiso[:,fidmag+6] = iso[:,fidmag+6]
		
	# Loop through colors multiple times to adjust all necessary filters
	print("\nAdjusting isochrone to fiducial sequence...", end='')
	sys.stdout.flush()
	colcomplete = np.zeros(len(fidcol))
	for l in range(3):
		# Loop through all colors specified in fiducial file
		for c in range(len(fidcol)):
			if colcomplete[c] == 1: continue
			# Check to see if one of the magnitudes necessary has already been solved for.
			goodmag = len(fiso[fiso[:,fidcol[c][0]+6]<80, 0])
			goodcol = len(fiso[fiso[:,fidcol[c][1]+6]<80, 0])
			
			if goodmag == 0 and goodcol == 0: continue	# Neither magnitude has data, skip it
			elif goodmag > 0 and goodcol > 0:			# Both magnitudes have been completely solved
				colcomplete[c] = 1
				continue
				
			# Compute interpolation for colors
			datmag = [float(f[0]) for f in fiddata if float(f[c+1]) > -1000]
			datcol = [float(f[c+1]) for f in fiddata if float(f[c+1]) > -1000]
			fit = interpolate.interp1d(datmag, datcol, kind=3)
			
			# Magnitude filter solved for, but not color
			if goodmag > 0:
				for s in range(fiso.shape[0]):
					if fiso[s,fidmag+6] < min(datmag) or fiso[s,fidmag+6] > max(datmag): continue
					fiso[s,fidcol[c][1]+6] = fiso[s,fidcol[c][0]+6] - fit(fiso[s,fidmag+6])
				
			# Color filter solved for, but not magnitude
			if goodcol > 0:
				for s in range(fiso.shape[0]):
					if fiso[s,fidmag+6] < min(datmag) or fiso[s,fidmag+6] > max(datmag): continue
					fiso[s,fidcol[c][0]+6] = fiso[s,fidcol[c][1]+6] + fit(fiso[s,fidmag+6])
					
	# Loop through filters and complete any missing entries
	for f in range(6, 23):
		# Find all values where this magnitude is already solved for
		goodmag = [i for i in range(fiso.shape[0]) if fiso[i,f] < 80]
		
		# No fiducial for this filter, make is the same as the original
		if len(goodmag) == 0: fiso[:,f] = iso[:,f]
		else:
			# From the last index on, fill values
			for i in range(max(goodmag)+1, fiso.shape[0]):
				orig_diff = iso[i,f] - iso[i-1,f]
				fiso[i,f] = orig_diff + fiso[i-1,f]
			# From the first index and below, fill values
			for i in range(min(goodmag)-1, -1, -1):
				orig_diff = iso[i,f] - iso[i+1,f]
				fiso[i,f] = orig_diff + fiso[i+1,f]
	print(" Done.")
	
	if file_output:
		ndirsplit = options['data'].split('/')
		if len(ndirsplit) == 1: 
			ndirsplit = os.path.realpath(options['data']).split('/')
			fidoutname = "iso_%s.fid.dat" % (ndirsplit[-2])
		else: fidoutname = "%s/iso_%s.fid.dat" % ('/'.join(ndirsplit[0:len(ndirsplit)-1]), ndirsplit[len(ndirsplit)-2])
		fio = open(fidoutname, "w")
		for s in range(fiso.shape[0]):
			outstr = "%6.3f " % options['age']
			for i in range(6): outstr += "%7.4f " % fiso[s,i]
			for i in range(6,23): outstr += "%6.3f " % fiso[s,i]
			print(outstr, file=fio)
		fio.close()
		print("    Adjusted isochrone written to '%s'" % fidoutname)
	
	return fiso
	