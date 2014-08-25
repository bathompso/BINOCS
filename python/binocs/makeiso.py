# MAKEISO subroutines
from __future__ import print_function, division
import sys, subprocess
import numpy as np


def padova(path, outpath):
	'''
	SUBROUTINE:			PADOVA
	DESCRIPTION: Converts files downloaded from Padova's CMD web interface [http://stev.oapd.inaf.it/cgi-bin/cmd] to a usable format
	INPUT:       path -- Path of folder containing the downloaded Padova web files
	             outpath -- Path to folder to hold output
	OUTPUT:      NONE
	FILE OUTPUT: '[outpath]/iso_[FeH].pv.syn.dat' -- File holding isochrone star information to be read into BINOCS
	                 0: log[Age] of isochrone
	                 1: Initial mass
	                 2: Actual mass (at specified age)
	                 3: log[Luminosity]
	                 4: log[g] (surface gravity)
	                 5: log[Temperature]
	                 6: Bolometric magnitude
	                 7-23: UBVRIugrizJHK[3][4][5][8] magnitudes
	'''
	# Detect what files are present in path directory
	webfiles = subprocess.check_output("ls "+path+"*.dat", shell=True).splitlines()
	# Loop through detected files and get [Fe/H]
	webfeh, nlines = [], []
	for f in range(len(webfiles)):
		df = open(webfiles[f], 'r')
		lines = df.read().splitlines()
		df.close()
		# Determine what line to read in (changes if SDSS filter file)
		if lines[3].find('SDSS') >= 0: tmp = lines[11].split()
		else: tmp = lines[10].split()
		# Save [Fe/H] value for this file, which is scaled from Z
		webfeh.append(np.log10(float(tmp[4]) / 0.01886))
		nlines.append(len(lines))
	# Find all unique [Fe/H] values, and print out formatted isochrone file
	for u in np.unique(webfeh):
		thisuni = [x for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015]
		thisnlines = max([nlines[x] for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015])
		# If we have all three types of files, we can print an output for this [Fe/H]
		if len < 3: continue
		# Determine output file name
		if u < 0: outname = "%s/iso_m%03d.pv.syn.dat" % (outpath, -1.0*u*100.0)
		else: outname = "%s/iso_p%03d.pv.syn.dat" % (outpath, u*100.0)
		print("Printing isochrone for [Fe/H] = %5.2f to '%s'" % (u, outname))
		# Loop through all webfiles for this [Fe/H] and read in data
		data = np.zeros([thisnlines, 24])
		for f in thisuni:
			df = open(webfiles[f], 'r')
			lines = df.read().splitlines()
			df.close()
			# Determine what file type this is
			if lines[3].find('SDSS') >= 0:
				print("    Reading SDSS+JHK file '%s'" % (webfiles[f]))
				adj = 1
				ftype = 2
			elif lines[11].find('V') >= 0:
				print("    Reading UBVRI file '%s'" % (webfiles[f]))
				adj = 0
				ftype = 1
			else:
				print("    Reading IRAC file '%s'" % (webfiles[f]))
				adj = 0
				ftype = 3
			for i in range(len(lines)):
				if lines[i].find('#') >= 0: continue
				tmp = lines[i].split()
				if len(tmp) == 0: continue
				# Save parameters to array
				for j in range(7): data[i-adj,j] = float(tmp[j])
				# Save magnitudes to array
				if ftype == 1:
					for j in range(7, 12): data[i-adj,j] = float(tmp[j])
				elif ftype == 2:
					for j in range(7, 15): data[i-adj,j+5] = float(tmp[j])
				else:
					for j in range(7, 11): data[i-adj,j+13] = float(tmp[j])
		# Print out newly matched file
		of = open(outname, 'w')
		for s in range(thisnlines):
			# Check to see whether all magnitudes exist
			badmag = [x for x in data[s,:] if x == 0 or x < -9.9]
			if len(badmag) > 0: continue
			# Print out star
			print("%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (data[s,0], data[s,1], data[s,2], data[s,3], data[s,4], data[s,5], data[s,6], data[s,7], data[s,8], data[s,9], data[s,10], data[s,11], data[s,12], data[s,13], data[s,14], data[s,15], data[s,16], data[s,17], data[s,18], data[s,19], data[s,20], data[s,21], data[s,22], data[s,23]), file=of)
		of.close()




def parsec(path, outpath):
	'''
	SUBROUTINE:			PARSEC
	DESCRIPTION: Converts files downloaded from PARSEC's CMD web interface [http://stev.oapd.inaf.it/cgi-bin/cmd] to a usable format
	INPUT:       path -- Path of folder containing the downloaded PARSEC web files
	             outpath -- Path to folder to hold output
	OUTPUT:      NONE
	FILE OUTPUT: '[outpath]/iso_[FeH].pc.syn.dat' -- File holding isochrone star information to be read into BINOCS
	                 0: log[Age] of isochrone
	                 1: Initial mass
	                 2: Actual mass (at specified age)
	                 3: log[Luminosity]
	                 4: log[g] (surface gravity)
	                 5: log[Temperature]
	                 6: Bolometric magnitude
	                 7-23: UBVRIugrizJHK[3][4][5][8] magnitudes
	'''
	# Detect what files are present in path directory
	webfiles = subprocess.check_output("ls "+path+"*.dat", shell=True).splitlines()
	# Loop through detected files and get [Fe/H]
	webfeh, nlines = [], []
	for f in range(len(webfiles)):
		df = open(webfiles[f], 'r')
		lines = df.read().splitlines()
		df.close()
		# Save [Fe/H] value for this file
		tmp = lines[11].split()
		webfeh.append(float(tmp[10]))
		nlines.append(len(lines))
	# Find all unique [Fe/H] values, and print out formatted isochrone file
	for u in np.unique(webfeh):
		thisuni = [x for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015]
		thisnlines = max([nlines[x] for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015])
		# If we have all three types of files, we can print an output for this [Fe/H]
		if len < 3: continue
		# Determine output file name
		if u < 0: outname = "%s/iso_m%03d.pc.syn.dat" % (outpath, -1.0*u*100.0)
		else: outname = "%s/iso_p%03d.pc.syn.dat" % (outpath, u*100.0)
		print("Printing isochrone for [Fe/H] = %5.2f to '%s'" % (u, outname))
		# Loop through all webfiles for this [Fe/H] and read in data
		data = np.zeros([thisnlines, 24])
		for f in thisuni:
			df = open(webfiles[f], 'r')
			lines = df.read().splitlines()
			df.close()
			# Determine what file type this is
			if lines[12].find('Ks') >= 0:
				print("    Reading SDSS+JHK file '%s'" % (webfiles[f]))
				ftype = 2
			elif lines[12].find('V') >= 0:
				print("    Reading UBVRI file '%s'" % (webfiles[f]))
				ftype = 1
			else:
				print("    Reading IRAC file '%s'" % (webfiles[f]))
				ftype = 3
			for i in range(len(lines)):
				if lines[i].find('#') >= 0: continue
				tmp = lines[i].split()
				if len(tmp) == 0: continue
				# Save parameters to array
				for j in range(7): data[i,j] = float(tmp[j+1])
				# Save magnitudes to array
				if ftype == 1:
					for j in range(7, 12): data[i,j] = float(tmp[j+1])
				elif ftype == 2:
					for j in range(7, 15): data[i,j+5] = float(tmp[j+1])
				else:
					for j in range(7, 11): data[i,j+13] = float(tmp[j+1])
		# Print out newly matched file
		of = open(outname, 'w')
		for s in range(thisnlines):
			# Check to see whether all magnitudes exist
			badmag = [x for x in data[s,:] if x == 0 or x < -9.9]
			if len(badmag) > 0: continue
			# Print out star
			print("%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (data[s,0], data[s,1], data[s,2], data[s,3], data[s,4], data[s,5], data[s,6], data[s,7], data[s,8], data[s,9], data[s,10], data[s,11], data[s,12], data[s,13], data[s,14], data[s,15], data[s,16], data[s,17], data[s,18], data[s,19], data[s,20], data[s,21], data[s,22], data[s,23]), file=of)
		of.close()




def dartmouth(path, outpath):
	'''
	SUBROUTINE:			DARTMOUTH
	DESCRIPTION: Converts files downloaded from Dartmouth's web interface [http://stellar.dartmouth.edu/models/isolf_new.html] to a usable format
	INPUT:       path -- Path of folder containing the downloaded Dartmouth web files
	             outpath -- Path to folder to hold output
	OUTPUT:      NONE
	FILE OUTPUT: '[outpath]/iso_[FeH].dm.syn.dat' -- File holding isochrone star information to be read into BINOCS
	                 0: log[Age] of isochrone
	                 1: Initial mass
	                 2: Actual mass (at specified age)
	                 3: log[Luminosity]
	                 4: log[g] (surface gravity)
	                 5: log[Temperature]
	                 6: Bolometric magnitude
	                 7-23: UBVRIugrizJHK[3][4][5][8] magnitudes
	'''
	# Detect what files are present in path directory
	webfiles = subprocess.check_output("ls "+path+"*.iso", shell=True).splitlines()
	# Loop through detected files and get [Fe/H]
	webfeh, nlines = [], []
	for f in range(len(webfiles)):
		df = open(webfiles[f], 'r')
		lines = df.read().splitlines()
		df.close()
		# Save [Fe/H] value for this file
		tmp = lines[3].split()
		webfeh.append(float(tmp[5]))
		nlines.append(len(lines))
	# Find all unique [Fe/H] values, and print out formatted isochrone file
	for u in np.unique(webfeh):
		thisuni = [x for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015]
		thisnlines = max([nlines[x] for x in range(len(webfeh)) if np.abs(webfeh[x] - u) < 0.015])
		# If we have all three types of files, we can print an output for this [Fe/H]
		if len < 3: continue
		# Determine output file name
		if u < 0: outname = "%s/iso_m%03d.dm.syn.dat" % (outpath, -1.0*u*100.0)
		else: outname = "%s/iso_p%03d.dm.syn.dat" % (outpath, u*100.0)
		print("Printing isochrone for [Fe/H] = %5.2f to '%s'" % (u, outname))
		# Loop through all webfiles for this [Fe/H] and read in data
		data = np.zeros([thisnlines, 24])
		for f in thisuni:
			df = open(webfiles[f], 'r')
			lines = df.read().splitlines()
			df.close()
			# Determine what file type this is
			if lines[5].find('SDSS') >= 0:
				print("    Reading SDSS file '%s'" % (webfiles[f]))
				ftype = 1
			elif lines[5].find('Bessel') >= 0:
				print("    Reading UBVRI+JHK file '%s'" % (webfiles[f]))
				ftype = 2
			else:
				print("    Reading IRAC file '%s'" % (webfiles[f]))
				ftype = 3
			for i in range(len(lines)):
				if lines[i].find('AGE') == 1: thisage = np.log10(float((lines[i])[5:11])*1E9)
				if lines[i].find('#') >= 0: continue
				tmp = lines[i].split()
				if len(tmp) == 0: continue
				# Save parameters to array 
				data[i,0] = thisage
				data[i,1], data[i,2] = float(tmp[1]), float(tmp[1])
				data[i,3], data[i,4], data[i,5] = float(tmp[4]), float(tmp[2]), float(tmp[3])		# LogL, LogT, LogG
				data[i,6] = -2.5 * data[i,3] + 4.75      											# Bolometric Magnitude
				# Save magnitudes to array
				if ftype == 1:
					for j in range(5, 10): data[i,j+7] = float(tmp[j])
				elif ftype == 2:
					for j in range(5, 10): data[i,j+2] = float(tmp[j])
					for j in range(10,13): data[i,j+7] = float(tmp[j])
				else:
					for j in range(5, 9): data[i,j+15] = float(tmp[j])
		# Print out newly matched file
		of = open(outname, 'w')
		for s in range(thisnlines):
			# Check to see whether all magnitudes exist
			badmag = [x for x in data[s,:] if x == 0 or x < -9.9]
			if len(badmag) > 0: continue
			# Print out star
			print("%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (data[s,0], data[s,1], data[s,2], data[s,3], data[s,4], data[s,5], data[s,6], data[s,7], data[s,8], data[s,9], data[s,10], data[s,11], data[s,12], data[s,13], data[s,14], data[s,15], data[s,16], data[s,17], data[s,18], data[s,19], data[s,20], data[s,21], data[s,22], data[s,23]), file=of)
		of.close()
		

		