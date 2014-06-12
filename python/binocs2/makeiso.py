# MAKEISO subroutines
import sys, subprocess
import numpy as np


def padova(path, outpath):
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
		print "Printing isochrone for [Fe/H] = %5.2f to '%s'" % (u, outname)
		# Loop through all webfiles for this [Fe/H] and read in data
		data = np.zeros([thisnlines, 24])
		for f in thisuni:
			df = open(webfiles[f], 'r')
			lines = df.read().splitlines()
			df.close()
			# Determine what file type this is
			if lines[3].find('SDSS') >= 0:
				print "    Reading SDSS+JHK file '%s'" % (webfiles[f])
				adj = 1
				ftype = 2
			elif lines[11].find('V') >= 0:
				print "    Reading UBVRI file '%s'" % (webfiles[f])
				adj = 0
				ftype = 1
			else:
				print "    Reading IRAC file '%s'" % (webfiles[f])
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
			print >>of, "%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (data[s,0], data[s,1], data[s,2], data[s,3], data[s,4], data[s,5], data[s,6], data[s,7], data[s,8], data[s,9], data[s,10], data[s,11], data[s,12], data[s,13], data[s,14], data[s,15], data[s,16], data[s,17], data[s,18], data[s,19], data[s,20], data[s,21], data[s,22], data[s,23])
		of.close()




def parsec(path, outpath):
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
		print "Printing isochrone for [Fe/H] = %5.2f to '%s'" % (u, outname)
		# Loop through all webfiles for this [Fe/H] and read in data
		data = np.zeros([thisnlines, 24])
		for f in thisuni:
			df = open(webfiles[f], 'r')
			lines = df.read().splitlines()
			df.close()
			# Determine what file type this is
			if lines[12].find('Ks') >= 0:
				print "    Reading SDSS+JHK file '%s'" % (webfiles[f])
				ftype = 2
			elif lines[12].find('V') >= 0:
				print "    Reading UBVRI file '%s'" % (webfiles[f])
				ftype = 1
			else:
				print "    Reading IRAC file '%s'" % (webfiles[f])
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
			print >>of, "%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (data[s,0], data[s,1], data[s,2], data[s,3], data[s,4], data[s,5], data[s,6], data[s,7], data[s,8], data[s,9], data[s,10], data[s,11], data[s,12], data[s,13], data[s,14], data[s,15], data[s,16], data[s,17], data[s,18], data[s,19], data[s,20], data[s,21], data[s,22], data[s,23])
		of.close()




def dartmouth(path, outpath):
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
		print "Printing isochrone for [Fe/H] = %5.2f to '%s'" % (u, outname)
		# Loop through all webfiles for this [Fe/H] and read in data
		data = np.zeros([thisnlines, 24])
		for f in thisuni:
			df = open(webfiles[f], 'r')
			lines = df.read().splitlines()
			df.close()
			# Determine what file type this is
			if lines[5].find('SDSS') >= 0:
				print "    Reading SDSS file '%s'" % (webfiles[f])
				ftype = 1
			elif lines[5].find('Bessel') >= 0:
				print "    Reading UBVRI+JHK file '%s'" % (webfiles[f])
				ftype = 2
			else:
				print "    Reading IRAC file '%s'" % (webfiles[f])
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
			print >>of, "%6.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (data[s,0], data[s,1], data[s,2], data[s,3], data[s,4], data[s,5], data[s,6], data[s,7], data[s,8], data[s,9], data[s,10], data[s,11], data[s,12], data[s,13], data[s,14], data[s,15], data[s,16], data[s,17], data[s,18], data[s,19], data[s,20], data[s,21], data[s,22], data[s,23])
		of.close()
		



def makeiso(isopath, outpath):
	tmp = [x for x in subprocess.check_output("ls "+isopath+"*", shell=True).splitlines() if x.find('.dat') >= 0]
	if len(tmp) == 0:
		print "\n!!! Dartmouth Isochrones Detected.\n"
		dartmouth(isopath, outpath)
	else:
		testfile = tmp[0]
		df = open(testfile, 'r')
		lines = df.read().splitlines()
		df.close()
		if lines[1].find('Marigo') >= 0:
			print "\n!!! Padova Isochrones Detected.\n"
			padova(isopath, outpath)
		elif lines[1].find('PARSEC') >= 0:
			print "\n!!! PARSEC Isochrones Detected.\n"
			parsec(isopath, outpath)
		