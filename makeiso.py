# The makeiso subroutines are contained within the binocs module
# This program is merely a wrapper that can be called directly from the command line
from __future__ import print_function, division
import sys, binocs, subprocess

# Print out a warning if the specified command line arguments aren't specified
if len(sys.argv) < 3:
	print("makeiso requires 2 arguments:")
	print("-- the path to the folder containing the downloaded isochrones")
	print("-- a folder to output the formatted isochrones")
	
# If necessary paths are specified, run makeiso
else: 
	isopath, outpath = sys.argv[1], sys.argv[2]
	tmp = [x for x in subprocess.check_output("ls "+isopath+"*", shell=True).splitlines() if x.find('.dat') >= 0]
	if len(tmp) == 0:
		print("\n!!! Dartmouth Isochrones Detected.\n")
		dartmouth(isopath, outpath)
	else:
		testfile = tmp[0]
		df = open(testfile, 'r')
		lines = df.read().splitlines()
		df.close()
		if lines[1].find('Marigo') >= 0:
			print("\n!!! Padova Isochrones Detected.\n")
			padova(isopath, outpath)
		elif lines[1].find('PARSEC') >= 0:
			print("\n!!! PARSEC Isochrones Detected.\n")
			parsec(isopath, outpath)
