# The makeiso routine is completely contained within the binocs module
# This program is merely a wrapper that can be called directly from the command line
from __future__ import print_function, division
import sys, binocs

# Print out a warning if the specified command line arguments aren't specified
if len(sys.argv) < 3:
	print("makeiso requires 2 arguments:")
	print("-- the path to the folder containing the downloaded isochrones")
	print("-- a folder to output the formatted isochrones")
	
# If necessary paths are specified, run makeiso
else: binocs.makeiso(sys.argv[1], sys.argv[2])

