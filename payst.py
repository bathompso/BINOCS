# The payst routine is completely contained within the binocs module
# This program is merely a wrapper that can be called directly from the command line
import sys, binocs

# Catalog matching mode
if len(sys.argv) >= 3: binocs.paystmatch(sys.argv[1], float(sys.argv[2])/3600.0)
		
# Catalog trimming mode
else: binocs.paysttrim(sys.argv[1])
		
		