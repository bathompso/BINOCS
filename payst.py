# The payst routine is completely contained within the binocs module
# This program is merely a wrapper that can be called directly from the command line
from __future__ import print_function
import sys, binocs

# Catalog matching mode
if len(sys.argv) >= 3: binocs.paystmatch(sys.argv[1], float(sys.argv[2])/3600.0)
		
# Catalog trimming mode
elif len(sys.argv) >= 2: binocs.paysttrim(sys.argv[1])

# No arguments provided, print documentation
else:
	print("==== PAYST ====")
	print("Run mode 1 [MERGE]:  payst [name of option file] [matching radius in arcsec]")
	print("Run mode 2 [TRIM]:   payst [name of merged data file]")
	print("\nOutput File Format:")
	print("1:     2MASS ID, or generated ID if no 2MASS data")
	print("2:     RA")
	print("3:     Dec")
	print("4-5:   U Magnitude + Uncertainty")
	print("6-7:   B Magnitude + Uncertainty")
	print("8-9:   V Magnitude + Uncertainty")
	print("10-11: R Magnitude + Uncertainty")
	print("12-13: I Magnitude + Uncertainty")
	print("14-15: u' Magnitude + Uncertainty")
	print("16-17: g' Magnitude + Uncertainty")
	print("18-19: r' Magnitude + Uncertainty")
	print("20-21: i' Magnitude + Uncertainty")
	print("22-23: z' Magnitude + Uncertainty")
	print("24-25: J Magnitude + Uncertainty")
	print("26-27: H Magnitude + Uncertainty")
	print("28-29: K Magnitude + Uncertainty")
	print("30-31: [3.6] Magnitude + Uncertainty")
	print("32-33: [4.5] Magnitude + Uncertainty")
	print("34-35: [5.8] Magnitude + Uncertainty")
	print("36-37: [8.0] Magnitude + Uncertainty")
	print("38-39: [12.0] Magnitude + Uncertainty")
	print("40-41: [22.0] Magnitude + Uncertainty")
	print("42-43: RV Value + Uncertainty")
	print("44-45: PM RA Value + Uncertainty")
	print("46-47: PM Dec Value + Uncertainty")
	print("48:    Membership Percentage")
	print("49:    Qualitative Membership Probability")
	print('')
		
		