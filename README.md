binaryfit
=========

binaryfit is a suite of python codes that facilitate binary star detection in open clusters. A paper detailing the inner workings of the main binaryfit routine, as well as initial results for several open clusters is in prep.

Dependencies
------------

Running the included python programs will require Python 2.7 or above, as well as numpy.

The binaryfit python routines are parallelized in order to run in a reasonable amount of time. Before running any of these routines, please make sure to install the python modules [ParallelPython](http://parallelpython.com) and [PyOpenCL](http://mathema.tician.de/software/pyopencl/)

Component Codes
---------------

### PAYST

PAYST (Recursive Acronym: Payst All Your Stuff Together) is a python code that handles the merging and trimming of individual photometry datasets. The binaryfit routine requires photometry over a large wavelength range (visual through mid-IR), but often studies only contain a specific region. PAYST matches stars between datasets by coordinates and outputs a matched dataset file.

##### Dataset Format

PAYST expects dataset files to be in a specific format.

- The first columns can be any number of star identifiers. These can be 2MASS designations, or other star numbers from the reference you are pulling the data from. Any number of identification columns can appear at the beginning of each line.

- After identifiers must be two columns of RA and Dec in decimal degrees.

- After RA / Dec must be all magnitudes. Each filter will have two columns: a magnitude and an uncertainty. For non-existant magnitudes, enter 99.999.

PAYST can also handle Proper Motion or Radial Velocity catalogs. These will have a slightly different format than photometry datafiles.

- Same as photometry, the first columns can be any number of identifiers.

- RA / Dec in decimal degrees.

- RV or Proper Motion values and uncertainties. If proper motion, pm\_RA is expected before pm\_Dec.

- Membership Probability (in %). Many studies specify cluster membership probability percentage for each star depending on its velocity. This column is expected to be an integer. If no known membership probability, enter -1 in this column.

- Qualitative Membership Flag. In RV studies, some stars are flagged as combinations of binary / single and member / likely-member / non-member. Flags include: BM, BLM, BN, SM, SLM, SN. If no membership flag is specified, use U in this column.

##### PAYST Option File

PAYST reads in an input option file that has lines of a similar format: "NGC188.2MASS.txt 1 J H K". Necessary columns are:

1. Dataset file name (relative path with respect to option file)

	- PAYST denotes stars with their 2MASS identification '00000000+0000000'. For PAYST to detect these, the 2MASS data file must have the string '2MASS' in its file name.
	
2. Number of identifier columns. This is simply the number of columns before RA / Dec in the data file. This lets PAYST know how many columns to skip when searching for magnitude information.

3. List of filters. Space-separated list of filters for which magnitudes are specified in the datafile. PAYST currently supports the following filters:

	- Johnson-Cousins Visual Filters: U B V R I
	
	- Thuan-Gunn SDSS Filters: SU SG SR SI SZ
	
	- 2MASS Filters: J H K
	
	- Spitzer IRAC Filters (\[3.6\]\[4.5\]\[5.8\]\[8.0\]): B1 B2 B3 B4
	
	- WISE Filters (\[3.4\]\[4.6\]\[12.0\]\[24.0\]): B1 B2 B5 B6  (note that B1 and B2 are shared with IRAC, as magnitudes are very similar)
	
	Make sure to use these filter names in the option file, as they are the specific strings PAYST is looking for.
	
4. Proper Motion and Radial Velocity Catalogs

	To specify a velocity catalog instead of photometry, simply leave the filter list blank. This will produse a line that only has the data file name, and number of identifier columns.
	
	PAYST determines what type of velocity catalog it has by looking at the file name. Make sure to have the string 'RV' or 'PM' in the file name of your velocity catalogs, depending on which type of data you have. If PAYST cannot find these strings in the file name, it will ask for input when running.

