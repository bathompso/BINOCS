PAYST
-----

### Dataset Format

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

### PAYST Option File

PAYST reads in an input option file that has lines of a similar format: "NGC188.2MASS.txt 1 J H K". Necessary columns are:

1. Dataset file name (relative path with respect to option file). PAYST denotes stars with their 2MASS identification '00000000+0000000'. For PAYST to detect these, the 2MASS data file must have the string '2MASS' in its file name.
	
2. Number of identifier columns. This is simply the number of columns before RA / Dec in the data file. This lets PAYST know how many columns to skip when searching for magnitude information.

3. List of filters. Space-separated list of filters for which magnitudes are specified in the datafile. PAYST currently supports the following filters:

	- Johnson-Cousins Visual Filters: U B V R I
	
	- Thuan-Gunn SDSS Filters: SU SG SR SI SZ
	
	- 2MASS Filters: J H K
	
	- Spitzer IRAC Filters (\[3.6\]\[4.5\]\[5.8\]\[8.0\]): B1 B2 B3 B4
	
	- WISE Filters (\[3.4\]\[4.6\]\[12.0\]\[24.0\]): B1 B2 B5 B6  (note that B1 and B2 are shared with IRAC, as filters are very similar)
	
4. Proper Motion and Radial Velocity Catalogs

	To specify a velocity catalog instead of photometry, simply leave the filter list blank. This will produse a line that only has the data file name, and number of identifier columns.
	
	PAYST determines what type of velocity catalog it has by looking at the file name. Make sure to have the string 'RV' or 'PM' in the file name of your velocity catalogs, depending on which type of data you have. If PAYST cannot find these strings in the file name, it will ask for input when running.
	
### Running PAYST

To run PAYST, simply run the command: `python payst2.py   /path/to/option/file   matching-radius`. Matching-radius is the maximum distance between matched sources, specified in arcseconds. A reasonable value is &le; 1.0.

### PAYST Output

PAYST will output a file named *.Merged.txt, based on the naming of your option file. The Merged dataset will have the following columns:

	1: 2MASS ID (if not matched to 2MASS, a desgination will be created with prefix 'BT')
	2: RA
	3: Dec
	4-5: U Magnitude + Uncertainty
	6-7: B Magnitude + Uncertainty
	8-9: V Magnitude + Uncertainty
	10-11: R Magnitude + Uncertainty
	12-13: I Magnitude + Uncertainty
	14-15: u' Magnitude + Uncertainty
	16-17: g' Magnitude + Uncertainty
	18-19: r' Magnitude + Uncertainty
	20-21: i' Magnitude + Uncertainty
	22-23: z' Magnitude + Uncertainty
	24-25: J Magnitude + Uncertainty
	26-27: H Magnitude + Uncertainty
	28-29: K Magnitude + Uncertainty
	30-31: [3.0] Magnitude + Uncertainty
	32-33: [4.5] Magnitude + Uncertainty
	34-35: [5.8] Magnitude + Uncertainty
	36-37: [8.0] Magnitude + Uncertainty
	38-39: [12.0] Magnitude + Uncertainty
	40-41: [22.0] Magnitude + Uncertainty
	42-43: RV Value + Uncertainty
	44-45: PM RA Value + Uncertainty
	46-47: PM Dec Value + Uncertainty
	48: Membership Percentage
	49: Qualitative Membership Probability
	
### PAYST Databse Trimming

While PAYST's main mode is to match separate datasets into a master catalog, PAYST is also capable of handling trimming the full dataset down to selections that are more interesting. To run PAYST in this mode, simply invoke the command: `python payst2.py /path/to/Merged.txt`, poiting PAYST to the output file from the previous run. After loading the database and asking for an initial CMD, PAYST will present you with several trimming options:

- _RA / Dec Trim_: Often, photometry datasets will cover much larger areas than the spatial extent of the cluster. To remove obvious background contamination, it helps to trim the dataset to a reasonable radius around the cluster. When choosing this option, PAYST will ask for the cluster's central RA and Dec, as well as a trimming radius in arcminutes. PAYST will display a spatial plot of all detected sources, along with your specified radius cut in order to confirm that you wish to trim out all outside stars.

- _Membership Trim_: If your dataset has RV or PM data included (along with membership percentages or qualitative flags), you can trim out only the stars you are interested in. This cut has two parts: one based on membership percentage and one on the qualitative flags. If you have a large amount of contamination in your CMD and good velocity catalogs, it may be helpful to select out only confirmed cluster members.

- _Reddening Trim_: To remove foreground or background contamination, you can select out cluster members based on reddening. Using the RJCE method [(Majewski et al. 2011)](http://adsabs.harvard.edu/abs/2011ApJ...739...25M), individual star reddenings can be determined from the _H_-[4.5] color. In the resulting plot, look for an overdensity in the digram and specify minimum and maximum _H_-[4.5] values to trim out obvious cluster non-members.

- _Full Photometry Cut_: For the main binaryfit routine to be effective, a star must contain 3 visual magnitudes, 3 near-IR magnitudes, and 2 mid-IR magnitudes. Any other stars will not be considered during the binaryfit routine. To select out only stars with sufficient photometry, choose this option. PAYST will list the number of stars with different numbers of valid magnitudes, as well as the number of stars with SED-fitting-worthy data. Selecting the '332' option will select only stars that are acceptable for future steps.

After accurately trimming your dataset, you can print out only those stars into a file named *.Trimmed.txt (it will have the same primary name as your *.Merged.txt file). The trimmed file will have the same format as the original merged file.

