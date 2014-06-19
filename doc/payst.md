PAYST
-----
Program handling matching photometry and membership data files togther into a master catalog.

### Input
*****

Individual data files must have a specific format to be read by `payst`: 

    [ID Columns] [RA] [DEC] [Magnitude + Uncertainty for each Filter]
    STAR001 08521310+1201377   133.054622000 12.027164000   14.402 0.029 13.725 0.031 13.626 0.034
    
In addition to photometry files, `payst` can handle matching membership data (radial velocity and proper motion) files as well, which should have a similar format to photometry files, with a few additional columns:

    RV:  [ID Columns] [RA] [DEC] [RV Magnitude + Uncertainty] [Membership %] [Qualitative Membership]
         16039  133.11575   11.99861   6.70   0.090  99  BLN
    PM:  [ID Columns] [RA] [DEC] [PM RA Magnitude + Uncertainty] [PM DEC Magnitude + Uncertainty] [Membership %] [Qualitative Membership]
         6   132.773239    11.578593  25.23   6.43 -29.33  28.44    -1   U
         
Qualitative membership designations are: `U` --- Unknown, `SM/SN` --- Single Member / Non-Member, `BM/BN` --- Binary Member / Non-Member, `BLM/BLN` --- Likely Binary Member / Non-Member. 

`payst` expects an option file with a line for each datafile to be matched, with the format:

    [File Name] [Number of ID Columns] [List of Filters] 
    NGC2682.2MASS.txt    2    J H K
    
Here, the 2nd element of each line corresponds to the number of columns in the datafile before RA + DEC. Acceptable filter characters are: 

* `U` `B` `V` `R` `I` for Johnson-Cousins *UBVRI* filters
* `SU` `SG` `SR` `SI` `SZ` for Thuan-Gunn SDSS *ugriz* filters
* `J` `H` `K` for 2MASS *JHK*<sub>S</sub> filters
* `B1` `B2` `B3` `B4` for IRAC \[3.6\] \[4.5\] \[5.8\] \[8.0\] filters
* `B1` `B2` `B5` `B6` for WISE \[3.4\] \[4.6\] \[12.0\] \[22.0\] filters.

If matching a PM or RV file, the list of filters is expected to be empty. `payst` determined which type of membership file is present by searching for `RV` or `PM` strings in the filename, so make sure to name membership data files accordingly.

### Execution
*****

`payst` can be run in two different modes: matching, and trimming. To match datafiles, you must specify two command line arguments to the `payst.py` program:

    python payst.py [option filename] [matching radius in arcsec]
    
The matching radius is used to determine how close together sources from different datafiles must be to be considered a match. By default, `payst` ignores any magnitude measurements with associated uncertainties greater than 0.1. To include more uncertain data, you can pass a third argument to payst which gives the upper-limit on the uncertainty which will be included in the matched file.

`payst`'s trimming mode can be started by specifying only a single command line argument when invoking the program:

    python payst [merged filename]
    
Here, `payst` is pointed to the previously-printed `*.Merged.txt` file, and there are several options for trimming the file down to a more useful size. Trimming to a specific radius around the cluster, or only stars with the necessary photometry, will speed up the SED-fitting process in the future.

### Output
*****

`payst` output files (from both modes) have the format:

| Column(s) | Data | Empty Value |
|:---------:|:-----|:------------|
| 1    | 2MASS ID | Made-Up ID `ID00000000+0000000` |
| 2    | RA  | |
| 3    | Dec | |
| 4-13  | UBVRI Magnitudes + Uncertainties | `99.999 9.999` |
| 14-23  | ugriz Magnitudes + Uncertainties | `99.999 9.999` |
| 24-29  | JHK<sub>S</sub> Magnitude + Uncertainty | `99.999 9.999` |
| 30-37  | [3.6] \[4.5\] [5.8] \[8.0\]<sub>IRAC</sub> Magnitudes + Uncertainties | `99.999 9.999` |
| 38-41  | [12.0] \[22.0\]<sub>WISE</sub> Magnitude + Uncertainty | `99.999 9.999` |
| 42-43  | RV Value + Uncertainty | `-9999 9.999` |
| 44-45  | PM RA Value + Uncertainty | `-9999 9.999` |
| 46-47  | PM Dec Value + Uncertainty | `-9999 9.999` |
| 48     | Membership Percentage | `-1` |
| 49     | Qualitative Membership Designation | `U` |
