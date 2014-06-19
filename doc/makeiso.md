MakeISO
-------
Program handling the formatting of isochrone files.

### Input
*****

The `makeiso` program takes web-downloaded isochrone tracks and converts them into a usable format for the SED-fitting routine. These downloaded isochrones can come from any of three sources:

* **Dartmouth**: The Dotter et al. (2008) isochrones can be downloaded from [here](http://stellar.dartmouth.edu/%7Emodels/isolf_new.html). For each desired metallicity, files for `UBV(RI)c + 2MASS + Kepler`, `SDSS ugriz`, and `Spitzer-IRAC` must be downloaded and placed in the same folder. Currently, there is no logic in place for files with varying [a/Fe] or Helium content.

* **Padova**: Marigo et al. (2008) isochrones can be downloaded from [here](http://stev.oapd.inaf.it/cgi-bin/cmd) (make sure to select Marigo et al. (2008) from the "previous set" section at the top). Select the second option in the "Ages/metallicities" section, to generate a sequence of isochrones. Change the Z value in this section to specify the isochrone's metallicity, but leave all other values on the page the defaults. Then, download the `UBVRIJHK`, `SDSS ugriz + 2MASS JHKs`, and `Spitzer IRAC+MIPS` isochrones files (from the "Photometric system" section) and place them in the same folder.

* **PARSEC**: Bressan et al. (2012) isochrones can be downloaded from the [same page as Padova](http://stev.oapd.inaf.it/cgi-bin/cmd). Download steps are the same as for Padova, only with keeping the "PARSEC version 1.1" evolutionary track selected.

### Execution
*****

The `makeiso` program can be run from the command line via:

    python makeiso.py [input folder] [output folder]
    
The input folder is the folder where all the downloaded isochrones were placed. `makeiso` assumes that all files within the specific folder are from the same isochrones system, so make sure to save Dartmouth, PARSEC, and Padova isochrones in separate folders.

### Output
*****

`makeiso` will search through all files in the input folder and output formatted files for all detected metallicities. Formatted isochrone files will have the filename:

    iso_[feh].[iso source].syn.dat
    [+0.25] Dartmouth isochrone: iso_p025.dm.syn.dat
    [-0.12] PARSEC isochrones:   iso_m012.pc.syn.dat
    
All isochrone files will have the same format:

| Column | Data |
|:------:|:-----|
| 1  | Log(Age) |
| 2  | Initial Mass (M<sub>sun</sub>) |
| 3  | Current Mass (M<sub>sun</sub>) |
| 4  | Log(Luminosity / L<sub>sun</sub>) |
| 5  | Log(T<sub>eff</sub>) |
| 6  | Log(g) |
| 7  | Bolometric Magnitude |
| 8-24| UBVRIugrizJHK<sub>S</sub>[3.6] \[4.5\] [5.8] \[8.0\] Magnitudes |