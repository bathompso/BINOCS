binaryfit
---------

Program executing the `binocs` SED-fitting procedure.

### Input
*****

Like `payst`, `binaryfit` expects an option file to be provided which specifies the files and parameters to be used in the SED fitting. This option file must have a line for each of the important parameters of the form:

    [parameter name] = [value]
    
| Parameter Name | Description |
|:--------------:|:------------|
| `data` | Path to `payst` output data file |
| `iso`  | Path to `makeiso` formatted isochrone file to be used |
| `dm`   | Isochrone mass interpolation grid spacing (in M<sub>sun</sub>), usually 0.01 |
| `age`  | Log(age) of isochrone to be used |
| `m-M`  | Distance modulus of cluster |
| `ebv`  | E(B-V) reddening of cluster |
| `nruns` | Number of fitting iterations, usually 300 |
| `fid`  | *(optional)* Path to fiducial ridgeline file |

### Execution
*****

To run the SED-fitting procedure:

    python binaryfit.py [option file]
    
### Output
*****

`binaryfit` will output several files, each within the same folder as the option file, and containing the filename of the datafile specified for the run. These output files are:

| Suffix | Description |
|:------:|:-----|
| `.fid.dat` | Fiducial-ridgeline adjusted isochrone. Same format as original `makeiso` isochrone. Only printed if `fid` parameter is located in option file. |
| `.bin` | Synthetic SED library containing all single and binary models to be compared to. Columns 1-2 are primary and secondary model masses. |
| `--binary.txt` | Initial results. Single / Binary cut-off set at q = 0.3 |
