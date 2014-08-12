BINOCS (Binary INformation from Open Clusters using SEDs)
=====

BINOCS is a suite of python codes that facilitate binary star detection in open clusters. This code, and the method it implements, is described in Thompson et al. (*in prep*)



Installation
------------
BINOCS code will run under Python versions 2.7.x, as well as Python 3.x

In addition to the default Python modules, BINOCS requires the following modules be installed: 

* NumPy
* SciPy
* PyOpenCL

All necessary modules can be installed via `pip`. Once pre-requisites are met, simply add the `binocs/python` folder to your `$PYTHONPATH` variable.



Available Routines
------

The BINOCS code consists of three separate routines, each of which is implemented by a program in the BINOCS root folder:

* `payst` --- handles matching of photometry data files, as well as membership data, into a formatted master catalog
* `makeiso` --- processes downloaded isochrone files into a readable format 
* `binaryfit` --- executes the BINOCS binary detection technique

Descriptions of input, execution, and output of each of the routines is located in files in the `doc` folder.

