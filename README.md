BINOCS (Binary INformation from Open Clusters using SEDs)
=====

BINOCS is a suite of python codes that facilitate binary star detection in open clusters. This code, and the method it implements, is described in Thompson et al. (*in prep*)



Installation
------------
BINOCS code will run under Python versions 2.7.x, as well as Python 3.x

In addition to the default Python modules, BINOCS requires the following modules be installed: 

* NumPy
* SciPy
* MatPlotLib
* PyOpenCL

All necessary modules can be installed via `pip`. Once pre-requisites are met, simply add the `binocs/python` folder to your `$PYTHONPATH` variable.



Available Routines
------

The BINOCS code consists of three separate routines, each of which is implemented by a program in the BINOCS root folder:

* `payst` --- handles matching of photometry data files, as well as membership data, into a formatted master catalog
* `makeiso` --- processes downloaded isochrone files into a readable format 
* `binaryfit` --- executes the BINOCS binary detection technique

Descriptions of input, execution, and output of each of the routines is located in files in the `doc` folder.



License
-------

Copyright 2014 Ben Thompson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

If you publish results using this software, please reference Thompson et al. (*in prep*).

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

