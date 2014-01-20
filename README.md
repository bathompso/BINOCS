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

PAYST (Recursive Acronym: Payst All Your Stuff Together) handles the merging and trimming of individual photometry datasets. The binaryfit routine requires photometry over a large wavelength range (visual through mid-IR), but often studies only contain a specific region. PAYST matches stars between datasets by coordinates and outputs a matched dataset file.

### MakeIso

MakeIso reads in base isochrone files and creates a full isochrone catalog that can be read by the main binaryfit program. Currently, makeiso can handle reading inputs from the Dartmouth, Padova, PARSEC and Y<sup>2</sup> isochrone systems.

### BinaryFit

BinaryFit compares a full cluster dataset (created through PAYST) and an isochrone model (created through MakeIso) to detect binary systems within the cluster. The output file will include determinations of binary, single and non-member stars, as well as mass determinations for member binaries and singles.


Specific details on how to run each component code is located within the individual code's readme file.


License
-------

Copyright 2014 Ben Thompson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

If you publish results using this software, please reference Thompson et al. (2014, in prep).

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

