# README #

SKPrime: A numerical implementation of the Schottky-Klein prime function in MATLAB.

This software is in the pre-release stage. Comprehensive help documentaiton for MATLAB is a work in progress. The theory behind the numerics is given in the paper listed in the citation section below.

### Citing the software ###

When citing this software please mention the URL of the master repository (https://github.com/ehkropf/SKPrime), and the paper

> D.G. Crowdy, C.C. Green, E.H. Kropf, M.M.S. Nasser. "The Schottky-Klein prime function: a theoretical and computational tool for applications." IMA Journal of Applied Mathematics, 2016, doi: 10.1093/imamat/hxw028.

### How do I get set up? ###

* Download a zip (or tar.gz) of the [current release](https://github.com/ACCA-Imperial/SKPrime/releases/tag/v0.1.5).
* Unpack the zip somewhere.
* Run the `install.m` function from MATLAB in the directory you've just unzipped.
    * This will put the SKPrime directory on your MATLAB search path.
* See the `example.m` file for usage examples.
* See the `example_accuracy.m` file for an accuracy test example.
* See the `example_slitmap.m` file for a slitmap with the prime function.

### FMM acceleration ###

If you have the [FMM2D software](http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html) installed, and the `zfmm2dpart.m` file is in your MATLAB search path, then the SKPrime software will take advantage of this for evaluation of large numbers of points in the domain.

### Contribution guidelines ###

* Ask Everett.
* Please see the [wiki](https://github.com/ehkropf/SKPrime/wiki).

### Contact ###

* Everett Kropf <e.kropf@imperial.ac.uk>

### License ###

This file is part of SKPrime.

SKPrime is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SKPrime is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SKPrime.  If not, see <http://www.gnu.org/licenses/>.
