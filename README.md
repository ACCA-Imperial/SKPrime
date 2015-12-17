# README #

SKPrime: A numerical implementation of the Schottky-Klein prime function in MATLAB.

This software is in the pre-release stage. Help documentaiton for MATLAB is a work in progress. Theory documentation is in-process, to be submitted very soon. Any quesions should be directed to Everett.

### How do I get set up? ###

* Download a zip of the master branch of the repository.
* Unpack the zip somewhere.
* Run the `install.m` function from MATLAB in the directory you've just unzipped.
    * This will put the SKPrime directory on your MATLAB search path.
* See the `example.m` file for usage examples.
* See the `example_accuracy.m` file for an accuracy test example.
* See the `example_slitmap.m` file for a slitmap with the prime function.

### FMM acceleration ###

If you have the [FMM2D software](http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html) installed, and the `zfmm2dpart.m` file is in your MATLAB search path, then the SKPrime software will take advantage of this for evaluation of points in the domain.

### Contribution guidelines ###

* Ask Everett.
* TBD

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
