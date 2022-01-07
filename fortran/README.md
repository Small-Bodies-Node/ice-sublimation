# cgifastrot.f


## Attribution

Please cite Cowan & A'Hearn (1979, M&P 21, 155; [NASA ADS 1979M&P....21..155C](https://ui.adsabs.harvard.edu/abs/1979M%26P....21..155C); [DOI 10.1007/BF00897085](https://dx.doi.org/10.1007/BF00897085)) for the concepts behind these tools.

If you find the FORTRAN code useful in your projects, also cite:

    J. J. Cowan, M. F. A'Hearn, and B. Prager. 2021.  Small-Bodies-Node/ice-sublimation. GitHub, Commit: ***.

Replace the year as appropriate for the version you are using.  Also replace *** with the first ~7 characters of the git commit hash.  This can be shown in GitHub on the [commits page](https://github.com/Small-Bodies-Node/ice-sublimation/commits/main).  The hash prefix to use is the unique string along the right hand side (e.g., c338c13).

## Usage

The rapidly rotating nucleus model is available as `cgifastrot.f`.  This code was used by the PDS Small Bodies Node for an online tool.  It requires `sublime.f`.  To compile:
```
gfortran cgifastrot.f sublime.f -o fastrot
```

Usage:
```
fastrot <species> <Av> <Air> <rh> <obliquity>
```
1. species - integer specifying the desired ice: (1) H2O, (2) H2O+CH4 clathrate, (3) CO2, (4) CO.
2. Av - visual albedo
3. Air - infrared albedo
4. rh - heliocentric distance in au
5. obliquity - 90 - angle between rotation axis and the solar direction

## Results

Tabulated results are available for the `fastrot` model and pole-on case, which is identical both to the non-rotating case and to the case of zero thermal inertia. The visual Bond albedo is 5% and the thermal emissivity is 100%. To calculate average sublimation for other obliquities (angles between equatorial plane and the comet-sun line), other values of the albedo and emissivity, or other heliocentric distances, this provided code may be used. The example output is given in the [data](data/) directory.

## Notice

Copyright 2010 The Ice Sublimation Tool Developers.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
