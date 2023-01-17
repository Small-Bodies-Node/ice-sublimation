# fastrot: Python version

## Notes

This version of fastrot will force the vapor pressure of CO and CO2 to 0 for very low temperatures (see code for limits) to avoid extrapolating polynomials well beyond their valid ranges.

## Attribution

Please cite Cowan & A'Hearn (1979, M&P 21, 155; [NASA ADS 1979M&P....21..155C](https://ui.adsabs.harvard.edu/abs/1979M%26P....21..155C); [DOI 10.1007/BF00897085](https://dx.doi.org/10.1007/BF00897085)) for the concepts behind these tools.

If you find the Python code useful in your projects, also cite:

    M. Van Selous and M. S. P. Kelley. 2021.  Small-Bodies-Node/ice-sublimation. GitHub, Commit: ***).

Replace the year as appropriate for the version you are using.  Also replace *** with the first ~7 characters of the git commit hash.  This can be shown in GitHub on the [commits page](https://github.com/Small-Bodies-Node/ice-sublimation/commits/main).  The hash prefix to use is the unique string along the right hand side (e.g., c338c13).

## Requirements

This script does not have any dependencies outside of the Python 3 standard library.

## Usage

```
fastrot.py [-h] --Av visual albedo --Air infrared albedo --rh
 heliocentric_distance --obl obliquity [--temp temperature]
        [--verbosity verbosity]
        species

example:
python fastrot.py 'H2O' --Av=0.00 --Air=0.05 --rh=3.98 --obl=90
```

1. species - string (or integer for backwards compatibility) specifying the desired ice: (1) H2O, (2) H2O-CH4 clathrate, (3) CO2, (4) CO.
2. Av - visual albedo
3. Air - infrared albedo
4. rh - heliocentric distance in au
5. obliquity - angle between the object's rotational axis and its orbital axis

## survey_fastrot.py

This script can be used to call `fastrot.py` over a parameter space with a single call.

### Requirements

`survey_fastrot.py` does not introduce any additional packages outside the standard library.
The only added requirement is that both `survey_fastrot.py` and `fastrot.py` are kept within the same directory.

### Usage

```
usage: survey_fastrot.py [-h] --species_list species [species ...] --Av_list visual albedo [visual albedo ...] --Air_list infrared albedo [infrared albedo ...]
                         --rh_list heliocentric_distance [heliocentric_distance ...] --obl_list obliquity [obliquity ...]
                        
example:
python survey_fastrot.py --species_list 'H2O' 'H2O-CH4' --Av_list 0.05 0.10 --Air_list 0 --rh_list 3.98 --obl_list 40 41 42
```

1. species_list species [species ...]
                        Desired ice species to be considered.
                        The valid inputs are: H2O, H2O-CH4, CO2, CO
2. Av_list - visual albedo [visual albedo ...]
3. Air_list - infrared albedo [infrared albedo ...]
4. rh_list - heliocentric_distance [heliocentric_distance ...]
5. obl_list - obliquity [obliquity ...]

## Tests

The script `test_fastrot.py` will compare `fastrot.py` output to a set of precomputed values from the FORTRAN code for a pole-on case and the same number of latitude steps.  All tests agree within 0.06%.

## Data

The [data](data/) directory provides the pole-on case, which is identical both to the non-rotating case and to the case of zero thermal inertia. The visual Bond albedo is 5% and the thermal emissivity is 100%. The parameters were chosen to provide a direct comparison between the updated model and the example output from the original FORTRAN code.

## Notice

Copyright 2021,2023 The Ice Sublimation Tool Developers.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
