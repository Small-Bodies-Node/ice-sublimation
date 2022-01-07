# ice-sublimation
Cometary ice sublimation model, based on Cowan & A'Hearn (1979, M&P 21, 155). An online version of this tool is at https://pds-smallbodies.astro.umd.edu/tools/ma-evap/index.shtml.  Two versions of the tool are available: the original FORTRAN version, and a Python version derived from the FORTRAN code.

## Attribution

Please cite Cowan & A'Hearn (1979, M&P 21, 155; [NASA ADS 1979M&P....21..155C](https://ui.adsabs.harvard.edu/abs/1979M%26P....21..155C); [DOI 10.1007/BF00897085](https://dx.doi.org/10.1007/BF00897085)) for the concepts behind these tools.

If you find the FORTRAN code useful in your projects, please cite:

    J. J. Cowan, M. F. A'Hearn, and B. Prager. 2021.  Small-Bodies-Node/ice-sublimation. GitHub, Commit: ***.

If you find the Python code useful in your projects:

    M. Van Selous and M. S. P. Kelley. 2021.  Small-Bodies-Node/ice-sublimation. GitHub, Commit: ***.

Replace the year as appropriate for the version you are using.  Also replace *** with the first ~7 characters of the git commit hash.  This can be shown in GitHub on the [commits page](https://github.com/Small-Bodies-Node/ice-sublimation/commits/main).  The hash prefix to use is the unique string along the right hand side (e.g., c338c13).

## Overview
This work, originally developed by J. Cowan and M. A'Hearn, provides a simple tool to calculate the sublimation of ices under various circumstances. The calculations are all based on the methods described by Cowan and A'Hearn (1979 Moon and Planets 21, 155-171), which in turn are based on earlier work by Delsemme and others as referenced in that paper. The calculations have not been substantially altered from the original publication and but updated vapor pressures and latent heats have been used. These changes to the input parameters lead to only small changes in the resultant sublimation rates. We note, as pointed out to us by W. Huebner and D. Boice, that the use of empirical results for heat of sublimation and vapor pressure lead the results to be inconsistent with the Clausius-Clapeyron equation, but this has negligible effect on the results for most (but not all) physical situations. Results are available for pure water, pure CO, pure CO2.

## fastrot
A rapidly rotating nucleus is one for which the thermal inertia is large enough that parallels of latitude become isotherms. The average sublimation over the nucleus is then a function of the obliquity, i.e., of the orientation of the rotation axis. The sublimation is relatively high if the nucleus is pole-on toward the sun (obliquity = 90°) and much smaller if the axis is perpendicular to the comet-sun line (obliquity = 0°). All calculations assume a spherical body. Note that a non-rotating comet is thus identical to a comet that is pole-on toward the sun. The approach given here calculates the sublimation separately at each latitude (and for a rapid rotator the sublimation is constant all the way around the parallel of latitude, even on the night side) and then calculates the appropriate average over the entire surface, i.e., the average over all 4πR^2 of the surface, including areas where the actual sublimation is zero.

The rapidly rotating nucleus model is available as [FORTRAN](fortran/cgifastrot.f) and [Python](python/fastrot.py) code. See the associated README files for instructions on usage.

