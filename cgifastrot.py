"""
    Description
    -----------
    This program calculates the average sublimation per unit area for
    a rapidly rotating cometary nucleus.  For a
    sufficiently rapid rotation, or equivalently for sufficiently high
    thermal inertia, a parallel of latitude is an isotherm and this
    is assumed by the program.

    PROGRAM ITERATES ENERGY BALANCE EQUATION BY NEWTON-RAPHSON METHOD
    TO GET EQUILIBRIUM TEMPERATURE.  SIMPSON'S RULE IS USED TO
    INTEGRATE OVER LATITUDE.

    The properties of the ice are handled in the separate subroutine
    sublime which provides the vapor pressure and latent heat, as well
    as their derivatives, for a given temperature.

    sb is the set of values of sin(b) at which the sublimation is calculated.
    b is latitude, frac is the effective value of cos(theta) for each latitude
    In this version, there are 181 steps in latitude (variable nb in the
    data statement).

   Modification History
   --------------------
   Mark Van Selous July 2021
   - Rewrote cgifastrot.f in python as cgifastrot.py
   - Increased the latitude step size from 41 to 181.
   - Thanks to this increased step size and having more accurate trig functions
   than were available in fortran77, there is very slight alteration in the final results.

    B. Prager 06/24
    - Original fastrot.f renamed to cgifastrot.f to make a distinction.
    - Removed file i/o from cgifastrot.f to avoid permissions issues with cgi scripts.
    - Added command line input of parameters. Parameter 1: Species. Parameter 2: Visual
    Albedo. Parameter 3: Infared Albedo. Parameter 4: Heliocentric Distance. Parameter 5:
    inclination.
    - Initialized the parameters such that they can accept nine characters of data.
"""

import math
from subprocess import *
import sys
import json
import logging

nb = 181
z = [0] * nb
frac = [0] * nb
q = 1

prec = 0
delsb = 2.0 / (nb - 1)

sb = [-1]
b = [math.asin(sb[0])]

for idx in range(1, nb):
    sb.append(sb[0] + idx * delsb)
    b.append(math.asin(sb[-1]))

speciesList = ['H2O', 'H2O_CH4', 'CO2', 'CO']

logging.basicConfig(level=logging.DEBUG)


def run_model(species, a_v, a_ir, rh, inc, temp=-1):
    """
    A call of this function replicates the behavior of the original cgifastrot.f script.
    After reading validating the input parameters, run_model() will iterate through mainloop().


    Parameters
    ----------

    species : str or int
        Desired ice species to be considered
        - 1: 'H2O'
        - 2: 'H20_CH4'
        - 3: 'CO2'
        - 4: 'CO'

    a_v : float (a_v > 0)
        Visual albedo

    a_ir : float
        Infrared albedo

    rh : float
        Heliocentric distance (in au)

    inc : float
        Obliquity - 90 - angle between rotation axis and the solar direction

    temp: float
        Initial temperature. If this parameter is not specified, a species
        dependent initial value will be used:
        - H2O: 190 K
        - H2O_CH4: 190 K
        - CO2: 100 K
        - CO: 60 K


    Returns
    -------

    species: str
        Inputted species

    inc: float
        Inputted obliquity

    rh: float
        Inputted heliocentric distance (in au)

    rlog: float
        rlog = log10(rh)

    a_v : float (a_v > 0)
        Inputted visual albedo

    a_ir : float
        Inputted infrared albedo

    zbar: float
        Average sublimation per unit area

    zlog: float
        zlog = log10(zbar)

    """

    if isinstance(species, int):
        if species in range(1, len(speciesList) + 1):
            # The species indexing if offset by +1 to make the inputs consistent with the original fortran script.
            species = speciesList[species - 1]
        else:
            logging.error(f"The inputted index of {species} is not currently supported \n"
                          f"Please input one of the following integers or strings: \n"
                          f"1: 'H2O', 2: 'H2O_CH4', 3: 'CO2', 4: 'CO'")
            sys.exit(1)
    if species not in speciesList:
        logging.error(f"The inputted species of {species} is not currently supported \n"
                      f"Please input one of the following integers or strings: \n"
                      f"1: 'H2O', 2: 'H2O_CH4', 3: 'CO2', 4: 'CO'")
        sys.exit(1)

    if a_v < 0:
        logging.error(
            f'A visual albedo of {a_v} is not a valid input.'
            ' Please input a value greater than 0.')
        sys.exit(1)

    logging.info("Input Parameters:")
    logging.info(
        f'Species = {species}, Avis = {a_v}, Air = {a_ir}, r_H = {rh}, Incl = {inc}')

    incl = (90 - inc) * math.pi / 180

    spec, mass, xlt, xltprim, press, pprim, temp = sublime(species, temp)
    root = 1 / math.sqrt(mass * 2 * math.pi * boltz)

    nflag = 1
    perc = 0
    gd = None
    for n in range(0, nb):
        temp, gd, perc, nflag = main_loop(
            n, species, a_v, a_ir, rh, incl, temp, root, nflag, perc, gd)

    zbar = 0.
    for nn in range(0, nb - 1):
        zbar = zbar + 0.5 * (z[nn] + z[nn + 1]) * delsb

    zbar = zbar / 2
    zlog = math.log10(zbar)
    rlog = math.log10(rh)
    logging.info("Providing Final Results")
    print(json.dumps({
        "Species": species,
        "Obliquity": inc,
        "r_H": rh,
        "rlog": "{:.3f}".format(rlog),
        "a_v": a_v,
        "a_ir": a_ir,
        "Zbar": "{:.3e}".format(zbar),
        "Zlog": "{:.3f}".format(zlog),
    }))

    return species, inc, rh, rlog, a_v, a_ir, zbar, zlog


def main_loop(n, species, a_v, a_ir, rh, incl, temp, root, nflag, perc, gd):
    """
    Parameters
    ----------
    n : int
        counter (ranging from 0 to nb-1) which tracks
    species: str
        Inputted species
    a_v : float (a_v > 0)
        Inputted visual albedo
    a_ir : float
        Inputted infrared albedo
    rh: float
        Inputted heliocentric distance (in au)
    incl : float
        Inputted obliquity expressed in radians.
    temp : float
        Temperature (in Kelvin). This is one of the three parameters which is updated by main_loop().
    root : float
        root is given by: root = 1 / sqrt(mass * 2 * pi * boltz)
        Below is an approximation for the four possible values:
        root(H2O) ~= 6.194e18
        root(H2O_CH4) ~= 6.194e18 (= root(H20))
        root(CO2) ~= 3.962e18
        root(CO) ~= 4.966e18
    nflag : int
        Counter for the total number of times main_loop() has been called.
    perc: float
        This is one of the three parameters which is updated by main_loop().
    gd : None / int
        Denotes the greatest index, n, at which calc_perc() was called.
        If gd = None, then calc_perc() has not been called yet.
        This is one of the three parameters which is updated by main_loop().

    Returns
    -------
    temp : float
        Temperature (in Kelvin). This is one of the three parameters which is updated by main_loop().
    gd : None / float
        gd is used in the case of nflag >= 10,000,000
    perc : float
    """
    root_t = math.sqrt(temp)
    if b[n] <= -incl:
        frac[n] = 0
        z[n] = 0
        return temp, gd, perc, nflag
    elif b[n] > incl:
        frac[n] = sb[n] * math.cos(incl)

    else:
        x1 = (
            math.cos(incl)
            * sb[n]
            * (math.acos(-math.tan(b[n]) * (1 / math.tan(incl))))
            / math.pi
        )
        x2 = (
            math.sin(incl)
            * math.cos(b[n])
            * math.sin(math.acos(-math.tan(b[n]) / math.tan(incl)))
            / math.pi
        )
        frac[n] = x1 + x2

    spec, mass, xlt, xltprim, press, pprim, temp = sublime(species, temp)
    sun = f0 * frac[n] * (1. - a_v) / rh ** 2
    radiat = (1 - a_ir) * sigma * temp ** 4
    evap = q * root / root_t * press * xlt
    phi = radiat + evap - sun
    z[n] = max(evap / xlt, 1e-30)

    drad = 4 * radiat / temp
    x1 = pprim * xlt
    x2 = press * xltprim

    devap = q * root / root_t * (x1 + x2)
    phipri = drad + devap

    dt = math.copysign(min(10, abs(phi / phipri / 2)), phi / phipri)
    tp = temp - dt
    temp = tp

    if abs(phi / sun) < 1e-4 or abs(phi) < 1e-4:
        gd, perc = calc_perc(n, rh, incl, perc)
        return temp, gd, perc, nflag
    if nflag >= 1000:
        subdis = z[gd] * 4 * math.pi * rh * 149.6e11
        if perc != 0:
            if (subdis / perc) < 1e-02:
                z[n] = 0
                gd, perc = calc_perc(n, rh, incl, perc)
                return temp, gd, perc, nflag
            else:
                logging.error('error calculating sublimation')
                sys.exit(1)

    nflag += 1
    temp, gd, perc, nflag = main_loop(
        n, species, a_v, a_ir, rh, incl, temp, root, nflag, perc, gd)

    return temp, gd, perc, nflag


def calc_perc(n, rh, incl, perc):
    logging.info('Incl: %f, Lat: %f', incl, b[n])
    gd = n
    perc = perc + z[n] * 4 * math.pi * rh * 149.6e11
    return gd, perc


############
a_v_in = 0.05
a_ir_in = 0.00
rh_in = 3.98
inc_in = 90

# speciesList = ['H2O', 'H2O_CH4', 'CO2', 'CO']
index_in = speciesList.index('H2O') + 1

run_model(index_in, a_v_in, a_ir_in, rh_in, inc_in, )
