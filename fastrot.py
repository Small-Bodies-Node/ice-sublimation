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
   - Rewrote cgifastrot.f and sublime.f in python as fastrot.py
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

# Constants
sigma = 5.67e-5
f0 = 1.39e6
dynmm = 1.333e3
boltz = 1.38e-16
ergcal = 6.953e-17
proton = 1.67e-24

tstart = {'H2O': 190,
          'H2O_CH4': 190,
          'CO2': 100,
          'CO': 60,
          }


def sublime(species, temp):
    """
    Description
    -----------
    Calculates the latent heat of sublimation and the vapor pressure of
    the solid for various ices and the derivatives thereof

    Parameters
    ----------
    species : str
        Desired ice species to be considered. The valid inputs are:
        - 'H2O'
        - 'H20_CH4'
        - 'CO2'
        - 'CO'
    temp: float
        Temperature (in Kelvin) If temp <= 0 K, then the code defaults to a species
        dependent initial value:
        - H2O: 190 K
        - H2O_CH4: 190 K
        - CO2: 100 K
        - CO: 60 K

    Returns
    -------
    species : str
    mass : float
    xlt : float
    xltprim : float
    press : float
    pprim : float
    temp : float
        The output temperature is either the input temperature
        or the species dependent initial temperature if the input
        temperature is less than 0 K.
    """
    mass, xlt, xltprim, press, pprim = None, None, None, None, None

    # If temp <=0, then use the species dependent starting point.
    if temp <= 0:
        temp = tstart[species]

    t = temp
    t2 = t * t
    t3 = t2 * t
    t4 = t2 * t2
    t5 = t4 * t
    t6 = t3 * t3

    if species == 'H2O':
        # 100
        mass = 18.
        xlt = 12420. - 4.8 * t
        xltprim = -4.8

        # from Marti & Mauersberger (1993 GRL 20, 363)
        p = -2663.5 / t + 12.537
        p = 10. * 10. ** p
        pp = (2663.5 / t2) * p

        # mass, xlt, xltprim, press, pprim = p5100(mass, xlt, xltprim, p, pp)

    elif species == 'H2O_CH4':
        mass = 18.
        xlt = 12160. + .5 * t - .033 * t2
        xltprim = 0.5 - 0.066 * t

        # from Marti & Mauersberger (1993 GRL 20, 363)
        p = -2663.5 / t + 12.537
        p = 10. * 10. ** p
        pp = (+2663.5 / t2) * p

        # mass, xlt, xltprim, press, pprim = p5100(mass, xlt, xltprim, p, pp)

    elif species == 'CO2':
        mass = 44.
        xlt = 6269. + 9.877 * t - .130997 * t2 + 6.2735e-4 * t3 - 1.2699e-6 * t4
        xltprim = 9.877 - .261994 * t + 1.88205e-3 * t2 - 5.0796e-6 * t3

        p = 21.3807649e0 - 2570.647e0/t - 7.78129489e4/t2 + 4.32506256e6/t3 - 1.20671368e8/t4 + 1.34966306e9/t5
        pp = 2570.647e0/t2 + 1.556258978e5/t3 - 12.97518768e6/t4 + 4.82685472e8/t5 - 6.7483153e9/t6
        p = dynmm * 10. ** p
        pp = pp * p

        # mass, xlt, xltprim, press, pprim = p5100(mass, xlt, xltprim, p, pp)

    elif species == 'CO':
        mass = 28
        if t > 68.127:
            sys.exit(f'error in CO temp, T = {t}')
        elif t > 61.544:
            xlt = 1855 + 3.253 * t - .06833 * t2
            xltprim = 3.253 - .13666 * t
            p = 16.8655152e0 - 748.151471e0/t - 5.84330795e0/t2 + 3.93853859e0/t3
            pp = 748.15147e0/t2 + 11.6866159e0/t3 - 11.81561577e0/t4
        elif t >= 14.0:
            xlt = 1893 + 7.331 * t + .01096 * t2 - .0060658 * t3 + 1.166e-4 * t4 - 7.8957e-7 * t5
            xltprim = 7.331 + .02192 * t - .0181974 * t2 + 4.664e-4 * t3 - 3.94785e-6 * t4

            p = 18.0741183e0 - 769.842078e0 / t - 12148.7759e0 / t2 + 2.7350095e5 / t3 - 2.9087467e6 / t4 + 1.20319418e7 / t5
            pp = 769.842078e0 / t2 + 24297.5518 / t3 - 820502.85e0 / t4 + 11634986.8e0 / t5 - 60159709.e0 / t6
            p = dynmm * 10. ** p
            pp = pp * p
        else:
            sys.exit(f'error in CO temp, T= {t}')

        # mass, xlt, xltprim, press, pprim = p5100(mass, xlt, xltprim, p, pp)

    mass = mass * proton
    xlt = xlt * ergcal
    xltprim = xltprim * ergcal
    press = p
    pprim = pp

    return species, mass, xlt, xltprim, press, pprim, temp


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
            n, species, a_v, a_ir, rh, inc, incl, temp, root, nflag, perc, gd)

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
        "rlog": rlog,
        "a_v": a_v,
        "a_ir": a_ir,
        "Zbar": zbar,
        "Zlog": zlog,
    }))

    return species, inc, rh, rlog, a_v, a_ir, zbar, zlog


def main_loop(n, species, a_v, a_ir, rh, inc, incl, temp, root, nflag, perc, gd):
    """
    Parameters
    ----------
    n : int
        counter (ranging from 0 to nb-1) which tracks
    species: str
        Inputted species
    a_v : float
        (a_v > 0)
        Inputted visual albedo
    a_ir : float
        Inputted infrared albedo
    rh: float
        Inputted heliocentric distance (in au)
    inc : float
        Inputted obliquity expressed in degrees
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
        gd is used in the case of nflag >= 5,000
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
        gd, perc = calc_perc(n, rh, inc, perc)
        return temp, gd, perc, nflag
    if nflag >= 5000:
        subdis = z[gd] * 4 * math.pi * rh * 149.6e11
        if perc != 0:
            if (subdis / perc) < 1e-02:
                z[n] = 0
                gd, perc = calc_perc(n, rh, inc, perc)
                return temp, gd, perc, nflag
            else:
                logging.error('error calculating sublimation')
                sys.exit(1)

    nflag += 1
    temp, gd, perc, nflag = main_loop(
        n, species, a_v, a_ir, rh, inc, incl, temp, root, nflag, perc, gd)

    return temp, gd, perc, nflag


def calc_perc(n, rh, inc, perc):
    logging.info('Incl: %f, Lat: %f', inc, b[n])
    gd = n
    perc = perc + z[n] * 4 * math.pi * rh * 149.6e11
    return gd, perc


############


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="NEED description", formatter_class=argparse.RawTextHelpFormatter,)
    parser.add_argument('species', type=str,
                        help="Desired ice species to be considered. "
                             "The valid inputs are: 'H2O', 'H20_CH4', 'CO2', and 'CO'")
    parser.add_argument('--a_v', metavar='visual albedo', type=float, required=True,)
    parser.add_argument('--a_ir', metavar='infrared albedo', type=float, required=True,)
    parser.add_argument('--rh', metavar='heliocentric_distance', type=float, required=True,)
    parser.add_argument('--inc', metavar='obliquity', type=float, required=True,)
    parser.add_argument('--temp', metavar='temperature', type=float, default=-1,
                        help="Not passing a starting temperature will default to a species dependent starting value:\n"
                             "H2O: 190 K\n"
                             "H2O_CH4: 190 K\n"
                             "CO2: 100 K\n"
                             "CO: 60 K\n")

    args = parser.parse_args()
    print(args)
    results = run_model(args.species, args.a_v, args.a_ir, args.rh, args.inc, args.temp)
