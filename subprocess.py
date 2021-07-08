'''
Description
-----------
Calculates the latent heat of sublimation and the vapor pressure of
the solid for various ices and the derivatives thereof


Modification History
--------------------
Mark Van Seluos July 2021
- Rewrote sublime.f in python as subprocess.py

Brian Prager 06/25/2010
- Commented out the C02 pressures derived by G.N. Brown Jr. and W.T. Ziegler
in "Vapor Pressure and the Heats of Vaporization and Sublimation of Liquids
and Solids of interest in cryogenics below 1-atm pressure" (Advances in
Cryogenic Engineering, 25, 1980) [Also used to derive C0]. The script could
end up looking at temperatures below the minimum accepted value for their
derivation, and would crash the script. The pressure was then replaced with
the Antoine Equation, which is justified in another paper. (Given in comments
above the equation).

'''

import sys

# Constants
sigma = 5.67e-5
f0 = 1.39e6
dynmm = 1.333e3
boltz = 1.38e-16
ergcal = 6.953e-17
proton = 1.67e-24


# If temp <=0, then use the species dependent starting point.

def sublime(species, temp):
    # ###Setup
    mass, xlt, xltprim, press, pprim = None, None, None, None, None

    tstart = {'H2O': 190,
              'H2O_CH4': 190,
              'CO2': 100,
              'CO': 60,
              }
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

        # There is an alternative formulation for pressure which is commented out in fortran

        mass, xlt, xltprim, press, pprim = p5100(mass, xlt, xltprim, p, pp)

    elif species == 'H2O_CH4':
        mass = 18.
        xlt = 12160. + .5 * t - .033 * t2
        xltprim = 0.5 - 0.066 * t

        # from Marti & Mauersberger (1993 GRL 20, 363)
        p = -2663.5 / t + 12.537
        p = 10. * 10. ** p
        pp = (+2663.5 / t2) * p

        # There is an alternative formulation for pressure which is commented out in fortran

        mass, xlt, xltprim, press, pprim = p5100(mass, xlt, xltprim, p, pp)

    elif species == 'CO2':
        mass = 44.
        xlt = 6269. + 9.877 * t - .130997 * t2 + 6.2735e-4 * t3 - 1.2699e-6 * t4
        xltprim = 9.877 - .261994 * t + 1.88205e-3 * t2 - 5.0796e-6 * t3

        p = 21.3807649e0 - 2570.647e0/t - 7.78129489e4/t2 + 4.32506256e6/t3 - 1.20671368e8/t4 + 1.34966306e9/t5
        pp = 2570.647e0/t2 + 1.556258978e5/t3 - 12.97518768e6/t4 + 4.82685472e8/t5 - 6.7483153e9/t6
        p = dynmm * 10. ** p
        pp = pp * p

        mass, xlt, xltprim, press, pprim = p5100(mass, xlt, xltprim, p, pp)

    elif species == 'CO':
        mass = 28
        if t > 68.127:
            sys.exit(f'error in CO temp, T = {t}')
            # f: stop

        elif t > 61.544:
            xlt = 1855 + 3.253 * t - .06833 * t2
            xltprim = 3.253 - .13666 * t

            p = 16.8655152e0 - 748.151471e0/t - 5.84330795e0/t2 + 3.93853859e0/t3
            pp = 748.15147e0/t2 + 11.6866159e0/t3 - 11.81561577e0/t4

            # Before I ran p5100 twice
            # I believe the return in Fortran ends the code so we only need 5100 once.
            # mass, xlt, xltprim, press, pprim = p5100(mass, xlt, xltprim, p, pp)

        elif t >= 14.0:
            xlt = 1893 + 7.331 * t + .01096 * t2 - .0060658 * t3 + 1.166e-4 * t4 - 7.8957e-7 * t5
            xltprim = 7.331 + .02192 * t - .0181974 * t2 + 4.664e-4 * t3 - 3.94785e-6 * t4

            p = 18.0741183e0 - 769.842078e0 / t - 12148.7759e0 / t2 + 2.7350095e5 / t3 - 2.9087467e6 / t4 + 1.20319418e7 / t5
            pp = 769.842078e0 / t2 + 24297.5518 / t3 - 820502.85e0 / t4 + 11634986.8e0 / t5 - 60159709.e0 / t6
            p = dynmm * 10. ** p
            pp = pp * p
        else:
            sys.exit(f'error in CO temp, T= {t}')
            # f: stop

        mass, xlt, xltprim, press, pprim = p5100(mass, xlt, xltprim, p, pp)

    return species, mass, xlt, xltprim, press, pprim, temp


def p5100(mass, xlt, xltprim, p, pp):
    mass = mass * proton
    xlt = xlt * ergcal
    xltprim = xltprim * ergcal
    press = p
    pprim = pp

    return mass, xlt, xltprim, press, pprim
