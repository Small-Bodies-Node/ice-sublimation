import math
from subprocess import *
import sys

nb = 41
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

tstart = [190., 190., 100., 60.]


def runModel(index, a0, a1, rh, inc, temp=-1, perc=0, gd=None):
    print("Input Parameters:")
    print(f'Species = {species[index]}, Avis = {a0}, Air = {a1}, r_H = {rh}, Incl = {inc}')
    p10(index, a0, a1, rh, inc, temp, perc, gd)


def p10(index, a0, a1, rh, inc, temp, perc, gd):
    if a0 < 0:
        p800()
    incl = (90 - inc) * math.pi / 180

    # If temp <=0, then use the species dependent starting point.
    if temp <= 0:
        temp = tstart[index]

    spec, mass, xlt, xltprim, press, pprim = sublime(index, temp)
    root = 1 / math.sqrt(mass * 2 * math.pi * boltz)

    nflag = 1
    for n in range(1, nb + 1):
        # do all up to p600
        temp, gd = p100(n, index, a0, a1, rh, inc, incl, temp, root, nflag, perc, gd)

    # now we run p600
    n = nb
    p600(n, index, a0, a1, rh, inc, incl, temp, root, nflag, perc, gd)


def p100(n, index, a0, a1, rh, inc, incl, temp, root, nflag, perc, gd):
    rootT = math.sqrt(temp)
    if b[n - 1] <= -incl:
        frac[n - 1] = 0
        z[n - 1] = 0
        # go to 600
        # which ends this iteration
        return temp, gd
    elif b[n - 1] > incl:
        frac[n - 1] = sb[n - 1] * math.cos(incl)
    else:
        x1 = math.cos(incl) * sb[n - 1] * (math.acos(-math.tan(b[n - 1]) * (1 / math.tan(incl)))) / math.pi
        x2 = math.sin(incl) * math.cos(b[n - 1]) * math.sin(math.acos(-math.tan(b[n - 1]) / math.tan(incl))) / math.pi
        frac[n - 1] = x1 + x2

    spec, mass, xlt, xltprim, press, pprim = sublime(index, temp)

    sun = f0 * frac[n - 1] * (1 - a0) / rh ** 2
    # print("a1: ",a1, "sigma: ",sigma, "temp: ", temp, "q: ", q, "root: ",root,"rootT: ",root,"press: ",press,"xlt: ",xlt)
    radiat = (1 - a1) * sigma * temp ** 4
    evap = q * root / rootT * press * xlt
    phi = radiat + evap - sun
    z[n - 1] = max(evap / xlt, 1e-30)

    # print("radiat", radiat,"evap: ",evap)
    drad = 4 * radiat / temp
    x1 = pprim * xlt
    x2 = press * xltprim

    devap = q * root / rootT * (x1 + x2)
    phipri = drad + devap

    dt = math.copysign(min(10, abs(phi / phipri / 2)), phi / phipri)
    # print("phi: ", phi, "dt: ", dt)
    tp = temp - dt
    temp = tp
    # print("phi:", phi, "sun", sun)
    # print("abs(phi/sun): ", abs(phi / sun), "abs(phi)", abs(phi))
    print("new temp: ", temp)
    if (abs(phi / sun) < 1e-4 or abs(phi) < 1e-4):
        # go to 300
        gd, perc = p300(n, rh, inc, perc)
    if nflag >= 10000000:
        subdis = z[gd - 1] * 4 * math.pi * rh * 149.6e11
        if perc != 0:
            if (subdis / perc) < 1e-02:
                z[n - 1] = 0
                # go to 300
                p300(n, rh, inc, perc)
            else:
                sys.exit('error calculating sublimation')

    nflag += 1
    # go to 100
    # start next iteration
    print("*" * 20)
    return temp, gd


def p300(n, rh, inc, perc):
    print(f"Incl: {inc}, Lat: {b[n - 1]}")
    gd = n
    perc = perc + z[n - 1] * 4 * math.pi * rh * 149.6e11
    return gd, perc


def p600(n, index, a0, a1, rh, inc, incl, temp, root, nflag, perc, gd):
    zbar = 0
    for nn in range(1, nb - 1 + 1):
        zbar = zbar + 0.5 * (z[nn - 1] + z[nn + 1 - 1]) * delsb

    p700(n, index, a0, a1, rh, inc, incl, temp, root, nflag, perc, gd, zbar)


def p700(n, index, a0, a1, rh, inc, incl, temp, root, nflag, perc, gd, zbar):
    zbar = zbar / 2
    zlog = math.log10(zbar)
    rlog = math.log10(rh)
    print("Final Result")
    print(
        f"Species: {species[index]}, Obliquity: {inc}, r_H: {rh}, rlog: {rlog:.3f}, a0: {a0}, a1: {a1}, Zbar: {zbar:.2e}, Zlog: {zlog:2.2f}")

    a0 = -1
    p10(index, a0, a1, rh, inc, temp, perc, gd)


def p800():
    sys.exit("Run complete")


############
a0_in = 0.05
a1_in = 0.00
rh_in = 3.98
inc_in = 90

species = ['H2O', 'H2O_CH4', 'CO2', 'CO']
index_in = species.index('H2O')

runModel(index_in, a0_in, a1_in, rh_in, inc_in,)
