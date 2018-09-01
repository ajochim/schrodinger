#!/usr/bin/env python3
"""Module that reads the input data, processed by the schrodinger solver"""

import sys
import numpy as np

def read_input(file):
    """Read given parameters and potential data of 'schrodinger.inp'"""

    # Read given parameters of "schrodinger.inp":
    try:
        schrodingerinp = open(file, "r")
        schrodingerlines = schrodingerinp.readlines()
        # print(schrodingerlines)
    except OSError:
        print("Could not open specified inputfile: {}".format(file))
        print("File is either not present or it exists but has wrong permissions\nExiting program")
        sys.exit(1)

    massentry = schrodingerlines[0].split()[0]
    mass = float(massentry)
    # print("\nmass m:\n", mass)

    xminentry = schrodingerlines[1].split()[0]
    xmin = float(xminentry)
    # print("\nxmin:\n", xmin)

    xmaxentry = schrodingerlines[1].split()[1]
    xmax = float(xmaxentry)
    # print("\nxmax:\n", xmax)

    npointentry = schrodingerlines[1].split()[2]
    npoint = int(npointentry)
    # print("\npoint:\n", npoint)

    firsteigventry = schrodingerlines[2].split()[0]
    firsteigv = int(firsteigventry)
    # print("\nfirst eigenvalue to print:\n", firsteigv)

    lasteigventry = schrodingerlines[2].split()[1]
    lasteigv = int(lasteigventry)
    # print("\nsecond eigenvalue to print:\n", secondEigv)

    intertype = schrodingerlines[3].split()[0]
    # print("\ninterpolation type:\n", intertype)

    potlen = len(schrodingerlines) - 5
    # print("\nnr. of pot-values:\n", potlen)
    xpot = np.zeros(potlen)
    ypot = np.zeros(potlen)
    for ii in range(5, potlen + 5):
        xpot[ii - 5] = float(schrodingerlines[ii].split()[0])
        ypot[ii - 5] = float(schrodingerlines[ii].split()[1])
    xypot = np.hstack((xpot.reshape((-1, 1)), ypot.reshape((-1, 1))))

    schrodingerinp.close()

    obtained_input = {"mass": mass, "xmin": xmin, "xmax": xmax, "npoint": npoint,
                      "firsteigv": firsteigv, "lasteigv": lasteigv, "intertype": intertype,
                      "potlen": potlen, "xypot": xypot}

    return obtained_input
