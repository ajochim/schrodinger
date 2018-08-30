#!/usr/bin/env python3
"""Module that solves the onedimensional Schrodinger equation for arbitrary potentials"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.interpolate import interp1d

def solve1d(file):
    """Solves the 1D schrodinger equation with potential/parameters given
       by the input file schrodinger.inp"""

    # Read given parameters of "schrodinger.inp":
    schrodingerinp = open(file, "r")
    schrodingerlines = schrodingerinp.readlines()
    # print(schrodingerlines)

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
    #print("\nx pot-values:\n", xpot, "\n\ny pot-values:\n", ypot)
    schrodingerinp.close()

    # Interpolation of the potential:
    if intertype == "linear":
        print("\nlinear interpolation selected")
        xinterp = np.linspace(xmin, xmax, npoint)
        yinterp = np.interp(xinterp, xpot, ypot)

    elif intertype == "cspline":
        print("\ncspline interpolation selected")
        fcspline = interp1d(xpot, ypot, kind='cubic')
        xinterp = np.linspace(xmin, xmax, npoint)
        yinterp = fcspline(xinterp)

    elif intertype == "polynomial":
        poldegree = 2
        print("\npolynomial interpolation of degree", poldegree, "selected")
        coeffpoly = np.polyfit(xpot, ypot, poldegree)  # minimizes the squared error
        fpoly = np.poly1d(coeffpoly)  # read coeffizients
        xinterp = np.linspace(xmin, xmax, npoint)
        yinterp = fpoly(xinterp)

    else:
        print("error: invalid interpolation method")

    # write obtained points to file "potential.dat":
    interppot = np.hstack((xinterp.reshape((-1, 1)), yinterp.reshape((-1, 1))))
    np.savetxt("potential.dat", interppot)


    #eigenvalues and eigenvectors:
    delta = abs(xinterp[1]-xinterp[0])
    # print("\nlength of xinterp:\n", len(xinterp), "\n\nobtained delta:\n", delta)
    aa = 1/(mass*(delta)**2)
    ludiag = np.zeros(npoint-1)
    ludiag[:] = (-aa)/2
    # print("\nlower and upper diagonal:\n", ludiag)
    maindiag = np.zeros(npoint)

    for jj in range(0, npoint):
        maindiag[jj] = aa + yinterp[jj]


    # eigenvalues in ascending order, the normalized eigenvector
    # corresponding to the eigenvalue eiva[i] is the column eive[:,i]
    eiva, eive = eigh_tridiagonal(maindiag, ludiag, select='a', select_range=None)
    # print("\neigenvalues:\n", eiva, "\neigenvectors\n:", eive)

    # write eigenvalues und eigenvectors to files:
    np.savetxt("energies.dat", eiva[(firsteigv-1):(lasteigv)])
    wavefunctions = np.hstack((xinterp.reshape((-1, 1)), eive))
    np.savetxt("wavefuncs.dat", wavefunctions[:, (firsteigv-1):(lasteigv +1)])

    #expectation values and uncertainties

    #check norm, norm factor is 1/sqrt(delta), use squared for less machine rounding errors
    '''
    for ii in range(lasteigv-firsteigv + 1):
        print("\n",ii +1, np.sum(eive[:,0])**2
    #print(eive[:,0])
    #print(np.abs(eive[:,0])**2)
    '''
    expectation_values = np.zeros(int(lasteigv-firsteigv + 1))
    expectation_values_squared = np.zeros(int(lasteigv-firsteigv + 1))
    deviation_x = np.zeros(int(lasteigv-firsteigv + 1))

    for ii in range(int(lasteigv-firsteigv + 1)):
        expectation_values[ii] = np.sum(xinterp[:]*((eive[:, ii])**2))
        expectation_values_squared[ii] = np.sum(((xinterp[:])**2)*(eive[:, ii])**2)
        deviation_x[ii] = np.sqrt(expectation_values_squared[ii] - (expectation_values[ii])**2)

    #print("\nexpectationvalues:\n",expectation_values, "\ndevitaion for x: ", deviation_x)
    #write expectationvalues and deviation to file
    expvalues = np.hstack((expectation_values.reshape((-1, 1)), deviation_x.reshape((-1, 1))))
    np.savetxt("expvalues.dat", expvalues)

    return eiva, xinterp, yinterp, xpot, ypot
