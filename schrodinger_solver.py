#!/usr/bin/env python3
"""Module that solves the onedimensional Schrodinger equation for arbitrary potentials"""

import sys
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.interpolate import interp1d


def interpolate(obtained_input):
    """Interpolation routine, that interpolates the given potential data with
       a specified interpolation method"""

    intertype = obtained_input["intertype"]
    xmin = obtained_input["xmin"]
    xmax = obtained_input["xmax"]
    npoint = obtained_input["npoint"]
    xypot = obtained_input["xypot"]
    xpot = xypot[:, 0]
    ypot = xypot[:, 1]

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
        coeffpoly = np.polyfit(xpot, ypot, poldegree)
        fpoly = np.poly1d(coeffpoly)
        xinterp = np.linspace(xmin, xmax, npoint)
        yinterp = fpoly(xinterp)

    else:
        print("error: invalid interpolation method")
        sys.exit(1)

    #format interpolated values
    interppot = np.hstack((xinterp.reshape((-1, 1)), yinterp.reshape((-1, 1))))

    return interppot

def solve1d(obtained_input, pot):
    """Solves the 1D schrodinger equation with potential/parameters given
       by the input file schrodinger.inp and returns dictionary"""

    mass = obtained_input["mass"]
    npoint = obtained_input["npoint"]
    firsteigv = obtained_input["firsteigv"]
    lasteigv = obtained_input["lasteigv"]

    xinterp = pot[:, 0]
    yinterp = pot[:, 1]

    #eigenvalues and eigenvectors:
    delta = abs(xinterp[1]-xinterp[0])
    aa = 1/(mass*(delta)**2)
    ludiag = np.zeros(npoint-1)
    ludiag[:] = (-aa)/2
    maindiag = np.zeros(npoint)

    for jj in range(0, npoint):
        maindiag[jj] = aa + yinterp[jj]

    # eigenvalues in ascending order, the normalized eigenvector
    # corresponding to the eigenvalue eiva[i] is the column eive[:,i]
    eiva, eive = eigh_tridiagonal(maindiag, ludiag, select='a', select_range=None)

    # reshape wavefunctions
    wavefunctions = np.hstack((xinterp.reshape((-1, 1)), eive))

    #expectation values and uncertainties
    expectation_values = np.zeros(int(lasteigv - firsteigv + 1))
    expectation_values_squared = np.zeros(int(lasteigv - firsteigv + 1))
    deviation_x = np.zeros(int(lasteigv - firsteigv + 1))

    for ii in range(int(lasteigv-firsteigv + 1)):
        expectation_values[ii] = np.sum(xinterp[:]*((eive[:, ii])**2))
        expectation_values_squared[ii] = np.sum(((xinterp[:])**2)*(eive[:, ii])**2)
        deviation_x[ii] = np.sqrt(expectation_values_squared[ii] - (expectation_values[ii])**2)

    #reshape expectationvalues and deviation
    expvalues = np.hstack((expectation_values.reshape((-1, 1)), deviation_x.reshape((-1, 1))))

    data = {} #dicitionary for all calculated values
    data.update({"potential": pot})
    data.update({"energies": eiva[(firsteigv - 1):(lasteigv)]})
    data.update({"wavefuncs": wavefunctions[:, (firsteigv - 1):(lasteigv + 1)]})
    data.update({"expvalues": expvalues})

    return data
