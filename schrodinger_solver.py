#!/usr/bin/env python3
"""Module that solves the onedimensional Schrödinger equation for arbitrary potentials"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.interpolate import interp1d

def solve1d(file):
    """"""
    # Read given parameters of "schrodinger.inp":
    schrodingerInp = open(file, "r")
    schrodingerLines = schrodingerInp.readlines()
    # print(schrodingerLines)

    massEntry = schrodingerLines[0].split()[0]
    mass = float(massEntry)
    # print("\nmass m:\n", mass)

    xMinEntry = schrodingerLines[1].split()[0]
    xMin = float(xMinEntry)
    # print("\nxMin:\n", xMin)

    xMaxEntry = schrodingerLines[1].split()[1]
    xMax = float(xMaxEntry)
    # print("\nxMax:\n", xMax)

    nPointEntry = schrodingerLines[1].split()[2]
    nPoint = int(nPointEntry)
    # print("\nPoint:\n", nPoint)

    firstEigvEntry = schrodingerLines[2].split()[0]
    firstEigv = int(firstEigvEntry)
    # print("\nfirst eigenvalue to print:\n", firstEigv)

    lastEigvEntry = schrodingerLines[2].split()[1]
    lastEigv = int(lastEigvEntry)
    # print("\nsecond eigenvalue to print:\n", secondEigv)

    interType = schrodingerLines[3].split()[0]
    # print("\ninterpolation type:\n", interType)

    potLen = len(schrodingerLines) - 5
    # print("\nnr. of pot-values:\n", potLen)
    xPot = np.zeros(potLen)
    yPot = np.zeros(potLen)
    for ii in range(5, potLen + 5):
        xPot[ii - 5] = float(schrodingerLines[ii].split()[0])
        yPot[ii - 5] = float(schrodingerLines[ii].split()[1])
    #print("\nx pot-values:\n", xPot, "\n\ny pot-values:\n", yPot)
    schrodingerInp.close()

    # Interpolation of the potential:
    if interType == "linear":
        print("\nlinear interpolation selected")
        xinterp = np.linspace(xMin, xMax, nPoint)
        yinterp = np.interp(xinterp, xPot, yPot)

    elif interType == "cspline":
        print("\ncspline interpolation selected")
        fCspline = interp1d(xPot, yPot, kind='cubic')
        xinterp = np.linspace(xMin, xMax, nPoint)
        yinterp = fCspline(xinterp)

    elif interType == "polynomial":
        polDegree = 2
        print("\npolynomial interpolation of degree", polDegree, "selected")
        coeffPoly = np.polyfit(xPot, yPot, polDegree)  # minimizes the squared error
        fPoly = np.poly1d(coeffPoly)  # read coeffizients
        xinterp = np.linspace(xMin, xMax, nPoint)
        yinterp = fPoly(xinterp)

    else:
        print("error: invalid interpolation method")

    # write obtained points to file "potential.dat":
    InterpPot = np.hstack((xinterp.reshape((-1, 1)), yinterp.reshape((-1, 1))))
    np.savetxt("potential.dat", InterpPot)


    #eigenvalues and eigenvectors:
    delta = abs(xinterp[1]-xinterp[0])
    # print("\nlength of xinterp:\n", len(xinterp), "\n\nobtained delta:\n", delta)
    a = 1/(mass*(delta)**2)
    ludiag = np.zeros(nPoint-1)
    ludiag[:] = (-a)/2
    # print("\nlower and upper diagonal:\n", ludiag)
    mainDiag = np.zeros(nPoint)

    for jj in range(0, nPoint):
        mainDiag[jj] = a + yinterp[jj]


    # eigenvalues in ascending order, the normalized eigenvector
    # corresponding to the eigenvalue eiva[i] is the column eive[:,i]
    eiva, eive = eigh_tridiagonal(mainDiag, ludiag, select='a', select_range=None)
    # print("\neigenvalues:\n", eiva, "\neigenvectors\n:", eive)

    # write eigenvalues und eigenvectors to files:
    np.savetxt("energies.dat", eiva[(firstEigv-1):(lastEigv)])
    wavefunctions = np.hstack((xinterp.reshape((-1, 1)), eive))
    np.savetxt("wavefuncs.dat", wavefunctions[:, (firstEigv-1):(lastEigv +1)])

    #expectation values and uncertainties

    #check norm, norm factor is 1/sqrt(delta), use squared for less machine rounding errors
    '''
    for ii in range(lastEigv-firstEigv + 1):
        print("\n",ii +1, np.sum(eive[:,0])**2
    #print(eive[:,0])
    #print(np.abs(eive[:,0])**2)
    '''
    expectation_values = np.zeros(int(lastEigv-firstEigv + 1))
    expectation_values_squared = np.zeros(int(lastEigv-firstEigv + 1))
    deviation_x = np.zeros(int(lastEigv-firstEigv + 1))

    for ii in range(int(lastEigv-firstEigv + 1)):
        expectation_values[ii] = np.sum(xinterp[:]*((eive[:, ii])**2))
        expectation_values_squared[ii] = np.sum(((xinterp[:])**2)*(eive[:, ii])**2)
        deviation_x[ii] = np.sqrt(expectation_values_squared[ii] - (expectation_values[ii])**2)

    #print("\nexpectationvalues:\n",expectation_values, "\ndevitaion for x: ", deviation_x)
    #write expectationvalues and deviation to file
    expvalues = np.hstack((expectation_values.reshape((-1, 1)), deviation_x.reshape((-1, 1))))
    np.savetxt("expvalues.dat", expvalues)

    return eiva, xinterp, yinterp, xPot, yPot
