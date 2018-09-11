#!/usr/bin/env python3
"""Program that solves the onedimensional Schrodinger equation for arbitrary potentials"""

import numpy as np

import schrodinger_in as io
import schrodinger_solver as solver
import schrodinger_visualize as visualize


#creates inputfile
def create_inputfile():
    """creates an input file via this script for faster use"""
    mass = 3
    xmin = 0
    xmax = 4*np.pi
    npoints = 1999
    firsteiva, lasteiva = 1, 6
    interpolation_type = "cspline"
    xx = np.linspace(xmin, xmax, npoints)

    potential = np.sin(xx)

    file = open("schrodinger_created.inp", "w")
    file.write(str(mass) + "\n")
    file.write(str(xmin) + " " + str(xmax) + " " + str(npoints) + "\n")
    file.write(str(firsteiva) + " " + str(lasteiva) + "\n")
    file.write(interpolation_type + "\n")
    file.write(str(npoints) + "\n")
    for ii in range(npoints):
        file.write(str(xx[ii]) + " " + str(potential[ii]) + "\n")

#create_inputfile()
#solves equation with given input file
#solver.solve1d(mass, npoint, firsteigv, lasteigv, interpot)
#schrodinger_solver.solve1d("schrodinger_created.inp")#if used create_inputfile()

OBTAINED_INPUT = io.read_input("/tests/input/schrodinger.inp")

#OBTAINED_INPUT = io.read_input("schrodinger_infwell.inp")
#OBTAINED_INPUT = io.read_input("schrodinger_potwell.inp")
#OBTAINED_INPUT = io.read_input("schrodinger_doublewell.inp")
#OBTAINED_INPUT = io.read_input("schrodinger_doublewellspline.inp")
#OBTAINED_INPUT = io.read_input("schrodinger_asymwell.inp")
#OBTAINED_INPUT = io.read_input("schrodinger_harmosz.inp")


INTERPOT = solver.interpolate(OBTAINED_INPUT)
solver.solve1d(OBTAINED_INPUT, INTERPOT)

#plot solutions, stretfactor to scale wavefunktions,
#split to move wavefunctions to eigenvalues
#markersize for expectation values and uncertainty
visualize.show(stretchfactor=3, split=True, markersize=15)
