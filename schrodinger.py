#!/usr/bin/env python3
"""Program that solves the onedimensional Schrodinger equation for arbitrary potentials"""

import schrodinger_io as io
import schrodinger_solver as solver
import schrodinger_visualize as visualize

OBTAINED_INPUT = io.read_input("schrodinger.inp")
INTERPOT = solver.interpolate(OBTAINED_INPUT)
CALCULATED = solver.solve1d(OBTAINED_INPUT, INTERPOT)
io.output(CALCULATED)

#plot solutions, stretfactor to scale wavefunktions,
#split to move wavefunctions to eigenvalues
#markersize for expectation values and uncertainty
visualize.show(stretchfactor=3, split=True, markersize=15)
