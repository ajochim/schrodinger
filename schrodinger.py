# -*- coding: utf-8 -*-
import schrodinger_solver
import schrodinger_visualize
import numpy as np

#creates inputfile
def create_inputfile():# creates an input file via this skript for faster use
    mass = 3
    xMin = 0
    xMax = 4*np.pi
    nPoints = 1999
    firstEiva, lastEiva = 1, 6
    interpolation_type = "cspline"
    xx = np.linspace(xMin, xMax, nPoints)
    
    potential = np.sin(xx)
    
    file = open("schrodinger_created.inp", "w")
    file.write(str(mass) + "\n")
    file.write(str(xMin) + " " + str(xMax) + " " + str(nPoints) + "\n")
    file.write(str(firstEiva) + " " + str(lastEiva) + "\n")
    file.write(interpolation_type + "\n")
    file.write(str(nPoints) + "\n")
    for ii in range(nPoints):
        file.write(str(xx[ii]) + " " + str(potential[ii]) + "\n")

#create_inputfile()

#solves equation with given input file
schrodinger_solver.solve1d("schrodinger.inp")
#schrodinger_solver.solve1d("schrodinger_created.inp")#if used create_inputfile()

#schrodinger_solver.solve1d("schrodinger_infwell.inp")
#schrodinger_solver.solve1d("schrodinger_potwell.inp")
#schrodinger_solver.solve1d("schrodinger_doublewell.inp")
#schrodinger_solver.solve1d("schrodinger_doublewellspline.inp")
#schrodinger_solver.solve1d("schrodinger_asymwell.inp")
#schrodinger_solver.solve1d("schrodinger_harmosz.inp")


#plot solutions, stretfactor to scale wavefunktions,
#split to move wavefunctions to eigenvalues
#markersize for expectation values and uncertainty
schrodinger_visualize.visualize(stretchfactor=2, split=True, markersize=15)