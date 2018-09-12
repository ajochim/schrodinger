#!/usr/bin/env python3
"""Program that solves the onedimensional Schrodinger equation for arbitrary potentials"""

import argparse
import schrodinger_io as io
import schrodinger_solver as solver
import schrodinger_visualize as visualize

# optional arguments
_DESCRIPTION = "Program that solves the onedimensional Schrodinger equation \
                for arbitrary potentials"
PARSER = argparse.ArgumentParser(description=_DESCRIPTION)
MSG = "Directory where schrodinger.inp is located: (default: .)"
PARSER.add_argument("-d", "--directory", default=".", help=MSG)
MSG = "Splitting wavefunctions in plot: (default: False)"
PARSER.add_argument("-s", "--split", default=False, action="store_true", help=MSG)
MSG = "Stretch wavefunctions in plot: (default: 1.0)"
PARSER.add_argument("-st", "--stretch", default=1.0, help=MSG)
MSG = "Change Markersize in plot: (default: 10.0)"
PARSER.add_argument("-m", "--markersize", default=10.0, help=MSG)
ARGS = PARSER.parse_args()
PATH = "{}/schrodinger.inp".format(ARGS.directory)
print("\nstarting parameters: ")
if ARGS.split:
    print("splitting wavefunctions in plot set to 'True'")
else:
    print("splitting wavefunctions in plot is currently set to 'False' \t(default)")
if float(ARGS.stretch) != float(PARSER.get_default("stretch")):
    print("changed stretchfactor to ", ARGS.stretch)
else:
    print("stretchfactor is currently set to ", PARSER.get_default("stretch"), "\t\t\t\t(default)")
if float(ARGS.markersize) != float(PARSER.get_default("markersize")):
    print("changed markersize to: ", str(ARGS.markersize))
else:
    print("markersize is currently set to ", PARSER.get_default("markersize"), "\t\t\t\t(default)")

# main program
OBTAINED_INPUT = io.read_input(PATH)
INTERPOT = solver.interpolate(OBTAINED_INPUT)
CALCULATED = solver.solve1d(OBTAINED_INPUT, INTERPOT)
io.output(CALCULATED)
visualize.show(stretchfactor=3, split=True, markersize=15)
