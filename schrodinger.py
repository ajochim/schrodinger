#!/usr/bin/env python3
"""Program that solves the onedimensional Schrodinger equation for arbitrary potentials"""

import argparse
import schrodinger_io as io
import schrodinger_solver as solver
import schrodinger_visualize as visualize

# optionsl arguments
_DESCRIPTION = "Optional arguments for schrodinger.py"
PARSER = argparse.ArgumentParser(description=_DESCRIPTION)
MSG = "Directory where schrodinger.inp is located: (default: .)"
PARSER.add_argument("-d", "--directory", default=".", help=MSG)
MSG = "Splitting wavefunctions in plot: (default: False)"
PARSER.add_argument("-s", "--split", default=False, action="store_true", help=MSG)
MSG = "Stretch wavefunctions in plot: (default: 3.0)"
PARSER.add_argument("-st", "--stretch", default=1.0, help=MSG)
MSG = "Change Markersize in plot: (default: 10.0)"
PARSER.add_argument("-m", "--markersize", default=10.0, help=MSG)
ARGS = PARSER.parse_args()
PATH = "{}/schrodinger.inp".format(ARGS.directory)
if ARGS.split:
    print("Splitting wavefunctions in plot set to True")
print("Stretchfactor: ", str(ARGS.stretch), "Markersize: ", str(ARGS.markersize))

# main program
OBTAINED_INPUT = io.read_input(PATH)
INTERPOT = solver.interpolate(OBTAINED_INPUT)
CALCULATED = solver.solve1d(OBTAINED_INPUT, INTERPOT)
io.output(CALCULATED)
visualize.show(stretchfactor=float(ARGS.stretch),
               split=ARGS.split, markersize=float(ARGS.markersize))
