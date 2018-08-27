#!/usr/bin/env python3
"""Contains routines to test the solvers modules"""

from os import path
import numpy as np
import pytest

import schrodinger_solver


TOLERANCE_ENERGIES = 1e-2
TOLERANCE_INTERP = 1e-1

INPUT_DIR = path.join("tests", "input")
EXPECTED_DIR = path.join("tests", "expected")

INTERP_TESTS = ["finite_pot_box_interp", "inf_pot_box_interp",
                "harm_osz_interp", "asym_harm_osz_interp", "double_pot_linear_interp",
                "double_pot_spline_interp"]
ENERGY_TESTS = ["finite_pot_box_energies", "inf_pot_box_energies",
                "harm_osz_energies"]


@pytest.mark.parametrize("testname", INTERP_TESTS)
def test_interpolation(testname):
    """Tests the interpolation of schrodinger_solver by comparing the
    interpolated potential with the given XY data"""
    input_files = "{}.{}".format(path.join(INPUT_DIR, testname), "inp")
  
    eiva, xinterp, yinterp, xPot, yPot = schrodinger_solver.solve1d(input_files)
    given_Pot = np.empty([len(xPot), 2])
    given_Pot[:, 0] = xPot
    given_Pot[:, 1] = yPot
    
    delta = abs(xinterp[1]-xinterp[0])
    test_interp = True
    
    for ii in given_Pot:
        x = ii[0]
        y = ii[1]
        # boolean array "bool_array"
        bool_array = np.abs(xinterp - x) <= delta / 2
        
        # calculate average, if x-value is exactly between two points
        # otherwise average is just the y value
        interval_average = np.sum(yinterp[bool_array]) / len(yinterp[bool_array])
        test_interp = test_interp and (np.abs(interval_average - y) < TOLERANCE_INTERP).all()

    assert test_interp

if __name__ == "__main__":
    pytest.main()
