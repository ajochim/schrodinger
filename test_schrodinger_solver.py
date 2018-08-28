#!/usr/bin/env python3
"""Contains routines to test the solvers modules"""

from os import path
import numpy as np
import pytest

import schrodinger_solver


TOLERANCE_ENERGIES = 8e-3
TOLERANCE_INTERP = 7e-2

INPUT_DIR = path.join("tests", "input")
EXPECTED_DIR = path.join("tests", "expected")

INPUT_INTERP = ["finite_pot_box", "inf_pot_box", "harm_osz", "asym_harm_osz",
         "double_pot_linear", "double_pot_spline"]
INPUT_ENERGY = ["finite_pot_box", "inf_pot_box", "harm_osz"]
ENERGY_TESTS = ["finite_pot_box_energies", "inf_pot_box_energies",
                "harm_osz_energies"]


@pytest.mark.parametrize("testname_interp", INPUT_INTERP)
def test_interpolation(testname_interp):
    """Tests the interpolation of schrodinger_solver by comparing the
    interpolated potential with the given XY data"""
    input_files_interp = "{}.{}".format(path.join(INPUT_DIR, testname_interp), "inp")
  
    eiva, xinterp, yinterp, xPot, yPot = schrodinger_solver.solve1d(input_files_interp)
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
        
        # extraordinary treatment of "finite_pot_box" discontinuities:
#        if testname == "finite_pot_box_interp":
#            index_array_true = np.where(bool_array)[0]
#            # wachsend = True; fallend = False
#            bool_asc_desc = (yinterp[index_array_true + 1] - yinterp[index_array_true - 1]) > 0
#            print("aufsteigend,absteigend", bool_asc_desc)
#            for jj in index_array_true:
#                if bool_asc_desc[jj] == True:
#                    buffer = bool_array[jj]
#                    bool_array[jj] = bool_array[jj + 1]
#                    bool_array[jj + 1] = buffer
#                elif bool_asc_desc[jj] == False:
#                    buffer = bool_array[jj]
#                    bool_array[jj] = bool_array[jj - 1]
#                    bool_array[jj - 1] = buffer
        
        # extraordinary treatment of "finite_pot_box" discontinuities:
        if testname_interp == "finite_pot_box":
            print("finite_pot_box discontinuities")
            bool_array_right = np.roll(bool_array, 1)
            bool_array_left = np.roll(bool_array, -1)
            interval_average_right = np.sum(yinterp[bool_array_right]) / len(yinterp[bool_array_right])
            interval_average_left = np.sum(yinterp[bool_array_left]) / len(yinterp[bool_array_left])
            interval_average = np.sum(yinterp[bool_array]) / len(yinterp[bool_array])
            
            inTol_leftright = ((np.abs(interval_average_right - y) < TOLERANCE_INTERP).all() or 
                      (np.abs(interval_average_left - y) < TOLERANCE_INTERP).all())
            inTol = (np.abs(interval_average - y) < TOLERANCE_INTERP).all() or inTol_leftright
            test_interp = test_interp and inTol
        else:
            interval_average = np.sum(yinterp[bool_array]) / len(yinterp[bool_array])
            test_interp = test_interp and (np.abs(interval_average - y) < TOLERANCE_INTERP).all()
            
#        # extraordinary treatment of discontinuities:
#        if test_interp == False:
#            bool_array_right = np.roll(bool_array, 1)
#            interval_average = np.sum(yinterp[bool_array_right]) / len(yinterp[bool_array_right])
#            test_interp = test_interp and (np.abs(interval_average - y) < TOLERANCE_INTERP).all()
#        if test_interp == False:
#            bool_array_left = np.roll(bool_array, -1)
#            interval_average = np.sum(yinterp[bool_array_left]) / len(yinterp[bool_array_left])
#            test_interp = test_interp and (np.abs(interval_average - y) < TOLERANCE_INTERP).all()
        
    assert test_interp

@pytest.mark.parametrize("testname_energie", INPUT_ENERGY)
def test_energies(testname_energie):
    """Tests the energy-levels of schrodinger_solver by comparing thoose
    with exact results or groundstates"""
    input_files_energies = "{}.{}".format(path.join(INPUT_DIR, testname_energie), "inp")
    
    reference_files = "{}{}.{}".format(path.join(EXPECTED_DIR, testname_energie), "_energies", "out")
    eiva_ref = np.loadtxt(reference_files)
    
    eiva, xinterp, yinterp, xPot, yPot = schrodinger_solver.solve1d(input_files_energies)
    eiva_restricted = eiva[0:len(eiva_ref)]
    
    test_energies = True
    test_energies = test_energies and (np.abs(eiva_restricted - eiva_ref) < TOLERANCE_ENERGIES).all()
    
    assert test_energies
    
    
    
if __name__ == "__main__":
    pytest.main()









