#!/usr/bin/env python3
"""Contains routines to test the solvers modules"""

from os import path
import numpy as np
import pytest

import schrodinger_io as io
import schrodinger_solver as solver


TOLERANCE_INTERP = 7e-2
TOLERANCE_ENERGIES = 3.6e-2
TOLERANCE_COMPARE = 1e-10

INPUT_DIR = path.join("tests", "input")
EXPECTED_DIR = path.join("tests", "expected")

INPUT_FULL = ["finite_pot_box", "inf_pot_box", "harm_osz", "asym_harm_osz",
              "double_pot_linear", "double_pot_spline"]
INPUT_ENERGY = ["finite_pot_box", "inf_pot_box", "harm_osz", "asym_harm_osz"]

ENERGY_TESTS = ["finite_pot_box_energies", "inf_pot_box_energies",
                "harm_osz_energies", "asym_harm_osz_energies"]


@pytest.mark.parametrize("testname_interp", INPUT_FULL)
def test_interpolation(testname_interp):
    """Tests the interpolation of schrodinger_solver by comparing the
    interpolated potential with the given XY data"""
    input_files_interp = "{}.{}".format(path.join(INPUT_DIR, testname_interp), "inp")

    obtained_input = io.read_input(input_files_interp)
    xypot = obtained_input["xypot"]
    xpot = xypot[:, 0]
    ypot = xypot[:, 1]

    interpot = solver.interpolate(obtained_input)
    xinterp = interpot[:, 0]
    yinterp = interpot[:, 1]

    given_pot = np.empty([len(xpot), 2])
    given_pot[:, 0] = xpot
    given_pot[:, 1] = ypot

    delta = abs(xinterp[1]-xinterp[0])
    test_interp = True

    for ii in given_pot:
        x_given = ii[0]
        y_given = ii[1]
        # boolean array "bool_array"
        bool_array = np.abs(xinterp - x_given) <= delta / 2

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
            interval_average_right = np.sum(yinterp[bool_array_right]) / \
            len(yinterp[bool_array_right])
            interval_average_left = np.sum(yinterp[bool_array_left]) / \
            len(yinterp[bool_array_left])
            interval_average = np.sum(yinterp[bool_array]) / \
            len(yinterp[bool_array])

            intol_lr = ((np.abs(interval_average_right - y_given) < TOLERANCE_INTERP).all() \
            or (np.abs(interval_average_left - y_given) < TOLERANCE_INTERP).all())

            intol_all = (np.abs(interval_average - y_given) < TOLERANCE_INTERP).all() \
            or intol_lr
            test_interp = test_interp and intol_all
        else:
            interval_average = np.sum(yinterp[bool_array]) / len(yinterp[bool_array])
            test_interp = test_interp \
            and (np.abs(interval_average - y_given) < TOLERANCE_INTERP).all()

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

    reference_files = "{}{}.{}" \
    .format(path.join(EXPECTED_DIR, testname_energie), "_energies", "out")
    eiva_ref = np.loadtxt(reference_files)

    obtained_input = io.read_input(input_files_energies)
    interpot = solver.interpolate(obtained_input)
    data = solver.solve1d(obtained_input, interpot)
    eiva = data["energies"]

    if np.size(eiva_ref) == 1:
        eiva_restricted = eiva[0:np.size(eiva_ref)]
    else:
        eiva_restricted = eiva[0:len(eiva_ref)]

    test_energies_assert = True
    test_energies_assert = test_energies_assert \
    and (np.abs(eiva_restricted - eiva_ref) < TOLERANCE_ENERGIES).all()

    assert test_energies_assert

@pytest.mark.parametrize("testname_compare", INPUT_FULL)
def test_compare(testname_compare):
    """Tests the constant functioning of the code by comparing previously
    calculated energy-levels and potential data with the current output of the solver"""
    input_comp = "{}.{}".format(path.join(INPUT_DIR, testname_compare), "inp")

    ref_en_comp = "{}{}.{}" \
    .format(path.join(EXPECTED_DIR, testname_compare), "_energies_compare", "out")
    eiva_ref_comp = np.loadtxt(ref_en_comp)

    ref_in_comp = "{}{}.{}" \
    .format(path.join(EXPECTED_DIR, testname_compare), "_interp_compare", "out")
    interp_ref_comp = np.loadtxt(ref_in_comp)
    xref = interp_ref_comp[:, 0]
    yref = interp_ref_comp[:, 1]

    obtained_input = io.read_input(input_comp)

    interpot = solver.interpolate(obtained_input)
    xinterp = interpot[:, 0]
    yinterp = interpot[:, 1]

    data = solver.solve1d(obtained_input, interpot)
    eiva = data["energies"]
    eiva_res = eiva[0:len(eiva_ref_comp)]

    test_compare_assert = True
    test_compare_assert = test_compare_assert and (np.abs(xinterp - xref) < TOLERANCE_INTERP).all()
    test_compare_assert = test_compare_assert and (np.abs(yinterp - yref) < TOLERANCE_INTERP).all()
    test_compare_assert = test_compare_assert and (np.abs(eiva_res - eiva_ref_comp) < TOLERANCE_ENERGIES).all()

    assert test_compare_assert


if __name__ == "__main__":
    pytest.main()
