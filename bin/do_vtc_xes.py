#!/usr/bin/env python3

import os
import shutil
import time
import subprocess

from NWChemScripting import *


INITIALCHARGE = int(input('Initial charge? '))
INITIALMULT = int(input('Initial multiplicity? '))
COMPOUND = input('Compound name? ')


def run_geometry_optimize():
    os.chdir(ROOTDIR + '/geometryoptimize')
    copyfile('template-input.nw', 'input.nw')
    replace_text_in_file('input.nw', 'COMPOUND', COMPOUND)
    replace_text_in_file('input.nw', 'charge INPUTCHARGE', 'charge {}'.format(INITIALCHARGE))
    replace_text_in_file('input.nw', 'mult INPUTMULT', 'mult {}'.format(INITIALMULT))

    assert start_job() == 0, 'Run did NOT start succesfully.'

    wait_for_calculation_completion('output.out')


def run_gnd_state():
    os.chdir(ROOTDIR + '/gndstate')

    copyfile('template-input-gnd.nw', 'input-gnd.nw')

    optimizedxyz = find_highest_number_xyz_file(ROOTDIR + '/geometryoptimize/xyzfiles/')
    copyfile(ROOTDIR + '/geometryoptimize/xyzfiles/{}'.format(optimizedxyz), ROOTDIR + '/gndstate/{}'.format(optimizedxyz))

    centeredfile = center_xyz(optimizedxyz, 3)

    replace_text_in_file('input-gnd.nw', 'COMPOUND', COMPOUND)
    replace_text_in_file('input-gnd.nw', 'charge INPUTCHARGE', 'charge {}'.format(INITIALCHARGE))
    replace_text_in_file('input-gnd.nw', 'mult INPUTMULT', 'mult {}'.format(INITIALMULT))

    replace_text_in_file('input-gnd.nw', 'load GEOMETRYFILE', 'load {}'.format(centeredfile))

    assert start_job() == 0, 'Run did NOT start succesfully.'

    wait_for_calculation_completion('output-gnd.out')


def run_xes_calc():

    os.chdir(ROOTDIR + '/xescalc')

    copyfile('template-input-vtc.nw', 'input-vtc.nw')

    replace_text_in_file('input-vtc.nw', 'COMPOUND', COMPOUND)
    replace_text_in_file('input-vtc.nw', 'charge INPUTCHARGE', 'charge {}'.format(INITIALCHARGE + 1))
    replace_text_in_file('input-vtc.nw', 'mult INPUTMULT', 'mult {}'.format(INITIALMULT + 1))

    centeredfile = find_highest_number_xyz_file(ROOTDIR + '/geometryoptimize/xyzfiles/').split('.xyz')[0] + '_centered.xyz'

    copyfile(ROOTDIR + '/gndstate/{}'.format(centeredfile), ROOTDIR + '/xescalc/{}'.format(centeredfile))

    copyfile(ROOTDIR + '/gndstate/{}.movecs'.format(COMPOUND), ROOTDIR + '/xescalc/{}.movecs'.format(COMPOUND))

    replace_text_in_file('input-vtc.nw', 'load GEOMETRYFILE', 'load {}'.format(centeredfile))

    replace_text_in_file('input-vtc.nw', 'HIGHESTOCCUPIEDBETA', str(get_highest_occupied_beta_movec(ROOTDIR + '/gndstate/output-gnd.out')))

    assert start_job() == 0, 'Run did NOT start succesfully.'

    wait_for_calculation_completion('output-vtc.out')


if __name__ == '__main__':
    # Check for geometry file
    molexists = os.path.exists('input_xyz/{}.mol'.format(COMPOUND))
    xyzexists = os.path.exists('input_xyz/{}.xyz'.format(COMPOUND))
    assert molexists or xyzexists, 'Error, specified compound molecular geometry file does not exist'

    # Make directory structure
    print('Making new directory named {}'.format(COMPOUND))
    assert not os.path.exists(COMPOUND), 'Error, folder already exists'
    os.mkdir(COMPOUND)
    print('Copying template folder to new directory')
    shutil.copytree('../template/geometryoptimize', '{}/geometryoptimize'.format(COMPOUND))
    shutil.copytree('../template/gndstate', '{}/gndstate'.format(COMPOUND))
    shutil.copytree('../template/xescalc', '{}/xescalc'.format(COMPOUND))

    # Convert mol to xyz if needed
    if not xyzexists:
        convert_mol_to_xyz('input_xyz/{}.mol'.format(COMPOUND))

    print('Copying xyz file to new folder.')
    shutil.copy('input_xyz/{}.xyz'.format(COMPOUND), '{}/geometryoptimize/{}.xyz'.format(COMPOUND, COMPOUND))

    # Set current working directory  and ROOTDIR needed for the 'run' functions to work
    os.chdir('{}'.format(COMPOUND))
    ROOTDIR = os.getcwd()

    print('Beginning calculation sequence...')
    run_geometry_optimize()
    run_gnd_state()
    run_xes_calc()

    # Process final output file to .dat file
    os.chdir(ROOTDIR)
    with open('xescalc/output-vtc.out', 'r') as xesoutput:
        with open('xescalc/{}.dat'.format(COMPOUND), 'w') as xesdat:
            proc = subprocess.Popen(['nw_spectrum_vtc_wespecmod.py', '-x'], stdin=xesoutput, stdout=xesdat)
