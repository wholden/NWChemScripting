#!/usr/bin/env python

import argparse
from NWChemScripting import center_xyz

parser = argparse.ArgumentParser(description='Script to center xyz file on a given atom')

parser.add_argument('infile', action='store')
parser.add_argument('targetline', action='store', type=int)

args = parser.parse_args()

infile = args.infile
targetline = args.targetline

assert infile[-3:] == 'xyz', 'File must be xyz file with .xyz extension'

center_xyz(infile, targetline)
