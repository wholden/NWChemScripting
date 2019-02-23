#!/usr/bin/env python

import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Script to generate bq charges for hydrogens in a cluster based on nearest neighbor identity.\nexample:\n make_bq_charges.py V2O3FinalCentered.xyz --atom_bqs O=-1.5,V=-0.5')

class StoreDictKeyPair(argparse.Action):
     def __call__(self, parser, namespace, values, option_string=None):
            my_dict = {}
            for kv in values.split(","):
                k,v = kv.split("=")
                my_dict[k] = v
            setattr(namespace, self.dest, my_dict)

parser.add_argument('infile', action='store')
parser.add_argument("--atom_bqs", dest="bqcharges", action=StoreDictKeyPair, metavar="KEY1=VAL1,KEY2=VAL2...")

args = parser.parse_args()

infile = args.infile
bqcharges = args.bqcharges

for k, v in bqcharges.items():
    bqcharges[k] = float(v)

assert infile[-3:] == 'xyz', 'File must be xyz file with .xyz extension'

def read_xyz(infile):
    atoms = []
    foundatoms = False

    with open(infile, 'r') as file:

        while not foundatoms:
            r = file.readline().split()
            if len(r) == 4:
                foundatoms=True
                for i in range(1, len(r)):
                    r[i] = float(r[i])
                atoms.append(r)
                break

        while True:
            r = file.readline().split()
            if len(r) != 4:
                break
            for i in range(1, len(r)):
                r[i] = float(r[i])
            atoms.append(r)
    return atoms

def get_nearest_atom_of_hydrogens(atomscoords):
    atomscoords = np.array(atomscoords)
    atoms = atomscoords[:, 0]
    coords = atomscoords[:, 1:].astype('float')
    hydrogens = np.where(atoms=='H')[0]
    nearests = []
    
    for h in hydrogens:
        dists = np.array([])
        for c in coords:
            dists = np.append(dists, np.linalg.norm(c - coords[h]))
        near = np.where(dists==np.min(dists[np.where(dists!=0)[0]]))[0][0]
        nearests.append(atoms[near])
    nearests = np.array(nearests)
    
#     return np.append(coords[hydrogens].T, nearests[np.newaxis, :], axis=0).T
    return coords[hydrogens], nearests

def write_bq_file(filename, nearestatoms, bqchargemap):
    with open(filename, 'w') as file:
#         for na in nearestatoms:
#             file.write('Bq {:>15.8f}{:>15.8f}{:>15.8f}{:>10.4f}\n'.format(
#                 float(na[0]), float(na[1]), float(na[2]), bqchargemap[na[3]]
#             ))
        for c, a in zip(*nearestatoms):
            file.write('Bq {:>15.8f}{:>15.8f}{:>15.8f}{:>10.4f}\n'.format(
                c[0], c[1], c[2], bqchargemap[a]
            ))

def make_bq_file_from_carved_cluster(outfile, infile, bqchargemap):
    write_bq_file(outfile, get_nearest_atom_of_hydrogens(read_xyz(infile)), bqchargemap)

def count_bq_charge(infile):
    with open(infile, 'r') as file:
        bqs = []
        while True:
            try:
                bqs.append(float(file.readline().split()[-1]))
            except IndexError:
                break
    return np.sum(bqs)
    
    
if __name__ == '__main__':
    make_bq_file_from_carved_cluster(infile[:-4]+'_bqcharges', infile, bqcharges)
    print('Total bq charge: {}'.format(count_bq_charge(infile[:-4]+'_bqcharges')))
