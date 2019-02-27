import numpy as np
import itertools
from scipy.optimize import minimize
from NWChemScripting import *

def make_constrain_ordered_xyz(file, neighbor0='P', neighbor1='H'):
    mol = mg.Molecule.from_file(file)
    assign_nearest_neighbors(mol)

    edge = mg.Molecule(['Zr'], [[0,0,0]]) #dummy Zr to let me instantiate the molecule
    edge.pop(0) #discard dummy

    poplist = []
    for i, site in enumerate(mol):
        if site.specie.name == neighbor0:
            for n in site.nearests:
                if n.specie.name == neighbor1:
                    poplist.append(i)
                    edge.append(site.specie, site.coords)
                    continue
    mol.remove_sites(poplist)

    poplist = []
    for i, site in enumerate(mol):
        if site.specie.name == neighbor1:
            for n in site.nearests:
                if n.specie.name == neighbor0:
                    poplist.append(i)
                    edge.append(site.specie, site.coords)
                    continue
    mol.remove_sites(poplist)

    filename = file.split('.xyz')[0].split('/')[-1] #oops only works on linux...
#     unconstrainf = '{}_Unconstrain.xyz'.format(filename)
#     mol.to('xyz', unconstrainf)
#     constrainf = '{}_Constrain.xyz'.format(filename)
#     edge.to('xyz', constrainf)

#     unca, uncc = read_xyz(unconstrainf)
#     ca, cc = read_xyz(constrainf)

    unca, uncc = mol.species, mol.cart_coords
    ca, cc = edge.species, edge.cart_coords

    print('Saving new file as: {}'.format('{}_Reordered'.format(filename)))
    write_xyz_from_atoms_coords('{}_Reordered'.format(filename), 
                                atoms=unca + ca, 
                                coords=np.concatenate([uncc, cc]), 
                                comment='Unconstrained {}; Constrained {}; Fix atoms {}:{}'.format(
                                    len(unca),
                                    len(ca),
                                    len(unca) + 1,
                                    len(unca) + len(ca)
                                ))


def get_neighbor_site_angles(mol, site, r=2.5):
    neighbors = [n[0] for n in mol.get_neighbors(site, r)]
    angles = []
    for comb in itertools.combinations(neighbors, 2):
        j = mol.index(site)
        i, k = [mol.index(s) for s in comb]
        angles.append(mol.get_angle(i, j, k))
    return sorted(np.array(angles))


def adjust_bond_length(site, neighbor, targetdist):
    directionvec = site.coords - neighbor.coords
    directionvec /= np.linalg.norm(directionvec)
    oldcoords = site.coords
    def func(c):
        site._coords = oldcoords + c * directionvec
        return (site.distance(neighbor) - targetdist) ** 2
    minimize(func, 0.1)


def min_sum_squared(vec0, vec1):
    sumsqs = []
    for p in itertools.permutations(vec1):
        sumsqs.append(sum((vec0 - np.array(p))**2))
    return np.min(sumsqs)


def avg_sum_squared(vec0, vec1):
    sumsqs = []
    for p in itertools.permutations(vec1):
        sumsqs.append(sum((vec0 - np.array(p))**2))
    return np.mean(sumsqs)


def total_sum_squared(vec0, vec1):
    sumsqs = []
    for p in itertools.permutations(vec1):
        sumsqs.append(sum((vec0 - np.array(p))**2))
    return np.sum(sumsqs)


def adjust_site_to_target_angles(mol, adjustsite, neighborsite, targetangles, r=2.5, minfunc=avg_sum_squared):
    def func(coords):
        adjustsite._coords = coords
        return minfunc(targetangles, get_neighbor_site_angles(mol, neighborsite, r))
    minimize(func, adjustsite._coords)
    return adjustsite._coords
