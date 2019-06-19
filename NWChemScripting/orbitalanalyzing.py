import pandas as pd
import numpy as np

def get_sticks(parsed, erange=(2440, 2475), oscillatorthresh=1e-10):
    x = []
    y = []
    for r, v in parsed.items():
        energy = abs(float(v['eV'])) + COMMONSHIFT
        if energy > erange[0] and energy < erange[1]:
            if abs(v['Total Oscillator Strength']) > oscillatorthresh:
                x.append(energy)
                y.append(abs(v['Total Oscillator Strength']))
    return np.array([x, y])

def get_transitions_erange_threshold(parsed, erange=(2440, 2475), oscillatorthresh=1e-10):
    trans = []
    for r, v in parsed.items():
        energy = abs(float(v['eV'])) + COMMONSHIFT
        if energy > erange[0] and energy < erange[1]:
            if abs(v['Total Oscillator Strength']) > oscillatorthresh:
                trans.append(r)
    new = {}
    for t in trans:
        new[t] = parsed[t]
    return new


def get_transitions_orbitals_beta_only(roots, movecs, rootthresh=0.25, bfnthresh=0.1):
    transorbs = []
    for rnum, root in roots.items():
        for transition in root['transitions']:
            if abs(transition['coeff']) > rootthresh:
                try:
                    occm = int(transition['occ (beta)'])
                except KeyError:
                    continue
                m = movecs[occm]
                for bfn in m['Bfns']:
                    if abs(float(bfn['Coefficient'])) > bfnthresh:
                        bfn = bfn.copy()
                        bfn['osc str'] = root['Total Oscillator Strength']
                        bfn['MOnum'] = occm
                        bfn['root'] = rnum
                        bfn['root energy'] = abs(float(root['eV'])) + COMMONSHIFT
                        bfn['root coeff abs.'] = abs(transition['coeff'])
                        transorbs.append(bfn)
    return pd.DataFrame(transorbs)


def calc_orbitals_fractions(orbs, erange, normalizerootcoeffs=True, normalizebfncoefficients=True):
    '''Calculate fractions of atomic orbitals contributing to transitions in a specified energy range.
    
    Normalize keywords do simple normalizations of the absolute values of the coefficients.
    '''
    orbs = orbs.copy()
    orbs = orbs[(orbs['root energy'] > erange[0]) & (orbs['root energy'] < erange[1])]
    
    # normalize root coefficients for each root
    if normalizerootcoeffs:
        for r in orbs['root'].unique():
            sumrootcoeff = 0
            for mo in orbs[orbs['root']==r]['MOnum'].unique():
                sumrootcoeff += orbs[(orbs['root']==r) & (orbs['MOnum']==mo)]['root coeff abs.'].iloc[0]

            orbs.loc[orbs['root']==r, 'root coeff abs.'] /= sumrootcoeff

    # normalize coefficient of bfn's for each MO
    if normalizebfncoefficients:
        for r in orbs['root'].unique():
            for mo in orbs[orbs['root']==r]['MOnum'].unique():
                subset = (orbs['root']==r) & (orbs['MOnum']==mo)
                orbs.loc[subset, 'Coefficient'] /= orbs.loc[subset, 'Coefficient'].abs().sum()
    
    fractions = {}
    for row in orbs.values:
        try:
            fractions['{} {}'.format(row[0], row[2][0])] += abs(float(row[4])) * abs(row[8]) * abs(row[6])
        except KeyError:
            fractions['{} {}'.format(row[0], row[2][0])] = 0.
            fractions['{} {}'.format(row[0], row[2][0])] += abs(float(row[4])) * abs(row[8]) * abs(row[6])
            
    normedfracs = {}
    for k, f in fractions.items():
        normedfracs[k] = f / np.sum(list(fractions.values()))
    
    return normedfracs

def check_vector_direction(vec, direction, threshold):
    return abs(np.dot(vec, direction)) / (abs(np.linalg.norm(vec)) * abs(np.linalg.norm(direction))) >= threshold

def filter_roots_by_tm_direction(roots, direction, threshold):
    filtered = {}
    for r, v in roots.items():
        if check_vector_direction(v['Transition Moments (XYZ)'], direction, threshold):
            filtered[r] = v
    return filtered