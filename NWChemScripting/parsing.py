import re
import numpy as np

def parse_roots_from_tddft_output(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        
    rootre = re.compile('Root\s+\d+')
    transitionre = re.compile(r'\s+Occ\.\s+(\d+)\s+(alpha|beta)\s+\w\s+-+\s+Virt.\s+(\d+)\s+(alpha|beta)\s+\w\s+(-?\d+\.\d+)\s+')
    
    rootlinenums = []
    for i, l in enumerate(lines):
        if rootre.search(l):
            rootlinenums.append(i)
            
    rootdict = {}
    for i in range(len(rootlinenums)):
        rootdict[i + 1] = {}
        
    for k, d in rootdict.items():
        d['transitions'] = []
        if k != len(rootdict):
            lower = rootlinenums[k-1]
            upper = rootlinenums[k]
        else:
            lower = rootlinenums[k-1]
            upper = None
        for l in lines[lower:upper]:
            m = transitionre.match(l)
            if m:
                trans = {}
                trans['occ ({})'.format(m.group(2))] = m.group(1)
                trans['virt ({})'.format(m.group(4))] = m.group(3)
                trans['coeff'] = float(m.group(5))
                d['transitions'].append(trans)

            dm = re.match(r'\s+Dipole Oscillator Strength\s+(-?\d+\.\d+)\s+', l)
            if dm:
                d['Dipole Oscillator Strength'] = float(dm.group(1))

            qm = re.match(r'\s+Electric Quadrupole\s+(-?\d+\.\d+)\s+', l)
            if qm:
                d['Electric Quadrupole'] = float(qm.group(1))

            mm = re.match(r'\s+Magnetic Dipole\s+(-?\d+\.\d+)\s+', l)
            if mm:
                d['Magnetic Dipole'] = float(mm.group(1))

            tm = re.match(r'\s+Total Oscillator Strength\s+(-?\d+\.\d+)\s+', l)
            if tm:
                d['Total Oscillator Strength'] = float(tm.group(1))

            auev = re.match(r'\s+Root\s+\d+\s+\w\s+(-?\d+\.\d+)\s+a.u.\s+(-?\d+\.\d+)\s+eV\s+', l)
            if auev:
                d['a.u.'] = auev.group(1)
                d['eV'] = auev.group(2)

            s2 = re.match(r'\s+<S2>\s+=\s+(-?\d+\.\d+)\s+', l)
            if s2:
                d['<S2>'] = s2.group(1)
                
    return rootdict


def parse_movec_occupation_energies(lines):
    vectorre = re.compile(r'\s+Vector\s+(\d+)\s+Occ=(\d+\.\d+)D(.\d+)\s+E=([\- ]?\d+.\d+)D(.\d+)')
    alphalinenums, betalinenums = parse_alpha_beta_vector_linenums(lines)
    
    avectors = []
    bvectors = []
    for a, b in zip(alphalinenums, betalinenums):
        for l in lines[a:b]:
            m = vectorre.match(l)
            if m:
                vec = {}
                vec['Vector Number'] = int(m.group(1))
                vec['Occupation'] = float(m.group(2)) * 10 ** int(m.group(3))
                vec['Energy'] = float(m.group(4)) * 10 ** int(m.group(5))
                avectors.append(vec)
        for l in lines[b:(b+b-a)]:
            m = vectorre.match(l)
            if m:
                vec = {}
                vec['Vector Number'] = int(m.group(1))
                vec['Occupation'] = float(m.group(2)) * 10 ** int(m.group(3))
                vec['Energy'] = float(m.group(4)) * 10 ** int(m.group(5))
                bvectors.append(vec)
    return avectors, bvectors


def parse_alpha_beta_vector_linenums(lines):
    alphalinere =  re.compile(r'\s+DFT Final Alpha Molecular Orbital Analysis\s*')

    alphalinenums = []
    for i, l in enumerate(lines):
        if alphalinere.search(l):
            alphalinenums.append(i)

    betalinere =  re.compile(r'\s+DFT Final Beta Molecular Orbital Analysis\s*')

    betalinenums = []
    for i, l in enumerate(lines):
        if betalinere.search(l):
            betalinenums.append(i)
    
    return alphalinenums, betalinenums


def parse_movec_info_all(lines):
    vectorre = re.compile(r'\s+Vector\s+(\d+)\s+Occ=(\d+\.\d+)D(.\d+)\s+E=([\- ]?\d+.\d+)D(.\d+)')
    mocenterre = re.compile(r'\s+MO Center=\s+(.\d+\.\d+)D(.\d+),\s+(.\d+\.\d+)D(.\d+),\s+(.\d+\.\d+)D(.\d+),\s+r\^2=\s+(\d+\.\d+)D(.\d+)')
    bfnre = re.compile(r'\s+(\d+)\s+(\-?\d+\.\d+)\s+(\d+)\s+(\w+)\s+(\w+)')
    alphalinenums, betalinenums = parse_alpha_beta_vector_linenums(lines)
    
    avectors = {}
    bvectors = {}

    for a, b in zip(alphalinenums, betalinenums):
        for l in lines[a:b]:
            m = vectorre.match(l)
            if m:
                # save previously parsed vector and instantiate new one
                try:
                    avectors[vecnum] = vec
                except UnboundLocalError:
                    pass
                
                vecnum = int(m.group(1))
                vec = {}
                vec['Occupation'] = float(m.group(2)) * 10 ** int(m.group(3))
                vec['Energy'] = float(m.group(4)) * 10 ** int(m.group(5))
                vec['Bfns'] = []

            m = mocenterre.match(l)
            if m:
                vec['MO Center'] = np.array([float(m.group(1)) * 10 ** int(m.group(2)),
                                            float(m.group(3)) * 10 ** int(m.group(4)),
                                            float(m.group(5)) * 10 ** int(m.group(6))])
                vec['r^s'] = float(m.group(7)) * 10 ** int(m.group(8))
                
            m = bfnre.findall(l)
            if m:
                for bfn in m:
                    bfndict = {}
                    bfndict['Bfn. #'] = bfn[0]
                    bfndict['Coefficient'] = bfn[1]
                    bfndict['Atom #'] = bfn[2]
                    bfndict['Atom'] = bfn[3]
                    bfndict['Atom Fn.'] = bfn[4]
                    vec['Bfns'].append(bfndict)
                continue
        
        # Assign last vector once line loop ends
        avectors[vecnum] = vec

                    
    for a, b in zip(alphalinenums, betalinenums):
        for l in lines[b:]:
            m = vectorre.match(l)
            if m:
                # save previously parsed vector and instantiate new one
                try:
                    bvectors[vecnum] = vec
                except UnboundLocalError:
                    pass
                
                vecnum = int(m.group(1))
                vec = {}
                vec['Occupation'] = float(m.group(2)) * 10 ** int(m.group(3))
                vec['Energy'] = float(m.group(4)) * 10 ** int(m.group(5))
                vec['Bfns'] = []

            m = mocenterre.match(l)
            if m:
                vec['MO Center'] = np.array([float(m.group(1)) * 10 ** int(m.group(2)),
                                            float(m.group(3)) * 10 ** int(m.group(4)),
                                            float(m.group(5)) * 10 ** int(m.group(6))])
                vec['r^s'] = float(m.group(7)) * 10 ** int(m.group(8))
                
            m = bfnre.findall(l)
            if m:
                for bfn in m:
                    bfndict = {}
                    bfndict['Bfn. #'] = bfn[0]
                    bfndict['Coefficient'] = bfn[1]
                    bfndict['Atom #'] = bfn[2]
                    bfndict['Atom'] = bfn[3]
                    bfndict['Atom Fn.'] = bfn[4]
                    vec['Bfns'].append(bfndict)
                continue
        
        # Assign last vector once line loop ends
        bvectors[vecnum] = vec
            
    return avectors, bvectors
