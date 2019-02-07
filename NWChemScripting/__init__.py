import re
import os
import time
from subprocess import call
import periodictable


def replace_text_in_file(infile, oldstr, newstr):
    # Read in the file
    with open(infile, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace(oldstr, newstr)

    # Write the file out again
    with open(infile, 'w') as file:
        file.write(filedata)


def find_highest_number_xyz_file(directory):
    numbers = []
    prefixes = []
    for filename in os.listdir(directory):
        match = re.match(r'(.*?)(\d*)(.xyz)', filename)
        prefixes.append(match.group(1))
        numbers.append(match.group(2))
    assert all([p==prefixes[0] for p in prefixes]), 'Some filename prefixes do not match'
    numbers = [int(n) for n in numbers]
    return '{}{:03d}{}'.format(match.group(1), max(numbers), match.group(3))


def read_last_line(file):
#https://stackoverflow.com/questions/3346430/what-is-the-most-efficient-way-to-get-first-and-last-line-of-a-text-file
    with open(file, "rb") as f:
        f.readline()         # Read the first line.
        f.seek(-2, os.SEEK_END)      # Jump to the second last byte.
        while f.read(1) != b"\n":    # Until EOL is found...
            f.seek(-2, os.SEEK_CUR)  # ...jump back the read byte plus one more.
        last = f.readline()          # Read last line.
    return last


def check_calculation_successful(outfile):
    try:
        res = read_last_line(outfile)[:12] == b' Total times'
        return res
    except IOError:
        return False


def get_highest_occupied_beta_movec(infile): #Deprecated, kept for backwards compatibility
    with open(infile, 'r') as f:
        content = f.read()
        betaorbitalsindex = content.index('DFT Final Beta Molecular Orbital Analysis')
        betaorbitals = content[betaorbitalsindex:]
        occ0index = betaorbitals.index('Occ=0')
        f.seek(betaorbitalsindex + occ0index)
        vectorindex = betaorbitals.index('Vector', occ0index - 14, occ0index)
        f.seek(betaorbitalsindex + vectorindex)
        r = f.readline()
    return int(r.split()[1]) - 1


def get_highest_occupied_movec(infile, channel='beta'):
    if channel == 'beta':
        channel = 'Beta'
    elif channel == 'alpha':
        channel = 'Alpha'
    else:
        rais RuntimeError('Channel must be \'alpha\' or \'beta\'')
    with open(infile, 'r') as f:
        content = f.read()
        orbitalsindex = content.index('DFT Final {} Molecular Orbital Analysis'.format(Channel))
        orbitals = content[orbitalsindex:]
        occ0index = orbitals.index('Occ=0')
        f.seek(orbitalsindex + occ0index)
        vectorindex = orbitals.index('Vector', occ0index - 14, occ0index)
        f.seek(orbitalsindex + vectorindex)
        r = f.readline()
    return int(r.split()[1]) - 1


def start_job():
    return call(['msub', 'job.run'])


def wait_for_calculation_completion(outfilename, maxwait=7200):
    w = 0
    while not os.path.isfile(outfilename):
        print('waiting for job to start: {}s'.format(w), end='\r')
        w = w + 1
        time.sleep(1)

    print('output file created. blocking for maxwait {}s'.format(maxwait))
    i = 0
    while not check_calculation_successful(outfilename):
        if i < maxwait:
            time.sleep(1)
            i = i + 1
        else:
            raise RuntimeError('Output file not created within {}s timeout period.'.format(maxwait))


def convert_mol_to_xyz(infile):
    name = infile.split('.mol')[0]
    with open(infile, 'r') as file:
        lines = file.readlines()

    istart = 0
    for l in lines:
        if 'V2000' in l.split() or 'V3000' in l.split():
            break
        else:
            istart = istart + 1

    atomre = re.compile(r'\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\w+)\s+')
    atoms = []
    for l in lines[istart + 1:]:
        m = atomre.match(l)
        if m is not None:
            atoms.append(m.groups())
        else:
            break

    with open(name + '.xyz', 'w') as outfile:
        outfile.write(str(len(atoms)) + '\n')
        outfile.write('File made automatically with Python script using input mol file. (v1)\n')
        for a in atoms:
            outfile.write('{}{:>15}{:>15}{:>15}\n'.format(a[-1], a[0], a[1], a[2]))


def center_xyz(infile, targetline):
    with open(infile, 'r') as file:
        with open(infile.split('.xyz')[0] + '_centered.xyz', 'w') as newfile:
            for i in range(targetline - 1):
                file.readline()
            target = file.readline().split()
            for i in range(1, len(target)):
                target[i] = float(target[i])

            file.seek(0)

            newfile.write(file.readline())
            newfile.write(file.readline())
            while True:
                r = file.readline().split()
                if len(r) != 4:
                    break
                for i in range(1, len(r)):
                    r[i] = float(r[i]) - target[i]
                newfile.write('{}  {:>15.5f}{:>15.5f}{:>15.5f}\n'.format(*r))
            newfile.write('\n')
    return infile.split('.xyz')[0] + '_centered.xyz'


def read_xyz(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        atoms = []
        coords = []
        for l in lines[2:]: #XYZ files must have two header rows
            split = l.split()
            atoms.append(split[0])
            coords.append([float(x) for x in split[1:]])
        assert len(atoms) == len(coords), 'Something went wrong, len of atoms doesnt equal length of coords'
    return atoms, coords


def basic_multiplicity_from_atoms(atoms):
    electrons = 0
    for a in atoms:
        electrons += periodictable.__getattribute__(a).number
    print('{} electrons, which means basic multiplicity {}'.format(electrons, electrons % 2 + 1))
    return electrons % 2 + 1


def make_xyz_animation(basename, directory=None):
    if directory is None:
        directory = os.getcwd() + '/'
    
    filere = re.compile(r'{}\d\d\d.xyz'.format(basename))
    
    xyzfiles = []
    for f in os.listdir(directory):
        if filere.match(f) is not None:
            xyzfiles.append(f)
    
    if not xyzfiles:
        raise FileNotFoundError('No files matching {}###.xyz were found.'.format(basename))
        
    assert not os.path.exists(directory+'{}animation.xyz'.format(basename)), 'File already exists'
        
    with open(directory+'{}animation.xyz'.format(basename), 'w') as outfile:
        for xyz in sorted(xyzfiles):
            with open(directory+xyz, 'r') as infile:
                outfile.write(infile.read())
