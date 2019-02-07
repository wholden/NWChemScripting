#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from NWChemScripting import *


if __name__ == '__main__':
    import argparse
    import os.path

    parser = argparse.ArgumentParser(description=('Script to generate spectra '
        'and/or plots from transitions calculated by NWChem.'))

    parser.add_argument('-f', action='store', dest='filename', type=str,
        help='File of transitions to generate spectrum from.')
    parser.add_argument('-l', action='store', dest='lorentz_width', 
        type=float, default=1.0,
        help='Set the Lorentz width of the transitions (eV).')
    parser.add_argument('-g', action='store', dest='gaussian_broadening', 
        type=float,
        help='Set Gaussian instrumental broadening (eV).')
    parser.add_argument('-p', action='store_true', default=False, dest='plot', 
        help='Make a png plot of the spectrum and transitions.')
    parser.add_argument('--lowerE', action='store', dest='lowerenergy', 
        type=float,
        help='Set lower bound of Erange.')
    parser.add_argument('--upperE', action='store', dest='upperenergy', 
        type=float,
        help='Set upper bound of Erange.')
    parser.add_argument('-s', action='store_true', default=False, dest='dosave', 
        help='Save the resulting processed spectrum.')
 
    result = parser.parse_args()

    erange = None
    if result.lowerenergy is not None and result.upperenergy is not None:
        erange = [result.lowerenergy, result.upperenergy]

    transitions = read_dft_transitions_file(result.filename)
    spectrum = spectrum_from_transitions(transitions, 
        lorentz_ev=result.lorentz_width, peaknorm=False, erange=erange)
    x, y = spectrum
    norm = np.max(y)
    spectrum = np.array([x, y / norm])  

    if result.gaussian_broadening is not None:
        spectrum = gaussian_broaden(spectrum, width_ev=result.gaussian_broadening)

    outputfile = result.filename.split('/')[-1].split('.')[0]
    if outputfile == '':
        outputfile = result.filename.split('\\')[-1].split('.')[0]
    if outputfile == '':
        print('WARNING: Filename could not be parsed. Saving as outspectrum.dat')
        outputfile = 'outspectrum'

    assert not os.path.exists(outputfile+'.processedspectrum'), 'File with outputname already exists!'
    assert not os.path.exists(outputfile+'.png'), 'File with outputname already exists!'

    if result.dosave:
        np.savetxt(outputfile+'.processedspectrum', spectrum.T)

    if result.plot:
        fig, ax = plt.subplots(dpi=300)
        ax.plot(*spectrum)

        #rescale so that stem matches spectral height
        rescale = Lorentzian(0, 0, result.lorentz_width/2)
        xs, ys = transitions
        if erange is not None:
            good = np.logical_and(xs >= erange[0], xs <= erange[1])
            xs, ys = xs[good], ys[good]
        markerline, stemlines, baseline = ax.stem(xs, ys * rescale / norm, 
                                                    basefmt='k', linefmt='b-')
        plt.setp(baseline, visible=False)
        plt.setp(stemlines, 'linewidth', 1)
        plt.setp(markerline, 'markersize', 3)

        plt.xlabel('Energy (eV)')

        plt.show()
        # fig.savefig(outputfile+'.png')

