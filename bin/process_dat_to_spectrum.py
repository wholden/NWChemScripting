#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


def read_dft_transitions_file(path):
    return np.loadtxt(path).T


#https://scipython.com/book/chapter-8-scipy/examples/the-voigt-profile/
def Lorentzian(x, xc, gamma):
    """ Return Lorentzian line shape at x with HWHM gamma """
    return gamma / np.pi / ((x-xc)**2 + gamma**2)


def spectrum_from_transitions(transitions, lorentz_ev=1, erange=None, numpoints=1000, peaknorm=True):
    x, y = transitions
    if erange is not None:
        good = np.logical_and(x >= erange[0], x <= erange[1])
        x, y = x[good], y[good]
        x_eval = np.linspace(erange[0], erange[1], numpoints)
    else:
        xmin = np.min(x)
        xmax = np.max(x)
        padding = (xmax - xmin) / 2
        x_eval = np.linspace(xmin - padding, xmax + padding, numpoints)
    
    spectrum = np.zeros_like(x_eval)
    for e, a in zip(x, y):
        spectrum += a * Lorentzian(x_eval, e, lorentz_ev/2)
    
    if peaknorm:
        spectrum = spectrum / np.max(spectrum)

    return np.array([x_eval, spectrum])


#Thank you stackoverflow
#https://stackoverflow.com/questions/24143320/gaussian-sum-filter-for-irregular-spaced-points
def gaussian_broaden(spectrum, width_ev=2, numpoints=1000, xmin=None, xmax=None):
    x, y = spectrum
    if xmin is None:
        xmin = np.min(x)
    if xmax is None:
        xmax = np.max(x)
    x_eval = np.linspace(xmin, xmax, numpoints)
    sigma = width_ev/(2*np.sqrt(2*np.log(2)))

    delta_x = x_eval[:, None] - x
    weights = np.exp(-delta_x*delta_x / (2*sigma*sigma)) / (np.sqrt(2*np.pi) * sigma)
    weights /= np.sum(weights, axis=1, keepdims=True)
    y_eval = np.dot(weights, y)

    return np.array([x_eval, y_eval])


def plot_spectrum_and_transitions(transitions, lorentz_ev=1, erange=None, 
                                numpoints=1000, gaussian_ev=None, show=True):
    spectrum = spectrum_from_transitions(transitions, lorentz_ev=lorentz_ev, 
                            erange=erange, numpoints=numpoints, peaknorm=False)
    x, y = spectrum
    norm = np.max(y)

    fig, ax = plt.subplots()
    ax.plot(x, y / norm)

    #rescale so that stem matches spectral height
    rescale = Lorentzian(0, 0, lorentz_ev/2)
    xs, ys = transitions
    if erange is not None:
        good = np.logical_and(xs >= erange[0], xs <= erange[1])
        xs, ys = xs[good], ys[good]
    markerline, stemlines, baseline = ax.stem(xs, ys * rescale / norm, 
                                                basefmt='k', linefmt='C0-')
    plt.setp(baseline, visible=False)
    plt.setp(stemlines, 'linewidth', 1)
    plt.setp(markerline, 'markersize', 3)
    plt.xlabel('Energy (eV)')

    if show:
        plt.show()

    return fig


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

