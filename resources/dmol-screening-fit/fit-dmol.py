# Modified from works of Jesper Byggmästar
# Intended for use on data from https://zenodo.org/records/14172633 (doi:0.5281/zenodo.14172632)

import numpy as np
import argparse
import scipy
import os
import json
from scipy.optimize import curve_fit

def screened_coulomb_ev(Z1, Z2, r, c1, c2, c3, c4, c5, c6):
    eps = 8.854187817e-12
    e = 1.60217657e-19

    a = 0.46848 / (Z1**0.23 + Z2**0.23)
    x = r / a
    phi = c1 * np.exp(-c2 * x) + c3 * np.exp(-c4 * x) + c5 * np.exp(-c6 * x)
    r = r * 1e-10  # Å --> m
    E = Z1 * Z2 * e * phi / (4.0 * np.pi * eps * r) # not e**2 => eV
    return E

def read_energies(dmolfile):

    # read dmol data
    r_dmol = []
    E_dmol = []

    # inconsistent number of columns, so cannot use np.loadtxt
    with open(dmolfile, 'r') as file:
        for line in file:
            line = line.split()
            r_dmol.append(float(line[0]))
            E_dmol.append(float(line[1]))

    # subtract last energy offset
    E_dmol = [e - E_dmol[-1] for e in E_dmol]

    # remove last two points (r=100, 1000)
    while r_dmol[-1] > 10:
        del r_dmol[-1]
        del E_dmol[-1]

    # sometimes first point is E=0 for some reason, remove it
    if E_dmol[0] < 1.0:
        del r_dmol[0]
        del E_dmol[0]

    return r_dmol, E_dmol

def read_energies_all(energies_dir):

    res = {}

    for fn in os.listdir(energies_dir):
        if 'energies.' not in fn: continue
        Z1 = int(fn.split('.')[1])
        Z2 = int(fn.split('.')[2])

        res[(Z1,Z2)] = read_energies(os.path.join(energies_dir, fn))

    return res

def fit(Z1, Z2, r, energies, fit_interval):

    for i in reversed(range(len(r))):
        if r[i] < fit_interval[0] or r[i] > fit_interval[1]:
            del r[i]
            del energies[i]

    p0 = [0.2, 1.0, 0.2, 1.0, 0.2, 1.0] # coeff guesses

    sigma = range(len(r)+1, 1, -1)
    coeff, cov = curve_fit(
        f=lambda _r, c1, c2, c3, c4, c5, c6: screened_coulomb_ev(Z1, Z2, _r, c1, c2, c3, c4, c5, c6),
        xdata=r, ydata=energies, p0=p0, sigma=sigma, maxfev=20000, bounds=(0, np.inf)
        )
    #coeff_err = np.sqrt(np.diag(cov))
    return coeff

def fit_all(energies_per_pair, fit_interval):
    res = {}

    for Zs, rs_energies in energies_per_pair.items():
        res[f'{Zs[0]},{Zs[1]}'] = fit(Zs[0], Zs[1], rs_energies[0], rs_energies[1], fit_interval).tolist()
        print(f'{len(res)/len(energies_per_pair)*100:.2f}%')

    return res

if __name__ == '__main__':
    #energies_dir = input('dmol energies.Z1.Z2 dir: ')
    energies_dir =  ''#'notes/zbl/nlh_potentials_opendata/dmol/original_data'
    #fit_interval = [float(x) for x in input('fit interval (e.g. "0.04,1.2" - no spaces!): ').split(',')]
    fit_interval = [0.04, 1.2]

    energies_per_pair = read_energies_all(energies_dir)
    coeffs_per_pair = fit_all(energies_per_pair, fit_interval)
    #print(coeffs_per_pair)

    with open('dmol-fit.json', 'w') as f_out:
        json.dump(coeffs_per_pair, f_out, indent=4)
