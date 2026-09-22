# Modified from works of Jesper Byggmästar

import numpy as np
import os
import urllib.request
import tarfile
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

HARTREE_TO_EV = 27.211386245988

ELEMENT_SYMBOLS = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca",
    21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr",
    41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd",
    61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb",
    71: "Lu", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg",
    81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
    91: "Pa", 92: "U"
}
SYMBOL_TO_Z = {v: k for k, v in ELEMENT_SYMBOLS.items()}

def screened_coulomb_ev(Z1, Z2, r, c1, c2, c3, c4, c5, c6):
    eps = 8.854187817e-12
    e = 1.60217657e-19

    phi = c1 * np.exp(-c2 * r) + c3 * np.exp(-c4 * r) + c5 * np.exp(-c6 * r)
    r_m = r * 1e-10  # Å --> m
    E = Z1 * Z2 * e * phi / (4.0 * np.pi * eps * r_m) # not e**2 => eV
    return E

def read_and_sort_energies(filepath, energy_col, hartree=False):
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            cols = line.split()
            try:
                dist = float(cols[0])
                energy = float(cols[energy_col])
                if hartree:
                    energy *= HARTREE_TO_EV
                data.append((dist, energy))
            except (ValueError, IndexError):
                continue # Skip malformed lines
    
    data.sort() # Sort by distance
    r = [item[0] for item in data]
    E = [item[1] for item in data]
    return r, E

def read_energies_dmol(filepath):
    r, E = read_and_sort_energies(filepath, 1)
    if not E: return [], []
    
    E = [e - E[-1] for e in E]

    while len(r) > 0 and r[-1] > 10:
        del r[-1]
        del E[-1]

    if len(r) > 0 and E[0] < 1.0:
        del r[0]
        del E[0]

    return r, E

def read_energies_mp2(filepath):
    return read_and_sort_energies(filepath, 2, hartree=True)

def read_energies_all(base_dir):
    res = {}

    # DMol
    dmol_dir = os.path.join(base_dir, 'nlh_potentials_opendata/dmol/original_data')
    if os.path.isdir(dmol_dir):
        for fn in os.listdir(dmol_dir):
            if not fn.startswith('energies.'): continue
            parts = fn.split('.')
            Z1, Z2 = int(parts[1]), int(parts[2])
            res[('dmol', Z1, Z2)] = read_energies_dmol(os.path.join(dmol_dir, fn))

    # MP2
    mp2_dir = os.path.join(base_dir, 'nlh_potentials_opendata/mp2/original_data')
    if os.path.isdir(mp2_dir):
        for fn in os.listdir(mp2_dir):
            if not fn.endswith('_unc-pc-2.dat'): continue
            parts = fn.replace('_unc-pc-2.dat', '').split('_')
            if len(parts) != 2: continue
            try:
                Z1, Z2 = SYMBOL_TO_Z[parts[0]], SYMBOL_TO_Z[parts[1]]
                res[('mp2', Z1, Z2)] = read_energies_mp2(os.path.join(mp2_dir, fn))
            except KeyError:
                print(f"Skipping {fn}, unknown element symbol.")

    return res

def fit(Z1, Z2, r, energies, fit_interval):
    r_fit, energies_fit = [], []
    for i in range(len(r)):
        if fit_interval[0] <= r[i] <= fit_interval[1]:
            r_fit.append(r[i])
            energies_fit.append(energies[i])
    
    if len(r_fit) < 6: return None

    p0 = [0.2, 1.0, 0.2, 1.0, 0.2, 1.0]
    sigma = range(len(r_fit) + 1, 1, -1)
    try:
        coeff, cov = curve_fit(
            f=lambda _r, c1, c2, c3, c4, c5, c6: screened_coulomb_ev(Z1, Z2, _r, c1, c2, c3, c4, c5, c6),
            xdata=r_fit, ydata=energies_fit, p0=p0, sigma=sigma, maxfev=20000, bounds=(0, np.inf)
        )
        return coeff
    except Exception as e:
        print(f"Fit failed for Z1={Z1}, Z2={Z2} with error: {e}")
        return None

def fit_all_and_save(energies_per_pair, fit_interval):
    results = {'dmol': [], 'mp2': []}

    total = len(energies_per_pair)
    count = 0
    last_printed_percent = -1
    for (source, Z1, Z2), (r, E) in energies_per_pair.items():
        count += 1
        coeffs = fit(Z1, Z2, r, E, fit_interval)
        if coeffs is not None:
            elem1 = ELEMENT_SYMBOLS.get(Z1, f"Z{Z1}")
            elem2 = ELEMENT_SYMBOLS.get(Z2, f"Z{Z2}")
            results[source].append(f"{elem1} {elem2} " + " ".join(map(str, coeffs)))
            
        current_percent = int(count / total * 100)
        if current_percent % 5 == 0 and current_percent != last_printed_percent:
            print(f'{current_percent}%')
            last_printed_percent = current_percent

    for source, lines in results.items():
        with open(f'{source}.dat', 'w') as f:
            f.write("\n".join(lines))

def download_nlh_data():
    url = "https://zenodo.org/records/17302337/files/nlh_potentials_opendata.tar.gz?download=1"
    filename = "nlh_potentials_opendata.tar.gz"
    if not os.path.exists(filename):
        print(f"Downloading {filename}...")
        urllib.request.urlretrieve(url, filename)
        print("Download complete.")
    
    print(f"Extracting {filename}...")
    with tarfile.open(filename, "r:gz") as tar:
        tar.extractall()
    print("Extraction complete.")

def plot_fit(elem1, elem2, fit_file, original_data_path):
    Z1, Z2 = SYMBOL_TO_Z[elem1], SYMBOL_TO_Z[elem2]
    
    coeffs = None
    with open(fit_file, 'r') as f:
        for line in f:
            parts = line.split()
            if (parts[0] == elem1 and parts[1] == elem2) or \
               (parts[0] == elem2 and parts[1] == elem1):
                coeffs = [float(c) for c in parts[2:]]
                break
    
    if coeffs is None:
        print(f"Coefficients for {elem1}-{elem2} not found in {fit_file}")
        return

    r_orig, E_orig = [], []
    if "dmol" in original_data_path:
        r_orig, E_orig = read_energies_dmol(original_data_path)
    elif "mp2" in original_data_path:
        r_orig, E_orig = read_energies_mp2(original_data_path)
    else:
        print("Could not determine original data type.")
        return

    r_dense = np.linspace(min(r_orig) if r_orig else 0, max(r_orig) if r_orig else 3, 1000)
    E_fit = screened_coulomb_ev(Z1, Z2, r_dense, *coeffs)

    # Log plot
    plt.figure()
    plt.plot(r_orig, E_orig, 'o--', label='Original Data')
    plt.plot(r_dense, E_fit, label='Fit')
    plt.yscale('log')
    plt.ylim([1, 1e8])
    plt.xlim([0, 3])
    plt.grid()
    plt.legend()
    plt.title(f'Log Scale Fit: {elem1}-{elem2}')

    # Linear plot
    plt.figure()
    plt.plot(r_orig, E_orig, 'o--', label='Original Data')
    plt.plot(r_dense, E_fit, label='Fit')
    plt.ylim([-10, 1000])
    plt.xlim([0, 5])
    plt.grid()
    plt.legend()
    plt.title(f'Linear Scale Fit: {elem1}-{elem2}')
    
    plt.show()

def main():
    download_nlh_data()

    base_dir = '.'
    fit_interval = [0.04, 1.2]

    energies_per_pair = read_energies_all(base_dir)
    fit_all_and_save(energies_per_pair, fit_interval)
    print("Fitting complete. Results saved to dmol.dat and mp2.dat")

if __name__ == '__main__':
    # To run fitting:
    main()

    # Testing specific energy
    c_fe_ni = [
        0.32818750683529224, 1.2763590338471473,
        0.5319935400634741,  0.4920242580140316,
        0.13562223656825356, 2.9359601944688674
        ]

    print(screened_coulomb_ev(26, 28, 0.5, *c_fe_ni))

    # To plot a fit:
    #plot_fit('Fe', 'Fe', 'dmol.dat', 'nlh_potentials_opendata/dmol/original_data/energies.26.26')
    #plot_fit('Al', 'Al', 'mp2.dat', 'nlh_potentials_opendata/mp2/original_data/Al_Al_unc-pc-2.dat')
    pass
