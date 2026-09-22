import os
import numpy as np
import ase.units as units
from ase.build import bulk
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from elastic import get_elastic_tensor, get_elementary_deformations
from jgap import JGAPCalculator

potential_file = "../coordination_fit/fe.jgap.h5"
if not os.path.exists(potential_file):
    print(f"Error: {potential_file} not found.")
    print("Please run the CoordinationFit example to generate this file first.")
    exit(1)

# 1. Minimize structure and find Lattice Constant
print("Building Bulk Fe (BCC)...")
atoms = bulk('Fe', 'bcc', a=2.8, cubic=True)
calc = JGAPCalculator([potential_file])
atoms.calc = calc

print("Minimizing bulk Fe to find the optimal lattice constant...")
ucf = FrechetCellFilter(atoms)
opt = BFGS(ucf, logfile=None)
opt.run(fmax=0.001)

a0 = atoms.cell[0, 0]
v0 = atoms.get_volume()
print(f"Relaxed lattice constant: a0 = {a0:.4f} Å")
print(f"Relaxed volume: V0 = {v0:.4f} Å^3\n")

# 2. Compute Elastic Constants
print("Calculating Elastic Constants...")
systems = get_elementary_deformations(atoms, n=5, d=0.33)

# Calculate stress for all deformed systems
for s in systems:
    s.calc = calc
    s.get_stress()

Cij, Bij = get_elastic_tensor(atoms, systems=systems)

C_GPa = Cij / units.GPa

print("\n--- Elastic Tensor (GPa) ---")
print(C_GPa)

# 3. Compute Vacancy Formation Energy
print("\nCalculating Vacancy Formation Energy...")
# Repeat minimized structure 10x in each direction
supercell = atoms.repeat((6, 6, 6))
supercell.calc = calc

N = len(supercell)
E_perfect = supercell.get_potential_energy()
e_bulk = E_perfect / N

# Create a vacancy by removing one atom
vacancy_supercell = supercell.copy()
del vacancy_supercell[0]
vacancy_supercell.calc = calc

E_vac_unrelaxed = vacancy_supercell.get_potential_energy()
E_vf_unrelaxed = E_vac_unrelaxed - (N - 1) * e_bulk

print("Minimizing vacancy supercell...")
opt_vac = BFGS(vacancy_supercell, logfile=None)
opt_vac.run(fmax=0.001)

E_vac_relaxed = vacancy_supercell.get_potential_energy()
E_vf_relaxed = E_vac_relaxed - (N - 1) * e_bulk

print(f"Supercell size: {N} atoms -> {N - 1} atoms with vacancy")
print(f"Bulk energy per atom: {e_bulk:.4f} eV")
print(f"Unrelaxed Vacancy Formation Energy: {E_vf_unrelaxed:.4f} eV")
print(f"Relaxed Vacancy Formation Energy:   {E_vf_relaxed:.4f} eV")
