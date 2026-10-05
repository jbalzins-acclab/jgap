#!/usr/bin/env python3
"""
test_potetnial_with_ase.py

Evaluates a fitted potential using the Atomic Simulation Environment (ASE):
  1. Per-config-type validation errors (Energy RMSE/MAE in meV/atom, Force RMSE in meV/Å,
     and Virial RMSE in meV/atom) against reference values in <test_xyz>.
     Can be disabled via --no-config-errors.
  2. Bulk physical properties: optimal lattice constant (a0), elastic constants (C11, C12, C44),
     bulk modulus (B), and relaxed vacancy formation energy (Evac_f) for BCC Fe.
     Optimization steps are capped at --max-steps (default: 100, as in validation).

Usage:
    python test_potetnial_with_ase.py <potential_file> <test_xyz> [--no-config-errors] [--skip-bulk] [--max-steps <N>]
"""

import os
import sys
import argparse
from collections import defaultdict
import numpy as np
import ase.io
import ase.units as units
from ase.build import bulk
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from elastic import get_elastic_tensor, get_elementary_deformations
import jgap
from jgap import JGAPCalculator


def evaluate_per_config_type_errors(test_xyz: str, calc: JGAPCalculator):
    """
    Computes and prints energy, force, and virial prediction errors grouped by config_type.
    """
    print(f"\nEvaluating per-config-type prediction errors on: {test_xyz}")
    frames = ase.io.read(test_xyz, ":")
    print(f"Read {len(frames)} frames from {test_xyz}")

    errors_by_type = defaultdict(
        lambda: {"e_diff_mev": [], "f_diff_mev": [], "v_diff_mev": [], "n_atoms": 0, "n_frames": 0}
    )

    for frame in frames:
        config_type = frame.info.get("config_type", "unspecified")

        # Extract reference values before assigning calc
        e_ref = frame.info.get("energy")
        if e_ref is None and frame.calc is not None:
            e_ref = frame.calc.results.get("energy")

        f_ref = frame.arrays.get("force")
        if f_ref is None:
            f_ref = frame.arrays.get("forces")
        if f_ref is None and frame.calc is not None:
            f_ref = frame.calc.results.get("forces")

        ref_v = frame.info.get("virial")
        if ref_v is None:
            ref_v = frame.info.get("virials")
        if ref_v is None:
            ref_v = frame.info.get("virial_fit")

        # Predict with potential
        frame.calc = calc
        e_pred = frame.get_potential_energy()
        f_pred = frame.get_forces()

        # Compute virials via JGAP potential if available
        v_pred = None
        try:
            jgap_atoms = jgap.Atoms.from_ase(frame)
            res = calc.potential.calculate_energy(jgap_atoms)
            v_pred = res.virials
        except Exception:
            pass

        n_atoms = len(frame)
        errors_by_type[config_type]["n_frames"] += 1
        errors_by_type[config_type]["n_atoms"] += n_atoms

        if e_ref is not None:
            e_diff = (e_pred - e_ref) / n_atoms * 1e3  # meV / atom
            errors_by_type[config_type]["e_diff_mev"].append(e_diff)

        if f_ref is not None:
            f_diff = (f_pred - f_ref).flatten() * 1e3  # meV / Å
            errors_by_type[config_type]["f_diff_mev"].extend(f_diff)

        if ref_v is not None and v_pred is not None:
            ref_v_arr = np.array(ref_v, dtype=float)
            if ref_v_arr.shape == (3, 3):
                ref_v_6 = np.array([
                    ref_v_arr[0, 0], ref_v_arr[0, 1], ref_v_arr[0, 2],
                    ref_v_arr[1, 1], ref_v_arr[1, 2], ref_v_arr[2, 2]
                ])
            elif ref_v_arr.size == 9:
                v_flat = ref_v_arr.flatten()
                ref_v_6 = np.array([
                    v_flat[0], v_flat[1], v_flat[2],
                    v_flat[4], v_flat[5], v_flat[8]
                ])
            elif ref_v_arr.size == 6:
                ref_v_6 = ref_v_arr.flatten()
            else:
                ref_v_6 = None

            if ref_v_6 is not None:
                pred_v_6 = np.array([v_pred.xx, v_pred.xy, v_pred.xz, v_pred.yy, v_pred.yz, v_pred.zz])
                v_diff = (pred_v_6 - ref_v_6) * 1e3 / n_atoms  # meV / atom
                errors_by_type[config_type]["v_diff_mev"].extend(v_diff)

    print("\n" + "=" * 115)
    print(
        f"{'Config Type':<24} {'Frames':>7} {'Atoms':>8} {'E RMSE (meV/at)':>16} {'E MAE (meV/at)':>16} {'F RMSE (meV/Å)':>16} {'V RMSE (meV/at)':>16}"
    )
    print("-" * 115)

    all_e_diff = []
    all_f_diff = []
    all_v_diff = []
    total_frames = 0
    total_atoms = 0

    for ctype, data in sorted(errors_by_type.items()):
        e_diffs = np.array(data["e_diff_mev"]) if data["e_diff_mev"] else np.array([])
        f_diffs = np.array(data["f_diff_mev"]) if data["f_diff_mev"] else np.array([])
        v_diffs = np.array(data["v_diff_mev"]) if data["v_diff_mev"] else np.array([])

        e_rmse_str = f"{np.sqrt(np.mean(e_diffs**2)):.3f}" if len(e_diffs) > 0 else "---"
        e_mae_str = f"{np.mean(np.abs(e_diffs)):.3f}" if len(e_diffs) > 0 else "---"
        f_rmse_str = f"{np.sqrt(np.mean(f_diffs**2)):.2f}" if len(f_diffs) > 0 else "---"
        v_rmse_str = f"{np.sqrt(np.mean(v_diffs**2)):.2f}" if len(v_diffs) > 0 else "---"

        if len(e_diffs) > 0:
            all_e_diff.extend(e_diffs)
        if len(f_diffs) > 0:
            all_f_diff.extend(f_diffs)
        if len(v_diffs) > 0:
            all_v_diff.extend(v_diffs)
        total_frames += data["n_frames"]
        total_atoms += data["n_atoms"]

        print(
            f"{ctype:<24} {data['n_frames']:>7d} {data['n_atoms']:>8d} {e_rmse_str:>16} {e_mae_str:>16} {f_rmse_str:>16} {v_rmse_str:>16}"
        )

    print("-" * 115)
    overall_e_rmse = f"{np.sqrt(np.mean(np.array(all_e_diff)**2)):.3f}" if all_e_diff else "---"
    overall_e_mae = f"{np.mean(np.abs(np.array(all_e_diff))):.3f}" if all_e_diff else "---"
    overall_f_rmse = f"{np.sqrt(np.mean(np.array(all_f_diff)**2)):.2f}" if all_f_diff else "---"
    overall_v_rmse = f"{np.sqrt(np.mean(np.array(all_v_diff)**2)):.2f}" if all_v_diff else "---"
    print(
        f"{'OVERALL':<24} {total_frames:>7d} {total_atoms:>8d} {overall_e_rmse:>16} {overall_e_mae:>16} {overall_f_rmse:>16} {overall_v_rmse:>16}"
    )
    print("=" * 115 + "\n")


def evaluate_bulk_properties(calc: JGAPCalculator, max_steps: int = 100):
    # 1. Minimize structure and find Lattice Constant
    print("Building Bulk Fe (BCC)...")
    atoms = bulk('Fe', 'bcc', a=2.8, cubic=True)
    atoms.calc = calc

    print(f"Minimizing bulk Fe to find optimal lattice constant (max {max_steps} steps)...")
    ucf = FrechetCellFilter(atoms)
    opt = BFGS(ucf, logfile=None)
    opt.run(fmax=0.005, steps=max_steps)

    a0 = atoms.cell[0, 0]
    v0 = atoms.get_volume()
    print(f"Relaxed lattice constant: a0 = {a0:.4f} Å (steps: {opt.nsteps})")
    print(f"Relaxed volume: V0 = {v0:.4f} Å^3\n")

    # 2. Compute Elastic Constants
    print("Calculating Elastic Constants...")
    systems = get_elementary_deformations(atoms, n=5, d=0.33)

    for s in systems:
        s.calc = calc
        s.get_stress()

    Cij, Bij = get_elastic_tensor(atoms, systems=systems)
    C_GPa = Cij / units.GPa

    print("\n--- Elastic Tensor (GPa) ---")
    print(C_GPa)
    if C_GPa.ndim == 1:
        c11 = float(C_GPa[0])
        c12 = float(C_GPa[1])
        c44 = float(C_GPa[2])
    else:
        c11 = float(C_GPa[0, 0])
        c12 = float(C_GPa[0, 1])
        c44 = float(C_GPa[3, 3])
    b_mod = (c11 + 2.0 * c12) / 3.0
    print(f"C11 = {c11:.1f} GPa, C12 = {c12:.1f} GPa, C44 = {c44:.1f} GPa, B = {b_mod:.1f} GPa")

    # 3. Compute Vacancy Formation Energy
    print("\nCalculating Vacancy Formation Energy...")
    supercell = atoms.repeat((6, 6, 6))
    supercell.calc = calc

    N = len(supercell)
    E_perfect = supercell.get_potential_energy()
    e_bulk = E_perfect / N

    vacancy_supercell = supercell.copy()
    del vacancy_supercell[0]
    vacancy_supercell.calc = calc

    E_vac_unrelaxed = vacancy_supercell.get_potential_energy()
    E_vf_unrelaxed = E_vac_unrelaxed - (N - 1) * e_bulk

    print(f"Minimizing vacancy supercell (max {max_steps} steps)...")
    opt_vac = BFGS(vacancy_supercell, logfile=None)
    opt_vac.run(fmax=0.005, steps=max_steps)

    E_vac_relaxed = vacancy_supercell.get_potential_energy()
    E_vf_relaxed = E_vac_relaxed - (N - 1) * e_bulk

    print(f"Supercell size: {N} atoms -> {N - 1} atoms with vacancy (steps: {opt_vac.nsteps})")
    print(f"Bulk energy per atom: {e_bulk:.4f} eV")
    print(f"Unrelaxed Vacancy Formation Energy: {E_vf_unrelaxed:.4f} eV")
    print(f"Relaxed Vacancy Formation Energy:   {E_vf_relaxed:.4f} eV\n")


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate a fitted potential with ASE: compute per-config-type errors and bulk physical properties."
    )
    parser.add_argument("potential_file", help="Path to .jgap.h5 or .tabgap.h5 potential")
    parser.add_argument("test_xyz", help="Path to test XYZ dataset containing reference energies, forces, and virials")
    parser.add_argument(
        "--no-config-errors",
        action="store_true",
        help="Disable per-config-type error summary table",
    )
    parser.add_argument(
        "--skip-bulk",
        action="store_true",
        help="Skip bulk BCC relaxation and physical property calculations",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=100,
        help="Maximum number of BFGS relaxation steps (default: 100, as in validation)",
    )

    args = parser.parse_args()

    if not os.path.exists(args.potential_file):
        print(f"Error: Potential file '{args.potential_file}' not found.")
        sys.exit(1)

    if not os.path.exists(args.test_xyz):
        print(f"Error: Test XYZ dataset '{args.test_xyz}' not found.")
        sys.exit(1)

    print(f"Using potential: {args.potential_file}")
    calc = JGAPCalculator(args.potential_file)

    # 1. Per config-type error printing (enabled by default, disabled with --no-config-errors)
    if not args.no_config_errors:
        evaluate_per_config_type_errors(args.test_xyz, calc)

    # 2. Bulk physical properties calculation
    if not args.skip_bulk:
        evaluate_bulk_properties(calc, max_steps=args.max_steps)


if __name__ == "__main__":
    main()
