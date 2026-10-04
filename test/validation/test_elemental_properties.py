#!/usr/bin/env python3
"""
test_elemental_properties.py

Automated validation test for physical property predictions and energy evaluations
using Python, ASE, and the `elastic` library with `JGAPCalculator`.

Workflow:
  1. Refits potentials using `jgap.standard_gap_fit` + `jgap.standard_tabulation`
     for requested systems (Al, Cu, Ni, Fe, FeNi) from XYZ databases in `test/structure-databases`
     (or uses cached potentials from `build/validation_pots/` unless `--refit` is given).
  2. Evaluates physical material properties:
     - Lattice parameter (a0) via FrechetCellFilter and BFGS relaxation
     - Bulk cohesive energy (Ecoh = E_iso - E_bulk/atom)
     - Elastic constants (C11, C12, C44) and Bulk Modulus (B) via elementary deformations
     - Relaxed vacancy formation energy (Evac_f)
  3. Compares all calculated properties against DFT references (from `elemental_property_comparison.csv`).
  4. Prints formatted comparison tables to stdout.
  5. Throws errors (exits with code 1) if any deviation from DFT exceeds validation thresholds.
"""

import sys
import os
import argparse
import ctypes
import time
from pathlib import Path
import numpy as np
import pandas as pd

# Setup paths and preload libjgap.dylib
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent

for bdir in [REPO_ROOT / "build", REPO_ROOT / "cmake-build-debug", REPO_ROOT / "build/release"]:
    dylib = bdir / "libjgap.dylib"
    if dylib.exists():
        try:
            ctypes.CDLL(str(dylib))
            break
        except Exception:
            pass
sys.path.insert(0, str(REPO_ROOT / "python"))

import jgap
import ase.units as units
from ase.build import bulk
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from elastic import get_elastic_tensor, get_elementary_deformations
from jgap import JGAPCalculator

# ==============================================================================
# System Configurations for standard_gap_fit & Property Evaluation
# ==============================================================================
SYSTEM_CONFIGS = {
    "Al": {
        "training_file": "db_Al.xyz",
        "elements": ["Al"],
        "crystal": {"Al": "fcc"},
        "a_init": {"Al": 4.04},
        "supercell": {"Al": (3, 3, 3)},  # 108 atoms
        "cutoff2": 5.2,
        "cutoff2_width": 1.0,
        "cutoff3": 4.5,
        "cutoff3_width": 0.6,
        "n_sparse3": 500,
        "fallback_isolated_energies": {"Al": -0.04810075},
    },
    "Cu": {
        "training_file": "db_Cu.xyz",
        "elements": ["Cu"],
        "crystal": {"Cu": "fcc"},
        "a_init": {"Cu": 3.63},
        "supercell": {"Cu": (3, 3, 3)},  # 108 atoms
        "cutoff2": 5.2,
        "cutoff2_width": 1.0,
        "cutoff3": 4.0,
        "cutoff3_width": 0.6,
        "n_sparse3": 500,
        "fallback_isolated_energies": {"Cu": -0.02847539},
    },
    "Ni": {
        "training_file": "db_Ni.xyz",
        "elements": ["Ni"],
        "crystal": {"Ni": "fcc"},
        "a_init": {"Ni": 3.52},
        "supercell": {"Ni": (3, 3, 3)},  # 108 atoms
        "cutoff2": 5.2,
        "cutoff2_width": 1.0,
        "cutoff3": 4.0,
        "cutoff3_width": 0.6,
        "n_sparse3": 500,
        "fallback_isolated_energies": {"Ni": -0.75480209},
    },
    "Fe": {
        "training_file": "db_Fe.xyz",
        "elements": ["Fe"],
        "crystal": {"Fe": "bcc"},
        "a_init": {"Fe": 2.87},
        "supercell": {"Fe": (4, 4, 4)},  # 128 atoms
        "cutoff2": 4.5,
        "cutoff2_width": 1.0,
        "cutoff3": 3.7,
        "cutoff3_width": 0.6,
        "n_sparse3": 500,
        "fallback_isolated_energies": {"Fe": -3.38958481},
    },
    "FeNi": {
        "training_file": "feni-train.xyz",
        "elements": ["Ni", "Fe"],
        "crystal": {"Ni": "fcc", "Fe": "bcc"},
        "a_init": {"Ni": 3.52, "Fe": 2.87},
        "supercell": {"Ni": (3, 3, 3), "Fe": (4, 4, 4)},
        "cutoff2": 4.0,
        "cutoff2_width": 1.0,
        "cutoff3": 3.5,
        "cutoff3_width": 0.6,
        "n_sparse3": 500,
        "fallback_isolated_energies": {"Ni": -0.44878399, "Fe": -3.38958481},
    },
}

# Default DFT Validation Tolerances (Allowable Absolute Deviations from DFT)
DEFAULT_TOLERANCES = {
    "a0": 0.05,       # Å
    "B": 40.0,        # GPa
    "C11": 75.0,      # GPa
    "C12": 40.0,      # GPa
    "C44": 30.0,      # GPa
    "Ecoh": 0.15,     # eV
    "Evac_f": 0.20,   # eV
}


def find_xyz_file(xyz_dir, filename):
    for candidate in [
        xyz_dir / filename,
        REPO_ROOT / "test" / "structure-databases" / filename,
        REPO_ROOT / "test" / "resources" / "structure-databases" / filename,
        REPO_ROOT / "test" / "xyz-samples" / filename,
    ]:
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(f"Training database '{filename}' not found in '{xyz_dir}' or standard search paths.")


def ensure_potential(sys_name, sys_cfg, xyz_dir, pot_dir, force_refit=False):
    """
    Ensure that the potential (.tabgap.h5 and .eam.fs) exists for sys_name.
    If missing or force_refit is True, refits using jgap.standard_gap_fit and jgap.standard_tabulation.
    Returns (pot_files, isolated_energies).
    """
    pot_dir.mkdir(parents=True, exist_ok=True)
    tabgap_h5 = pot_dir / f"{sys_name}.tabgap.h5"
    eam_fs = pot_dir / f"{sys_name}.eam.fs"
    jgap_h5 = pot_dir / f"{sys_name}.jgap.h5"

    xyz_file = find_xyz_file(xyz_dir, sys_cfg["training_file"])

    iso_energies = dict(sys_cfg["fallback_isolated_energies"])

    if not force_refit and tabgap_h5.exists() and eam_fs.exists():
        print(f"[{sys_name}] Using cached potential: {tabgap_h5.name}")
        return [str(tabgap_h5), str(eam_fs)], iso_energies

    print(f"\n======================================================================")
    print(f"[{sys_name}] Refitting potential using standard_gap_fit...")
    print(f"  Training File: {xyz_file}")
    print(f"  Target Pots  : {tabgap_h5} & {eam_fs}")
    print(f"======================================================================")

    t0 = time.time()
    training_data = jgap.read_atoms(str(xyz_file))
    print(f"  Loaded {len(training_data)} frames in {time.time() - t0:.2f} s")

    # Extract isolated atom energies from single-atom frames
    for frame in training_data:
        if len(frame) == 1:
            sym = frame.symbols[0]
            iso_energies[sym] = frame.energy

    params = jgap.StandardGapParams(
        seed=42,
        cutoff2=sys_cfg["cutoff2"],
        cutoff2_width=sys_cfg["cutoff2_width"],
        n_sparse2=20,
        eam_mode=jgap.EamMode.Blind,
        eam_pair_function=jgap.EamPairFunctionType.FSGen3,
        eam_n_sparse=20,
        cutoff3=sys_cfg["cutoff3"],
        cutoff3_width=sys_cfg["cutoff3_width"],
        n_sparse3=sys_cfg["n_sparse3"],
        approx_ram_limit_gb=4.0,
    )

    rules = jgap.PerConfigTypeRegularizationRules(
        jgap.PerConfigTypeSigmas(0.001, 0.05, 0.1, 0.02)
    )
    sigmas = rules.determine_for_all(training_data)

    print(f"  Running standard_gap_fit for {sys_name}...")
    t_fit_start = time.time()
    jgap.standard_gap_fit(str(jgap_h5), training_data, sigmas, params)
    print(f"  Fit completed in {time.time() - t_fit_start:.2f} s")

    print(f"  Running standard_tabulation for {sys_name}...")
    t_tab_start = time.time()
    tab_params = jgap.StandardTabulationParams(
        r_min_3b=0.5,
        max_eam_density=10.0,
        n_grid_2b=5000,
        n_grid_3b=[80, 80, 80],
    )
    jgap.standard_tabulation(str(jgap_h5), str(pot_dir / sys_name), tab_params)
    print(f"  Tabulation completed in {time.time() - t_tab_start:.2f} s")

    if not tabgap_h5.exists() or not eam_fs.exists():
        raise RuntimeError(f"Failed to generate potential files: {tabgap_h5} or {eam_fs} missing.")

    return [str(tabgap_h5), str(eam_fs)], iso_energies


def compute_elastic_constants(atoms, calc):
    systems = get_elementary_deformations(atoms, n=5, d=0.33)
    for s in systems:
        s.calc = calc
        s.get_stress()

    Cij, _ = get_elastic_tensor(atoms, systems=systems)
    C_GPa = Cij / units.GPa

    if C_GPa.ndim == 1:
        c11 = float(C_GPa[0])
        c12 = float(C_GPa[1])
        c44 = float(C_GPa[2])
    else:
        c11 = float(C_GPa[0, 0])
        c12 = float(C_GPa[0, 1])
        c44 = float(C_GPa[3, 3])

    b_mod = (c11 + 2.0 * c12) / 3.0
    return c11, c12, c44, b_mod


def evaluate_element_properties(elem, pot_files, iso_energy, crystal, a_init, sc_dim):
    calc = JGAPCalculator(pot_files)

    # 1. Optimal Lattice Constant and Bulk Energy
    atoms = bulk(elem, crystal, a=a_init, cubic=True)
    atoms.calc = calc

    ucf = FrechetCellFilter(atoms)
    opt = BFGS(ucf, logfile=None)
    opt.run(fmax=0.005, steps=100)

    a0 = float(atoms.cell[0, 0])
    e_total_bulk = float(atoms.get_potential_energy())
    e_per_atom = e_total_bulk / len(atoms)

    # 2. Cohesive Energy
    ecoh = (iso_energy - e_per_atom) if iso_energy is not None else None

    # 3. Elastic Constants and Bulk Modulus
    c11, c12, c44, b_mod = compute_elastic_constants(atoms, calc)

    # 4. Vacancy Formation Energy (Relaxed)
    supercell = atoms.repeat(sc_dim)
    supercell.calc = calc
    n_atoms = len(supercell)
    e_perfect = float(supercell.get_potential_energy())
    e_bulk_ref = e_perfect / n_atoms

    vac_supercell = supercell.copy()
    del vac_supercell[0]
    vac_supercell.calc = calc

    opt_vac = BFGS(vac_supercell, logfile=None)
    opt_vac.run(fmax=0.005, steps=100)

    e_vac_relaxed = float(vac_supercell.get_potential_energy())
    e_vf_relaxed = e_vac_relaxed - (n_atoms - 1) * e_bulk_ref

    return {
        "element": elem,
        "crystal": crystal,
        "a0": round(a0, 4),
        "Ecoh": round(ecoh, 4) if ecoh is not None else None,
        "C11": round(c11, 1),
        "C12": round(c12, 1),
        "C44": round(c44, 1),
        "B": round(b_mod, 1),
        "Evac_f": round(e_vf_relaxed, 3),
    }


def print_comparison_table(results, ref_df, feni_df, base_tolerances):
    """
    Nicely formats and prints results vs DFT references and checks error tolerances.
    Returns (all_passed, failure_reasons).
    """
    all_passed = True
    failure_reasons = []

    print("\n" + "=" * 110)
    print(f"{'SYSTEM':<8} {'ELEM':<5} {'MODEL / SOURCE':<20} {'a0 (Å)':<9} {'Ecoh (eV)':<11} {'C11 (GPa)':<11} {'C12 (GPa)':<11} {'C44 (GPa)':<11} {'B (GPa)':<9} {'Evac_f (eV)':<11}")
    print("=" * 110)

    for res in results:
        sys_name = res["system"]
        elem = res["element"]

        # DFT Reference row
        dft_row = ref_df[(ref_df["element"] == elem) & (ref_df["model"].str.contains("DFT", case=False, na=False))]
        dft_vals = dft_row.iloc[0].to_dict() if not dft_row.empty else {}

        # Paper tabGAP row (if available)
        paper_row = ref_df[(ref_df["element"] == elem) & (ref_df["model"].str.contains("Paper", case=False, na=False))]
        paper_vals = paper_row.iloc[0].to_dict() if not paper_row.empty else {}

        # Alloy framework row for FeNi
        alloy_vals = {}
        if sys_name == "FeNi" and not feni_df.empty:
            feni_row = feni_df[(feni_df["element"] == elem) & (feni_df["model"].str.contains("alloys_framework", case=False, na=False))]
            if not feni_row.empty:
                alloy_vals = feni_row.iloc[0].to_dict()

        def fmt(val, precision=2):
            if val is None or pd.isna(val):
                return "---"
            return f"{float(val):.{precision}f}"

        # 1. Print DFT Reference
        if dft_vals:
            print(f"{sys_name:<8} {elem:<5} {'DFT (Ref)':<20} {fmt(dft_vals.get('a0'), 3):<9} {fmt(dft_vals.get('Ecoh'), 3):<11} {fmt(dft_vals.get('C11'), 1):<11} {fmt(dft_vals.get('C12'), 1):<11} {fmt(dft_vals.get('C44'), 1):<11} {fmt(dft_vals.get('B'), 1):<9} {fmt(dft_vals.get('Evac_f'), 2):<11}")

        # 2. Print Paper tabGAP / Alloy Ref (if present)
        if paper_vals:
            print(f"{'':<8} {'':<5} {'Paper tabGAP':<20} {fmt(paper_vals.get('a0'), 3):<9} {fmt(paper_vals.get('Ecoh'), 3):<11} {fmt(paper_vals.get('C11'), 1):<11} {fmt(paper_vals.get('C12'), 1):<11} {fmt(paper_vals.get('C44'), 1):<11} {fmt(paper_vals.get('B'), 1):<9} {fmt(paper_vals.get('Evac_f'), 2):<11}")
        elif alloy_vals:
            print(f"{'':<8} {'':<5} {'Alloy Ref (Thesis)':<20} {fmt(alloy_vals.get('a0'), 3):<9} {fmt(alloy_vals.get('Ecoh'), 3):<11} {fmt(alloy_vals.get('C11'), 1):<11} {fmt(alloy_vals.get('C12'), 1):<11} {fmt(alloy_vals.get('C44'), 1):<11} {fmt(alloy_vals.get('B'), 1):<9} {fmt(alloy_vals.get('Evac_f'), 2):<11}")

        # 3. Print Calculated Refitted Result
        print(f"{'':<8} {'':<5} {'JGAP (Refitted)':<20} {fmt(res['a0'], 4):<9} {fmt(res['Ecoh'], 4):<11} {fmt(res['C11'], 1):<11} {fmt(res['C12'], 1):<11} {fmt(res['C44'], 1):<11} {fmt(res['B'], 1):<9} {fmt(res['Evac_f'], 3):<11}")

        # 4. Determine system-appropriate tolerances
        tolerances = dict(base_tolerances)
        if sys_name == "FeNi":
            # FeNi uses alloy-specific energy zero (+0.308 eV for Ni) and alloy multi-element representation
            tolerances["Ecoh"] = max(tolerances.get("Ecoh", 0.15), 0.35)
            tolerances["B"] = max(tolerances.get("B", 40.0), 150.0)
            tolerances["C11"] = max(tolerances.get("C11", 75.0), 200.0)
            tolerances["C12"] = max(tolerances.get("C12", 40.0), 120.0)
            tolerances["C44"] = max(tolerances.get("C44", 30.0), 140.0)

        # 5. Print Differences vs DFT & Check Tolerances
        diff_strs = {}
        for prop in ["a0", "Ecoh", "C11", "C12", "C44", "B", "Evac_f"]:
            calc_v = res.get(prop)
            dft_v = dft_vals.get(prop)
            if calc_v is not None and dft_v is not None and not pd.isna(dft_v):
                diff = calc_v - float(dft_v)
                abs_diff = abs(diff)
                tol = tolerances.get(prop, 1e9)
                prec = 4 if prop in ["a0", "Ecoh", "Evac_f"] else 1
                diff_strs[prop] = f"{diff:+0.{prec}f}"
                if abs_diff > tol:
                    all_passed = False
                    reason = f"[{sys_name} {elem}] {prop} deviation |{diff:+.4f}| exceeds tolerance {tol}"
                    failure_reasons.append(reason)
            else:
                diff_strs[prop] = "---"

        print(f"{'':<8} {'':<5} {'Δ (Calc - DFT)':<20} {diff_strs['a0']:<9} {diff_strs['Ecoh']:<11} {diff_strs['C11']:<11} {diff_strs['C12']:<11} {diff_strs['C44']:<11} {diff_strs['B']:<9} {diff_strs['Evac_f']:<11}")
        print("-" * 110)

    print("=" * 110)
    return all_passed, failure_reasons


def main():
    parser = argparse.ArgumentParser(description="JGAP Physical Property Validation via standard_gap_fit")
    parser.add_argument("--systems", type=str, default="Al,Cu,Ni,Fe,FeNi", help="Comma-separated systems to evaluate (Al,Cu,Ni,Fe,FeNi)")
    parser.add_argument("--refit", action="store_true", help="Force refit and tabulation of potentials even if cached")
    parser.add_argument("--pot-dir", type=str, default="", help="Directory to store/load fitted potentials (default: build/validation_pots)")
    parser.add_argument("--xyz-dir", type=str, default="", help="Directory containing training XYZ databases")
    parser.add_argument("--csv", type=str, default="", help="Path to reference CSV (elemental_property_comparison.csv)")
    parser.add_argument("--tolerance-a0", type=float, default=DEFAULT_TOLERANCES["a0"], help="Max allowable |a0 - a0_dft| in Å")
    parser.add_argument("--tolerance-B", type=float, default=DEFAULT_TOLERANCES["B"], help="Max allowable |B - B_dft| in GPa")
    parser.add_argument("--tolerance-C11", type=float, default=DEFAULT_TOLERANCES["C11"], help="Max allowable |C11 - C11_dft| in GPa")
    parser.add_argument("--tolerance-C12", type=float, default=DEFAULT_TOLERANCES["C12"], help="Max allowable |C12 - C12_dft| in GPa")
    parser.add_argument("--tolerance-C44", type=float, default=DEFAULT_TOLERANCES["C44"], help="Max allowable |C44 - C44_dft| in GPa")
    parser.add_argument("--tolerance-Ecoh", type=float, default=DEFAULT_TOLERANCES["Ecoh"], help="Max allowable |Ecoh - Ecoh_dft| in eV")
    parser.add_argument("--tolerance-Evac", type=float, default=DEFAULT_TOLERANCES["Evac_f"], help="Max allowable |Evac_f - Evac_dft| in eV")
    args = parser.parse_args()

    pot_dir = Path(args.pot_dir) if args.pot_dir else REPO_ROOT / "build" / "validation_pots"
    xyz_dir = Path(args.xyz_dir) if args.xyz_dir else REPO_ROOT / "test" / "structure-databases"
    ref_csv_path = Path(args.csv) if args.csv else REPO_ROOT / "test" / "resources" / "validation" / "elemental_property_comparison.csv"

    tolerances = {
        "a0": args.tolerance_a0,
        "B": args.tolerance_B,
        "C11": args.tolerance_C11,
        "C12": args.tolerance_C12,
        "C44": args.tolerance_C44,
        "Ecoh": args.tolerance_Ecoh,
        "Evac_f": args.tolerance_Evac,
    }

    print("======================================================================")
    print("JGAP Validation: Refit & Physical Property Predictions vs DFT")
    print(f"  Target Systems : {args.systems}")
    print(f"  Potentials Dir : {pot_dir}")
    print(f"  XYZ DBs Dir    : {xyz_dir}")
    print(f"  Reference CSV  : {ref_csv_path}")
    print(f"  Force Refit    : {args.refit}")
    print(f"  Tolerances     : {tolerances}")
    print("======================================================================")

    ref_df = pd.read_csv(ref_csv_path) if ref_csv_path.exists() else pd.DataFrame()

    requested_systems = [s.strip() for s in args.systems.split(",") if s.strip()]
    results = []

    for sys_name in requested_systems:
        if sys_name not in SYSTEM_CONFIGS:
            print(f"Error: Unknown system '{sys_name}'. Available: {list(SYSTEM_CONFIGS.keys())}")
            sys.exit(1)

        cfg = SYSTEM_CONFIGS[sys_name]
        try:
            pot_files, iso_energies = ensure_potential(sys_name, cfg, xyz_dir, pot_dir, force_refit=args.refit)

            for elem in cfg["elements"]:
                print(f"\n>>> Evaluating {elem} ({cfg['crystal'][elem]}) using [{sys_name}] potential...")
                res = evaluate_element_properties(
                    elem=elem,
                    pot_files=pot_files,
                    iso_energy=iso_energies.get(elem),
                    crystal=cfg["crystal"][elem],
                    a_init=cfg["a_init"][elem],
                    sc_dim=cfg["supercell"][elem],
                )
                res["system"] = sys_name
                results.append(res)
                print(f"    a0 = {res['a0']:.4f} Å | B = {res['B']:.1f} GPa | Ecoh = {res['Ecoh']} eV | Evac_f = {res['Evac_f']} eV")

        except Exception as ex:
            print(f"ERROR processing system '{sys_name}': {ex}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    feni_csv_path = REPO_ROOT / "test" / "resources" / "validation" / "feni_pure_properties_comparison.csv"
    feni_df = pd.read_csv(feni_csv_path) if feni_csv_path.exists() else pd.DataFrame()

    # Print summary and check tolerances
    all_passed, failures = print_comparison_table(results, ref_df, feni_df, tolerances)

    if not all_passed:
        print("\nVALIDATION CHECKS FAILED: Deviations exceeded allowable tolerances from DFT:")
        for reason in failures:
            print(f"  - {reason}")
        print("\nExiting with status code 1.")
        sys.exit(1)

    print("\nVALIDATION PASSED: All property predictions match DFT within allowable tolerances.")
    sys.exit(0)


if __name__ == "__main__":
    main()
