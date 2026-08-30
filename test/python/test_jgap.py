import os
import sys
import tempfile
import numpy as np
import unittest

# Ensure python path has the jgap module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../python")))

import jgap
from jgap.ase import JGAPCalculator


class TestJGAP(unittest.TestCase):

    def test_species(self):
        s_fe = jgap.Species("Fe")
        self.assertEqual(s_fe.symbol, "Fe")
        self.assertEqual(s_fe.atomic_number, 26)
        self.assertAlmostEqual(s_fe.mass, 55.845, places=2)
        self.assertEqual(repr(s_fe), "<Species 'Fe'>")
        self.assertEqual(str(s_fe), "Fe")

        s_fe_num = jgap.Species.from_atomic_number(26)
        self.assertEqual(s_fe_num, s_fe)

    def test_vector3(self):
        v1 = jgap.Vector3(1.0, 2.0, 3.0)
        v2 = jgap.Vector3(4.0, 5.0, 6.0)

        self.assertEqual(v1.x, 1.0)
        self.assertEqual(v1.y, 2.0)
        self.assertEqual(v1.z, 3.0)

        v_sum = v1 + v2
        self.assertEqual(v_sum.x, 5.0)
        self.assertEqual(v_sum.y, 7.0)
        self.assertEqual(v_sum.z, 9.0)

        v_scaled = v1 * 2.0
        self.assertEqual(v_scaled.x, 2.0)
        self.assertEqual(v_scaled.y, 4.0)
        self.assertEqual(v_scaled.z, 6.0)

        np_arr = v1.to_numpy()
        np.testing.assert_allclose(np_arr, [1.0, 2.0, 3.0])

        v_from_np = jgap.Vector3.from_numpy(np.array([7.0, 8.0, 9.0]))
        self.assertEqual(v_from_np.x, 7.0)
        self.assertEqual(v_from_np.y, 8.0)
        self.assertEqual(v_from_np.z, 9.0)

    def test_lattice(self):
        v1 = jgap.Vector3(5.0, 0.0, 0.0)
        v2 = jgap.Vector3(0.0, 5.0, 0.0)
        v3 = jgap.Vector3(0.0, 0.5, 5.0)

        lat = jgap.Lattice(v1, v2, v3)
        self.assertAlmostEqual(lat.volume(), 125.0, places=5)

        mat = lat.to_numpy()
        self.assertEqual(mat.shape, (3, 3))
        np.testing.assert_allclose(mat[0], [5.0, 0.0, 0.0])

        lat2 = jgap.Lattice(mat)
        self.assertAlmostEqual(lat2.volume(), 125.0, places=5)

    def test_virials(self):
        v = jgap.Virials(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
        self.assertEqual(v.xx, 1.0)
        self.assertEqual(v.xy, 2.0)
        self.assertEqual(v.xz, 3.0)
        self.assertEqual(v.yy, 4.0)
        self.assertEqual(v.yz, 5.0)
        self.assertEqual(v.zz, 6.0)

        voigt = v.to_voigt()
        np.testing.assert_allclose(voigt, [1.0, 4.0, 6.0, 5.0, 3.0, 2.0])

        mat = v.to_matrix()
        self.assertEqual(mat.shape, (3, 3))
        self.assertEqual(mat[0, 1], 2.0)
        self.assertEqual(mat[1, 0], 2.0)
        self.assertEqual(mat[0, 2], 3.0)
        self.assertEqual(mat[2, 0], 3.0)
        self.assertEqual(mat[1, 2], 5.0)
        self.assertEqual(mat[2, 1], 5.0)

        v2 = jgap.Virials(voigt)
        np.testing.assert_allclose(v2.to_voigt(), voigt)

    def test_atoms_properties_and_ase_interop(self):
        pos = np.array([[0.0, 0.0, 0.0], [1.4, 1.4, 1.4]])
        symbols = ["Fe", "Ni"]
        lat = jgap.Lattice(np.eye(3) * 10.0)
        pbc = [True, True, True]

        atoms = jgap.Atoms(positions=pos, symbols=symbols, lattice=lat, pbc=pbc)
        self.assertEqual(len(atoms), 2)
        self.assertEqual(atoms.n_atoms(), 2)
        self.assertEqual(atoms.symbols, ["Fe", "Ni"])
        np.testing.assert_allclose(atoms.positions, pos)
        self.assertEqual(list(atoms.pbc), [True, True, True])

        # Modification via properties
        new_pos = np.array([[0.1, 0.1, 0.1], [1.5, 1.5, 1.5]])
        atoms.positions = new_pos
        np.testing.assert_allclose(atoms.positions, new_pos)

        # Optional properties
        atoms.energy = -15.5
        self.assertEqual(atoms.energy, -15.5)
        forces = np.array([[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]])
        atoms.forces = forces
        np.testing.assert_allclose(atoms.forces, forces)

        # ASE round-trip
        ase_atoms = atoms.to_ase()
        self.assertEqual(len(ase_atoms), 2)
        self.assertEqual(ase_atoms.get_chemical_symbols(), ["Fe", "Ni"])
        np.testing.assert_allclose(ase_atoms.get_positions(), new_pos)
        np.testing.assert_allclose(ase_atoms.cell.array, np.eye(3) * 10.0)
        self.assertEqual(list(ase_atoms.pbc), [True, True, True])

        # Convert back from ASE
        atoms_from_ase = jgap.Atoms.from_ase(ase_atoms)
        self.assertEqual(len(atoms_from_ase), 2)
        self.assertEqual(atoms_from_ase.symbols, ["Fe", "Ni"])
        np.testing.assert_allclose(atoms_from_ase.positions, new_pos)

    def test_regularization_rules(self):
        defaults = jgap.PerConfigTypeSigmas(0.001, 0.05, 0.1, 0.02)
        rules = jgap.PerConfigTypeRegularizationRules(
            default_sigmas=defaults,
            exact_config_type_sigmas={
                "liquid": jgap.PerConfigTypeSigmas(0.01, 0.2, 0.5, 0.1)
            },
            config_type_contains_sigmas={
                "melt": jgap.PerConfigTypeSigmas(0.02, 0.3, 0.6, 0.15),
                "melt_high_temp": jgap.PerConfigTypeSigmas(0.05, 0.4, 0.8, 0.2),
            }
        )

        pos = np.array([[0.0, 0.0, 0.0], [1.4, 1.4, 1.4]])
        atoms1 = jgap.Atoms(positions=pos, symbols=["Fe", "Fe"])
        atoms1.config_type = "bulk" # fallback to default

        atoms2 = jgap.Atoms(positions=pos, symbols=["Fe", "Fe"])
        atoms2.config_type = "liquid" # exact match

        atoms3 = jgap.Atoms(positions=pos, symbols=["Fe", "Fe"])
        atoms3.config_type = "fast_melt_traj" # contains "melt"

        atoms4 = jgap.Atoms(positions=pos, symbols=["Fe", "Fe"])
        atoms4.config_type = "fast_melt_high_temp_traj" # contains "melt" and "melt_high_temp" -> longest match

        sigmas1 = rules.determine(atoms1)
        sigmas2 = rules.determine(atoms2)
        sigmas3 = rules.determine(atoms3)
        sigmas4 = rules.determine(atoms4)

        self.assertAlmostEqual(sigmas1.energy, 0.001)
        self.assertAlmostEqual(sigmas2.energy, 0.01)
        self.assertAlmostEqual(sigmas3.energy, 0.02)
        self.assertAlmostEqual(sigmas4.energy, 0.05)

        all_sigmas = rules.determine_for_all([atoms1, atoms2, atoms3, atoms4])
        self.assertEqual(len(all_sigmas), 4)
        self.assertAlmostEqual(all_sigmas[0].energy, sigmas1.energy)
        self.assertAlmostEqual(all_sigmas[1].energy, sigmas2.energy)
        self.assertAlmostEqual(all_sigmas[2].energy, sigmas3.energy)
        self.assertAlmostEqual(all_sigmas[3].energy, sigmas4.energy)

    def test_read_xyz(self):
        xyz_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../resources/xyz-samples/feni-test.xyz")
        )
        if os.path.exists(xyz_path):
            frames = jgap.read_atoms(xyz_path)
            self.assertGreater(len(frames), 0)
            first_frame = frames[0]
            self.assertGreater(len(first_frame), 0)
            self.assertEqual(first_frame.positions.ndim, 2)
            self.assertEqual(first_frame.positions.shape[1], 3)

    def test_standard_gap_fit_and_tabulate(self):
        xyz_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../resources/xyz-samples/feni-test.xyz")
        )
        if not os.path.exists(xyz_path):
            return

        frames = jgap.read_atoms(xyz_path)
        training_frames = frames[:3]

        params = jgap.StandardGapParams(
            seed=42,
            cutoff2=3.5,
            n_sparse2=5,
            eam_pair_function=jgap.EamPairFunctionType.FSGen3,
            eam_n_sparse=5,
            eam_min_density=0.05,
            cutoff3=3.0,
            n_sparse3=20,
        )

        rules = jgap.SimpleRegularizationRules()
        sigmas = rules.determine_for_all(training_frames)

        with tempfile.TemporaryDirectory() as tmpdir:
            pot_file = os.path.join(tmpdir, "model.jgap.h5")

            # 1. Fit potential and save directly to file
            jgap.standard_gap_fit(pot_file, training_frames, sigmas, params)
            self.assertTrue(os.path.exists(pot_file))

            # 2. Load potential & evaluate
            loaded_pot = jgap.load_potential(pot_file)
            test_frame = training_frames[0]
            res = loaded_pot.calculate_energy(test_frame)
            self.assertIsInstance(res.energy, float)
            self.assertEqual(res.forces.shape, (len(test_frame), 3))

            # 3. Tabulate
            tab_prefix = os.path.join(tmpdir, "model_tab")
            tab_params = jgap.StandardTabulationParams(
                r_min_3b=1.0,
                max_eam_density=8.0,
                n_grid_2b=500,
                n_grid_3b=[20, 20, 20],
            )
            jgap.standard_tabulation(pot_file, tab_prefix, tab_params)
            self.assertTrue(os.path.exists(f"{tab_prefix}.tabgap.h5"))

            # 4. Load TabGap & Use with ASE Calculator
            loaded_tabgap = jgap.load_potential(f"{tab_prefix}.tabgap.h5")
            ase_atoms = test_frame.to_ase()
            ase_atoms.calc = JGAPCalculator(loaded_tabgap)
            e = ase_atoms.get_potential_energy()
            f = ase_atoms.get_forces()
            s = ase_atoms.get_stress()
            self.assertIsInstance(e, float)
            self.assertEqual(f.shape, (len(test_frame), 3))
            self.assertEqual(len(s), 6)

    def test_approx_ram_limit_gb(self):
        xyz_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../resources/xyz-samples/feni-test.xyz")
        )
        if not os.path.exists(xyz_path):
            return

        frames = jgap.read_atoms(xyz_path)[:3]
        rules = jgap.SimpleRegularizationRules()
        sigmas = rules.determine_for_all(frames)

        with tempfile.TemporaryDirectory() as tmpdir:
            # 1. Tiny RAM limit (< M*M) -> throws error
            params_tiny = jgap.StandardGapParams(
                seed=42,
                n_sparse2=10,
                n_sparse3=50,
                approx_ram_limit_gb=1e-12,
            )
            pot_file_tiny = os.path.join(tmpdir, "tiny_ram.jgap.h5")
            with self.assertRaises(RuntimeError):
                jgap.standard_gap_fit(pot_file_tiny, frames, sigmas, params_tiny)

            # 2. Tight RAM limit between M*M and 2*M*M -> triggers StreamingQrGapFit with small chunk rows
            # M is ~370 sparse points -> M*M*8 bytes ~ 1.1 MB (0.00102 GB). 2*M*M*8 ~ 2.19 MB (0.00204 GB).
            params_inplace = jgap.StandardGapParams(
                seed=42,
                n_sparse2=10,
                n_sparse3=50,
                approx_ram_limit_gb=0.0015, # ~1.6 MB (between M*M and 2*M*M)
            )
            pot_file_inplace = os.path.join(tmpdir, "inplace.jgap.h5")
            jgap.standard_gap_fit(pot_file_inplace, frames, sigmas, params_inplace)
            self.assertTrue(os.path.exists(pot_file_inplace))

            # 3. Medium RAM limit between 2*M*M and (N+M)*M -> triggers StreamingQrGapFit with larger chunk rows
            # (N+M)*M*8 ~ 4.46 MB (0.00415 GB).
            params_streaming = jgap.StandardGapParams(
                seed=42,
                n_sparse2=10,
                n_sparse3=50,
                approx_ram_limit_gb=0.0030, # ~3.2 MB (between 2*M*M and (N+M)*M)
            )
            pot_file_streaming = os.path.join(tmpdir, "streaming.jgap.h5")
            jgap.standard_gap_fit(pot_file_streaming, frames, sigmas, params_streaming)
            self.assertTrue(os.path.exists(pot_file_streaming))

            # 4. StandardGapParams with split_sets -> triggers SplitQRGapFit
            params_split = jgap.StandardGapParams(
                seed=42,
                n_sparse2=10,
                n_sparse3=50,
                split_sets=[["Fe"], ["Ni"]],
                approx_ram_limit_gb=0.0030,
            )
            pot_file_split = os.path.join(tmpdir, "split.jgap.h5")
            jgap.standard_gap_fit(pot_file_split, frames, sigmas, params_split)
            self.assertTrue(os.path.exists(pot_file_split))


if __name__ == "__main__":
    unittest.main()
