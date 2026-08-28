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

    def test_standard_gap_params(self):
        rules = jgap.SimpleRegularizationRules(0.002, 0.04, 0.08, 0.01)
        pf = jgap.FSGenPairFunction(4.5, 3.0)

        params = jgap.StandardGapParams(
            seed=123,
            cutoff2=4.0,
            n_sparse2=15,
            eam_mode=jgap.EamMode.Blind,
            eam_pf=pf,
            eam_n_sparse=10,
            eam_min_density=0.05,
            cutoff3=3.5,
            n_sparse3=100,
            regularization_rules=rules,
        )

        self.assertEqual(params.seed, 123)
        self.assertEqual(params.cutoff2, 4.0)
        self.assertEqual(params.n_sparse2, 15)
        self.assertEqual(params.eam_n_sparse, 10)
        self.assertEqual(params.cutoff3, 3.5)
        self.assertEqual(params.n_sparse3, 100)
        self.assertIn("StandardGapParams", repr(params))

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
            eam_n_sparse=5,
            eam_min_density=0.05,
            cutoff3=3.0,
            n_sparse3=20,
        )

        # 1. Fit potential
        gap_pot = jgap.standard_gap_fit(training_frames, params)
        self.assertIsInstance(gap_pot, jgap.GapPotential)
        self.assertGreater(gap_pot.num_components(), 0)

        # 2. Evaluate energy & forces
        test_frame = training_frames[0]
        res = gap_pot.calculate_energy(test_frame)
        self.assertIsInstance(res.energy, float)
        self.assertEqual(res.forces.shape, (len(test_frame), 3))

        with tempfile.TemporaryDirectory() as tmpdir:
            # 3. Save & Load GapPotential
            pot_file = os.path.join(tmpdir, "model.jgap.h5")
            gap_pot.save(pot_file)

            loaded_pot = jgap.load_potential(pot_file)
            res_loaded = loaded_pot.calculate_energy(test_frame)
            self.assertAlmostEqual(res.energy, res_loaded.energy, places=6)
            np.testing.assert_allclose(res.forces, res_loaded.forces, rtol=1e-5, atol=1e-5)

            # 4. Tabulate
            tab_params = jgap.TabulationParams(
                max_cutoffs=gap_pot.get_cutoffs(),
                r_min_3b=1.0,
                max_eam_density=8.0,
                n_grid_2b=500,
                n_grid_3b=[20, 20, 20],
            )
            tabgap_pot = gap_pot.tabulate(tab_params)
            self.assertIsInstance(tabgap_pot, jgap.TabGapPotential)

            res_tab = tabgap_pot.calculate_energy(test_frame)
            self.assertIsInstance(res_tab.energy, float)

            # 5. Save & Load TabGapPotential
            tab_prefix = os.path.join(tmpdir, "model_tab")
            written_files = jgap.save_tabgap(tabgap_pot, tab_prefix)
            self.assertGreater(len(written_files), 0)

            loaded_tabgap = jgap.load_potential(written_files)
            res_tab_loaded = loaded_tabgap.calculate_energy(test_frame)
            self.assertAlmostEqual(res_tab.energy, res_tab_loaded.energy, places=6)

            # 6. Use with ASE Calculator
            ase_atoms = test_frame.to_ase()
            ase_atoms.calc = JGAPCalculator(loaded_tabgap)
            e = ase_atoms.get_potential_energy()
            f = ase_atoms.get_forces()
            s = ase_atoms.get_stress()
            self.assertAlmostEqual(e, res_tab_loaded.energy, places=6)
            np.testing.assert_allclose(f, res_tab_loaded.forces, rtol=1e-5, atol=1e-5)
            self.assertEqual(len(s), 6)


if __name__ == "__main__":
    unittest.main()
