#
# Copyright (C) 2026 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import unittest as ut
import tempfile

import espressomd
import pyMBE
import pyMBE.lib.nanoparticle_tools as nanoparticle_tools
from pyMBE.storage.instances.nanoparticle import NanoparticleInstance
from pyMBE.storage.templates.particle import ParticleTemplate
from pyMBE.storage.pint_quantity import PintQuantity


# ESPResSo allows only one System instance per process.
espresso_system = espressomd.System(box_l=[40, 40, 40])


class TestNanoparticleCreation(ut.TestCase):

    def setUp(self):
        """
        Reset ESPResSo state before each unit test.

        Notes:
            - ESPResSo allows one ``System`` instance per process in this test
            file, so particles are cleared between tests.
        """
        espresso_system.part.clear()

    def _build_pmb_with_particles(self):
        """
        Build a pyMBE object with particle templates used by nanoparticle tests.

        Returns:
            ('pyMBE.pymbe_library'):
                Configured pyMBE object containing templates for ``core``,
                ``A``, and ``B`` particles.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.set_reduced_units(unit_length=0.4 * pmb.units.nm,
                              Kw=1e-14)

        pmb.define_particle(name="core",
                            z=0,
                            sigma=4.0 * pmb.units("reduced_length"),
                            epsilon=1.0 * pmb.units("reduced_energy"))
        pmb.define_particle(name="A",
                            z=-1,
                            sigma=1.0 * pmb.units("reduced_length"),
                            epsilon=1.0 * pmb.units("reduced_energy"))
        pmb.define_particle(name="B",
                            z=1,
                            sigma=1.0 * pmb.units("reduced_length"),
                            epsilon=1.0 * pmb.units("reduced_energy"))
        return pmb

    def test_create_nanoparticle_registers_instances(self):
        """
        Unit test: verify nanoparticle creation and database bookkeeping.

        Notes:
            - Checks nanoparticle instance IDs, particle-to-molecule linkage, and
            particle count consistency with template-derived properties.
        """
        pmb = self._build_pmb_with_particles()
        nanoparticle_name = "np"
        pmb.define_nanoparticle(name=nanoparticle_name,
                                core_particle_name="core",
                                surface_density_of_sites=0.2 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=0.5,
                                number_of_patches_of_primary_sites=2,
                                secondary_site_particle_name="B")

        core_positions = [[5.0, 5.0, 5.0],
                          [25.0, 25.0, 25.0]]
        created_np_ids = pmb.create_nanoparticle(name=nanoparticle_name,
                                                 number_of_nanoparticles=2,
                                                 espresso_system=espresso_system,
                                                 list_core_particle_positions=core_positions,
                                                 fix=True)

        self.assertEqual(created_np_ids, [0, 1])

        np_instances = pmb.db.get_instances(pmb_type="nanoparticle")
        self.assertEqual(sorted(np_instances.keys()), [0, 1])
        self.assertEqual(np_instances[0].molecule_id, 0)
        self.assertEqual(np_instances[1].molecule_id, 1)

        tpl = pmb.db.get_template(pmb_type="nanoparticle",
                                  name=nanoparticle_name)
        properties = tpl.calculate_nanoparticle_properties(pmb)
        expected_particles_per_np = 1 + properties["number_of_primary_sites"] + properties["number_of_secondary_sites"]

        for nanoparticle_id, expected_core_pos in enumerate(core_positions):
            particles_in_np = pmb.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                                     attribute="molecule_id",
                                                                     value=nanoparticle_id)
            self.assertEqual(len(particles_in_np), expected_particles_per_np)

            core_candidates = [pid for pid in particles_in_np
                               if pmb.db.get_instance(pmb_type="particle", instance_id=pid).name == "core"]
            self.assertEqual(len(core_candidates), 1)
            core_pid = core_candidates[0]

            core_pos = list(espresso_system.part.by_id(core_pid).pos)
            self.assertListEqual(core_pos, expected_core_pos)
            self.assertListEqual(list(espresso_system.part.by_id(core_pid).fix), [True, True, True])

            for pid in particles_in_np:
                self.assertEqual(pmb.db.get_instance(pmb_type="particle", instance_id=pid).molecule_id,
                                 nanoparticle_id)
                self.assertListEqual(list(espresso_system.part.by_id(pid).fix), [True, True, True])

        particle_id_map = pmb.get_particle_id_map(object_name=nanoparticle_name)
        self.assertEqual(len(particle_id_map["all"]), expected_particles_per_np * 2)

    def test_create_nanoparticle_input_validation_and_empty(self):
        """
        Unit test: verify input validation and empty-creation behavior.

        Notes:
            - Covers zero requested nanoparticles and invalid core-position inputs.
        """
        pmb = self._build_pmb_with_particles()
        pmb.define_nanoparticle(name="np",
                                core_particle_name="core",
                                surface_density_of_sites=0.2 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=0.5,
                                number_of_patches_of_primary_sites=2,
                                secondary_site_particle_name="B")

        self.assertEqual(pmb.create_nanoparticle(name="np",
                                                 number_of_nanoparticles=0,
                                                 espresso_system=espresso_system),
                         [])

        with self.assertRaises(ValueError):
            pmb.create_nanoparticle(name="np",
                                    number_of_nanoparticles=2,
                                    espresso_system=espresso_system,
                                    list_core_particle_positions=[[1.0, 2.0, 3.0]])

        with self.assertRaises(ValueError):
            pmb.create_nanoparticle(name="np",
                                    number_of_nanoparticles=1,
                                    espresso_system=espresso_system,
                                    list_core_particle_positions=[[1.0, 2.0]])

    def test_create_nanoparticle_calls_enable_motion_when_not_fixed(self):
        """
        Unit test: verify rigid-body motion hook is called when ``fix=False``.
        """
        pmb = self._build_pmb_with_particles()
        pmb.define_nanoparticle(name="np",
                                core_particle_name="core",
                                surface_density_of_sites=0.2 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=0.5,
                                number_of_patches_of_primary_sites=2,
                                secondary_site_particle_name="B")

        calls = []

        def fake_enable_motion_of_rigid_object(instance_id, pmb_type, espresso_system):
            _ = espresso_system
            calls.append((instance_id, pmb_type))

        pmb.enable_motion_of_rigid_object = fake_enable_motion_of_rigid_object

        created_np_ids = pmb.create_nanoparticle(name="np",
                                                 number_of_nanoparticles=1,
                                                 espresso_system=espresso_system,
                                                 fix=False)

        self.assertEqual(created_np_ids, [0])
        self.assertEqual(calls, [(0, "nanoparticle")])

    def test_create_nanoparticle_sites_positions_variants(self):
        """
        Unit test: verify site-position generation across template variants.

        Notes:
            - Covers two-patch, multi-patch, and zero-site configurations.
        """
        pmb = self._build_pmb_with_particles()

        # Two primary patches + secondary sites
        pmb.define_nanoparticle(name="np_two",
                                core_particle_name="core",
                                surface_density_of_sites=0.2 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=0.5,
                                number_of_patches_of_primary_sites=2,
                                secondary_site_particle_name="B")
        tpl_two = pmb.db.get_template(pmb_type="nanoparticle", name="np_two")
        properties_two = tpl_two.calculate_nanoparticle_properties(pmb)
        sites_two = pmb._create_nanoparticle_sites_positions(nanoparticle_tpl=tpl_two)
        self.assertEqual(len(sites_two), 3)
        self.assertEqual(sum(item["number_of_sites"] for item in sites_two), properties_two["total_number_of_sites"])

        # More than two primary patches, no secondary type.
        # Overlap detection can be sensitive to geometric degeneracy for small systems;
        # here we stub overlap validation to exercise this branch deterministically.
        pmb.define_nanoparticle(name="np_three",
                                core_particle_name="core",
                                surface_density_of_sites=0.2 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=1.0,
                                number_of_patches_of_primary_sites=3,
                                secondary_site_particle_name=None)
        tpl_three = pmb.db.get_template(pmb_type="nanoparticle", name="np_three")
        original_check_patch_overlaps = nanoparticle_tools.check_patch_overlaps
        nanoparticle_tools.check_patch_overlaps = lambda sites_positions, number_patches: 0
        try:
            sites_three = pmb._create_nanoparticle_sites_positions(nanoparticle_tpl=tpl_three)
        finally:
            nanoparticle_tools.check_patch_overlaps = original_check_patch_overlaps
        self.assertEqual(len(sites_three), 3)
        self.assertTrue(all(item["particle_name"] == "A" for item in sites_three))

        # Zero sites (density = 0) should return an empty specification
        pmb.define_nanoparticle(name="np_zero",
                                core_particle_name="core",
                                surface_density_of_sites=0.0 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=1.0,
                                number_of_patches_of_primary_sites=1,
                                secondary_site_particle_name=None)
        tpl_zero = pmb.db.get_template(pmb_type="nanoparticle", name="np_zero")
        self.assertEqual(pmb._create_nanoparticle_sites_positions(nanoparticle_tpl=tpl_zero), [])

    def test_create_nanoparticle_zero_site_patch_branch(self):
        """
        Unit test: cover branch where a generated patch has zero sites.
        """
        pmb = self._build_pmb_with_particles()
        pmb.define_nanoparticle(name="np",
                                core_particle_name="core",
                                surface_density_of_sites=0.2 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=0.5,
                                number_of_patches_of_primary_sites=2,
                                secondary_site_particle_name="B")

        original_helper = pmb._create_nanoparticle_sites_positions
        pmb._create_nanoparticle_sites_positions = lambda nanoparticle_tpl: [
            {"particle_name": "A", "positions": [], "number_of_sites": 0},
            {"particle_name": "B", "positions": [[0.0, 0.0, 0.0]], "number_of_sites": 1},
        ]
        try:
            created_ids = pmb.create_nanoparticle(name="np",
                                                  number_of_nanoparticles=1,
                                                  espresso_system=espresso_system,
                                                  list_core_particle_positions=[[1.0, 1.0, 1.0]],
                                                  fix=True)
        finally:
            pmb._create_nanoparticle_sites_positions = original_helper
        self.assertEqual(created_ids, [0])

    def test_nanoparticle_tools_standalone_functions(self):
        """
        Unit test: verify helper functions in ``nanoparticle_tools``.
        """
        points = nanoparticle_tools.uniform_distribution_sites_on_sphere(number_of_edges=1, tolerance=1e-6)
        self.assertEqual(len(points), 1)
        self.assertEqual(len(points[0]), 3)

        distances = nanoparticle_tools.calculate_distance_vector_point([[0, 0, 0], [1, 0, 0]], [0, 0, 0])
        self.assertEqual(distances[0], 0.0)
        self.assertEqual(distances[1], 1.0)

        _, patch = nanoparticle_tools.define_patch(points=[[0, 0, 0], [1, 0, 0], [0, 1, 0]],
                                                   central_point=[0, 0, 0],
                                                   patch_size=2)
        self.assertEqual(len(patch), 2)

        self.assertEqual(nanoparticle_tools.check_patch_overlaps(sites_positions=[[(0, 0, 0)], [(1, 0, 0)]],
                                                                 number_patches=2), 0)
        with self.assertRaises(ValueError):
            nanoparticle_tools.check_patch_overlaps(sites_positions=[[(0, 0, 0)], [(0, 0, 0)]],
                                                    number_patches=2)

        mean_d, std_d, err_d = nanoparticle_tools.calculate_distance_between_points_on_sphere(
            points=[[(1.0, 0.0, 0.0),
                     (0.0, 1.0, 0.0),
                     (0.0, 0.0, 1.0),
                     (-1.0, 0.0, 0.0),
                     (0.0, -1.0, 0.0),
                     (0.0, 0.0, -1.0),
                     (0.707, 0.707, 0.0)]]
        )
        self.assertGreaterEqual(mean_d, 0.0)
        self.assertGreaterEqual(std_d, 0.0)
        self.assertGreaterEqual(err_d, 0.0)

        dipole_vec, dipole_mag = nanoparticle_tools.calculate_dipole_moment(charges=[1.0, -1.0],
                                                                             positions=[[1.0, 0.0, 0.0],
                                                                                        [0.0, 1.0, 0.0]])
        self.assertEqual(len(dipole_vec), 3)
        self.assertGreaterEqual(dipole_mag, 0.0)

        q_tensor, q_mag, eigenvalues = nanoparticle_tools.calculate_quadrupole_moment(
            charges=[1.0, -1.0],
            positions=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        )
        self.assertEqual(q_tensor.shape, (3, 3))
        self.assertEqual(len(eigenvalues), 3)
        self.assertGreaterEqual(q_mag, 0.0)

    def test_nanoparticle_template_edge_cases_and_manager_paths(self):
        """
        Unit test: verify nanoparticle-template edge cases and manager paths.

        Notes:
            - Covers template validation errors and nanoparticle-specific manager
            collection logic.
        """
        pmb = self._build_pmb_with_particles()
        pmb.define_nanoparticle(name="np",
                                core_particle_name="core",
                                surface_density_of_sites=0.2 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=0.5,
                                number_of_patches_of_primary_sites=2,
                                secondary_site_particle_name="B")

        # Manager nanoparticle branch in _collect_particle_templates and _get_templates_df
        particle_counts = pmb.db.get_particle_templates_under(template_name="np",
                                                              pmb_type="nanoparticle",
                                                              return_counts=True)
        self.assertEqual(set(particle_counts.keys()), {"core", "A", "B"})
        self.assertFalse(pmb.get_templates_df("nanoparticle").empty)

        # NanoparticleTemplate validations (fraction and number of patches)
        tpl = pmb.db.get_template(pmb_type="nanoparticle", name="np")
        with self.assertRaises(ValueError):
            tpl.copy(update={"fraction_primary_sites": -0.1}).calculate_nanoparticle_properties(pmb)
        with self.assertRaises(ValueError):
            tpl.copy(update={"number_of_patches_of_primary_sites": 0}).calculate_nanoparticle_properties(pmb)

        # Cover state_name is None path for _get_initial_state_charge_number
        pmb.db._register_template(ParticleTemplate(name="A_no_init",
                                                   sigma=PintQuantity.from_quantity(1.0 * pmb.units.reduced_length, "length", pmb.units),
                                                   epsilon=PintQuantity.from_quantity(1.0 * pmb.units.reduced_energy, "energy", pmb.units),
                                                   cutoff=PintQuantity.from_quantity(1.2 * pmb.units.reduced_length, "length", pmb.units),
                                                   offset=PintQuantity.from_quantity(0.0 * pmb.units.reduced_length, "length", pmb.units),
                                                   initial_state=None))
        pmb.define_particle_states(particle_name="A_no_init",
                                   states=[{"name": "A_no_init_state", "z": 0}])
        pmb.define_nanoparticle(name="np_no_init_site",
                                core_particle_name="core",
                                surface_density_of_sites=0.1 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A_no_init",
                                fraction_primary_sites=1.0,
                                number_of_patches_of_primary_sites=1,
                                secondary_site_particle_name=None)
        tpl_no_init = pmb.db.get_template(pmb_type="nanoparticle", name="np_no_init_site")
        self.assertIn("total_number_of_sites", tpl_no_init.calculate_nanoparticle_properties(pmb))

        # Cover core initial_state is None path
        pmb.db._register_template(ParticleTemplate(name="core_no_init",
                                                   sigma=PintQuantity.from_quantity(4.0 * pmb.units.reduced_length, "length", pmb.units),
                                                   epsilon=PintQuantity.from_quantity(1.0 * pmb.units.reduced_energy, "energy", pmb.units),
                                                   cutoff=PintQuantity.from_quantity(4.2 * pmb.units.reduced_length, "length", pmb.units),
                                                   offset=PintQuantity.from_quantity(0.0 * pmb.units.reduced_length, "length", pmb.units),
                                                   initial_state=None))
        pmb.define_nanoparticle(name="np_bad_core",
                                core_particle_name="core_no_init",
                                surface_density_of_sites=0.1 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=1.0,
                                number_of_patches_of_primary_sites=1,
                                secondary_site_particle_name=None)
        tpl_bad_core = pmb.db.get_template(pmb_type="nanoparticle", name="np_bad_core")
        with self.assertRaises(ValueError):
            tpl_bad_core.calculate_nanoparticle_properties(pmb)

        # Cover no-particle-state branch for _get_initial_state_charge_number
        pmb.db._register_template(ParticleTemplate(name="A_no_states",
                                                   sigma=PintQuantity.from_quantity(1.0 * pmb.units.reduced_length, "length", pmb.units),
                                                   epsilon=PintQuantity.from_quantity(1.0 * pmb.units.reduced_energy, "energy", pmb.units),
                                                   cutoff=PintQuantity.from_quantity(1.2 * pmb.units.reduced_length, "length", pmb.units),
                                                   offset=PintQuantity.from_quantity(0.0 * pmb.units.reduced_length, "length", pmb.units),
                                                   initial_state=None))
        pmb.define_nanoparticle(name="np_no_states_site",
                                core_particle_name="core",
                                surface_density_of_sites=0.1 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A_no_states",
                                fraction_primary_sites=1.0,
                                number_of_patches_of_primary_sites=1,
                                secondary_site_particle_name=None)
        tpl_no_states = pmb.db.get_template(pmb_type="nanoparticle", name="np_no_states_site")
        with self.assertRaises(ValueError):
            tpl_no_states.calculate_nanoparticle_properties(pmb)

    def test_nanoparticle_instance_and_io_roundtrip(self):
        """
        Unit test: verify nanoparticle instance validation and I/O roundtrip.

        Notes:
            - Confirms serialized nanoparticle templates and instances load back
            correctly from disk.
        """
        pmb = self._build_pmb_with_particles()
        pmb.define_nanoparticle(name="np",
                                core_particle_name="core",
                                surface_density_of_sites=0.2 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=0.5,
                                number_of_patches_of_primary_sites=2,
                                secondary_site_particle_name="B")

        with self.assertRaises(ValueError):
            NanoparticleInstance(name="np", molecule_id=-1)

        pmb.create_nanoparticle(name="np",
                                number_of_nanoparticles=1,
                                espresso_system=espresso_system,
                                fix=True)

        new_pmb = pyMBE.pymbe_library(seed=43)
        with tempfile.TemporaryDirectory() as tmp_directory:
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)

        self.assertFalse(new_pmb.get_templates_df(pmb_type="nanoparticle").empty)
        self.assertFalse(new_pmb.get_instances_df(pmb_type="nanoparticle").empty)


if __name__ == "__main__":
    ut.main()
