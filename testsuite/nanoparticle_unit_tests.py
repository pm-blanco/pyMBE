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
import math
from unittest.mock import patch

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
                                total_number_of_sites=10,
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
                                total_number_of_sites=10,
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
                                total_number_of_sites=10,
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
                                total_number_of_sites=10,
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
                                total_number_of_sites=10,
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
                                total_number_of_sites=0,
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
                                total_number_of_sites=10,
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
                                total_number_of_sites=10,
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
                                total_number_of_sites=5,
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
                                total_number_of_sites=5,
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
                                total_number_of_sites=5,
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
                                total_number_of_sites=10,
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

    def test_get_nanoparticle_properties(self):
        """
        Unit test: verify ``get_nanoparticle_properties`` returns the correct properties.

        Notes:
            - Checks that the function returns the same result as calling
              ``calculate_nanoparticle_properties`` directly on the template.
            - Verifies that all expected keys are present in the returned dict.
            - Verifies that a non-existent name raises an error.
        """
        pmb = self._build_pmb_with_particles()
        pmb.define_nanoparticle(name="np",
                                core_particle_name="core",
                                total_number_of_sites=10,
                                primary_site_particle_name="A",
                                fraction_primary_sites=0.5,
                                number_of_patches_of_primary_sites=2,
                                secondary_site_particle_name="B")

        properties = nanoparticle_tools.get_nanoparticle_properties(pmb, "np")

        expected_keys = ["nanoparticle_surface_area",
                         "nanoparticle_volume",
                         "total_number_of_sites",
                         "real_surface_density_of_sites",
                         "number_of_primary_sites",
                         "number_of_primary_sites_per_patch",
                         "number_of_secondary_sites",
                         "real_fraction_primary_sites",
                         "primary_site_charge_number",
                         "secondary_site_charge_number",
                         "total_charge",
                         "surface_charge_density",
                         "volume_charge_density"]
        for key in expected_keys:
            self.assertIn(key, properties)

        self.assertEqual(properties["total_number_of_sites"], 10)

        tpl = pmb.db.get_template(pmb_type="nanoparticle", name="np")
        direct = tpl.calculate_nanoparticle_properties(pmb)
        self.assertEqual(properties, direct)

        with self.assertRaises(Exception):
            nanoparticle_tools.get_nanoparticle_properties(pmb, "nonexistent_np")

    def test_print_nanoparticle_properties(self):
        """
        Unit test: verify ``print_nanoparticle_properties`` prints all keys.

        Notes:
            - Checks that every key in the properties dict appears in stdout.
            - Verifies that the nanoparticle name appears in the header.
            - Verifies that the function works without a name argument.
        """
        pmb = self._build_pmb_with_particles()
        pmb.define_nanoparticle(name="np",
                                core_particle_name="core",
                                total_number_of_sites=10,
                                primary_site_particle_name="A",
                                fraction_primary_sites=0.5,
                                number_of_patches_of_primary_sites=2,
                                secondary_site_particle_name="B")

        properties = nanoparticle_tools.get_nanoparticle_properties(pmb, "np")

        import io
        import sys
        captured = io.StringIO()
        sys.stdout = captured
        try:
            nanoparticle_tools.print_nanoparticle_properties(properties, name="np")
        finally:
            sys.stdout = sys.__stdout__
        output = captured.getvalue()

        self.assertIn("np", output)
        for key in properties:
            self.assertIn(key, output)

        # Also verify it runs without a name (no error, generic header)
        captured2 = io.StringIO()
        sys.stdout = captured2
        try:
            nanoparticle_tools.print_nanoparticle_properties(properties)
        finally:
            sys.stdout = sys.__stdout__
        self.assertGreater(len(captured2.getvalue()), 0)


class TestRelaxNanoparticleOverlaps(ut.TestCase):
    """Tests for relax_nanoparticle_overlaps in nanoparticle_tools."""

    def setUp(self):
        espresso_system.part.clear()

    def _build_pmb(self):
        """
        Build a pyMBE object for overlap-relaxation tests.

        Particle setup:
            - ``np_core``: sigma=1, offset=3 (mimics nanoparticle core with
              diameter=4, i.e. offset = diameter - sigma = 4 - 1).
            - ``np_ion``:  sigma=1, offset=0 (small salt ion).
            - ``np_site``: sigma=0, acidic (excluded from LJ automatically).

        Lorentz-Berthelot target offsets derived at test time from
        ``pmb.get_lj_parameters``; no expected value is hardcoded in the tests.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_particle(
            name="np_core",
            z=0,
            sigma=1.0 * pmb.units("reduced_length"),
            epsilon=1.0 * pmb.units("reduced_energy"),
            offset=3.0 * pmb.units("reduced_length"),
            cutoff=2 ** (1 / 6) * pmb.units("reduced_length"),
        )
        pmb.define_particle(
            name="np_ion",
            z=1,
            sigma=1.0 * pmb.units("reduced_length"),
            epsilon=1.0 * pmb.units("reduced_energy"),
        )
        pmb.define_particle(
            name="np_site",
            acidity="acidic",
            pka=4.0,
            sigma=0.0 * pmb.units("reduced_length"),
            epsilon=1.0 * pmb.units("reduced_energy"),
        )
        pmb.define_nanoparticle(
            name="np",
            core_particle_name="np_core",
            total_number_of_sites=4,
            primary_site_particle_name="np_site",
            fraction_primary_sites=1.0,
            number_of_patches_of_primary_sites=1,
            secondary_site_particle_name=None,
        )
        return pmb

    def _get_max_target_offset(self, pmb, core_name):
        """Return the largest Lorentz-Berthelot offset involving the core particle."""
        max_target = 0.0
        for name in pmb.db.get_templates("particle"):
            lj = pmb.get_lj_parameters(core_name, name)
            if lj:
                val = lj["offset"].to("reduced_length").magnitude
                if val > max_target:
                    max_target = val
        return max_target

    def test_final_offsets_match_lj_targets(self):
        """
        Unit test: after relax_nanoparticle_overlaps the ESPResSo LJ offset for
        core-core and core-ion pairs must equal the Lorentz-Berthelot targets
        returned by pmb.get_lj_parameters.

        Notes:
            - relax_espresso_system is mocked to isolate the LJ-setting logic
              from actual MD integration.
        """
        pmb = self._build_pmb()
        pmb.setup_lj_interactions(espresso_system=espresso_system)
        type_map  = pmb.get_type_map()
        core_type = type_map["np_core"]
        ion_type  = type_map["np_ion"]
        target_cc = pmb.get_lj_parameters("np_core", "np_core")["offset"].to("reduced_length").magnitude
        target_ci = pmb.get_lj_parameters("np_core", "np_ion")["offset"].to("reduced_length").magnitude

        with patch("pyMBE.lib.handy_functions.relax_espresso_system"):
            nanoparticle_tools.relax_nanoparticle_overlaps(
                espresso_system=espresso_system,
                pmb=pmb,
                nanoparticle_name="np",
                seed=42,
            )

        params_cc = espresso_system.non_bonded_inter[core_type, core_type].lennard_jones.get_params()
        params_ci = espresso_system.non_bonded_inter[core_type, ion_type].lennard_jones.get_params()
        self.assertAlmostEqual(params_cc["offset"], target_cc)
        self.assertAlmostEqual(params_ci["offset"], target_ci)

    def test_core_ion_offset_never_exceeds_per_pair_target(self):
        """
        Unit test: during the growth loop the core-ion offset must be clamped at
        its per-pair Lorentz-Berthelot target even while the global counter climbs
        toward the core-core target.  The last relax call (final pass) must receive
        the exact per-pair target.

        Notes:
            - relax_espresso_system is mocked with a side_effect that records the
              core-ion LJ offset in ESPResSo at each call.
            - The per-pair target is derived from pmb.get_lj_parameters, not
              hardcoded.
        """
        pmb = self._build_pmb()
        pmb.setup_lj_interactions(espresso_system=espresso_system)
        type_map  = pmb.get_type_map()
        core_type = type_map["np_core"]
        ion_type  = type_map["np_ion"]
        target_ci = pmb.get_lj_parameters("np_core", "np_ion")["offset"].to("reduced_length").magnitude

        captured_ci = []

        def capture(*_args, **_kwargs):
            p = espresso_system.non_bonded_inter[core_type, ion_type].lennard_jones.get_params()
            captured_ci.append(p["offset"])

        with patch("pyMBE.lib.handy_functions.relax_espresso_system", side_effect=capture):
            nanoparticle_tools.relax_nanoparticle_overlaps(
                espresso_system=espresso_system,
                pmb=pmb,
                nanoparticle_name="np",
                seed=42,
                delta_offset=0.5,
            )

        self.assertGreater(len(captured_ci), 0)
        for offset in captured_ci:
            self.assertLessEqual(offset, target_ci + 1e-10,
                                 msg=f"Core-ion offset {offset:.4f} exceeded per-pair target {target_ci:.4f}")
        self.assertAlmostEqual(captured_ci[-1], target_ci,
                               msg="Final relax call must receive the exact per-pair target offset")

    def test_sigma_zero_site_particles_not_touched(self):
        """
        Unit test: particles with sigma=0 must be excluded from the LJ growth
        loop.  After the call, the ESPResSo LJ sigma for every core-site type
        pair must remain 0.

        Notes:
            - Both protonated and deprotonated state types of ``np_site`` are
              found via the type_map (keys containing ``np_site``).
            - relax_espresso_system is mocked to isolate the LJ-setting logic.
        """
        pmb = self._build_pmb()
        pmb.setup_lj_interactions(espresso_system=espresso_system)
        type_map   = pmb.get_type_map()
        core_type  = type_map["np_core"]
        site_types = [t for k, t in type_map.items() if "np_site" in k]

        with patch("pyMBE.lib.handy_functions.relax_espresso_system"):
            nanoparticle_tools.relax_nanoparticle_overlaps(
                espresso_system=espresso_system,
                pmb=pmb,
                nanoparticle_name="np",
                seed=42,
            )

        for site_type in site_types:
            params = espresso_system.non_bonded_inter[core_type, site_type].lennard_jones.get_params()
            self.assertAlmostEqual(params["sigma"], 0.0,
                                   msg=f"Core-site LJ sigma must stay 0 for site ESPResSo type {site_type}")

    def test_thermostat_off_after_completion(self):
        """
        Unit test: the Langevin thermostat must be OFF when
        relax_nanoparticle_overlaps returns, consistent with the documented
        behavior of relax_espresso_system.

        Notes:
            - relax_espresso_system is NOT mocked: the real function must run so
              the thermostat state can be inspected.
            - delta_offset is set to the full max_target so a single loop
              iteration plus the final pass (2 MD runs) keeps the test fast.
        """
        pmb = self._build_pmb()
        espresso_system.time_step = 0.001
        espresso_system.cell_system.skin = 0.4
        pmb.setup_lj_interactions(espresso_system=espresso_system)
        espresso_system.thermostat.set_langevin(kT=1.0, gamma=0.1, seed=42)
        max_target = self._get_max_target_offset(pmb, "np_core")

        nanoparticle_tools.relax_nanoparticle_overlaps(
            espresso_system=espresso_system,
            pmb=pmb,
            nanoparticle_name="np",
            seed=42,
            delta_offset=max_target,   # single loop iteration + final pass
        )

        self.assertIsNone(espresso_system.thermostat.kT,
                          "Thermostat must be OFF after relax_nanoparticle_overlaps")

    def test_no_relax_calls_when_no_valid_lj_pairs(self):
        """
        Unit test: when every particle (including the core) has sigma=0,
        target_lj is empty and the function must return early without ever
        calling relax_espresso_system.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_particle(
            name="core_zero",
            z=0,
            sigma=0.0 * pmb.units("reduced_length"),
            epsilon=1.0 * pmb.units("reduced_energy"),
        )
        pmb.define_particle(
            name="site_zero",
            acidity="acidic",
            pka=4.0,
            sigma=0.0 * pmb.units("reduced_length"),
            epsilon=1.0 * pmb.units("reduced_energy"),
        )
        pmb.define_nanoparticle(
            name="np_zero",
            core_particle_name="core_zero",
            total_number_of_sites=2,
            primary_site_particle_name="site_zero",
            fraction_primary_sites=1.0,
            number_of_patches_of_primary_sites=1,
            secondary_site_particle_name=None,
        )
        pmb.setup_lj_interactions(espresso_system=espresso_system)

        with patch("pyMBE.lib.handy_functions.relax_espresso_system") as mock_relax:
            nanoparticle_tools.relax_nanoparticle_overlaps(
                espresso_system=espresso_system,
                pmb=pmb,
                nanoparticle_name="np_zero",
                seed=42,
            )

        mock_relax.assert_not_called()

    def test_verbose_prints_one_line_per_loop_iteration(self):
        """
        Unit test: with verbose=True a progress line containing "offset" must be
        printed for each iteration of the growth loop, no more and no less.

        Notes:
            - The expected line count is derived from max_target and delta_offset
              so the test stays correct if _build_pmb changes the particle offset.
            - relax_espresso_system is mocked to prevent MD logging from
              contaminating the captured output.
        """
        import io, sys
        pmb = self._build_pmb()
        pmb.setup_lj_interactions(espresso_system=espresso_system)
        delta_offset = 0.5
        max_target = self._get_max_target_offset(pmb, "np_core")
        expected_iterations = math.ceil(max_target / delta_offset)

        captured = io.StringIO()
        with patch("pyMBE.lib.handy_functions.relax_espresso_system"):
            sys.stdout = captured
            try:
                nanoparticle_tools.relax_nanoparticle_overlaps(
                    espresso_system=espresso_system,
                    pmb=pmb,
                    nanoparticle_name="np",
                    seed=42,
                    delta_offset=delta_offset,
                    verbose=True,
                )
            finally:
                sys.stdout = sys.__stdout__

        lines = [l for l in captured.getvalue().splitlines() if l.strip()]
        self.assertEqual(len(lines), expected_iterations,
                         f"Expected {expected_iterations} progress lines, got {len(lines)}")
        for line in lines:
            self.assertIn("offset", line.lower())


if __name__ == "__main__":
    ut.main()
