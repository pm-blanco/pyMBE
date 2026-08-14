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

import unittest as ut
from unittest.mock import Mock
from pyMBE import pymbe_library
from pyMBE.simulation_builder.lammps_engine import LammpsSimulation
import espressomd

box_l=[10]*3
espresso_system=espressomd.System(box_l=box_l)
pmb=pymbe_library(seed=32)
lammps_engine=LammpsSimulation()

class Test(ut.TestCase):
    """Unit tests for simulation-engine selection and dispatch."""
    def test_dummy_engine(self):
        """Checks that operations fail when no simulation engine is configured."""
        with self.assertRaises(RuntimeError):
            pmb.setup_lj_interactions()
    def test_missing_simulation_engine(self):
        """Checks that selecting ``None`` as the simulation engine is rejected."""
        with self.assertRaises(ValueError):
            pmb.set_simulation_engine(simulation_engine=None,box_l=box_l)

    def test_lammps_engine_initialization(self):
        """Checks that a LAMMPS engine stores its initialization parameters."""
        database = object()
        backend = object()
        units = object()
        engine = LammpsSimulation(
            box_l=box_l,
            db=database,
            lammps=backend,
            units=units,
            kT=2.0,
            Kw=3.0,
            seed=4,
        )

        self.assertIs(engine.db, database)
        self.assertIs(engine.lammps, backend)
        self.assertIs(engine.units, units)
        self.assertEqual(engine.box_l, box_l)
        self.assertEqual(engine.kT, 2.0)
        self.assertEqual(engine.Kw, 3.0)
        self.assertEqual(engine.seed, 4)

    def test_lammps_engine_not_implemented(self):
        """Checks that unimplemented LAMMPS operations raise ``NotImplementedError``."""
        methods = [
            lammps_engine.display_commands,
            lammps_engine.setup_lj_interactions,
            lammps_engine.setup_langevin,
            lammps_engine.setup_cpH,
            lammps_engine.run_simulation,
        ]
        for method in methods:
            with self.subTest(method=method.__name__):
                with self.assertRaises(NotImplementedError):
                    method()
        with self.assertRaises(NotImplementedError):
            lammps_engine.some_future_method()


    def test_lammps_engine_selection(self):
        """Checks that a LAMMPS-compatible backend is wrapped by ``LammpsSimulation``."""
        class FakeLammps:
            def display_commands(self):
                """Provides the command-display method required by the protocol."""
                pass # pragma: no cover
            def setup_lj_interactions(self):
                """Provides the Lennard-Jones setup method required by the protocol."""
                pass # pragma: no cover
            def setup_langevin(self):
                """Provides the Langevin setup method required by the protocol."""
                pass # pragma: no cover
            def setup_cpH(self):
                """Provides the constant-pH setup method required by the protocol."""
                pass # pragma: no cover
            def run_simulation(self):
                """Provides the simulation-run method required by the protocol."""
                pass # pragma: no cover

        local_pmb = pymbe_library(seed=32)
        fake_lammps = FakeLammps()
        local_pmb.set_simulation_engine(fake_lammps, box_l=box_l)

        self.assertIsInstance(local_pmb.simulation_engine, LammpsSimulation)
        self.assertIs(local_pmb.simulation_engine.lammps, fake_lammps)
        self.assertEqual(local_pmb.simulation_engine.box_l, box_l)

    def test_espresso_engine_helpers_and_volume_validation(self):
        """Checks ESPResSo helper methods and volume-validation behavior."""
        local_pmb = pymbe_library(seed=32)
        local_pmb.set_simulation_engine(espresso_system)
        engine = local_pmb.simulation_engine
        particle = espresso_system.part.add(pos=[1, 2, 3])
        self.assertEqual(list(engine._get_particle_pos_espresso(particle.id)), [1, 2, 3])
        self.assertListEqual(list(engine.get_box_side_length()), box_l)
        espresso_system.part.clear()

        engine.db = Mock()
        engine.db._get_instances_df.return_value.index.size = 0
        engine.espresso_system = Mock()
        with self.assertRaises(ValueError):
            engine.change_volume_and_rescale_particles(0)
        engine.change_volume_and_rescale_particles(5, dir="xy")
        engine.espresso_system.change_volume_and_rescale_particles.assert_called_once_with(
            d_new=5, dir="xy")

    def test_espresso_engine_update_particle_id_without_engine(self):
        """Checks that updating IDs without an attached engine raises an error."""
        local_pmb = pymbe_library(seed=32)
        local_pmb.set_simulation_engine(espresso_system)
        local_pmb.simulation_engine.espresso_system = None
        with self.assertRaises(ValueError):
            local_pmb.simulation_engine.update_particle_id(old_pid=1, new_pid=2)

    def test_espresso_engine_update_particle_id(self):
        """Checks that particle IDs are propagated to both bond endpoints."""
        local_pmb = pymbe_library(seed=32)
        local_pmb.set_simulation_engine(espresso_system)
        engine = local_pmb.simulation_engine
        engine.db = Mock()
        engine.db._find_instance_ids_by_attribute.side_effect = [[10], [11]]

        engine.update_particle_id(old_pid=2, new_pid=8)

        self.assertEqual(engine.db._update_instance.call_count, 3)
        calls = engine.db._update_instance.call_args_list
        self.assertEqual(calls[0].kwargs["value"], 8)
        self.assertEqual(calls[1].kwargs["instance_id"], 10)
        self.assertEqual(calls[2].kwargs["instance_id"], 11)

    def test_espresso_engine_add_instances_id_collision(self):
        """Checks that colliding particle IDs are remapped before insertion."""
        local_pmb = pymbe_library(seed=32)
        local_pmb.set_simulation_engine(espresso_system)
        engine = local_pmb.simulation_engine
        engine.db = Mock()
        engine.db._find_instance_ids_by_attribute.side_effect = [[1], [2], [], []]
        engine._get_particle_ids_in_espresso = Mock(return_value=[1, 2])
        engine.update_particle_id = Mock()
        engine._add_particle = Mock()

        with self.assertWarns(UserWarning):
            engine.add_instances_to_engine()

        engine.update_particle_id.assert_called_once_with(old_pid=1, new_pid=3)
        engine._add_particle.assert_called_once_with(3)

    def test_espresso_engine_no_particles_to_add(self):
        """Checks that adding an empty particle set raises a runtime error."""
        local_pmb = pymbe_library(seed=32)
        local_pmb.set_simulation_engine(espresso_system)
        local_pmb.simulation_engine.db = Mock()
        local_pmb.simulation_engine.db._find_instance_ids_by_attribute.return_value = []
        with self.assertRaises(RuntimeError):
            local_pmb.simulation_engine.add_instances_to_engine()

if __name__=='__main__':
    ut.main()
