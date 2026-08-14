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

# Import pyMBE and other libraries
import pyMBE
import numpy as np
import unittest as ut
import espressomd

# Create an instance of the ESPResSo system
box_l=[10]*3
espresso_system = espressomd.System(box_l=box_l)


class Test(ut.TestCase):
    def setUp(self):
        espresso_system.part.clear()
        espresso_system.bonded_inter.clear()

    def define_templates(self, pmb):
        pmb.define_particle(name='A',
            z=0,
            sigma=0.4*pmb.units.nm,
            epsilon=1*pmb.units('reduced_energy'))

        pmb.define_particle(name='B',
            z=0,
            sigma=0.4*pmb.units.nm,
            epsilon=1*pmb.units('reduced_energy'))

        pmb.define_particle(name='C',
            z=0,
            sigma=0.4*pmb.units.nm,
            epsilon=1*pmb.units('reduced_energy'))

        self.harmonic_angle_params = {
            'k': 50 * pmb.units('reduced_energy'),
            'phi_0': np.pi * pmb.units(''),
        }

        self.cosine_angle_params = {
            'k': 30 * pmb.units('reduced_energy'),
            'phi_0': np.pi / 2 * pmb.units(''),
        }

        self.harmonic_cosine_angle_params = {
            'k': 40 * pmb.units('reduced_energy'),
            'phi_0': np.pi / 3 * pmb.units(''),
        }

        self.harmonic_bond_params = {
            'r_0': 0.4 * pmb.units.nm,
            'k': 400 * pmb.units('reduced_energy / reduced_length**2'),
        }

    def get_angle_object(self, central_particle_id):
        """
        Returns the angle interaction object stored in ESPResSo for a given central particle.
        Angle bonds are added to the central particle in ESPResSo.
        """
        bonds = espresso_system.part.by_id(central_particle_id).bonds
        for bond_tuple in bonds:
            bond_obj = bond_tuple[0] # Testing only one angle on the central particle
            # Angle interactions involve 2 partner particles (the two sides)
            if len(bond_tuple) == 3:
                return bond_obj
        return None

    def check_angle_setup(self, angle_object, input_parameters, angle_type):
        """
        Checks that pyMBE sets up an angle interaction object correctly.

        Args:
            angle_object: instance of an ESPResSo angle interaction object.
            input_parameters('dict'): dictionary with the parameters for the angle potential.
            angle_type('str'): label identifying the angle type.
        """
        # Check that pyMBE stores the correct type of angle
        type_map = {
            "harmonic": "AngleHarmonic",
            "cosine": "AngleCosine",
            "harmonic_cosine": "AngleCossquare",
        }
        self.assertIn(type_map[angle_type].lower(), str(type(angle_object)).lower(),
                      msg=f"pyMBE does not store the correct type of angle interaction, expected {type_map[angle_type]}")
        # Check that pyMBE defines the angle with the correct parameters
        angle_params = angle_object.get_params()
        np.testing.assert_almost_equal(actual=angle_params['bend'],
                                       desired=input_parameters['k'].m_as('reduced_energy'),
                                       verbose=True)
        np.testing.assert_almost_equal(actual=angle_params['phi0'],
                                       desired=float(input_parameters['phi_0'].magnitude),
                                       verbose=True)

    def test_angle_setup_harmonic(self):
        """
        Unit test to check the setup of harmonic angle potentials in pyMBE.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        # Define bonds between A-B and B-C
        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C']])

        # Define harmonic angle for triplet A-B-C
        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=self.harmonic_angle_params,
                         particle_triplets=[('A', 'B', 'C')])

        # Create three particles: A, B, C
        pid_A = pmb.create_particle(name="A",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_B = pmb.create_particle(name="B",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_C = pmb.create_particle(name="C",
                                    box_l=box_l,
                                    number_of_particles=1)

        # Create bonds: A-B and B-C
        pmb.create_bond(particle_id1=pid_A[0],
                        particle_id2=pid_B[0])
        pmb.create_bond(particle_id1=pid_B[0],
                        particle_id2=pid_C[0])

        # Create the angle: A-B-C (B is central)
        
        pmb.create_angular_potential(particle_id1=pid_A[0],
                         particle_id2=pid_B[0],
                         particle_id3=pid_C[0])
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        angle_object = self.get_angle_object(central_particle_id=pid_B[0])
        self.assertIsNotNone(angle_object, "No angle object found on the central particle")

        self.check_angle_setup(angle_object=angle_object,
                               input_parameters=self.harmonic_angle_params,
                               angle_type="harmonic")

        # Clean-up
        for inst_id in pid_A + pid_B + pid_C:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="angle")
        pmb.db.delete_templates(pmb_type="bond")

    def test_angle_setup_cosine(self):
        """
        Unit test to check the setup of cosine angle potentials in pyMBE.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        # Define bonds between A-B
        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B']])

        # Define cosine angle for triplet A-B-A
        pmb.define_angular_potential(angle_type="cosine",
                         angle_parameters=self.cosine_angle_params,
                         particle_triplets=[('A', 'B', 'A')])

        # Create three particles
        pid_A1 = pmb.create_particle(name="A",
                                     box_l=box_l,
                                     number_of_particles=1)
        pid_B = pmb.create_particle(name="B",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_A2 = pmb.create_particle(name="A",
                                     box_l=box_l,
                                     number_of_particles=1)

        # Create bonds: A1-B and B-A2
        pmb.create_bond(particle_id1=pid_A1[0],
                        particle_id2=pid_B[0],
                       )
        pmb.create_bond(particle_id1=pid_B[0],
                        particle_id2=pid_A2[0],
                        )

        pmb.create_angular_potential(particle_id1=pid_A1[0],
                         particle_id2=pid_B[0],
                         particle_id3=pid_A2[0],
                         )
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()

        angle_object = self.get_angle_object(central_particle_id=pid_B[0])
        self.assertIsNotNone(angle_object, "No angle object found on the central particle")

        self.check_angle_setup(angle_object=angle_object,
                               input_parameters=self.cosine_angle_params,
                               angle_type="cosine")

        # Clean-up
        for inst_id in pid_A1 + pid_B + pid_A2:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="angle")
        pmb.db.delete_templates(pmb_type="bond")

    def test_angle_setup_harmonic_cosine(self):
        """
        Unit test to check the setup of harmonic-cosine angle potentials in pyMBE.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C']])

        pmb.define_angular_potential(angle_type="harmonic_cosine",
                         angle_parameters=self.harmonic_cosine_angle_params,
                         particle_triplets=[('A', 'B', 'C')])

        pid_A = pmb.create_particle(name="A",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_B = pmb.create_particle(name="B",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_C = pmb.create_particle(name="C",
                                    box_l=box_l,
                                    number_of_particles=1)

        pmb.create_bond(particle_id1=pid_A[0],
                        particle_id2=pid_B[0])
        pmb.create_bond(particle_id1=pid_B[0],
                        particle_id2=pid_C[0])

        pmb.create_angular_potential(particle_id1=pid_A[0],
                         particle_id2=pid_B[0],
                         particle_id3=pid_C[0])
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        angle_object = self.get_angle_object(central_particle_id=pid_B[0])
        self.assertIsNotNone(angle_object, "No angle object found on the central particle")

        self.check_angle_setup(angle_object=angle_object,
                               input_parameters=self.harmonic_cosine_angle_params,
                               angle_type="harmonic_cosine")

        for inst_id in pid_A + pid_B + pid_C:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="angle")
        pmb.db.delete_templates(pmb_type="bond")

    def test_angle_instance_reuses_cached_object(self):
        """
        Test that cached ESPResSo angle interactions are reused for the same template.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=self.harmonic_angle_params,
                         particle_triplets=[('A', 'B', 'C')])

        angle_template = pmb.get_angle_template(side_name1="A",
                                                central_name="B",
                                                side_name2="C")

        pmb.set_simulation_engine(espresso_system)
        first_angle_object = pmb.simulation_engine._get_angle_instance(angle_template=angle_template)
        second_angle_object = pmb.simulation_engine._get_angle_instance(angle_template=angle_template)

        self.assertIs(first_angle_object, second_angle_object)
        self.assertEqual(len(pmb.db.espresso_angle_instances), 1)

        pmb.db.delete_templates(pmb_type="angle")

    def test_default_angle(self):
        """
        Unit test to check the setup of default angle potentials.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        # Define bonds between A-B and B-C
        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C']])

        # Define only a default angle (no specific triplet)
        pmb.define_default_angular_potential(angle_type="harmonic",
                                 angle_parameters=self.harmonic_angle_params)

        # Create three particles
        pid_A = pmb.create_particle(name="A",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_B = pmb.create_particle(name="B",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_C = pmb.create_particle(name="C",
                                    box_l=box_l,
                                    number_of_particles=1)

        # Create bonds: A-B and B-C
        pmb.create_bond(particle_id1=pid_A[0],
                        particle_id2=pid_B[0])
        pmb.create_bond(particle_id1=pid_B[0],
                        particle_id2=pid_C[0])

        # Create angle using default (no specific A-B-C template exists)
        pmb.create_angular_potential(particle_id1=pid_A[0],
                         particle_id2=pid_B[0],
                         particle_id3=pid_C[0],
                         use_default_angle=True)
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        angle_object = self.get_angle_object(central_particle_id=pid_B[0])
        self.assertIsNotNone(angle_object, "No angle object found on the central particle")

        self.check_angle_setup(angle_object=angle_object,
                               input_parameters=self.harmonic_angle_params,
                               angle_type="harmonic")

        # Clean-up
        for inst_id in pid_A + pid_B + pid_C:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="angle")
        pmb.db.delete_templates(pmb_type="bond")

    def test_specific_angle_over_default(self):
        """
        Test that a specific angle template is preferred over the default.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        # Define bonds between A-B and B-C
        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C']])

        # Define default angle with cosine parameters
        pmb.define_default_angular_potential(angle_type="cosine",
                                 angle_parameters=self.cosine_angle_params)

        # Define specific angle for A-B-C with harmonic parameters
        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=self.harmonic_angle_params,
                         particle_triplets=[('A', 'B', 'C')])

        pid_A = pmb.create_particle(name="A",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_B = pmb.create_particle(name="B",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_C = pmb.create_particle(name="C",
                                    box_l=box_l,
                                    number_of_particles=1)

        # Create bonds: A-B and B-C
        pmb.create_bond(particle_id1=pid_A[0],
                        particle_id2=pid_B[0])
        pmb.create_bond(particle_id1=pid_B[0],
                        particle_id2=pid_C[0])

        # Should use the specific A-B-C template, not the default
        pmb.create_angular_potential(particle_id1=pid_A[0],
                         particle_id2=pid_B[0],
                         particle_id3=pid_C[0],
                         use_default_angle=True)
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        angle_object = self.get_angle_object(central_particle_id=pid_B[0])
        self.assertIsNotNone(angle_object)

        # Should be harmonic (specific), not cosine (default)
        self.check_angle_setup(angle_object=angle_object,
                               input_parameters=self.harmonic_angle_params,
                               angle_type="harmonic")

        # Clean-up
        for inst_id in pid_A + pid_B + pid_C:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="angle")
        pmb.db.delete_templates(pmb_type="bond")

    def test_angle_raised_exceptions(self):
        """
        Unit test to check that angle-related methods raise the correct exceptions.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        for callback in [pmb.define_angular_potential, pmb.define_default_angular_potential]:
            with self.subTest(msg=f'using method {callback.__qualname__}()'):
                self.check_angle_exceptions(callback, pmb)

    def check_angle_exceptions(self, callback, pmb):
        # Check exception for unknown angle type
        angle_type = 'quartic'
        angle_params = {'k': 50 * pmb.units('reduced_energy'),
                        'phi_0': np.pi * pmb.units('')}

        input_parameters = {"angle_type": angle_type, "angle_parameters": angle_params}
        if callback == pmb.define_angular_potential:
            input_parameters["particle_triplets"] = [('A', 'B', 'C')]

        np.testing.assert_raises(NotImplementedError, callback, **input_parameters)

        # Check exception for missing k
        angle_type = 'harmonic'
        angle_params = {'phi_0': np.pi * pmb.units('')}

        input_parameters = {"angle_type": angle_type, "angle_parameters": angle_params}
        if callback == pmb.define_angular_potential:
            input_parameters["particle_triplets"] = [('A', 'B', 'C')]

        np.testing.assert_raises(ValueError, callback, **input_parameters)

        # Check exception for missing phi_0
        angle_type = 'harmonic'
        angle_params = {'k': 50 * pmb.units('reduced_energy')}

        input_parameters = {"angle_type": angle_type, "angle_parameters": angle_params}
        if callback == pmb.define_angular_potential:
            input_parameters["particle_triplets"] = [('A', 'B', 'C')]

        np.testing.assert_raises(ValueError, callback, **input_parameters)

        # Check exception for duplicate triplet in define_angular_potential
        if callback == pmb.define_angular_potential:
            test = {"angle_type": "harmonic",
                    "angle_parameters": self.harmonic_angle_params,
                    "particle_triplets": [("A", "B", "A"), ("A", "B", "A")]}

            np.testing.assert_raises(RuntimeError,
                                     pmb.define_angular_potential,
                                     **test)

    def test_sanity_get_angle_template(self):
        """
        Tests that get_angle_template raises ValueError when no template is found.
        """
        pmb = pyMBE.pymbe_library(seed=51)
        inputs = {"side_name1": "X",
                  "central_name": "Y",
                  "side_name2": "Z"}
        np.testing.assert_raises(ValueError,
                                 pmb.get_angle_template,
                                 **inputs)

    def test_canonical_angle_name(self):
        """
        Test that angle names are canonical (side particles sorted alphabetically).
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        # Define angle with C-B-A order — should produce same template as A-B-C
        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=self.harmonic_angle_params,
                         particle_triplets=[('C', 'B', 'A')])

        # The canonical name should be A-B-C (sides sorted, central stays in middle)
        tpl = pmb.get_angle_template(side_name1="A",
                                     central_name="B",
                                     side_name2="C")
        self.assertEqual(tpl.name, "A-B-C")

        # Also retrievable with reversed side order
        tpl2 = pmb.get_angle_template(side_name1="C",
                                      central_name="B",
                                      side_name2="A")
        self.assertEqual(tpl2.name, "A-B-C")

        pmb.db.delete_templates(pmb_type="angle")

    def test_angle_without_bonds_raises_error(self):
        """
        Test that creating an angle without bonds between the particles raises a ValueError.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        # Define angle template but no bond templates
        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=self.harmonic_angle_params,
                         particle_triplets=[('A', 'B', 'C')])

        # Create three particles without any bonds between them
        pid_A = pmb.create_particle(name="A",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_B = pmb.create_particle(name="B",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_C = pmb.create_particle(name="C",
                                    box_l=box_l,
                                    number_of_particles=1)
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        # Attempting to create angle without bonds should raise ValueError
        with self.assertRaises(ValueError, msg="create_angular_potential should raise ValueError when no bonds exist"):
            pmb.create_angular_potential(particle_id1=pid_A[0],
                             particle_id2=pid_B[0],
                             particle_id3=pid_C[0])

        # Clean-up
        for inst_id in pid_A + pid_B + pid_C:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="angle")

    def test_angle_with_partial_bonds_raises_error(self):
        """
        Test that creating an angle with only one bond (missing the other) raises a ValueError.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        # Define bond only between A-B, not B-C
        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B']])

        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=self.harmonic_angle_params,
                         particle_triplets=[('A', 'B', 'C')])

        pid_A = pmb.create_particle(name="A",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_B = pmb.create_particle(name="B",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_C = pmb.create_particle(name="C",
                                    box_l=box_l,
                                    number_of_particles=1)

        # Create only the A-B bond, not B-C
        pmb.create_bond(particle_id1=pid_A[0],
                        particle_id2=pid_B[0])
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        # Should raise ValueError because B-C bond is missing
        with self.assertRaises(ValueError, msg="create_angular_potential should raise ValueError when a bond is missing"):
            pmb.create_angular_potential(particle_id1=pid_A[0],
                             particle_id2=pid_B[0],
                             particle_id3=pid_C[0])

        # Clean-up
        for inst_id in pid_A + pid_B + pid_C:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="angle")
        pmb.db.delete_templates(pmb_type="bond")

    def test_generate_angles_for_entity_without_particles_returns(self):
        """
        Test that auto-generating angles for an empty entity returns without changes.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb._generate_angles_for_entity(
                                        entity_id=999,
                                        entity_id_col="residue_id")

        self.assertEqual(len(pmb.get_instances_df(pmb_type="angle")), 0)

    def test_generate_angles_for_entity_creates_angle_from_bonds(self):
        """
        Test that auto-generating angles creates one angle from a bonded triplet.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C']])
        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=self.harmonic_angle_params,
                         particle_triplets=[('A', 'B', 'C')])

        pid_A = pmb.create_particle(name="A",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_B = pmb.create_particle(name="B",
                                     box_l=box_l,
                                    number_of_particles=1)
        pid_C = pmb.create_particle(name="C",
                                   box_l=box_l,
                                    number_of_particles=1)

        for particle_id in (pid_A[0], pid_B[0], pid_C[0]):
            pmb.db._update_instance(instance_id=particle_id,
                                    pmb_type="particle",
                                    attribute="residue_id",
                                    value=0)

        pmb.create_bond(particle_id1=pid_A[0],
                        particle_id2=pid_B[0])
        pmb.create_bond(particle_id1=pid_B[0],
                        particle_id2=pid_C[0])
        

        pmb._generate_angles_for_entity(
                                        entity_id=0,
                                        entity_id_col="residue_id")
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        angle_object = self.get_angle_object(central_particle_id=pid_B[0])
        self.assertIsNotNone(angle_object, "No angle object found on the central particle")
        self.assertEqual(len(pmb.get_instances_df(pmb_type="angle")), 1)

        for inst_id in pid_A + pid_B + pid_C:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="angle")
        pmb.db.delete_templates(pmb_type="bond")

    def test_generate_angles_for_entity_skips_missing_templates(self):
        """
        Test that auto-generating angles skips bonded triplets without a matching angle template.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C']])

        pid_A = pmb.create_particle(name="A",
                                    box_l=box_l,
                                    number_of_particles=1)
        pid_B = pmb.create_particle(name="B",
                                     box_l=box_l,
                                    number_of_particles=1)
        pid_C = pmb.create_particle(name="C",
                                     box_l=box_l,
                                    number_of_particles=1)

        for particle_id in (pid_A[0], pid_B[0], pid_C[0]):
            pmb.db._update_instance(instance_id=particle_id,
                                    pmb_type="particle",
                                    attribute="residue_id",
                                    value=0)

        pmb.create_bond(particle_id1=pid_A[0],
                        particle_id2=pid_B[0])
        pmb.create_bond(particle_id1=pid_B[0],
                        particle_id2=pid_C[0])
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        pmb._generate_angles_for_entity(
                                        entity_id=0,
                                        entity_id_col="residue_id")

        self.assertIsNone(self.get_angle_object(central_particle_id=pid_B[0]))
        self.assertEqual(len(pmb.get_instances_df(pmb_type="angle")), 0)

        for inst_id in pid_A + pid_B + pid_C:
            pmb.delete_instances_in_system(instance_id=inst_id,
                                           pmb_type="particle")
        pmb.db.delete_templates(pmb_type="bond")

    def test_create_residue_with_gen_angle_generates_angles(self):
        """
        Test that create_residue generates angles automatically when gen_angle is enabled.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C']])
        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=self.harmonic_angle_params,
                         particle_triplets=[('A', 'B', 'C')])
        pmb.define_residue(name="R_angle",
                           central_bead="B",
                           side_chains=["A", "C"])

        residue_id = pmb.create_residue(name="R_angle",
                                        box_l=box_l,
                                        gen_angle=True)
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        particle_ids = pmb.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                              attribute="residue_id",
                                                              value=residue_id)
        central_particle_id = next(pid for pid in particle_ids
                                   if pmb.db.get_instance(pmb_type="particle",
                                                          instance_id=pid).name == "B")

        angle_object = self.get_angle_object(central_particle_id=central_particle_id)
        self.assertIsNotNone(angle_object, "No angle object found on the central particle")
        self.assertEqual(len(pmb.get_instances_df(pmb_type="angle")), 1)

        pmb.delete_instances_in_system(instance_id=residue_id,
                                       pmb_type="residue")
        pmb.db.delete_templates(pmb_type="angle")
        pmb.db.delete_templates(pmb_type="bond")

    def test_create_molecule_with_gen_angle_generates_angles(self):
        """
        Test that create_molecule generates backbone angles automatically when gen_angle is enabled.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        self.define_templates(pmb)

        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=self.harmonic_bond_params,
                        particle_pairs=[['A', 'B'], ['B', 'C']])
        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=self.harmonic_angle_params,
                         particle_triplets=[('A', 'B', 'C')])
        pmb.define_residue(name="R_A",
                           central_bead="A",
                           side_chains=[])
        pmb.define_residue(name="R_B",
                           central_bead="B",
                           side_chains=[])
        pmb.define_residue(name="R_C",
                           central_bead="C",
                           side_chains=[])
        pmb.define_molecule(name="M_angle",
                            residue_list=["R_A", "R_B", "R_C"])

        molecule_ids = pmb.create_molecule(name="M_angle",
                                           number_of_molecules=1,
                                            box_l=box_l,
                                           gen_angle=True)
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        particle_ids = pmb.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                              attribute="molecule_id",
                                                              value=molecule_ids[0])
        central_particle_id = next(pid for pid in particle_ids
                                   if pmb.db.get_instance(pmb_type="particle",
                                                          instance_id=pid).name == "B")

        angle_object = self.get_angle_object(central_particle_id=central_particle_id)
        self.assertIsNotNone(angle_object, "No angle object found on the central particle")
        self.assertEqual(len(pmb.get_instances_df(pmb_type="angle")), 1)

        pmb.delete_instances_in_system(instance_id=molecule_ids[0],
                                       pmb_type="molecule")
        pmb.db.delete_templates(pmb_type="angle")
        pmb.db.delete_templates(pmb_type="bond")


if __name__ == '__main__':
    ut.main()
