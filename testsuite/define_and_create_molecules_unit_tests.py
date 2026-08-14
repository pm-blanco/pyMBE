#
# Copyright (C) 2024-2026 pyMBE-dev team
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

import pyMBE
import numpy as np
import espressomd
import unittest as ut


# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

# The unit tests for define_particle are in lj_tests.py and set_particle_acidity

# Define particles, residues, and molecules for testing
particle_parameters={"S1":{"name":"S1",
                           "sigma":1*pmb.units.nm,
                           "offset":0.5*pmb.units.nm,
                           "epsilon":1.0*pmb.units.reduced_energy,
                           "z":1},
                     "S2":{"name":"S2",
                           "sigma":2*pmb.units.nm,
                           "epsilon":1.0*pmb.units.reduced_energy,
                           "offset":1.5*pmb.units.nm,
                           "z": 1},
                     "S3":{"name":"S3",
                           "sigma":3*pmb.units.nm,
                           "epsilon":1.0*pmb.units.reduced_energy,
                           "offset":2.5*pmb.units.nm,
                           "z":2}}
for particle_set in particle_parameters.values():
    pmb.define_particle(**particle_set)

residue_parameters={"R1":{"name": "R1",
                        "central_bead": "S1",
                        "side_chains": []},
                    "R2":{"name": "R2",
                        "central_bead": "S1",
                        "side_chains": ["S2","S3"]},
                    "R3":{"name": "R3",
                        "central_bead": "S2",
                        "side_chains": ["R2"]}}

for parameter_set in residue_parameters.values():
    pmb.define_residue(**parameter_set)

molecule_parameters={"M1":{"name": "M1",
                     "residue_list": []},
                    "M2":{"name": "M2",
                    "residue_list": ["R1","R2","R3"]}}

for parameter_set in molecule_parameters.values():
    pmb.define_molecule(**parameter_set)

# Create an instance of an espresso system
box_l=[10]*3
espresso_system=espressomd.System(box_l = box_l)
particle_positions=[[0,0,0],[1,1,1]]

bond_type = 'harmonic'
bond = {'r_0'    : 0.4*pmb.units.nm,
        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

pmb.define_default_bond(bond_type = bond_type,
                        bond_parameters = bond)
type_map=pmb.get_type_map()

class Test(ut.TestCase):
    def test_residue_definition(self):
        """
        Unit test: check that define_residue() stores the parameters correctly in the pyMBE database
        """
        for residue_name in residue_parameters.keys():
            input_parameters=residue_parameters[residue_name]
            res_tpl = pmb.db.get_template(pmb_type="residue", 
                                          name=residue_name)
            self.assertEqual(res_tpl.pmb_type, 
                             "residue")
            self.assertEqual(res_tpl.central_bead, 
                             input_parameters["central_bead"])
            self.assertEqual(res_tpl.side_chains, 
                             input_parameters["side_chains"])

    def test_molecule_definition(self):
        """
        Unit test: check that define_molecule() stores the parameters correctly in the pyMBE database
        """
        for molecule_name in molecule_parameters.keys():
            input_parameters=molecule_parameters[molecule_name]
            mol_tpl = pmb.db.get_template(pmb_type="molecule", 
                                          name=molecule_name)
            self.assertEqual(mol_tpl.pmb_type, 
                             "molecule")
            self.assertEqual(mol_tpl.residue_list, 
                             input_parameters["residue_list"])
            
    def test_create_and_delete_particles(self):
        """
        Docstring for test_create_and_delete_particles_residues_molecules
    
        """
        retval = pmb.create_particle(name="S1",
                                    box_l=box_l,
                                    number_of_particles=2,
                                    fix=[True,True,True],
                                    position=particle_positions)
        
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()

        self.assertListEqual(retval, [0, 1])

        particle_ids=pmb.get_particle_id_map(object_name="S1")["all"]
        
        for pid in particle_ids:
            particle=espresso_system.part.by_id(pid)
            self.assertEqual(first=particle.type, 
                             second=type_map["S1"])
            self.assertEqual(first=particle.q, 
                                        second=particle_parameters["S1"]["z"])
            self.assertListEqual(list1=list(particle.fix), 
                                 list2=[True]*3)
            self.assertListEqual(list1=list(particle.pos), 
                                 list2=particle_positions[pid])
        starting_number_of_particles=len(espresso_system.part.all())

        ### Tests that when number of particles is 0 or -1 no particles is created
        for number_of_particles in [0, -1]:
            retval = pmb.create_particle(name="S1",
                                        box_l=box_l,
                                        number_of_particles=number_of_particles)
            self.assertEqual(len(retval), 0)

        
        # If no particles have been created, only two particles should be in the system (from the previous test)
        self.assertEqual(first=len(espresso_system.part.all()), 
                                second=starting_number_of_particles)
        
        ### Returns ValueError as no particle with name S23 has been previously defined.
        with self.assertRaises(ValueError):
            pmb.create_particle(name="S23",
                                box_l=box_l,
                                number_of_particles=1)
            
        # If no particles have been created, only two particles should be in the system (from the previous test)
        self.assertEqual(first=len(espresso_system.part.all()), 
                                second=starting_number_of_particles)

        # Unit tests for delete particle
        starting_number_of_particles=len(espresso_system.part.all())
        starting_number_of_rows=len(pmb.get_instances_df(pmb_type="particle"))
        # This should delete one particle instance
        pmb.delete_instances_in_system(instance_id=0,
                                    pmb_type="particle")
        self.assertEqual(first=len(espresso_system.part.all()), 
                                second=starting_number_of_particles-1)
        particle_df = pmb.get_instances_df(pmb_type="particle")
        self.assertEqual(first=len(particle_df), 
                                second=starting_number_of_rows-1)

        # Delete the other particle instance to simplify the rest of the tests
        pmb.delete_instances_in_system(instance_id=1,
                                    pmb_type="particle")
        

        central_bead_position=[[0,0,0]]
        backbone_vector=np.array([1.,2.,3.])

        pmb.create_residue(name="R2",
                        box_l=box_l,
                        central_bead_position=central_bead_position,
                        backbone_vector=backbone_vector,
                        use_default_bond=True)

        particle_ids=pmb.get_particle_id_map(object_name="R2")["all"]

        pmb.add_instances_to_engine()
        # Check that the particle properties are correct
        for pid in particle_ids:
            particle=espresso_system.part.by_id(pid)
            particle_tpl = pmb.db.get_instance(pmb_type="particle",
                                            instance_id=pid)
            particle_name = particle_tpl.name
            self.assertEqual(first=particle.type, 
                            second=type_map[particle_name])
            self.assertEqual(first=particle.q, 
                            second=particle_parameters[particle_name]["z"])
            # Check that the position are correct
            # Central bead
            if particle_name == "S1":
                self.assertListEqual(list1=list(particle.pos), 
                                     list2=central_bead_position[0])
            else: # Side chains should be in positions perpendicular to the backbone vector
                self.assertAlmostEqual(first=np.dot(particle.pos,backbone_vector), 
                                        second=0, 
                                        places=10)
            # Check that particles have the correct residue id
            residue_id = particle_tpl.residue_id
            self.assertEqual(first=residue_id, 
                            second=0)

        # Check that particles are correctly bonded
        # Central bead S1 (id 0) should be bonded to S2 (id 1) and S3 (id 2)
        bonded_pairs=[]
        bond_df = pmb.get_instances_df(pmb_type="bond")
        for bond_index in bond_df.index:
            particle_id1= bond_df.loc[bond_index,"particle_id1"]
            particle_id2= bond_df.loc[bond_index,"particle_id2"]
            bonded_pair=frozenset([particle_id1,particle_id2])
            bonded_pairs.append(bonded_pair)
            bonded_in_espresso = False
            for pid in bonded_pair:
                for bond in espresso_system.part.by_id(pid).bonds[:]:
                    partner_id  = bond[1]
                    if partner_id in bonded_pair:
                        bonded_in_espresso=True
            self.assertEqual(first=bonded_in_espresso, 
                            second=True)

        self.assertEqual(first=frozenset(bonded_pairs), 
                        second=frozenset([frozenset([0,1]),frozenset([0,2])]))

        pmb.create_residue(name="R3",
                            box_l=box_l,
                            use_default_bond=True)

        particle_ids=pmb.get_particle_id_map(object_name="R3")["all"]
        pmb.add_instances_to_engine()
        # Check that the particle properties are correct
        for pid in particle_ids:
            particle=espresso_system.part.by_id(pid)
            particle_tpl = pmb.db.get_instance(pmb_type="particle",
                                                instance_id=pid)
            particle_name = particle_tpl.name
            self.assertEqual(first=particle.type, 
                             second=type_map[particle_name])
            self.assertEqual(first=particle.q, 
                             second=particle_parameters[particle_name]["z"])
            # Check that particles have the correct residue id
            residue_id = particle_tpl.residue_id
            self.assertEqual(first=residue_id, 
                             second=1)

        # Check that particles are correctly bonded, new bonds are:
        # Central bead S2 (id 3) should be bonded to R2 central bead S1 (id 4)
        # Central bead S1 (id 4) should be bonded to side chains S2 (id 5) and S3 (id 6)

        bonded_pairs=[]
        bond_df = pmb.get_instances_df(pmb_type="bond")

        for bond_index in bond_df.index:
            particle_id1= bond_df.loc[bond_index,"particle_id1"]
            particle_id2= bond_df.loc[bond_index,"particle_id2"]
            bonded_pair=frozenset([particle_id1,particle_id2])
            bonded_pairs.append(bonded_pair)
            bonded_in_espresso = False
            for pid in bonded_pair:
                for bond in espresso_system.part.by_id(pid).bonds[:]:
                    partner_id  = bond[1]
                    if partner_id in bonded_pair:
                        bonded_in_espresso=True
                        
            self.assertEqual(first=bonded_in_espresso, 
                            second=True)
        self.assertEqual(first=frozenset(bonded_pairs), 
                        second=frozenset([frozenset([0,1]),
                                            frozenset([0,2]),
                                            frozenset([3,4]),
                                            frozenset([4,5]),
                                            frozenset([4,6])]))
        starting_number_of_particles=len(espresso_system.part.all())
        with self.assertRaises(ValueError):
            pmb.create_residue(name="R51",
                               box_l=box_l,
                               use_default_bond=True)
        # If no particles have been created, the number of particles should be the same as before
        self.assertEqual(first=len(espresso_system.part.all()), 
                        second=starting_number_of_particles)

        # Tests for delete_residue
        # This should delete 3 particles (residue 0 is a R2 residue)
        starting_number_of_particles=len(espresso_system.part.all())
        pmb.delete_instances_in_system(instance_id=0,
                                    pmb_type="residue")
        self.assertEqual(first=len(espresso_system.part.all()), 
                         second=starting_number_of_particles-3)
        # There should be only one residue instance now in the pyMBE database
        self.assertEqual(first=len(pmb.get_instances_df(pmb_type="residue")), 
                                second=1)
        # And there should be only 4 particles (central bead + 2 side chains + central bead of R3)
        self.assertEqual(first=len(pmb.get_instances_df(pmb_type="particle")), 
                                second=4)
        # Delete the other residue instance to simplify the rest of the tests
        pmb.delete_instances_in_system(instance_id=1,
                                    pmb_type="residue")
        
        backbone_vector = np.array([1,3,-4])
        magnitude = np.linalg.norm(backbone_vector)
        backbone_vector = backbone_vector/magnitude
        pmb.create_molecule(name="M2",
                            number_of_molecules=1,
                            box_l=box_l,
                            backbone_vector = backbone_vector,
                            use_default_bond=True)
        pmb.add_instances_to_engine()
        particle_ids=pmb.get_particle_id_map(object_name="M2")["all"]

        # Residue and molecule IDs expected
        # For the  M2 molecule created, the residue and molecule IDs should be as follows:
        # R1 (residue_id=0, molecule_id=0), R2 (residue_id=1, molecule_id=0), R3 (residue_id=2, molecule_id=0)

        residue_ids={0: 0, 1: 1, 2: 1, 3: 1, 4: 2, 5: 2, 6: 2, 7: 2}
        molecule_ids={0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0}

        # Check that the particle properties are correct
        for pid in particle_ids:
            particle=espresso_system.part.by_id(pid)
            particle_tpl = pmb.db.get_instance(pmb_type="particle",
                                                instance_id=pid)
            particle_name = particle_tpl.name
            self.assertEqual(first=particle.type, 
                            second=type_map[particle_name])
            self.assertEqual(first=particle.q, 
                            second=particle_parameters[particle_name]["z"])
            # Check that particles have the correct residue id
            residue_id = particle_tpl.residue_id
            self.assertEqual(first=residue_id, 
                            second=residue_ids[pid])
            # Check that particles have the correct molecule id
            molecule_id = particle_tpl.molecule_id
            self.assertEqual(first=molecule_id, 
                            second=molecule_ids[pid])

        # Check that the molecule has the right residues
        residue_list=[]
        residue_df = pmb.get_instances_df(pmb_type="residue")
        for res_index in residue_df[residue_df['molecule_id']==0].index:
            resname = residue_df.loc[res_index,"name"]
            residue_list.append(resname)
        self.assertEqual(first=frozenset(residue_list), 
                         second=frozenset(molecule_parameters["M2"]["residue_list"]))

        # Expected bonded pairs for the molecule
        # Molecule 0:
        # S1(0)-S1(1) (R1-R2)
        # S1(1)-S2(2) (R2)
        # S1(1)-S3(3) (R2)
        # S2(1)-S2(4) (R2-R3)
        # S2(4)-S1(5) (R3)
        # S1(5)-S2(6) (R3)
        # S1(5)-S3(7) (R3)

        bonded_pairs_ref=[frozenset([0,1]),
                        frozenset([1,2]),
                        frozenset([1,3]),
                        frozenset([1,4]),
                        frozenset([4,5]),
                        frozenset([5,6]),
                        frozenset([5,7])]

        # Check that particles are correctly bonded
        bonded_pairs=[]
        bond_df = pmb.get_instances_df(pmb_type="bond")
        for bond_index in bond_df.index:
            particle_id1= bond_df.loc[bond_index,"particle_id1"]
            particle_id2= bond_df.loc[bond_index,"particle_id2"]
            bonded_pair=frozenset([particle_id1,particle_id2])
            bonded_pairs.append(bonded_pair)
            bonded_in_espresso = False
            for pid in bonded_pair:
                for bond in espresso_system.part.by_id(pid).bonds[:]:
                    partner_id  = bond[1]
                    if partner_id in bonded_pair:
                        bonded_in_espresso=True
            self.assertEqual(first=bonded_in_espresso, 
                            second=True)
        self.assertEqual(first  = frozenset(bonded_pairs), 
                        second = frozenset(bonded_pairs_ref))
        central_bead_positions = []
        residue_map=pmb.get_particle_id_map(object_name="M2")["residue_map"]
        for res_id in residue_map.keys():
            central_bead_id = min(residue_map[res_id])
            central_bead_pos = espresso_system.part.by_id(central_bead_id).pos
            central_bead_positions.append(central_bead_pos)

        # Here one expects 3 central bead positions for residues R1, R2, and R3
        self.assertEqual(len(central_bead_positions),len(molecule_parameters["M2"]["residue_list"]))
        backbone_direction_1 = central_bead_positions[1] - central_bead_positions[0]
        backbone_direction_2 = central_bead_positions[2] - central_bead_positions[1]
        backbone_direction_1 /= np.linalg.norm(backbone_direction_1)
        backbone_direction_2 /= np.linalg.norm(backbone_direction_2)
        np.testing.assert_almost_equal(actual=  list(backbone_direction_1),
                                       desired= list(backbone_vector))
        np.testing.assert_almost_equal(actual= list(backbone_direction_2),
                                       desired= list(backbone_vector))
        starting_number_of_particles=len(espresso_system.part.all())
        pmb.create_molecule(name="M2",
                            number_of_molecules=0,
                            box_l=box_l,
                            use_default_bond=True)
        pmb.create_molecule(name="M2",
                            number_of_molecules=-1,
                            box_l=box_l,
                            use_default_bond=True)
        ### No particle instances are added thus no need to add them to the engine 
        # If no particles have been created, only two particles should be in the system (from the previous test)
        self.assertEqual(first=len(espresso_system.part.all()), 
                        second=starting_number_of_particles)
        # Check that providing the wrong molecule name raises a ValueError
        self.assertRaises(ValueError, pmb.create_molecule,
                        name="M3",
                        number_of_molecules=1,
                        box_l=box_l,
                        use_default_bond=True)
        
        # Tests for delete_molecule
        # create another molecule just to have two molecules in the system
        pmb.create_molecule(name="M2",
                            number_of_molecules=1,
                            box_l=box_l,
                            backbone_vector = backbone_vector,
                            use_default_bond=True)
        # This should delete 8 particles (molecule 0 is a M2 molecule)
        pmb.add_instances_to_engine()
        starting_number_of_particles=len(espresso_system.part.all())
        pmb.delete_instances_in_system(instance_id=0,
                                    pmb_type="molecule")
        self.assertEqual(first=len(espresso_system.part.all()), 
                        second=starting_number_of_particles-8)
        # There should only one molecule instance now in the pyMBE database
        self.assertEqual(first=len(pmb.get_instances_df(pmb_type="molecule")), 
                        second=1)
        # There should be only 3 residues (from the remaining M2 molecule)
        self.assertEqual(first=len(pmb.get_instances_df(pmb_type="residue")), 
                        second=3)
        # There should be only 8 particles (from the remaining M2 molecule)
        self.assertEqual(first=len(pmb.get_instances_df(pmb_type="particle")), 
                        second=8)
        
    def test_set_particle_initial_state(self):
        """
        Unit tests for set_particle_initial_state
        """
        pmb = pyMBE.pymbe_library(23)
        pmb.define_particle(name="A",
                            sigma=1*pmb.units.nm,
                            epsilon=1*pmb.units.reduced_energy,
                            z=1)
        state_1 = pmb.db.get_template(pmb_type="particle", name="A").initial_state
        pmb.set_particle_initial_state(particle_name="A",
                                       state_name="random")
        state_2 = pmb.db.get_template(pmb_type="particle", name="A").initial_state
        self.assertIsNot(state_1,
                         state_2)
        self.assertEqual(state_2,
                         "random")

    def test_get_radius_map(self):
        """
        Tests for get_radius_map        
        """

        self.assertEqual(first=len(pmb.get_radius_map()),
                         second=len(particle_parameters.values()))

        
        second_radii=[]
        for particle in particle_parameters.values():
            second_radii.append((particle['sigma'].m_as("reduced_length")+particle['offset'].m_as("reduced_length"))/2)

        first_radii=[pmb.get_radius_map()[0],
                    pmb.get_radius_map()[1],
                    pmb.get_radius_map()[2],]

        for index in range(len(first_radii)):
            self.assertAlmostEqual(first=first_radii[index],
                                   second=second_radii[index])

        self.assertEqual(first=isinstance(pmb.get_radius_map()[0],float),
                         second=True)
        self.assertEqual(first=pmb.get_radius_map(dimensionless=False)[0].dimensionality,
                         second=pmb.units.nm.dimensionality)  

        # Test the sanity test
        pmb2 = pyMBE.pymbe_library(24)
        empty_map = pmb2.get_radius_map()
        self.assertEqual(empty_map,
                         {})  

if __name__ == "__main__":
    ut.main()
