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
import pyMBE
import espressomd
from pyMBE.storage.templates.hydrogel import HydrogelTemplate
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.instances.residue import ResidueInstance
from pyMBE.storage.instances.molecule import MoleculeInstance
from pyMBE.storage.instances.peptide import PeptideInstance
from pyMBE.storage.instances.protein import ProteinInstance
from pyMBE.storage.instances.bond import BondInstance
from pyMBE.storage.instances.angle import AngleInstance
from pyMBE.storage.instances.hydrogel import HydrogelInstance
from pyMBE.storage.templates.bond import BondTemplate
from pyMBE.storage.templates.angle import AngleTemplate
from pyMBE.storage.templates.hydrogel import HydrogelNode
from pyMBE.storage.pint_quantity import PintQuantity
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
from pydantic.error_wrappers import ValidationError
import pint
import numpy as np
box_l=[10]*3
espresso_system=espressomd.System(box_l =box_l )

class Test(ut.TestCase):
    def test_label_id_map(self):
        pmb = pyMBE.pymbe_library(seed=23)

        self.assertEqual(pmb.db._get_label_id_map("particle"), "particle_map")
        self.assertEqual(pmb.db._get_label_id_map("molecule"), "molecule_map")
        self.assertEqual(pmb.db._get_label_id_map("hydrogel"), "assembly_map")

    def test_particle_instance_position_is_nparray(self):
        """Checks that is not possible to instantiate a particle instance with a list it is only possible with a numpy array.
        """
        
        with self.assertRaises(ValidationError):
            ParticleInstance(name="B",
                                    particle_id=1,
                                    initial_state="B",
                                    position=[box_l[0]*0.5,box_l[1]*0.5,box_l[2]*0.5])
        

    def test_sanity_hydrogel_node_template(self):
        """
        Sanity test for HydrogelNode template validator
        """
        inputs={"particle_name": "A",
                "lattice_index": 1}
        self.assertRaises(ValueError,
                          HydrogelNode,
                          **inputs)

    def test_sanity_database_methods(self):
        """
        Sanity tests for exceptions in:
             _register_instance
             _update_instance
             get_instances
             delete_instance
             delete_instances
             delete_templates
             _register_reaction
             get_reaction
             delete_reaction
        """
        pmb = pyMBE.pymbe_library(23)
        # Unit tests for _register_instance
        class DummyInstance():
            pass

        inputs = {"instance": DummyInstance()}
        self.assertRaises(TypeError,
                          pmb.db._register_instance,
                          **inputs)
        pmb.define_particle(name="A",
                            sigma=1*pmb.units.nm,
                            epsilon=1*pmb.units.reduced_energy,
                            pka=9,
                            acidity="acidic")
        part_inst = ParticleInstance(name="A",
                                    particle_id=0,
                                    initial_state="A",
                                    position=np.array([box_l[0]*0.5,box_l[1]*0.5,box_l[2]*0.5]))
        pmb.db._register_instance(part_inst)
        inputs = {"instance": part_inst}
        self.assertRaises(ValueError,
                          pmb.db._register_instance,
                          **inputs)
        templateless_part_inst = ParticleInstance(name="B",
                                    particle_id=1,
                                    initial_state="B",
                                    position=np.array([box_l[0]*0.5,box_l[1]*0.5,box_l[2]*0.5]))
        inputs = {"instance": templateless_part_inst}
        self.assertRaises(ValueError,
                          pmb.db._register_instance,
                          **inputs)
        # Unit test for get_instances
        self.assertEqual(pmb.db._instances["particle"],
                         pmb.db.get_instances(pmb_type="particle"))

        # Unit tests for _update_instance
        inputs = {"instance_id": 2,
                  "pmb_type": "particle",
                  "attribute": "particle_id",
                  "value": 0}
        self.assertRaises(ValueError,
                          pmb.db._update_instance,
                          **inputs)
        
        pmb.db._register_template(HydrogelTemplate(name="test",
                                                    node_map=[],
                                                    chain_map=[]))
        
        pmb.db._register_instance(HydrogelInstance(name="test",
                                                  assembly_id=0))
        inputs = {"instance_id": 0,
                  "pmb_type": "hydrogel",
                  "attribute": "assembly_id",
                  "value": 1}
        self.assertRaises(ValueError,
                          pmb.db._update_instance,
                          **inputs)
        
        # Unit test for _register_reaction
        inputs = {"participants":[ReactionParticipant(particle_name="A",
                                                   state_name="A",
                                                   coefficient=-1),
                                  ReactionParticipant(particle_name="B",
                                                   state_name="B",
                                                   coefficient=1)],
                 "pK":1,
                 "reaction_type":"test"}
        reaction = Reaction(**inputs)
        inputs = {"reaction": reaction}
        pmb.db._register_reaction(reaction)
        self.assertRaises(ValueError,
                          pmb.db._register_reaction,
                          **inputs)
        # Unit tests for get_reaction:
        ## Test that one gets back the right reaction
        self.assertEqual(reaction,
                         pmb.db.get_reaction(name=reaction.name))
        ## Sanity test, giving an unknown reaction name triggers a ValueError
        inputs = {"name" : "test"}
        self.assertRaises(ValueError,
                          pmb.db.get_reaction,
                          **inputs)
        # Sanity test for delete_reaction
        inputs = {"reaction_name": "test"}
        self.assertRaises(ValueError,
                          pmb.db.delete_reaction,
                          **inputs)
        
        # Sanity Unit test for delete_instance
        inputs = {"pmb_type": "molecule",
                  "instance_id": 0}
        self.assertRaises(ValueError,
                          pmb.db.delete_instance,
                          **inputs)
        inputs = {"pmb_type": "particle",
                  "instance_id": 3}
        self.assertRaises(ValueError,
                          pmb.db.delete_instance,
                          **inputs)
        # Sanity tests for delete_template
        ## Triggers a ValueError because no molecule template has been defined
        inputs = {"pmb_type": "molecule",
                  "name": "test"}
        self.assertRaises(ValueError,
                          pmb.db.delete_template,
                          **inputs)
        ## Triggers a ValueError because no particle of this name has been defined
        inputs = {"pmb_type": "particle",
                  "name": "test"}
        self.assertRaises(ValueError,
                          pmb.db.delete_template,
                          **inputs)
        ## Triggers a ValueError because particle instances of this template have been created
        inputs = {"pmb_type": "particle",
                  "name": "A"}
        self.assertRaises(ValueError,
                          pmb.db.delete_template,
                          **inputs)
        # Unit tests for delete_instances
        ## Trying to delete instances from an empty category does nothing
        previous_instances = pmb.db._instances.copy()
        pmb.db.delete_instances(pmb_type="molecule")
        self.assertEqual(previous_instances,
                         pmb.db._instances)
        
        ## Calling the function deletes all instances of a given pmb_type
        part_inst = ParticleInstance(name="A",
                                    particle_id=1,
                                    initial_state="A",
                                    position=np.array([box_l[0]*0.5,box_l[1]*0.5,box_l[2]*0.5]))
        pmb.db._register_instance(part_inst)
        pmb.db.delete_instances(pmb_type="particle")
        assert "particle" not in pmb.db._instances.keys()
    
    def test_find_instance_ids(self):
        """
        Sanity test for `_find_instance_ids_by_attribute`
        and `_find_instance_ids_by_name`
        """
        pmb = pyMBE.pymbe_library(23)
        pmb.define_particle(name="A",
                            sigma=1*pmb.units.nm,
                            epsilon=1*pmb.units.reduced_energy,
                            pka=9,
                            acidity="acidic")
        pmb.define_particle(name="B",
                            sigma=1*pmb.units.nm,
                            epsilon=1*pmb.units.reduced_energy)
        pmb.define_residue(name="R1",
                           central_bead="A",
                           side_chains=["B"])
        bond_type = 'harmonic'
        bond = {'r_0'    : 0.4*pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        pmb.define_default_bond(bond_type = bond_type,
                                bond_parameters = bond)
        pmb.define_molecule(name="M1",
                            residue_list=["R1"]*2)
        pmb.create_molecule(name="M1",
                            number_of_molecules=1,
                            use_default_bond=True,
                            box_l=box_l)
        instance_ids_r1 = pmb.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                              attribute="residue_id",
                                                              value=0)
        instance_ids_r2 = pmb.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                              attribute="residue_id",
                                                              value=1)
        self.assertEqual(instance_ids_r1,
                         [0,1])
        self.assertEqual(instance_ids_r2,
                         [2,3])
        instance_ids_m1 = pmb.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                              attribute="molecule_id",
                                                              value=0)
        self.assertEqual(instance_ids_m1,
                         [0,1,2,3])
        instance_ids_by_name_A = pmb.db.find_instance_ids_by_name(pmb_type="particle",
                                                                   name="A")
        instance_ids_by_name_B = pmb.db.find_instance_ids_by_name(pmb_type="particle",
                                                                   name="B")
        self.assertEqual(instance_ids_by_name_A,
                         [0,2])
        self.assertEqual(instance_ids_by_name_B,
                         [1,3])
        # Sanity test, no ids are returned if the instance does not exist   
        instance_ids_test = pmb.db.find_instance_ids_by_name(pmb_type="peptide",
                                                                   name="B")
        self.assertEqual(instance_ids_test,
                         [])
          
        # Check that the pyMBE database finds a specific instance
        self.assertEqual(pmb.db._has_instance(pmb_type="particle",
                                              instance_id=3),
                         True)
        self.assertEqual(pmb.db._has_instance(pmb_type="particle",
                                              instance_id=4),
                         False)
        # Sanity test, unexisting pyMBE type
        inputs = {"pmb_type": "unknown",
                  "instance_id": 0}
        self.assertRaises(ValueError,
                          pmb.db._has_instance,
                          **inputs)
        
    def test_count_templates(self):
        """
        Sanity test for `_collect_particle_templates`
        """
        pmb = pyMBE.pymbe_library(23)
        pmb.define_particle(name="A",
                            sigma=1*pmb.units.nm,
                            epsilon=1*pmb.units.reduced_energy,
                            pka=9,
                            acidity="acidic")
        pmb.define_particle(name="B",
                            sigma=1*pmb.units.nm,
                            epsilon=1*pmb.units.reduced_energy)
        pmb.define_residue(name="R1",
                           central_bead="A",
                           side_chains=["B"])
        pmb.define_residue(name="R2",
                           central_bead="A",
                           side_chains=["R1"])
        pmb.define_molecule(name="M1",
                            residue_list=["R1"]*2)
        A_states = pmb.db._collect_particle_templates(name="A",
                                                      pmb_type="particle_state")
        self.assertEqual(A_states,
                         {"A":1})
        AH_states = pmb.db._collect_particle_templates(name="AH",
                                                      pmb_type="particle_state")
        self.assertEqual(AH_states,
                         {"A":1})
        A_particles = pmb.db._collect_particle_templates(name="A",
                                                      pmb_type="particle")
        B_particles = pmb.db._collect_particle_templates(name="B",
                                                      pmb_type="particle")
        self.assertEqual(A_particles,
                         {"A":1})
        self.assertEqual(B_particles,
                         {"B":1})
        R1_counts = pmb.db._collect_particle_templates(name="R1",
                                                      pmb_type="residue")
        self.assertEqual(R1_counts,
                         {"A":1,
                          "B":1})
        R2_counts = pmb.db._collect_particle_templates(name="R2",
                                                      pmb_type="residue")
        self.assertEqual(R2_counts,
                         {"A":2,
                          "B":1})
        M1_counts = pmb.db._collect_particle_templates(name="M1",
                                                      pmb_type="molecule")
        self.assertEqual(M1_counts,
                         {"A":2,
                          "B":2})
        inputs={"name": "test",
                "pmb_type": "unknown"}
        self.assertRaises(NotImplementedError,
                          pmb.db._collect_particle_templates,
                          **inputs)
        # Sanity test for unknown types in _has_template
        inputs = {"pmb_type": "unknown",
                  "name": "A"}
        self.assertRaises(ValueError,
                          pmb.db._has_template,
                          **inputs)
        # Sanity tests for get_particle_templates_under
        templates_R1 = pmb.db.get_particle_templates_under(template_name="R1")
        self.assertEqual(templates_R1,
                         {"A","B"})
        # Sanity tests, raise ValueError when pmb_type cannot be safely infered
        pmb.define_residue(name="A",
                           central_bead="A",
                           side_chains=["B"])
        inputs = {"template_name": "A"}
        self.assertRaises(ValueError,
                          pmb.db.get_particle_templates_under,
                          **inputs)

    def test_sanity_db(self):
        """
        Sanity tests for the pyMBE database
        """
        pmb = pyMBE.pymbe_library(23)
        pmb.define_molecule(name ="test",
                            residue_list=[])
        pmb.define_peptide(name="test",
                           sequence="",
                           model="1beadAA")
        inputs = {"name": "test",
                  "allowed_types": {"molecule", "peptide"}}
        self.assertRaises(ValueError,
                          pmb._get_template_type,
                          **inputs)
        
    def test_instance_id_validators(self):
        """
        Tests that negative values of instances raise a ValueError in the pyMBE database
        """
        inputs = {"name":"A",
                  "particle_id":-1,
                   "initial_state":"A",
                   "position":np.array([box_l[0]*0.5,box_l[1]*0.5,box_l[2]*0.5])}
        self.assertRaises(ValueError,
                          ParticleInstance,
                          **inputs)
        inputs = {"name":"A",
                  "residue_id":-1}
        self.assertRaises(ValueError,
                          ResidueInstance,
                          **inputs)
        inputs = {"name":"A",
                  "molecule_id":-1}
        self.assertRaises(ValueError,
                          MoleculeInstance,
                          **inputs)
        inputs = {"name":"A",
                  "molecule_id":-1}
        self.assertRaises(ValueError,
                          PeptideInstance,
                          **inputs)
        inputs = {"name":"A",
                  "molecule_id":-1}
        self.assertRaises(ValueError,
                          ProteinInstance,
                          **inputs)
        inputs = {"name":"A",
                  "assembly_id":-1}
        self.assertRaises(ValueError,
                          HydrogelInstance,
                          **inputs)
        inputs = {"name":"A",
                  "bond_id":-1,
                  "particle_id1":1,
                  "particle_id2":2}
        self.assertRaises(ValueError,
                          BondInstance,
                          **inputs)
        inputs = {"name":"A",
                  "bond_id":1,
                  "particle_id1":-1,
                  "particle_id2":2}
        self.assertRaises(ValueError,
                          BondInstance,
                          **inputs)
        inputs = {"name":"A",
                  "bond_id":1,
                  "particle_id1":1,
                  "particle_id2":-2}
        self.assertRaises(ValueError,
                          BondInstance,
                          **inputs)
        inputs = {"name":"A",
                  "angle_id":-1,
                  "particle_id1":1,
                  "particle_id2":2,
                  "particle_id3":3}
        self.assertRaises(ValueError,
                          AngleInstance,
                          **inputs)
        
    def test_make_name_bond_template(self):
        inputs = {"bond_type": "harmonic",
                  "parameters": {"r": PintQuantity(magnitude=1,
                                                   units="nm",
                                                   dimension="length")}}
        bond_tpl = BondTemplate(**inputs)
        self.assertRaises(RuntimeError,
                          bond_tpl._make_name)

    def test_make_name_angle_template(self):
        inputs = {"angle_type": "harmonic",
                  "parameters": {"k": PintQuantity(magnitude=1,
                                                   units="reduced_energy",
                                                   dimension="energy"),
                                 "phi_0": PintQuantity(magnitude=1,
                                                        units="dimensionless",
                                                        dimension="dimensionless")}}
        angle_tpl = AngleTemplate(**inputs)
        self.assertRaises(RuntimeError,
                          angle_tpl._make_name)
        
    def test_exceptions_pint_quantity(self):
        units = pint.UnitRegistry()
        inputs = {"q":1,
                  "expected_dimension": "length",
                  "ureg": units}
        self.assertRaises(TypeError,
                          PintQuantity.from_quantity,
                          **inputs)
        inputs = {"q":1*units.nm,
                  "expected_dimension": "unknown",
                  "ureg": units}
        self.assertRaises(ValueError,
                          PintQuantity.from_quantity,
                          **inputs)
        inputs = {"q":1*units.nm**2,
                  "expected_dimension": "length",
                  "ureg": units}
        self.assertRaises(ValueError,
                          PintQuantity.from_quantity,
                          **inputs)
        
    def test_exceptions_reaction_template(self):
        """
        Tests sanity of the Reaction template
        """
        # Reactions with less than 2 participants trigger a value error
        inputs = {"participants":[ReactionParticipant(particle_name="A",
                                                   state_name="A",
                                                   coefficient=1)],
                 "pK":1,
                 "reaction_type":"test"}
        
        self.assertRaises(ValueError,
                          Reaction,
                          **inputs)
        # Reactions with a participant with a 0 stechiometric coeff. trigger a value error
        inputs = {"particle_name":"A",
                   "state_name":"A",
                   "coefficient":0} 
        self.assertRaises(ValueError,
                          ReactionParticipant,
                          **inputs)
        
        # Adding a new participant with a 0 stechiometric coeff. triggers a ValueError
        react_tpl =Reaction(participants=[ReactionParticipant(particle_name="A",
                                                   state_name="A",
                                                   coefficient=-1),
                                            ReactionParticipant(particle_name="B",
                                                            state_name="B",
                                                            coefficient=1)],
                            pK=1,
                            reaction_type="test")
        inputs={"particle_name": "C",
                "state_name":"C",
                "coefficient":0}
        self.assertRaises(ValueError,
                          react_tpl.add_participant,
                          **inputs)

if __name__ == '__main__':
    ut.main()
