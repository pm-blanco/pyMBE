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

import numpy as np 
import espressomd
import unittest as ut
import tempfile
import os
import json
import pathlib
import pyMBE
import pyMBE.lib.handy_functions as hf

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)
Box_L = 100 * pmb.units.reduced_length
box_l=[Box_L.to('reduced_length').magnitude] * 3
espresso_system=espressomd.System(box_l = box_l)
pmb.set_simulation_engine(espresso_system)
protein_pdb = '1f6s'
path_to_parfile = pathlib.Path(__file__).parent / "tests_data" / "protein_topology_dict.json"
path_to_cg=pmb.root / "parameters" / "globular_proteins" / f"{protein_pdb}.vtf"
topology_dict, sequence = pmb.read_protein_vtf (filename=path_to_cg)
ref_sequence = "nQLTKCEVFRELKDLKGYGGVSLPEWVCTTFHTSGYDTQAIVQNNDSTEYGLFQINNKIWCKDDQNPHSSNICNISCDKFLDDDLTDDIMCVKKILDKVGINYWLAHKALCSEKLDQWLCEKc" # VTF file is missing a terminal L aminoacid on its sequence
ref_residue_list = []
for key in ref_sequence:
    ref_residue_list.append(f"AA-{key}")


class Test(ut.TestCase):
    def test_protein_setup(self):
        """
        Unit tests for setting up globular proteins in pyMBE.
        """
        def custom_deserializer(dct):
            if "value" in dct and "unit" in dct:
                return pmb.units.Quantity(dct["value"], dct["unit"])  
            return dct  
        with open (path_to_parfile, "r") as file:
            load_json = json.load(file,object_hook=custom_deserializer)
        np.testing.assert_equal(actual= topology_dict, 
                                desired= load_json,
                                verbose = True)
        protein_model = '2beadAA'
        hf.define_protein_AA_particles(topology_dict=topology_dict,
                                       pmb=pmb,
                                       pka_set={})
        residue_list = hf.define_protein_AA_residues(sequence=sequence,
                                                    model=protein_model,
                                                    pmb=pmb)
        # Sanity test for unvalid lj_setups
        input_params = {"topology_dict": topology_dict,
                        "pmb": pmb,
                        "pka_set": {},
                        "lj_setup_mode": "random"}
        self.assertRaises(ValueError,
                          hf.define_protein_AA_particles,
                          **input_params,)


        # Define a residue for the metal ion
        pmb.define_residue(name="AA-Ca",
                           central_bead="Ca",
                           side_chains=[])
        pmb.define_protein (name=protein_pdb,
                            sequence=sequence, 
                            model = protein_model)
        # Check residue templates
        for residue_name in residue_list:
            residue_template = pmb.db.get_template(name=residue_name,
                                                pmb_type="residue")
            assert residue_template is not None
            assert residue_template.pmb_type == "residue"
            assert residue_template.name == residue_name
        # Check protein template
        protein_template = pmb.db.get_template(name=protein_pdb,
                                            pmb_type="protein")
        assert protein_template is not None
        assert protein_template.name == protein_pdb
        np.testing.assert_equal(actual=protein_template.sequence,
                                desired=ref_sequence,
                                verbose=True)
        np.testing.assert_equal(actual=protein_template.residue_list,
                                desired=ref_residue_list,
                                verbose=True)
        input_parameters={"name": protein_pdb,
                        "sequence": sequence,
                        "model" : "3beadAA"}
        np.testing.assert_raises(ValueError, 
                                 pmb.define_protein, 
                                 **input_parameters)
        with self.assertRaises(ValueError):
            pmb.create_protein(name="missing_protein_template",
                               number_of_proteins=1,
                               box_l=box_l,
                               topology_dict=topology_dict)
        molecule_id = pmb.create_protein(name=protein_pdb,
                                        number_of_proteins=1,
                                        box_l=box_l,
                                        topology_dict=topology_dict)[0]
        pmb.add_instances_to_engine()
        particle_id_list = pmb.get_particle_id_map(object_name=protein_pdb)["all"]
        center_of_mass_es = pmb.calculate_center_of_mass(instance_id=molecule_id,
                                                        pmb_type="protein")
        center_of_mass = np.zeros(3)
        axis_list = [0,1,2]
        for aminoacid in topology_dict.keys():
            initial_pos = topology_dict[aminoacid]['initial_pos']
            for axis in axis_list:
                center_of_mass[axis] +=  initial_pos[axis]
        center_of_mass = center_of_mass/ len(topology_dict.keys())
        distance_es = np.zeros(3)
        distance_topology = np.zeros(3)
        for id in particle_id_list:
            initial_pos_es = espresso_system.part.by_id(id).pos
            part_inst = pmb.db.get_instance(instance_id=id,
                                       pmb_type="particle")
            part_tpl = pmb.db.get_template(name=part_inst.name,
                                           pmb_type="particle")
            part_state =  pmb.db.get_template(name=part_tpl.initial_state,
                                           pmb_type="particle_state")
            charge = espresso_system.part.by_id(id).q
            es_type = espresso_system.part.by_id(id).type
            np.testing.assert_equal(part_state.z,
                                    charge)
            np.testing.assert_equal(part_state.es_type,
                                    es_type)
            residue_id = part_inst.residue_id
            res_inst = pmb.db.get_instance(instance_id=residue_id,
                                            pmb_type="residue")
            residue_name = res_inst.name
            if "G" in residue_name:
                continue
            initial_pos = topology_dict[f"{residue_name[3:]}{residue_id}"]['initial_pos']
            for axis in axis_list:
                distance_es[axis] = (initial_pos_es[axis] - center_of_mass_es[axis])**2
                distance_topology[axis] = (initial_pos[axis] - center_of_mass[axis])**2
            relative_distance_es = np.sqrt(np.sum(distance_es))
            relative_distance = np.sqrt(np.sum(distance_es))
            np.testing.assert_equal(actual=relative_distance_es, 
                                    desired=relative_distance, 
                                    verbose=True)
        starting_number_of_particles=len(espresso_system.part.all())
        pmb.create_protein(name=protein_pdb,
                            number_of_proteins=0,
                            box_l=box_l,
                            topology_dict=topology_dict)
        pmb.create_protein(name=protein_pdb,
                            number_of_proteins=-1,
                            box_l=box_l,
                            topology_dict=topology_dict)
        ### No particles instances are created, thus it is not needed to call add_instances_to_engine
        np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                                desired=starting_number_of_particles, 
                                verbose=True)
        positions = []
        for pid in particle_id_list:
            positions.append(espresso_system.part.by_id(pid).pos)
        pmb.enable_motion_of_rigid_object(instance_id=molecule_id,
                                        pmb_type="protein")

        momI = 0
        center_of_mass = pmb.calculate_center_of_mass(instance_id=molecule_id, 
                                                      pmb_type="protein")
        for p in espresso_system.part:
            if p.mass > 1: 
                rigid_object_id = p.id 
                rigid_object_mass = espresso_system.part.by_id(rigid_object_id).mass
                rigid_object_rotation = espresso_system.part.by_id(rigid_object_id).rotation
                rigid_object_intertia  = np.copy(espresso_system.part.by_id(rigid_object_id).rinertia)
                np.testing.assert_equal(actual=rigid_object_mass, 
                                        desired=len(particle_id_list), 
                                        verbose=True)
                np.testing.assert_equal(actual=rigid_object_rotation, 
                                        desired=[1, 1, 1], 
                                        verbose=True)
                for pid in particle_id_list:
                    momI += np.power(np.linalg.norm(center_of_mass - espresso_system.part.by_id(pid).pos), 2)
                rinertia = np.ones(3) * momI
                np.testing.assert_array_almost_equal(rinertia, rigid_object_intertia)

    def test_define_peptide_1beadAA(self):
        """
        Test that define_peptide_AA_residues correctly defines
        residue templates for a short peptide using the 1beadAA model.
        """
        pmb2 = pyMBE.pymbe_library(seed=123)
        # Define particles needed by the residues
        for bead in ["A", "G", "L"]:
            pmb2.define_particle(name=bead,
                                 sigma=0.355 * pmb2.units.nm,
                                 epsilon=1 * pmb2.units("reduced_energy"))
        sequence = ["A", "G", "L", "A"]  # include repetition on purpose
        model = "1beadAA"
        hf.define_peptide_AA_residues(sequence=sequence, 
                                      model=model, 
                                      pmb=pmb2)
        # Expected residue template names
        expected_residues = {"AA-A", "AA-G", "AA-L"}
        # Retrieve all residue templates from the DB
        stored_residues = {tpl.name for tpl in pmb2.db.get_templates(pmb_type="residue").values()}
        self.assertSetEqual(stored_residues,
                            expected_residues,
                            "Residue templates stored in DB do not match expected 1beadAA residues")
        # Check residue definitions
        for res_name in expected_residues:
            residue = pmb2.db.get_template(pmb_type="residue",
                                           name=res_name)
            aa = res_name.split("-")[1]
            self.assertEqual(residue.central_bead,
                             aa,
                             f"Central bead for {res_name} should be '{aa}'")
            self.assertEqual(residue.side_chains,
                             [],
                             f"Residue {res_name} should have no side chains in 1beadAA model")

    def test_define_peptide_residues_2beadAA(self):
        """
        Test residue definition for the 2beadAA model:
        - standard residues use CA + side-chain
        - G, c, n are single-bead residues
        - residues are defined only once
        """
        pmb2 = pyMBE.pymbe_library(seed=123)
        sequence = ["A", "G", "L", "c", "n", "A"]
        hf.define_peptide_AA_residues(sequence, 
                                      model="2beadAA", 
                                      pmb=pmb2)
        # Expected residue templates
        expected = {"AA-A":  {"central_bead": "CA", "side_chains": ["A"]},
                    "AA-L":  {"central_bead": "CA", "side_chains": ["L"]},
                    "AA-G":  {"central_bead": "G",  "side_chains": []},
                    "AA-c":  {"central_bead": "c",  "side_chains": []},
                    "AA-n":  {"central_bead": "n",  "side_chains": []},}
        for resname, props in expected.items():
            tpl = pmb2.db.get_template(pmb_type="residue", name=resname)
            self.assertEqual(tpl.central_bead, props["central_bead"])
            self.assertEqual(tpl.side_chains, props["side_chains"])
        # Ensure residues were defined only once
        residue_templates = pmb2.db.get_templates(pmb_type="residue").values()
        residue_names = [tpl.name for tpl in residue_templates]
        self.assertEqual(len(residue_names), len(set(residue_names)))
        self.assertEqual(set(residue_names), set(expected.keys()))

    def test_protein_parser(self):
        """
        Unit tests for protein_sequence_parser
        """
        def test_sequence(input,output):
            """
            Tests that the pyMBE parses correctly the input sequence.

            Args:
                input(`str` or `lst` of `str`): input protein sequence.
                ouput(`lst` of `str`): expected ouput protein sequence.
            """
            clean_sequence= hf.protein_sequence_parser(sequence = input)
            np.testing.assert_equal(actual=clean_sequence, 
                                    desired=output, 
                                    verbose=True)
            # check that  correctly returns que protein sequence 
        test_sequence(input="REKH",
                    output=["R", "E", "K", "H"])
        test_sequence(input="rekh",
                    output=["R", "E", "K", "H"])
        test_sequence(input="R-E-K-H",
                    output=["R", "E", "K", "H"])
        test_sequence(input="r-e-k-h",
                    output=["R", "E", "K", "H"])
        test_sequence(input="ARG-GLU-LYS-HIS",
                    output=["R", "E", "K", "H"])
        test_sequence(input="arg-glu-lys-his",
                    output=["R", "E", "K", "H"])
        test_sequence(input=["R","E", "K", "H"],
                    output=["R", "E", "K", "H"])
        test_sequence(input=["r","e", "k", "h"],
                    output=["R", "E", "K", "H"])
        test_sequence(input=["ARG","GLU", "LYS", "HIS"],
                    output=["R", "E", "K", "H"])
        test_sequence(input=["arg","glu", "lys", "his"],
                    output=["R", "E", "K", "H"])

        # check that protein_sequence_parser() raises a ValueError if a wrong residue key is provided
        input_parameters = {"sequence":"rekx"}
        np.testing.assert_raises(ValueError, 
                                 hf.protein_sequence_parser, 
                                 **input_parameters)
        input_parameters = {"sequence":"ARG-GLU-TUR-HIS"}
        np.testing.assert_raises(ValueError, 
                                 hf.protein_sequence_parser, 
                                 **input_parameters)
        input_parameters = {"sequence":"A-E-E-X"}
        np.testing.assert_raises(ValueError, 
                                 hf.protein_sequence_parser, 
                                 **input_parameters)
        input_parameters = {"sequence":"a-e-e-x"}
        np.testing.assert_raises(ValueError, 
                                 hf.protein_sequence_parser,
                                 **input_parameters)
        input_parameters = {"sequence":["A", "E","X"]}
        np.testing.assert_raises(ValueError, 
                                 hf.protein_sequence_parser, 
                                 **input_parameters)
    def test_check_aminoacid_key(self):
        """
        Unit tests for check_aminoacid_key()
        """
        # Check  that check_aminoacid_key returns True for any latter valid in the one letter amino acid code
        valid_AA_keys=['V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C']
        for key in valid_AA_keys:
            np.testing.assert_equal(actual=hf.check_aminoacid_key(key=key), 
                                desired=True, 
                                verbose=True)
        
        # Check  that check_aminoacid_key returns False for a key not valid in the one letter amino acid code 
        np.testing.assert_equal(actual=hf.check_aminoacid_key(key="B"), 
                                desired=False, 
                                verbose=True)
        

    def test_metal_functions(self):
        """
        Unit tests for the helpers for the protein metal ions
        """
        # Check  that check_if_metal_ion returns True for any key corresponding to a supported metal ion 
        for key in hf.get_metal_ions_charge_number_map().keys():
            np.testing.assert_equal(actual=hf.check_if_metal_ion(key=key), 
                                    desired=True, 
                                    verbose=True)
        # Check  that check_if_metal_ion returns False for a key not corresponding to a supported metal ion 
        np.testing.assert_equal(actual=hf.check_if_metal_ion(key="B"), 
                                desired=False, 
                                verbose=True)
        # Check  that get_metal_ions_charge_number_map returns the correct charge map for metals 
        metal_charge_number_map = {"Ca": 2}
        pmb_metal_charge_number_map = hf.get_metal_ions_charge_number_map()

        np.testing.assert_equal(actual=pmb_metal_charge_number_map, 
                                desired=metal_charge_number_map, 
                                verbose=True)

    def test_define_protein_AA_residues(self):
        """
        Unit test for define_protein_AA_residues
        """
        valid_protein_model = ['1beadAA', 
                               '2beadAA']
        test_sequence = ['c','n', 'G','V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C' ]

        valid_protein_model = ['1beadAA', '2beadAA']

        output =['AA-c', 'AA-n', 'AA-G', 'AA-V', 'AA-I', 'AA-L', 'AA-E', 'AA-Q', 'AA-D', 'AA-N', 'AA-H', 'AA-W', 'AA-F', 'AA-Y', 'AA-R', 'AA-K', 'AA-S', 'AA-T', 'AA-M', 'AA-A', 'AA-G', 'AA-P', 'AA-C']
        test_pmb = pyMBE.pymbe_library(23)
        for protein_model in valid_protein_model:
            pmb_residue_list = hf.define_protein_AA_residues(sequence=test_sequence,
                                                             model = protein_model,
                                                             pmb=test_pmb)
            np.testing.assert_equal(actual=pmb_residue_list, 
                                    desired=output, 
                                    verbose=True)    
            test_pmb.db.delete_templates(pmb_type="residue")
        

    def test_define_peptide_sanity(self):
        """
        Sanity tests for define_peptide
        """
        # check that define_peptide() raises a ValueError if a wrong model key is provided
        input_parameters = {"name": "generic_peptide",
                            "sequence": "EEEEEEE",
                            "model": "3beadAA" }
        np.testing.assert_raises(ValueError, 
                                 pmb.define_peptide, 
                                 **input_parameters)
        input_parameters = {"name": "generic_peptide",
                            "sequence": "EEEEEEE",
                            "model": "beadAA" }
        np.testing.assert_raises(ValueError, 
                                 pmb.define_peptide, 
                                 **input_parameters)
    
    def test_read_protein_vtf_unknown_residue(self):
        pmb = pyMBE.pymbe_library(23)
        vtf_content = """\
atom 1 name CA resname XXX resid 1 chain A radius 1.0
1 0.0 0.0 0.0
"""
        with tempfile.NamedTemporaryFile("w", delete=False) as f:
            f.write(vtf_content)
            filename = f.name
        try:
            with self.assertRaisesRegex(ValueError,"Unknown residue name 'XXX' in VTF file"):
                pmb.read_protein_vtf(filename)
        finally:
            os.remove(filename)

    def test_read_protein_vtf_duplicate_particle_label(self):
        pmb = pyMBE.pymbe_library(23)
        vtf_content = """\
atom 1 name CA resname ALA resid 1 chain A radius 1.0
atom 2 name CA resname ALA resid 1 chain A radius 1.0
1 0.0 0.0 0.0
2 1.0 0.0 0.0
"""
        with tempfile.NamedTemporaryFile("w", delete=False) as f:
            f.write(vtf_content)
            filename = f.name
        try:
            with self.assertRaisesRegex(ValueError,"Duplicate particle label 'CA1'"):
                pmb.read_protein_vtf(filename)
        finally:
            os.remove(filename)

class TestGetResiduesFromTopologyDict(ut.TestCase):

    def test_1beadAA_basic(self):
        """
        Single-bead-per-residue model:
        residue name is taken directly from the bead prefix.
        """
        topology = {"A1": {},
                    "G2": {},
                    "L3": {}}
        residues = hf.get_residues_from_topology_dict(topology, model="1beadAA")
        self.assertEqual(len(residues), 3)
        self.assertEqual(residues["1"]["resname"], "A")
        self.assertEqual(residues["2"]["resname"], "G")
        self.assertEqual(residues["3"]["resname"], "L")
        self.assertEqual(residues["1"]["beads"], ["A1"])
        self.assertEqual(residues["2"]["beads"], ["G2"])
        self.assertEqual(residues["3"]["beads"], ["L3"])

    def test_2beadAA_basic(self):
        """
        Two-bead-per-residue model:
        CA beads are ignored when determining residue name.
        """
        topology = {"CA1": {},
                    "L1": {},
                    "CA2": {},
                    "V2": {}}
        residues = hf.get_residues_from_topology_dict(topology, model="2beadAA")
        self.assertEqual(len(residues), 2)
        self.assertEqual(residues["1"]["resname"], "L")
        self.assertEqual(residues["2"]["resname"], "V")
        self.assertCountEqual(residues["1"]["beads"], ["CA1", "L1"])
        self.assertCountEqual(residues["2"]["beads"], ["CA2", "V2"])

    def test_2beadAA_excludes_CA(self):
        """
        CA beads must not overwrite the residue name.
        """
        topology = {"L1": {},
                    "CA1": {}}
        residues = hf.get_residues_from_topology_dict(topology, model="2beadAA")
        self.assertEqual(residues["1"]["resname"], "L")

    def test_2beadAA_only_CA_is_glycine(self):
        """
        Residues containing only CA beads are assigned glycine ('G').
        """
        topology = {"CA1": {},
                    "CA2": {}}
        residues = hf.get_residues_from_topology_dict(topology, model="2beadAA")
        self.assertEqual(residues["1"]["resname"], "G")
        self.assertEqual(residues["2"]["resname"], "G")
        self.assertEqual(residues["1"]["beads"], ["CA1"])
        self.assertEqual(residues["2"]["beads"], ["CA2"])

    def test_invalid_model_raises(self):
        """
        Unknown protein model must raise ValueError.
        """
        topology = {"A1": {}}
        with self.assertRaises(ValueError):
            hf.get_residues_from_topology_dict(topology, model="3beadAA")

    def test_invalid_bead_id_raises(self):
        """
        Bead identifiers without a numeric residue index must raise ValueError.
        """
        topology = {"CA": {},   # no index
                    "L1": {}}
        with self.assertRaises(ValueError):
            hf.get_residues_from_topology_dict(topology, model="1beadAA")

if __name__ == "__main__":
    ut.main()
