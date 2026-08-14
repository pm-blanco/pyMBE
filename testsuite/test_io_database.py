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

import tempfile
import espressomd
import pandas as pd
import unittest as ut
import json
import os
import pyMBE
from pyMBE.lib.lattice import DiamondLattice
import pyMBE.lib.handy_functions as hf
from pyMBE.storage.io import _decode, _encode, _load_database_csv, _save_database_csv
from pyMBE.storage.pint_quantity import PintQuantity
from pyMBE.storage.instances.bond import BondInstance
from pathlib import Path
import csv
import numpy as np

box_l=[100]*3
espresso_system=espressomd.System (box_l = box_l)

class DummyDB:
    def __init__(self):
        self._templates = {}
        self._instances = {}
        self._reactions = {}

class Test(ut.TestCase):

    def test_instance_fallback_model_dump_failure(self):
        """
        Tests database behavior in failure of instance model dump
        """
        class BadInstance:
            name = "bad_inst"

        db = DummyDB()
        db._templates = {}
        db._instances["weird"] = {"x": BadInstance()}
        db._reactions = {}
        with tempfile.TemporaryDirectory() as tmp:
            _save_database_csv(db, tmp)

            text = Path(tmp, "instances_weird.csv").read_text()
            self.assertIn("bad_inst", text)

    def test_template_fallback_model_dump_failure(self):
        """
        Tests database behavior in failure of template model dump
        """
        class BadTemplate:
            name = "bad"

        db = DummyDB()
        db._templates["weird"] = {"bad": BadTemplate()}
        db._instances = {}
        db._reactions = {}
        with tempfile.TemporaryDirectory() as tmp:
            _save_database_csv(db, tmp)
            text = Path(tmp, "templates_weird.csv").read_text()
            self.assertIn("bad", text)

    def test_invalid_metadata_json(self):
        """
        Tests that invalid metadata files in the database are ignored
        """
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "metadata.json").write_text('"this is a string"')
            db = DummyDB()
            metadata = _load_database_csv(db, folder)
            self.assertEqual(metadata, {})

    def test_lj_shift_as_pint_quantity(self):
        """
        Tests a LJ shift as a pint quantity 
        """
        sigma = {"magnitude": 1.0,
                "units": "nm",
                "dimension": "[length]"}
        cutoff = {"magnitude": 1.0,
                "units": "nm",
                "dimension": "[length]"}
        offset = {"magnitude": 1.0,
                "units": "nm",
                "dimension": "[length]"}
        epsilon = {"magnitude": 1.0,
                "units": "J",
                "dimension": "[energy]"}
        shift = {"magnitude": 1.0,
                "units": "nm",
                "dimension": "[length]"}

        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            with open(folder / "templates_lj.csv", "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    ["name", "state1", "state2", "sigma", "epsilon", "cutoff", "offset", "shift"]
                )
                writer.writerow([
                    "lj1", "A", "B",
                    json.dumps(sigma),
                    json.dumps(epsilon),
                    json.dumps(cutoff),
                    json.dumps(offset),
                    json.dumps(shift),
                ])


            db = DummyDB()
            _load_database_csv(db, folder)

            lj = db._templates["lj"]["A-B"]
            self.assertIsInstance(lj.shift, PintQuantity)

    def test_bond_template_scalar_parameter(self):
        """
        Tests Bond template with non-PintQuantity parameter
        """
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "templates_bond.csv").write_text(
                "name,bond_type,parameters\n"
                'b1,harmonic,"{""k"":{""magnitude"":24,""units"":""kilojoule / nm**2"",""dimension"":""energy/length**2""}}"\n'
            )

            db = DummyDB()
            _load_database_csv(db, folder)

            bond = db._templates["bond"]["b1"]
            self.assertEqual(bond.parameters["k"].magnitude, 24)


    def test_lists_string_coerced_to_list(self):
        """
        Tests that string lists returned by the encoder are parsed back to list properly
        """
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)

            (folder / "templates_residue.csv").write_text(
                "name,central_bead,side_chains\n"
                "RES1,BB,XYZ\n"
            )

            db = DummyDB()
            _load_database_csv(db, folder)

            tpl = db._templates["residue"]["RES1"]
            self.assertEqual(tpl.side_chains, ["X", "Y", "Z"])

        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "templates_molecule.csv").write_text("name,residue_list\n"
                                                            "MOL1,ABC\n")

            db = DummyDB()
            _load_database_csv(db, folder)

            tpl = db._templates["molecule"]["MOL1"]
            self.assertEqual(tpl.residue_list, ["A", "B", "C"])

        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)

            (folder / "templates_peptide.csv").write_text(
                "name,model,residue_list,sequence\n"
                "PEP1,CG,XYZ,XYZ\n"
            )

            db = DummyDB()
            _load_database_csv(db, folder)

            tpl = db._templates["peptide"]["PEP1"]
            self.assertEqual(tpl.residue_list, ["X", "Y", "Z"])

        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)

            (folder / "templates_protein.csv").write_text(
                "name,model,residue_list,sequence\n"
                "PROT1,CG,DEF,DEF\n"
            )

            db = DummyDB()
            _load_database_csv(db, folder)

            tpl = db._templates["protein"]["PROT1"]
            self.assertEqual(tpl.residue_list, ["D", "E", "F"])

    def test_io_particles_and_particle_states_templates(self):
        """
        Checks that information in the pyMBE database about 
        particle and particle_state templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        units = pmb.units
        pmb.define_particle(name="Z",
                    sigma=3.5 * units.reduced_length,
                    cutoff=4 * units.reduced_length,
                    offset=0 * units.reduced_length,
                    epsilon=0.2 * units.reduced_energy,
                    acidity="acidic",
                    pka=4.25)
        pmb.define_particle(name="X",
                        sigma=3.5 * units.reduced_length,
                        cutoff=4 * units.reduced_length,
                        offset=0 * units.reduced_length,
                        epsilon=0.2 * units.reduced_energy,
                        z=1)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded particle_state templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="particle_state"),
                                      new_pmb.get_templates_df(pmb_type="particle_state"))
        # Test that the loaded particle templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="particle"),
                                      new_pmb.get_templates_df(pmb_type="particle"))
        
    def test_io_lj_templates(self):
        """
        Checks that information in the pyMBE database about 
        LennardJOnes templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        units = pmb.units
        pmb.define_particle(name="Z",
                            sigma=3.5 * units.reduced_length,
                            cutoff=4 * units.reduced_length,
                            offset=0 * units.reduced_length,
                            epsilon=0.2 * units.reduced_energy,
                            acidity="acidic",
                            pka=4.25)
        pmb.define_particle(name="X",
                            sigma=3.5 * units.reduced_length,
                            cutoff=4 * units.reduced_length,
                            offset=0 * units.reduced_length,
                            epsilon=0.2 * units.reduced_energy,
                            z=1)
        pmb.set_simulation_engine(espresso_system)
        pmb.setup_lj_interactions()
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded LJ templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="lj"),
                                      new_pmb.get_templates_df(pmb_type="lj"))
    
    def test_io_bond_templates(self):
        """
        Checks that information in the pyMBE database about 
        residue templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        parameters1 = {"k":  100.0 * pmb.units.reduced_energy / (pmb.units.reduced_length**2),
                       "r_0": 1.0  * pmb.units.reduced_length}
        parameters2 = {'r_0'    : 0.4 * pmb.units.nm,
                        'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                        'd_r_max': 0.8 * pmb.units.nm}
        pmb.define_bond(bond_type="harmonic",
                    bond_parameters=parameters1,
                    particle_pairs=[["Z","Z"], 
                                    ["Z","X"],
                                    ["X","X"]])
        pmb.define_default_bond(bond_type="FENE",
                                bond_parameters=parameters2)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded bond templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="bond"),
                                      new_pmb.get_templates_df(pmb_type="bond"))

    def test_io_angle_templates(self):
        """
        Checks that information in the pyMBE database about
        angle templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        parameters1 = {"k": 50.0 * pmb.units.reduced_energy,
                       "phi_0": 1.0 * pmb.units("")}
        parameters2 = {"k": 30.0 * pmb.units.reduced_energy,
                       "phi_0": 0.5 * pmb.units("")}
        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=parameters1,
                         particle_triplets=[("Z", "X", "Z"),
                                            ("X", "Z", "Y")])
        pmb.define_default_angular_potential(angle_type="cosine",
                                 angle_parameters=parameters2)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded angle templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="angle"),
                                      new_pmb.get_templates_df(pmb_type="angle"))

    def test_io_residues_templates(self):
        """
        Checks that information in the pyMBE database about 
        residue templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_residue(name="R1", 
                           central_bead="Z", 
                           side_chains=["X","Z"])
        pmb.define_residue(name="R2", 
                           central_bead="X", 
                           side_chains=[])
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded residue templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="residue"),
                                      new_pmb.get_templates_df(pmb_type="residue"))
        
    def test_io_molecule_templates(self):
        """
        Checks that information in the pyMBE database about 
        molecule templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_molecule(name="M1", 
                            residue_list=["R1","R2"]*2)
        pmb.define_molecule(name="M2", 
                            residue_list=["R2","R1"]*20)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded molecule templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="molecule"),
                                      new_pmb.get_templates_df(pmb_type="molecule"))
        
    def test_io_peptide_templates(self):
        """
        Checks that information in the pyMBE database about 
        peptide templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_peptide(name="Peptide1",
                           model="1beadAA",
                           sequence="KKKKDDDD")
        pmb.define_peptide(name="Peptide2",
                           model="2beadAA",
                           sequence="RRRREEEE")
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded peptide templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="peptide"),
                                      new_pmb.get_templates_df(pmb_type="peptide"))
        
    def test_io_protein_templates(self):
        """
        Checks that information in the pyMBE database about 
        protein templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.define_protein(name="1beb",
                            model="2beadAA",
                            sequence="KKKKKK")
        pmb.define_protein(name="Avidin",
                            model="1beadAA",
                            sequence="KKKKKK")
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded protein templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="protein"),
                                      new_pmb.get_templates_df(pmb_type="protein"))
        
    def test_io_hydrogel_templates(self):
        """
        Checks that information in the pyMBE database about 
        hydrogel templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        diamond_lattice = DiamondLattice(4, 3.5 * pmb.units.reduced_length)
        lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)
        # Setting up node topology
        indices = diamond_lattice.indices
        node_topology = []
        for index in range(len(indices)):
            node_topology.append({"particle_name": "X",
                                "lattice_index": indices[index]})
        # Setting up chain topology
        node_labels = lattice_builder.node_labels
        chain_labels = lattice_builder.chain_labels
        reverse_node_labels = {v: k for k, v in node_labels.items()}
        chain_topology = []

        for chain_data in chain_labels.items():
            node_label_pair = chain_data[0]
            node_label_s, node_label_e = [int(x) for x in node_label_pair.strip("()").split(",")]
            chain_topology.append({'node_start':reverse_node_labels[node_label_s],
                                    'node_end': reverse_node_labels[node_label_e],
                                    'molecule_name':"M1"})

        pmb.define_hydrogel(name="my_hydrogel", 
                            node_map=node_topology, 
                            chain_map=chain_topology)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded hydrogel templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_templates_df(pmb_type="hydrogel"),
                                      new_pmb.get_templates_df(pmb_type="hydrogel"))

    def test_io_reaction_templates(self):
        """
        Checks that information in the pyMBE database about 
        reaction templates is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        path_to_pka=pmb.root / "parameters" / "pka_sets" / "Nozaki1967.json"
        pmb.load_pka_set(filename=path_to_pka)
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded protein templates are equal to the original
        pd.testing.assert_frame_equal(pmb.get_reactions_df(),
                                      new_pmb.get_reactions_df())
        
    def test_io_instances(self):
        """
        Checks that information in the pyMBE database about 
        instances created into espresso is stored to file and loaded back properly.
        """
        pmb = pyMBE.pymbe_library(seed=42)
        # Test instances of a hydrogel plus a standalone residue with an angle
        # (tests hydrogel, molecule, residue, bond, angle and particle instances)
        pmb.define_particle(name="Z",
                    sigma=3.5 * pmb.units.reduced_length,
                    cutoff=4  * pmb.units.reduced_length,
                    offset=0  * pmb.units.reduced_length,
                    epsilon=0.2 * pmb.units.reduced_energy,
                    acidity="acidic",
                    pka=4.25)
        pmb.define_particle(name="X",
                            sigma=3.5 * pmb.units.reduced_length,
                            cutoff=4 * pmb.units.reduced_length,
                            offset=0 * pmb.units.reduced_length,
                            epsilon=0.2 * pmb.units.reduced_energy,
                            z=1)
        pmb.define_residue(name="R1", 
                   central_bead="Z", 
                   side_chains=["X","Z"])
        parameters = {"k":  100.0 * pmb.units.reduced_energy / (pmb.units.reduced_length**2),
                      "r_0": 1.0  * pmb.units.reduced_length}
        pmb.define_bond(bond_type="harmonic",
                        bond_parameters=parameters,
                        particle_pairs=[["Z","Z"], 
                                        ["Z","X"],
                                        ["X","X"]])
        angle_parameters = {"k": 50.0 * pmb.units.reduced_energy,
                            "phi_0": 1.0 * pmb.units("")}
        pmb.define_angular_potential(angle_type="harmonic",
                         angle_parameters=angle_parameters,
                         particle_triplets=[("X", "Z", "Z")])
        pmb.define_molecule(name="M1", 
                            residue_list=["R1"]*1)
        angle_residue_id = pmb.create_residue(name="R1",
                                              box_l=box_l,
                                              gen_angle=True)
        diamond_lattice = DiamondLattice(4, 3.5 * pmb.units.reduced_length)
        lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)
        # Setting up node topology --> Nodes are particles of type "X"
        indices = diamond_lattice.indices
        node_topology = []
        for index in range(len(indices)):
            node_topology.append({"particle_name": "X",
                                "lattice_index": indices[index]})
        # Setting up chain topology --> Chains are molecules of type "M1"
        node_labels = lattice_builder.node_labels
        chain_labels = lattice_builder.chain_labels
        reverse_node_labels = {v: k for k, v in node_labels.items()}
        chain_topology = []
        for chain_data in chain_labels.items():
            node_label_pair = chain_data[0]
            node_label_s, node_label_e = [int(x) for x in node_label_pair.strip("()").split(",")]
            chain_topology.append({'node_start': reverse_node_labels[node_label_s],
                                    'node_end': reverse_node_labels[node_label_e],
                                    'molecule_name': "M1"})
        pmb.define_hydrogel("my_hydrogel", 
                            node_topology, 
                            chain_topology)
        assembly_id = pmb.create_hydrogel(name="my_hydrogel",
                                          box_l=box_l)

        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded instances are equal to the original
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="hydrogel"),
                                      new_pmb.get_instances_df(pmb_type="hydrogel"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="molecule"),
                                      new_pmb.get_instances_df(pmb_type="molecule"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="residue"),
                                      new_pmb.get_instances_df(pmb_type="residue"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="bond"),
                                      new_pmb.get_instances_df(pmb_type="bond"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="angle"),
                                      new_pmb.get_instances_df(pmb_type="angle"))
        instances_before=pmb.get_instances_df(pmb_type="particle").copy()
        positions_before=instances_before['position'].to_numpy()
        instances_after=new_pmb.get_instances_df(pmb_type="particle").copy()
        positions_after=instances_after['position'].to_numpy()
        pd.testing.assert_frame_equal(instances_before.drop(columns=['position']),
                                      instances_after.drop(columns=['position']))
        for i in range(len(positions_after)):
            self.assertTrue(np.allclose(positions_after[i],positions_before[i],atol=1e-12))
      
        # Clean up before the next test
        pmb.delete_instances_in_system(
                                       instance_id=assembly_id,
                                       pmb_type="hydrogel")
        pmb.delete_instances_in_system(
                                       instance_id=angle_residue_id,
                                       pmb_type="residue")
        pmb.db.delete_templates(pmb_type="angle")
        # Test instances of a peptide (tests peptide, residue, bond and particle instances)
        path_to_interactions=pmb.root / "parameters" / "peptides" / "Lunkad2021"
        path_to_pka=pmb.root / "parameters" / "pka_sets" / "Hass2015.json"
        ### pmb.load_database(folder): Returns metadata, but not used ?
        pmb.load_database (folder=path_to_interactions) # Defines particles  
        pmb.load_pka_set(filename=path_to_pka) 
        pka_set = pmb.get_pka_set()
        for particle_name in pka_set.keys():
            pmb.define_monoprototic_particle_states(particle_name=particle_name,
                                                    acidity=pka_set[particle_name]["acidity"])
        hf.define_peptide_AA_residues(sequence="KKKDDD",
                                      model="1beadAA",
                                      pmb=pmb)
        parameters2 = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                'd_r_max': 0.8 * pmb.units.nm}
        pmb.define_default_bond(bond_type="FENE",
                                bond_parameters=parameters2)
        pmb.define_peptide(name="Peptide1",
                            model="1beadAA",
                            sequence="KKKKDDDD")
     
        pep_ids = pmb.create_molecule(name="Peptide1",
                            number_of_molecules=2,
                            box_l=box_l, ###set box_l
                            use_default_bond=True)
        
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()

        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded instances are equal to the original
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="peptide"),
                                      new_pmb.get_instances_df(pmb_type="peptide"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="residue"),
                                      new_pmb.get_instances_df(pmb_type="residue"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="bond"),
                                      new_pmb.get_instances_df(pmb_type="bond"))
        
        instances_before=pmb.get_instances_df(pmb_type="particle")
        positions_before=instances_before['position'].to_numpy()
        instances_after=new_pmb.get_instances_df(pmb_type="particle")
        positions_after=instances_after['position'].to_numpy()

        pd.testing.assert_frame_equal(instances_before.drop(columns=['position']),
                                      instances_after.drop(columns=['position']))
        for i in range(len(positions_after)):
            self.assertTrue(np.allclose(positions_after[i],positions_before[i],atol=1e-12))
        # Clean up before the next test
        for pepid in pep_ids:
            pmb.delete_instances_in_system(
                                        instance_id=pepid,
                                        pmb_type="peptide")
        pmb.db.delete_templates(pmb_type="particle")
        pmb.db.delete_templates(pmb_type="particle_state")
        pmb.db.delete_templates(pmb_type="residue")
        pmb.db.delete_reactions()
        # Test instances of a protein (tests protein, residue and particle instances)
        path_to_protein_structure = pmb.root / "parameters" / "globular_proteins" / "1beb.vtf"    
        topology_dict, sequence = pmb.read_protein_vtf (filename=path_to_protein_structure)
        pmb.load_pka_set(filename=path_to_pka)
        # Define AA particles and residues
        hf.define_protein_AA_particles(topology_dict=topology_dict,
                                       pmb=pmb,
                                       pka_set=pka_set)
        hf.define_protein_AA_residues(sequence=sequence,
                                      model="2beadAA",
                                      pmb=pmb)
        pmb.define_protein(name="1beb",
                           model="2beadAA",
                           sequence="KKKKKK")
        prot_ids = pmb.create_protein(name="1beb",
                                    number_of_proteins=1,
                                    box_l=box_l,
                                    topology_dict=topology_dict)
        
        pmb.set_simulation_engine(espresso_system)
        pmb.add_instances_to_engine()
        
        new_pmb = pyMBE.pymbe_library(23)
        with tempfile.TemporaryDirectory() as tmp_directory:
            # Save and load the database
            pmb.save_database(tmp_directory)
            new_pmb.load_database(tmp_directory)
        # Test that the loaded instances are equal to the original
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="protein"),
                                      new_pmb.get_instances_df(pmb_type="protein"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="residue"),
                                      new_pmb.get_instances_df(pmb_type="residue"))
        pd.testing.assert_frame_equal(pmb.get_instances_df(pmb_type="bond"),
                                      new_pmb.get_instances_df(pmb_type="bond"))
        instances_before=pmb.get_instances_df(pmb_type="particle")
        positions_before=instances_before['position'].to_numpy()
        instances_after=new_pmb.get_instances_df(pmb_type="particle")
        positions_after=instances_after['position'].to_numpy()
        pd.testing.assert_frame_equal(instances_before.drop(columns=['position']),
                                      instances_after.drop(columns=['position']))
        
        for i in range(len(positions_after)):
            self.assertTrue(np.allclose(positions_after[i],positions_before[i],atol=1e-12))
        # Clean up 
        for protid in prot_ids:
            pmb.delete_instances_in_system(
                                        instance_id=protid,
                                        pmb_type="protein")
            
    def test_database_io_exceptions(self):
        """
        Unit test to check exceptions in the io of the pyMBE database
        """
        pmb = pyMBE.pymbe_library(51)
        inputs = {"folder": "test",
                  "format": "random"}
        self.assertRaises(ValueError,
                          pmb.load_database,
                          **inputs)
        self.assertRaises(ValueError,
                          pmb.save_database,
                          **inputs)

    def test_decode_edge_cases(self):
        """
        Tests the edge cases in the IO decoder of the database
        """
        self.assertIsNone(_decode(None))
        self.assertIsNone(_decode(float("nan")))
        self.assertIsNone(_decode(""))
        self.assertIsNone(_decode("nan"))

        # malformed JSON → fallback to raw string
        self.assertEqual(_decode("{not:json}"), "{not:json}")

        # already-native types
        self.assertEqual(_decode({"a": 1}), {"a": 1})
        self.assertEqual(_decode([1, 2]), [1, 2])
        self.assertEqual(_decode(3), 3)

        value = 3.14159
        result = _decode(value)

        self.assertIsInstance(result, float)
        self.assertEqual(result, value)

        value = (1, 2, 3)   # tuple is not dict, list, int, bool, float, or str
        result = _decode(value)
        self.assertIsNone(result)

    def test_encode_edge_cases(self):
        """
        Tests the edge cases in the IO encoder of the database
        """
        

        self.assertEqual(_encode(None), "")

        pq = PintQuantity(magnitude=3.0, units="nm", dimension="length")
        encoded = _encode(pq)
        self.assertIsInstance(encoded, str)
        self.assertIn("magnitude", encoded)

        class Dummy:
            def __str__(self):
                return "dummy"

        self.assertEqual(_encode(Dummy()), json.dumps("dummy"))

    def test_load_empty_database_folder(self):
        """
        Tests that an empty folder does not populate the pyMBE database
        """
        new_pmb = pyMBE.pymbe_library(2)
        with tempfile.TemporaryDirectory() as tmp:
            new_pmb.load_database(tmp)
        # database should remain empty
        self.assertEqual(len(new_pmb.db._templates), 0)
        self.assertEqual(len(new_pmb.db._instances), 0)
        self.assertEqual(len(new_pmb.db._reactions), 0)

    def test_partial_database_files(self):
        """
        Test that the database does not break if a file is missing
        """
        pmb = pyMBE.pymbe_library(1)
        pmb.define_particle(name="X", 
                            sigma=1*pmb.units.reduced_length,
                            epsilon=1*pmb.units.reduced_energy)
        with tempfile.TemporaryDirectory() as tmp:
            pmb.save_database(tmp)
            # manually delete one CSV
            os.remove(os.path.join(tmp, "templates_particle.csv"))
            new_pmb = pyMBE.pymbe_library(2)
            new_pmb.load_database(tmp)
        # particle templates missing, but no crash
        self.assertTrue(new_pmb.get_templates_df("particle").empty)

    def test_metadata_roundtrip(self):
        """
        Test that covers:
            - happy path of metadata reading
            - return value of load_database
        """
        pmb = pyMBE.pymbe_library(1)
        with tempfile.TemporaryDirectory() as tmp:
            pmb.save_database(tmp)
            meta = {"creator": "test", "version": 1}
            with open(os.path.join(tmp, "metadata.json"), "w") as f:
                json.dump(meta, f)
            new_pmb = pyMBE.pymbe_library(2)
            loaded_meta = new_pmb.load_database(tmp)
        self.assertEqual(loaded_meta, meta)

    def test_invalid_metadata_file(self):
        """
        Covers handling of broken metadata files
        """
        pmb = pyMBE.pymbe_library(1)

        with tempfile.TemporaryDirectory() as tmp:
            pmb.save_database(tmp)
            with open(os.path.join(tmp, "metadata.json"), "w") as f:
                f.write("not json")

            new_pmb = pyMBE.pymbe_library(2)
            meta = new_pmb.load_database(tmp)

        self.assertEqual(meta, {})

    def test_default_bond_particle_names(self):
        """
        Test io for default bonds
        """
        pmb = pyMBE.pymbe_library(1)
        pmb.define_default_bond(bond_type="FENE", 
                                bond_parameters={'r_0'    : 0.5 * pmb.units.nm,
                                                'k'      : 500 * pmb.units('reduced_energy / reduced_length**2'),
                                                'd_r_max': 0.5 * pmb.units.nm})

        with tempfile.TemporaryDirectory() as tmp:
            pmb.save_database(tmp)
            new_pmb = pyMBE.pymbe_library(2)
            new_pmb.load_database(tmp)

        df = new_pmb.get_templates_df("bond")
        self.assertTrue(df["particle_name1"].isna().any())

    def test_non_sequential_instance_ids(self):
        """
        Tests that  IDs are not implicitly re-indexed.
        """
        pmb = pyMBE.pymbe_library(1)
        pmb.db._instances["bond"] = {10: BondInstance(name="b", bond_id=10, particle_id1=1, particle_id2=2),
                                     42: BondInstance(name="b", bond_id=42, particle_id1=3, particle_id2=4)}
        with tempfile.TemporaryDirectory() as tmp:
            pmb.save_database(tmp)
            new_pmb = pyMBE.pymbe_library(2)
            new_pmb.load_database(tmp)
        self.assertSetEqual(set(new_pmb.db._instances["bond"].keys()), {10, 42},)

    def test_load_database_missing_folder(self):
        """
        Tests that that the io raises an error if the folder does not exist
        """
        pmb = pyMBE.pymbe_library(1)
        with self.assertRaises(FileNotFoundError):
            pmb.load_database("does_not_exist")

    

if __name__ == '__main__':
    ut.main()
