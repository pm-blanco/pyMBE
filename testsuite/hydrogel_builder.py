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
import unittest as ut
import pyMBE
from pyMBE.lib.lattice import DiamondLattice
import espressomd

pmb = pyMBE.pymbe_library(seed=42)

# Define node particle
NodeType = "node_type"
pmb.define_particle(name=NodeType, 
                    sigma=0.355*pmb.units.nm, 
                    epsilon=1*pmb.units('reduced_energy'))

# define monomers
BeadType1 = "C"
pmb.define_particle(name=BeadType1, 
                    sigma=0.355*pmb.units.nm, 
                    epsilon=1*pmb.units('reduced_energy'))
BeadType2 = "M"
pmb.define_particle(name=BeadType2, 
                    sigma=0.355*pmb.units.nm, 
                    epsilon=1*pmb.units('reduced_energy'))

Res1 = "res_1"
pmb.define_residue(
    name=Res1,  # Name of the residue
    central_bead=BeadType1,  # Define the central bead name
    side_chains=[]  # Assuming no side chains for the monomer
)

Res2 = "res_2"
pmb.define_residue(
    name=Res2,  # Name of the residue
    central_bead=BeadType2,  # Define the central bead name
    side_chains=[]  # Assuming no side chains for the monomer
)

molecule_name = 'alternating_residue'
mpc=8
residue_list = [Res1]*(mpc//2) + [Res2]*(mpc//2)
pmb.define_molecule(name=molecule_name, residue_list = residue_list)

# define bond parameters
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
generic_bond_length = 0.355*pmb.units.nm
HARMONIC_parameters = {'r_0'    : generic_bond_length,
                      'k'      : generic_harmonic_constant}

pmb.define_bond(bond_type = 'harmonic',
                bond_parameters = HARMONIC_parameters, particle_pairs = [[BeadType1, BeadType1],
                                                                        [BeadType1, BeadType2],
                                                                        [BeadType2, BeadType2],
                                                                        [NodeType, BeadType1],
                                                                        [NodeType, BeadType2]])

diamond_lattice = DiamondLattice(mpc, generic_bond_length)
box_l = [diamond_lattice.box_l]*3
espresso_system = espressomd.System(box_l = box_l)
lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)


pmb.create_molecule(name=molecule_name,
                    number_of_molecules=1,
                    box_l=box_l,
                    use_default_bond=False,
                    list_of_first_residue_positions = [[np.random.uniform(0,diamond_lattice.box_l)]*3])

# Setting up node topology
indices = diamond_lattice.indices
node_topology = []

for index in range(len(indices)):
    node_topology.append({"particle_name": NodeType,
                          "lattice_index": indices[index]})

# Setting up chain topology
chain_topology = []
for node_conectivity in diamond_lattice.connectivity:
    node_start = str(diamond_lattice.indices[node_conectivity[0]])
    node_end = str(diamond_lattice.indices[node_conectivity[1]])
    chain_topology.append({'node_start':node_start,
                           'node_end': node_end,
                           'molecule_name':molecule_name})

#######################################################
hydrogel_name="my_hydrogel"
pmb.define_hydrogel(hydrogel_name,node_topology, chain_topology)

# Creating hydrogel
hydrogel_id= pmb.create_hydrogel(hydrogel_name,box_l=box_l)
hydrogel_tpl = pmb.db.get_template(pmb_type="hydrogel", 
                                   name=hydrogel_name)
hydrogel_inst = pmb.db.get_instance(pmb_type="hydrogel", 
                                    instance_id=hydrogel_id)

pmb.set_simulation_engine(espresso_system)
pmb.add_instances_to_engine()

class Test(ut.TestCase):
    def test_create_hydrogel_missing_template(self):
        """
        Unit test that create_hydrogel raises if the template is missing.
        """
        with self.assertRaises(ValueError):
            pmb.create_hydrogel("missing_hydrogel_template",box_l=box_l)

    def test_hydrogel_template_storage(self):
        """
        Unit test that checks that the hydrogel input information
        (node_map and chain_map) is correctly stored in the pyMBE database.
        """

        hydrogel_tpl = pmb.db.get_template(pmb_type="hydrogel",
                                           name=hydrogel_name)
        # --- Test node_map storage ---
        self.assertEqual(len(hydrogel_tpl.node_map), 
                         len(node_topology))
        # Convert both representations to comparable sets
        expected_nodes = {(node["particle_name"], tuple(node["lattice_index"])) for node in node_topology}
        stored_nodes = {(node.particle_name, tuple(node.lattice_index)) for node in hydrogel_tpl.node_map}
        self.assertSetEqual(stored_nodes,
                            expected_nodes,
                            "Stored hydrogel node_map does not match input definition")
        # --- Test chain_map storage ---
        self.assertEqual(len(hydrogel_tpl.chain_map), len(chain_topology))
        expected_chains = {(chain["node_start"], chain["node_end"], chain["molecule_name"]) for chain in chain_topology}
        stored_chains = {(chain.node_start, chain.node_end, chain.molecule_name) for chain in hydrogel_tpl.chain_map}
        self.assertSetEqual(stored_chains, expected_chains, "Stored hydrogel chain_map does not match input definition")
   
    def test_hydrogel_instance_info(self):
        """
        Unit test to check that hydrogel instance store information properly
        """
        self.assertEqual(hydrogel_inst.name, hydrogel_name)
        self.assertEqual(hydrogel_inst.assembly_id, hydrogel_id)  

    def test_node_positions(self):
        """
        Unit test that checks that nodes are created in the right position
        """
        hydrogel_tpl = pmb.db.get_template(pmb_type="hydrogel", name=hydrogel_name)
        # Get all particles belonging to this hydrogel
        particle_ids = pmb.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                            attribute="assembly_id",
                                                            value=hydrogel_id)
        node_particles = {pid: pmb.db.get_instance("particle", pid) for pid in particle_ids  if pmb.db.get_instance("particle", pid).name == NodeType}
        self.assertEqual(len(node_particles), len(hydrogel_tpl.node_map))
        for node_tpl in hydrogel_tpl.node_map:
            node_index = np.array(node_tpl.lattice_index)
            expected_pos = node_index * 0.25 * diamond_lattice.box_l
            found = False
            for pid, inst in node_particles.items():
                pos = espresso_system.part.by_id(pid).pos
                if np.allclose(pos, expected_pos, atol=1e-7):
                    self.assertEqual(inst.name, node_tpl.particle_name)
                    found = True
                    break

            self.assertTrue(found, f"Node at {node_index} not found")

    def test_chain_creation(self):
        """
        Unit test that checks that the chains are created as defined in the hydrogel template in the database.
        """
        hydrogel_tpl = pmb.db.get_template(pmb_type="hydrogel", name=hydrogel_name)
        molecule_ids = pmb.db._find_instance_ids_by_attribute(pmb_type="molecule",
                                                            attribute="assembly_id",
                                                            value=hydrogel_id)
        self.assertEqual(len(molecule_ids), len(hydrogel_tpl.chain_map))
        for chain_tpl in hydrogel_tpl.chain_map:
            expected_name = chain_tpl.molecule_name
            found = False
            for mol_id in molecule_ids:
                mol = pmb.db.get_instance("molecule", mol_id)
                if mol.name == expected_name:
                    found = True
                    break
            self.assertTrue(found, f"Chain {expected_name} not found")
    
    def test_chain_length(self):
        """
        Unit test to test that chains are created in the right position
        """
        molecule_ids = pmb.db._find_instance_ids_by_attribute(pmb_type="molecule",
                                                            attribute="assembly_id",
                                                            value=hydrogel_id)
        expected = (diamond_lattice.mpc -1) * generic_bond_length.m_as('reduced_length')       
        for mol_id in molecule_ids:
            particle_ids = pmb.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                                attribute="molecule_id",
                                                                value=mol_id)
            positions = np.array([espresso_system.part.by_id(pid).pos for pid in particle_ids])
            contour = np.sum(np.linalg.norm(np.diff(positions, axis=0), axis=1))*pmb.units.nm
            np.testing.assert_allclose(contour.m_as('reduced_length'), expected, atol=1e-7)
    
    def test_exceptions(self):
        """
        Unit tests for the sanity tests
        """
        #  check that only non-negative values of monomers per chain are allowed 
        np.testing.assert_raises(ValueError, DiamondLattice, 0, generic_bond_length)
        np.testing.assert_raises(ValueError, DiamondLattice, "invalid", generic_bond_length)
        np.testing.assert_raises(ValueError, DiamondLattice, -5, generic_bond_length)
        # check that any objects are other than DiamondLattice passed to initialize_lattice_builder raises a TypeError 
        np.testing.assert_raises(TypeError, pmb.initialize_lattice_builder, None)
        # Check exceptions when the node and chain maps are incomplete
        incomplete_node_map = [{"particle_name": NodeType, "lattice_index": [0, 0, 0]},{"particle_name": NodeType, "lattice_index": [1, 1, 1]}]
        incomplete_chain_map = [{"node_start": "[0 0 0]", "node_end":"[1 1 1]" , "residue_list": residue_list}]
        np.testing.assert_raises(ValueError, pmb.define_hydrogel, "test_hydrogel", incomplete_node_map, chain_topology)
        np.testing.assert_raises(ValueError, pmb.define_hydrogel, "test_hydrogel", node_topology, incomplete_chain_map)

if __name__ == "__main__":
    ut.main()
