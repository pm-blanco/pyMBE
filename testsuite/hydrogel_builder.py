#
# Copyright (C) 2024-2025 pyMBE-dev team
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

CounterIon = "counter_ion"
pmb.define_particle(name=CounterIon, 
                    sigma=0.5*pmb.units.nm, 
                    epsilon=1.5*pmb.units("reduced_energy"))

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
pmb.define_molecule(name=molecule_name, residue_list = [Res1, Res2, Res1, Res2, Res1, Res2])

# define bond parameters
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
generic_bond_length = 0.355*pmb.units.nm
HARMONIC_parameters = {'r_0'    : generic_bond_length,
                      'k'      : generic_harmonic_constant}

pmb.define_bond(bond_type = 'harmonic',
                bond_parameters = HARMONIC_parameters, particle_pairs = [[BeadType1, BeadType1],
                                                                        [BeadType1, BeadType2],
                                                                        [BeadType2, BeadType2]])
pmb.define_bond(bond_type = 'harmonic',
                bond_parameters = HARMONIC_parameters, particle_pairs = [[NodeType, BeadType1],
                                                                        [NodeType, BeadType2]])

mpc=8
diamond_lattice = DiamondLattice(mpc, generic_bond_length)
box_l = diamond_lattice.box_l
espresso_system = espressomd.System(box_l = [box_l]*3)
pmb.add_bonds_to_espresso(espresso_system = espresso_system)
lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)

pmb.create_particle(name=CounterIon,
                    espresso_system=espresso_system,
                    number_of_particles=1,
                    position=[[np.random.uniform(0,box_l)]*3])

pmb.create_molecule(name=molecule_name,
                    number_of_molecules=1,
                    espresso_system=espresso_system,
                    use_default_bond=False,
                    list_of_first_residue_positions = [[np.random.uniform(0,box_l)]*3])

# Setting up node topology
indices = diamond_lattice.indices
node_topology = []

for index in range(len(indices)):
    node_topology.append({"particle_name": NodeType,
                          "lattice_index": indices[index]})

# Setting up chain topology
connectivity = diamond_lattice.connectivity
node_labels = lattice_builder.node_labels
reverse_node_labels = {v: k for k, v in node_labels.items()}
connectivity_with_labels = {(reverse_node_labels[i], reverse_node_labels[j]) for i, j in connectivity}
chain_topology = []
residue_list = [Res1]*(mpc//2) + [Res2]*(mpc//2)
for node_s, node_e in connectivity_with_labels:
    chain_topology.append({'node_start':node_s,
                              'node_end': node_e,
                              'residue_list':residue_list})
#######################################################
hydrogel_name="my_hydrogel"
pmb.define_hydrogel(hydrogel_name,node_topology, chain_topology)

# Creating hydrogel
hydrogel_info = pmb.create_hydrogel(hydrogel_name, espresso_system)

################################################

def compare_node_maps(map1, map2):
    # Ensure lengths are the same
    np.testing.assert_equal(len(map1), len(map2))
    
    # Sort lists by lattice_index to ensure correct comparison
    map1_sorted = sorted(map1, key=lambda x: tuple(x["lattice_index"]))
    map2_sorted = sorted(map2, key=lambda x: tuple(x["lattice_index"]))
    
    # Compare each node's details
    for node1, node2 in zip(map1_sorted, map2_sorted):
        np.testing.assert_equal(node1["particle_name"], 
                                node2["particle_name"])
        np.testing.assert_equal(node1["lattice_index"],
                                node2["lattice_index"])
        
    return 

def parse_string_to_array(string):
    """
    Convert a string representation of a list (e.g., '[3 1 3]') into a numpy array.
    """
    string = string.strip("[]")  # Remove brackets
    elements = map(int, string.split())  # Split by spaces and convert to integers
    return np.array(list(elements))

def compare_chain_maps(chain_topology_1, chain_topology_2):
    """
    Compare two chain topology maps by checking if they have the same set of edges with corresponding residue lists.
    """
    np.testing.assert_equal(len(chain_topology_1), len(chain_topology_2))

    # Convert string coordinates to arrays and sort lists by (node_start, node_end)
    def preprocess_chain(chain_topology):
        processed = []
        for edge in chain_topology:
            processed.append({
                'node_start': parse_string_to_array(edge['node_start']),
                'node_end': parse_string_to_array(edge['node_end']),
                'residue_list': edge['residue_list']  # Keep as is
            })
        return sorted(processed, key=lambda x: (x['node_start'].tolist(), x['node_end'].tolist()))

    chain_topology_1 = preprocess_chain(chain_topology_1)
    chain_topology_2 = preprocess_chain(chain_topology_2)

    # Compare edges one by one
    for edge1, edge2 in zip(chain_topology_1, chain_topology_2):
        np.testing.assert_equal(edge1['node_start'].tolist(), 
                                edge2['node_start'].tolist())
        np.testing.assert_equal(edge1['node_end'].tolist(), 
                                edge2['node_end'].tolist())
        # Check if the residue lists are the same
        np.testing.assert_equal(edge1['residue_list'], 
                                edge2['residue_list'])

    return # All edges match

class Test(ut.TestCase):
 
    def test_format_node(self):
        assert pmb.format_node([1, 2, 3]) == "[1 2 3]"
        assert pmb.format_node([4, 5, 6]) == "[4 5 6]"

    def test_hydrogel_info(self):
        assert hydrogel_info["name"] == hydrogel_name

    def test_node_positions(self):
        for _, node_id in hydrogel_info["nodes"].items():
            node_pos = espresso_system.part.by_id(int(node_id[0])).pos
            node_name_in_espresso = pmb.df[(pmb.df["pmb_type"] == "particle") & (pmb.df["particle_id"] == node_id[0])]["name"].values[0]
            node_label = node_labels[pmb.format_node(list((node_pos*(4/lattice_builder.box_l)).astype(int)))]
            node_data = node_topology[node_label]
            node_name = node_data["particle_name"] 
            # Assert node's name and position are correctly set
            np.testing.assert_equal(node_name_in_espresso, node_name)
            np.testing.assert_allclose(np.copy(node_pos), np.array(node_data["lattice_index"]) * 0.25 * diamond_lattice.box_l, atol=1e-7)

    def test_chain_placement_and_connectivity(self):
        for molecule_id, molecule_data in hydrogel_info["chains"].items():
            # Ensure that chain's node_start and node_end are correctly set
            node_start = molecule_data["node_start"]
            node_end = molecule_data["node_end"]
            chain_name_in_espresso = pmb.df[(pmb.df["pmb_type"] == "molecule") & (pmb.df["molecule_id"] == molecule_id)]["name"].values[0]
            # Assert chain's node_start and node_end
            np.testing.assert_equal(chain_name_in_espresso, f"chain_{node_start}_{node_end}")
            # Check if chain is connected in the espresso system (e.g., check bond or distance between node_start and node_end)
            node_start_id = hydrogel_info["nodes"][node_start][0]
            node_end_id = hydrogel_info["nodes"][node_end][0]
            start_pos = espresso_system.part.by_id(int(node_start_id)).pos
            end_pos = espresso_system.part.by_id(int(node_end_id)).pos
            vec_between_nodes = end_pos - start_pos
            # Ensure that the chain is connected (check distance, should be within acceptable bond length range)
            vec_between_nodes = vec_between_nodes -  diamond_lattice.box_l * np.round(vec_between_nodes / diamond_lattice.box_l)
            distance_between_nodes = np.linalg.norm(vec_between_nodes)
            np.testing.assert_allclose(distance_between_nodes, (diamond_lattice.mpc+1)*generic_bond_length.magnitude, atol=0.0000001)
    
    def test_all_residue_placement(self):
        def get_residue_list(chain_topology, node_start, node_end):
            node_start_array = parse_string_to_array(node_start)
            node_end_array = parse_string_to_array(node_end)
            for edge in chain_topology:
                if (np.array_equal(parse_string_to_array(edge['node_start']), node_start_array) and
                   np.array_equal(parse_string_to_array(edge['node_end']), node_end_array)):
                    return edge['residue_list']
        
        for _, chain_data in hydrogel_info["chains"].items():
            residue_in_chain = chain_data.copy()
            node_start = residue_in_chain.pop("node_start")
            node_end = residue_in_chain.pop("node_end")
            node_start_label = lattice_builder.node_labels[node_start]
            node_end_label = lattice_builder.node_labels[node_end]
            vec_between_nodes = (np.array([float(x) for x in node_end.strip('[]').split()]) -
                                     np.array([float(x) for x in node_start.strip('[]').split()])) * 0.25 * lattice_builder.box_l
            vec_between_nodes = vec_between_nodes - lattice_builder.box_l * np.round(vec_between_nodes / lattice_builder.box_l)
            backbone_vector = vec_between_nodes / (diamond_lattice.mpc + 1)
   
            for (res_id, res_data) in residue_in_chain.items():
                central_bead_id = res_data["central_bead_id"]
                
                # Get the position of the central bead from the espresso system 
                central_bead_pos = espresso_system.part.by_id(central_bead_id).pos

                # Calculate the expected position of the residue's central bead
                residue_index = list(residue_in_chain.keys()) .index(res_id)
                expected_position = np.array([float(x) for x in node_start.strip('[]').split()]) * 0.25 * diamond_lattice.box_l + (residue_index + 1) * backbone_vector
                
                # Validate that the central bead's position matches the expected position
                np.testing.assert_allclose(np.copy(central_bead_pos), expected_position, atol=1e-7)
                expected_node_start = reverse_node_labels[node_start_label]
                expected_node_end = reverse_node_labels[node_end_label]
                expected_res_name = get_residue_list(chain_topology, expected_node_start, expected_node_end)[residue_index]
                residue_name = pmb.df[(pmb.df["pmb_type"]=="residue") & (pmb.df["residue_id"]==res_id)]["name"].values[0]
                np.testing.assert_equal(node_start, expected_node_start)
                np.testing.assert_equal(node_end, expected_node_end)
                np.testing.assert_equal(residue_name, expected_res_name)

    def test_exceptions(self):
        print("*** Unit Test: check that only non-negative values of monomers per chain are allowed ***")
        np.testing.assert_raises(ValueError, DiamondLattice, 0, generic_bond_length)
        np.testing.assert_raises(ValueError, DiamondLattice, "invalid", generic_bond_length)
        np.testing.assert_raises(ValueError, DiamondLattice, -5, generic_bond_length)
        print("*** Unit Test passed ***")
        print("*** Unit test: check that any objects are other than DiamondLattice passed to initialize_lattice_builder raises a TypeError ***")
        np.testing.assert_raises(TypeError, pmb.initialize_lattice_builder, None)
        print("*** Unit test passed ***")
        # Check exceptions when the node and chain maps are incomplete
        incomplete_node_map = [{"particle_name": NodeType, "lattice_index": [0, 0, 0]},{"particle_name": NodeType, "lattice_index": [1, 1, 1]}]
        incomplete_chain_map = [{"node_start": "[0 0 0]", "node_end":"[1 1 1]" , "residue_list": residue_list}]
        np.testing.assert_raises(ValueError, pmb.define_hydrogel, "test_hydrogel", incomplete_node_map, chain_topology)
        np.testing.assert_raises(ValueError, pmb.define_hydrogel, "test_hydrogel", node_topology, incomplete_chain_map)
        # Check that two hydrogels with the same name can be defined in the dataframe
        pmb.define_hydrogel(hydrogel_name,node_topology, chain_topology)
        hydrogel_count = len(pmb.df[pmb.df["name"] == hydrogel_name])
        assert hydrogel_count == 2, f"Hydrogel '{hydrogel_name}' should be redefined."
        assert hydrogel_name in pmb.df["name"].values
        assert pmb.df.loc[pmb.df["name"] == hydrogel_name, "pmb_type"].values[0] == "hydrogel"

    def test_hydrogel_definitions_in_df(self):
        # Verify node_map and chain_map are correctly added
        compare_node_maps(pmb.df.loc[pmb.df["name"] == hydrogel_name, "node_map"].values[0], node_topology)
        compare_chain_maps(pmb.df.loc[pmb.df["name"] == hydrogel_name, "chain_map"].values[0], chain_topology)
        for chain_id in chain_topology:
            molecule_name = f"chain_{chain_id['node_start']}_{chain_id['node_end']}"
            assert molecule_name in pmb.df["name"].values
        #####-- Invalid hydrogel name --#####
        # Test if create_hydrogel raises an exception when provided with invalid data
        print("*** Unit Test: Check invalid inputs for create_hydrogel ***")
        with self.assertLogs(level='WARNING') as cm:
            pmb.create_hydrogel("invalid_hydrogel", espresso_system)
            self.assertEqual(cm.output, ["WARNING:root:Hydrogel with name 'invalid_hydrogel' is not defined in the DataFrame, no hydrogel will be created."])
        print("*** Invalid Input Test Passed ***")
        # Check if the molecules (chains) are correctly stored in the hydrogel data
        for ((molecule_id, molecule_data),_) in zip(hydrogel_info["chains"].items(),chain_topology):
            molecule_name_in_espresso = pmb.df[(pmb.df["pmb_type"] == "molecule") & (pmb.df["molecule_id"] == molecule_id)]["name"].values[0]
            np.testing.assert_equal(molecule_name_in_espresso, f"chain_{molecule_data['node_start']}_{molecule_data['node_end']}")

        print("*** Checking if the ends of an arbitrarly chosen chain is connected to node_start and node_end ***")

        molecule = hydrogel_info["chains"][1]
        Res_node_start = list(molecule.values())[0]
        Res_node_end   = list(molecule.values())[-3]
        central_bead_near_node_start = Res_node_start["central_bead_id"]
        central_bead_near_node_end = Res_node_end["central_bead_id"]

        node_ids = []
        for indice in node_labels.keys():
            index_pos = np.array(list(int(x) for x in indice.strip('[]').split()))*0.25*lattice_builder.box_l
            node_id = espresso_system.part.select(lambda p: (p.pos == index_pos).all()).id[0]
            node_ids.append(node_id)

        bead_ids_in_random_molecule = [i for i in range(central_bead_near_node_start, central_bead_near_node_end+1)]
        particle_ids = pmb.df["particle_id"].fillna(-1).to_numpy()
        particle_ids2 = pmb.df["particle_id2"].fillna(-1).to_numpy()

        mask = np.isin(particle_ids, node_ids) & np.isin(particle_ids2, bead_ids_in_random_molecule)
        filtered_df = pmb.df[mask]
        
        # Extract scalar values for central_bead_node_start and central_bead_node_end
        central_bead_node_start = filtered_df[filtered_df["particle_id2"] == central_bead_near_node_start]["particle_id"].iloc[0]
        central_bead_node_end = filtered_df[filtered_df["particle_id2"] == central_bead_near_node_end]["particle_id"].iloc[0]

        bond_name_node_start = filtered_df[
            (filtered_df["particle_id"] == central_bead_node_start) & 
            (filtered_df["particle_id2"] == central_bead_near_node_start)
        ]["name"].iloc[0]

        bond_name_node_end = filtered_df[
            (filtered_df["particle_id"] == central_bead_node_end) & 
            (filtered_df["particle_id2"] == central_bead_near_node_end)
        ]["name"].iloc[0]
        
        all_not_na = filtered_df['bond_object'].notna().all()
        
        assert all_not_na, "Bond object is not defined near nodes"
        
        central_bead_name_near_node_start = pmb.df[pmb.df["particle_id"]==central_bead_near_node_start]["name"].values[0]
        central_bead_name_near_node_end = pmb.df[pmb.df["particle_id"]==central_bead_near_node_end]["name"].values[0]

        if central_bead_name_near_node_start == BeadType1:
            possible_bond_names = [NodeType+"-"+BeadType1, BeadType1+"-"+NodeType]
            assert bond_name_node_start in possible_bond_names

        if central_bead_name_near_node_end == BeadType2:
            possible_bond_names = [NodeType+"-"+BeadType2, BeadType2+"-"+NodeType]
            assert bond_name_node_end in possible_bond_names
        
        print("*** Unit Test passed ***")

if __name__ == "__main__":
    ut.main()
