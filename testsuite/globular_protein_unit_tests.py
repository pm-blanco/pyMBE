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
import espressomd
import pyMBE
import re
import json
import pathlib
from pint import UnitRegistry

ureg = UnitRegistry()

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

c_salt    =  0.01  * pmb.units.mol / pmb.units.L  
c_protein =  2e-4 * pmb.units.mol / pmb.units.L 
Box_V =  1. / (pmb.N_A*c_protein)
Box_L = Box_V**(1./3.) 

def custom_deserializer(dct):
    if "value" in dct and "unit" in dct:
        return ureg.Quantity(dct["value"], dct["unit"])  
    return dct  

protein_pdb = '1f6s'

print("*** Unit test: check that read_protein_vtf_in_df() loads the protein topology correctly ***")

path_to_parfile = pathlib.Path(__file__).parent / "tests_data" / "protein_topology_dict.json"
path_to_cg=pmb.root / "parameters" / "globular_proteins" / f"{protein_pdb}.vtf"
topology_dict = pmb.read_protein_vtf_in_df (filename=path_to_cg)


with open (path_to_parfile, "r") as file:
    load_json = json.load(file,object_hook=custom_deserializer)

np.testing.assert_equal(actual= topology_dict, 
                        desired= load_json,
                        verbose = True)

print("*** Unit test passed ***")


print("*** Unit test: check that define_protein() defines the aminoacids in the protein correctly ***")

protein_model = '2beadAA'

pmb.define_protein (name=protein_pdb, 
                    topology_dict=topology_dict, 
                    model = protein_model,
                    lj_setup_mode = "wca")
sequence = []
clean_sequence= []

for aminoacid in topology_dict.keys():
    
    input_parameters=topology_dict[aminoacid]
    residue_name = re.split(r'\d+', aminoacid)[0]
    sequence.append(residue_name)
    
    if residue_name not in ['CA', 'n', 'c','Ca']:
        clean_sequence.append(residue_name)

residue_list = pmb.define_AA_residues(sequence= clean_sequence,
                                      model = protein_model)

for residue in residue_list:
    for index in pmb.df[pmb.df['name']==residue].index:
        np.testing.assert_equal(actual=str(pmb.df.loc[index, "pmb_type"].values[0]), 
                        desired="residue", 
                        verbose=True)

protein_index = pmb.df[pmb.df['name']==protein_pdb].index

np.testing.assert_equal(actual=str(pmb.df.loc[protein_index, "name"].values[0]), 
                            desired=protein_pdb, 
                            verbose=True)        

np.testing.assert_equal(actual=pmb.df.loc[protein_index, ('sequence','')].values[0], 
                    desired=clean_sequence, 
                    verbose=True)


np.testing.assert_equal(actual=pmb.df.loc[protein_index, ('residue_list','')].values[0], 
                    desired=residue_list, 
                    verbose=True)  

print("*** Unit test passed ***")

print("*** Unit test: check that define_protein() raises a ValueError if the user provides a wrong model")

input_parameters={"name": protein_pdb,
                 "topology_dict": topology_dict,
                 "model" : "3beadAA",
                "lj_setup_mode": "wca"}

np.testing.assert_raises(ValueError, pmb.define_protein, **input_parameters)

input_parameters={"name": protein_pdb,
                 "topology_dict": topology_dict,
                 "model" : protein_model,
                "lj_setup_mode": "awc"}

np.testing.assert_raises(ValueError, pmb.define_protein, **input_parameters)

print("*** Unit test passed ***")


print("*** Unit test: check that create_protein() creates all the particles in the protein into the espresso_system with the properties defined in pmb.df  ***")

espresso_system=espressomd.System(box_l = [Box_L.to('reduced_length').magnitude] * 3)

# Here we upload the pka set from the reference_parameters folder
path_to_pka=pmb.root / "parameters" / "pka_sets" / "Nozaki1967.json"
pmb.load_pka_set(filename=path_to_pka)

pmb.create_protein(name=protein_pdb,
                    number_of_proteins=1,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)

residue_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].residue_id.dropna().to_list()

particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()

molecule_id = pmb.df.loc[pmb.df['name']==protein_pdb].molecule_id.values[0]

center_of_mass_es = pmb.calculate_center_of_mass_of_molecule ( molecule_id=molecule_id,espresso_system=espresso_system)

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
    charge = espresso_system.part.by_id(id).q
    es_type = espresso_system.part.by_id(id).type

    residue_id = pmb.df.loc[pmb.df['particle_id']==id].residue_id.values[0]
    residue_name = pmb.df.loc[pmb.df['particle_id']==id].name.values[0]
    
    initial_pos = topology_dict[f"{residue_name}{residue_id}"]['initial_pos']
    index = pmb.df.loc[pmb.df['particle_id']==id].index

    for axis in axis_list:
        distance_es[axis] = (initial_pos_es[axis] - center_of_mass_es[axis])**2
        distance_topology[axis] = (initial_pos[axis] - center_of_mass[axis])**2

    relative_distance_es = np.sqrt(np.sum(distance_es))
    relative_distance = np.sqrt(np.sum(distance_es))
    
    np.testing.assert_equal(actual=relative_distance_es, 
                        desired=relative_distance, 
                        verbose=True)

    np.testing.assert_equal(actual=charge, 
                        desired=pmb.df.loc[index, ("state_one","z")].values[0], 
                        verbose=True)
      
    np.testing.assert_equal(actual=es_type, 
                    desired=pmb.df.loc[index, ("state_one","es_type")].values[0], 
                    verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that create_protein() does not create any protein for number_of_proteins <= 0  ***")

starting_number_of_particles=len(espresso_system.part.all())

pmb.create_protein(name=protein_pdb,
                    number_of_proteins=0,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)

pmb.create_protein(name=protein_pdb,
                    number_of_proteins=-1,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)

np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles, 
                        verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that enable_motion_of_rigid_object does the setup correctly ***")

positions = []
for pid in particle_id_list:
    positions.append(espresso_system.part.by_id(pid).pos)

pmb.enable_motion_of_rigid_object(espresso_system=espresso_system,
                                  name=protein_pdb)

momI = 0
molecule_id = pmb.df.loc[pmb.df['name']==protein_pdb].molecule_id.values[0]
for p in espresso_system.part:
    center_of_mass = pmb.calculate_center_of_mass_of_molecule ( molecule_id=molecule_id,espresso_system=espresso_system)
    if p.mass > 1: 
        rigid_object_id = p.id 
        rigid_object_mass = espresso_system.part.by_id(rigid_object_id).mass
        rigid_object_rotation = espresso_system.part.by_id(rigid_object_id).rotation
        rigid_object_intertia  = np.copy(espresso_system.part.by_id(rigid_object_id).rinertia)

        np.testing.assert_equal(actual=rigid_object_mass, 
                    desired=len(particle_id_list), 
                    verbose=True)
        print ('mass passed ')

        np.testing.assert_equal(actual=rigid_object_rotation, 
        desired=[1, 1, 1], 
        verbose=True)
        print ('rotation  passed ')

        for pid in particle_id_list:
            momI += np.power(np.linalg.norm(center_of_mass - espresso_system.part.by_id(pid).pos), 2)
        rinertia = np.ones(3) * momI

        np.testing.assert_array_almost_equal(rinertia, rigid_object_intertia)

print("*** Unit test passed ***")

print("*** Unit test: check that enable_motion_of_rigid_object() raises a ValueError if a wrong pmb_type is provided***")

input_parameters = {"espresso_system":espresso_system,
                    "name": "CA"}

np.testing.assert_raises(ValueError, pmb.enable_motion_of_rigid_object, **input_parameters)

print("*** Unit test passed ***")

print("*** Unit test: check that protein_sequence_parser() correctly returns que protein sequence ***")

def test_sequence(input,output):
    """
    Tests that the pyMBE parses correctly the input sequence.

    Args:
        input(`str` or `lst` of `str`): input protein sequence.
        ouput(`lst` of `str`): expected ouput protein sequence.
    """
    clean_sequence= pmb.protein_sequence_parser(sequence = input)
    np.testing.assert_equal(actual=clean_sequence, 
                            desired=output, 
                            verbose=True)

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

print("*** Unit test: check that protein_sequence_parser() raises a ValueError if a wrong residue key is provided***")

input_parameters = {"sequence":"rekx"}
np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

input_parameters = {"sequence":"ARG-GLU-TUR-HIS"}
np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

input_parameters = {"sequence":"A-E-E-X"}
np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

input_parameters = {"sequence":"a-e-e-x"}
np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

input_parameters = {"sequence":["A", "E","X"]}
np.testing.assert_raises(ValueError, pmb.protein_sequence_parser, **input_parameters)

print("*** Unit test passed ***")

print("*** Unit test: Check  that check_aminoacid_key returns True for any latter valid in the one letter amino acid code***")
valid_AA_keys=['V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C']
for key in valid_AA_keys:
    np.testing.assert_equal(actual=pmb.check_aminoacid_key(key=key), 
                        desired=True, 
                        verbose=True)
print("*** Unit test passed ***\n")
print("*** Unit test: Check  that check_aminoacid_key returns False for a key not valid in the one letter amino acid code ***")
np.testing.assert_equal(actual=pmb.check_aminoacid_key(key="B"), 
                        desired=False, 
                        verbose=True)
print("*** Unit test passed ***\n")

print("*** Unit test: Check  that check_if_metal_ion returns True for any key corresponding to a supported metal ion ***")
for key in pmb.get_metal_ions_charge_number_map().keys():
    np.testing.assert_equal(actual=pmb.check_if_metal_ion(key=key), 
                        desired=True, 
                        verbose=True)
print("*** Unit test passed ***\n")
print("*** Unit test: Check  that check_if_metal_ion returns False for a key not corresponding to a supported metal ion ***")
np.testing.assert_equal(actual=pmb.check_if_metal_ion(key="B"), 
                        desired=False, 
                        verbose=True)
print("*** Unit test passed ***\n")

print("*** Unit test: Check  that get_metal_ions_charge_number_map returns the correct charge map for metals ***")
metal_charge_number_map = {"Ca": 2}
pmb_metal_charge_number_map = pmb.get_metal_ions_charge_number_map()

np.testing.assert_equal(actual=pmb_metal_charge_number_map, 
                    desired=metal_charge_number_map, 
                    verbose=True)
print("*** Unit test passed ***\n")


print("*** Unit test: check that define_AA_residues()")

test_sequence = ['c','n', 'G','V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C' ]

valid_protein_model = ['1beadAA', '2beadAA']

output =['AA-c', 'AA-n', 'AA-G', 'AA-V', 'AA-I', 'AA-L', 'AA-E', 'AA-Q', 'AA-D', 'AA-N', 'AA-H', 'AA-W', 'AA-F', 'AA-Y', 'AA-R', 'AA-K', 'AA-S', 'AA-T', 'AA-M', 'AA-A', 'AA-G', 'AA-P', 'AA-C']

for protein_model in valid_protein_model:

    pmb_residue_list = pmb.define_AA_residues(sequence= test_sequence,
                                        model = protein_model)

    np.testing.assert_equal(actual=pmb_residue_list, 
                            desired=output, 
                            verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: check that define_peptide() raises a ValueError if a wrong model key is provided")

input_parameters = {"name": "generic_peptide",
                    "sequence": "EEEEEEE",
                    "model": "3beadAA" }

np.testing.assert_raises(ValueError, pmb.define_peptide, **input_parameters)

input_parameters = {"name": "generic_peptide",
                    "sequence": "EEEEEEE",
                    "model": "beadAA" }

np.testing.assert_raises(ValueError, pmb.define_peptide, **input_parameters)

print("*** Unit test passed ***")

print("*** Unit test: check that search_particles_in_residue() returns the correct list of residues")

list_of_residues = ['AA-c', 'AA-n', 'AA-G', 'AA-V', 'AA-I', 'AA-L', 'AA-E', 'AA-Q', 'AA-D', 'AA-N', 'AA-H', 
                    'AA-W', 'AA-F', 'AA-Y', 'AA-R', 'AA-K', 'AA-S', 'AA-T', 'AA-M', 'AA-A','AA-P', 'AA-C']

for residue_name in list_of_residues:

    residue = residue_name.replace('AA-','')

    list_of_particles_in_residue= pmb.search_particles_in_residue(residue_name = residue_name)
    if residue in ['c', 'n']:        
        np.testing.assert_equal(actual=list_of_particles_in_residue, 
                                desired=[residue], 
                                verbose=True)
    elif residue == "G":
        np.testing.assert_equal(actual=list_of_particles_in_residue, 
                                desired=[], 
                                verbose=True)
    else:
        np.testing.assert_equal(actual=list_of_particles_in_residue, 
                                desired=['CA', residue], 
                                verbose=True)
        
                               
print("*** Unit test passed ***")

print("*** Unit test: check that search_particles_in_residue() returns the correct list of residues for nested residues case")

pmb.define_particle(
    name = "I",
    sigma = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))
    
# Acidic particle
pmb.define_particle(
    name = "A",
    sigma = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))
    
# Basic particle
pmb.define_particle(
    name = "B",
    sigma = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))

pmb.define_residue(
    name = "Res_1",
    central_bead = "I",
    side_chains = ["A","B"])
    
pmb.define_residue(
    name = "Res_2",
    central_bead = "I",
    side_chains = ["Res_1"])

list_of_particles_in_residue= pmb.search_particles_in_residue(residue_name = "Res_2")

np.testing.assert_equal(actual=list_of_particles_in_residue, 
                                desired=['I', 'I', 'A', 'B'], 
                                verbose=True)

print("*** Unit test passed ***")

print("*** Unit test: Check that create_protein() does not create any protein for an undefined protein name ***")
starting_number_of_particles=len(espresso_system.part.all())
pmb.create_protein(name="undefined_protein",
                    number_of_proteins=1,
                    espresso_system=espresso_system,
                    topology_dict=topology_dict)
np.testing.assert_equal(actual=len(espresso_system.part.all()), 
                        desired=starting_number_of_particles, 
                        verbose=True)
print("*** Unit test passed ***")
