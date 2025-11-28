#
# Copyright (C) 2024 pyMBE-dev team
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
# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(seed=42)

print("*** Create_molecule with input position list unit test ***")
print("*** Unit test: Check that the positions of the central bead of the first residue in the generated molecules are equal to the input positions ***")
# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)
solvent_permitivity = 78.3
N_molecules = 3
chain_length = 5
molecule_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L

pos_list = [[10,10,10], [20,20,20], [30,30,30]]
pmb.define_particle(name='central_mon',
                        sigma=0.35*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name='side_mon',
                        sigma=0.35*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))

pmb.define_residue(
    name = 'res1',
    central_bead = 'central_mon',
    side_chains = ['side_mon', 'side_mon']
    )

bond_type = 'harmonic'
generic_bond_length=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

harmonic_bond = {'r_0'    : generic_bond_length,
                 'k'      : generic_harmonic_constant,
                 }

pmb.define_default_bond(bond_type = bond_type, bond_parameters = harmonic_bond)

# Defines the peptide in the pyMBE data frame
molecule_name = 'generic_molecule'
pmb.define_molecule(name=molecule_name, residue_list = ['res1']*chain_length)

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name, z=1, sigma=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  z=-1, sigma=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

# System parameters
volume = N_molecules/(pmb.N_A*molecule_concentration)
L = volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration = N_molecules/(volume*pmb.N_A)

# Create an instance of an espresso system
espresso_system=espressomd.System(box_l = [L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

pmb.write_pmb_df("df1.csv")
# Create your molecules into the espresso system
molecules = pmb.create_molecule(name=molecule_name,
                        number_of_molecules= N_molecules,
                        espresso_system=espresso_system,
                        use_default_bond=True,
                        list_of_first_residue_positions = pos_list)
pmb.write_pmb_df("df2.csv")
# Running unit test here. Use np.testing.assert_almost_equal of the input position list and the central_bead_pos list under here.
central_bead_pos = []
for molecule_id in molecules:
    info = next(iter(molecules[molecule_id].values()))
    central_bead_id = info['central_bead_id']
    side_chain_ids = info['side_chain_ids']
    central_bead_pos.append(espresso_system.part.by_id(central_bead_id).pos.tolist())

np.testing.assert_almost_equal(pos_list, central_bead_pos)

print("*** Unit test passed ***\n")

print("*** Unit test: Check that create_molecule raises a ValueError if the user does not provide a nested list for list_of_first_residue_positions***")
input_parameters={"name": "generic_molecule",
                 "number_of_molecules": 1,
                 "espresso_system": espresso_system,
                 "list_of_first_residue_positions": [1,2,3]}
np.testing.assert_raises(ValueError, pmb.create_molecule, **input_parameters)
print("*** Unit test passed ***\n")

print("*** Unit test: Check that create_molecule raises a ValueError if the user does not provide a nested list with three coordinates***")
input_parameters={"name": "generic_molecule",
                 "number_of_molecules": 1,
                 "espresso_system": espresso_system,
                 "list_of_first_residue_positions": [[1,2]]}
np.testing.assert_raises(ValueError, pmb.create_molecule, **input_parameters)
print("*** Unit test passed ***\n")

print("*** Unit test: Check that create_molecule raises a ValueError if the user does not provide a the same number of first_residue_positions as number_of_molecules***")
input_parameters={"name": "generic_molecule",
                 "number_of_molecules": 2,
                 "espresso_system": espresso_system,
                 "list_of_first_residue_positions": [[1,2,3]]}
np.testing.assert_raises(ValueError, pmb.create_molecule, **input_parameters)
print("*** Unit test passed ***\n")

print("*** Unit test: Check that center_molecule_in_simulation_box works correctly for cubic boxes***")

molecule_id = pmb.df.loc[pmb.df['name']==molecule_name].molecule_id.values[0]
pmb.center_molecule_in_simulation_box(molecule_id=molecule_id, 
                                      espresso_system=espresso_system)
center_of_mass = pmb.calculate_center_of_mass_of_molecule(molecule_id=molecule_id, 
                                                          espresso_system=espresso_system)
center_of_mass_ref = [L.to('reduced_length').magnitude/2]*3

np.testing.assert_almost_equal(center_of_mass, center_of_mass_ref)

print("*** Unit test passed ***\n")

print("*** Unit test: Check that center_molecule_in_simulation_box works correctly for non-cubic boxes***")

espresso_system.change_volume_and_rescale_particles(d_new=3*L.to('reduced_length').magnitude, dir="z")
molecule_id = pmb.df.loc[pmb.df['name']==molecule_name].molecule_id.values[2]
pmb.center_molecule_in_simulation_box(molecule_id=molecule_id, espresso_system=espresso_system)
center_of_mass = pmb.calculate_center_of_mass_of_molecule(molecule_id=molecule_id, espresso_system=espresso_system)
center_of_mass_ref = [L.to('reduced_length').magnitude/2, L.to('reduced_length').magnitude/2, 1.5*L.to('reduced_length').magnitude]

np.testing.assert_almost_equal(center_of_mass, center_of_mass_ref)

print("*** Unit test passed ***")

print("*** Unit test: Check that center_molecule_in_simulation_box raises a Value Error if a wrong molecule_id is provided***")

input_parameters = {"molecule_id": 20 ,
                    "espresso_system":espresso_system}

np.testing.assert_raises(ValueError, pmb.center_molecule_in_simulation_box, **input_parameters)

print("*** Unit test passed ***")
