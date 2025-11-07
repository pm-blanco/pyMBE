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

import tempfile
import espressomd
import pandas as pd
import numpy as np
import logging
import io
import pyMBE.storage.df_management as df_management

# Create an in-memory log stream
log_stream = io.StringIO()
logging.basicConfig(level=logging.INFO, 
                    format="%(levelname)s: %(message)s",
                    handlers=[logging.StreamHandler(log_stream)])
# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(seed=42)

print ('*** Unit tests: read and write from pyMBE dataframe ***')

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)

# Define particles
pmb.define_particle(
    name = "I",
    sigma = 0.3*pmb.units.nm,
    epsilon = 1*pmb.units('reduced_energy'),
    z = 0)

pmb.define_particle(
    name = "A",
    acidity = "acidic",
    pka = 4,
    sigma = 0.3*pmb.units.nm,
    epsilon = 1*pmb.units('reduced_energy'),)
    
pmb.define_particle(
    name = "B",
    acidity = "basic",
    pka = 9,
    sigma = 0.3*pmb.units.nm,
    epsilon = 1*pmb.units('reduced_energy'),)

# Define residues
pmb.define_residue(
    name = "Res_1",
    central_bead = "I",
    side_chains = ["A","B"])

pmb.define_residue(
    name = "Res_2",
    central_bead = "I",
    side_chains = ["Res_1"])   

# Define peptide
peptide_name = 'generic_peptide'
peptide_sequence = 'EEEEEEE'
peptide_model = '2beadAA'
pmb.define_peptide(name=peptide_name, sequence=peptide_sequence, model=peptide_model)

# Define a molecule
molecule_name = "A_molecule"
n_molecules = 1

pmb.define_molecule(
    name = molecule_name,
    residue_list = ["Res_1", "Res_1",
                    "Res_2", "Res_1",
                    "Res_1", "Res_2",
                    "Res_2"])

# Define a bond
bond_type = 'harmonic'
generic_bond_length=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

harmonic_bond = {'r_0'    : generic_bond_length,
                 'k'      : generic_harmonic_constant,
                 }


pmb.define_default_bond(bond_type = bond_type, bond_parameters = harmonic_bond)
pmb.define_bond(bond_type = bond_type, 
                bond_parameters = harmonic_bond,
                particle_pairs=[["A","A"],["B","B"]])
bond_type = 'FENE'
FENE_bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                'd_r_max': 0.8 * pmb.units.nm}

pmb.define_bond(bond_type = bond_type, 
                bond_parameters = FENE_bond,
                particle_pairs=[["A","B"]])

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name, z=1, sigma=0.35*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  z=-1, sigma=0.35*pmb.units.nm,  epsilon=1*pmb.units('reduced_energy'))

# System parameters
molecule_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L
volume = n_molecules/(pmb.N_A*molecule_concentration)
L = volume ** (1./3.) # Side of the simulation box
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)
pmb.add_bonds_to_espresso(espresso_system = espresso_system)

# Setup potential energy

pmb.setup_lj_interactions (espresso_system=espresso_system)
pd.options.display.max_colwidth = 10

# Copy the pmb.df into a new DF for the unit test 
stored_df = pmb.df.copy()

with tempfile.TemporaryDirectory() as tmp_directory:
    # Write the pymbe DF to a csv file
    df_filename = f'{tmp_directory}/df-example_molecule.csv'
    pmb.write_pmb_df (filename = df_filename)
    # Read the same pyMBE df from a csv a load it in pyMBE
    read_df = pmb.read_pmb_df(filename = df_filename)
    # Write the pyMBE DF to a txt file
    df_filename_test = f"{tmp_directory}/df-example_molecule.txt"
    pmb.write_pmb_df (filename = df_filename_test)
    np.testing.assert_raises(ValueError, pmb.read_pmb_df, df_filename_test)

stored_df['node_map'] = stored_df['node_map'].astype(object)
stored_df['chain_map'] = stored_df['chain_map'].astype(object)
stored_df['l0'] = stored_df['l0'].astype(object)

read_df['node_map'] = read_df['node_map'].astype(object)
read_df['chain_map'] = read_df['chain_map'].astype(object)
read_df['l0'] = read_df['l0'].astype(object)

# Preprocess data for the Unit Test
# The espresso bond object must be converted to a dict in order to compare them using assert_frame_equal
stored_df['bond_object']  = stored_df['bond_object'].apply(lambda x: (x.name(), x.get_params(), x._bond_id) if pd.notnull(x) else x)
read_df['bond_object']  = read_df['bond_object'].apply(lambda x: (x.name(), x.get_params(), x._bond_id) if pd.notnull(x) else x)
print("*** Unit test: check that the dataframe stored by pyMBE to file is the same as the one read from the file (same values and variable types) ***")

# One needs to replace the pd.NA by np.nan otherwise the comparison between pint objects fails
stored_df = stored_df.replace({pd.NA: np.nan})
read_df = read_df.replace({pd.NA: np.nan})

pd.testing.assert_frame_equal(stored_df, 
                                read_df,
                                rtol=1e-5)
print("*** Unit test passed***")

# Test that copy_df_entry raises an error if one provides a non-valid column name
print("*** Unit test: check that copy_df_entry raises an error if the entry does not exist ***")
np.testing.assert_raises(ValueError, df_management._DFManagement._copy_df_entry, df = pmb.df, name='test', column_name='non_existing_column',number_of_copies=1)
print("*** Unit test passed***")