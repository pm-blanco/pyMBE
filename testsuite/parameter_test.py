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

import pathlib
import pyMBE
import pandas as pd
import numpy as np
import pyMBE.storage.df_management as df_management

pmb = pyMBE.pymbe_library(seed=42)

print("*** Unit test: check that the different pKa sets are correctly formatted ***")

data_root = pathlib.Path(__file__).parent / "test_parameters"
params_root = pathlib.Path(pyMBE.__file__).parent / "parameters"
pka_root = params_root / "pka_sets"
peptides_root = params_root / "peptides"

for path in pka_root.glob("*.json"):
    print(f"Checking {path.stem}")
    pmb.load_pka_set(path)

print("*** Test passed ***")

print("*** Unit test: check that the order to execute load_pka_set() and load_interaction_parameters does not change the resulting parameters in pmb.df ***")
path_to_interactions=pmb.root / "parameters" / "peptides" / "Lunkad2021.json"
path_to_pka=pmb.root / "parameters" / "pka_sets" / "Hass2015.json"

# First order of loading parameters
pmb.df = df_management._DFManagement._setup_df() # clear the pmb_df
pmb.load_interaction_parameters (filename=peptides_root / "Lunkad2021.json")
pmb.load_pka_set(filename=pka_root / "Hass2015.json")
df_1 = pmb.df.copy()
df_1 = df_1.sort_values(by="name").reset_index(drop=True)
# Drop espresso types (they depend on the order of loading)
df_1 = df_1.drop(labels=('state_one', 'es_type'), axis=1).drop(labels=('state_two', 'es_type'), axis=1)
# Drop bond_object  (assert_frame_equal does not process it well)
df_1 = df_1.sort_index(axis=1).drop(labels="bond_object", axis=1)
# Second order of loading parameters
pmb.df = df_management._DFManagement._setup_df() # clear the pmb_df
pmb.load_pka_set (filename=path_to_pka)
#print(pmb.df["acidity"])
pmb.load_interaction_parameters(filename=path_to_interactions) 
#print(pmb.df["acidity"])
df_2 = pmb.df.copy()
df_2 = df_2.sort_values(by="name").reset_index(drop=True)
# Drop espresso types (they depend on the order of loading)
df_2 = df_2.drop(labels=('state_one', 'es_type'), axis=1).drop(labels=('state_two', 'es_type'), axis=1)
# Drop bond_object  (assert_frame_equal does not process it well)
df_2 = df_2.sort_index(axis=1).drop(labels="bond_object", axis=1)

df_1 = df_1.replace({pd.NA: np.nan})
df_2 = df_2.replace({pd.NA: np.nan})
pd.testing.assert_frame_equal(df_1,df_2)

print("*** Test passed ***")

print("*** Unit test: check that  load_interaction_parameters loads FENE bonds correctly ***")
pmb.df = df_management._DFManagement._setup_df() # clear the pmb_df
pmb.load_interaction_parameters (filename=data_root / "test_FENE.json")

expected_parameters = {'r_0' : 0.4*pmb.units.nm,
                        'k'  : 400 * pmb.units('reduced_energy / reduced_length**2'),
                        'd_r_max': 0.8 * pmb.units.nm}
reduced_units = {'r_0'    : 'reduced_length',
                     'k'      : 'reduced_energy / reduced_length**2',
                     'd_r_max': 'reduced_length'}
parameters_in_df = pmb.df[pmb.df.pmb_type == "bond"].parameters_of_the_potential.values[0]

for key in expected_parameters.keys():
    np.testing.assert_equal(actual=parameters_in_df[key],
                desired=expected_parameters[key].m_as(reduced_units[key]),
                verbose=True)

print("*** Test passed ***")
print("*** Unit test: check that  load_interaction_parameters loads residue, molecule and peptide objects correctly ***")

pmb.df = df_management._DFManagement._setup_df() # clear the pmb_df
pmb.load_interaction_parameters (filename=data_root / "test_molecules.json")

expected_residue_parameters={"central_bead":  "A", "side_chains": ["B","C"] }
expected_molecule_parameters={"residue_list":   ["R1","R1", "R1"]}
expected_peptide_parameters= {"sequence":   ['K', 'K', 'K', 'K', 'K', 'D', 'D', 'D', 'D', 'D'], "model": "1beadAA" }

# Check residue
np.testing.assert_equal(actual=pmb.df[pmb.df.name == "R1"].central_bead.values[0],
                desired=expected_residue_parameters["central_bead"],
                verbose=True)

np.testing.assert_equal(actual=frozenset(pmb.df[pmb.df.name == "R1"].side_chains.values[0]),
                desired=frozenset(expected_residue_parameters["side_chains"]),
                verbose=True)
# Check molecule
np.testing.assert_equal(actual=frozenset(pmb.df[pmb.df.name == "M1"].residue_list.values[0]),
                desired=frozenset(expected_molecule_parameters["residue_list"]),
                verbose=True)
# Check peptide
np.testing.assert_equal(actual=pmb.df[pmb.df.name == "P1"].sequence.values[0],
                desired=expected_peptide_parameters["sequence"],
                verbose=True)
np.testing.assert_equal(actual=frozenset(pmb.df[pmb.df.name == "P1"].model.values[0]),
                desired=frozenset(expected_peptide_parameters["model"]),
                verbose=True)
print("*** Test passed ***")
print("*** Unit test: check that  load_interaction_parameters raises a ValueError if one loads a data set with an unknown pmb_type ***")
pmb.df = df_management._DFManagement._setup_df() # clear the pmb_df
input_parameters={"filename": data_root / "test_non_valid_object.json"}
np.testing.assert_raises(ValueError, pmb.load_interaction_parameters, **input_parameters)
print("*** Test passed ***")
print("*** Unit test: check that  load_interaction_parameters raises a ValueError if one loads a bond not supported by pyMBE ***")
pmb.df = df_management._DFManagement._setup_df() # clear the pmb_df
input_parameters={"filename": data_root / "test_non_valid_bond.json"}
np.testing.assert_raises(ValueError, pmb.load_interaction_parameters, **input_parameters)
print("*** Test passed ***")
print("*** Unit test: check that  check_pka_set raises a ValueError if data is missing important fields ***")
np.testing.assert_raises(ValueError, pmb.check_pka_set, {"name" : {}})
np.testing.assert_raises(ValueError, pmb.check_pka_set, {"name" : {"pka_value": 1.}})
np.testing.assert_raises(ValueError, pmb.check_pka_set, {"name" : {"acidity": 1.}})
print("*** Test passed ***")
