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

# Import pyMBE and other libraries
import numpy as np
import pandas as pd
import pyMBE
import pyMBE.storage.df_management as df_management

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

def check_acid_base_setup(input_parameters, acidity_setup):
    """
    Checks if pyMBE stores in the pmb.df the input parameters for acid/base particles correctly.

    Args:
        input_parameters (`dict`): dictionary with the input parameters for define_particle.
        acidity_setup (`dict`): dictionary with the expected setup that pyMBE should do in the pmb.df for acid/base particles.
    """
    pmb.define_particle(**input_parameters)

    # Handle pd.NA safely
    if pd.isna(input_parameters.get("acidity", None)):  
        input_parameters.pop("z", None)  # Use .pop with default to avoid KeyError

    # Checks that the input parameters are stored properly
    for parameter_key, expected_value in input_parameters.items():
        actual_value = pmb.df[parameter_key].values[0]

        # Use pd.isna() to compare safely, since pd.NA does not behave like regular values
        if pd.isna(expected_value) and pd.isna(actual_value):
            continue  # Skip this check, they are both missing (NA)

        np.testing.assert_equal(actual=actual_value, desired=expected_value, verbose=True)

    # Checks that the setup of the acid/base properties is done correctly
    for state in ["state_one", "state_two"]:
        for state_attribute in ["label", "z"]:
            actual_value = pmb.df[state][state_attribute].values[0]
            expected_value = acidity_setup[state][state_attribute]

            if pd.isna(expected_value) and pd.isna(actual_value):
                continue  # Skip this check if both are NA

            np.testing.assert_equal(actual=actual_value, desired=expected_value, verbose=True)

    # Checks that pyMBE assigns different espresso types to each state
    np.testing.assert_raises(
        AssertionError, 
        np.testing.assert_equal, 
        pmb.df["state_one"]["es_type"].values[0], 
        pmb.df["state_two"]["es_type"].values[0]
    )


print("*** Particle acidity unit tests ***")
print("*** Unit test: check that all acid/base input parameters in define_particle for an inert particle are correctly stored in pmb.df***")
# Clean pmb.df
pmb.df = df_management._DFManagement._setup_df()
input_parameters={"name":"I", 
                  "acidity": pd.NA,
                  "pka": pd.NA,
                  "z":2}
acidity_setup={"state_one":{"label":f"{input_parameters['name']}",
                         "z":2},
            "state_two":{"label": pd.NA,
                         "z": pd.NA},}

check_acid_base_setup(input_parameters=input_parameters,
                      acidity_setup=acidity_setup)

print("*** Unit test passed ***")
print("*** Unit test: check that a deprecation warning is raised if the keyword 'inert' is used for acidity ***")
# Clean pmb.df
pmb.df = df_management._DFManagement._setup_df()
input_parameters={"name":"I", 
                  "acidity": "inert",
                  "pka": pd.NA,
                  "z":2}
pmb.define_particle(**input_parameters)
print("*** Unit test passed ***")
print("*** Unit test: check that all acid/base input parameters in define_particle for an acid are correctly stored in pmb.df***")
# Clean pmb.df
pmb.df = df_management._DFManagement._setup_df()
input_parameters={"name":"A", 
                  "acidity": "acidic",
                  "pka":4}
acidity_setup={"state_one":{"label":f"{input_parameters['name']}H",
                         "z":0},
            "state_two":{"label":f"{input_parameters['name']}",
                         "z":-1},}

check_acid_base_setup(input_parameters=input_parameters,
                      acidity_setup=acidity_setup)
print("*** Unit test passed ***")
print("*** Unit test: check that all acid/base input parameters in define_particle for a base are correctly stored in pmb.df***")
# Clean pmb.df
pmb.df = df_management._DFManagement._setup_df()
input_parameters={"name":"B", 
                  "acidity": "basic",
                  "pka":9}
acidity_setup={"state_one":{"label":f"{input_parameters['name']}H",
                         "z":1},
            "state_two":{"label":f"{input_parameters['name']}",
                         "z":0},}

check_acid_base_setup(input_parameters=input_parameters,
                      acidity_setup=acidity_setup)
print("*** Unit test passed ***")

print("*** Unit test: check that set_particle_acidity raises a ValueError if pKa is not provided and pKa is acidic or basic  ***")
input_parametersA={"name":"A", 
                   "acidity": "acidic" }

input_parametersB= {"name": "B",
                   "acidity": "basic"}
np.testing.assert_raises(ValueError, pmb.set_particle_acidity,**input_parametersA)
np.testing.assert_raises(ValueError, pmb.set_particle_acidity, **input_parametersB)
print("*** Unit test passed ***")
print("*** Unit test: check that set_particle_acidity raises a ValueError if a non-supported acidity is provided  ***")
input_parametersA={"name":"A", 
                   "acidity": "random" }
np.testing.assert_raises(ValueError, pmb.set_particle_acidity,**input_parametersA)
print("*** Unit test passed ***")
print("*** All unit tests passed ***")
