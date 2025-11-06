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

# Import pyMBE and other libraries
import pyMBE
import numpy as np
import pyMBE.storage.df_management as df_management

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

def check_charge_number_map(input_parameters):
    """
    Checks if pyMBE stores in the pmb.df the input parameters for acid/base particles correctly.

    Args:
        input_parameters(`dict`): dictionary with the input parameters for define_particle.

    """
    pmb.define_particle(**input_parameters)

    if input_parameters["acidity"] == "inert":
        np.testing.assert_equal(actual=pmb.get_charge_number_map(),
                                desired={0: input_parameters["z"]},
                                verbose=True)
    elif input_parameters["acidity"] == "acidic":
        np.testing.assert_equal(actual=pmb.get_charge_number_map(),
                                desired={0: 0, 1: -1},
                                verbose=True)
    elif input_parameters["acidity"] == "basic":
        np.testing.assert_equal(actual=pmb.get_charge_number_map(),
                                desired={0: 1, 1: 0},
                                verbose=True)

print("*** get_charge_number_map unit tests ***")
print("*** Unit test: check that get_charge_number_map works correctly for inert particles***")
# Clean pmb.df
pmb.df = df_management._DFManagement._setup_df()
input_parameters={"name":"I", 
                  "acidity": "inert",
                  "pka": np.nan,
                  "z":5}

check_charge_number_map(input_parameters)

print("*** Unit test passed ***")
print("*** Unit test: check that get_charge_number_map works correctly for acidic particles***")
# Clean pmb.df
pmb.df = df_management._DFManagement._setup_df()
input_parameters={"name":"A", 
                  "acidity": "acidic",
                  "pka":4}

check_charge_number_map(input_parameters)

print("*** Unit test passed ***")
print("*** Unit test: check that get_charge_number_map works correctly for basic particles***")
# Clean pmb.df
pmb.df = df_management._DFManagement._setup_df()
input_parameters={"name":"B", 
                  "acidity": "basic",
                  "pka":4}

check_charge_number_map(input_parameters)

print("*** Unit test passed ***")
print("*** All unit tests passed ***")
