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
import pandas as pd
from scipy import interpolate

def determine_reservoir_concentrations_test_ideal(pH_res, c_salt_res):

    # Create an instance of the pyMBE library
    pmb = pyMBE.pymbe_library(seed=42)

    # Determine the reservoir composition using pyMBE
    cH_res_pyMBE, cOH_res_pyMBE, cNa_res_pyMBE, cCl_res_pyMBE = pmb.determine_reservoir_concentrations(
            pH_res,
            c_salt_res * pmb.units.mol/ pmb.units.L,
            lambda x: 1.0)

    # Determine the reservoir concentrations without pyMBE
    cH_res = 10**(-pH_res) * pmb.units.mol/ pmb.units.L
    cOH_res = 10**(pH_res-14) * pmb.units.mol/ pmb.units.L
    c_salt_res = c_salt_res * pmb.units.mol/ pmb.units.L
    cNa_res = max(c_salt_res, c_salt_res + cOH_res - cH_res) 
    cCl_res = max(c_salt_res, c_salt_res + cH_res - cOH_res) 

    # Check the determined concentrations 
    np.testing.assert_allclose(cH_res, cH_res_pyMBE)
    np.testing.assert_allclose(cOH_res, cOH_res_pyMBE)
    np.testing.assert_allclose(cNa_res, cNa_res_pyMBE)
    np.testing.assert_allclose(cCl_res, cCl_res_pyMBE)

def determine_reservoir_concentrations_test_interacting(pH_res, c_salt_res, reference_data):

    # Create an instance of the pyMBE library
    pmb = pyMBE.pymbe_library(seed=42)

    # Load the excess chemical potential data
    monovalent_salt_ref_data=pd.read_csv(pmb.root / "parameters" / "salt" / "excess_chemical_potential_excess_pressure.csv")
    ionic_strength = pmb.units.Quantity(monovalent_salt_ref_data["cs_bulk_[1/sigma^3]"].values, "1/reduced_length**3")
    excess_chemical_potential = pmb.units.Quantity(monovalent_salt_ref_data["excess_chemical_potential_[kbT]"].values, "reduced_energy")
    excess_chemical_potential_interpolated = interpolate.interp1d(ionic_strength.m_as("1/reduced_length**3"), 
                                                                                    excess_chemical_potential.m_as("reduced_energy"))
    activity_coefficient_monovalent_pair = lambda x: np.exp(excess_chemical_potential_interpolated(x.to('1/(reduced_length**3 * N_A)').magnitude))

    # Determine the reservoir composition using pyMBE
    cH_res_pyMBE, cOH_res_pyMBE, cNa_res_pyMBE, cCl_res_pyMBE = pmb.determine_reservoir_concentrations(
            pH_res,
            c_salt_res * pmb.units.mol/ pmb.units.L,
            activity_coefficient_monovalent_pair)

    # Check the results
    row = reference_data[(reference_data['pH_res'] == pH_res) & (reference_data['c_salt_res'] == c_salt_res)]
    np.testing.assert_allclose(cH_res_pyMBE.magnitude, row["cH_res"])
    np.testing.assert_allclose(cOH_res_pyMBE.magnitude, row["cOH_res"])
    np.testing.assert_allclose(cNa_res_pyMBE.magnitude, row["cNa_res"])
    np.testing.assert_allclose(cCl_res_pyMBE.magnitude, row["cCl_res"])


print("*** Unit test: check that pyMBE correctly calculates the reservoir composition when setting up the grand-reaction method (ideal case). ***")
for pH_res in np.linspace(1.0, 13.0, 5):
    for c_salt_res in np.logspace(-5, 0, 5):
        determine_reservoir_concentrations_test_ideal(pH_res, c_salt_res)
print("*** Unit test passed ***")

print("*** Unit test: check that pyMBE correctly calculates the reservoir composition when setting up the grand-reaction method (interacting case). ***")
# Load reference data
reference_data = pd.read_csv('./determine_reservoir_concentrations_test_data/data.csv')
determine_reservoir_concentrations_test_interacting(10.0, 0.1, reference_data)
print("*** Unit test passed ***")
