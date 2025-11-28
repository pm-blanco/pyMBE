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

import unittest as ut
import numpy as np
import pathlib
import pyMBE

mode="short" # Supported modes: "short", "long"
pH_samples=25 # If more through testing is needed, set to 200

class Test(ut.TestCase):
    data_root = pathlib.Path(__file__).parent / "henderson_hasselbalch_tests_data"

    def test(self):
        pmb = pyMBE.pymbe_library(seed=42)
        print("*** Running HH tests ***\n")

        # Peptide parameters
        sequence1 = 5 * "D" + 8 * "H"
        sequence2 = 3 * "E" + 7 * "R"
        pep1_concentration = 1e-2 * pmb.units.mol/pmb.units.L
        pep2_concentration = 1e-2 * pmb.units.mol/pmb.units.L
        c_salt=5e-3 * pmb.units.mol/ pmb.units.L
        model = '1beadAA'

        # Load pKa-values
        path_to_pka=pmb.root / "parameters" / "pka_sets" / "Nozaki1967.json"
        pmb.load_pka_set(path_to_pka)

        # Define the peptides in the pyMBE data frame
        pmb.define_peptide(name = "peptide_1",
                sequence = sequence1,
                model = model)

        pmb.define_peptide(name = "peptide_2",
                sequence = sequence2,
                model = model)

        with self.subTest(msg="Test edge cases"):
            # reference data
            ref_data_HH = np.loadtxt(self.data_root / "HH_no_pH_list.csv", delimiter=",")
            
            # Test that the function returns a list of None when no residues are defined 
            pH_values = [0, 14]
            pmb.define_molecule(name = "test",
                                residue_list = [])
            Z_HH = pmb.calculate_HH(molecule_name = "test",
                                      pH_list = pH_values)
            np.testing.assert_array_equal(Z_HH, 
                                          [None]*len(pH_values))
            
            pmb.define_molecule(name = "mol1",
                                residue_list=["TT"])
            Z_HH = pmb.calculate_HH(molecule_name = "mol1",
                                      pH_list = pH_values)
            np.testing.assert_array_equal(Z_HH, 
                                          [None]*len(pH_values))

            # Test that the function ignores residues with undefined particles
            pmb.define_residue(name = "RT",
                               central_bead="T",
                               side_chains=["TT"])
            pmb.define_molecule(name = "mol2",
                                residue_list=["RT"])
            Z_HH = pmb.calculate_HH(molecule_name = "mol2",
                                      pH_list = pH_values)
            np.testing.assert_array_equal(Z_HH, 
                                          [None]*len(pH_values))

            # Test that the function ignores undefined residues when other residues are defined
            pmb.define_peptide(name = "peptide_3",
                                sequence =sequence1+"T",
                                model= model)
            Z_HH_1 = pmb.calculate_HH(molecule_name = "peptide_3")
            np.testing.assert_allclose(Z_HH_1, ref_data_HH[0,:])

            

        with self.subTest(msg="Check Henderson-Hasselbalch equation"):
            # Check case where no pH_list is provided
            Z_HH_1 = pmb.calculate_HH(molecule_name = "peptide_1")
            Z_HH_2 = pmb.calculate_HH(molecule_name = "peptide_2")

            ref_data_HH = np.loadtxt(self.data_root / "HH_no_pH_list.csv", delimiter=",")
            np.testing.assert_allclose(Z_HH_1, ref_data_HH[0,:])
            np.testing.assert_allclose(Z_HH_2, ref_data_HH[1,:])


            # Check case where pH_list is provided
            pH_range = np.linspace(2, 12, num=200)[::200//pH_samples]
            Z_HH_1 = pmb.calculate_HH(molecule_name = "peptide_1",
                                      pH_list = pH_range)
            Z_HH_2 = pmb.calculate_HH(molecule_name = "peptide_2",
                                      pH_list = pH_range)

            ref_data_HH = np.loadtxt(self.data_root / "HH.csv", delimiter=",")
            np.testing.assert_allclose(Z_HH_1, ref_data_HH[0,::200//pH_samples])
            np.testing.assert_allclose(Z_HH_2, ref_data_HH[1,::200//pH_samples])

        with self.subTest(msg="Check Henderson-Hasselbalch equation with non-ionizable groups"):

            # Define additional non-ionizable groups
            pmb.define_particle(name = "N0",
                                z=0,
                                )
            pmb.define_particle(name = "N1",
                                z=1,
                                )
            path_to_pka=pmb.root / "parameters" / "pka_sets" / "Nozaki1967.json"
            pmb.load_pka_set(path_to_pka)
            pmb.define_residue(name = "RD",
                               central_bead="D",
                               side_chains=[])
            pmb.define_residue(name = "RH",
                               central_bead="H",
                               side_chains=[])
            pmb.define_residue(name = "RN0",
                               central_bead="N0",
                               side_chains=[])
            pmb.define_residue(name = "RN1",
                               central_bead="N1",
                               side_chains=[])

            # Load the reference data
            ref_data_HH = np.loadtxt(self.data_root / "HH_no_pH_list.csv", delimiter=",")
        

            # Check the case with non-ionizable groups without charge
            pmb.define_molecule(name = "mol_1",
                                residue_list = 5*["RD"] + 8*["RH"] + 3*["RN0"])
            Z_HH_1 = pmb.calculate_HH(molecule_name = "mol_1")
            np.testing.assert_allclose(Z_HH_1, 
                                       ref_data_HH[0,:])
           
            # Check the case with non-ionizable groups with charge
            pmb.define_molecule(name = "mol_2",
                                residue_list = 5*["RD"] + 8*["RH"] + 3*["RN1"])
            Z_HH_2 = pmb.calculate_HH(molecule_name = "mol_2")
            np.testing.assert_allclose(Z_HH_2, 
                                       ref_data_HH[0,:]+3)

        with self.subTest(msg="Check Henderson-Hasselbalch equation + Donnan"):
            # Check case where no pH_list is provided
            HH_Donnan_dict = pmb.calculate_HH_Donnan(
                    c_macro = {"peptide_1": pep1_concentration,
                               "peptide_2": pep2_concentration},
                    c_salt = c_salt)

            ref_data_HH_Donnan = np.loadtxt(self.data_root / "HH_Donnan_no_pH_list.csv", delimiter=",")
            np.testing.assert_allclose(HH_Donnan_dict["charges_dict"]["peptide_1"], ref_data_HH_Donnan[0,:])
            np.testing.assert_allclose(HH_Donnan_dict["charges_dict"]["peptide_2"], ref_data_HH_Donnan[1,:])


            # Check case where pH_list is provided
            HH_Donnan_dict = pmb.calculate_HH_Donnan(
                    c_macro = {"peptide_1": pep1_concentration,
                               "peptide_2": pep2_concentration},
                    c_salt = c_salt,
                    pH_list = pH_range)

            ref_data_HH_Donnan = np.loadtxt(self.data_root / "HH_Donnan.csv", delimiter=",")
            np.testing.assert_allclose(HH_Donnan_dict["charges_dict"]["peptide_1"], ref_data_HH_Donnan[0,::200//pH_samples])
            np.testing.assert_allclose(HH_Donnan_dict["charges_dict"]["peptide_2"], ref_data_HH_Donnan[1,::200//pH_samples])

        with self.subTest(msg="Check that HH and HH_Don are consistent"):
            Z_HH_1 = pmb.calculate_HH(molecule_name = "peptide_1",
                                      pH_list = HH_Donnan_dict["pH_system_list"])
            Z_HH_2 = pmb.calculate_HH(molecule_name = "peptide_2",
                                      pH_list = HH_Donnan_dict["pH_system_list"])

            np.testing.assert_allclose(Z_HH_1, HH_Donnan_dict["charges_dict"]["peptide_1"])
            np.testing.assert_allclose(Z_HH_2, HH_Donnan_dict["charges_dict"]["peptide_2"])


if __name__ == "__main__":
    ut.main()
