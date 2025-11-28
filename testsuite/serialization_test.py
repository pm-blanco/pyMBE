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

import json
import unittest as ut
import numpy as np
import pandas as pd
import pyMBE
import pyMBE.lib.analysis
import scipy.constants
import pyMBE.storage.df_management as df_management

class Serialization(ut.TestCase):

    def test_json_encoder(self):
        encoder = df_management._DFManagement._NumpyEncoder
        # Python types
        self.assertEqual(json.dumps(1, cls=encoder), "1")
        self.assertEqual(json.dumps([1, 2], cls=encoder), "[1, 2]")
        self.assertEqual(json.dumps((1, 2), cls=encoder), "[1, 2]")
        self.assertEqual(json.dumps({1: 2}, cls=encoder), """{"1": 2}""")
        # NumPy types
        self.assertEqual(json.dumps(np.array([1, 2]), cls=encoder), "[1, 2]")
        self.assertEqual(json.dumps(np.array(1), cls=encoder), "1")
        self.assertEqual(json.dumps(np.int32(1), cls=encoder), "1")
        # Pandas types
        with self.assertRaisesRegex(TypeError, "Object of type Series is not JSON serializable"):
            json.dumps(pd.Series([1, 2]), cls=encoder)

    def test_parameters_to_path(self):
        params = {"kT": 2., "phi": -np.pi, "n": 3, "fene": True, "name": "pep"}
        name = pyMBE.lib.analysis.built_output_name(params)
        self.assertEqual(name, "kT_2_phi_-3.14_n_3_fene_True_name_pep")
        params_out = pyMBE.lib.analysis.get_params_from_file_name(name)
        params_ref = {"kT": "2", "phi": "-3.14", "n": "3",
                      "fene": "True", "name": "pep"}
        self.assertEqual(params_out, params_ref)

    def test_pint_units(self):
        ref_output =  "\n".join(["Current set of reduced units:",
                                 "0.355 nanometer = 1 reduced_length",
                                 "4.1164e-21 joule = 1 reduced_energy",
                                 "1.6022e-19 coulomb = 1 reduced_charge",
                                 "Temperature: 298.15 kelvin"
                                ]) 
        pmb = pyMBE.pymbe_library(seed=42)
        reduced_units = pmb.get_reduced_units()
        self.assertEqual(reduced_units, ref_output)
        np.testing.assert_allclose(
            [pmb.kB.magnitude, pmb.N_A.magnitude, pmb.e.magnitude],
            [scipy.constants.k, scipy.constants.N_A, scipy.constants.e],
            rtol=1e-8, atol=0.)
        self.assertAlmostEqual((pmb.kT / pmb.kB).magnitude, 298.15, delta=1e-7)
        self.assertAlmostEqual((pmb.kT / scipy.constants.k).magnitude, 298.15,
                               delta=1e-7)


if __name__ == "__main__":
    ut.main()
