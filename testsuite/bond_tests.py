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
import pyMBE
import numpy as np
import unittest as ut
import json.decoder
import json
import io
import logging
import pyMBE.storage.df_management as df_management

# Create an in-memory log stream
log_stream = io.StringIO()
logging.basicConfig(level=logging.INFO, 
                    format="%(levelname)s: %(message)s",
                    handlers=[logging.StreamHandler(log_stream)] )

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

class Test(ut.TestCase):

    def setUp(self):
        pmb.df = df_management._DFManagement._setup_df()

    def check_bond_setup(self, bond_object, input_parameters, bond_type):
        """
        Checks that pyMBE sets up a bond object correctly.

        Args:
                bond_object(`espressomd.interactions`): instance of a espresso bond object.
                input_parameters(`dict`): dictionary with the parameters for the bond potential.
                bond_type(`str`): label identifying the bond type
        """
        # Check that pyMBE stores the correct type of bond
        self.assertIn(bond_type.lower(), str(bond_object).lower(),
                      msg="pyMBE does not store the correct type of bond")
        # Check that pyMBE defines the bond with the correct parameters
        bond_params = bond_object.get_params()
        reduced_units = {'r_0'    : 'reduced_length',
                         'k'      : 'reduced_energy / reduced_length**2',
                         'd_r_max': 'reduced_length'}
        for key in input_parameters.keys():
            np.testing.assert_equal(actual=bond_params[key],
                    desired=input_parameters[key].m_as(reduced_units[key]),
                    verbose=True)
        
    def test_bond_harmonic(self):
        pmb.define_particle(name='A', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

        bond_type = 'harmonic'
        bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        # check default bond
        pmb.define_default_bond(bond_type = bond_type,
                                bond_parameters = bond)

        bond_object = pmb.filter_df(pmb_type='bond')['bond_object'].values[0]
        self.check_bond_setup(bond_object=bond_object,
                        input_parameters=bond,
                        bond_type=bond_type)

        # check particle bond
        pmb.define_bond(bond_type = bond_type,
                        bond_parameters = bond,
                        particle_pairs = [['A', 'A']])

        bond_object = pmb.filter_df(pmb_type='bond')['bond_object'].values[1]
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=bond,
                              bond_type=bond_type)

        # check bond deserialization
        bond_params =  bond_object.get_params()
        bond_params["bond_id"] = bond_object._bond_id
        deserialized = df_management._DFManagement._convert_str_to_bond_object(
            f'{bond_object.__class__.__name__}({json.dumps(bond_params)})')
        self.check_bond_setup(bond_object=deserialized,
                              input_parameters=bond,
                              bond_type=bond_type)

    def test_bond_fene(self):
        pmb.define_particle(name='A',
                             z=0, 
                             sigma=0.4*pmb.units.nm, 
                             epsilon=1*pmb.units('reduced_energy'))

        bond_type = 'FENE'
        bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                'd_r_max': 0.8 * pmb.units.nm}

        # check default bond
        pmb.define_default_bond(bond_type = bond_type,
                                bond_parameters = bond)

        bond_object = pmb.filter_df(pmb_type='bond')['bond_object'].values[0]
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=bond,
                              bond_type=bond_type)

        # check particle bond
        pmb.define_bond(bond_type = bond_type,
                        bond_parameters = bond,
                        particle_pairs = [['A', 'A']])

        bond_object = pmb.filter_df(pmb_type='bond')['bond_object'].values[1]
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=bond,
                              bond_type=bond_type)

        # check bond deserialization
        bond_params =  bond_object.get_params()
        bond_params["bond_id"] = bond_object._bond_id
        deserialized = df_management._DFManagement._convert_str_to_bond_object(
            f'{bond_object.__class__.__name__}({json.dumps(bond_params)})')
        self.check_bond_setup(bond_object=deserialized,
                              input_parameters=bond,
                              bond_type=bond_type)

        # check bond default equilibrium length
        bond_type = 'FENE'
        bond = {'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                'd_r_max': 0.8 * pmb.units.nm}

        pmb.define_bond(bond_type = bond_type,
                            bond_parameters = bond,
                            particle_pairs = [['A', 'A']])
        log_contents = log_stream.getvalue()
        self.assertIn("no value provided for r_0. Defaulting to r_0 = 0", log_contents)  
        bond['r_0'] = 0. * pmb.units.nm
        bond_object = pmb.filter_df(pmb_type='bond')['bond_object'].values[2]
        self.check_bond_setup(bond_object=bond_object,
                              input_parameters=bond,
                              bond_type=bond_type)

    def test_bond_harmonic_and_fene(self):
        pmb.define_particle(name='A', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
        pmb.define_particle(name='B', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

        bond_type_1 = 'harmonic'
        bond_1 = {'r_0'    : 0.4 * pmb.units.nm,
                  'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}
        pmb.define_bond(bond_type = bond_type_1,
                        bond_parameters = bond_1,
                        particle_pairs = [['A', 'A']])

        bond_type_2 = 'FENE'
        bond_2 = {'r_0'    : 0.4 * pmb.units.nm,
                  'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                  'd_r_max': 0.8 * pmb.units.nm}

        pmb.define_bond(bond_type = bond_type_2,
                        bond_parameters = bond_2,
                        particle_pairs = [['B', 'B']])

        bond_object_1 = pmb.filter_df(pmb_type='bond')['bond_object'][2]
        bond_object_2 = pmb.filter_df(pmb_type='bond')['bond_object'][3]

        self.check_bond_setup(bond_object=bond_object_1,
                              input_parameters=bond_1,
                              bond_type=bond_type_1)
        self.check_bond_setup(bond_object=bond_object_2,
                              input_parameters=bond_2,
                              bond_type=bond_type_2)

    def test_bond_raised_exceptions(self):
        pmb.define_particle(name='A', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
        for callback in [pmb.define_bond, pmb.define_default_bond]:
            with self.subTest(msg=f'using method {callback.__qualname__}()'):
                self.check_bond_exceptions(callback)

    def check_bond_exceptions(self, callback):
        # check exceptions for unknown bond types
        bond_type = 'Quartic'
        bond = {'r_0'    : 0.4 * pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type, "bond_parameters" : bond}
        if callback == pmb.define_bond:
            input_parameters["particle_pairs"] = [['A', 'A']]

        np.testing.assert_raises(NotImplementedError, callback, **input_parameters)

        # check exceptions for missing bond equilibrium length
        bond_type = 'harmonic'
        bond = {'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type, "bond_parameters" : bond}
        if callback == pmb.define_bond:
            input_parameters["particle_pairs"] = [['A', 'A']]

        np.testing.assert_raises(ValueError, callback, **input_parameters)

        # check exceptions for missing bond force constant
        bond_type = 'harmonic'
        bond = {'r_0'    : 0.4 * pmb.units.nm}

        input_parameters={"bond_type": bond_type, "bond_parameters" : bond}
        if callback == pmb.define_bond:
            input_parameters["particle_pairs"] = [['A', 'A']]

        np.testing.assert_raises(ValueError, callback, **input_parameters)

        # check exceptions for missing bond force constant
        bond_type = 'FENE'
        bond = {'r_0'    : 0.4*pmb.units.nm,
                'd_r_max': 0.8 * pmb.units.nm}

        input_parameters={"bond_type": bond_type, "bond_parameters" : bond}
        if callback == pmb.define_bond:
            input_parameters["particle_pairs"] = [['A', 'A']]

        np.testing.assert_raises(ValueError, callback, **input_parameters)

        # check exceptions for missing bond maximal length
        bond_type = 'FENE'
        bond = {'r_0'    : 0.4*pmb.units.nm,
                'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}

        input_parameters={"bond_type": bond_type, "bond_parameters" : bond}
        if callback == pmb.define_bond:
            input_parameters["particle_pairs"] = [['A', 'A']]

        np.testing.assert_raises(ValueError, callback, **input_parameters)

    def test_bond_framework(self):
        pmb.define_particle(name='A', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
        pmb.define_particle(name='B', z=0, sigma=0.4*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

        pmb.add_bonds_to_espresso(None)
        log_contents = log_stream.getvalue()
        assert "there are no bonds defined in pymbe.df" in log_contents

        bond_type_1 = 'harmonic'
        bond_1 = {'r_0'    : 0.4 * pmb.units.nm,
                  'k'      : 400 * pmb.units('reduced_energy / reduced_length**2')}
        pmb.define_bond(bond_type = bond_type_1,
                        bond_parameters = bond_1,
                        particle_pairs = [['A', 'A']])

        bond_type_2 = 'FENE'
        bond_2 = {'r_0'    : 0.4 * pmb.units.nm,
                  'k'      : 400 * pmb.units('reduced_energy / reduced_length**2'),
                  'd_r_max': 0.8 * pmb.units.nm}

        pmb.define_bond(bond_type = bond_type_2,
                        bond_parameters = bond_2,
                        particle_pairs = [['B', 'B']])

        bond_object_1 = pmb.filter_df(pmb_type='bond')['bond_object'][2]
        bond_object_2 = pmb.filter_df(pmb_type='bond')['bond_object'][3]

        self.check_bond_setup(bond_object=bond_object_1,
                              input_parameters=bond_1,
                              bond_type=bond_type_1)
        self.check_bond_setup(bond_object=bond_object_2,
                              input_parameters=bond_2,
                              bond_type=bond_type_2)

        # check deserialization exceptions
        with self.assertRaises(ValueError):
            df_management._DFManagement._convert_str_to_bond_object('Not_A_Bond()')
        with self.assertRaises(json.decoder.JSONDecodeError):
            df_management._DFManagement._convert_str_to_bond_object('HarmonicBond({invalid_json})')
        with self.assertRaises(NotImplementedError):
            df_management._DFManagement._convert_str_to_bond_object('QuarticBond({"r_0": 1., "k": 1.})')

        # check bond keys
        self.assertEqual(pmb.find_bond_key('A', 'A'), 'A-A')
        self.assertEqual(pmb.find_bond_key('B', 'B'), 'B-B')
        self.assertEqual(pmb.find_bond_key('A', 'A', use_default_bond=True), 'A-A')
        self.assertEqual(pmb.find_bond_key('Z', 'Z', use_default_bond=True), 'default')
        self.assertIsNone(pmb.find_bond_key('A', 'B'))
        self.assertIsNone(pmb.find_bond_key('B', 'A'))
        self.assertIsNone(pmb.find_bond_key('Z', 'Z'))
        self.assertEqual(pmb.find_bond_key('A', 'B', use_default_bond=True), 'default')
        
        self.assertIsNone(pmb.search_bond('A', 'B', hard_check=False))
        log_contents = log_stream.getvalue()
        self.assertIn("Bond not defined between particles A and B", log_contents)
        
        with self.assertRaises(ValueError):
            pmb.search_bond('A', 'B', use_default_bond=True)

        with self.assertRaises(ValueError):
            pmb.search_bond('A', 'B' , hard_check=True)

        # check invalid bond index
        pmb.add_value_to_df(key=('particle_id',''), new_value=10,
                            index=np.where(pmb.df['name']=='A')[0][0])
        pmb.add_value_to_df(key=('particle_id',''), new_value=20,
                            index=np.where(pmb.df['name']=='B')[0][0])
        self.assertIsNone(pmb.add_bond_in_df(10, 20, use_default_bond=False))
        self.assertIsNone(pmb.add_bond_in_df(10, 20, use_default_bond=True))

        # check bond lengths
        self.assertAlmostEqual(pmb.get_bond_length('A', 'A'),
                               bond_object_1.r_0, delta=1e-7)
        self.assertAlmostEqual(pmb.get_bond_length('B', 'B'),
                               bond_object_2.r_0, delta=1e-7)
        self.assertIsNone(pmb.get_bond_length('A', 'B'))
        log_contents = log_stream.getvalue()
        self.assertIn("Bond not defined between particles A and B", log_contents)
        with self.assertRaises(ValueError):
            pmb.get_bond_length('A', 'B', hard_check=True)


if __name__ == '__main__':
    ut.main()
