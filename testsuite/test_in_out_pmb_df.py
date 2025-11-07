#
# Copyright (C) 2025 pyMBE-dev team
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

import pyMBE
import io
import pandas as pd
import unittest as ut
import logging
import re
import pyMBE.storage.df_management as df_management

# Create an in-memory log stream
log_stream = io.StringIO()
logging.basicConfig(level=logging.DEBUG, 
                    format="%(levelname)s: %(message)s",
                    handlers=[logging.StreamHandler(log_stream)])

pmb = pyMBE.pymbe_library(seed=42)
particle_parameters={"S1":{"name": "S1",
                                "sigma":0.355*pmb.units.nm, 
                                "epsilon":1*pmb.units('reduced_energy'),
                                    "z":0},
                            "S2":{"name": "S2",
                                "sigma":0.355*pmb.units.nm, 
                                "epsilon":1*pmb.units('reduced_energy'),
                                    "z":1},
                            "S3":{"name": "S3",
                                "sigma":0.355*pmb.units.nm, 
                                "epsilon":1*pmb.units('reduced_energy'),
                                    "z":2}}

pmb.define_particles(parameters=particle_parameters)

generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
generic_bond_length = 0.355*pmb.units.nm
HARMONIC_parameters = {'r_0'    : generic_bond_length,
                    'k'      : generic_harmonic_constant}
pmb.define_bond(bond_type = 'harmonic',
                bond_parameters = HARMONIC_parameters, particle_pairs = [["S1", "S2"]])

class Serialization(ut.TestCase):
    def test_add_to_df(self):
        print("*** Unit test: test that check if the add_value_to_df() raises two warnings in default mode (overwrite=False)  ***")
        index=0
        key = ('name','')
        old_value = pmb.df.loc[index,pd.IndexSlice[key]]
        new_value='T2'
        name=pmb.df.loc[index,key]
        pmb_type=pmb.df.loc[index,('pmb_type','')]
        df_management._DFManagement._add_value_to_df(df = pmb.df,
                                                     index = index,
                                                     key = key,
                                                     new_value = new_value)
        log_contents = log_stream.getvalue()
        warning_msg1=f"You are attempting to redefine the properties of {name} of pmb_type {pmb_type}"
        warning_msg2=f"pyMBE has preserved of the entry `{key}`: old_value = {old_value}. If you want to overwrite it with new_value = {new_value}, activate the switch overwrite = True"
        # Check if the warnings are in the log
        assert re.search(re.escape(warning_msg1), log_contents)
        assert re.search(re.escape(warning_msg2), log_contents)
        print("*** Unit passed ***")
        print("*** Unit test: test that check if the add_value_to_df() raises one warning when overwrite=True ***")
        index=0
        key = ('name','')
        old_value = pmb.df.loc[index,pd.IndexSlice[key]]
        new_value='T2'
        name=pmb.df.loc[index,key]
        pmb_type=pmb.df.loc[index,('pmb_type','')]
        df_management._DFManagement._add_value_to_df(df = pmb.df,
                                                     index = index,
                                                     key = key,
                                                     new_value = new_value,
                                                     overwrite = True)
        log_contents = log_stream.getvalue()
        warning_msg1=f"You are attempting to redefine the properties of {name} of pmb_type {pmb_type}"
        warning_msg2=f"Overwritting the value of the entry `{key}`: old_value = {old_value} new_value = {new_value}"
        # Check if the warnings are in the log
        assert re.search(re.escape(warning_msg1), log_contents)
        assert re.search(re.escape(warning_msg2), log_contents)
        print("*** Unit passed ***")
        
    def test_delete_entries_df(self):
        print("*** Unit test: test that entries in df are deleted properly using `delete_entries_in_df`  ***")
        
        pmb.df = df_management._DFManagement._delete_entries_in_df(df=pmb.df,
                                                                    entry_name="S1-S2")
        assert pmb.df[pmb.df["name"]=="S1-S2"].empty
        pmb.df = df_management._DFManagement._delete_entries_in_df(df=pmb.df,
                                                                   entry_name="S1")
        assert pmb.df[pmb.df["name"]=="S1"].empty

        residue_parameters={"R1":{"name": "R1",
                                "central_bead": "S2",
                                "side_chains": []},
                            "R2":{"name": "R2",
                                "central_bead": "S2",
                                "side_chains": ["S2","S3"]}}

        for parameter_set in residue_parameters.values():
            pmb.define_residue(**parameter_set)

        pmb.df = df_management._DFManagement._delete_entries_in_df(df=pmb.df,
                                                                   entry_name="R1")
        assert pmb.df[pmb.df["name"]=="R1"].empty

        molecule_parameters={"M1":{"name": "M1",
                            "residue_list": ["R2","R2","R2"]}}

        for parameter_set in molecule_parameters.values():
            pmb.define_molecule(**parameter_set)

        pmb.df = df_management._DFManagement._delete_entries_in_df(df=pmb.df,
                                                                   entry_name="M1")
        assert pmb.df[pmb.df["name"]=="M1"].empty
        print("*** Unit passed ***")

if __name__ == "__main__":
    ut.main()

