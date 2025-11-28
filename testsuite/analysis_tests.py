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

import unittest as ut
import pandas as pd
import pyMBE.lib.analysis as ana
import pathlib


class Serialization(ut.TestCase):
    data_root = pathlib.Path(__file__).parent / "tests_data"

    def test_analyze_time_series(self):
        print("*** Unit test: test that analysis.analyze_time_series analyzes all data in a folder correctly ***")
        analyzed_data = ana.analyze_time_series(path_to_datafolder=self.data_root,
                                       filename_extension="_time_series.csv",
                                       minus_separator=True)
        analyzed_data[["Dens","eps"]] = analyzed_data[["Dens","eps"]].apply(pd.to_numeric)
        reference_data = pd.read_csv(self.data_root / "average_data.csv", header=[0,1])
        analyzed_data.columns = analyzed_data.sort_index(axis=1,level=[0,1],ascending=[True,True]).columns
        reference_data.columns = reference_data.sort_index(axis=1,level=[0,1],ascending=[True,True]).columns
        pd.testing.assert_frame_equal(analyzed_data.dropna(),reference_data.dropna(), check_column_type=False, check_dtype=False)
        print("*** Unit passed ***")
        print("*** Unit test: test that analysis.analyze_time_series ignores files that should not be analyzed ***")
        analyzed_data = ana.analyze_time_series(path_to_datafolder=self.data_root,
                                       ignore_files=["average_data.csv","N-064_Solvent-good_Init-coil_time_series_analyzed.csv"],
                                       minus_separator=True,
                                       filename_extension="_time_series.csv")
        analyzed_data[["Dens","eps"]] = analyzed_data[["Dens","eps"]].apply(pd.to_numeric)
        reference_data = pd.read_csv(self.data_root / "average_data.csv", header=[0,1])
        analyzed_data.columns = analyzed_data.sort_index(axis=1,level=[0,1],ascending=[True,True]).columns
        reference_data.columns = reference_data.sort_index(axis=1,level=[0,1],ascending=[True,True]).columns
        pd.testing.assert_frame_equal(analyzed_data.dropna(),reference_data.dropna(), check_column_type=False, check_dtype=False)
        print("*** Unit passed ***")
    
    def test_get_dt(self):
        print("*** Unit test: test that analysis.get_dt returns the right time step ***")
        data = pd.DataFrame.from_dict( {'time': [0, 1, 2,], 'obs': ['1.0', '2.0', '4.0']} )
        dt, n_warnings = ana.get_dt(data)
        self.assertAlmostEqual(dt, 1.0, delta = 1e-7)
        self.assertEqual(n_warnings, 0)
        print("*** Unit passed ***")       
        
        print("*** Unit test: test that analysis.get_dt prints a warning if there are values with different time steps ***")
        data = pd.DataFrame.from_dict( {'time': [0, 1, 1,], 'obs': ['1.0', '2.0', '4.0']} )
        dt, n_warnings = ana.get_dt(data,verbose=True)
        self.assertAlmostEqual(dt, 1.0, delta = 1e-7)
        self.assertEqual(n_warnings, 1)
        print("*** Unit passed ***")

        print("*** Unit test: test that analysis.get_dt raises a ValueError if the two first values have the same time step***")
        data = pd.DataFrame.from_dict( {'time': [0, 0, 1,], 'obs': ['1.0', '2.0', '4.0']} )
        inputs = {"data": data}
        self.assertRaises(ValueError, ana.get_dt, **inputs)
        print("*** Unit passed ***")

        print("*** Unit test: test that analysis.get_dt raises a ValueError if the column with the time is not found ***")
        data = pd.DataFrame.from_dict( {'ns': [0, 1, 2,], 'obs': ['1.0', '2.0', '4.0']} )
        inputs = {"data": data}
        self.assertRaises(ValueError, ana.get_dt, **inputs)
        
        print("*** Unit passed ***")

        print("*** Unit test: test that analysis.get_dt raises a ValueError if the time is stored at uneven intervals ***")
        data = pd.DataFrame.from_dict( {'time': [0, 1, 4,], 'obs': ['1.0', '2.0', '4.0']} )
        inputs = {"data": data}
        self.assertRaises(ValueError, ana.get_dt, **inputs)
        
        print("*** Unit passed ***")

    def test_add_data_to_df(self):
        print("*** Unit test: test that analysis.add_data_to_df creates a Pandas Dataframe from a dictionary correctly ***")
        data = {'A': [2], 
                'B': ['1.0']}
        reference_df = pd.DataFrame(data, 
                                    index=[0])
        
        analysis_df = ana.add_data_to_df(df=pd.DataFrame(),
                                         data_dict=data,
                                         index=[0])
        pd.testing.assert_frame_equal(reference_df,analysis_df)
        print("*** Unit passed ***")
        print("*** Unit test: test that analysis.add_data_to_df updates a Pandas Dataframe with new data from dictionary correctly ***")
        data ["C"] = False
        reference_df =  pd.concat([reference_df, pd.DataFrame(data,index=[len(analysis_df)])])
        analysis_df = ana.add_data_to_df(df=analysis_df,
                                         data_dict=data,
                                         index=[len(analysis_df)])
        print("*** Unit passed ***")


    def test_get_params_from_file_name(self):
        print("*** Unit test: test that get_params_from_file_name parses a filename without minus separator ***")
        filename = 'density_0.001_N_1000_T_2.00.csv'
        correct_params = {'density': '0.001', 'N': '1000', 'T': '2.00'}
        params = ana.get_params_from_file_name(filename, 
                                                minus_separator = False)
        self.assertEqual(correct_params,params)
        print("*** Unit passed ***")

        print("*** Unit test: test that get_params_from_file_name parses a filename with minus separator ***")
        filename = 'N-064_Solvent-good_Init-coil_time_series.csv'
        correct_params = {'N': 64, 'Solvent': 'good', 'Init': 'coil'}
        params = ana.get_params_from_file_name(filename, 
                                                minus_separator = True,
                                                filename_extension="_time_series.csv")
        self.assertEqual(correct_params,params)
        print("*** Unit passed ***")

        print("*** Unit test: test that get_params_from_file_name parses a filename with a different extension ***")
        filename = 'density_0.001_N_1000_T_2.00_time_series.txt'
        correct_params = {'density': '0.001', 'N': '1000', 'T': '2.00'}
        params = ana.get_params_from_file_name(filename, 
                                                minus_separator = False,
                                                filename_extension="_time_series.txt")
        self.assertEqual(correct_params,params)
        print("*** Unit passed ***")

        print("*** Unit test: test that get_params_from_file_name raises a ValueError if a filename with a wrong formating is provided ***")
        inputs = {"file_name": 'density_0.001_N_1000_T_f_2.00_time_series.txt',
                  "filename_extension": "_time_series.txt"}
        self.assertRaises(ValueError, ana.get_params_from_file_name, **inputs)
        print("*** Unit passed ***")
        
    def test_block_analyze(self):
        print("*** Unit test: test that block_analyze yields the expected outputs and reports the number of blocks and the block size. It should print that it encountered 1 repeated time value. ***")
        data = pd.read_csv(self.data_root / "N-064_Solvent-good_Init-coil_time_series.csv")
        analyzed_data = ana.block_analyze(full_data=data, verbose=True)
        analyzed_data = ana.add_data_to_df(df=pd.DataFrame(),
                            data_dict=analyzed_data.to_dict(),
                            index=[0])
        reference_data = pd.read_csv(self.data_root / "N-064_Solvent-good_Init-coil_time_series_analyzed.csv", header=[0,1])
        pd.testing.assert_frame_equal(analyzed_data.dropna(),reference_data.dropna(), check_column_type=False)
        print("*** Unit passed ***")
        
        print("*** Unit test: test that block_analyze analyzes correcly a subset of the data ***")
        analyzed_data = ana.block_analyze(full_data=data, columns_to_analyze="Rg")
        analyzed_data = ana.add_data_to_df(df=pd.DataFrame(),
                            data_dict=analyzed_data.to_dict(),
                            index=[0])
        reference_data = pd.read_csv(self.data_root / "N-064_Solvent-good_Init-coil_time_series_analyzed.csv", header=[0,1])
        reference_data = reference_data[[("mean","Rg"),("err_mean","Rg"),("n_eff","Rg"),("tau_int","Rg")]]
        pd.testing.assert_frame_equal(analyzed_data.dropna(),reference_data.dropna(), check_column_type=False)
        print("*** Unit passed ***")

        print("*** Unit test: test that block_analyze raises a ValueError if there is no time column ***")
        data = pd.DataFrame.from_dict( {'ns': [0, 1, 2,], 'obs': ['1.0', '2.0', '4.0']} )
        inputs = {"full_data": data, "verbose": False, "dt": 1}
        self.assertRaises(ValueError, ana.block_analyze, **inputs)
        
        print("*** Unit passed ***")
        

if __name__ == "__main__":
    print("*** lib.analysis unit tests ***")
    ut.main()
