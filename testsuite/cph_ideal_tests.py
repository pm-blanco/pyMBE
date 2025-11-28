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

import sys
import pathlib
import tempfile
import subprocess
import multiprocessing
import pandas as pd
import unittest as ut
import numpy as np

root = pathlib.Path(__file__).parent.parent

def kernel(pH_value,temp_dir_path,sample_path):
    """
    Runs a set of tests for a given peptide sequence.

    Args:
        pH_value(`float`): pH of the media.
        temp_dir_path(`str`): path to the folder were to output the data.
        sample_path(`str`): path to the sample script to produce the time series.
    """
    run_command=[sys.executable, sample_path,
                    "--pH", str(pH_value), "--output", temp_dir_path, "--test"]
    print(subprocess.list2cmdline(run_command))
    subprocess.check_output(run_command)
    return 

def analyze_and_test_data(temp_dir_path, analysis_path, plot_path, rtol, atol):
    """
    Analyzes the data and checks that is consistent against the HH analytical result.

    Args:
        temp_dir_path(`str`): path of the folder were to output the data.
        analysis_path(`str`): path to the script to analyze the time series.
        plot_path(`str`): path to the plotting script, used to calculate the analytical solution.
        atol(`float`): absolute tolerance for the numerical test.
        rtol(`float`): relative tolerance for the numerical test.
    """
    # Analyze the data
    run_command=[sys.executable, analysis_path, "--data_folder", temp_dir_path]
    subprocess.check_output(run_command)
    analyzed_data=pd.read_csv(temp_dir_path + "/analyzed_data.csv", header=[0, 1])
    # Produce the ideal data
    run_command=[sys.executable, plot_path, "--output", temp_dir_path, "--mode", "store_HH"]
    subprocess.check_output(run_command)
    HH_data=pd.read_csv(temp_dir_path + "/HH_data.csv") 
    np.testing.assert_allclose(np.sort(analyzed_data["mean"]["charge"].to_numpy()), 
                               np.sort(HH_data["Z_HH"].to_numpy()), 
                               rtol=rtol, 
                               atol=atol)
    
class Test(ut.TestCase):
    def test_peptide(self):
        sample_path = root / "samples" / "peptide_cpH.py"
        analysis_path = root / "samples" / "analyze_time_series.py"
        plot_path = root / "samples" / "plot_peptide_cpH.py"
        test_pH_values=[2,4,5,6]
        rtol=0.01 # relative tolerance
        atol=0.05 # absolute tolerance
        temp_dir=tempfile.TemporaryDirectory()
        N_cases=len(test_pH_values)
        with multiprocessing.Pool(processes=2) as pool:
            pool.starmap(kernel, zip(test_pH_values,[temp_dir.name]*N_cases,[sample_path]*N_cases), chunksize=2)[0]
        analyze_and_test_data(temp_dir_path=temp_dir.name,
                              analysis_path=analysis_path, 
                              plot_path=plot_path, 
                              rtol=rtol, 
                              atol=atol)
        temp_dir.cleanup()

    def test_branched_polyampholyte(self):
        sample_path = root / "samples" / "branched_polyampholyte.py"
        analysis_path = root / "samples" / "analyze_time_series.py"
        plot_path = root / "samples" / "plot_branched_polyampholyte.py"
        test_pH_values=[3.5,4.5,8.5,9.5]
        rtol=0.01 # relative tolerance
        atol=0.05 # absolute tolerance
        temp_dir=tempfile.TemporaryDirectory()
        N_cases=len(test_pH_values)
        with multiprocessing.Pool(processes=2) as pool:
            pool.starmap(kernel, zip(test_pH_values,[temp_dir.name]*N_cases,[sample_path]*N_cases), chunksize=2)[0]
        analyze_and_test_data(temp_dir_path=temp_dir.name,
                              analysis_path=analysis_path, 
                              plot_path=plot_path, 
                              rtol=rtol, 
                              atol=atol)
        temp_dir.cleanup()
    
if __name__ == "__main__":
    ut.main()
