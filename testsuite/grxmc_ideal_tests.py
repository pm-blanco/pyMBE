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

# Tests that samples/peptide_mixture_grxmc_ideal.py, analysis.py, plot_peptide_mixture_grxmc_ideal.py work properly
# Import pyMBE and other libraries
import sys
import subprocess
import numpy as np
import sys
import pathlib
import tempfile
import subprocess
import multiprocessing
import pandas as pd
import unittest as ut

root = pathlib.Path(__file__).parent.parent
data_root = root / "samples" / "time_series" / "peptide_mixture_grxmc_ideal"
sample_path = root / "samples" / "peptide_mixture_grxmc_ideal.py"
analysis_path = root / "samples" / "analyze_time_series.py"
ideal_path = root / "samples" / "plot_peptide_mixture_grxmc_ideal.py"
test_pH_values=[2,5,7,10,12]
rtol=0.01 # relative tolerance
atol=0.05 # absolute tolerance

def kernel(pH_value,temp_dir_path,mode):
    """
    Runs a set of tests for a given peptide sequence.

    Args:
        pH_value(`float`): pH of the media.
        temp_dir_path(`str`): path of the folder were to output the data.
        mode(`str`): mode for the setup of the grxmc mode, supported modes are "standard" and "unified"
    """
    run_command=[sys.executable, sample_path, "--mode", mode,
                    "--pH", str(pH_value), "--output", temp_dir_path, "--test"]
    print(subprocess.list2cmdline(run_command))
    subprocess.check_output(run_command)
    return 

def analyze_and_test_data(temp_dir_path):
    """
    Analyzes the data and checks that is consistent against the HH analytical result.

    Args:
        temp_dir_path(`str`): path of the folder were to output the data.
    """
    # Analyze the data
    run_command=[sys.executable, analysis_path, "--data_folder", temp_dir_path]
    subprocess.check_output(run_command)
    analyzed_data=pd.read_csv(temp_dir_path + "/analyzed_data.csv", header=[0, 1])
    # Produce the ideal data
    run_command=[sys.executable, ideal_path, "--output", temp_dir_path, "--mode", "store_HH"]
    subprocess.check_output(run_command)
    HH_data=pd.read_csv(temp_dir_path + "/HH_data.csv") 
    np.testing.assert_allclose(np.sort(analyzed_data["mean"]["charge_peptide1"].to_numpy()), 
                               np.sort(HH_data["Z_HH_peptide1"].to_numpy()), 
                               rtol=rtol, 
                               atol=atol)
    np.testing.assert_allclose(np.sort(analyzed_data["mean"]["charge_peptide2"].to_numpy()), 
                               np.sort(HH_data["Z_HH_peptide2"].to_numpy()), 
                               rtol=rtol, 
                               atol=atol)
    np.testing.assert_allclose(np.sort(analyzed_data["mean"]["xi_plus"].to_numpy()), 
                               np.sort(HH_data["xi_HH"].to_numpy()), 
                               rtol=rtol*5, 
                               atol=atol*5)    

class Test(ut.TestCase):
    def test_standard_grxmc(self):
        temp_dir=tempfile.TemporaryDirectory()
        mode="standard"
        N_cases=len(test_pH_values)
        with multiprocessing.Pool(processes=2) as pool:
            pool.starmap(kernel, zip(test_pH_values,[temp_dir.name]*N_cases,[mode]*N_cases), chunksize=2)[0]
        analyze_and_test_data(temp_dir_path=temp_dir.name)
        temp_dir.cleanup()
    def test_unified_grxmc(self):
        temp_dir=tempfile.TemporaryDirectory()
        mode="unified"
        N_cases=len(test_pH_values)
        with multiprocessing.Pool(processes=2) as pool:
            pool.starmap(kernel, zip(test_pH_values,[temp_dir.name]*N_cases,[mode]*N_cases), chunksize=2)[0]
        analyze_and_test_data(temp_dir_path=temp_dir.name)
        temp_dir.cleanup()

if __name__ == "__main__":
    ut.main()


