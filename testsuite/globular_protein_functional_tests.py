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

import pyMBE
from pyMBE.lib import analysis
import sys
import pathlib
import tempfile
import subprocess
import multiprocessing
import numpy as np
import pandas as pd
import unittest as ut
import glob 
import importlib.resources

root = pathlib.Path(__file__).parent.parent
data_root = root / "testsuite" / "globular_protein_tests_data"
script_path = root / "samples" / "Beyer2024" / "globular_protein.py"
test_pH_values = [2, 5, 7]
tasks = ["1beb", "1f6s"]
mode = "test"

def kernel(protein_pdb):
    """
    Runs a set of tests for a given protein pdb.

    Args:
        protein_pdb(`str`): PDB code of the protein.
    """
    with tempfile.TemporaryDirectory() as time_series_path:
        for pH in test_pH_values:
            print(f"pH = {pH}")
            run_command=[sys.executable, script_path, "--pdb", protein_pdb, "--pH", str(pH),
                         "--path_to_cg", importlib.resources.files(pyMBE) / "parameters" / "globular_proteins" / f"{protein_pdb}.vtf",
                         "--mode", "test", "--no_verbose", "--output", time_series_path]
            print(subprocess.list2cmdline(run_command))
            subprocess.check_output(run_command)
        # Analyze all time series
        data=analysis.analyze_time_series(path_to_datafolder=time_series_path,
                                          filename_extension="_time_series.csv")
    return (protein_pdb, data)


def kernel_move (protein_pdb):

    with tempfile.TemporaryDirectory() as time_series_path:

        run_command=[sys.executable, script_path, "--pdb", protein_pdb, "--pH", str(2),
                        "--path_to_cg", importlib.resources.files(pyMBE) / "parameters" / "globular_proteins" / f"{protein_pdb}.vtf",
                        "--mode", "test", "--output",  time_series_path, "--move_protein", "--no_verbose", "--ideal" ]
        print(subprocess.list2cmdline(run_command))
        subprocess.check_output(run_command)
        frame_folder = f"{time_series_path}/frames"
        list_files = glob.glob(f"{frame_folder}/*.vtf")
        coords={"first_frame": [],
                "last_frame": []}
        # Read the first and the last trajectory to check that the protein has moved
        for frame in list_files:
            num = int(frame[-5])
            with open(frame) as f:
                for line in f:
                    line_clean = line.split()
                    if line_clean: 
                        header = line_clean[0]
                        if header.isnumeric():
                            coord_part = line_clean[1:]
                            if int (num) == 0:
                                coords["first_frame"].append(coord_part)
                            elif int (num) == (len(list_files)-1):
                                coords["last_frame"].append(coord_part)
    return coords



class Test(ut.TestCase):
    
    def test_globular_protein(self):
        with multiprocessing.Pool(processes=2) as pool:
            results = dict(pool.map(kernel, tasks, chunksize=1))

        rtol=0.1 # relative tolerance
        atol=0.5 # absolute tolerance
        for protein_pdb, data in results.items():
            # Save data for future testing
            if mode == "save":
                data.to_csv(data_root / f"{protein_pdb}.csv", index=False)
                continue
            assert mode == "test", f"Mode {mode} not supported, valid modes: ['save', 'test']"
            with self.subTest(msg=f"Protein {protein_pdb}"):
                # Get reference test data
                ref_data=pd.read_csv(data_root / f"{protein_pdb}.csv", header=[0, 1])
                # Check charge
                test_charge=np.sort(data["mean","charge"].to_numpy())
                ref_charge=np.sort(ref_data["mean","charge"].to_numpy())
                np.testing.assert_allclose(
                    test_charge, ref_charge, rtol=rtol, atol=atol)
    

    def test_globular_protein_enable_motion(self):        

        coords=kernel_move("1beb")
        np.testing.assert_raises(AssertionError, 
                                 np.testing.assert_array_equal, 
                                 coords["first_frame"], 
                                 coords["last_frame"])

if __name__ == "__main__":
    ut.main()
