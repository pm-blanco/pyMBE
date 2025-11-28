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
from pyMBE.lib import analysis
import numpy as np
import pandas as pd
import unittest as ut

root = pathlib.Path(__file__).parent.parent
data_root = root / "testsuite" / "peptide_tests_data"
script_path = root / "samples" / "Beyer2024" / "peptide.py"
test_pH_values = [3, 7, 11]
tasks = [
  "K"*5+"D"*5, # K_5-D_5 case
  "E"*5+"H"*5, # E_5-H_5 case
  "nDSHAKRHHGYKRKFHEKHHSHRGYc", # histatin-5 case, slow simulation
]
mode = "test"

def kernel(sequence):
    """
    Runs a set of tests for a given peptide sequence.

    Args:
        sequence(`str`): Amino acid sequence of the peptide.
    """
    with tempfile.TemporaryDirectory() as time_series_path:
        for pH in test_pH_values:
            print(f"pH = {pH}")
            run_command=[sys.executable, script_path, "--sequence", sequence,
                         "--pH", str(pH), "--mode", "test", "--no_verbose",
                         "--output", time_series_path]
            print(subprocess.list2cmdline(run_command))
            subprocess.check_output(run_command)
        # Analyze all time series
        data=analysis.analyze_time_series(path_to_datafolder=time_series_path)
    return (sequence, data)


class Test(ut.TestCase):

    def test_peptide(self):
        with multiprocessing.Pool(processes=2) as pool:
            results = dict(pool.map(kernel, tasks, chunksize=2))

        rtol=0.1 # relative tolerance
        atol=0.5 # absolute tolerance
        for sequence, data in results.items():
            # Save data for future testing
            if mode == "save":
                data.to_csv(data_root / f"{sequence}.csv", index=False)
                continue
            assert mode == "test", f"Mode {mode} not supported, valid modes: ['save', 'test']"
            with self.subTest(msg=f"Sequence {sequence}"):
                # Get reference test data
                ref_data=pd.read_csv(data_root / f"{sequence}.csv", header=[0, 1])
                # Check charge
                test_charge=np.sort(data["mean","charge"].to_numpy())
                ref_charge=np.sort(ref_data["mean","charge"].to_numpy())
                np.testing.assert_allclose(test_charge, ref_charge, rtol=rtol, atol=atol)
                # Check rg
                test_rg=np.sort(data["mean","rg"].to_numpy())
                ref_rg=np.sort(ref_data["mean","rg"].to_numpy())
                np.testing.assert_allclose(test_rg, ref_rg, rtol=rtol, atol=atol)

if __name__ == "__main__":
    ut.main()
