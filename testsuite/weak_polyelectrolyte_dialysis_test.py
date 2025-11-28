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

from pyMBE.lib import analysis
import sys
import pathlib
import tempfile
import subprocess
import numpy as np
import pandas as pd
import unittest as ut

root = pathlib.Path(__file__).parent.parent
data_root = root / "testsuite" / "weak_polyelectrolyte_dialysis_test_data"
script_path = root / "samples" / "Beyer2024" / "weak_polyelectrolyte_dialysis.py"

test_pH_values=[3,5,7,9]
c_salt_res=0.01
c_mon_sys=0.435
pKa_value=4.0


def kernel():
    with tempfile.TemporaryDirectory() as time_series_path:
        for pH in test_pH_values:
            print(f"pH = {pH}")
            run_command=[sys.executable, script_path, "--c_salt_res", str(c_salt_res),
                         "--c_mon_sys", str(c_mon_sys), "--pKa_value", str(pKa_value),
                         "--pH_res", str(pH), "--mode", "test", "--output",
                         time_series_path, "--no_verbose"]
            print(subprocess.list2cmdline(run_command))
            subprocess.check_output(run_command)
        # Analyze all time series
        data=analysis.analyze_time_series(path_to_datafolder=time_series_path,
                                          filename_extension="_time_series.csv")
    return data


class Test(ut.TestCase):

    def test_polyelectrolyte_dialysis(self):
        """
        Test weak polyelectrolyte dialysis with G-RxMC (interacting).
        """
        rtol=0.1 # relative tolerance
        atol=0.05 # absolute tolerance
        data = kernel()
        # Get reference test data
        ref_data=pd.read_csv(data_root / "data.csv", header=[0, 1])
        # Check charge
        test_charge=np.sort(data["mean","alpha"].to_numpy())
        ref_charge=np.sort(ref_data["mean","alpha"].to_numpy())
        np.testing.assert_allclose(test_charge, ref_charge, rtol=rtol, atol=atol)


if __name__ == "__main__":
    ut.main()
