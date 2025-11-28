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
#

import subprocess
import pathlib
import pyMBE
import sys

pmb = pyMBE.pymbe_library(seed=23)

# Runs `samples/weak_polyacid_hydrogel_grxmc.py` to reproduce selected cases from Landsgesell2022
samples_path = pathlib.Path(__file__).parent.parent
script_path=samples_path / "weak_polyacid_hydrogel_grxmc.py"
data_path=samples_path / "Landsgesell2022" / "time_series"
mode="short-run"
mpc=40
pKa_value=4
# Cases for the PV curve
c_salt_res=0.05
L_fractions=[0.30259595,0.50963529]
pH_value = 5

for L_fraction in L_fractions:
    run_command=[sys.executable, script_path, 
                 "--mpc", str(mpc), 
                 "--pH_res", str(pH_value),
                 "--csalt_res", str(c_salt_res), 
                 "--pKa", str(pKa_value),
                 "--L_fraction", str(L_fraction),
                 "--mode", mode,
                 "--output", data_path]
    print(subprocess.list2cmdline(run_command))
    subprocess.check_output(run_command)

# Cases for the alpha vs pH curve
c_salt_res=0.01
#key = pH-value , value = equilibrium L fraction (swelling equilibrium of the gel at that pH value)
swelling_eq={4: 0.36630036630074175,
             6: 0.5574136008913612}
for pH_value in swelling_eq.keys():
    run_command=[sys.executable, script_path, 
                 "--mpc", str(mpc), 
                 "--pH_res", str(pH_value),
                 "--csalt_res", str(c_salt_res), 
                 "--pKa", str(pKa_value),
                 "--L_fraction", str(swelling_eq[pH_value]),
                 "--mode", mode,
                 "--output", data_path]
    print(subprocess.list2cmdline(run_command))
    subprocess.check_output(run_command)

# Analyze the data
script_path=samples_path / "analyze_time_series.py"
run_command=[sys.executable, script_path, 
                 "--data_folder", data_path]
print(subprocess.list2cmdline(run_command))
subprocess.check_output(run_command)
