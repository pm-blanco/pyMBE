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

import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import argparse
import pathlib
import pyMBE
from scipy import interpolate

parser = argparse.ArgumentParser(description='Plots alpha vs pH  from weak_gel.py and the corresponding reference data from Landsgesell2022.')
parser.add_argument('--path_to_data',
                    type=pathlib.Path,
                    required= False,
                    default=pathlib.Path(__file__).parent / "time_series" / "analyzed_data.csv",
                    help='path to the analyzed data')
args = parser.parse_args()

pmb = pyMBE.pymbe_library(seed=42)

# Values from Landsgesell2022
unit_length=0.355*pmb.units.nm
temperature=300*pmb.units.K
pmb.set_reduced_units(unit_length=unit_length,
                      temperature=temperature,
                      verbose=False)
ref_cs = pmb.units.Quantity(0.00134770889, "1/reduced_length**3")
ref_pH = 5
ref_max_L = 89.235257606948  # Maximum chain length for mpc=39
# Read the reference data
data_path = pathlib.Path(__file__).parent.parent.parent / "testsuite" / "data"
data_ref = pd.read_csv(data_path / "Landsgesell2022a.csv")
pressures_ref = pmb.units.Quantity((data_ref[(data_ref["pH"] == 5) & (data_ref["cs_bulk"] == 0.00134770889)])["total_isotropic_pressures"].values,"reduced_energy/reduced_length**3")
box_l_ref     = data_ref[(data_ref["pH"] == 5) & (data_ref["cs_bulk"] == 0.00134770889)]["#final_box_l"].values

# Read the analyzed data
analyzed_data = pd.read_csv(args.path_to_data, header=[0,1])
analyzed_data = analyzed_data[(analyzed_data["csalt"]["value"] == round((ref_cs/pmb.N_A).m_as("M"),2)) & (analyzed_data["pH"]["value"] == ref_pH)]
analyzed_preassure = pmb.units.Quantity(analyzed_data["mean"]["pressure"].values, "reduced_energy/reduced_length**3")
analyzed_preassure_err = pmb.units.Quantity(analyzed_data["err_mean"]["pressure"].values, "reduced_energy/reduced_length**3")


# Calculate pressure in the reservour
## load the monovalent salt reference data
data_path = pmb.root / "parameters" / "salt"
monovalent_salt_ref_data=pd.read_csv(data_path / "excess_chemical_potential_excess_pressure.csv")
cs_bulk = pmb.units.Quantity(monovalent_salt_ref_data['cs_bulk_[1/sigma^3]'].values,"1/reduced_length**3")
excess_press = pmb.units.Quantity(monovalent_salt_ref_data['excess_pressure_[kT/sigma^3]'].values, "reduced_energy/reduced_length**3")
excess_press = interpolate.interp1d(cs_bulk.m_as("1/L"), 
                                    excess_press.m_as("bar"))
p_id = 2*ref_cs*pmb.kT
p_res = p_id + excess_press(ref_cs.m_as("1/L"))*pmb.units.bar

# Plot
plt.rcParams["font.family"] = "serif"
plt.tight_layout()
mpl.rc('axes', linewidth=1)
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['lines.linewidth'] = 1.0

width = 8
height = 8
plt.figure(figsize=(width,height), dpi=300)
plt.grid(which='major', 
            color='#CCCCCC', 
            linestyle='--', 
            linewidth=0.6)
plt.errorbar(analyzed_data["Lfraction"]["value"], 
             (analyzed_preassure-p_res).m_as("bar"), 
             yerr=analyzed_preassure_err.m_as("bar"), 
             fmt="o", 
             capsize=3, 
             elinewidth=2,
             ecolor='tab:blue', 
             label="pyMBE")

plt.plot(box_l_ref/ref_max_L, 
         (pressures_ref-p_res).m_as("bar"), 
         marker="*",
         label="Landsgesell et al. 2022")
plt.xlabel(r"$L/L_{max}$", fontsize=17)
plt.ylabel(r"$P_{sys}-P_{res}$ [bar]", fontsize=17)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(fontsize=17)
plt.savefig("./PV_curve.pdf", 
            bbox_inches='tight')
plt.close()