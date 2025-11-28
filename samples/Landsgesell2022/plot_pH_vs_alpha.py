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

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pyMBE
import pandas as pd
import argparse
import pathlib

parser = argparse.ArgumentParser(description='Plots alpha vs pH  from weak_gel.py and the corresponding reference data from Landsgesell2022.')
parser.add_argument('--path_to_data',
                    type=pathlib.Path,
                    required= False,
                    default=pathlib.Path(__file__).parent / "time_series" / "analyzed_data.csv",
                    help='path to the analyzed data')
args = parser.parse_args()

pmb = pyMBE.pymbe_library(seed=42)

ref_cs = 0.01

# Read the reference data
data_path = pathlib.Path(__file__).parent.parent.parent / "testsuite" / "data"
data_ref = pd.read_csv(data_path / "Landsgesell2022a.csv")
data_ref['cs'] = pd.to_numeric(data_ref['cs'], errors='coerce')
data_ref = data_ref[np.isclose(data_ref['cs'], ref_cs)].sort_values(by="pH")
# Read the analyzed time series
analyzed_data = pd.read_csv(args.path_to_data, header=[0,1])
analyzed_data_to_plot={"pH" : [],
                       "alpha": [],
                       "alpha_err": []    
                       }
for pH_value in data_ref["pH"]:
    sorted_ref_df = data_ref[data_ref["pH"] == pH_value]
    reference_L_fraction = sorted_ref_df["best_L"].values[0]/sorted_ref_df["max_L"].values[0]
    mask=(analyzed_data["pH"]["value"] == pH_value) & (abs(analyzed_data["Lfraction"]["value"] - reference_L_fraction) < 1e3) & (analyzed_data["csalt"]["value"] == ref_cs)
    sorted_analyzed_df = analyzed_data[mask]
    if not sorted_analyzed_df.empty:
        analyzed_data_to_plot["pH"].append(sorted_analyzed_df["pH"]["value"].values[0])
        analyzed_data_to_plot["alpha"].append(sorted_analyzed_df["mean"]["alpha"].values[0])
        analyzed_data_to_plot["alpha_err"].append(sorted_analyzed_df["err_mean"]["alpha"].values[0])


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
plt.xlabel("pH in the reservoir",fontsize=17)
plt.ylabel(r"degree of ionization $\alpha$", fontsize=17)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.errorbar(analyzed_data_to_plot["pH"], 
             analyzed_data_to_plot["alpha"], 
             yerr = analyzed_data_to_plot["alpha_err"], 
             fmt='o', 
             capsize=3, 
             elinewidth=2, 
             ecolor='tab:blue', 
             label=r"pyMBE")
plt.plot(data_ref["pH"].drop_duplicates(),
         data_ref["alpha"].drop_duplicates(),
         marker="*",
         label="Landsgesell et al. 2022")

plt.legend(fontsize=17, loc="lower right")
plt.savefig("./alpha_pH_curve.pdf", 
            bbox_inches='tight')
plt.close()
