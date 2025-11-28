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

# Plots the charge of a peptide with a given sequence using either a custom or reference set of pKa values

import numpy as np
import matplotlib.pyplot as plt
import argparse
import pathlib
import pandas as pd
# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(seed=42)

parser = argparse.ArgumentParser(description='Plots the titration data from peptide_cpH.py and the corresponding analytical solution.')
parser.add_argument('--sequence',
                    type=str,
                    default= 'EEEEDDDD', 
                    help='sequence of the peptide')
parser.add_argument('--path_to_data',
                    type=pathlib.Path,
                    required= False,
                    default=pathlib.Path(__file__).parent / "time_series" / "peptide_cpH" / "analyzed_data.csv",
                    help='path to the analyzed data')
parser.add_argument('--output',
                    type=pathlib.Path,
                    required= False,
                    default=pathlib.Path(__file__).parent / "time_series" / "peptide_cpH",
                    help='output directory')
parser.add_argument('--mode',
                    type=str,
                    required= False,
                    default="plot",
                    choices=["plot", "store_HH"],
                    help='mode to execute the script; available options are: store_HH (stores the analytical HH solution, used for testing) and plot (produces the plots)')
args = parser.parse_args()

# Define peptide parameters
sequence = args.sequence
# Define the peptide in the pyMBE dataframe and load the pka set
# This is necessary to calculate the analytical solution from the Henderson-Hasselbach equation
peptide = 'generic_peptide'
pmb.define_peptide (name=peptide, 
                    sequence=sequence,
                   model="1beadAA") # not really relevant for plotting
path_to_pka=pmb.root / "parameters" / "pka_sets" / "Hass2015.json"
pmb.load_pka_set(path_to_pka)

# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation
if args.mode == "plot":
    pH_range_HH = np.linspace(2, 12, num=100)
elif args.mode == "store_HH":
    pH_range_HH = [2,4,5,6]
Z_HH = pmb.calculate_HH(molecule_name=peptide, 
                                  pH_list=pH_range_HH) 

if args.mode == "plot":
    # Read the analyzed data produced with peptide_mixture_grxmc_ideal
    analyzed_data = pd.read_csv(args.path_to_data, header=[0,1])

    # Plot the results
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.errorbar(analyzed_data["pH"]["value"], 
                analyzed_data["mean"]["charge"],
                yerr=analyzed_data["err_mean"]["charge"],
                fmt = 'o', 
                capsize=3,
                color="red", 
                label='cpH simulation')
    
    ax.plot(pH_range_HH, 
            Z_HH, 
            "--r", 
            label='HH analytical solution')
    
    plt.legend()
    plt.xlabel('pH')
    plt.ylabel('Net charge of the peptide / e')
    plt.show()
    plt.close()

elif args.mode == "store_HH":
    HH_data=pd.DataFrame({"pH": pH_range_HH,
                         "Z_HH": Z_HH,})
    HH_data.to_csv(args.output / "HH_data.csv",
                        index=False)
