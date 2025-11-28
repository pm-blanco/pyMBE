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

parser = argparse.ArgumentParser(description='Plots the titration data from branched_polyampholyte.py and the corresponding analytical solution.')
parser.add_argument('--path_to_data',
                    type=pathlib.Path,
                    required= False,
                    default=pathlib.Path(__file__).parent / "time_series" / "branched_polyampholyte" / "analyzed_data.csv",
                    help='path to the analyzed data')
parser.add_argument('--output',
                    type=pathlib.Path,
                    required= False,
                    default=pathlib.Path(__file__).parent / "time_series" / "branched_polyampholyte",
                    help='output directory')
parser.add_argument('--mode',
                    type=str,
                    required= False,
                    default="plot",
                    choices=["plot", "store_HH"],
                    help='mode to execute the script; available options are: store_HH (stores the analytical HH solution, used for testing) and plot (produces the plots)')
args = parser.parse_args()

# Define the molecule (necessary to calculate the HH analytical solution with pyMBE)

# Acidic particle
pmb.define_particle(
    name = "A",
    acidity = "acidic",
    pka = 4,
    sigma = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))
    
# Basic particle
pmb.define_particle(
    name = "B",
    acidity = "basic",
    pka = 9,
    sigma = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))

# Define different residues
pmb.define_residue(
    name = "Res_1",
    central_bead = "I",
    side_chains = ["A","B"])
    
pmb.define_residue(
    name = "Res_2",
    central_bead = "I",
    side_chains = ["Res_1"])

# Define the molecule
pmb.define_molecule(
    name = "polyampholyte",
    residue_list = 2*["Res_1"] + ["Res_2"] + 2*["Res_1"] + 2*["Res_2"])

# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation
if args.mode == "plot":
    pH_range_HH = np.linspace(2, 12, num=100)
elif args.mode == "store_HH":
    pH_range_HH = [3.5,4.5,8.5,9.5]
Z_HH = pmb.calculate_HH(molecule_name="polyampholyte",
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
    plt.ylabel('Net charge  of the polyampholyte / e')
    plt.show()
    plt.close()

elif args.mode == "store_HH":
    HH_data=pd.DataFrame({"pH": pH_range_HH,
                         "Z_HH": Z_HH,})
    HH_data.to_csv(args.output / "HH_data.csv",
                        index=False)
