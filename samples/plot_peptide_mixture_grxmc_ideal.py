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
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import pathlib

parser = argparse.ArgumentParser(description='Plots the titration data from peptide.py and the corresponding analytical solution.')
parser.add_argument('--sequence1',
                    type=str,
                    default= 'nEHc', 
                    help='sequence of the first peptide')
parser.add_argument('--pep1_conc',
                    type=float,
                    default= 1e-2, 
                    help='concentration of the first peptide, in mol/L')
parser.add_argument('--pep2_conc',
                    type=float,
                    default= 1e-2, 
                    help='concentration of the second peptide, in mol/L')
parser.add_argument('--sequence2',
                    type=str,
                    default= 'nEEHHc', 
                    help='sequence of the second peptide')
parser.add_argument('--csalt',
                    type=float,
                    default= 5e-3, 
                    help='Added salt concentration in the system, in mol/L')
parser.add_argument('--path_to_data',
                    type=pathlib.Path,
                    required= False,
                    default=pathlib.Path(__file__).parent / "time_series" / "peptide_mixture_grxmc_ideal" / "analyzed_data.csv",
                    help='path to the analyzed data')
parser.add_argument('--output',
                    type=pathlib.Path,
                    required= False,
                    default=pathlib.Path(__file__).parent / "time_series" / "peptide_mixture_grxmc_ideal",
                    help='output directory')
parser.add_argument('--mode',
                    type=str,
                    required= False,
                    default="plot",
                    choices=["plot", "store_HH"],
                    help='mode to execute the script; available options are: store_HH (stores the analytical HH solution, used for testing) and plot (produces the plots)')
args = parser.parse_args()

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)
c_salt=args.csalt * pmb.units.mol/ pmb.units.L

# Define peptide parameters
sequence1 = args.sequence1
pep1_concentration = args.pep1_conc *pmb.units.mol/pmb.units.L

sequence2 = args.sequence2
pep2_concentration = args.pep2_conc *pmb.units.mol/pmb.units.L

# Define the peptides in the pyMBE data frame and load the pka set
# This is necesary to calculate the analytical solution from the Henderson-Hasselbach equation
peptide1 = 'generic_peptide1'
pmb.define_peptide (name=peptide1, 
                    sequence=sequence1,
                   model="1beadAA") # not really relevant for plotting
peptide2 = 'generic_peptide2'
pmb.define_peptide (name=peptide2, 
                    sequence=sequence2, 
                    model="1beadAA") # not really relevant for plotting
path_to_pka=pmb.root / "parameters" / "pka_sets" / "Hass2015.json"
pmb.load_pka_set(path_to_pka)

# Calculate the ideal titration curve of the peptide with Henderson-Hasselbach equation
if args.mode == "plot":
    pH_range_HH = np.linspace(2, 12, num=100)
elif args.mode == "store_HH":
    pH_range_HH = [2,5,7,10,12]
HH_charge_dict = pmb.calculate_HH_Donnan(c_macro={peptide1: pep1_concentration, peptide2: pep2_concentration}, 
                                            c_salt=c_salt, 
                                            pH_list=pH_range_HH)
Z_HH_Donnan = HH_charge_dict["charges_dict"]
xi_HH = HH_charge_dict["partition_coefficients"]

if args.mode == "plot":
    # Read the analyzed data produced with peptide_mixture_grxmc_ideal
    analyzed_data = pd.read_csv(args.path_to_data, header=[0,1])

    # Plot the results
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.errorbar(analyzed_data["pH"]["value"], 
                analyzed_data["mean"]["charge_peptide1"],
                yerr=analyzed_data["err_mean"]["charge_peptide1"],
                fmt = 'o', 
                capsize=3,
                color="red", 
                label='peptide 1 (GRxmC simulation)')
    ax.errorbar(analyzed_data["pH"]["value"], 
                analyzed_data["mean"]["charge_peptide2"],
                yerr=analyzed_data["err_mean"]["charge_peptide2"],
                fmt = 'o', 
                capsize=3,
                color="green", 
                label='peptide 2 (GRxMC simulation)')
    ax.plot(pH_range_HH, 
            Z_HH_Donnan[peptide1], 
            "--r", 
            label='peptide 1 (Theory: HH+Donnan)')
    ax.plot(pH_range_HH, 
            Z_HH_Donnan[peptide2], 
            "--g", 
            label='peptide 2 (Theory: HH+Donnan)')

    plt.legend()
    plt.xlabel('pH')
    plt.ylabel('Net charge of the peptide / e')
    plt.show()
    plt.close()

    fig, ax = plt.subplots(figsize=(10, 7))
    plt.errorbar(analyzed_data["pH"]["value"], 
                analyzed_data["mean"]["xi_plus"], 
                yerr=analyzed_data["err_mean"]["xi_plus"], 
                fmt = 'o', 
                capsize=3, 
                label='Simulation')
    ax.plot(pH_range_HH, 
            np.asarray(xi_HH), 
            "-k", 
            label='Theory (HH+Donnan)')
    plt.legend()
    plt.xlabel('pH')
    plt.ylabel(r'partition coefficient $\xi_+$')
    plt.show()
    plt.close()
elif args.mode == "store_HH":
    HH_data=pd.DataFrame({"pH": pH_range_HH,
                 "Z_HH_peptide1": HH_charge_dict["charges_dict"][peptide1],
                 "Z_HH_peptide2": HH_charge_dict["charges_dict"][peptide2],
                 "xi_HH": HH_charge_dict["partition_coefficients"]})
    HH_data.to_csv(args.output / "HH_data.csv", index=False)
