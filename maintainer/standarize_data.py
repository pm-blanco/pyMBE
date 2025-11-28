#
# Copyright (C) 2023-2024 pyMBE-dev team
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

# Import pyMBE and other libraries
import pyMBE
import numpy as np
import pandas as pd
import argparse 
import pathlib

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

# Outputs
output_filenames={"data_landsgesell.csv": "Landsgesell2020a.csv",
                  "Lys-AspMSDE.csv": "Lunkad2021a.csv",
                  "Glu-HisMSDE.csv": "Lunkad2021b.csv",
                  "histatin5_SoftMatter.txt": "Blanco2020a.csv",
                  "1beb-10mM-torres.dat": "Torres2017.csv",
                  "1f6s-10mM-torres.dat": "Torres2022.csv",
                  "landsgesell2022": "Landsgesell2022a.csv"}

parser = argparse.ArgumentParser(description='Script to standarize the data from various authors')
parser.add_argument('--src_filename', 
                    type=str, 
                    required= True,  
                    choices=output_filenames.keys(),
                    help='File name from testsuite/data/src')
args = parser.parse_args()

# Inputs
filename=args.src_filename

files_lunkad=["Glu-HisMSDE.csv","Lys-AspMSDE.csv"]
files_blanco=["histatin5_SoftMatter.txt"]
files_torres = ["1f6s-10mM-torres.dat","1beb-10mM-torres.dat" ]
files_landsgesell2020=["data_landsgesell.csv"]
files_landsgesell2022=["equilibrium_values_gel_MD.csv","weak-gel_total_data.csv"]

data_path = pathlib.Path(__file__).parent.parent / "testsuite" / "data"

# Get path to the data from publication
if filename == "landsgesell2022":
    ref_paths=[]
    for name in files_landsgesell2022:
        ref_paths.append(data_path / "src" / name)
else:
    ref_path=data_path / "src" / filename


if filename in files_lunkad:
    data=pd.read_csv(ref_path)
    Z_ref = 5*-1*data['aaa']+5*data['aab']
    # Error propagation calculation
    # 1/4 factor added to correct for bug in the original calculation of the error reported by the authors       
    Z_ref_err = 5/4*np.sqrt((data['eaa'])**2+(data['eab'])**2)   
    pH_range = np.linspace(2, 12, num=21)
    
elif filename in files_blanco:
    data=np.loadtxt(ref_path, delimiter=",")
    Z_ref=data[:,1]         
    Z_ref_err=data[:,2]
    pH_range = np.linspace(2, 12, num=21)
    
elif filename in files_torres:
    Z_ref = []
    Z_ref_err = []
    pH_range = []
    with open (ref_path,'r') as file:
        for line in file: 
            line_split = line.split ()
            pH = float (line_split[0])
            pH_range.append(pH)
            Z_ref.append (float(line_split[1]))
            Z_ref_err.append(float(line_split[2]))

elif filename in files_landsgesell2020:
    data = pd.read_csv(ref_path, sep="\t", index_col=False)
elif filename == "landsgesell2022":
    dfs = []
    for path in ref_paths:
        dfs.append(pd.read_csv(path))
    data = pd.concat(dfs, axis=0, ignore_index=True)

else:
    raise RuntimeError()

if filename in files_lunkad+files_blanco+files_torres:
    # Create the pandas DataFrame
    data=pd.DataFrame({"pH": pH_range,
                      "charge": Z_ref,
                      "charge_error": Z_ref_err})

# Store the data
data.to_csv(data_path / output_filenames[filename],
            index=False)
