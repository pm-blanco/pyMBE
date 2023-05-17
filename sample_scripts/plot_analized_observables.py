from matplotlib import colors
from matplotlib.axis import YAxis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse 
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import StrMethodFormatter

import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import pyMBE
pmb = pyMBE.pymbe_library()

filename = '../reference_data/alac-10mM-torres2022.dat'

# filename_fortran = 'results/fortrancode/analyzed_observables_fortran_code.csv'

filename_espresso = 'analyzed_observables.csv'

protein_name = '1f6s'
protein_filename = os.path.join (parentdir,'sample_scripts/coarse_grain_model_of_1f6s.vtf' )
pmb.load_protein_vtf_in_df (name=protein_name,filename=protein_filename)
protein_sequence = pmb.df.loc[pmb.df['name']== protein_name].sequence.values[0]

pmb.load_pka_set (filename=os.path.join(parentdir,'reference_parameters/pka_sets/Nozaki1967.txt'))

pH_range = np.linspace(2, 7, num=31)
Z_HH = pmb.calculate_HH (sequence = protein_sequence, pH_list =pH_range )
Z_HH_Ca = []
# Here we have to add +2 for the Calcium in the protein 

for value in Z_HH:
    new_value= value +2 
    Z_HH_Ca.append (new_value)

pH_list = []
znet = []
sigma = []

pH_ideal = []
zideal = []

with open (filename,'r') as file:
    for line in file: 
        line_split = line.split ()
        pH = float (line_split[0])
        pH_list.append(pH)
        znet.append (float(line_split[1]))
        sigma.append(float(line_split[2]))


# fortran_df = pd.read_csv(filename_fortran)
# cols = fortran_df.columns.to_list ()


full_data = pd.read_csv(filename_espresso)
columns = full_data.columns.to_list()
full_data.sort_values(by=['pH'],inplace=True)

fig = plt.figure(figsize = (20,20))
ax1 = fig.add_subplot (111)

ax1.hlines(y=0,xmin=2,xmax=7,lw = 3, colors='grey',linestyles='dashed')

ax1.plot(
    pH_range,
    Z_HH_Ca,
    linewidth = 3,
    #marker = 'o',
    #markersize = 10,
    # markerfacecolor = 'none',
    #alpha = 0.8,
    color = 'black',
    label = 'Ideal' 
    )

ax1.errorbar (
    pH_list,
    znet,
    yerr=sigma,
    linewidth = 2,
    elinewidth = 3,
    marker = 'o',
    markersize = 10,
    markerfacecolor = 'none',
    alpha = 0.8, #changes the line opacity
    color = 'red',
    label = 'Simulations Data Arg ')

ax1.errorbar (
    full_data['pH'],
    full_data['Znet'],
    yerr = full_data['ZnetErr'],
    linewidth = 2,
    elinewidth = 3,
    marker = 's',
    markersize = 10,
    markerfacecolor = 'none',
    alpha = 0.8, #changes the line opacity
    color = 'blue',
    label = 'ESPResSo')
print (full_data['pH'], full_data['Znet'])


# ax1.plot(
#     pH_list,
#     znet,
#     linewidth = 0,
#     marker = 'o',
#     markersize = 10,
#     # markerfacecolor = 'none',
#     alpha = 0.8, #changes the line opacity
#     color = 'red',
#     label = 'MC Sim Arg'
    
#     )



ax1.set_title ('Net charge vs pH. $c_{salt}$ = 0.01 M',fontsize ='40',pad = 30)
ax1.set_xlabel('$\it{pH}$',size = 45)
ax1.set_ylabel('$\it{Z}$', size = 45)

ax1.invert_yaxis ()
ax1.set_xlim ([2.0,7.0])
ax1.set_ylim ([-10.0,20.0])

y = np.linspace(-10,20,7)
minor_ticks = np.linspace(-10,20,31)

# pH_list = fortran_df ['pH'].tolist()

ax1.set_xticks(pH_list, minor=True)
ax1.set_yticks (y,minor=True)
ax1.set_yticks (minor_ticks,minor=True)

# ax1.set_xticks(5)

ax1.tick_params(axis='x',which='major',length=30,direction='in',width=3,colors='black',pad=15,labelsize=40)
ax1.tick_params(axis='x',which='minor',length=15,direction='in',width=3)

ax1.tick_params(axis='y',which='major',length=30,direction='in',width=3,colors='black',pad=10,labelsize=40) 
ax1.tick_params(axis='y',which='minor',length=15,direction='in',width=3,colors='black')

# # ax2.tick_params(axis='y',which='major',length=20,direction='in',width=3,colors='black',pad=10,labelsize=30)
# # ax2.tick_params(axis='y',which='minor',length=10,direction='in',width=3)

# ax1.yaxis.label.set_color('black')
ax1.spines['left'].set_color('black')
ax1.spines['left'].set_lw(3)
ax1.spines['top'].set_lw(3)
ax1.spines['right'].set_lw(3)
ax1.spines['bottom'].set_lw(3)

ax1.legend(frameon =False)

# 
# ax1.set_xscale('log')
# ax1.set_yscale('log')


plt.legend(prop={'size': 35})
pdf_name = 'analyzed_observables.pdf'
plt.savefig(pdf_name)
plt.show()


