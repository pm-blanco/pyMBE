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

# Load espresso, pyMBE and other necessary libraries
import espressomd
import pandas as pd
import tqdm
from espressomd.io.writer import vtf
from pathlib import Path
import pyMBE
import argparse

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

# Load some functions from the handy_scripts library for convenience
from pyMBE.lib.handy_functions import setup_electrostatic_interactions
from pyMBE.lib.handy_functions import relax_espresso_system
from pyMBE.lib.handy_functions import setup_langevin_dynamics
from pyMBE.lib.handy_functions import do_reaction
from pyMBE.lib.analysis import built_output_name

parser = argparse.ArgumentParser(description='Sample script to run the pre-made peptide models with pyMBE')
parser.add_argument('--sequence',
                    type=str,
                    default= 'EEEEDDDD', # 8 ionizable side-chains
                    help='sequence of the peptide')
parser.add_argument('--pH',
                    type=float,
                    default=7,
                    help='pH of the solution')
parser.add_argument('--output',
                    type=Path,
                    required= False,
                    default=Path(__file__).parent / "time_series" / "peptide_cpH",
                    help='output directory')
parser.add_argument('--test', 
                    default=False, 
                    action='store_true',
                    help='to run a short simulation for testing the script')
parser.add_argument('--no_verbose', action='store_false', help="Switch to deactivate verbosity",default=True)
args = parser.parse_args()

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
frames_path = args.output / "frames"
frames_path.mkdir(parents=True, exist_ok=True)

# Simulation parameters
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)
pH_value = args.pH
N_samples = 1000 # to make the demonstration quick, we set this to a very low value
MD_steps_per_sample = 100 # to make the demonstration quick, we set this to a very low value
N_samples_print = 10  # Write the full trajectory data every X samples
langevin_seed = 100
dt = 0.01
solvent_permitivity = 78.3
verbose=args.no_verbose
ideal=False

if args.test:
    MD_steps_per_sample = 1
    ideal=True

# Peptide parameters
sequence = args.sequence
model = '2beadAA'  # Model with 2 beads per each aminoacid
pep_concentration = 5e-4 *pmb.units.mol/pmb.units.L # mols of peptide chains
N_peptide_chains = 1

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

# Other system parameters, derived from those above
volume = N_peptide_chains/(pmb.N_A*pep_concentration)
L = volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration = N_peptide_chains/(volume*pmb.N_A)

# Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.

path_to_interactions=pmb.root / "parameters" / "peptides" / "Lunkad2021.json"
path_to_pka=pmb.root / "parameters" / "pka_sets" / "Hass2015.json"
pmb.load_interaction_parameters (filename=path_to_interactions) 
pmb.load_pka_set (path_to_pka)

generic_bond_length=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

HARMONIC_parameters = {'r_0'    : generic_bond_length,
                       'k'      : generic_harmonic_constant,
                      }

pmb.define_default_bond(bond_type = 'harmonic',
                        bond_parameters = HARMONIC_parameters)


# Defines the peptide in the pyMBE data frame
peptide_name = 'generic_peptide'
pmb.define_peptide (name=peptide_name, 
                    sequence=sequence, 
                    model=model)

pmb.define_particle(name=cation_name, 
                    z=1, 
                    sigma=0.35*pmb.units.nm, 
                    epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=anion_name,  
                    z=-1, 
                    sigma=0.35*pmb.units.nm,  
                    epsilon=1*pmb.units('reduced_energy'))

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)
espresso_system.time_step=dt
espresso_system.cell_system.skin=0.4
# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
pmb.create_molecule(name=peptide_name, 
                    number_of_molecules=N_peptide_chains,
                    espresso_system=espresso_system, 
                    use_default_bond=True)
# Create counterions for the peptide chains
pmb.create_counterions(object_name=peptide_name,
                        cation_name=cation_name,
                        anion_name=anion_name,
                        espresso_system=espresso_system) 

# check what is the actual salt concentration in the box
# if the number of salt ions is a small integer, then the actual and desired salt concentration may significantly differ
c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,
                                        cation_name=cation_name,
                                        anion_name=anion_name,
                                        c_salt=c_salt)

with open(frames_path / "trajectory0.vtf", mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

#List of ionisable groups
basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisable_groups = basic_groups + acidic_groups
total_ionisable_groups = len(list_ionisable_groups)

if verbose:
    print(f"The box length of your system is {L.to('reduced_length')} {L.to('nm')}")
    print(f"The peptide concentration in your system is {calculated_peptide_concentration.to('mol/L')} with {N_peptide_chains} peptides")
    print(f"The ionisable groups in your peptide are {list_ionisable_groups}")

cpH, labels = pmb.setup_cpH(counter_ion=cation_name, constant_pH=pH_value)
if verbose:
    print(f"The acid-base reaction has been successfully setup for {labels}")

# Setup espresso to track the ionization of the acid/basic groups in peptide
type_map =pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
cpH.set_non_interacting_type (type=non_interacting_type)
if verbose:
    print('The non-interacting type is set to ', non_interacting_type)
if not ideal:
    ##Setup the potential energy
    if verbose:
        print('Setup LJ interaction (this can take a few seconds)')
    pmb.setup_lj_interactions (espresso_system=espresso_system)
    if verbose:
        print('Minimize energy before adding electrostatics')
    relax_espresso_system(espresso_system=espresso_system,
                          seed=langevin_seed)
    if verbose:
        print('Setup and tune electrostatics (this can take a few seconds)')
    setup_electrostatic_interactions(units=pmb.units,
                                    espresso_system=espresso_system,
                                    kT=pmb.kT,
                                    verbose=verbose)
    if verbose:
        print('Minimize energy after adding electrostatics')
    relax_espresso_system(espresso_system=espresso_system,
                          seed=langevin_seed)

if verbose:
    print('Setup Langevin dynamics')
setup_langevin_dynamics(espresso_system=espresso_system, 
                        kT = pmb.kT, 
                        seed = langevin_seed,
                        time_step=dt,
                        tune_skin=False)
# for this example, we use a hard-coded skin value; In general it should be optimized by tuning
espresso_system.cell_system.skin=0.4

#Save the pyMBE dataframe in a CSV file
pmb.write_pmb_df(filename='df.csv')

# Initialize the time series with arbitrary values at time = 0
time_series={} # for convenience, here we save the whole time series in a python dictionary
time_series["time"] = [0.0] 
time_series["charge"] = [0.0] 
# In production, one should save the time series to a file, then read and analyze it after the simulation

# Main loop for performing simulations at different pH-values
N_frame=0
for sample in tqdm.trange(N_samples):

    # LD sampling of the configuration space
    espresso_system.integrator.run(steps=MD_steps_per_sample)        
    # cpH sampling of the reaction space
    do_reaction(cpH, steps=total_ionisable_groups) # rule of thumb: one reaction step per titratable group (on average)
    
    # Get peptide net charge
    charge_dict=pmb.calculate_net_charge(espresso_system=espresso_system, 
                                            molecule_name=peptide_name,
                                            dimensionless=True)
    time_series["time"].append(espresso_system.time)
    time_series["charge"].append(charge_dict["mean"])
    if sample % N_samples_print == 0:
        N_frame+=1
        with open(frames_path / f"trajectory{N_frame}.vtf", mode='w+t') as coordinates:
            vtf.writevsf(espresso_system, coordinates)
            vtf.writevcf(espresso_system, coordinates)
   
# Store time series

data_path=args.output
data_path.mkdir(parents=True, exist_ok=True)
time_series=pd.DataFrame(time_series)
filename=built_output_name(input_dict={"sequence":sequence,"pH":pH_value})

time_series.to_csv(data_path / f"{filename}_time_series.csv", index=False)


