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

#Load espresso, pyMBE and other necessary libraries
import espressomd
from pathlib import Path
import pandas as pd
import argparse
from espressomd.io.writer import vtf
import pyMBE
from pyMBE.lib.analysis import built_output_name
from pyMBE.lib.handy_functions import do_reaction

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

# Command line arguments

parser = argparse.ArgumentParser(description='Script that runs a simulation of an ideal peptide mixture in the grand-reaction ensemble using pyMBE and ESPResSo.')
parser.add_argument('--mode',
                    type=str,
                    default= "standard",
                    choices=["standard", "unified"],
                    help='Set if the grand-reaction method is used with unified ions or not')
parser.add_argument('--test', 
                    default=False, 
                    action='store_true',
                    help='to run a short simulation for testing the script')
parser.add_argument('--sequence1',
                    type=str,
                    default= 'nEHc', 
                    help='sequence of the first peptide')
parser.add_argument('--sequence2',
                    type=str,
                    default= 'nEEHHc', 
                    help='sequence of the second peptide')
parser.add_argument('--pH',
                    type=float,
                    default=7,
                    help='pH of the solution')
parser.add_argument('--output',
                    type=Path,
                    required= False,
                    default=Path(__file__).parent / "time_series" / "peptide_mixture_grxmc_ideal",
                    help='output directory')
parser.add_argument('--no_verbose', action='store_false', help="Switch to deactivate verbose",default=True)
args = parser.parse_args()

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
frames_path = args.output / "frames"
frames_path.mkdir(parents=True, exist_ok=True)

# Simulation parameters
verbose=args.no_verbose
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm, 
                      Kw=1e-14)
N_samples = 1000 # to make the demonstration quick, we set this to a very low value
MD_steps_per_sample = 1000
N_samples_print = 1000  # Write the trajectory every 100 samples
LANGEVIN_SEED = 42
dt = 0.001
solvent_permitivity = 78.3
pH_value= args.pH
# Peptide parameters
sequence1 = args.sequence1
model = '1beadAA'  # Model with 2 beads per each aminoacid
pep1_concentration = 1e-2 *pmb.units.mol/pmb.units.L
N_peptide1_chains = 10

sequence2 = args.sequence2
pep2_concentration = 1e-2 *pmb.units.mol/pmb.units.L
N_peptide2_chains = 10

if args.test:
    MD_steps_per_sample = 1
    N_peptide1_chains = 5
    N_peptide2_chains = 5
    
# Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.
# Note that this parameterization only includes some of the natural aminoacids
# For the other aminoacids the user needs to use  a parametrization including all the aminoacids in the peptide sequence
path_to_pka=pmb.root / "parameters" / "pka_sets" / "Hass2015.json"
path_to_interactions=pmb.root / "parameters" / "peptides" / "Lunkad2021.json"

pmb.load_interaction_parameters(filename=path_to_interactions) 
pmb.load_pka_set(path_to_pka)

# Defines the bonds
bond_type = 'harmonic'
generic_bond_length=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

harmonic_bond = {'r_0'    : generic_bond_length,
                 'k'      : generic_harmonic_constant}

pmb.define_default_bond(bond_type = bond_type, 
                        bond_parameters = harmonic_bond)

# Defines the peptides in the pyMBE data frame
peptide1 = 'generic_peptide1'
pmb.define_peptide (name=peptide1, 
                    sequence=sequence1, 
                    model=model)
peptide2 = 'generic_peptide2'
pmb.define_peptide (name=peptide2, 
                    sequence=sequence2, 
                    model=model)

# Solution parameters
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

if args.mode == 'standard':
    proton_name = 'Hplus'
    hydroxide_name = 'OHminus'
    sodium_name = 'Na'
    chloride_name = 'Cl'

    pmb.define_particle(name=proton_name, 
                        z=1, 
                        sigma=0.35*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))
    pmb.define_particle(name=hydroxide_name,  
                        z=-1, 
                        sigma=0.35*pmb.units.nm,  
                        epsilon=1*pmb.units('reduced_energy'))
    pmb.define_particle(name=sodium_name, 
                        z=1, 
                        sigma=0.35*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))
    pmb.define_particle(name=chloride_name,  
                        z=-1, 
                        sigma=0.35*pmb.units.nm,  
                        epsilon=1*pmb.units('reduced_energy'))

elif args.mode == 'unified':
    cation_name = 'Na'
    anion_name = 'Cl'

    pmb.define_particle(name=cation_name, 
                        z=1, 
                        sigma=0.35*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))
    pmb.define_particle(name=anion_name,  
                        z=-1, 
                        sigma=0.35*pmb.units.nm,  
                        epsilon=1*pmb.units('reduced_energy'))


# System parameters
volume = N_peptide1_chains/(pmb.N_A*pep1_concentration)
L = volume ** (1./3.) # Side of the simulation box
calculated_peptide_concentration = N_peptide1_chains/(volume*pmb.N_A)

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
pmb.create_molecule(name=peptide1, 
                    number_of_molecules=N_peptide1_chains,
                    espresso_system=espresso_system, 
                    use_default_bond=True)
pmb.create_molecule(name=peptide2, 
                    number_of_molecules=N_peptide2_chains,
                    espresso_system=espresso_system, 
                    use_default_bond=True)

if args.mode == 'standard':
    pmb.create_counterions(object_name=peptide1,
                           cation_name=proton_name,
                           anion_name=hydroxide_name,
                           espresso_system=espresso_system) # Create counterions for the peptide chains with sequence 1
    pmb.create_counterions(object_name=peptide2,
                           cation_name=proton_name,
                           anion_name=hydroxide_name,
                           espresso_system=espresso_system) # Create counterions for the peptide chains with sequence 2

    c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,
                                              cation_name=sodium_name,
                                              anion_name=chloride_name,
                                              c_salt=c_salt)
elif args.mode == 'unified':
    pmb.create_counterions(object_name=peptide1, 
                           cation_name=cation_name,
                           anion_name=anion_name,
                           espresso_system=espresso_system) # Create counterions for the peptide chains with sequence 1
    pmb.create_counterions(object_name=peptide2, 
                           cation_name=cation_name,
                           anion_name=anion_name,
                           espresso_system=espresso_system) # Create counterions for the peptide chains with sequence 2

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
total_ionisable_groups = len (list_ionisable_groups)
# Get peptide net charge
if verbose:
    print("The box length of your system is", L.to('reduced_length'), L.to('nm'))

if args.mode == 'standard':
    grxmc, sucessful_reactions_labels, ionic_strength_res = pmb.setup_grxmc_reactions(pH_res=pH_value, 
                                                                                   c_salt_res=c_salt, 
                                                                                   proton_name=proton_name, 
                                                                                   hydroxide_name=hydroxide_name, 
                                                                                   salt_cation_name=sodium_name, 
                                                                                   salt_anion_name=chloride_name,
                                                                                   activity_coefficient=lambda x: 1.0)
elif args.mode == 'unified':
    grxmc, sucessful_reactions_labels, ionic_strength_res = pmb.setup_grxmc_unified(pH_res=pH_value, 
                                                                                 c_salt_res=c_salt, 
                                                                                 cation_name=cation_name, 
                                                                                 anion_name=anion_name,
                                                                                 activity_coefficient=lambda x: 1.0)
if verbose:
    print('The acid-base reaction has been sucessfully setup for ', sucessful_reactions_labels)

# Setup espresso to track the ionization of the acid/basic groups in peptide
type_map =pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map(type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
grxmc.set_non_interacting_type (type=non_interacting_type)
if verbose:
    print('The non interacting type is set to ', non_interacting_type)

espresso_system.time_step = dt

#Save the initial state
with open(frames_path / "trajectory1.vtf", mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)

# Setup espresso to do langevin dynamics
espresso_system.time_step= dt 
espresso_system.integrator.set_vv()
espresso_system.thermostat.set_langevin(kT=pmb.kT.to('reduced_energy').magnitude, gamma=0.1, seed=LANGEVIN_SEED)
espresso_system.cell_system.skin=0.4

#Save the pyMBE dataframe in a CSV file
pmb.write_pmb_df (filename='df.csv')
time_series={}
for label in ["time","charge_peptide1","charge_peptide2","num_plus","xi_plus"]:
    time_series[label]=[] 

# Main simulation loop
N_frame=0
for step in range(N_samples):
    espresso_system.integrator.run(steps=MD_steps_per_sample)        
    do_reaction(grxmc, steps=total_ionisable_groups)
    time_series["time"].append(espresso_system.time)
    # Get net charge of peptide1 and peptide2
    charge_dict_peptide1=pmb.calculate_net_charge(espresso_system=espresso_system, 
                                            molecule_name=peptide1,
                                            dimensionless=True)
    charge_dict_peptide2=pmb.calculate_net_charge(espresso_system=espresso_system, 
                                            molecule_name=peptide2,
                                            dimensionless=True)
    time_series["charge_peptide1"].append(charge_dict_peptide1["mean"])
    time_series["charge_peptide2"].append(charge_dict_peptide2["mean"])
    if args.mode == 'standard':
        num_plus = espresso_system.number_of_particles(type=type_map["Na"])+espresso_system.number_of_particles(type=type_map["Hplus"])
    elif args.mode == 'unified':
        num_plus = espresso_system.number_of_particles(type=type_map["Na"])      
    time_series["num_plus"].append(num_plus)
    concentration_plus = (num_plus/(pmb.N_A * L**3)).to("mol/L")
    xi_plus = (concentration_plus/ionic_strength_res).magnitude
    time_series["xi_plus"].append(xi_plus)
    if step % N_samples_print == 0:
        N_frame+=1
        with open(frames_path / f"trajectory{N_frame}.vtf", mode='w+t') as coordinates:
            vtf.writevsf(espresso_system, coordinates)
            vtf.writevcf(espresso_system, coordinates)

# Store time series
data_path=args.output
data_path.mkdir(parents=True, exist_ok=True)
time_series=pd.DataFrame(time_series)

filename=built_output_name(input_dict={"mode":args.mode,
                                       "sequence1":sequence1,
                                       "sequence2": sequence2,
                                       "pH":pH_value})

time_series.to_csv(data_path / f"{filename}_time_series.csv",
                    index=False)
