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

#Load espresso, pyMBE and other necessary libraries
import espressomd
from pathlib import Path
import pandas as pd
import argparse
from espressomd.io.writer import vtf
import pyMBE
from pyMBE.lib.analysis import built_output_name
from pyMBE.lib.handy_functions import do_reaction
from pyMBE.storage.df_management import _DFManagement as _DFm
import numpy as np

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
verbose = args.no_verbose
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm, 
                      Kw=1e-14)
N_samples           = 1000	# to make the demonstration quick, we set this to a very low value
MD_steps_per_sample = 1000
N_samples_print     = 1000	# Write the trajectory every 100 samples
LANGEVIN_SEED 	    = 42
dt                  = 0.001
solvent_permitivity = 78.3
pH_value= args.pH

# Nanoparticle parameters
phi_np          = 0.1			# Volume fraction of the nanoparticle
np_diameter     = 4			# Diameter of the nanoparticle in reduced units
surf_den_sites  = 0.2			# Surface density of sites in sites/reduced units^2
pka_acidic_site = 4.0
pka_basic_site  = 10.0

if args.test: 
    MD_steps_per_sample = 1
    phi_np              = 0.1
    np_diameter         = 4
    surf_den_sites      = 0.2	


# Defines the components of the nanoparticle (core particle, acidic and basic sites) in the pyMBE data frame

core_particle = "core_particle"
pmb.define_particle(
    name    = "core_particle",
    z       = 0,
    sigma   = np_diameter*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))

acidic_site = "acidic_site"
pmb.define_particle(
    name    = "acidic_site",
    acidity = "acidic",
    pka     = pka_acidic_site,
    sigma   = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))

basic_site = "basic_site"
pmb.define_particle(
    name    = "basic_site",
    acidity = "basic",
    pka     = pka_basic_site,
    sigma   = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))

# Define nanoparticle

def define_nanoparticle(name, core_particle_name, surface_density_of_sites, sites_distribution):
        """
        Defines a pyMBE object of type `nanoparticle` in `pymbe.df`.

        Args:
            name(`str`): Unique label that identifies the `nanoparticle`.
            core_particle_name(`str`): `name` of the `particle` to be placed as the core_particle of the `nanoparticle`.
            surface_density_of_sites(`pint.Quantity`): surface density of sites on the surface of the core_particle, should be expressed in reduced units^-2. Together with the radius of the core_particle, this parameter is used to calculate the total_number_of_sites. 
            sites_distribution(`dict`): Dictionary containing the distribution of ionizable sites on the surface of the nanoparticle. Currently, nanoparticles can have a maximum of two different kind of sites and only pmb_objects of type `particle` are supported. 
            The dictionary should have the structure: {"main"     : {"particle_name"     = acidic_site,
                                                                     "fraction"          = fraction,
                                                                     "number_of_patches" = n_patches},
                                                       "secondary": {"particle_name"     = basic_site}} 
            Here, the "main" sites are located in n="number_of_patches" approximately circular patches over the surface of the core_particle, while the "secondary" ones correspond to the remaining sites around these patches. The parameter "fraction" determines the proportion of the total_number_of_sites that are classified as "main" sites, consequently, (1-fraction) defines the proportion of "secondary" sites.
            Currently, pyMBE supports nanoparticles with uniformily distributed sites, and the patches are positioned approxiamtely equidistant from one another.  

	"""
        # Sanity checks

        _DFm._check_if_multiple_pmb_types_for_name(name=name,
                                                   pmb_type_to_be_defined='nanoparticle',
                                                   df=pmb.df)
        
        pmb.check_dimensionality(surface_density_of_sites,"[length]**-2")

        radius                      = pmb.get_radius_map()
        types                       = pmb.get_type_map()
        S_sphere                    = 4 * np.pi * radius[types[core_particle_name]]**2
        N_sites                     = int(S_sphere  * surface_density_of_sites.magnitude)
        real_surface_charge_density = N_sites / S_sphere
        number_main_sites           = int(N_sites * sites_distribution['main']['fraction'])
        real_fraction               = number_main_sites/N_sites

        index = len(pmb.df)
        pmb.df.at [index,'name']                  = name
        pmb.df.at [index,'pmb_type']              = 'nanoparticle'
        pmb.df.at [index,'core_particle']         = core_particle_name,
        pmb.df.at [index,'surface_density_sites'] = real_surface_charge_density, 
        pmb.df.at [index,'main_site']             = sites_distribution['main']['particle_name']
        pmb.df.at [index,'fraction_main_site']    = real_fraction
        pmb.df.at [index,'number_main_patches']   = sites_distribution['main']['number_of_patches']
        pmb.df.at [index,'secondary_site']        = sites_distribution['secondary']['particle_name']
        pmb.df.fillna(pd.NA, inplace=True)
        return

define_nanoparticle(    name                     = "nanoparticle",
                        core_particle_name       = core_particle,
			surface_density_of_sites = surf_den_sites*pmb.units('reduced_length^-2'),
                        sites_distribution   = {"main"     : {"particle_name"     : acidic_site,
                                                            "fraction"          : 0.5,
                                                            "number_of_patches" : 1},
                                              "secondary": {"particle_name"     : basic_site}},
)

#Save the pyMBE dataframe in a CSV file
pmb.write_pmb_df (filename='df.csv')

exit()

#def  create_nanoparticle(name=nanoparticle_name)
    # Get info about the NP properties from the df

    # Call specific nanoparticle builder depending on the NP properties

#pmb.create_nanoparticle(name=nanoparticle_name)


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
