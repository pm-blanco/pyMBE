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
import numpy as np
import pyMBE.lib.np_aux as np_aux

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
                    default=Path(__file__).parent / "time_series" / "nanoparticle_mixture_grxmc_ideal",
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
vol_frac_of_nanoparticles = 0.1		# Volume fraction of the nanoparticle
number_of_nanoparticles   = 10      # Total number of the nanoparticles
nanoparticle_diameter     = 4		# Diameter of the nanoparticle in reduced units
surface_denstity_of_sites = 0.2  	# Surface density of sites in sites/reduced units^2
pka_A_site                = 4.0
pka_B_site                = 10.0

# Names for the componentes of the nanoparticles

core_particle = "core_particle"
A_site        = "A_site"
B_site        = "B_site"

# Patchy distribution of sites A and B

sites_distribution = {"main"     : {"particle_name"     : A_site,
                                    "fraction"          : 0.5,
                                    "number_of_patches" : 2},
                      "secondary": {"particle_name"     : B_site}}

# Short simulation setup for testing

if args.test: 
    MD_steps_per_sample = 1
    phi_np              = 0.1
    np_diameter         = 4
    surf_den_sites      = 0.2	


# Defines the components of the nanoparticle (core particle, A and B type of sites) in the pyMBE data frame

pmb.define_particle(
    name    = core_particle,
    z       = 0,
    sigma   = nanoparticle_diameter*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))

pmb.define_particle(
    name    = A_site,
    acidity = "acidic",
    pka     = pka_A_site,
    sigma   = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))

pmb.define_particle(
    name    = B_site,
    acidity = "basic",
    pka     = pka_B_site,
    sigma   = 1*pmb.units('reduced_length'),
    epsilon = 1*pmb.units('reduced_energy'))

nanoparticle_name = "nanoparticle"
pmb.define_nanoparticle(name                     = nanoparticle_name,
                        core_particle_name       = core_particle,
	                    surface_density_of_sites = surface_denstity_of_sites*pmb.units('reduced_length^-2'),
                        primary_site_particle_name = A_site,
                        fraction_primary_sites = sites_distribution["main"]["fraction"],
                        number_of_patches_of_primary_sites = sites_distribution["main"]["number_of_patches"],
                        secondary_site_particle_name = B_site)

# Saline solution parameters

c_salt = 5e-3 * pmb.units.mol/ pmb.units.L

if args.mode == 'standard':
    proton_name    = 'Hplus'
    hydroxide_name = 'OHminus'
    sodium_name    = 'Na'
    chloride_name  = 'Cl'

    pmb.define_particle(name    = proton_name, 
                        z       = 1, 
                        sigma   = 0.35*pmb.units.nm, 
                        epsilon = 1*pmb.units('reduced_energy'))
    pmb.define_particle(name    = hydroxide_name,  
                        z       = -1,  
                        sigma   = 0.35*pmb.units.nm,  
                        epsilon = 1*pmb.units('reduced_energy'))
    pmb.define_particle(name    = sodium_name, 
                        z       = 1, 
                        sigma   = 0.35*pmb.units.nm, 
                        epsilon = 1*pmb.units('reduced_energy'))
    pmb.define_particle(name    = chloride_name,  
                        z       = -1, 
                        sigma   = 0.35*pmb.units.nm,  
                        epsilon = 1*pmb.units('reduced_energy'))

elif args.mode == 'unified':
    cation_name = 'Na'
    anion_name  = 'Cl'

    pmb.define_particle(name    = cation_name, 
                        z       = 1, 
                        sigma   = 0.35*pmb.units.nm, 
                        epsilon = 1*pmb.units('reduced_energy'))
    pmb.define_particle(name    = anion_name,  
                        z       = -1, 
                        sigma   = 0.35*pmb.units.nm,  
                        epsilon = 1*pmb.units('reduced_energy'))

# System parameters

nanoparticle_tpl      = pmb.db.get_template(name=nanoparticle_name, pmb_type="nanoparticle")
properties            = nanoparticle_tpl.calculate_nanoparticle_properties(pmb)
nanoparticle_volume   = properties["nanoparticle_volume"].to(pmb.units('reduced_length**3'))
volume                = number_of_nanoparticles * nanoparticle_volume / vol_frac_of_nanoparticles
L                     = volume ** (1./3.) # Length of the simulation box

# Create an instance of an espresso system

espresso_system = espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

# Create main and secondary patches
 
def create_patches(nanoparticle_template, nanoparticle_radius, tolerance=1e-6, angle_between_patches = 180):
        """
        Creates a list with the `number_main_patches` lists of the main sites positions and the list of positions for the secondary sites. The coordinates are generated using [0,0,0] as the center of the nanoparticle. Prints a final report with the len of each patch. Calculates the distance between sites and the dipole and quadrupole moments. 

        Args:
            nanoparticle_template(`pint.Quantity`): Template for the nanoparticle.
            nanoparticle_radius(`pint.Quantity`): Radius of the nanoparticle expressed in reduced units.
            tolerance(`float`): Set the tolerance of the numerical method for distributing `total_number_of_sites` in a sphere with radius `nanoparticle_radius`. Defaults = 1e-6.
            angle_between_patches(`float`): If `number_main_patches` = 2, this parameter corresponds to the angle between the vectors formed from the center of the nanoparticle and the center of the two patches. Defaults = 180 (the patches are in the poles). If `number_main_patches` > 2, the patches are distributed approximately uniformly over the surface. 
        
        Returns:
            sites_positions_per_patch(`list` of `list`): List with the list of the positions of the main and secondary sites in the form: [[sites_positions_main_patch_1],[sites_positions_main_patch_2],...,[sites_positions_secondary_patch]].
            distance_between_sites(`dict`): Dictionary containing the average, standard deviation and standard error of the distances between sites on the surface of the nanoparticle.  
            dipole_moment(`dict`): Dictionary containing the vector and magnitude of the dipole moment.
            quadrupole_moment(`dict`): Dictionary containing the vector, magnitude and eigenvalues of the quadrupole moment.

        Note:
        The sites are created 1/2 of the reduced unit inside of the nanoparticle surface to avoid overlapping of charges due to electrostatic attractions in abcense of excluded volume. 
        """
        main_site_name              = nanoparticle_template.primary_site_particle_name
        secondary_site_name         = nanoparticle_template.secondary_site_particle_name
        total_number_of_sites       = properties["total_number_of_sites"]
        number_main_patches         = nanoparticle_template.number_of_patches_of_primary_sites
        number_main_sites_per_patch = properties["number_of_primary_sites_per_patch"]
        number_main_sites           = properties["number_of_primary_sites"]
        number_secondary_sites      = properties["number_of_secondary_sites"]

        radius_sites                = (nanoparticle_radius - (1/2)*pmb.units('reduced_length')).magnitude 
        root_edges                  = np_aux.uniform_distribution_sites_on_sphere(number_of_edges = total_number_of_sites, tolerance=tolerance)
        nanoparticle_edges          = np.multiply(root_edges, radius_sites)
        
        if number_main_patches <= 2:
            initial_edge                = [nanoparticle_edges[0]]
            sites_positions_per_patch   = []
            distances_to_center_site_patch_1, sites_positions_main_patch_1 = np_aux.define_patch(points = nanoparticle_edges, central_point = initial_edge[0], patch_size    = number_main_sites_per_patch)
            sites_positions_per_patch.append(sites_positions_main_patch_1)
            
            if number_main_patches == 2:
                distance_omega            = (2 * radius_sites**2 - 2 * radius_sites**2 *np.cos(np.radians(angle_between_patches)))**(1/2)
                distance_omega_to_patch_1 = list(np.abs(distances_to_center_site_patch_1 - distance_omega))
                initial_edge.append(nanoparticle_edges[distance_omega_to_patch_1.index(min(distance_omega_to_patch_1))])
                distances_to_center_site_patch_2, sites_positions_main_patch_2 = np_aux.define_patch(points = nanoparticle_edges, central_point = initial_edge[1], patch_size    = number_main_sites_per_patch)
                sites_positions_per_patch.append(sites_positions_main_patch_2)
                
                # Checking that there is no overlapping between patches
                np_aux.check_patch_overlaps(sites_positions=sites_positions_per_patch,number_patches=number_main_patches)
        
        elif number_main_patches > 2:
            sites_positions_per_patch  = []
            initial_patch_edges        = np_aux.uniform_distribution_sites_on_sphere(number_of_edges = number_main_patches, tolerance=tolerance)
            initial_patch_scaled_edges = np.multiply(initial_patch_edges, radius_sites)
            initial_edge = []
            for scaled_edge in initial_patch_scaled_edges:
                comparison_patch_to_nanoparticle = np_aux.calculate_distance_vector_point(nanoparticle_edges,scaled_edge)
                initial_edge.append(nanoparticle_edges[comparison_patch_to_nanoparticle.index(min(comparison_patch_to_nanoparticle))])
            for main_patch in range(number_main_patches):
                distances_to_center_site_patch, sites_positions_main_patch = np_aux.define_patch(points = nanoparticle_edges, central_point = initial_edge[main_patch], patch_size    = number_main_sites_per_patch)
                sites_positions_per_patch.append(sites_positions_main_patch)

            # Checking that there is no overlapping between patches
            np_aux.check_patch_overlaps(sites_positions=sites_positions_per_patch,number_patches=number_main_patches)
        
        else:
            sites_positions_per_patch = []
        
        # Creating the secondary patch with the remaining sites
        
        sites_positions_secondary_patch = nanoparticle_edges
        for site_position in np.vstack(sites_positions_per_patch):
            sites_positions_secondary_patch = set(map(tuple, sites_positions_secondary_patch)).difference({tuple(site_position)})
        sites_positions_per_patch.append(list(sites_positions_secondary_patch))
	
        # Final report

        counted_total_number_of_sites = 0
        for patch_index, patch in enumerate(sites_positions_per_patch):
            if patch_index < (len(sites_positions_per_patch)-1):
                print('Sites in main patch ',patch_index+1,'     : ', len(patch))
                counted_total_number_of_sites += len(patch)
            else: 
                print('Sites in secondary patch    : ', len(patch))
                counted_total_number_of_sites += len(patch)

        print('Total sites                 : ', counted_total_number_of_sites)
        
        # Calculating the distances between sites        

        avg_distance_between_sites, standard_deviation, standard_error = np_aux.calculate_distance_between_points_on_sphere(points=sites_positions_per_patch)        
        print('Mean spacing between sites  : ', avg_distance_between_sites)
        print('Standard deviation          : ', standard_deviation,)
        print('Standard error              : ', standard_error)
        distance_between_sites = {"average_distance_between_sites" : avg_distance_between_sites,
                                  "standard_deviation"             : standard_deviation,
                                  "standard_error"                 : standard_error}

        # Calculating the dipole and quadrupole moments

        charges                      = pmb.get_charge_number_map()        
        types                        = pmb.get_type_map()
        main_sites_charge            = charges[types[main_site_name]]
        secondary_site_charges       = charges[types[secondary_site_name]]
        positions_map                = np.vstack(sites_positions_per_patch)
        charges_map                  = np.concatenate((np.ones(number_main_sites)*main_sites_charge,np.ones(number_secondary_sites)*secondary_site_charges))
        dipole_vector, dipole_magnitude = np_aux.calculate_dipole_moment(charges_map, positions_map)
        print('Dipole moment magnitude     : ', dipole_magnitude,'in e * reduced length')
        dipole_moment = {"dipole_vector"    : dipole_vector,
                         "dipole_magnitude" : dipole_magnitude}
        quadrupole_matrix, quadrupole_magnitude, quadrupole_eigenvalues = np_aux.calculate_quadrupole_moment(charges_map, positions_map)
        print('Quadrupole moment magnitude : ', quadrupole_magnitude, 'in e * reduced length**2' )
        quadrupole_moment = {"quadrupole_matrix"      : quadrupole_matrix,
                             "quadrupole_magnitude"   : quadrupole_magnitude,
                             "quadrupole_eigenvalues" : quadrupole_eigenvalues}

        return sites_positions_per_patch, distance_between_sites, dipole_moment, quadrupole_moment

# Create nanoparticles

def create_nanoparticle(name, espresso_system, number_of_nanoparticles, list_core_particle_positions=None, fix=False):
        """
        Creates `number_of_nanoparticles` nanoparticles of type `name` into `espresso_system` and bookkeeps them into `pymbe.df`.
        
        Args:
            name(`str`): Label of the nanoparticle type to be created. `name` must be a `nanoparticle` defined in `pmb_df`.  
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            number_of_nanoparticles(`int`): Number of nanoparticles to be created.
            list_core_particle_positions(list of [`float`,`float`,`float`], optional): Initial positions of the nanoparticles' core. If not given, core particles are created in random positions. Defaults to None.
            fix(`bool`, optional): Controls if the nanoparticle motion is frozen in the integrator, it is used to create rigid objects. Defaults to False.
        Returns:
            nanoparticles_info(`dict`): Dictionary with the ids of the nanoparticle cores and their corresponding sites created into `espresso_system`.
        """ 

        if number_of_nanoparticles <=0:
            return []
        
        # Get information from the nanoparticle type `name` from the df
        nanoparticle_tpl = pmb.db.get_template(name=nanoparticle_name, pmb_type="nanoparticle")
        properties = nanoparticle_tpl.calculate_nanoparticle_properties(pmb)

        radius                      = pmb.get_radius_map()
        types                       = pmb.get_type_map()
        core_particle_name          = nanoparticle_tpl.core_particle_name
        nanoparticle_radius         = radius[types[core_particle_name]]*pmb.units('reduced_length')
        main_site_name              = nanoparticle_tpl.primary_site_particle_name
        secondary_site_name         = nanoparticle_tpl.secondary_site_particle_name
        total_number_of_sites       = properties["total_number_of_sites"]
        number_main_patches         = nanoparticle_tpl.number_of_patches_of_primary_sites
        number_main_sites_per_patch = properties["number_of_primary_sites_per_patch"]
        number_main_sites           = properties["number_of_primary_sites"]
        number_secondary_sites      = properties["number_of_secondary_sites"]
        sites_types                 = [main_site_name]*number_main_patches + [secondary_site_name]*1	
        number_sites_per_type       = [number_main_sites_per_patch]*number_main_patches + [number_secondary_sites]*1
        sites_positions_per_patch, distance_between_sites, dipole_magnitude, quadrupole_moment = create_patches(
                                      nanoparticle_template = nanoparticle_tpl,
                                      nanoparticle_radius   = nanoparticle_radius, 
                                      angle_between_patches = 180)
        print(distance_between_sites)
        print(dipole_magnitude)
        print(quadrupole_moment)
        
        nanoparticles_info      = {}
        for nanoparticle_index in range(number_of_nanoparticles):
            # create the principal bead
            if not list_core_particle_positions:
                core_particle_id = pmb.create_particle(name                 = core_particle_name,
                                                        espresso_system     = espresso_system,
                                                        number_of_particles = 1)[0]
            else:
                core_particle_id = pmb.create_particle(name                 = core_particle_name,
                                                        espresso_system     = espresso_system,
                                                        position            = [list_core_particle_positions[nanoparticle_index]],
                                                        number_of_particles = 1)[0]
            core_particle_position = espresso_system.part.by_id(core_particle_id).pos
            
            # Internal bookkeeping of the core_particle_id
            nanoparticles_info[nanoparticle_index] = {}
            nanoparticles_info[nanoparticle_index]['core_particle_id'] = core_particle_id 
            
            # create the main and secondary sites
            sites_ids = []
            for index_patch, sites_type in enumerate(sites_types):
                number_sites    = number_sites_per_type[index_patch]
                sites_positions = sites_positions_per_patch[index_patch]
                
                sites_id = pmb.create_particle(name                = sites_type,
                                               espresso_system     = espresso_system,
                                               position            = sites_positions,
                                               number_of_particles = number_sites)
                sites_ids.append(sites_id)
            nanoparticles_info[nanoparticle_index]['sites_ids'] = sites_ids

        espresso_system.integrator.set_vv()
        espresso_system.virtual_sites  = espressomd.virtual_sites.VirtualSitesRelative()
        espresso_system.min_global_cut = (nanoparticle_radius - (1/2)*pmb.units('reduced_length')).magnitude
        
        print(nanoparticles_info)
        for nanoparticle_index, nanoparticle_info in enumerate(nanoparticles_info):
            com     = []
            print(nanoparticle_info)
            for patch in nanoparticle_info['sites_ids']:
                print(patch)
                exit()
                com.append(np.average(patch.pos, 0))
            com.append(np.average(surf_parts_per_np_base[i].pos, 0))
            com_avg = sum(com)/(N_patchy_acid+N_patchy_base)  
            
            momI    = 0
            for patchy in range(N_patchy_acid):
                for par in surf_parts_per_np_acid[patchy][i]:
                    momI += np.power(np.linalg.norm(com_avg - par.pos), 2)
            for par in surf_parts_per_np_base[i]:
                momI += np.power(np.linalg.norm(com_avg - par.pos), 2)
            
            central_parts.fix[i] = [False, False, False]
            central_parts.pos[i] = com_avg
            central_parts.mass[i] = N_s                   # ASK
            central_parts.rinertia[i] = np.ones(3) * momI
            
            for patchy in range(N_patchy_acid):
                for par in surf_parts_per_np_acid[patchy][i]:
                    par.vs_auto_relate_to(central_parts.id[i])
            for par in surf_parts_per_np_base[i]:
                par.vs_auto_relate_to(central_parts.id[i])
            
        print("Relaxation of the raspberry surface particles done")

        return nanoparticles_info
        
nanoparticle_info = create_nanoparticle(name=nanoparticle_name, espresso_system=espresso_system, number_of_nanoparticles=number_of_nanoparticles,list_core_particle_positions=None)
print(nanoparticle_info)
exit()

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
