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
vol_frac_of_nanoparticles = 0.1		# Volume fraction of the nanoparticle
number_of_nanoparticles   = 10          # Total number of the nanoparticles
nanoparticle_diameter     = 4		# Diameter of the nanoparticle in reduced units
surface_denstity_of_sites = 0.2		# Surface density of sites in sites/reduced units^2
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

# Define nanoparticle

def define_nanoparticle(name, core_particle_name, surface_density_of_sites, sites_distribution):
        """
        Defines a pyMBE object of type `nanoparticle` in `pymbe.df`.

        Args:
            name(`str`): Unique label that identifies the `nanoparticle`.
            core_particle_name(`str`): `name` of the `particle` to be placed as the core_particle of the `nanoparticle`.
            surface_density_of_sites(`pint.Quantity`): surface density of sites on the surface of the core_particle, should be expressed in reduced units^-2. Together with the radius of the core_particle, this parameter is used to calculate the total_number_of_sites. 
            sites_distribution(`dict`): Dictionary containing the distribution of ionizable sites on the surface of the nanoparticle. Currently, nanoparticles can have a maximum of two different kind of sites and only pmb_objects of type `particle` are supported. 
            The dictionary should have the structure: {"main"     : {"particle_name"     = A_site,
                                                                     "fraction"          = fraction,
                                                                     "number_of_patches" = n_patches},
                                                       "secondary": {"particle_name"     = B_site}} 
            Here, the "main" sites are located in n="number_of_patches" approximately circular patches over the surface of the core_particle, while the "secondary" ones correspond to the remaining sites around these patches. The parameter "fraction" determines the proportion of the total_number_of_sites that are classified as "main" sites, consequently, (1-fraction) defines the proportion of "secondary" sites.
            Currently, pyMBE supports nanoparticles with uniformily distributed sites, and the patches are positioned approxiamtely equidistant from one another.  

	"""
        # Sanity checks
	
	# Check if there is an existent pyMBE object using the requested name

        _DFm._check_if_multiple_pmb_types_for_name(name=name,
                                                   pmb_type_to_be_defined='nanoparticle',
                                                   df=pmb.df)
        
	# Check if the dimensionality if the surface density of sites is correct

        pmb.check_dimensionality(surface_density_of_sites,"[length]**-2")

	# Calculatin the total_number_of_sites, number_of_main_sites and number_of_secondary_sites in the nanoparticle surface, and recalculate both the surface density of
	
        radius                        = pmb.get_radius_map()
        types                         = pmb.get_type_map()
        nanoparticle_surface_area     = 4 * np.pi * (radius[types[core_particle_name]]*pmb.units('reduced_length'))**2
        nanoparticle_volume           = 4 / 3 * np.pi * (radius[types[core_particle_name]]*pmb.units('reduced_length'))**3
        total_number_of_sites         = int(nanoparticle_surface_area  * surface_density_of_sites)
        print(total_number_of_sites)
        real_surface_density_of_sites = total_number_of_sites / nanoparticle_surface_area
        number_main_sites             = int(total_number_of_sites * sites_distribution['main']['fraction'])
        number_secondary_sites        = total_number_of_sites - number_main_sites
        real_fraction                 = number_main_sites/total_number_of_sites

        index = len(pmb.df)
        pmb.df.at [index,'name']                      = name
        pmb.df.at [index,'pmb_type']                  = 'nanoparticle'
        pmb.df.at [index,'core_particle']             = core_particle_name,
        pmb.df.at [index,'nanoparticle_surface_area'] = nanoparticle_surface_area.magnitude,     #ASK
        pmb.df.at [index,'nanoparticle_volume']       = nanoparticle_volume.magnitude,           #ASK
        pmb.df.at [index,'surface_density_sites']     = real_surface_density_of_sites.magnitude, #ASK
        pmb.df.at [index,'total_number_of_sites']     = total_number_of_sites,
        pmb.df.at [index,'main_site']                 = sites_distribution['main']['particle_name']
        pmb.df.at [index,'fraction_main_site']        = real_fraction
        pmb.df.at [index,'number_main_patches']       = sites_distribution['main']['number_of_patches']
        pmb.df.at [index,'number_main_sites']         = number_main_sites
        pmb.df.at [index,'secondary_site']            = sites_distribution['secondary']['particle_name']
        pmb.df.at [index,'number_secondary_sites']    = number_secondary_sites
        pmb.df.fillna(pd.NA, inplace=True)
        return


nanoparticle_name = "nanoparticle"
define_nanoparticle(    name                     = nanoparticle_name,
                        core_particle_name       = core_particle,
			surface_density_of_sites = surface_denstity_of_sites*pmb.units('reduced_length^-2'),
                        sites_distribution       = sites_distribution,
                    )

# Save the pyMBE dataframe in a CSV file

pmb.write_pmb_df (filename='df_before.csv')

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

nanoparticle_index    = np.where(pmb.df['name']==nanoparticle_name)
nanoparticle_volume   = pmb.df.loc[pmb.df.index[nanoparticle_index]].nanoparticle_volume.values[0]*pmb.units('reduced_length**3')
volume                = number_of_nanoparticles * nanoparticle_volume / vol_frac_of_nanoparticles
L                     = volume ** (1./3.) # Side of the simulation box

# Create an instance of an espresso system

espresso_system = espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)

# Distributing points evenly on the surface of a sphere

def uniform_distribution_sites_on_sphere(number_of_edges=2, tolerance=1e-6):
    """
    This algorithm is based on iterative force‑based relaxation for distributing points on a sphere, conceptually similar to the Thomson problem (J.J. Thomson, 1904) for minimizing repulsive potential energy of charges on a sphere. See also related uniform sphere point distribution techniques in computational geometry.
    References:
    – Thomson problem — Wikipedia (overview of the physics problem), Wikipedia.
    – Simple schemes for uniform point distribution. Cheng Guan Koay, J Comput Sci. 2011 Dec;2(4):377–381. doi: 10.1016/j.jocs.2011.06.007.
    
    Args:
        number_of_edges(`int`): Number of total points to distribute on the surface of sphere with radius 1 and origin in [0,0,0]. Defaults = 2.
        tolerance(`float`): Set the tolerance of the numerical method. Defaults = 1e-6.
    Returns:
        edges(`list` of `float`): List with the minimized distribution of points.
    """

    # Generates initial configuration
    
    edges = []
    for i in range(number_of_edges):
        theta = pmb.rng.random() * 2*np.pi        
        phi   = np.arcsin(pmb.rng.random() * 2 - 1)
        edges.append((np.cos(theta)*np.cos(phi), np.sin(theta)*np.cos(phi), np.sin(phi)))

    # Iterates until fullfil the tolerance

    while 1:
        # Determine the total force acting on each point.
        forces = []
        for i in range(len(edges)):
            p = edges[i]
            f = (0,0,0)
            ftotal = 0
            for j in range(len(edges)):
                if j == i: continue
                q = edges[j]

                # Find the distance vector, and its length.
                dv = (p[0]-q[0], p[1]-q[1], p[2]-q[2])
                dl = np.sqrt(dv[0]**2 + dv[1]**2 + dv[2]**2)

                # The force vector is dv divided by dl^3. (We divide by dl once to make dv a unit vector, then by dl^2 to make its length correspond to the force.)
                dl3 = dl ** 3
                fv = (dv[0]/dl3, dv[1]/dl3, dv[2]/dl3)

                # Add to the total force on the point p.
                f = (f[0]+fv[0], f[1]+fv[1], f[2]+fv[2])

            # Add in the forces array.
            forces.append(f)

            # Add to the running sum of the total forces/distances.
            ftotal = ftotal + np.sqrt(f[0]**2 + f[1]**2 + f[2]**2)

        # Scale the forces to ensure the points do not move too far in one go. Otherwise there will be chaotic jumping around and never any convergence.
        
        if ftotal > 0.25:
            fscale = 0.25 / ftotal
        else:
            fscale = 1

        # Move each point, and normalise. While we do this, also track the distance each point ends up moving.
        
        dist = 0
        for i in range(len(edges)):
            p = edges[i]
            f = forces[i]
            p2 = (p[0] + f[0]*fscale, p[1] + f[1]*fscale, p[2] + f[2]*fscale)
            pl = np.sqrt(p2[0]**2 + p2[1]**2 + p2[2]**2)
            p2 = (p2[0] / pl, p2[1] / pl, p2[2] / pl)
            dv = (p[0]-p2[0], p[1]-p2[1], p[2]-p2[2])
            dl = np.sqrt(dv[0]**2 + dv[1]**2 + dv[2]**2)
            dist = dist + dl
            edges[i] = p2
      
        # Check for convergence and finish.
        
        if dist < tolerance:
            break

    return edges

# Auxiliary functions to calculate the patchy distribution 

def calculate_distance_vector_point(A,p):
    C = []
    for a in A:
        C.append(((a[0] - p[0])**2 + (a[1] - p[1])**2 + (a[2] - p[2])**2)**(1/2))
    return C

def calculate_patch(points,central_point,patch_size):
    site_positions            = []
    distance_to_central_point = calculate_distance_vector_point(points,central_point)
    points_index              = sorted(range(len(distance_to_central_point)), key=lambda sub: distance_to_central_point[sub])[:patch_size]
    for index in points_index:
        site_positions.append((points[index][0],points[index][1],points[index][2]))
    return distance_to_central_point, site_positions

# Create main and secondary patches
 
def create_patches(nanoparticle_radius, total_number_of_sites, number_main_patches, number_main_sites,tolerance=1e-6, angle_between_patches = 180):
        """
        Creates a list with the lists of positions for `number_main_patches` of the main sites and the list of positions for the secondary sites, using as origin the coordinates [0,0,0].

        Args:
            nanoparticle_radius(`pint.Quantity`): Radius of the nanoparticle expressed in reduced units.
            total_number_of_sites(`int`): Number of total sites on the nanoparticle surface.
            number_main_sites(`int`): Number of only main sites on the nanoparticle surface.
            number_main_patches(`int`): Number of main sites patches.
            tolerance(`float`): Set the tolerance of the numerical method for distributing `total_number_of_sites` in a sphere with radius `nanoparticle_radius`. Defaults = 1e-6.
            angle_between_patches(`float`): If only 2 main patches are defined, this parameter corresponds to the angle between the vectors formed from the center of the nanoparticle and the center of the two patches. Defaults = 180 (the patches are in the poles).
        Returns:
            sites_positions(`list` of `list`): List with the list of the positions of the main and secondary sites in the form: [[sites_positions_main_patch_1],[sites_positions_main_patch_2],...,[sites_positions_secondary_patch]].
        Note:
        The sites are created 1/2 of the reduced unit inside of the nanoparticle surface to avoid overlapping of charges due to electrostatic attractions in abcense of excluded volume. 
        """

        radius_sites                = (nanoparticle_radius - (1/2)*pmb.units('reduced_length')).magnitude 
        root_edges                  = uniform_distribution_sites_on_sphere(number_of_edges = total_number_of_sites, tolerance=tolerance)
        nanoparticle_edges          = np.multiply(root_edges, radius_sites)
        number_main_sites_per_patch = int(number_main_sites/number_main_patches,)
        
        if number_main_patches <= 2:
            initial_edge                = [root_edges[0]]
            sites_positions             = []
            distances_to_center_site_patch_1, sites_positions_main_patch_1 = calculate_patch(points = root_edges, central_point = initial_edge[0], patch_size    = number_main_sites_per_patch)
            sites_positions.append(sites_positions_main_patch_1)
            
            if number_main_patches == 2:
                distance_omega            = (2 * radius_sites**2 - 2 * radius_sites**2 *np.cos(np.radians(angle_between_patches)))**(1/2)
                distance_omega_to_patch_1 = list(np.abs(distances_to_center_site_patch_1 - distance_omega))
                initial_edge.append(root_edges[distance_omega_to_patch_1.index(min(distance_omega_to_patch_1))])
                distances_to_center_site_patch_2, sites_positions_main_patch_2 = calculate_patch(points = root_edges, central_point = initial_edge[1], patch_size    = number_main_sites_per_patch)
                
                # Checking that there is no overlapping
                print(sites_positions_main_patch_1)
                print(sites_positions_main_patch_2)
                overlapped_sites = set(sites_positions_main_patch_1) & set(sites_positions_main_patch_2)
                print(overlapped_sites)
                if len(overlapped_sites) != 0:
                    raise ValueError("The patchies are overlapping in {} sites. Please adjust the angle between them.\n".format(len(overlapped_sites))); 
                else:
                    sites_positions.append(sites_positions_main_patch_2)

        if number_main_patches > 2:
            sites_positions            = []
            initial_patch_edges        = uniform_distribution_sites_on_sphere(number_of_edges = number_main_patches, tolerance=tolerance)
            initial_patch_scaled_edges = np.multiply(initial_patch_edges, radius_sites)
            initial_edge = []
            for scaled_edge in initial_patch_scaled_edges:
                comparison_patch_to_nanoparticle = calculate_distance_vector_point(root_edges,scaled_edge)
                initial_edge.append(root_edges[comparison_patch_to_nanoparticle.index(min(comparison_patch_to_nanoparticle))])
            for main_patch in range(number_main_patches):
                distances_to_center_site_patch, sites_positions_main_patch = calculate_patch(points = root_edges, central_point = initial_edge[main_patch], patch_size    = number_main_sites_per_patch)
                sites_positions.append(sites_positions_main_patch)

            # Checking that there is no overlapping
            overlapped_sites = []
            for i in range(number_main_patches-1):  
                for j in range(i + 1, number_main_patches): 
                    overlapped_sites.append(set(sites_positions[i]) & set(sites_positions[j]))
            for overlap in overlapped_sites: 
                if len(overlap) != 0:
                    raise ValueError("The patchies are overlapping in {} sites. Please adjust the angle between them.\n".format(len(overlap)));
                else:
                    0
        else:
            sites_positions = []

        return sites_positions

'''

edges_base = edges_seed
for patchy in edges_acid:
    edges_base = set(map(tuple, edges_base)).difference(set(patchy))
edges_base = list(edges_base)

total = 0
for i, edge_acid in enumerate(edges_acid):
    print('Sites in patch ',i+1,' : ', len(edge_acid))
    total += len(edge_acid)
print('Total acid sites : ', total)
print('Total basic sites: ', len(edges_base))
print('Angle separation: ', omega)

########## Calculating the distances between charges ########

max_dis = (A_np/N_s)**(1/2)
dis = []
for i in range(len(edges_seed)):
    if i==0:
        dis_ind = ((edges_seed[0][0]-edges_seed[-1][0])**2+(edges_seed[0][1]-edges_seed[-1][1])**2+(edges_seed[0][2]-edges_seed[-1][2])**2)**(1/2)
    if i>=1:
        dis_ind = ((edges_seed[i-1][0]-edges_seed[i][0])**2+(edges_seed[i-1][1]-edges_seed[i][1])**2+(edges_seed[i-1][2]-edges_seed[i][2])**2)**(1/2)
    if dis_ind > max_dis.magnitude:
        dis_ind= dis_ind/2
    dis.append(dis_ind)
avg_dis = np.mean(dis)
dev_dis = np.std(dis)
err_dis = dev_dis / (len(dis))**(1/2)

######## Calculating the dipole and quadrupole moments ######

positions_map = edges_base
for i in range (len(edges_acid)):
    positions_map = np.concatenate((positions_map,edges_acid[i]))

charges_map = np.concatenate((np.ones(N_base),-np.ones(N_acid)))

dp_mnt, dp_mag         = aux.calculate_dipole_moment(charges_map, positions_map)

qq_mnt, qd_mag, qd_eig = aux.calculate_quadrupole_moment(charges_map, positions_map)

print('dipole moment magnitude: (sim)'   , dp_mag*ureg('e*sigma'))
print('dipole moment magnitude: '   , dp_mag*ureg('e*sigma').to('D'))

print('quadrupole moment magnitude: (sim)', qd_mag*ureg('e*sigma**2'))
print('quadrupole moment magnitude: ', qd_mag*ureg('e*sigma**2').to('D*angstrom'))

'''

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
            created_pid_list(`list` of `float`): List with the ids of the particles created into `espresso_system`.
        """ 

        if number_of_nanoparticles <=0:
            return []
        if not _DFm._check_if_name_is_defined_in_df(name=name, df=pmb.df):
            logging.warning(f"Nanoparticle with name '{name}' is not defined in the pyMBE DataFrame, no nanoparticle will be created.")
            return []
        pmb._check_if_name_has_right_type(name=name,
                                           expected_pmb_type="nanoparticle")
        
        # Get information from the nanoparticle type `name` from the df
        
        radius                    = pmb.get_radius_map()
        types                     = pmb.get_type_map()
        core_particle_name        = pmb.df.loc[pmb.df['name'] == name].core_particle.values[0]
        nanoparticle_radius       = radius[types[core_particle_name]]*pmb.units('reduced_length')
        main_site                 = pmb.df.loc[pmb.df['name'] == name].main_site.values[0]
        secondary_site            = pmb.df.loc[pmb.df['name'] == name].secondary_site.values[0]
        total_number_of_sites     = int(pmb.df.loc[pmb.df['name'] == name].total_number_of_sites.values[0])
        number_main_patches       = int(pmb.df.loc[pmb.df['name'] == name].number_main_patches.values[0])
        number_main_sites         = int(pmb.df.loc[pmb.df['name'] == name].number_main_sites.values[0])
        number_secondary_sites    = int(pmb.df.loc[pmb.df['name'] == name].number_secondary_sites.values[0])
        nanoparticle_types        = [core_particle_name, main_site, secondary_site]	
        number_particles_per_type = [number_of_nanoparticles, number_main_sites, number_secondary_sites]
        nanoparticle_edges = create_patches(nanoparticle_radius   = nanoparticle_radius, 
                                            total_number_of_sites = total_number_of_sites, 
                                            number_main_patches   = number_main_patches,
                                            number_main_sites     = number_main_sites,
                                            angle_between_patches = 80)
	
        print(nanoparticle_edges)
        exit()
        # Copy the data of the nanoparticle `number_of_nanoparticles` times in the `df`

        pmb.df = _DFm._copy_df_entry(df                = pmb.df,
                                      name             = name,
                                      column_name      = 'molecule_id',
                                      number_of_copies = number_of_nanoparticles)

        # Get a list of the index in `df` corresponding to the new nanoparticles to be created
        
        nanoparticles_info      = {}
        nanoparticle_index      = np.where(pmb.df['name'] == name)
        nanoparticle_index_list = list(nanoparticle_index[0])[-number_of_nanoparticles:]
        
        for core_particle_position_index, nanoparticle_index in enumerate(nanoparticle_index_list):     
            nanoparticle_id      = _DFm._assign_molecule_id(df             = pmb.df,   
                                                            molecule_index = nanoparticle_index)
            nanoparticles_info[nanoparticle_id] = {}
            
            if list_core_particle_positions is None:
                core_particle_position = None
            else:
                for item in list_core_particle_positions:
                    core_particle_position = [np.array(list_core_particle_positions[core_particle_position_index])]
            '''
            for index_type, nanoparticle_type in enumerate(nanoparticle_types):
                particles_info = pmb.create_particle(name                = nanoparticle_type,
                                                     espresso_system     = espresso_system,
                                                     number_of_particles = 1,
                                                     position            = core_particle_position,
                                                     )
            main_sites_particles_info = pmb.create_particle(name                = core_particle_name,
                                                                 espresso_system     = espresso_system,
                                                                 number_of_particles = 1,
                                                                 position            = core_particle_position,
                                                                 )
            secondary_sites_particles_info = pmb.create_particle(name                = core_particle_name,
                                                                 espresso_system     = espresso_system,
                                                                 number_of_particles = 1,
                                                                 position            = core_particle_position,
                                                                 )
                    # Add the correct molecule_id to all particles in the residue
                    for index in self.df[self.df['residue_id']==residue_id].index:
                        _DFm._add_value_to_df(df = self.df,
                                              key = ('molecule_id',''),
                                              index = int (index),
                                              new_value = molecule_id,
                                              overwrite = True)
                    central_bead_id = residues_info[residue_id]['central_bead_id']
                    previous_residue = residue
                    residue_position = espresso_system.part.by_id(central_bead_id).pos
                    previous_residue_id = central_bead_id
                    first_residue = False
            '''
 
create_nanoparticle(name=nanoparticle_name, espresso_system=espresso_system, number_of_nanoparticles=number_of_nanoparticles)

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
