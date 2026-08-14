#
# Copyright (C) 2024-2026 pyMBE-dev team
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

from pathlib import Path
import tqdm
import espressomd
import argparse
import numpy as np
import pandas as pd 
from espressomd.io.writer import vtf
# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(seed=42)

#Import functions from handy_functions script 
from pyMBE.lib.handy_functions import define_protein_AA_particles, define_protein_AA_residues
from pyMBE.lib import analysis
# Here you can adjust the width of the panda columns displayed when running the code 
pd.options.display.max_colwidth = 10

#This line allows you to see the complete amount of rows in the dataframe
pd.set_option('display.max_rows', None)

parser = argparse.ArgumentParser(description='Script to run globular protein simulation in espressomd')

parser.add_argument('--pdb',
                    type=str, 
                    required= True,  
                    help='PDB code of the protein')
parser.add_argument('--pH', 
                    type=float, 
                    required= True,  
                    help='pH value')
parser.add_argument('--path_to_cg', 
                    type=Path,
                    required= True,  
                    help='Path to the CG structure of the protein')
parser.add_argument('--move_protein', 
                    action="store_true",
                    default=False,  
                    help='Activates the motion of the protein')
parser.add_argument('--ideal', 
                    action="store_true",
                    default=False,  
                    help='Sets up an ideal system without steric and electrostatic interactions ')
parser.add_argument('--mode',
                    type=str,
                    default= "short-run",
                    choices=["short-run","long-run", "test"],
                    help='sets for how long the simulation runs')
parser.add_argument('--output',
                    type=Path,
                    required= False,
                    default=Path(__file__).parent / "time_series" / "globular_protein",
                    help='output directory')
parser.add_argument('--no_verbose', 
                    action='store_false', 
                    help="Switch to deactivate verbose",default=True)
args = parser.parse_args ()
mode=args.mode
verbose=args.no_verbose
protein_name = args.pdb
pH_value = args.pH 
model = '2beadAA'
inputs={"pH": args.pH,
        "pdb": args.pdb}

#System Parameters 
langevin_seed = 77 

c_salt    =  0.01  * pmb.units.mol / pmb.units.L  
c_protein =  2e-4 * pmb.units.mol / pmb.units.L 
Box_V =  1. / (pmb.N_A*c_protein)
Box_L = Box_V**(1./3.) 
solvent_permitivity = 78.3
epsilon = 1*pmb.units('reduced_energy')
sigma = 1*pmb.units("reduced_length")
ion_size = 0.4*pmb.units.nm

#Simulation Parameters
 #  in LJ units of time
dt = 0.01
stride_traj = 100 # in LJ units of time

if mode == 'short-run':
    stride_obs = 10 #  in LJ units of time
    integ_steps = int (stride_obs/dt)
    t_max = 5e3
    N_samples = int (t_max / stride_obs)
    
elif mode == 'long-run':
    stride_obs = 10 #  in LJ units of time
    integ_steps = int (stride_obs/dt)
    t_max = 1e5
    N_samples = int (t_max / stride_obs)
    
elif mode == 'test':
    t_max = 1e3
    stride_obs = 10 #  in LJ units of time
    integ_steps = int (stride_obs/dt)
    N_samples = int (t_max / stride_obs)
else: 
    raise RuntimeError()    

#Switch for Electrostatics and WCA interactions 
WCA = True
Electrostatics = True

if args.ideal:
    WCA=False
    Electrostatics=False

data_path = args.output

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
frames_path = args.output / "frames"
frames_path.mkdir(parents=True, exist_ok=True)
box_l=[Box_L.to('reduced_length').magnitude] * 3
espresso_system = espressomd.System(box_l=box_l)
espresso_system.time_step=dt
espresso_system.cell_system.skin=0.4
#Reads the VTF file of the protein model
topology_dict, sequence = pmb.read_protein_vtf (filename=args.path_to_cg)
# Here we upload the pka set from the reference_parameters folder
path_to_pka=pmb.root / "parameters" / "pka_sets" / "Nozaki1967.json"
pmb.load_pka_set(filename=path_to_pka)
pka_set = pmb.get_pka_set()

#Defines the protein in the pyMBE database
define_protein_AA_particles(topology_dict=topology_dict,
                            pmb=pmb,
                            pka_set=pka_set)
residue_list = define_protein_AA_residues(sequence=sequence,
                                        model=model,
                                        pmb=pmb)

# Define a residue for the metal ion
if args.pdb == "1f6s":
    pmb.define_residue(name="AA-Ca",
                        central_bead="Ca",
                        side_chains=[])
        
pmb.define_protein(name=protein_name, 
                   sequence=sequence, 
                   model = model)

# Here we define the solution particles in the pyMBE database
cation_name = 'Na'
anion_name = 'Cl'

pmb.define_particle(name = cation_name, 
                    z = 1, 
                    sigma=sigma, 
                    epsilon=epsilon,
                    offset=ion_size-sigma)

pmb.define_particle(name = anion_name,  
                    z =-1, 
                    sigma=sigma, 
                    epsilon=epsilon,
                    offset=ion_size-sigma)

#We create the protein in espresso 
protein_id = pmb.create_protein(name=protein_name,
                                number_of_proteins=1,
                                box_l=box_l,
                                topology_dict=topology_dict)[0]

pmb.set_simulation_engine(espresso_system)
pmb.add_instances_to_engine()
#Here we activate the motion of the protein 
if args.move_protein:
    pmb.enable_motion_of_rigid_object(instance_id=protein_id,
                                      pmb_type="protein")

# Here we put the protein on the center of the simulation box
pmb.center_object_in_simulation_box(instance_id=protein_id,
                                    box_l=box_l,
                                    pmb_type="protein")

if not args.ideal:
    # Estimate the radius of the protein
    protein_ids = pmb.get_particle_id_map(object_name=protein_name)["all"]
    protein_radius=0
    protein_center=espresso_system.box_l/2
    for pid in protein_ids:
        protein_part = espresso_system.part.by_id(pid)
        dist = protein_part.pos - protein_center
        dist = np.linalg.norm(dist)
        if dist > protein_radius:
            protein_radius = dist
    # Create counter-ions 
    protein_net_charge = pmb.calculate_net_charge(
                                                object_name=protein_name,
                                                pmb_type="protein",
                                                dimensionless=True)["mean"]
    ## Get coordinates outside the volume occupied by the protein
    counter_ion_coords=pmb.generate_coordinates_outside_sphere(center=protein_center,
                                                                radius=protein_radius,
                                                                max_dist=Box_L.m_as("reduced_length"),
                                                                n_samples=abs(protein_net_charge))
    if protein_net_charge > 0:
        counter_ion_name=anion_name
    elif protein_net_charge < 0:
        counter_ion_name=cation_name

    pmb.create_particle(name=counter_ion_name,
                        box_l=box_l,
                        number_of_particles=int(abs(protein_net_charge)),
                        position=counter_ion_coords)

    # Create added salt ions
    N_ions= int((Box_V.to("reduced_length**3")*c_salt.to('mol/reduced_length**3')*pmb.N_A).magnitude)

    ## Get coordinates outside the volume occupied by the protein
    added_salt_ions_coords=pmb.generate_coordinates_outside_sphere(center=protein_center,
                                                                radius=protein_radius,
                                                                max_dist=Box_L.m_as("reduced_length"),
                                                                n_samples=N_ions*2)
    ## Create cations
    pmb.create_particle(name=cation_name,
                        box_l=box_l,
                        number_of_particles=N_ions,
                        position=added_salt_ions_coords[:N_ions])
    ## Create anions
    pmb.create_particle(name=anion_name,
                        box_l=box_l,
                        number_of_particles=N_ions,
                        position=added_salt_ions_coords[N_ions:])
    pmb.add_instances_to_engine()

#Here we calculated the ionisable groups
acid_base_ids = []
list_ionisable_groups = []
for name in pka_set.keys():
    part_ids = pmb.db.find_instance_ids_by_name(pmb_type="particle",
                                                name=name)
    if part_ids:
        acid_base_ids+=part_ids
        list_ionisable_groups+=[name]  
total_ionisable_groups = len(acid_base_ids)

if verbose:
    print(f"The box length of the system is {Box_L.to('reduced_length')} {Box_L.to('nm')}")
    print(f"The ionisable groups in the protein are {list_ionisable_groups}")
    print(f"The total amount of ionisable groups is {total_ionisable_groups}")

#Setup of the reactions in espresso 
cpH = pmb.setup_cpH(counter_ion=cation_name, 
                    constant_pH= pH_value)
if verbose:
    print("The acid-base reaction has been successfully set up for:")
    print(pmb.get_reactions_df())

type_map = pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
cpH.set_non_interacting_type (type=non_interacting_type)
if verbose:
    print(f"The non interacting type is set to {non_interacting_type}")

#Save the initial state 
n_frame = 0
with open(frames_path / f"trajectory{n_frame}.vtf", mode='w+t') as coordinates:
    vtf.writevsf(espresso_system, coordinates)
    vtf.writevcf(espresso_system, coordinates)
# Setup the potential energy

if WCA:
    pmb.setup_lj_interactions()
    pmb.simulation_engine.relax_espresso_system(
                          seed=langevin_seed)
    if Electrostatics:
        pmb.simulation_engine.setup_electrostatic_interactions(units=pmb.units,
                                        kT=pmb.kT)
        
pmb.simulation_engine.setup_langevin_dynamics(
                        kT = pmb.kT, 
                        seed = langevin_seed)

observables_df = pd.DataFrame()
time_step = []
net_charge_list = []

Z_sim=[]

pmb.save_database (folder=data_path/"database")

#Here we start the main loop over the Nsamples 

labels_obs=["time","charge"]
time_series={}

for label in labels_obs:
    time_series[label]=[]

AA_label_list = []
for amino in list_ionisable_groups:
    label = f'AA-{amino}'
    time_series[f"charge_{label}"] = []
    AA_label_list.append(label)

for step in tqdm.trange(N_samples, disable=not verbose):
    espresso_system.integrator.run (steps = integ_steps)
    pmb.simulation_engine.do_reaction(cpH, steps=total_ionisable_groups)
    protein_net_charge = pmb.calculate_net_charge(
                                                object_name=protein_name,
                                                pmb_type="protein",
                                                dimensionless=True)["mean"]
    # Store observables
    time_series["time"].append(espresso_system.time)
    time_series["charge"].append(protein_net_charge)
    charge_residues_per_type = {}
    for label in AA_label_list:
        charge_residues_per_type[label]=[]
        charge_res=pmb.calculate_net_charge (
                                            object_name=label,
                                            pmb_type="residue",
                                            dimensionless=True)["mean"]
        time_series[f"charge_{label}"].append(charge_res)
    if step % stride_traj == 0 :
        n_frame +=1
        with open(frames_path / f"trajectory{n_frame}.vtf", mode='w+t') as coordinates:
            vtf.writevsf(espresso_system, coordinates)
            vtf.writevcf(espresso_system, coordinates)

data_path.mkdir(parents=True, exist_ok=True)
    
time_series=pd.DataFrame(time_series)
filename=analysis.built_output_name(input_dict=inputs)

time_series.to_csv(data_path / f"{filename}_time_series.csv", index=False)

if verbose:
    print("*** DONE ***")
