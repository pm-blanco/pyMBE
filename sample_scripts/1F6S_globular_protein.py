import os
import sys
import inspect
from tqdm import tqdm
import espressomd
import argparse

from espressomd import interactions
from espressomd.io.writer import vtf
from espressomd import electrostatics 

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library()

from handy_scripts.handy_functions import calculate_net_charge_in_molecule

# Here you can adjust the width of the panda columns displayed when running the code 
pmb.pd.options.display.max_colwidth = 10

#This line allows you to see the complete amount of rows in the dataframe
pmb.pd.set_option('display.max_rows', None)

# The trajectories of the simulations will be stored using espresso built-up functions in separed files in the folder 'frames'
if not os.path.exists('./frames'):
    os.makedirs('./frames')

parser = argparse.ArgumentParser(description='Script to run globular protein simulation in espressomd')
parser.add_argument('-pH', type=float, required= True,  help='Expected pH value')
args = parser.parse_args ()

#System Parameters 
pH_value = args.pH 
print (pH_value)

#agregar check para cuando no se da un valor pH 

# print (f"Running simulation for: {pH_value}")
pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)

c_salt    =  0.01  * pmb.units.mol / pmb.units.L  
c_protein =  2e-4 * pmb.units.mol / pmb.units.L 
Box_V =  1. / (pmb.N_A*c_protein)
Box_L = Box_V**(1./3.) 

SEED = 77 
solvent_permitivity = 78.3

#Simulation Parameters
t_max = 1e3 #  in LJ units of time
stride_obs = 10 #  in LJ units of time
stride_traj = 100 # in LJ units of time
dt = 0.01
N_samples = int (t_max / stride_obs)
integ_steps = int (stride_obs/dt)

probability_reaction = 0.5 

bead_size = 0.355*pmb.units.nm
epsilon = 1*pmb.units('reduced_energy')

WCA = True
Electrostatics = True

espresso_system = espressomd.System(box_l=[Box_L.to('reduced_length').magnitude] * 3)
espresso_system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()

#Directory of the protein model 
protein_name = '1f6s'
protein_filename = os.path.join (parentdir,'sample_scripts/coarse_grain_model_of_1f6s.vtf' )


#Reads the VTF file of the protein model
protein_positions = pmb.load_protein_vtf_in_df (name=protein_name,filename=protein_filename)

#Create dictionary with the value of epsilon for each residue
protein_sequence = pmb.df.loc[pmb.df['name']== protein_name].sequence.values[0]
clean_sequence = pmb.protein_sequence_parser(sequence=protein_sequence)

epsilon_dict = {}
for residue in clean_sequence:
    if residue not in epsilon_dict.keys():
        epsilon_dict [residue] = epsilon

    epsilon_dict ['CA'] = epsilon

#Define epsilon for each particle
pmb.define_epsilon_value_of_particles (eps_dict = epsilon_dict)

# print (pmb.filter_df(pmb_type='particle'))

#Metal Ion calcium definition
pmb.define_particle(name='Ca',q=+2,epsilon=epsilon)

# Solution 
cation_name = 'Na'
anion_name = 'Cl'
pmb.define_particle(name=cation_name,  q=1, diameter=0.2*pmb.units.nm, epsilon=epsilon)
pmb.define_particle(name=anion_name,  q=-1, diameter=0.2*pmb.units.nm,  epsilon=epsilon)

# Here we upload the pka set from the reference_parameters folder 

pmb.load_pka_set (filename=os.path.join(parentdir,'reference_parameters/pka_sets/Nozaki1967.txt'))

#We create the protein in espresso 
pmb.create_protein_in_espresso(name=protein_name,
                               number_of_proteins=1,
                               espresso_system=espresso_system,
                               positions=protein_positions)

#Here we activate the motion of the protein
# pmb.activate_motion_of_rigid_object(espresso_system=espresso_system,name=protein_name)

# If we want to center the first protein created
protein_id = pmb.df.loc[pmb.df['name']==protein_name].molecule_id.values[0]
pmb.center_molecule_in_simulation_box (molecule_id=protein_id,espresso_system=espresso_system)

pmb.create_counterions_in_espresso (pmb_object_name='particle',cation_name=cation_name,anion_name=anion_name,espresso_system=espresso_system)

c_salt_calculated = pmb.create_added_salt_in_espresso (espresso_system=espresso_system,cation_name=cation_name,anion_name=anion_name,c_salt=c_salt)

basic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='basic')].name.to_list()
acidic_groups = pmb.df.loc[(~pmb.df['particle_id'].isna()) & (pmb.df['acidity']=='acidic')].name.to_list()
list_ionisible_groups = basic_groups + acidic_groups
total_ionisible_groups = len (list_ionisible_groups)

print('The box length of the system is', Box_L.to('reduced_length'), Box_L.to('nm'))
print('The ionisable groups in the protein are ', list_ionisible_groups)
print ('The total amount of ionizable groups are:',total_ionisible_groups)

RE, sucessfull_reactions_labels = pmb.setup_constantpH_reactions_in_espresso (counter_ion=cation_name, constant_pH= pH_value, SEED = SEED )
print('The acid-base reaction has been sucessfully setup for ', sucessfull_reactions_labels)

type_map = pmb.get_type_map()
types = list (type_map.values())
espresso_system.setup_type_map( type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
RE.set_non_interacting_type (type=non_interacting_type)
print('The non interacting type is set to ', non_interacting_type)

# Setup the potential energy

if (WCA):

    pmb.setup_lj_interactions_in_espresso (espresso_system=espresso_system)
    
    print('\nMinimazing system energy\n')
    espresso_system.cell_system.skin = 0.4
    espresso_system.time_step = dt 
    print('steepest descent')
    espresso_system.integrator.set_steepest_descent(f_max=0, gamma=0.1, max_displacement=0.1)
    espresso_system.integrator.run(1000)
    print('velocity verlet')
    espresso_system.integrator.set_vv()  # to switch back to velocity Verlet
    espresso_system.integrator.run(1000)
    espresso_system.thermostat.turn_off()
    print('\nMinimization finished \n')

    if (Electrostatics):
        #Electrostatic energy setup 
        print ('Electrostatics setup')
        bjerrum_length = pmb.e.to('reduced_charge')**2 / (4 *pmb.units.pi*pmb.units.eps0* solvent_permitivity * pmb.kT.to('reduced_energy'))
        coulomb_prefactor = bjerrum_length.to('reduced_length') * pmb.kT.to('reduced_energy')
        coulomb = espressomd.electrostatics.P3M ( prefactor = coulomb_prefactor.magnitude, accuracy=1e-3)

        print('\nBjerrum length ', bjerrum_length.to('nm'), '=', bjerrum_length.to('reduced_length'))

        espresso_system.time_step = dt
        espresso_system.actors.add(coulomb)

        # save the optimal parameters and add them by hand
        p3m_params = coulomb.get_params()
        espresso_system.actors.remove(coulomb)

        coulomb = espressomd.electrostatics.P3M (
                                    prefactor = coulomb_prefactor.magnitude,
                                    accuracy = 1e-3,
                                    mesh = p3m_params['mesh'],
                                    alpha = p3m_params['alpha'] ,
                                    cao = p3m_params['cao'],
                                    r_cut = p3m_params['r_cut'],
                                    tune = False,
                                        )

        espresso_system.actors.add(coulomb)

        print("\nElectrostatics successfully added to the system \n")

#Save the initial state 

n_frame = 0
pmb.write_output_vtf_file(espresso_system=espresso_system,filename=f"frames/trajectory{n_frame}.vtf")

print (f'Optimizing skin\n')
espresso_system.time_step = dt 
espresso_system.integrator.set_vv()
espresso_system.thermostat.set_langevin(kT=pmb.kT.to('reduced_energy').magnitude, gamma=0.1, seed=SEED)

espresso_system.cell_system.tune_skin ( min_skin = 10, 
                                        max_skin = espresso_system.box_l[0]/2., tol=1e-3, 
                                        int_steps=1000, adjust_max_skin=True)

print('Optimized skin value: ', espresso_system.cell_system.skin, '\n')

observables_df = pmb.pd.DataFrame()
time_step = []
net_charge_list = []
net_charge_amino_save = {}

Z_sim=[]
particle_id_list = pmb.df.loc[~pmb.df['molecule_id'].isna()].particle_id.dropna().to_list()

pmb.df.to_csv('df.csv')

for step in tqdm(range(N_samples)):
        
        if pmb.np.random.random() > probability_reaction:
            espresso_system.integrator.run (steps = integ_steps)
        else:
            RE.reaction( reaction_steps = total_ionisible_groups)

        calculated_net_charge = calculate_net_charge_in_molecule (espresso_system=espresso_system,pmb_df = pmb.df, name = protein_name)

        net_charge = calculated_net_charge['net_charge']
        net_charge_residues = calculated_net_charge ['net_charge_residues']

        time_step.append (str(espresso_system.time))
        net_charge_list.append (net_charge)

        if len(net_charge_amino_save.keys()) == 0:
            for amino in net_charge_residues.keys():
                net_charge_amino_save [amino] = []
        for amino in net_charge_residues.keys():            
            net_charge_amino_save [amino].append (net_charge_residues[amino])

        if (step % stride_traj == 0  ):
            n_frame +=1
            pmb.write_output_vtf_file(espresso_system=espresso_system,filename=f"frames/trajectory{n_frame}.vtf")

observables_df['time'] = time_step 
observables_df['Znet'] = net_charge_list

for amino in net_charge_amino_save.keys():
    observables_df[amino] = net_charge_amino_save[amino]

observables_df.to_csv(f'pH-{pH_value}_observables.csv',index=False)

print(pmb.df)
print (observables_df)

# from espressomd import visualization

# espresso_system.time_step=1e-3
# espresso_system.thermostat.set_langevin(kT=1, gamma=1.0, seed=24)
# espresso_system.cell_system.tune_skin(min_skin=0.01, max_skin =4.0, tol=0.1, int_steps=1000)
# visualizer = espressomd.visualization.openGLLive(espresso_system, bond_type_radius=[0.3], background_color=[1, 1, 1],particle_type_colors=[[1.02,0.51,0], # Brown
#                 [1,1,1],  # Grey
#                 [2.55,0,0], # Red
#                 [0,0,2.05],  # Blue
#                 [0,0,2.05],  # Blue
#                 [2.55,0,0], # Red
#                 [2.05,1.02,0]])     
# visualizer.run(1)

# filename = 'protein.png'

# visualizer = visualization.openGLLive(
#     espresso_system, bond_type_radius=[0.3], particle_coloring='type', draw_axis=False, background_color=[1, 1, 1],
# particle_type_colors=[[1.02,0.51,0], # Brown
#                 [1,1,1],  # Grey
#                 [2.55,0,0], # Red
#                 [0,0,2.05],  # Blue
#                 [0,0,2.05],  # Blue
#                 [2.55,0,0], # Red
#                 [2.05,1.02,0]]) # Orange
# visualizer.screenshot(filename)

# from PIL import Image
# img = Image.open(filename)
# img.show()
