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

#######################################################
# Loading modules 
#######################################################

# Load python modules
import espressomd
from pathlib import Path
import numpy as np
import pandas as pd
import tqdm
from scipy import interpolate
import argparse

# Import pyMBE
import pyMBE
from pyMBE.lib import analysis

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

# Load some functions from the handy_scripts library for convenience
from pyMBE.lib.handy_functions import setup_electrostatic_interactions
from pyMBE.lib.handy_functions import relax_espresso_system
from pyMBE.lib.handy_functions import setup_langevin_dynamics
from pyMBE.lib.handy_functions import do_reaction


#######################################################
# Setting parameters for the simulation
#######################################################

##### Read command-line arguments using argparse
parser = argparse.ArgumentParser()
parser.add_argument("--c_salt_res", 
                    type=float, 
                    help="Concentration of NaCl in the reservoir in mol/l.",
                    default=0.01)
parser.add_argument("--c_mon_sys", 
                    type=float, 
                    help="Concentration of acid monomers in the system in mol/l.",
                    default=0.01)
parser.add_argument("--pH_res", 
                    type=float, 
                    help="pH-value in the reservoir.",
                    default=7)
parser.add_argument("--pKa_value", 
                    type=float, 
                    help="pKa-value of the polyacid monomers.",
                    default=4)
parser.add_argument('--mode',
                    type=str,
                    default= "short-run",
                    choices=["short-run","long-run", "test"],
                    help='sets for how long the simulation runs')
parser.add_argument('--output',
                    type=Path,
                    required= False,
                    default=Path(__file__).parent / "time_series" / "grxmc",
                    help='output directory')
parser.add_argument('--no_verbose', 
                    action='store_false', 
                    help="Switch to deactivate verbose",
                    default=True)
args = parser.parse_args()

# Inputs
inputs={"csalt": args.c_salt_res,
        "cmon": args.c_mon_sys,
        "pH": args.pH_res,
        "pKa": args.pKa_value}
mode=args.mode
verbose=args.no_verbose

# Units and general parameters
pmb.set_reduced_units(unit_length=0.355*pmb.units.nm)
solvent_permittivity = 78.9

# Integration parameters
dt = 0.01
langevin_seed = 42

# Parameters of the polyacid model
Chain_length = 50
if mode == "test":
    N_chains = 2
else:
    N_chains = 16
polyacid_name = 'polyacid'

# Add the polyacid to the pmb instance
pmb.define_particle(name='A', 
                    acidity='acidic', 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'), 
                    pka=args.pKa_value)
pmb.define_residue(name='rA', 
                   central_bead="A", 
                   side_chains=[])
pmb.define_molecule(name=polyacid_name, 
                    residue_list=['rA']*Chain_length)


bond_type = 'FENE'
fene_spring_constant = 30 * pmb.units('reduced_energy / reduced_length**2')
fene_r_max = 1.5 * pmb.units('reduced_length')

fene_bond = {'k'      : fene_spring_constant,
             'd_r_max': fene_r_max, 
            }

pmb.define_bond(bond_type = bond_type, 
                bond_parameters = fene_bond, 
                particle_pairs = [['A','A']])

# Parameters of the small ions
proton_name = 'Hplus'
hydroxide_name = 'OHminus'
sodium_name = 'Na'
chloride_name = 'Cl'

pmb.define_particle(name=proton_name, 
                    z=1, 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=hydroxide_name,  
                    z=-1, 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=sodium_name, 
                    z=1, 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'))
pmb.define_particle(name=chloride_name, 
                    z=-1, 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'))

# System parameters (some are read from the command line)
c_mon_sys = args.c_mon_sys * pmb.units.mol/ pmb.units.L
c_salt_res = args.c_salt_res * pmb.units.mol/ pmb.units.L
pH_res = args.pH_res 
pka_set = {'A': {"pka_value": args.pKa_value, "acidity": "acidic"}}
volume = N_chains * Chain_length/(pmb.N_A*c_mon_sys)
L = volume ** (1./3.)
if verbose:
    print("Box length:", L.to('reduced_length').magnitude)


#######################################################
# Setting up the espresso system
#######################################################

# Create an instance of an espresso system
espresso_system = espressomd.System(box_l = [L.to('reduced_length').magnitude]*3)
espresso_system.time_step=dt
espresso_system.cell_system.skin=0.4
if verbose:
    print("Created espresso object")

# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)
if verbose:
    print("Added bonds")

# Create molecules and ions in the espresso system
pmb.create_molecule(name=polyacid_name, 
                    number_of_molecules=N_chains, 
                    espresso_system=espresso_system)
pmb.create_counterions(object_name=polyacid_name, 
                       cation_name=proton_name, 
                       anion_name=hydroxide_name, 
                       espresso_system=espresso_system)
c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system, 
                                          cation_name=sodium_name, 
                                          anion_name=chloride_name, 
                                          c_salt=c_salt_res)
if verbose:
    print("Created molecules")

# Set up the reactions
path_to_ex_pot=pmb.root / "parameters" / "salt"
monovalent_salt_ref_data=pd.read_csv(f"{path_to_ex_pot}/excess_chemical_potential_excess_pressure.csv")
ionic_strength = pmb.units.Quantity(monovalent_salt_ref_data["cs_bulk_[1/sigma^3]"].values, "1/reduced_length**3")
excess_chemical_potential = pmb.units.Quantity(monovalent_salt_ref_data["excess_chemical_potential_[kbT]"].values, "reduced_energy")
excess_chemical_potential_interpolated = interpolate.interp1d(ionic_strength.m_as("1/reduced_length**3"), 
                                                                                excess_chemical_potential.m_as("reduced_energy"))
activity_coefficient_monovalent_pair = lambda x: np.exp(excess_chemical_potential_interpolated(x.to('1/(reduced_length**3 * N_A)').magnitude))
if verbose:
    print("Setting up reactions...")
grxmc, labels, ionic_strength_res = pmb.setup_grxmc_reactions(pH_res=pH_res, 
                                                              c_salt_res=c_salt_res, 
                                                              proton_name=proton_name, 
                                                              hydroxide_name=hydroxide_name, 
                                                              salt_cation_name=sodium_name, 
                                                              salt_anion_name=chloride_name, 
                                                              activity_coefficient=activity_coefficient_monovalent_pair, 
                                                              pka_set=pka_set)
if verbose:
    print('The acid-base reaction has been sucessfully set up for ', labels)

# Setup espresso to track the ionization of the acid groups
type_map = pmb.get_type_map()
types = list(type_map.values())
espresso_system.setup_type_map(type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
grxmc.set_non_interacting_type (type=non_interacting_type)

#Set up the interactions
pmb.setup_lj_interactions(espresso_system=espresso_system)

# Relax the system
relax_espresso_system(espresso_system=espresso_system, 
                      seed=langevin_seed)
setup_langevin_dynamics(espresso_system=espresso_system, 
                        kT = pmb.kT, 
                        seed = langevin_seed,
                        time_step=dt,
                        tune_skin=False)
if verbose:
    print("Running warmup without electrostatics")
for i in tqdm.trange(100, disable=not verbose):
    espresso_system.integrator.run(steps=1000)
    do_reaction(grxmc, steps=1000)

setup_electrostatic_interactions(units=pmb.units,
                                espresso_system=espresso_system,
                                kT=pmb.kT,
                                solvent_permittivity=solvent_permittivity,
                                verbose=verbose)
espresso_system.thermostat.turn_off()
relax_espresso_system(espresso_system=espresso_system,
                      seed=langevin_seed,
                      max_displacement=0.01)
setup_langevin_dynamics(espresso_system=espresso_system, 
                        kT = pmb.kT, 
                        seed = langevin_seed,
                        time_step=dt,
                        tune_skin=False)

if verbose:
    print("Running warmup with electrostatics")
if mode == "long-run":
    N_warmup_loops = 1000
else:
    N_warmup_loops = 100
for i in tqdm.trange(N_warmup_loops, disable=not verbose):
    espresso_system.integrator.run(steps=1000)
    do_reaction(grxmc, steps=100)


# Main loop
print("Started production run.")

labels_obs=["time", "alpha"]
time_series={}

for label in labels_obs:
    time_series[label]=[]

if mode == "long-run":
    N_production_loops = 5000
else:
    N_production_loops = 100
for i in tqdm.trange(N_production_loops, disable=not verbose):
    espresso_system.integrator.run(steps=1000)
    do_reaction(grxmc, steps=100)

    # Measure time
    time_series["time"].append(espresso_system.time)

    # Measure degree of ionization
    charge_dict=pmb.calculate_net_charge(espresso_system=espresso_system, molecule_name=polyacid_name, dimensionless=True)
    time_series["alpha"].append(np.abs(charge_dict["mean"])/Chain_length)

data_path = args.output
data_path.mkdir(parents=True, exist_ok=True)

time_series=pd.DataFrame(time_series)
filename=analysis.built_output_name(input_dict=inputs)

time_series.to_csv(data_path / f"{filename}_time_series.csv", index=False)
