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
#

import espressomd
from pathlib import Path
import numpy as np
import pandas as pd
import tqdm
from scipy import interpolate
import argparse
from pyMBE.lib.lattice import DiamondLattice
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
parser.add_argument("--mpc",
                    type=int,
                    help="Number of monomers per chain ",
                    default=40)
parser.add_argument("--L_fraction",
                    type=float,
                    help="Target fraction of the maximum box length to achieve during the volume compression of the gel",
                    default=1)
parser.add_argument("--csalt_res",
                    type=float,
                    help="Concentration of NaCl in the reservoir in mol/l.",
                    default=0.01)
parser.add_argument("--pH_res",
                    type=float,
                    help="pH-value in the reservoir.",
                    default=7)
parser.add_argument("--pKa",
                    type=float,
                    help="pKa-value of the hydrogel monomers",
                    default=4)
parser.add_argument('--mode',
                    type=str,
                    default= "short-run",
                    choices=["long-run", "short-run", "test"],
                    help='sets for how long the simulation runs')
parser.add_argument('--output',
                    type=Path,
                    required= False,
                    default=Path(__file__).parent / "time_series" / "weak_polyacid_hydrogel_grxmc",
                    help='output directory')


args = parser.parse_args()
mode=args.mode
c_salt_res = args.csalt_res * pmb.units.mol/ pmb.units.L
solvent_permittivity = 78.9
seed=42
NodeType = "N"
pmb.define_particle(name=NodeType, 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'))

BeadType = "C"
pmb.define_particle(name=BeadType, 
                    acidity="acidic", 
                    sigma=1*pmb.units('reduced_length'), 
                    epsilon=1*pmb.units('reduced_energy'), 
                    pka=args.pKa)

pmb.define_residue(name='Res', 
                   central_bead=BeadType, 
                   side_chains=[])

bond_type = 'FENE'
bond_length = 0.966 * pmb.units("reduced_length")
fene_spring_constant = 30 * pmb.units('reduced_energy / reduced_length**2')
fene_r_max = 1.5 * pmb.units('reduced_length')
fene_r0 = 0 * pmb.units('reduced_length')

fene_bond = {'k'      : fene_spring_constant,
             'd_r_max': fene_r_max,
             "r_0": fene_r0}

pmb.define_bond(bond_type = bond_type, 
                bond_parameters = fene_bond, 
                particle_pairs = [[BeadType,NodeType],[BeadType, BeadType]])

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

diamond_lattice = DiamondLattice(args.mpc, bond_length)
espresso_system = espressomd.System(box_l = [diamond_lattice.box_l]*3)
pmb.add_bonds_to_espresso(espresso_system = espresso_system)
lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)

# Setting up node topology
indices = diamond_lattice.indices
node_topology = []

for index in range(len(indices)):
    node_topology.append({"particle_name": NodeType,
                          "lattice_index": indices[index]})

# Setting up chain topology
node_labels = lattice_builder.node_labels
chain_labels = lattice_builder.chain_labels
reverse_node_labels = {v: k for k, v in node_labels.items()}
chain_topology = []
residue_list = ["Res"]*args.mpc

for chain_data in chain_labels.items():
    chain_label = chain_data[1]
    node_label_pair = chain_data[0]
    node_label_s, node_label_e = [int(x) for x in node_label_pair.strip("()").split(",")]
    chain_topology.append({'node_start':reverse_node_labels[node_label_s],
                              'node_end': reverse_node_labels[node_label_e],
                              'residue_list':residue_list})

pmb.define_hydrogel("my_hydrogel", node_topology, chain_topology)
hydrogel_info = pmb.create_hydrogel("my_hydrogel", espresso_system)

c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system, 
                                          cation_name=sodium_name, 
                                          anion_name=chloride_name, 
                                          c_salt=c_salt_res)

print(f"Salt concentration {c_salt_calculated.to('mol/L')} added")

print("*** Setting LJ interactions ***")

dt = 0.01  # Timestep 
espresso_system.time_step = dt
pmb.setup_lj_interactions(espresso_system=espresso_system)

print("*** Relaxing the system... ***")
relax_espresso_system(espresso_system=espresso_system,
                      seed=seed,
                      Nsteps_iter_relax=100)

setup_langevin_dynamics(espresso_system=espresso_system,
                        kT = pmb.kT,
                        seed = seed,
                        time_step=dt,
                        tune_skin=False)

L_max = diamond_lattice.box_l
L_target = args.L_fraction * L_max
steps_size = 0.1
steps_needed = int(abs((L_max - L_target)/steps_size))

print("*** Box progressive rescaling... ***")

for j in tqdm.trange(steps_needed):
    espresso_system.change_volume_and_rescale_particles(
        d_new=espresso_system.box_l[0] - steps_size, dir="xyz")
    espresso_system.integrator.run(500)
# Just to make sure that the system has the target size
espresso_system.change_volume_and_rescale_particles(
    d_new=L_target, dir="xyz")
relax_espresso_system(espresso_system=espresso_system,
                      seed=seed,
                      Nsteps_iter_relax=100,
                      max_displacement=0.01)

print("*** Setting up the GRxMC method and electrostatics... ***")

# Set up the reactions
path_to_ex_pot=pmb.root / "parameters" / "salt"
monovalent_salt_ref_data=pd.read_csv(f"{path_to_ex_pot}/excess_chemical_potential_excess_pressure.csv")
ionic_strength = pmb.units.Quantity(monovalent_salt_ref_data["cs_bulk_[1/sigma^3]"].values, "1/reduced_length**3")
excess_chemical_potential = pmb.units.Quantity(monovalent_salt_ref_data["excess_chemical_potential_[kbT]"].values, "reduced_energy")
excess_chemical_potential_interpolated = interpolate.interp1d(ionic_strength.m_as("1/reduced_length**3"), 
                                                                                excess_chemical_potential.m_as("reduced_energy"))
activity_coefficient_monovalent_pair = lambda x: np.exp(excess_chemical_potential_interpolated(x.to('1/(reduced_length**3 * N_A)').magnitude))
pka_set = {BeadType: {"pka_value": args.pKa, 
                      "acidity": "acidic"}}

grxmc, labels, ionic_strength_res = pmb.setup_grxmc_reactions(pH_res=args.pH_res, 
                                                              c_salt_res=c_salt_res, 
                                                              proton_name=proton_name, 
                                                              hydroxide_name=hydroxide_name, 
                                                              salt_cation_name=sodium_name, 
                                                              salt_anion_name=chloride_name, 
                                                              activity_coefficient=activity_coefficient_monovalent_pair, 
                                                              pka_set=pka_set)

# Setup espresso to track the ionization of the acid groups
type_map = pmb.get_type_map()
types = list(type_map.values())
espresso_system.setup_type_map(type_list = types)

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
grxmc.set_non_interacting_type (type=non_interacting_type)

for i in tqdm.trange(100):
    espresso_system.integrator.run(steps=1000)
    do_reaction(grxmc,1000)

setup_electrostatic_interactions(units=pmb.units,
                                 espresso_system=espresso_system,
                                 kT=pmb.kT,
                                 solvent_permittivity=solvent_permittivity)

if mode == "long-run":
    N_warmup_loops = 1000
else:
    N_warmup_loops = 100

print("*** Running warmup with electrostatics... ***")
for i in tqdm.trange(N_warmup_loops):
    espresso_system.integrator.run(steps=1000)
    do_reaction(grxmc,100)

# Main loop
print("*** Starting production run... ***")

labels_obs=["time", "alpha", "pressure"]
time_series={}

for label in labels_obs:
    time_series[label]=[]

if mode == "long-run":
    N_production_loops = 5000
elif mode == "test":
    N_production_loops = 500
else:
    N_production_loops = 100

for i in tqdm.trange(N_production_loops):
    espresso_system.integrator.run(steps=1000)
    do_reaction(grxmc,200)

    # Measure time
    time_series["time"].append(espresso_system.time)

    # Measure degree of ionization
    time_series["alpha"].append(len(espresso_system.part.select(type=2,q=-1))/(16*args.mpc))
    time_series["pressure"].append(espresso_system.analysis.pressure()["total"])

inputs={"csalt": args.csalt_res,
        "Lfraction": args.L_fraction,
        "pH": args.pH_res,
        "pKa": args.pKa}

data_path = args.output
data_path.mkdir(parents=True, exist_ok=True)

time_series=pd.DataFrame(time_series)
filename=analysis.built_output_name(input_dict=inputs)

time_series.to_csv(data_path / f"{filename}_time_series.csv", index=False)

