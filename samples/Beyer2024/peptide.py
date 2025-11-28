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
from pathlib import Path
import espressomd
import pandas as pd
import argparse
import tqdm

# Import pyMBE
import pyMBE
from pyMBE.lib import analysis
from pyMBE.lib import handy_functions as hf
from pyMBE.lib.handy_functions import do_reaction

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)

parser = argparse.ArgumentParser(description='Script to run the peptide test cases for pyMBE')
parser.add_argument('--sequence',
                    type=str,
                    required= True,
                    help='sequence of the peptide')
parser.add_argument('--pH',
                    type=float,
                    required= True,
                    help='pH of the solution')
parser.add_argument('--mode',
                    type=str,
                    default= "short-run",
                    choices=["short-run","long-run", "test"],
                    help='sets for how long the simulation runs')
parser.add_argument('--output',
                    type=Path,
                    required= False,
                    default=Path(__file__).parent / "Beyer2024" / "time_series" / "peptides",
                    help='output directory')
parser.add_argument('--no_verbose', action='store_false', help="Switch to deactivate verbose",default=True)
args = parser.parse_args()

# Inputs
sequence=args.sequence
pH=args.pH
inputs={"pH": args.pH,
        "sequence": args.sequence}
mode=args.mode
verbose=args.no_verbose


langevin_seed = 100
dt = 0.01
solvent_permitivity = 78.3

# Sanity check
Lunkad_test_sequences=["E"*5+"H"*5,"K"*5+"D"*5]
Blanco_test_sequence=["nDSHAKRHHGYKRKFHEKHHSHRGYc"]

valid_sequences=Lunkad_test_sequences+Blanco_test_sequence

if sequence not in valid_sequences:
    raise ValueError(f"ERROR: the only valid peptide sequence for this test script are {valid_sequences}")

if sequence in Lunkad_test_sequences:
    path_to_interactions=pmb.root / "parameters" / "peptides" / "Lunkad2021.json"
    path_to_pka=pmb.root / "parameters" / "pka_sets" / "CRC1991.json"
    pmb.load_interaction_parameters(filename=path_to_interactions)
    pmb.load_pka_set(filename=path_to_pka)
    model = '2beadAA'  # Model with 2 beads per each aminoacid
    N_peptide_chains = 4
    sigma=1*pmb.units.Quantity("reduced_length")
    offset_cation=0*pmb.units.Quantity("reduced_length")
    offset_anion=0*pmb.units.Quantity("reduced_length")
    c_salt=1e-2 * pmb.units.mol/ pmb.units.L
    chain_length=len(sequence)*2

elif sequence in Blanco_test_sequence:
    pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)
    pmb.load_interaction_parameters (pmb.root / "parameters" / "peptides" / "Blanco2021.json")
    pmb.load_pka_set (pmb.root / "parameters" / "pka_sets" / "Nozaki1967.json")
    model = '1beadAA'
    N_peptide_chains = 1
    c_salt = 5e-3 * pmb.units.mol/ pmb.units.L
    sigma=1*pmb.units.Quantity("reduced_length")
    offset_cation=0.2*pmb.units.nm-sigma
    offset_anion=0.36*pmb.units.nm-sigma
    chain_length=len(sequence)

pep_concentration = 5.56e-4 *pmb.units.mol/pmb.units.L 

# Simulation parameters
if mode == "short-run":
    Nsamples = 1000
    MD_steps_per_sample = 1000
elif mode == "long-run":
    Nsamples = 5000
    MD_steps_per_sample = 5000
elif mode == "test":
    Nsamples = 500
    MD_steps_per_sample = 700
    c_salt = 5e-3 * pmb.units.mol/ pmb.units.L
    N_peptide_chains = 1
else:
    raise RuntimeError()


pmb.define_peptide (name=sequence, 
                    sequence=sequence, 
                    model=model)

# Solution parameters
cation_name = 'Na'
anion_name = 'Cl'
c_salt=5e-3 * pmb.units.mol/ pmb.units.L

pmb.define_particle(name=cation_name,
                    z=1,
                    sigma=sigma,
                    epsilon=1*pmb.units('reduced_energy'),
                    offset=offset_cation)

pmb.define_particle(name=anion_name,
                    z=-1,
                    sigma=sigma,
                    epsilon=1*pmb.units('reduced_energy'),
                    offset=offset_anion)

# System parameters
volume = N_peptide_chains/(pmb.N_A*pep_concentration)
L = volume ** (1./3.) # Side of the simulation box

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)
espresso_system.time_step=dt
espresso_system.cell_system.skin=0.4
# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)

# Create your molecules into the espresso system
pmb.create_molecule(name=sequence,
                    number_of_molecules=N_peptide_chains,
                    espresso_system=espresso_system)

# Create counterions for the peptide chains
pmb.create_counterions(object_name=sequence,
                    cation_name=cation_name,
                    anion_name=anion_name,
                    espresso_system=espresso_system)

c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,
                     cation_name=cation_name,
                     anion_name=anion_name,
                     c_salt=c_salt)

cpH, labels = pmb.setup_cpH(counter_ion=cation_name,
                                                constant_pH=pH)

if verbose:
    print(f"The box length of your system is {L.to('reduced_length')} = {L.to('nm')}")
    print(f"The acid-base reaction has been successfully setup for {labels}")

# Setup espresso to track the ionization of the acid/basic groups in peptide
type_map =pmb.get_type_map()
espresso_system.setup_type_map(type_list = list(type_map.values()))

# Setup the non-interacting type for speeding up the sampling of the reactions
non_interacting_type = max(type_map.values())+1
cpH.set_non_interacting_type (type=non_interacting_type)
if verbose:
    print(f"The non-interacting type is set to {non_interacting_type}")

#Setup the potential energy
pmb.setup_lj_interactions (espresso_system=espresso_system)
hf.relax_espresso_system(espresso_system=espresso_system,
                          seed=langevin_seed)
hf.setup_electrostatic_interactions(units=pmb.units,
                                    espresso_system=espresso_system,
                                    kT=pmb.kT,
                                    verbose=verbose)
hf.relax_espresso_system(espresso_system=espresso_system,
                          seed=langevin_seed)


hf.setup_langevin_dynamics(espresso_system=espresso_system,
                            kT = pmb.kT,
                            seed = langevin_seed,
                            time_step=dt,
                            tune_skin=False)

espresso_system.cell_system.skin=0.4

# Main loop

labels_obs=["time","charge","rg"]
time_series={}

for label in labels_obs:
    time_series[label]=[]

for sample in tqdm.trange(Nsamples,disable=not verbose):
    # Run LD
    espresso_system.integrator.run(steps=MD_steps_per_sample)
    # Run MC
    do_reaction(cpH, steps=len(sequence))
    # Sample observables
    charge_dict=pmb.calculate_net_charge(espresso_system=espresso_system,
                                            molecule_name=sequence,
                                            dimensionless=True)

    Rg = espresso_system.analysis.calc_rg(chain_start=0,
                                        number_of_chains=N_peptide_chains,
                                        chain_length=chain_length)
    # Store observables
    time_series["time"].append(espresso_system.time)
    time_series["charge"].append(charge_dict["mean"])
    time_series["rg"].append(Rg[0])

data_path=args.output
data_path.mkdir(parents=True, exist_ok=True)

time_series=pd.DataFrame(time_series)
filename=analysis.built_output_name(input_dict=inputs)

time_series.to_csv(data_path / f"{filename}_time_series.csv", index=False)

if verbose:
    print("*** DONE ***")
