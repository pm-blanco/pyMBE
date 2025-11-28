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

import numpy as np 
import espressomd
import pyMBE
import logging
import io

# Create an in-memory log stream
log_stream = io.StringIO()
logging.basicConfig(level=logging.INFO, 
                    format="%(levelname)s: %(message)s",
                    handlers=[logging.StreamHandler(log_stream)])
espresso_system = espressomd.System(box_l = [100]*3)

def build_peptide_in_espresso(seed):
    pmb = pyMBE.pymbe_library(seed=seed)

    # Simulation parameters
    pmb.set_reduced_units(unit_length=0.4*pmb.units.nm)

    # Peptide parameters
    sequence = 'EEEEEEE'
    model = '2beadAA'  # Model with 2 beads per each aminoacid

    # Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.
    path_to_interactions=pmb.root / "parameters" / "peptides" / "Lunkad2021.json"
    path_to_pka=pmb.root / "parameters" / "pka_sets" / "CRC1991.json"
    pmb.load_interaction_parameters(filename=path_to_interactions) 
    pmb.load_pka_set(path_to_pka)

    # Defines the peptide in the pyMBE data frame
    peptide_name = 'generic_peptide'
    pmb.define_peptide(name=peptide_name, 
                       sequence=sequence, 
                       model=model)

    # Bond parameters
    generic_bond_length=0.4 * pmb.units.nm
    generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

    HARMONIC_parameters = {'r_0'    : generic_bond_length,
                           'k'      : generic_harmonic_constant}

    pmb.define_default_bond(bond_type = 'harmonic',
                            bond_parameters = HARMONIC_parameters)

    # Add all bonds to espresso system
    pmb.add_bonds_to_espresso(espresso_system=espresso_system)

    # Create molecule in the espresso system
    pmb.create_molecule(name=peptide_name, 
                        number_of_molecules=1, 
                        espresso_system=espresso_system, 
                        use_default_bond=True)

    # Extract positions of particles in the peptide
    positions = []
    molecule_id = pmb.df.loc[pmb.df['name']==peptide_name].molecule_id.values[0]
    particle_id_list = pmb.df.loc[pmb.df['molecule_id']==molecule_id].particle_id.dropna().to_list()
    for pid in particle_id_list:
        positions.append(espresso_system.part.by_id(pid).pos)

    return np.asarray(positions)


print("*** Check that the using the same seed results in the same initial particle positions***")
positions1 = build_peptide_in_espresso(42)
positions2 = build_peptide_in_espresso(42)

np.testing.assert_almost_equal(positions1, positions2)

print("*** Test passed ***")
