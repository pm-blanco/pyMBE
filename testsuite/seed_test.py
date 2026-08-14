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

import numpy as np 
import espressomd
import pyMBE
from pyMBE.lib import handy_functions as hf
import unittest as ut

box_l=[100]*3
espresso_system = espressomd.System(box_l = box_l)

def build_peptide_in_espresso(seed):
    pmb = pyMBE.pymbe_library(seed=seed)
    # Peptide parameters
    sequence = 'EEEEEEE'
    model = '2beadAA'  # Model with 2 beads per each aminoacid
    # Load peptide parametrization from Lunkad, R. et al.  Molecular Systems Design & Engineering (2021), 6(2), 122-131.
    path_to_interactions=pmb.root / "parameters" / "peptides" / "Lunkad2021"
    path_to_pka=pmb.root / "parameters" / "pka_sets" / "CRC1991.json"
    pmb.load_database(folder=path_to_interactions) 
    pmb.load_pka_set(path_to_pka)
    pka_set = pmb.get_pka_set()
    for particle_name in pka_set.keys():
        pmb.define_monoprototic_particle_states(acidity=pka_set[particle_name]["acidity"],
                                             particle_name=particle_name)
    # define residues
    hf.define_peptide_AA_residues(sequence=sequence,
                                  model=model, 
                                  pmb=pmb)
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
    # Create molecule in the espresso system
    pmb.create_molecule(name=peptide_name, 
                        number_of_molecules=1, 
                        box_l=box_l, 
                        use_default_bond=True)
    
    pmb.set_simulation_engine(espresso_system)
    pmb.add_instances_to_engine()
    # Extract positions of particles in the peptide
    particle_id_list = pmb.get_particle_id_map("generic_peptide")["all"]
    positions = []
    for pid in particle_id_list:
        positions.append(espresso_system.part.by_id(pid).pos)
    pmb.delete_instances_in_system(instance_id=0,
                                   pmb_type="peptide")
    return np.asarray(positions)

class Test(ut.TestCase):
    def test_deterministic_build_pyMBE(self):
        """
        Check that the using the same seed results in the same initial particle positions
        """
        positions1 = build_peptide_in_espresso(42)
        positions2 = build_peptide_in_espresso(42)
        np.testing.assert_equal(positions1, 
                                positions2)

if __name__ == "__main__":
    ut.main()