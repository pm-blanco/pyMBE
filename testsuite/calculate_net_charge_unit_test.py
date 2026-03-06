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
import unittest as ut
from pyMBE.lib.lattice import DiamondLattice

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(seed=42)

pmb.define_particle(name='0P',
                    z=0,
                    sigma=1*pmb.units.reduced_length,
                    epsilon=1*pmb.units.reduced_energy)

pmb.define_particle(name='+1p',
                    z=+1,
                    sigma=1*pmb.units.reduced_length,
                    epsilon=1*pmb.units.reduced_energy)

pmb.define_particle(name='-1p',
                    z=-1,
                    sigma=1*pmb.units.reduced_length,
                    epsilon=1*pmb.units.reduced_energy)

pmb.define_residue(
    name = 'R1',
    central_bead = '+1p',
    side_chains = ['0P']
    )

pmb.define_residue(
    name = 'R2',
    central_bead = '-1p',
    side_chains = ['R1']
    )

bond_type = 'harmonic'
generic_bond_length=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

harmonic_bond = {'r_0'    : generic_bond_length,
                 'k'      : generic_harmonic_constant,
                }

pmb.define_default_bond(bond_type = bond_type, 
                        bond_parameters = harmonic_bond)

molecule_name = 'generic_molecule'
pmb.define_molecule(name=molecule_name, 
                    residue_list = ['R1']*2+['R2']*3)

# Create an instance of an espresso system
espresso_system=espressomd.System(box_l = [10]*3)

# Create your molecules into the espresso system

diamond_lattice = DiamondLattice(4, 3.5 * pmb.units.reduced_length)
lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)
indices = diamond_lattice.indices
node_topology = []
for index in range(len(indices)):
    node_topology.append({"particle_name": "+1p",
                        "lattice_index": indices[index]})
node_labels = lattice_builder.node_labels
chain_labels = lattice_builder.chain_labels
reverse_node_labels = {v: k for k, v in node_labels.items()}
chain_topology = []
for chain_data in chain_labels.items():
    node_label_pair = chain_data[0]
    node_label_s, node_label_e = [int(x) for x in node_label_pair.strip("()").split(",")]
    chain_topology.append({'node_start': reverse_node_labels[node_label_s],
                            'node_end': reverse_node_labels[node_label_e],
                            'molecule_name': molecule_name})
pmb.define_hydrogel("my_hydrogel", 
                    node_topology, 
                    chain_topology)
pmb.create_hydrogel(name="my_hydrogel",
                    espresso_system=espresso_system,
                    use_default_bond=True)

class Test(ut.TestCase):
    def test_calculate_net_charge_with_units(self):
        """
        *** Unit test: check that calculate_net_charge calculates the charge in a hydrogel properly with units
        """
        
        # Check that it calculates properly the charge of the whole hydrogel
        charge_map = pmb.calculate_net_charge(object_name="my_hydrogel",
                                        pmb_type="hydrogel",
                                        espresso_system=espresso_system)

        np.testing.assert_equal(charge_map["mean"], 40.0*pmb.units.Quantity(1,'reduced_charge'))
        np.testing.assert_equal(charge_map["instances"], {0: 40.0*pmb.units.Quantity(1,'reduced_charge')})

        # Check that it calculates properly the charge of the chains in the hydrogel
        charge_map = pmb.calculate_net_charge(object_name=molecule_name,
                                              pmb_type="molecule",
                                              espresso_system=espresso_system)
        # Check mean charge
        np.testing.assert_equal(charge_map["mean"], 2.0*pmb.units.Quantity(1,'reduced_charge'))
        # Check molecule charge map
        np.testing.assert_equal(charge_map["instances"],
                                {0: 2.0*pmb.units.Quantity(1,'reduced_charge'), 
                                 1: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                 2: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                 3: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                 4: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                 5: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                 6: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                 7: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                 8: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                 9: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                10: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                11: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                12: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                13: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                14: 2.0*pmb.units.Quantity(1,'reduced_charge'),
                                15: 2.0*pmb.units.Quantity(1,'reduced_charge')})

        # Check that it calculates properly the charge of the residues in the hydrogel
        charge_map_r1 = pmb.calculate_net_charge(object_name="R1",
                                                 pmb_type="residue",
                                                 espresso_system=espresso_system)
        charge_map_r2 = pmb.calculate_net_charge(object_name="R2",
                                                 pmb_type="residue",
                                                 espresso_system=espresso_system)
        res_charge_map = charge_map_r1["instances"] | charge_map_r2["instances"]
        np.testing.assert_equal(res_charge_map[0], 1.0*pmb.units.Quantity(1,'reduced_charge'))
        np.testing.assert_equal(res_charge_map[1], 1.0*pmb.units.Quantity(1,'reduced_charge'))
        np.testing.assert_equal(res_charge_map[2], 0.0*pmb.units.Quantity(1,'reduced_charge'))
        np.testing.assert_equal(res_charge_map[3], 0.0*pmb.units.Quantity(1,'reduced_charge'))
        np.testing.assert_equal(res_charge_map[4], 0.0*pmb.units.Quantity(1,'reduced_charge'))
        # Check that particle-level charge is computed per particle instance
        charge_map_particle = pmb.calculate_net_charge(object_name="+1p",
                                                       pmb_type="particle",
                                                       espresso_system=espresso_system)
        np.testing.assert_equal(charge_map_particle["mean"], 1.0*pmb.units.Quantity(1,'reduced_charge'))
        particle_ids = pmb.get_particle_id_map(object_name="+1p")["all"]
        expected_particle_map = {pid: 1.0*pmb.units.Quantity(1,'reduced_charge') for pid in particle_ids}
        np.testing.assert_equal(charge_map_particle["instances"], expected_particle_map)

        
    def test_calculate_net_charge_without_units(self):
        """
        *** Unit test: check that calculate_net_charge calculates the charge in a molecule properly without units
        """
        # Check that it calculates properly the charge of the whole hydrogel
        charge_map = pmb.calculate_net_charge(object_name="my_hydrogel",
                                            pmb_type="hydrogel",
                                            espresso_system=espresso_system,
                                            dimensionless=True)
        np.testing.assert_equal(charge_map["mean"], 40.0)
        np.testing.assert_equal(charge_map["instances"], {0: 40.0})
        # Check the case where the returned charge does not have a dimension
        charge_map = pmb.calculate_net_charge(object_name=molecule_name,
                                              pmb_type="molecule",
                                              espresso_system=espresso_system,
                                              dimensionless=True)
        # Check mean charge
        np.testing.assert_equal(charge_map["mean"], 2.0)
        # Check molecule charge map
        np.testing.assert_equal(charge_map["instances"],
                                {0: 2.0, 
                                 1: 2.0,
                                 2: 2.0,
                                 3: 2.0,
                                 4: 2.0,
                                 5: 2.0,
                                 6: 2.0,
                                 7: 2.0,
                                 8: 2.0,
                                 9: 2.0,
                                10: 2.0,
                                11: 2.0,
                                12: 2.0,
                                13: 2.0,
                                14: 2.0,
                                15: 2.0})
        charge_map_r1 = pmb.calculate_net_charge(object_name="R1",
                                                 pmb_type="residue",
                                                 espresso_system=espresso_system,
                                                 dimensionless=True)
        charge_map_r2 = pmb.calculate_net_charge(object_name="R2",
                                                 pmb_type="residue",
                                                 espresso_system=espresso_system,
                                                 dimensionless=True)
        res_charge_map = charge_map_r1["instances"] | charge_map_r2["instances"]
        np.testing.assert_equal(res_charge_map[0], 1.0)
        np.testing.assert_equal(res_charge_map[1], 1.0)
        np.testing.assert_equal(res_charge_map[2], 0.0)
        np.testing.assert_equal(res_charge_map[3], 0.0)
        np.testing.assert_equal(res_charge_map[4], 0.0)
        # Check that particle-level charge is computed per particle instance
        charge_map_particle = pmb.calculate_net_charge(object_name="+1p",
                                                       pmb_type="particle",
                                                       espresso_system=espresso_system,
                                                       dimensionless=True)
        np.testing.assert_equal(charge_map_particle["mean"], 1.0)
        particle_ids = pmb.get_particle_id_map(object_name="+1p")["all"]
        expected_particle_map = {pid: 1.0 for pid in particle_ids}
        np.testing.assert_equal(charge_map_particle["instances"], expected_particle_map)

if __name__ == '__main__':
    ut.main()
