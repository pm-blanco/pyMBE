#
# Copyright (C) 2024-2025 pyMBE-dev team
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

import pyMBE
import espressomd
import matplotlib.pyplot as plt
from pyMBE.lib.lattice import DiamondLattice

pmb = pyMBE.pymbe_library(seed=42)
reduced_unit_set = pmb.get_reduced_units()
# Monomers per chain
mpc = 40
# Define node particle
NodeType = "node_type"
pmb.define_particle(name=NodeType, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
# define monomers
BeadType1 = "C"
pmb.define_particle(name=BeadType1, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))
BeadType2 = "M"
pmb.define_particle(name=BeadType2, sigma=0.355*pmb.units.nm, epsilon=1*pmb.units('reduced_energy'))

Res1 = "res_1"
pmb.define_residue(
    name=Res1,  # Name of the residue
    central_bead=BeadType1,  # Define the central bead name
    side_chains=[]  # Assuming no side chains for the monomer
)

Res2 = "res_2"
pmb.define_residue(
    name=Res2,  # Name of the residue
    central_bead=BeadType2,  # Define the central bead name
    side_chains=[]  # Assuming no side chains for the monomer
)

residue_list = [Res1]*(mpc//2) + [Res2]*(mpc//2)

# Defining bonds in the hydrogel for all different pairs
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
generic_bond_l = 0.355*pmb.units.nm
HARMONIC_parameters = {'r_0'    : generic_bond_l,
                       'k'      : generic_harmonic_constant}
pmb.define_bond(bond_type = 'harmonic',
                        bond_parameters = HARMONIC_parameters, particle_pairs = [[BeadType1, BeadType1],
                                                                                [BeadType1, BeadType2],
                                                                                [BeadType2, BeadType2],
                                                                                [NodeType, BeadType1],
                                                                                [NodeType, BeadType2]])
# Provide mpc and bond_l to Diamond Lattice
diamond_lattice = DiamondLattice(mpc, generic_bond_l)
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
connectivity = diamond_lattice.connectivity
node_labels = lattice_builder.node_labels
reverse_node_labels = {v: k for k, v in node_labels.items()}
connectivity_with_labels = {(reverse_node_labels[i], reverse_node_labels[j]) for i, j in connectivity}
chain_topology = []

for node_s, node_e in connectivity_with_labels:
    chain_topology.append({'node_start':node_s,
                              'node_end': node_e,
                              'residue_list':residue_list})


pmb.define_hydrogel("my_hydrogel",node_topology, chain_topology)
hydrogel_info = pmb.create_hydrogel("my_hydrogel", espresso_system)

fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")
lattice_builder.draw_lattice(ax)
lattice_builder.draw_simulation_box(ax)
plt.legend(fontsize=12)
plt.show()
