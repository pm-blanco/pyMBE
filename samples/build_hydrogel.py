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
#
# Builds a diamond-lattice hydrogel (bead types A, C, D) with angular potentials.
#
# By default only the chain-internal angles A-A-A and A-A-C is defined.
# Use --include_crosslinker_angles to also add A-C-D and C-D-C templates.
# Use --visualize to display a 3D plot of the lattice.
# Use --mpc to set the number of monomers per chain (default: 5).

import argparse
import numpy as np
import espressomd
import matplotlib.pyplot as plt

import pyMBE
from pyMBE.lib.lattice import DiamondLattice


def define_hydrogel_with_angular_potential(espresso_system, include_crosslinker_angles, mpc=5):
    """
    Build a hydrogel with angular potentials on the diamond lattice.

    Defines three bead types: crosslinker nodes (D), chain-end beads (C), and
    internal chain beads (A). Always defines the chain-internal A-A-A and A-A-C
    angular potential templates. When include_crosslinker_angles is True, also
    defines the D-C-A (canonical: A-C-D) and C-D-C templates.

    Args:
        espresso_system (espressomd.system.System): ESPResSo system object in
            which the hydrogel particles, bonds, and angles will be created.
        include_crosslinker_angles (bool): If True, also define angular potential
            templates for the crosslinker-adjacent triplets D-C-A (canonical:
            A-C-D) and C-D-C in addition to the chain-internal A-A-A and A-A-C
            templates.
        mpc (int): Number of monomers per chain. Must be >= 5 (two terminal C
            beads plus at least three internal A beads to generate A-A-A angles).
            Defaults to 5.

    Returns:
        tuple:
            - pmb (pyMBE.pymbe_library): pyMBE instance holding all templates
              and instances for this hydrogel.
            - hydrogel_id (int): Assembly ID of the created hydrogel, as
              returned by pmb.create_hydrogel.
            - lattice_builder (pyMBE.lib.lattice.LatticeBuilder): Lattice
              builder instance, usable for drawing the lattice.
    """
    pmb = pyMBE.pymbe_library(seed=42)

    node_type = "D"
    internal_bead_type = "A"
    chain_end_bead_type = "C"
    terminal_residue_name = "TerminalRes"
    internal_residue_name = "InternalRes"
    molecule_name = "hydrogel_chain"

    bond_length = 0.355 * pmb.units.nm

    pmb.define_particle(name=node_type,
                        sigma=bond_length,
                        epsilon=1*pmb.units("reduced_energy"))
    pmb.define_particle(name=internal_bead_type,
                        sigma=bond_length,
                        epsilon=1*pmb.units("reduced_energy"))
    pmb.define_particle(name=chain_end_bead_type,
                        sigma=bond_length,
                        epsilon=1*pmb.units("reduced_energy"))

    pmb.define_residue(name=terminal_residue_name,
                       central_bead=chain_end_bead_type,
                       side_chains=[])
    pmb.define_residue(name=internal_residue_name,
                       central_bead=internal_bead_type,
                       side_chains=[])
    pmb.define_molecule(name=molecule_name,
                        residue_list=(
                            [terminal_residue_name]
                            + [internal_residue_name] * (mpc - 2)
                            + [terminal_residue_name]
                        ))

    harmonic_bond_parameters = {
        "r_0": bond_length,
        "k": 400 * pmb.units("reduced_energy / reduced_length**2"),
    }
    pmb.define_bond(bond_type="harmonic",
                    bond_parameters=harmonic_bond_parameters,
                    particle_pairs=[
                        [internal_bead_type, internal_bead_type],
                        [internal_bead_type, chain_end_bead_type],
                        [node_type, chain_end_bead_type],
                    ])

    angle_parameters = {
        "k": 25 * pmb.units("reduced_energy"),
        "phi_0": np.pi * pmb.units(""),
    }
    pmb.define_angular_potential(angle_type="harmonic",
                                 angle_parameters=angle_parameters,
                                 particle_triplets=[
                                     (chain_end_bead_type, internal_bead_type, internal_bead_type),  # A-A-C
                                     (internal_bead_type, internal_bead_type, internal_bead_type),   # A-A-A
                                 ])

    if include_crosslinker_angles:
        pmb.define_angular_potential(angle_type="harmonic",
                                     angle_parameters=angle_parameters,
                                     particle_triplets=[
                                         (node_type, chain_end_bead_type, internal_bead_type),  # canonical: A-C-D
                                         (chain_end_bead_type, node_type, chain_end_bead_type),  # canonical: C-D-C
                                     ])

    diamond_lattice = DiamondLattice(mpc, bond_length)
    espresso_system.box_l = [diamond_lattice.box_l] * 3
    lattice_builder = pmb.initialize_lattice_builder(diamond_lattice)

    node_topology = [
        {"particle_name": node_type, "lattice_index": index}
        for index in diamond_lattice.indices
    ]
    chain_topology = []
    for node_connectivity in diamond_lattice.connectivity:
        node_start = str(diamond_lattice.indices[node_connectivity[0]])
        node_end = str(diamond_lattice.indices[node_connectivity[1]])
        chain_topology.append({
            "node_start": node_start,
            "node_end": node_end,
            "molecule_name": molecule_name,
        })

    lattice_builder.chains = chain_topology
    pmb.define_hydrogel(name="simple_hydrogel",
                        node_map=node_topology,
                        chain_map=chain_topology)
    hydrogel_id = pmb.create_hydrogel(name="simple_hydrogel",
                                      box_l=espresso_system.box_l,
                                      gen_angle=True)
    return pmb, hydrogel_id, lattice_builder


def summarize_angles(pmb):
    """
    Return angle templates, instances, and per-name counts from a pyMBE instance.

    Args:
        pmb (pyMBE.pymbe_library): pyMBE instance whose database will be queried.

    Returns:
        tuple:
            - angle_templates_df (pandas.DataFrame): DataFrame of all defined
              angle templates.
            - angle_instances_df (pandas.DataFrame): DataFrame of all angle
              instances currently in the system.
            - angle_counts (pandas.Series): Count of angle instances grouped
              by canonical angle name, sorted alphabetically by name.
    """
    angle_templates_df = pmb.get_templates_df(pmb_type="angle")
    angle_instances_df = pmb.get_instances_df(pmb_type="angle")
    angle_counts = angle_instances_df["name"].value_counts().sort_index()
    return angle_templates_df, angle_instances_df, angle_counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build a diamond-lattice hydrogel (A/C/D beads) with angular potentials."
    )
    parser.add_argument(
        "--mpc",
        type=int,
        default=5,
        help="Number of monomers per chain (default: 5).",
    )
    parser.add_argument(
        "--include_crosslinker_angles",
        action="store_true",
        help="Also define crosslinker-adjacent angle templates (A-C-D and C-D-C).",
    )
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="Display a 3D plot of the hydrogel lattice after creation.",
    )
    args = parser.parse_args()

    espresso_system = espressomd.System(box_l=[1.0] * 3)

    pmb, hydrogel_id, lattice_builder = define_hydrogel_with_angular_potential(
        espresso_system, args.include_crosslinker_angles, mpc=args.mpc
    )

    if args.visualize:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        lattice_builder.draw_lattice(ax=ax, pmb=pmb)
        lattice_builder.draw_simulation_box(ax)
        plt.legend(fontsize=12)
        plt.show()

    angle_templates_df, angle_instances_df, angle_counts = summarize_angles(pmb)
    print(f"\nHydrogel assembly id: {hydrogel_id}")
    print("\nAngle templates:")
    print(angle_templates_df)
    print("\nAngle instance dataframe:")
    print(angle_instances_df)
    print("\nAngle instance counts:")
    print(angle_counts)

    expected_angle_names = {"A-A-A", "A-A-C", "A-C-D", "C-D-C"} if args.include_crosslinker_angles else {"A-A-A", "A-A-C"}
    actual_angle_names = set(angle_counts.index)
    if actual_angle_names != expected_angle_names:
        raise AssertionError(
            f"Expected angle names {sorted(expected_angle_names)}, "
            f"got {sorted(actual_angle_names)}"
        )
    print("\nOK: all expected angular potential triplets generated.")
