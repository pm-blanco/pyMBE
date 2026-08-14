#
# Copyright (C) 2026 pyMBE-dev team
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

from collections import Counter
import numpy as np
import unittest as ut
import pyMBE
from pyMBE.lib.lattice import DiamondLattice
import espressomd
box_l=[10] * 3
espresso_system = espressomd.System(box_l=box_l)


def build_simple_hydrogel_with_optional_angles(junction_angle_mode="none", mpc_local=4):
    """
    Build a simple hydrogel used to validate hydrogel-specific angle support.

    Defines two particle types: chain beads (C) and crosslinker nodes (N),
    connected in a diamond lattice. Always defines the chain-internal C-C-C
    angular potential template. Crosslinker-adjacent templates are added
    depending on junction_angle_mode.

    Args:
        junction_angle_mode (str): Controls which crosslinker angular potential
            templates are defined. Accepted values:
            - ``"none"``: no crosslinker angles defined (default).
            - ``"partial"``: only C-N-C defined (missing N-C-C, triggers error
              on hydrogel creation with gen_angle=True).
            - ``"full"``: both C-N-C and C-C-N defined.
        mpc_local (int): Number of monomers per chain in the diamond lattice.
            Defaults to 4.

    Returns:
        tuple:
            - pmb_local (pyMBE.pymbe_library): pyMBE instance with all
              templates and topology defined but hydrogel not yet created.
            - espresso_system (espressomd.system.System): ESPResSo system
              where the hydrogel will be created.
            - hydrogel_name_local (str): Name under which the hydrogel is
              registered in pmb_local (``"simple_hydrogel"``).
            - diamond_lattice_local (DiamondLattice): The diamond lattice
              instance used to build the topology.
    """
    pmb_local = pyMBE.pymbe_library(seed=7)
    node_type = "N"
    bead_type = "C"
    residue_name = "Res"
    molecule_name = "simple_chain"
    bond_length = 0.355 * pmb_local.units.nm
    harmonic_bond_parameters = {
        "r_0": bond_length,
        "k": 400 * pmb_local.units("reduced_energy / reduced_length**2"),
    }

    pmb_local.define_particle(name=node_type,
                              sigma=bond_length,
                              epsilon=1 * pmb_local.units("reduced_energy"))
    pmb_local.define_particle(name=bead_type,
                              sigma=bond_length,
                              epsilon=1 * pmb_local.units("reduced_energy"))
    pmb_local.define_residue(name=residue_name,
                             central_bead=bead_type,
                             side_chains=[])
    pmb_local.define_molecule(name=molecule_name,
                              residue_list=[residue_name] * mpc_local)
    pmb_local.define_bond(bond_type="harmonic",
                          bond_parameters=harmonic_bond_parameters,
                          particle_pairs=[[bead_type, bead_type],
                                          [bead_type, node_type]])

    angle_parameters = {
        "k": 5 * pmb_local.units("reduced_energy"),
        "phi_0": np.pi * pmb_local.units(""),
    }
    pmb_local.define_angular_potential(angle_type="harmonic",
                           angle_parameters=angle_parameters,
                           particle_triplets=[(bead_type, bead_type, bead_type)])

    if junction_angle_mode in {"partial", "full"}:
        particle_triplets = [(bead_type, node_type, bead_type)]
        if junction_angle_mode == "full":
            particle_triplets.append((node_type, bead_type, bead_type))
        pmb_local.define_angular_potential(angle_type="harmonic",
                               angle_parameters=angle_parameters,
                               particle_triplets=particle_triplets)

    diamond_lattice_local = DiamondLattice(mpc_local, bond_length)
    espresso_system.box_l = [diamond_lattice_local.box_l] * 3
    pmb_local.initialize_lattice_builder(diamond_lattice_local)

    node_topology_local = [
        {"particle_name": node_type, "lattice_index": index}
        for index in diamond_lattice_local.indices
    ]
    chain_topology_local = []
    for node_connectivity in diamond_lattice_local.connectivity:
        node_start = str(diamond_lattice_local.indices[node_connectivity[0]])
        node_end = str(diamond_lattice_local.indices[node_connectivity[1]])
        chain_topology_local.append({
            "node_start": node_start,
            "node_end": node_end,
            "molecule_name": molecule_name,
        })

    hydrogel_name_local = "simple_hydrogel"
    pmb_local.define_hydrogel(hydrogel_name_local,
                              node_topology_local,
                              chain_topology_local)
    return pmb_local, espresso_system, hydrogel_name_local, diamond_lattice_local


def get_angle_counts(pmb_local):
    """
    Return counts of angular potential instances grouped by canonical angle name.

    Args:
        pmb_local (pyMBE.pymbe_library): pyMBE instance whose database will
            be queried for angular potential instances.

    Returns:
        collections.Counter: Maps each canonical angle name (str) to the
            number of instances of that angular potential in the system.
    """
    return Counter(
        angle.name
        for angle in pmb_local.db.get_instances(pmb_type="angle").values()
    )


def expected_crosslinker_angle_counts(diamond_lattice_local):
    """
    Return the expected angular potential instance counts for the simple hydrogel
    when all crosslinker templates are defined (junction_angle_mode="full").

    Args:
        diamond_lattice_local (DiamondLattice): The diamond lattice instance
            used to build the hydrogel, providing mpc and connectivity info.

    Returns:
        dict: Maps canonical angle name (str) to the expected instance count
            (int). Keys are ``"C-C-C"``, ``"C-C-N"``, and ``"C-N-C"``.
    """
    number_of_nodes = 8
    number_of_chains = 16
    crosslinker_functionality = 4

    return {
        "C-C-C": number_of_chains * (diamond_lattice_local.mpc - 2),
        "C-C-N": 2 * number_of_chains,
        "C-N-C": number_of_nodes * (crosslinker_functionality * (crosslinker_functionality - 1) // 2),
    }


class Test(ut.TestCase):
    def setUp(self):
        espresso_system.part.clear()
        espresso_system.bonded_inter.clear()
        espresso_system.non_bonded_inter.reset()

    def test_hydrogel_gen_angle_defaults_to_false(self):
        """
        Hydrogel creation should not generate angles unless gen_angle is requested.
        """
        pmb_local, _, hydrogel_name_local, _ = build_simple_hydrogel_with_optional_angles(
            junction_angle_mode="full"
        )
        pmb_local.create_hydrogel(hydrogel_name_local, box_l)
        self.assertEqual(len(pmb_local.db.get_instances(pmb_type="angle")), 0)

    def test_hydrogel_chain_angles_are_created_without_crosslinker_templates(self):
        """
        If only chain angles are defined, hydrogel creation should still succeed
        and create only the chain-internal angles.
        """
        pmb_local, espresso_system_local, hydrogel_name_local, diamond_lattice_local = build_simple_hydrogel_with_optional_angles(
            junction_angle_mode="none"
        )
        with self.assertLogs(level="WARNING") as log_context:
            pmb_local.create_hydrogel(hydrogel_name_local, box_l, gen_angle=True)

        angle_counts = get_angle_counts(pmb_local)
        pmb_local.set_simulation_engine(espresso_system_local)
        pmb_local.add_instances_to_engine()
        self.assertEqual(
            angle_counts,
            {"C-C-C": len(diamond_lattice_local.connectivity) * (diamond_lattice_local.mpc - 2)},
        )
        self.assertTrue(
            any("No angle templates defined for hydrogel crosslinkers" in message for message in log_context.output)
        )

    def test_hydrogel_junction_angle_counts_are_correct_when_defined(self):
        """
        Hydrogel creation should add the correct number of junction angles when
        the explicit templates exist.
        """
        pmb_local, espresso_system_local, hydrogel_name_local, diamond_lattice_local = build_simple_hydrogel_with_optional_angles(
            junction_angle_mode="full"
        )
        pmb_local.create_hydrogel(hydrogel_name_local, box_l, gen_angle=True)
        pmb_local.set_simulation_engine(espresso_system_local)
        pmb_local.add_instances_to_engine()
        self.assertEqual(get_angle_counts(pmb_local),
                         expected_crosslinker_angle_counts(diamond_lattice_local))

    def test_hydrogel_partial_crosslinker_angle_definitions_raise(self):
        """
        Defining only a subset of required crosslinker-adjacent angles should fail.
        """
        pmb_local, _, hydrogel_name_local, _ = build_simple_hydrogel_with_optional_angles(
            junction_angle_mode="partial"
        )
        with self.assertRaises(ValueError):
            pmb_local.create_hydrogel(hydrogel_name_local, box_l, gen_angle=True)


if __name__ == "__main__":
    ut.main()
