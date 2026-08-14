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
#
# Example: demonstrating angle potential support in pyMBE.
#
# This script shows three usages of `gen_angle=True`:
#   1. create_residue with gen_angle=True
#        - builds a single residue with auto-generated angles
#   2. create_residue with a nested residue as side chain
#        - builds a residue whose side chain is itself a residue; verifies
#          that all bonded triplets (intra- and cross-level) are generated
#   3. create_molecule with gen_angle=True
#        - builds a chain of residues; angles are auto-generated across the
#          full molecule, including angles that span residue boundaries


import pyMBE
import numpy as np
import espressomd

# ----------------------------------------------------------------------
# Set up the ESPResSo system and pyMBE
# ----------------------------------------------------------------------
box_l=[20] * 3
espresso_system = espressomd.System(box_l=box_l)
pmb = pyMBE.pymbe_library(seed=42)

# ----------------------------------------------------------------------
# Particle definitions
# ----------------------------------------------------------------------
pmb.define_particle(name='A',
                    z=0,
                    sigma=0.4 * pmb.units.nm,
                    epsilon=1 * pmb.units('reduced_energy'))

pmb.define_particle(name='B',
                    z=0,
                    sigma=0.4 * pmb.units.nm,
                    epsilon=1 * pmb.units('reduced_energy'))

pmb.define_particle(name='C',
                    z=0,
                    sigma=0.4 * pmb.units.nm,
                    epsilon=1 * pmb.units('reduced_energy'))

pmb.define_particle(name='D',
                    z=0,
                    sigma=0.4 * pmb.units.nm,
                    epsilon=1 * pmb.units('reduced_energy'))

# ----------------------------------------------------------------------
# Default bond used everywhere
# ----------------------------------------------------------------------
HARMONIC_bond_parameters = {
    'r_0': 0.4 * pmb.units.nm,
    'k': 40 * pmb.units('reduced_energy / reduced_length**2'),
}
pmb.define_default_bond(bond_type='harmonic',
                        bond_parameters=HARMONIC_bond_parameters)

# ----------------------------------------------------------------------
# Default angle used by both demos
# ----------------------------------------------------------------------
default_angle_params = {
    'k': 50 * pmb.units('reduced_energy'),
    'phi_0': np.pi * pmb.units(''),
}
pmb.define_default_angular_potential(angle_type="harmonic",
                         angle_parameters=default_angle_params)


def show_database(label):
    print(f"\n############ {label} ############")
    for pmb_type in ["particle", "residue", "molecule", "bond", "angle"]:
        print(f"\n=== {pmb_type} templates ===")
        print(pmb.get_templates_df(pmb_type=pmb_type))
        print(f"\n=== {pmb_type} instances ===")
        print(pmb.get_instances_df(pmb_type=pmb_type))


# ======================================================================
# Demo 1: create_residue with gen_angle=True
# ======================================================================
# Topology: central A bonded to side particles A, B, C.
# Auto-generated angles (central A has 3 neighbors): A-A-B, A-A-C, B-A-C
pmb.define_residue(name="Res_1",
                   central_bead="A",
                   side_chains=["A", "B", "C"])

residue_id = pmb.create_residue(name="Res_1",
                                box_l=box_l,
                                use_default_bond=True,
                                gen_angle=True)

show_database("Demo 1: single residue with gen_angle=True")

# ----------------------------------------------------------------------
# Clean up: remove the residue (and its particles, bonds, angles) from
# both ESPResSo and the pyMBE database before running Demo 2.
# ----------------------------------------------------------------------
pmb.delete_instances_in_system(instance_id=residue_id,
                               pmb_type="residue")


# ======================================================================
# Demo 2: create_residue with a nested residue as side chain
# ======================================================================
# Topology:  D₁       where SubRes_3 = B bonded to C and D₂
#             \
#              A
#              |
#              B -- C
#              |
#              D₂
#
# SubRes_3 (central B, side chains C and D) is passed as a side chain of
# the outer residue NestedRes_3 (central A, direct side chain D).
# D₁ and D₂ are distinct instances of particle type D.
#
# All bonded triplets and their canonical angle names:
#   centered at A:  D₁ -- A -- B  →  "B-A-D"  (outer {B,D} sorted)
#   centered at B:  A  -- B -- C  →  "A-B-C"
#   centered at B:  A  -- B -- D₂ →  "A-B-D"
#   centered at B:  C  -- B -- D₂ →  "C-B-D"
pmb.define_residue(name="SubRes_3",
                   central_bead="B",
                   side_chains=["C", "D"])

pmb.define_residue(name="NestedRes_3",
                   central_bead="A",
                   side_chains=["D", "SubRes_3"])

nested_residue_id = pmb.create_residue(name="NestedRes_3",
                                       box_l=box_l,
                                       use_default_bond=True,
                                       gen_angle=True)

show_database("Demo 2: nested residue with gen_angle=True")

angle_instances_df = pmb.get_instances_df(pmb_type="angle")
particle_instances_df = pmb.get_instances_df(pmb_type="particle")
print(f"\nAngle instance count: {len(angle_instances_df)}")

id_to_type = dict(zip(particle_instances_df["particle_id"],
                      particle_instances_df["name"]))
actual_triplets = set()
for _, row in angle_instances_df.iterrows():
    p1 = id_to_type[row["particle_id1"]]
    center = id_to_type[row["particle_id2"]]
    p3 = id_to_type[row["particle_id3"]]
    outer = sorted([p1, p3])
    actual_triplets.add(f"{outer[0]}-{center}-{outer[1]}")

print("Resolved angle triplets:", sorted(actual_triplets))
expected_triplets = {"A-B-C", "A-B-D", "B-A-D", "C-B-D"}
assert actual_triplets == expected_triplets, (
    f"Demo 2: expected {sorted(expected_triplets)}, got {sorted(actual_triplets)}"
)
print("Demo 2: all expected angle triplets generated OK")

# ======================================================================
# Demo 3: create_molecule with gen_angle=True
# ======================================================================
# Topology: a chain of 4 Res_2 residues, each with central A and one
# side chain B. After auto-generation, angles include both intra-residue
# (A-A-B) and inter-residue (A-A-A) triplets.
pmb.define_residue(name="Res_2",
                   central_bead="A",
                   side_chains=["B"])

pmb.define_molecule(name="Mol_1",
                    residue_list=["Res_2"] * 4)

molecule_ids = pmb.create_molecule(name="Mol_1",
                                   number_of_molecules=1,
                                   box_l=box_l,
                                   backbone_vector=[1.0, 0.0, 0.0],
                                   use_default_bond=True,
                                   gen_angle=True)

show_database("Demo 3: linear molecule with gen_angle=True")

pmb.delete_instances_in_system(instance_id=molecule_ids[0],
                               pmb_type="molecule"
                               )

