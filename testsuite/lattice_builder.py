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
import numpy as np
import pyMBE.lib.lattice
import unittest as ut
import matplotlib
import matplotlib.pyplot as plt
import pyMBE
import espressomd
import pint

matplotlib.use("Agg") # use a non-graphic backend


units = pint.UnitRegistry()

mpc = 4
bond_l = 0.355 * units.nm



diamond = pyMBE.lib.lattice.DiamondLattice(mpc, bond_l)
box_l=[diamond.box_l] * 3
espresso_system = espressomd.System(box_l=box_l)


# Define node particle
NodeType1 = "node_type1"
NodeType2 = "node_type2"


# define monomers
BeadType1 = "carbon"
BeadType2 = "oxygen"
BeadType3 = "nitrogen"

Res1 = "res_1"
Res2 = "res_2"
Res3 = "res_3"


def define_templates(pmb):
    # Defining bonds in the hydrogel for all different pairs
    generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')
    generic_bond_l = 0.355*pmb.units.nm
    HARMONIC_parameters = {'r_0'    : generic_bond_l,
                        'k'      : generic_harmonic_constant}
    pmb.define_particle(name=NodeType1, 
                        sigma=0.355*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))
    pmb.define_particle(name=NodeType2, 
                        sigma=0.355*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))
    pmb.define_particle(name=BeadType1, 
                        sigma=0.355*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))
    pmb.define_particle(name=BeadType2, 
                        sigma=0.355*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))
    pmb.define_particle(name=BeadType3, 
                        sigma=0.355*pmb.units.nm, 
                        epsilon=1*pmb.units('reduced_energy'))
    pmb.define_residue(name=Res1,
                       central_bead=BeadType1,
                       side_chains=[])
    pmb.define_residue(name=Res2,
                       central_bead=BeadType2,  
                       side_chains=[])
    pmb.define_residue(name=Res3,
                       central_bead=BeadType3,
                       side_chains=[])
    pmb.define_bond(bond_type = 'harmonic',
                    bond_parameters = HARMONIC_parameters, particle_pairs = [[BeadType1, BeadType1],
                                                                            [BeadType1, BeadType2],
                                                                            [BeadType1, BeadType3],
                                                                            [BeadType2, BeadType2],
                                                                            [BeadType2, BeadType3],
                                                                            [BeadType3, BeadType3],
                                                                            [BeadType1, NodeType1],
                                                                            [BeadType1, NodeType2],
                                                                            [BeadType2, NodeType1],
                                                                            [BeadType2, NodeType2],
                                                                            [BeadType3, NodeType1],
                                                                            [BeadType3, NodeType2]])



class Test(ut.TestCase):
    colormap = {"default_linker":"green",
                "default_monomer":"blue",
                BeadType1: "pink",
                BeadType2: "purple",
                BeadType3: "black",
                Res3: "red",
                NodeType2: "orange",
                NodeType1: "cyan",
                Res1: "yellow",
                Res2: "magenta"}
    
    def test_lattice_setup(self):
        """
        Unit tests for the lattice builder module
        """        
        pmb = pyMBE.pymbe_library(42)
        define_templates(pmb=pmb)
        # --- Invalid low-level operations ---
        with self.assertRaises(ValueError):
            pmb._create_hydrogel_node("[1 1 1]", NodeType1, box_l=box_l)

        with self.assertRaises(ValueError):
            pmb._create_hydrogel_chain(
                "[0 0 0]", "[1 1 1]",
                {0: [0, 0, 0], 1: diamond.box_l / 4.0 * np.ones(3)},
                [diamond.box_l]*3
            )

        # --- Lattice initialization ---
        lattice = pmb.initialize_lattice_builder(diamond)

        assert len(lattice.nodes) == len(diamond.indices)
        assert len(lattice.chains) == 0

        # --- Default chains ---
        lattice.add_default_chains(mpc=2)
        assert len(lattice.chains) == len(diamond.connectivity)

        # --- Default node types ---
        assert lattice.get_node("[1 1 1]") == "default_linker"
        assert lattice.get_node("[0 0 0]") == "default_linker"

        # --- Custom node assignment ---
        lattice.set_node(node="[1 1 1]", residue=NodeType1)
        lattice.set_node(node="[0 0 0]", residue=NodeType2)

        np.testing.assert_equal(lattice.get_node("[1 1 1]"), NodeType1)
        np.testing.assert_equal(lattice.get_node("[0 0 0]"), NodeType2)

        # untouched nodes remain default
        np.testing.assert_equal(lattice.get_node("[2 2 0]"), "default_linker")
        np.testing.assert_equal(lattice.get_node("[3 1 3]"), "default_linker")

        pmb2 = pyMBE.pymbe_library(23)
        define_templates(pmb=pmb2)
        diamond2 = pyMBE.lib.lattice.DiamondLattice(mpc, bond_l)
        lattice = pmb2.initialize_lattice_builder(diamond2)

        sequence = [Res3, Res1, Res2, Res1]
        node_a = "[0 0 0]"
        node_b = "[1 1 1]"

        # 1. Define chain in forward direction
        lattice.set_chain(
            node_start=node_a,
            node_end=node_b,
            sequence=sequence
        )

        np.testing.assert_equal(
            actual=lattice.get_chain(node_a, node_b),
            desired=sequence,
            verbose=True
        )

        # 2. Define chain explicitly in reverse direction
        lattice.set_chain(
            node_start=node_b,
            node_end=node_a,
            sequence=sequence
        )

        np.testing.assert_equal(
            actual=lattice.get_chain(node_b, node_a),
            desired=sequence,
            verbose=True
        )

        # 3. Geometry-safe reversal:
        #    forward lookup returns reversed sequence
        np.testing.assert_equal(
            actual=lattice.get_chain(node_a, node_b),
            desired=sequence[::-1],
            verbose=True
        )

        # 4. Invalid chain lookup
        with self.assertRaises(RuntimeError):
            lattice.get_chain("[1 1 1]", "[2 2 0]")

        # 5. Non-strict mode reverse detection
        lattice.strict = False
        key, reverse = lattice._get_node_vector_pair("[1 1 1]", "[3 3 1]")

        assert not reverse, "Expected reverse=False in non-strict mode"
        np.testing.assert_equal(actual=key, desired=(1, 5))

        # 6. Coverage for reverse branch in _create_hydrogel_chain
        # --------------------------------------------------------


        # Register molecule template and default_linker particle
        pmb2.define_molecule(name="test_chain",
                            residue_list=sequence)
        pmb2.define_particle(name="default_linker",
                             sigma=1*pmb2.units.reduced_length,
                             epsilon=1*pmb2.units.reduced_energy)
        pmb2.define_default_bond(bond_type="harmonic",
                                 bond_parameters={"r_0": 1*pmb2.units.reduced_length,
                                                  "k": 400 * pmb2.units('reduced_energy / reduced_length**2')})
        
        # Define nodes dictionary as expected by _create_hydrogel_chain
        # using actual hydrogel node particles registered in pyMBE.
        nodes = {}
        for label, index in lattice.node_labels.items():
            node_index = lattice.lattice.indices[index]
            node_name = lattice.get_node(label)
            node_pos, node_id = pmb2._create_hydrogel_node(box_l=box_l,
                                                           node_index=node_index,
                                                           node_name=node_name)
            nodes[label] = {"name": node_name,
                            "pos": node_pos,
                            "id": node_id}
        from pyMBE.storage.templates.hydrogel import HydrogelChain
        # Define hydrogel chain template (reverse geometry)
        hydrogel_chain = HydrogelChain(node_start=node_b,   # reversed on purpose
                                      node_end=node_a,
                                      molecule_name="test_chain")
        mol_id = pmb2._create_hydrogel_chain(hydrogel_chain=hydrogel_chain,
                                            nodes=nodes,
                                            box_l=box_l,
                                            use_default_bond=True)
        # Extract created particle IDs
        chain_pids = pmb2.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                            attribute="molecule_id",
                                                            value=mol_id)
        # Extract residue sequence from particle instances
        created_residues_id = [pmb2.db.get_instance("particle", pid).residue_id  for pid in chain_pids]
        created_residues    = [pmb2.db.get_instance("residue",  rid).name  for rid in created_residues_id]
        # Reverse branch MUST reverse the residue list
        np.testing.assert_equal(actual=created_residues, desired=sequence[::-1], verbose=True)

    def test_non_palindromic_self_loop_chain_raises(self):
        """
        A self-loop hydrogel chain with a non-palindromic residue sequence is
        ambiguous and must be rejected by _create_hydrogel_chain.
        """
        espresso_system.part.clear()
        espresso_system.bonded_inter.clear()
        espresso_system.non_bonded_inter.reset()

        pmb = pyMBE.pymbe_library(seed=11)
        define_templates(pmb=pmb)
        lattice = pmb.initialize_lattice_builder(diamond)
        lattice.strict = False

        non_palindromic_sequence = [Res1, Res2, Res3, Res1]
        pmb.define_molecule(name="self_loop_chain",
                            residue_list=non_palindromic_sequence)

        from pyMBE.storage.templates.hydrogel import HydrogelChain
        hydrogel_chain = HydrogelChain(node_start="[0 0 0]",
                                       node_end="[0 0 0]",
                                       molecule_name="self_loop_chain")

        with self.assertRaisesRegex(
            ValueError,
            "could not resolve a unique topology for that chain",
        ):
            pmb._create_hydrogel_chain(hydrogel_chain=hydrogel_chain,
                                       nodes={},
                                       box_l=box_l,
                                       use_default_bond=False)

    def test_plot(self):
        pmb = pyMBE.pymbe_library(seed=42)
        diamond = pyMBE.lib.lattice.DiamondLattice(mpc, bond_l)
        lattice = pmb.initialize_lattice_builder(diamond)
        define_templates(pmb)
        pmb.define_molecule(name="test",
                            residue_list=[Res1])
        # Setting up chain topology
        connectivity = diamond.connectivity
        node_labels = lattice.node_labels
        reverse_node_labels = {v: k for k, v in node_labels.items()}
        connectivity_with_labels = {(reverse_node_labels[i], reverse_node_labels[j]) for i, j in connectivity}
        chain_topology = []

        for node_s, node_e in connectivity_with_labels:
            chain_topology.append({'node_start':node_s,
                                'node_end': node_e,
                                'molecule_name':"test"})
        # --- Colormap ---
        lattice.set_colormap(self.colormap)
        for index, (label, color) in enumerate(self.colormap.items()):
            np.testing.assert_equal(lattice.get_monomer_color(label), color)
            np.testing.assert_equal(lattice.get_monomer_color_index(label), index)

        # --- Invalid colormap access ---
        with self.assertRaisesRegex(
            RuntimeError, "monomer 'unknown' has no associated color"
        ):
            lattice.get_monomer_color("unknown")

        with self.assertRaises(AssertionError):
            lattice.set_colormap("red")

        # --- Invalid node access ---
        with self.assertRaisesRegex(
            AssertionError, r"node '\[0 5 13\]' doesn't exist in a diamond lattice"
        ):
            lattice.get_node("[0 5 13]")

        # --- Plot smoke tests ---
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(projection="3d", computed_zorder=False)
        lattice.chains= chain_topology
        lattice.draw_lattice(ax,
                             pmb=pmb)
        lattice.draw_simulation_box(ax)
        plt.close(fig)

   

if __name__ == "__main__":
    ut.main()
