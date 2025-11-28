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

import itertools
import numpy as np

class LatticeBuilder:
    """
    Generic lattice builder.
    Args:
        lattice(`object`): Data class that represents a lattice, for example a :class:`DiamondLattice`.
        strict(`bool`): If set to `True`, the lattice connectivity cannot be altered.
            For example, a chain cannot be created between two nodes that are
            not explicitly connected according to the definition in `lattice`.
    Attributes:
        kwargs_node_labels(`dict`): Keyword arguments passed to the matplotlib
            `Axes3D.text() <https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.text.html>`_
            function to draw lattice node labels.
        kwargs_monomers(`dict`): Keyword arguments passed to the matplotlib
            `Axes3D.scatter() <https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.scatter.html>`_
            function to draw chain monomers.
        kwargs_bonds(`dict`): Keyword arguments passed to the matplotlib
            `Axes3D.plot() <https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.plot.html>`_
            function to draw chain bonds.
        kwargs_box(`dict`): Keyword arguments passed to the matplotlib
            `Axes3D.plot() <https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.mplot3d.axes3d.Axes3D.plot.html>`_
            function to draw the simulation box.
    """

    def __init__(self, lattice, strict=True):
        self.lattice = lattice
        self.strict = strict
        self.node_labels = {str([int(x) for x in indices]).replace(",", ""): i
                            for i, indices in enumerate(self.lattice.indices)}
        self.chain_labels = {f"({start}, {end})": idx for idx, (start, end) in enumerate(lattice.connectivity)}
        self.nodes = {label: "default_linker" for label in self.node_labels}
        self.colormap = {}
        self.chains = {}
        self.kwargs_node_labels = {}
        self.kwargs_monomers = {}
        self.kwargs_bonds = {}
        self.kwargs_box = {}
        self.mpc = lattice.mpc
        self.box_l = lattice.box_l
        

    def _get_node_by_label(self, node):
        assert node in self.node_labels, f"node '{node}' doesn't exist in a {self.lattice.name} lattice"
        return node

    def _get_node_vector_pair(self, node_start, node_end):
        node1 = self._get_node_by_label(node_start)
        node2 = self._get_node_by_label(node_end)
        reverse = False
        key = None
        key_a = (self.node_labels[node1], self.node_labels[node2])
        key_b = key_a[::-1]
        if key_a in self.lattice.connectivity:
            key = key_a
        elif key_b in self.lattice.connectivity:
            key = key_b
        else:
            assert not self.strict, f"there is no chain between '{node_start}' and '{node_end}' in a {self.lattice.name} lattice (strict mode is enabled)"
            key = key_b if key_b in self.chains else key_a
        if key == key_b:
            reverse = True
        return key, reverse

    @classmethod
    def _make_sphere(cls, radius, resolution):
        u = np.linspace(0, 2 * np.pi, resolution[0])
        v = np.linspace(0, np.pi, resolution[1])
        x = radius * np.outer(np.cos(u), np.sin(v))
        y = radius * np.outer(np.sin(u), np.sin(v))
        z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
        return (x, y, z)

    def add_default_chains(self, mpc):
        """
        Build a default lattice network. Skip node pairs that already have a chain defined.
        Args:
            mpc(`int`): Number of monomers per chain.
        """
        for key in self.lattice.connectivity:
            if key not in self.chains:
                self.chains[key] = mpc * ["default_monomer"]

    def draw_lattice(self, ax):
        """
        Draw the lattice in an `Axes3D <https://matplotlib.org/stable/api/toolkits/mplot3d/axes3d.html>`_ canvas.
        Args:
            ax: Axes.
        """
        kwargs_node_labels = {"zdir": (1., 1., 1.), "horizontalalignment": "left", "verticalalignment": "bottom", **self.kwargs_node_labels}
        kwargs_bonds = {"linestyle": "-", "marker": None, "color": "gray", **self.kwargs_bonds}
        kwargs_monomers = {**self.kwargs_monomers}
        scatter_data = {}
        # gather monomers at lattice nodes
        for node_label, node_type in self.nodes.items():
            node_id = self.node_labels[node_label]
            for image_box in itertools.product((0, 4), repeat=3):
                image_indices = self.lattice.indices[node_id] + np.array(image_box)
                if np.max(image_indices) <= 4:
                    image_label = str([int(x) for x in image_indices]).replace(",", "")
                    ax.text(*image_indices + np.array([-0.15, 0., 0.]), image_label, **kwargs_node_labels)
                    if node_type not in scatter_data:
                        scatter_data[node_type] = []
                    scatter_data[node_type].append(image_indices)
        # gather monomers from the chains
        for (start_node, end_node), sequence in self.chains.items():
            node_connection_vec = (self.lattice.indices[end_node, :] - self.lattice.indices[start_node, :]) / 4.
            node_connection_vec -= np.rint(node_connection_vec)
            node_connection_vec *= 4.
            bond_vector = node_connection_vec / (len(sequence) + 1)
            for j in range(len(sequence) + 1):
                pos = self.lattice.indices[start_node, :]
                vec = np.vstack((pos + (j + 0) * bond_vector,
                                 pos + (j + 1) * bond_vector))
                ax.plot(vec[:, 0], vec[:, 1], zs=vec[:, 2], zorder=1, **kwargs_bonds)
            # draw bonds
            for j, node_type in enumerate(sequence):
                pos = self.lattice.indices[start_node, :] + (j + 1) * bond_vector
                if node_type not in scatter_data:
                    scatter_data[node_type] = []
                scatter_data[node_type].append(pos)
        # draw monomers
        resolution = (16, 8)
        self.sphere = self._make_sphere(radius=0.1, resolution=resolution)
        node_types = scatter_data.keys()
        if self.colormap:
            node_types = sorted(node_types, key=lambda x: self.get_monomer_color(x))
        for node_type in node_types:
            if self.colormap:
                color = self.colormap[node_type]
                kwargs_monomers["c"] = color
                
            node_positions = scatter_data[node_type]
            pos = np.array(node_positions)
            # plotting nodes and monomers
            ax_data = ax.scatter(pos[:,0], pos[:,1], pos[:,2], edgecolor="none",
                                 zorder=2, label=node_type, s=12**2, **kwargs_monomers)
            color = ax_data.get_facecolors()[0]
            facecolors = np.tile(color, resolution).reshape((*resolution, len(color)))
            for x, y, z in node_positions:
                ax.plot_surface(x + self.sphere[0], y + self.sphere[1], z + self.sphere[2], zorder=3,
                                shade=False, facecolors=facecolors)

    def draw_simulation_box(self, ax):
        """
        Draw the simulation box in an `Axes3D <https://matplotlib.org/stable/api/toolkits/mplot3d/axes3d.html>`_ canvas.
        Args:
            ax: Axes.
        """
        kwargs_box = {"linestyle": "-", "marker": None, "color": "black", **self.kwargs_box}
        for x in (0, 1):
            ax.plot(np.array([x, x, x, x, x]) * 4.,
                    np.array([0, 0, 1, 1, 0]) * 4.,
                    zs=np.array([0, 1, 1, 0, 0]) * 4.,
                    **kwargs_box)
        for y, z in itertools.product((0, 1), repeat=2):
            ax.plot(np.array([0, 1]) * 4.,
                    np.array([y, y]) * 4.,
                    zs=np.array([z, z]) * 4.,
                    **kwargs_box)

    def get_chain(self, node_start, node_end):
        """
        Get a chain between two nodes.
        Args:
            node_start(`str`): Start node label, e.g. ``[0 0 0]`` for the node at the origin.
            node_end(`str`): End node label, e.g. ``[1 1 1]`` for the first node along the diagonal.
        Returns:
            `list`: Sequence of residue labels.
        """
        key, reverse = self._get_node_vector_pair(node_start, node_end)
        if key not in self.chains:
            raise RuntimeError(f"no chain has been defined between '{node_start}' and '{node_end}' yet")
        sequence = self.chains[key]
        if reverse:
            sequence = sequence[::-1]
        return sequence

    def get_monomer_color(self, label):
        """
        Get the color corresponding to a monomer label.
        Only works if a custom color map was set via :meth:`set_colormap`!
        Args:
            label(`str`): Monomer label.
        Returns:
            The color.
        """
        self.get_monomer_color_index(label) # to run assertions
        return self.colormap[label]

    def get_monomer_color_index(self, label):
        """
        Get the color wheel index corresponding to a monomer label.
        Only works if a custom color map was set via :meth:`set_colormap`!
        Args:
            label(`str`): Monomer label.
        Returns:
            `int`: The color index.
        """
        if label not in self.colormap:
            raise RuntimeError(f"monomer '{label}' has no associated color in the colormap")
        return self.colormap_sorted_names.index(label)

    def get_node(self, node):
        """
        Get a node residue label.
        Args:
            node(`str`): Node label, e.g. ``[0 0 0]`` for the node at the origin.
        Returns:
            `str`: Residue label.
        """
        key = self._get_node_by_label(node)
        return self.nodes[key]

    def set_colormap(self, colormap):
        """
        Set a discrete color map. By default, the standard matplotlib color wheel will be used.
        Args:
            colormap(`dict`): Discrete color wheel that maps monomer labels to colors.
        """
        assert isinstance(colormap, dict)
        self.colormap = colormap.copy()
        self.colormap_sorted_names = list(self.colormap.keys())
    
    def set_chain(self, node_start, node_end, sequence):
        """
        Set a chain between two nodes.
        Args:
            node_start(`str`): Start node label, e.g. ``[0 0 0]`` for the node at the origin.
            node_end(`str`): End node label, e.g. ``[1 1 1]`` for the first node along the diagonal.
            sequence(`list`): Sequence of residue labels.
        """
        assert len(sequence) != 0 and not isinstance(sequence, str)
        key, reverse = self._get_node_vector_pair(node_start, node_end)
        assert node_start != node_end or sequence == sequence[::-1], \
            (f"chain cannot be defined between '{node_start}' and '{node_end}' since it "
            "would form a loop with a non-symmetric sequence (under-defined stereocenter)")
        if reverse:
            sequence = sequence[::-1]
        self.chains[key] = sequence
    
    def set_node(self, node, residue):
        """
        Set a node residue type.
        Args:
            node(`str`): Node label, e.g. ``[0 0 0]`` for the node at the origin.
            residue(`str`): Node residue label.
        """
        key = self._get_node_by_label(node)
        self.nodes[key] = residue

class DiamondLattice:
    """
    Representation of the diamond lattice.
    """
    name = "diamond"
    indices = np.array([[0, 0, 0], [1, 1, 1],
                        [2, 2, 0], [0, 2, 2],
                        [2, 0, 2], [3, 3, 1],
                        [1, 3, 3], [3, 1, 3]])
    connectivity = {(0, 1), (1, 2), (1, 3), (1, 4),
                    (2, 5), (3, 6), (4, 7), (5, 0),
                    (5, 3), (5, 4), (6, 0), (6, 2),
                    (6, 4), (7, 0), (7, 2), (7, 3)}
    def __init__(self,mpc,bond_l):
        if not isinstance(mpc, int) or mpc <= 0:
            raise ValueError("mpc must be a non-zero positive integer.")
        self.mpc = mpc
        self.bond_l = bond_l
        self.box_l = (self.mpc+1)*self.bond_l.magnitude / (np.sqrt(3)*0.25)

