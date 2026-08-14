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

from pyMBE.simulation_builder.base_engine import SimulationEngine


class LammpsSimulation(SimulationEngine):
    def __init__(self, box_l=None, db=None, lammps=None, units=None,
                 kT=None, Kw=None, seed=None):
        self.box_l = box_l
        self.db = db
        self.lammps = lammps
        self.units = units
        self.kT = kT
        self.Kw = Kw
        self.seed = seed

    def display_commands(self, *args, **kwargs):
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')

    def setup_lj_interactions(self, *args, **kwargs):
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')

    def setup_langevin(self, *args, **kwargs):
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')

    def setup_cpH(self, *args, **kwargs):
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')

    def run_simulation(self, *args, **kwargs):
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')

    def __getattr__(self, attr):
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')
    def _add_angle(self): # pragma: no cover
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')
    def _check_bond_inputs(self): # pragma: no cover
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')
    
    def _create_bond_instance(self): # pragma: no cover
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')
    def _get_bond_instance(self): # pragma: no cover
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')
    def add_instances_to_engine(self): # pragma: no cover
        raise NotImplementedError('Lammps Simulation Engine is not yet implemented')

