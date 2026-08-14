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

from typing import Literal, Optional,List
from numpy.typing import NDArray
from pydantic import validator
from ..base_type import PMBBaseModel

class ParticleInstance(PMBBaseModel):
    """
    Concrete instance of a particle placed in the simulation.

    Attributes:
        pmb_type ('Literal["particle"]'):
            Fixed string identifying this object as a particle instance. Always ``"particle"``.

        name ('str'):
            Name of the particle template from which this instance is derived.

        particle_id ('int'):
            Unique non-negative integer identifying the particle within the database. Assigned sequentially by the database manager.

        initial_state ('str'):
            Name of the particle state at creation time. 

        residue_id ('Optional[int]'):
            Optional identifier of the ``ResidueInstance`` this particle belongs to. Particles that are not part of a residue should  leave this field as ``None``.

        molecule_id ('Optional[int]'):
            Optional identifier of the ``MoleculeInstance`` this particle  belongs to. Particles not belonging to any molecule should keep this as ``None``.

        assembly_id ('Optional[int]'):
            Identifier of the super-parent assembly (e.g. hydrogel) to which this particle instance belongs. ``None`` indicates that the particle is not assigned to any assembly.

    Notes:
        - ``initial_state`` is stored as a plain string to ensure clean serialization and avoid engine-specific objects.
        - Connectivity, bonding, and spatial ordering are external to this class and handled by the database or simulation backend.
    """
    pmb_type: Literal["particle"] = "particle"
    name: str 
    particle_id: int
    initial_state: str
    position: NDArray
    fix: List[bool]=[False,False,False]
    residue_id: Optional[int] = None
    molecule_id: Optional[int] = None
    assembly_id: Optional[int] = None
    added_to_engine: bool = False
    
    class Config:
        arbitrary_types_allowed=True
        validate_assignment = True
        extra = "forbid"
        
    @validator("particle_id")
    def validate_particle_id(cls, pid):
        if pid < 0:
            raise ValueError("particle_id must be a non-negative integer.")
        return pid
