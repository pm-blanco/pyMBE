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

from pyMBE.storage.base_type import PMBBaseModel
from pydantic import validator
from typing import Optional

class NanoparticleInstance(PMBBaseModel):
    """
    Persistent instance representation of a nanoparticle.

    Attributes:
        pmb_type ('str'):
            Fixed string identifying this object as a nanoparticle instance.
            Always ``"nanoparticle"``.
        
        name ('str'):
            Name of the nanoparticle **template** from which this instance was
            created. This must correspond to an existing
            ``NanoparticleTemplate`` in the database.

        molecule_id ('int'):
            Unique non-negative integer identifying the nanoparticle instance within the database.

        assembly_id (int | None):
            Identifier of the super-parent assembly (e.g. hydrogel) to which this nanoparticle belongs. ``None`` indicates that the nanoparticle is not assigned to any assembly.
    """
    pmb_type: str = "nanoparticle"
    name: str
    molecule_id: int 
    assembly_id: Optional[int] = None
    
    @validator("molecule_id")
    def validate_molecule_id(cls, mid):
        if mid < 0:
            raise ValueError("molecule_id must be a non-negative integer.")
        return mid
