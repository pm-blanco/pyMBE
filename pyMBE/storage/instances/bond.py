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

from typing import Literal
from pyMBE.storage.base_type import PMBBaseModel
from pydantic import validator

class BondInstance(PMBBaseModel):
    """
    Instance representation of a bond between two particles.

    Attributes:
        pmb_type ('Literal["bond"]'):
            Fixed identifier set to ``"bond"`` for all bond instances.

        bond_id ('int'):
            Unique non-negative integer identifying this bond instance.

        name ('str'):
            Name of the bond template from which this instance was created.

        particle_id1 ('int'):
            ID of the first particle involved in the bond.

        particle_id2 ('int'):
            ID of the second particle involved in the bond.

    Notes:
        - ``particle_id1`` and ``particle_id2`` must correspond to
          particle instance IDs already registered in the database.
        - This class does **not** store simulation engine–specific
          objects (e.g., Espresso bond handles). Those should be created
          by a runtime builder separate from the persistent database.
    """
    pmb_type: Literal["bond"] = "bond"
    bond_id: int
    name : str            # bond template name
    particle_id1: int
    particle_id2: int
    added_to_engine: bool = False

    @validator("bond_id", "particle_id1", "particle_id2")
    def validate_non_negative_int(cls, value, field):
        if value < 0:
            raise ValueError(f"{field.name} must be a non-negative integer.")
        return value
    
