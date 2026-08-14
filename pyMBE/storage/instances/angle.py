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

class AngleInstance(PMBBaseModel):
    """
    Instance representation of an angle between three particles.

    Attributes:
        pmb_type ('Literal["angle"]'):
            Fixed identifier set to ``"angle"`` for all angle instances.

        angle_id ('int'):
            Unique non-negative integer identifying this angle instance.

        name ('str'):
            Name of the angle template from which this instance was created.

        particle_id1 ('int'):
            ID of the first side particle.

        particle_id2 ('int'):
            ID of the central particle.

        particle_id3 ('int'):
            ID of the second side particle.
    """
    pmb_type: Literal["angle"] = "angle"
    angle_id: int
    name: str
    particle_id1: int
    particle_id2: int
    particle_id3: int
    added_to_engine: bool = False

    @validator("angle_id", "particle_id1", "particle_id2", "particle_id3")
    def validate_non_negative_int(cls, value, field):
        if value < 0:
            raise ValueError(f"{field.name} must be a non-negative integer.")
        return value
