from typing import Optional
from pydantic import Field, field_validator
from ..base_type import PMBBaseModel


class ParticleInstance(PMBBaseModel):
    """
    A placed particle within the simulation.
    """
    pmb_type: str = "particle"
    particle_id: int
    state_name: str

    @field_validator("particle_id")
    def validate_particle_id(cls, pid):
        if pid < 0:
            raise ValueError("particle_id must be a non-negative integer.")
        return pid
