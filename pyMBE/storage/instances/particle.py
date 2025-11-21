# instances/particle.py

from pydantic import Field, field_validator
from ..base_type import PMBBaseModel


class ParticleInstance(PMBBaseModel):
    """
    Instantiated particle in the system:
    - particle_id is unique
    - template_name links to ParticleTemplate
    - can override active_state if needed
    """

    pmb_type: str = Field(default="particle", frozen=True)

    particle_id: int
    name: str
    active_state: str | None = None

    @field_validator("particle_id")
    def validate_particle_id(cls, pid):
        if pid < 0:
            raise ValueError("particle_id must be a non-negative integer.")
        return pid
