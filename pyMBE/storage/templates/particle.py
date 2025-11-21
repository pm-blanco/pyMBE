from typing import Dict
from pydantic import Field, field_validator

from ..base_type import PMBBaseModel


class ParticleState(PMBBaseModel):
    pmb_type: str = Field(default="particle", frozen=True)

    label: str
    es_type: int
    charge: int


class ParticleTemplate(PMBBaseModel):
    """
    Template describing the type of particle:
    - sigma, epsilon
    - allowed states
    - template_name = unique string identifier
    """

    pmb_type: str = Field(default="particle", frozen=True)

    name: str
    sigma: float
    epsilon: float

    states: Dict[str, ParticleState] = Field(default_factory=dict)
    default_state: str | None = None

    # ---------------- Validators -----------------

    @field_validator("default_state")
    def validate_default_state(cls, v, values):
        if v is None:
            return v
        if "states" in values and v not in values["states"]:
            raise ValueError(
                f"default_state '{v}' not found in states "
                f"({list(values['states'].keys())})"
            )
        return v

    # ---------------- Helpers -----------------

    def add_state(self, state: ParticleState):
        if state.label in self.states:
            raise ValueError(f"State '{state.label}' already exists.")
        self.states[state.label] = state
