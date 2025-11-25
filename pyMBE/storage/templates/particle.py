from typing import Dict, Literal
from pydantic import Field, field_validator

from ..base_type import PMBBaseModel


class ParticleState(PMBBaseModel):
    pmb_type: Literal["particle_state"] = "particle_state"
    name: str                      # e.g. "HA", "A-", "H+"
    charge: int
    es_type: float                  # label in espresso


class ParticleTemplate(PMBBaseModel):
    """
    Template describing the type of particle:
    - sigma, epsilon
    - allowed states
    - template_name = unique string identifier
    """

    pmb_type: str = Field(default="particle", frozen=True)
    sigma: float
    epsilon: float
    states: Dict[str, ParticleState] = {}

    # ---------------- Validators -----------------

    def add_state(self, state: ParticleState):
        if state.name in self.states:
            raise ValueError(f"State {state.name} already exists in template {self.name}")
        self.states[state.name] = state

    @classmethod
    def single_state(cls, name: str, charge: int, es_type: str, epsilon: float = 1.0):
        """
        Convenience constructor for particles such as H+ that only need one state.
        """
        state = ParticleState(name=name, charge=charge, es_type=es_type)
        return cls(name=name, epsilon=epsilon, states={name: state})

