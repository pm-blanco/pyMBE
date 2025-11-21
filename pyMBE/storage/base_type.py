from pydantic import BaseModel, Field

class PMBBaseModel(BaseModel):
    """
    Base class for all pyMBE models:
    - Hard-coded pmb_type in subclasses
    """

    pmb_type: str = Field(frozen=True)

    class Config:
        validate_assignment = True
        extra = "forbid"
