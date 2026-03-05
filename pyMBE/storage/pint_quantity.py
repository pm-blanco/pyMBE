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
 
from dataclasses import dataclass
import pint
# dimension -> representative unit used to check dimensionality
_DIMENSION_REPRESENTATIVE = {"length": "nm",
                             "energy": "meV",
                             "energy/length**2": "meV/nm**2",
                             "length**-2": "nm**-2",
                             "dimensionless": "dimensionless",} # extend as needed

@dataclass
class PintQuantity:
    """
    Internal representation of a Pint quantity for pyMBE storage.

    Attributes:
        magnitude ('float'): 
            Numeric value of the quantity in the stored units.

        units ('str'): 
            String representation of the units (e.g., "nm", "meV", "meV/nm**2").

        dimension ('str'): 
            Logical dimension of the quantity, e.g., "length", "energy", etc.
    
    Notes:
        - Stores the magnitude and units of a quantity in a base/SI-like format along with its logical physical dimension.
    """

    magnitude: float
    units: str        
    dimension: str    

    @classmethod
    def from_quantity(cls, q, expected_dimension, ureg):
        """
        Create a PintQuantity from a Pint Quantity, validating its dimension.

        Args:
            q ('pint.Quantity'): 
                Pint Quantity to store.

            expected_dimension ('str'): 
                Expected logical dimension ("length", "energy", etc.).

            ureg ('pint.UnitRegistry'): 
                Pint UnitRegistry used for unit conversion.

        Returns:
            'PintQuantity': 
                Internal pyMBE representation in SI units.
        """
        if not isinstance(q, pint.Quantity):
            raise TypeError("from_quantity expects a pint.Quantity")

        # Build a representative unit for the dimension using the provided registry
        if expected_dimension not in _DIMENSION_REPRESENTATIVE:
            raise ValueError(f"Unknown expected_dimension '{expected_dimension}'")

        rep_unit = ureg(_DIMENSION_REPRESENTATIVE[expected_dimension])

        # Use pint's dimensionality check
        try:
            if not q.check(rep_unit):
                raise ValueError(f"Quantity {q} does not have expected dimension '{expected_dimension}'")
        except Exception as e:
            # If check fails because registries differ
            raise

        # Use the dimension representative unit
        rep_unit_str = _DIMENSION_REPRESENTATIVE[expected_dimension]
        rep_unit = ureg(rep_unit_str)

        # Convert to the representative SI unit
        q_base = q.to(rep_unit)

        # Store magnitude and unit name
        mag = float(q_base.magnitude)
        unit_str = rep_unit_str
        return cls(magnitude=mag, units=unit_str, dimension=expected_dimension)

    def to_quantity(self, ureg):
        """
        Convert the stored PintQuantity back into a Pint Quantity.

        Args:
            ureg ('pint.UnitRegistry'): 
                Pint UnitRegistry used to construct the Quantity.

        Returns:
            'pint.Quantity': 
                Pint Quantity with the stored magnitude and units.
        """
        return self.magnitude * ureg(self.units)

    def to_dict(self):
        """
        Serialize the PintQuantity to a dictionary.

        Returns:
            'dict': 
                Dictionary with keys "magnitude", "units", and "dimension".
        """
        return {"magnitude": self.magnitude, "units": self.units, "dimension": self.dimension}

    @classmethod
    def from_dict(cls, d):
        """
        Deserialize a PintQuantity from a dictionary.

        Args:
            d ('dict'): 
                Dictionary containing "magnitude", "units", and "dimension".

        Returns:
            'pint.PintQuantity': 
                Reconstructed PintQuantity object.
        """
        return cls(magnitude=d["magnitude"], units=d["units"], dimension=d["dimension"])

