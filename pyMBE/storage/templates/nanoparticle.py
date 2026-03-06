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
from pydantic import Field
from pyMBE.storage.pint_quantity import PintQuantity
import numpy as np

class NanoparticleTemplate(PMBBaseModel):
    """
    Template defining a nanoparticle in the pyMBE database.

    Attributes:
        pmb_type ('str'): 
            Fixed type identifier. Always "nanoparticle".

        name ('str'): 
            Unique name of the nanoparticle template.

        core_particle_name ('str'):
            Name of the particle template used as the nanoparticle core.

        surface_density_of_sites ('PintQuantity'):
            Surface density of grafting/interaction sites on the nanoparticle core.

        primary_site_particle_name ('str'):
            Name of the particle template used for the primary site type.

        fraction_primary_sites ('float'):
            Fraction of all surface sites assigned to the primary site type.
            Expected range is typically between 0 and 1.

        number_of_patches_of_primary_sites ('int'):
            Number of surface patches that contain the primary site type.

        secondary_site_particle_name ('str | None'):
            Optional particle template name for a secondary site type.
            If not provided, only a single site type is used.

    """
    pmb_type: str = Field(default="nanoparticle", frozen=True)
    name: str
    core_particle_name: str
    surface_density_of_sites: PintQuantity
    primary_site_particle_name: str
    fraction_primary_sites: float
    number_of_patches_of_primary_sites: int 
    secondary_site_particle_name: str | None = None

    def calculate_nanoparticle_properties(self, pmb):
        """
        Compute various nanoparticle properties from template parameters.

        The method uses:
        - Core radius from the core particle template.
        - Site counts inferred from surface area and site surface density.
        - Site charges taken from each site's particle initial state.

        Args:
            pmb ('pyMBE.pymbe_library'):
                Active pyMBE object with a populated template database and unit registry.

        Returns:
            ('dict'): Dictionary with geometric, compositional, and electrostatic properties:
                - ``nanoparticle_surface_area`` ('pint.Quantity')
                - ``nanoparticle_volume`` ('pint.Quantity')
                - ``total_number_of_sites`` ('int')
                - ``real_surface_density_of_sites`` ('pint.Quantity')
                - ``number_of_primary_sites`` ('int')
                - ``number_of_primary_sites_per_patch`` ('int')
                - ``number_of_secondary_sites`` ('int')
                - ``real_fraction_primary_sites`` ('float')
                - ``primary_site_charge_number`` ('int')
                - ``secondary_site_charge_number`` ('int')
                - ``total_charge`` ('pint.Quantity')
                - ``surface_charge_density`` ('pint.Quantity')
                - ``volume_charge_density`` ('pint.Quantity')

        Notes:
            - If ``secondary_site_particle_name`` is not set, all sites are assigned
              to the primary site type.
            - Primary-site counts are rounded to ensure an integer number of sites
              per patch and exact patch occupancy.
        """
        if not (0.0 <= self.fraction_primary_sites <= 1.0):
            raise ValueError("fraction_primary_sites must be between 0 and 1.")
        if self.number_of_patches_of_primary_sites <= 0:
            raise ValueError("number_of_patches_of_primary_sites must be > 0.")

        def _get_initial_state_charge_number(particle_name: str) -> int:
            particle_tpl = pmb.db.get_template(name=particle_name, pmb_type="particle")
            state_name = particle_tpl.initial_state
            if state_name is None:
                particle_states = pmb.db.get_particle_states_templates(particle_name=particle_name)
                if not particle_states:
                    raise ValueError(f"Particle '{particle_name}' has no defined states.")
                state_name = next(iter(particle_states.keys()))
            state_tpl = pmb.db.get_template(name=state_name, pmb_type="particle_state")
            return state_tpl.z

        core_particle_tpl = pmb.db.get_template(name=self.core_particle_name, pmb_type="particle")
        core_initial_state_name = core_particle_tpl.initial_state
        if core_initial_state_name is None:
            raise ValueError(
                f"Core particle '{self.core_particle_name}' has no initial_state set."
            )
        core_initial_state = pmb.db.get_template(name=core_initial_state_name, pmb_type="particle_state")
        core_radius = pmb.get_radius_map(dimensionless=False)[core_initial_state.es_type]

        nanoparticle_surface_area = 4.0 * np.pi * core_radius**2
        nanoparticle_volume = (4.0 / 3.0) * np.pi * core_radius**3

        surface_density_of_sites = self.surface_density_of_sites.to_quantity(pmb.units)
        
        total_number_of_sites = int(
            np.round((nanoparticle_surface_area * surface_density_of_sites).to("dimensionless").magnitude)
        )
        real_surface_density_of_sites = total_number_of_sites / nanoparticle_surface_area

        if self.secondary_site_particle_name is None:
            number_of_primary_sites = total_number_of_sites
            number_of_primary_sites_per_patch = int(
                np.round(number_of_primary_sites / self.number_of_patches_of_primary_sites)
            )
            real_number_of_primary_sites = (
                number_of_primary_sites_per_patch * self.number_of_patches_of_primary_sites
            )
            number_of_secondary_sites = 0
            secondary_site_charge_number = 0
        else:
            number_of_primary_sites = int(np.round(total_number_of_sites * self.fraction_primary_sites))
            number_of_primary_sites_per_patch = int(
                np.round(number_of_primary_sites / self.number_of_patches_of_primary_sites)
            )
            real_number_of_primary_sites = (
                number_of_primary_sites_per_patch * self.number_of_patches_of_primary_sites
            )
            number_of_secondary_sites = total_number_of_sites - real_number_of_primary_sites
            secondary_site_charge_number = _get_initial_state_charge_number(
                self.secondary_site_particle_name
            )

        if total_number_of_sites > 0:
            real_fraction_primary_sites = real_number_of_primary_sites / total_number_of_sites
        else:
            real_fraction_primary_sites = 0.0

        primary_site_charge_number = _get_initial_state_charge_number(self.primary_site_particle_name)
        total_charge_number = (
            real_number_of_primary_sites * primary_site_charge_number
            + number_of_secondary_sites * secondary_site_charge_number
        )
        total_charge = total_charge_number * pmb.units("reduced_charge")
        surface_charge_density = total_charge / nanoparticle_surface_area
        volume_charge_density = total_charge / nanoparticle_volume

        return {"nanoparticle_surface_area": nanoparticle_surface_area,
                "nanoparticle_volume": nanoparticle_volume,
                "total_number_of_sites": total_number_of_sites,
                "real_surface_density_of_sites": real_surface_density_of_sites,
                "number_of_primary_sites": real_number_of_primary_sites,
                "number_of_primary_sites_per_patch": number_of_primary_sites_per_patch,
                "number_of_secondary_sites": number_of_secondary_sites,
                "real_fraction_primary_sites": real_fraction_primary_sites,
                "primary_site_charge_number": primary_site_charge_number,
                "secondary_site_charge_number": secondary_site_charge_number,
                "total_charge": total_charge,
                "surface_charge_density": surface_charge_density,
                "volume_charge_density": volume_charge_density,}
