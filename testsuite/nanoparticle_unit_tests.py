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

import unittest as ut

import espressomd
import pyMBE


class TestNanoparticleCreation(ut.TestCase):

    def test_create_nanoparticle_registers_instances(self):
        pmb = pyMBE.pymbe_library(seed=42)
        pmb.set_reduced_units(unit_length=0.4 * pmb.units.nm,
                              Kw=1e-14)

        pmb.define_particle(name="core",
                            z=0,
                            sigma=4.0 * pmb.units("reduced_length"),
                            epsilon=1.0 * pmb.units("reduced_energy"))
        pmb.define_particle(name="A",
                            z=-1,
                            sigma=1.0 * pmb.units("reduced_length"),
                            epsilon=1.0 * pmb.units("reduced_energy"))
        pmb.define_particle(name="B",
                            z=1,
                            sigma=1.0 * pmb.units("reduced_length"),
                            epsilon=1.0 * pmb.units("reduced_energy"))

        nanoparticle_name = "np"
        pmb.define_nanoparticle(name=nanoparticle_name,
                                core_particle_name="core",
                                surface_density_of_sites=0.2 * pmb.units("reduced_length^-2"),
                                primary_site_particle_name="A",
                                fraction_primary_sites=0.5,
                                number_of_patches_of_primary_sites=2,
                                secondary_site_particle_name="B")

        espresso_system = espressomd.System(box_l=[40, 40, 40])
        core_positions = [[5.0, 5.0, 5.0],
                          [25.0, 25.0, 25.0]]
        created = pmb.create_nanoparticle(name=nanoparticle_name,
                                          number_of_nanoparticles=2,
                                          espresso_system=espresso_system,
                                          list_core_particle_positions=core_positions,
                                          fix=True)

        self.assertEqual(sorted(created.keys()), [0, 1])

        np_instances = pmb.db.get_instances(pmb_type="nanoparticle")
        self.assertEqual(sorted(np_instances.keys()), [0, 1])
        self.assertEqual(np_instances[0].molecule_id, 0)
        self.assertEqual(np_instances[1].molecule_id, 1)

        nanoparticle_tpl = pmb.db.get_template(pmb_type="nanoparticle",
                                               name=nanoparticle_name)
        properties = nanoparticle_tpl.calculate_nanoparticle_properties(pmb)
        expected_particles_per_np = 1 + properties["number_of_primary_sites"] + properties["number_of_secondary_sites"]

        for nanoparticle_id, nanoparticle_info in created.items():
            core_particle_id = nanoparticle_info["core_particle_id"]
            core_instance = pmb.db.get_instance(pmb_type="particle",
                                                instance_id=core_particle_id)
            self.assertEqual(core_instance.name, "core")
            self.assertEqual(core_instance.molecule_id, nanoparticle_id)

            core_pos = list(espresso_system.part.by_id(core_particle_id).pos)
            self.assertListEqual(core_pos, core_positions[nanoparticle_id])
            self.assertListEqual(list(espresso_system.part.by_id(core_particle_id).fix), [True, True, True])

            all_sites_ids = nanoparticle_info["all_sites_ids"]
            for site_id in all_sites_ids:
                site_instance = pmb.db.get_instance(pmb_type="particle",
                                                    instance_id=site_id)
                self.assertEqual(site_instance.molecule_id, nanoparticle_id)
                self.assertListEqual(list(espresso_system.part.by_id(site_id).fix), [True, True, True])

            self.assertEqual(1 + len(all_sites_ids), expected_particles_per_np)

        particle_id_map = pmb.get_particle_id_map(object_name=nanoparticle_name)
        self.assertEqual(len(particle_id_map["all"]), expected_particles_per_np * 2)


if __name__ == "__main__":
    ut.main()
