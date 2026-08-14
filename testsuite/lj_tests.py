#
# Copyright (C) 2024-2026 pyMBE-dev team
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

# Import pyMBE and other libraries
import pyMBE
import numpy as np
import unittest as ut

# Create an instance of pyMBE library
pmb = pyMBE.pymbe_library(seed=42)
import espressomd
box_l=[50]*3
espresso_system=espressomd.System(box_l = box_l)

pmb.set_simulation_engine(espresso_system,box_l)


class Test(ut.TestCase):
    def test_particle_definition(self):
        """
        Unit test to check that define_particle stores correctly all LJ input parameters in the pyMBE database.
        """    
        input_parameters={"name":"D", 
                            "sigma":1*pmb.units.nm, 
                            "epsilon":pmb.units.Quantity(1,"reduced_energy"), 
                            "cutoff":2*pmb.units.nm, 
                            "offset":3*pmb.units.nm}

        pmb.define_particle(**input_parameters)
        part_tpl = pmb.db.get_template(name="D", 
                                    pmb_type="particle")
        for parameter_key in input_parameters.keys():
            atr = getattr(part_tpl, parameter_key)
            if isinstance(atr, str):
                self.assertEqual(first=atr,
                                 second=input_parameters[parameter_key])
            else:
                if parameter_key == "epsilon":    
                    self.assertAlmostEqual(first=atr.to_quantity(pmb.units).to("reduced_energy").magnitude,
                                    second=input_parameters[parameter_key].to("reduced_energy").magnitude)
                else:
                    self.assertEqual(first=atr.to_quantity(pmb.units).to("reduced_length").magnitude,
                                    second=input_parameters[parameter_key].to("reduced_length").magnitude)
            # Clean template from the database
        pmb.db.delete_template(name="D", 
                               pmb_type="particle")
        pmb.db.delete_template(name="D", 
                               pmb_type="particle_state")

        input_parameters={"name":"D", 
                            "sigma":1*pmb.units.nm, 
                            "epsilon":pmb.units.Quantity(1,"reduced_energy")}

        pmb.define_particle(**input_parameters)
        part_tpl = pmb.db.get_template(name="D", 
                                    pmb_type="particle")
        self.assertEqual(first=part_tpl.offset.to_quantity(pmb.units),
                            second=pmb.units.Quantity(0,"reduced_length"))
        self.assertEqual(first=part_tpl.cutoff.to_quantity(pmb.units), 
                            second=pmb.units.Quantity(2**(1./6.),"reduced_length"))
        # Clean template from the database
        pmb.db.delete_template(name="D", 
                            pmb_type="particle")
        # check that define_particle raises a ValueError if sigma is provided with the wrong dimensionality
        input_parameters={"name":"E", 
                        "sigma":1*pmb.units.ns, 
                            "epsilon":pmb.units.Quantity(1,"reduced_energy") }
        self.assertRaises(ValueError, pmb.define_particle, **input_parameters)
        # Unit test: check that define_particle raises a ValueError if offset is provided with the wrong dimensionality
        input_parameters={"name":"E", 
                        "offset":1*pmb.units.ns, 
                            "sigma":1*pmb.units.nm, 
                            "epsilon":pmb.units.Quantity(1,"reduced_energy") }
        self.assertRaises(ValueError, pmb.define_particle, **input_parameters)
        # Unit test: check that define_particle raises a ValueError if cutoff is provided with the wrong dimensionality
        input_parameters={"name":"E", 
                        "cutoff":1*pmb.units.ns, 
                            "sigma":1*pmb.units.nm, 
                            "epsilon":pmb.units.Quantity(1,"reduced_energy") }
        self.assertRaises(ValueError, pmb.define_particle, **input_parameters)
        # Unit test: check that define_particle raises a ValueError if epsilon is provided with the wrong dimensionality
        input_parameters={"name":"E", 
                        "epsilon":1*pmb.units.ns, 
                            "sigma":1*pmb.units.nm, }
        self.assertRaises(ValueError, pmb.define_particle, **input_parameters)
        
    def test_lj_interaction_setup(self):
        """
        Unit test to check that setup_lj_interactions sets up correctly LJ interactions between acid/base particles.
        """
        # Define particles
        A_input_parameters={"name":"A", 
                            "sigma":1*pmb.units.nm, 
                            "epsilon":pmb.units.Quantity(1,"reduced_energy"), 
                            "cutoff":2**(1./6.)*pmb.units.nm, 
                            "offset":1*pmb.units.nm}

        B_input_parameters={"name":"B", 
                            "sigma":2*pmb.units.nm, 
                            "epsilon":pmb.units.Quantity(2,"reduced_energy"), 
                            "cutoff":2*2**(1./6.)*pmb.units.nm, 
                            "offset":2*pmb.units.nm,
                            "acidity": "acidic",
                            "pka": 3}
        C_input_parameters={"name":"C", 
                        "sigma":0*pmb.units.nm, 
                        "epsilon":pmb.units.Quantity(2,"reduced_energy"), 
                        "cutoff":2*2**(1./6.)*pmb.units.nm, 
                        "offset":2*pmb.units.nm}
        pmb.define_particle(**A_input_parameters)
        pmb.define_particle(**B_input_parameters)
        pmb.define_particle(**C_input_parameters)
        # Setup LJ interactions shift="auto"
        pmb.setup_lj_interactions()
        # Check A-A LJ setup
        lj_templates = pmb.db.get_templates(pmb_type="lj")
        # Check B-B, B-BH, BH-BH setup
        labels=["A-A", "B-B", "B-BH", "BH-BH"]
        for label in labels:
            lj_template = lj_templates[label]
            if label == "A-A":
                input_params = A_input_parameters
            else:
                input_params = B_input_parameters
            for parameter_key in ["sigma","offset","cutoff"]:
                value_in_pyMBE = getattr(lj_template, parameter_key).to_quantity(pmb.units)
                self.assertEqual(first=value_in_pyMBE.to("reduced_length").magnitude, 
                                second=input_params[parameter_key].to("reduced_length").magnitude)
            self.assertAlmostEqual(first=lj_template.epsilon.to_quantity(pmb.units).to("reduced_energy").magnitude, 
                                second=input_params["epsilon"].to("reduced_energy").magnitude)
        # Clean LJ interactions
        pmb.db.delete_templates(pmb_type="lj")
        # ValueError if combining-rule other than Lorentz_-Berthelot is used
        input_params = { "combining_rule": "Geometric"}
        self.assertRaises(ValueError, pmb.setup_lj_interactions, **input_params)
        # Check initialization with shift=0
        pmb.setup_lj_interactions( shift_potential=False)
        # Calculate the reference parameters using Lorentz-Berthelot combining rule
        # Check A-BH, A-B, setup
        labels=["A-BH", "A-B"]
        ref_lj_parameters={}
        for parameter_key in ["sigma","offset","cutoff"]:
            ref_lj_parameters[parameter_key]=(A_input_parameters[parameter_key]+B_input_parameters[parameter_key])/2
        ref_lj_parameters["epsilon"]=np.sqrt(A_input_parameters["epsilon"]*B_input_parameters["epsilon"])
        lj_templates = pmb.db.get_templates(pmb_type="lj")
        lj_df = pmb.db._get_templates_df(pmb_type="lj")
        self.assertEqual(lj_df[lj_df.name == "A-A"]["shift"].values[0].m_as("nanometer"),
                         0)
        for label in labels:
            lj_template = lj_templates[label]
            for parameter_key in ["sigma","offset","cutoff"]:
                value_in_pyMBE = getattr(lj_template, parameter_key).to_quantity(pmb.units)
                self.assertEqual(first=value_in_pyMBE.to("reduced_length").magnitude, 
                                second=ref_lj_parameters[parameter_key].to("reduced_length").magnitude)
            self.assertAlmostEqual(first=lj_template.epsilon.to_quantity(pmb.units).to("reduced_energy").magnitude, 
                                second=ref_lj_parameters["epsilon"].to("reduced_energy").magnitude)
        # Check that no interaction between particle C and any other particle has been set up
        # Particle C has sigma = 0 (ideally behaving particle)
        for label in lj_templates.keys():
            self.assertFalse("C" in label)
        input_params = {"particle_name1":"A", 
                        "particle_name2":"B", 
                        "combining_rule":"Geometric"}
        self.assertRaises(ValueError, pmb.get_lj_parameters, **input_params)

if __name__ == "__main__":
    ut.main()