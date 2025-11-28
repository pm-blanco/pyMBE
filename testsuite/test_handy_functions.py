#
# Copyright (C) 2025 pyMBE-dev team
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

"""Test that the utility functions in lib/handy_functions work properly."""

import unittest as ut
import espressomd
import espressomd.version
import pyMBE
import pyMBE.lib.handy_functions as hf
import logging
import io
import re
import numpy as np
# Create an in-memory log stream
log_stream = io.StringIO()
logging.basicConfig(level=logging.INFO, 
                    format="%(levelname)s: %(message)s",
                    handlers=[logging.StreamHandler(log_stream)] )

# Create instances of espresso and pyMBE
espresso_system = espressomd.System(box_l=[60,60,60])
seed = 23
pmb = pyMBE.pymbe_library(seed=seed)
kT = pmb.kT

langevin_inputs={"espresso_system":espresso_system, 
                "kT" : kT, 
                "seed" : seed,
                "time_step" : 1e-2, 
                "gamma": 1, 
                "tune_skin" : False, 
                "min_skin": 1,
                "max_skin": None,
                "tolerance": 1e-3,
                "int_steps": 200,
                "adjust_max_skin": True}

relax_inputs={"espresso_system":espresso_system, 
              "gamma":0.01,
              "Nsteps_steepest_descent":5000, 
              "max_displacement":0.1, 
              "Nsteps_iter_relax":500,
              "seed": seed}

electrostatics_inputs={"units": pmb.units, 
                       "espresso_system": espresso_system, 
                       "kT": pmb.kT, 
                       "c_salt": None, 
                       "solvent_permittivity":78.5, 
                       "method": 'p3m', 
                       "tune_p3m":True, 
                       "accuracy":1e-3,
                       "verbose":False}


class Test(ut.TestCase):

    def setUp(self):
        espresso_system.time_step = 1e-2
        espresso_system.cell_system.skin = 5.

    def tearDown(self):
        espresso_system.part.clear()
        espresso_system.non_bonded_inter.reset()

    def test_exceptions_langevin_setup(self):
        """Test exceptions in :func:`lib.handy_functions.setup_langevin_dynamics`"""
        broken_inputs  = langevin_inputs.copy()
        broken_inputs["seed"] = "pyMBE"
        self.assertRaises(TypeError, hf.setup_langevin_dynamics, **broken_inputs)
        broken_inputs  = langevin_inputs.copy()
        broken_inputs["time_step"] = -1
        self.assertRaises(ValueError, hf.setup_langevin_dynamics, **broken_inputs)
        broken_inputs  = langevin_inputs.copy()
        broken_inputs["gamma"] = -1
        self.assertRaises(ValueError, hf.setup_langevin_dynamics, **broken_inputs)
        broken_inputs  = langevin_inputs.copy()
        broken_inputs["min_skin"] = 10
        broken_inputs["max_skin"] = 1
        self.assertRaises(ValueError, hf.setup_langevin_dynamics, **broken_inputs)

    def test_langevin_setup(self):
        """Test :func:`lib.handy_functions.setup_langevin_dynamics`"""
        hf.setup_langevin_dynamics(**langevin_inputs)
        ## Test setup of the integrator
        self.assertEqual(first=langevin_inputs["time_step"],
                         second=espresso_system.time_step,
                         msg="The input time step in `lib.handy_functions.setup_langevin_dynamics` is not consistent with the one in the espresso simulation System")
        ## Test setup of the thermostat
        if espressomd.version.friendly() == "4.2":
            thermostat_setup=espresso_system.thermostat.get_state()[0]
        else:
            thermostat_setup=espresso_system.thermostat.langevin.get_params()
            thermostat_setup["kT"] = espresso_system.thermostat.kT
            thermostat_setup["type"] = "LANGEVIN"
        self.assertEqual(first="LANGEVIN",
                         second=thermostat_setup["type"],
                         msg="`lib.handy_functions.setup_langevin_dynamics` is setting a different thermostat than Langevin")
        self.assertEqual(first=langevin_inputs["kT"].m_as("reduced_energy"),
                         second=thermostat_setup["kT"],
                         msg="`lib.handy_functions.setup_langevin_dynamics` is setting a different kT than the input one")
        self.assertEqual(first=langevin_inputs["seed"],
                         second=thermostat_setup["seed"],
                         msg="`lib.handy_functions.setup_langevin_dynamics` is setting a different seed than the input one")
        np.testing.assert_allclose(
             [langevin_inputs["gamma"]]*3,
             np.copy(thermostat_setup["gamma"]),
             err_msg="`lib.handy_functions.setup_langevin_dynamics` is setting a different seed than the input one")

        espresso_system.thermostat.turn_off()
        espresso_system.part.add(pos=[1,1,1])
        espresso_system.part.add(pos=[2,2,2])
        espresso_system.non_bonded_inter[0,0].lennard_jones.set_params(
            epsilon = 1, sigma = 1, cutoff = 3, shift = "auto")
        langevin_inputs_opt_skin  = langevin_inputs.copy() 
        langevin_inputs_opt_skin["tune_skin"] = True
        hf.setup_langevin_dynamics(**langevin_inputs_opt_skin)
        log_contents = log_stream.getvalue()
        optimized_skin = float(re.search(r"Optimized skin value: ([\d.]+)", log_contents).group(1))
        self.assertEqual(first=optimized_skin,
                         second=espresso_system.cell_system.skin,
                         msg="The optimized skin has not been set in espresso")

    def test_exceptions_relax_espresso_system(self):
        """Test exceptions in :func:`lib.handy_functions.relax_espresso_system`"""
        broken_inputs  = relax_inputs.copy()
        broken_inputs["gamma"] = -1
        self.assertRaises(ValueError, hf.relax_espresso_system, **broken_inputs)
        broken_inputs  = relax_inputs.copy()
        broken_inputs["Nsteps_steepest_descent"] = -1
        self.assertRaises(ValueError, hf.relax_espresso_system, **broken_inputs)
        broken_inputs  = relax_inputs.copy()
        broken_inputs["Nsteps_iter_relax"] = -1
        self.assertRaises(ValueError, hf.relax_espresso_system, **broken_inputs)
        broken_inputs  = relax_inputs.copy()
        broken_inputs["max_displacement"] = -1
        self.assertRaises(ValueError, hf.relax_espresso_system, **broken_inputs)
        broken_inputs  = relax_inputs.copy()

    def test_relax_espresso_system(self):
        """Test :func:`lib.handy_functions.relax_espresso_system`"""
        positions = [[1,1,1], 
                     [1.15,1.15,1.15], 
                     [1.5,1.5,1.5],
                     [2,2,2],
                     [2.15,2.15,2.15],
                     [1,1,1.5]]
        for pos in positions:
            espresso_system.part.add(pos=pos)
        espresso_system.non_bonded_inter[0,0].lennard_jones.set_params(
            epsilon = 1, sigma = 1, cutoff = 2**(1./6.), shift = "auto")
        min_dist = hf.relax_espresso_system(**relax_inputs)
        self.assertGreater(a=min_dist,
                           b=0.9,
                           msg="lib.handy_functions.relax_espresso_system is unable to relax a simple lj system")
        espresso_system.part.clear()
        # Test that the LD run is actually performed
        for pos in positions:
            espresso_system.part.add(pos=pos)
        test_inputs  = relax_inputs.copy()
        test_inputs["Nsteps_steepest_descent"] = 1
        test_inputs["max_displacement"] = 0
        hf.relax_espresso_system(**relax_inputs)
        new_positions = espresso_system.part.all().pos
        distances = np.linalg.norm(np.array(positions) - new_positions, 
                                   axis=1)
        # Check that the positions have changed
        self.assertTrue(np.all(distances > 0.1),
                        msg="lib.handy_functions.relax_espresso_system does not change the positions of the particles during relaxation")



    def test_exceptions_electrostatics(self):
        """Test exceptions in :func:`lib.handy_functions.setup_electrostatic_interactions`"""
        broken_inputs  = electrostatics_inputs.copy()
        broken_inputs["units"] = "pyMBE"
        self.assertRaises(TypeError, hf.setup_electrostatic_interactions, **broken_inputs)
        broken_inputs  = electrostatics_inputs.copy()
        broken_inputs["method"] = "dH"
        self.assertRaises(ValueError, hf.setup_electrostatic_interactions, **broken_inputs)
        broken_inputs  = electrostatics_inputs.copy()
        broken_inputs["method"] = "dh"
        self.assertRaises(ValueError, hf.setup_electrostatic_interactions, **broken_inputs)
        broken_inputs  = electrostatics_inputs.copy()
        broken_inputs["c_salt"] = 10*pmb.units.nm
        self.assertRaises(ValueError, hf.setup_electrostatic_interactions, **broken_inputs)

    def test_setup_electrostatics(self):
        """Test :func:`lib.handy_functions.setup_electrostatic_interactions`"""
        espresso_system.part.add(pos=[1,1,1], q=1)
        espresso_system.part.add(pos=[5.15,5.15,5.15], q=-1)
        Bjerrum_length = pmb.e.to('reduced_charge')**2 / (4 * pmb.units.pi * pmb.units.eps0 * electrostatics_inputs["solvent_permittivity"] * electrostatics_inputs["kT"].to('reduced_energy'))
        coulomb_prefactor=Bjerrum_length*electrostatics_inputs["kT"]
        # Test the P3M setup
        hf.setup_electrostatic_interactions(**electrostatics_inputs)
        if espressomd.version.friendly() == "4.2":
            coulomb = espresso_system.actors.active_actors.copy()[0]
        else:
            coulomb = espresso_system.electrostatics.solver
        coulomb_params = coulomb.get_params()
        self.assertEqual(first=coulomb.name(),
                         second='Coulomb::CoulombP3M',
                         msg="lib.handy_functions.setup_electrostatic_interactions sets up the wrong electrostatic method")
        self.assertGreaterEqual(a=electrostatics_inputs["accuracy"],
                                b=coulomb.accuracy,
                                msg="lib.handy_functions.setup_electrostatic_interactions sets up the P3M method with the wrong accuracy")
        self.assertAlmostEqual(first=coulomb_params["prefactor"],
                                second=coulomb_prefactor.m_as("reduced_length * reduced_energy"),
                                msg="lib.handy_functions.setup_electrostatic_interactions sets up the wrong coulomb prefactor for the P3M method")
        self.assertEqual(first=electrostatics_inputs["tune_p3m"],
                         second=coulomb_params["is_tuned"],
                         msg="lib.handy_functions.setup_electrostatic_interactions does not tune the P3M method")
        if espressomd.version.friendly() == "4.2":
            espresso_system.actors.remove(coulomb)
        else:
            coulomb = espresso_system.electrostatics.solver = None
        ## Test the setup of the P3M method without tuning it with some input parameters
        electrostatics_inputs["tune_p3m"] = False
        electrostatics_inputs["params"] = {"mesh": [8, 8, 8], 
                                           "cao": 5, 
                                           "alpha": 1.1265e+01,
                                           "r_cut": 1}
        hf.setup_electrostatic_interactions(**electrostatics_inputs)
        if espressomd.version.friendly() == "4.2":
            coulomb = espresso_system.actors.active_actors.copy()[0]
        else:
            coulomb = espresso_system.electrostatics.solver
        coulomb_params = coulomb.get_params()
        for param in electrostatics_inputs["params"]:
            np.testing.assert_allclose(
                 np.copy(electrostatics_inputs["params"][param]),
                 np.copy(coulomb_params[param]),
                 err_msg="lib.handy_functions.setup_electrostatic_interactions sets up the wrong P3M parameters")
        if espressomd.version.friendly() == "4.2":
            espresso_system.actors.remove(coulomb)
        else:
            coulomb = espresso_system.electrostatics.solver = None
        electrostatics_inputs["params"] = None
        # Test the Debye–Hückel setup
        electrostatics_inputs["method"] = "dh"
        electrostatics_inputs["c_salt"] = pmb.units.Quantity(1, "mol/L")
        kappa=1./np.sqrt(8*pmb.units.pi*Bjerrum_length*pmb.N_A*electrostatics_inputs["c_salt"])
        hf.setup_electrostatic_interactions(**electrostatics_inputs)
        if espressomd.version.friendly() == "4.2":
            dh = espresso_system.actors.active_actors.copy()[0]
        else:
            dh = espresso_system.electrostatics.solver
        dh_params = dh.get_params()
        self.assertEqual(first=dh.name(),
                         second='Coulomb::DebyeHueckel',
                         msg="lib.handy_functions.setup_electrostatic_interactions sets up the wrong electrostatic method")
        self.assertAlmostEqual(first=dh_params["prefactor"],
                                second=coulomb_prefactor.m_as("reduced_length * reduced_energy"),
                                msg="lib.handy_functions.setup_electrostatic_interactions sets up the wrong coulomb prefactor for the DH method")
        self.assertAlmostEqual(first=dh_params["kappa"],
                                second=(1./kappa).m_as('1/ reduced_length'),
                                msg="lib.handy_functions.setup_electrostatic_interactions sets up the wrong Debye screening length for the DH method")
        self.assertAlmostEqual(first=dh_params["r_cut"],
                                second=3*kappa.m_as('reduced_length'),
                                msg="lib.handy_functions.setup_electrostatic_interactions sets up the wrong cut-off for the DH method")
        if espressomd.version.friendly() == "4.2":
            espresso_system.actors.remove(dh)
        else:
            coulomb = espresso_system.electrostatics.solver = None
        electrostatics_inputs["c_salt"] = pmb.units.Quantity(1, "mol/L")*pmb.N_A
        hf.setup_electrostatic_interactions(**electrostatics_inputs)
        if espressomd.version.friendly() == "4.2":
            dh = espresso_system.actors.active_actors.copy()[0]
        else:
            dh = espresso_system.electrostatics.solver
        dh_params = dh.get_params()
        self.assertAlmostEqual(first=dh_params["kappa"],
                                second=(1./kappa).m_as('1/ reduced_length'),
                                msg="lib.handy_functions.setup_electrostatic_interactions sets up the wrong Debye screening length for the DH method")
        self.assertAlmostEqual(first=dh_params["r_cut"],
                                second=3*kappa.m_as('reduced_length'),
                                msg="lib.handy_functions.setup_electrostatic_interactions sets up the wrong cut-off for the DH method")
        if espressomd.version.friendly() == "4.2":
            espresso_system.actors.remove(dh)
        else:
            coulomb = espresso_system.electrostatics.solver = None
        # Test a non-default cut-off
        electrostatics_inputs["params"] = {"r_cut": 3}
        hf.setup_electrostatic_interactions(**electrostatics_inputs)
        if espressomd.version.friendly() == "4.2":
            dh = espresso_system.actors.active_actors.copy()[0]
        else:
            dh = espresso_system.electrostatics.solver
        dh_params = dh.get_params()
        self.assertAlmostEqual(first=dh_params["r_cut"],
                                second=electrostatics_inputs["params"]["r_cut"],
                                msg="lib.handy_functions.setup_electrostatic_interactions sets up the wrong cut-off for the DH method")

if __name__ == "__main__":
    ut.main()