#
# Copyright (C) 2024-2025 pyMBE-dev team
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

import logging

def setup_electrostatic_interactions(units, espresso_system, kT, c_salt=None, solvent_permittivity=78.5, method='p3m', tune_p3m=True, accuracy=1e-3, params=None, verbose=False):
    """
    Sets up electrostatic interactions in an ESPResSo system.

    Args:
        units(`pint.UnitRegistry`): Unit registry for handling physical units.
        espresso_system(`espressomd.system.System`): system object of espressomd library.
        kT(`pint.Quantity`): Thermal energy.
        c_salt(`pint.Quantity`): Added salt concentration. If provided, the program outputs the debye screening length. It is a mandatory parameter for the Debye-Hückel method.
        solvent_permittivity (`float`): Solvent relative permittivity. Defaults to 78.5, correspoding to its value in water at 298.15 K.
        method(`str`): Method for computing electrostatic interactions. Defaults to "p3m". 
        tune_p3m(`bool`): If True, tunes P3M parameters for efficiency. Defaults to True. 
        accuracy(`float`): Desired accuracy for electrostatics. Defaults to 1e-3.
        params(`dict`): Additional parameters for the electrostatic method. For P3M, it can include 'mesh', 'alpha', 'cao' and `r_cut`. For Debye-Hückel, it can include 'r_cut'.
        verbose(`bool`): If True, enables verbose output for P3M tuning. Defaults to False.

    Note:
        `c_salt` is a mandatory argument for setting up the Debye-Hückel electrostatic potential.
        The calculated Bjerrum length is ouput to the log. If `c_salt` is provided, the calculated Debye screening length is also output to the log.
        Currently, the only supported electrostatic methods are P3M ("p3m") and Debye-Hückel ("dh").
    """
    import espressomd.electrostatics
    import espressomd.version
    import numpy as np
    import scipy.constants
    logging.debug("*** Starting electrostatic interactions setup... ***")
    # Initial sanity checks
    if not hasattr(units, 'Quantity'):
        raise TypeError("Invalid 'units' argument: Expected a pint.UnitRegistry object")
    valid_methods_list=['p3m', 'dh']
    if method not in valid_methods_list:
        raise ValueError('Method not supported, supported methods are', valid_methods_list)
    if c_salt is None and method == 'dh':
        raise ValueError('Please provide the added salt concentration c_salt to setup the Debye-Huckel potential')
    e = scipy.constants.e * units.C
    N_A = scipy.constants.N_A / units.mol
    BJERRUM_LENGTH = e**2 / (4 * units.pi * units.eps0 * solvent_permittivity * kT)
    logging.info(f" Bjerrum length {BJERRUM_LENGTH.to('nm')} = {BJERRUM_LENGTH.to('reduced_length')}")
    COULOMB_PREFACTOR=BJERRUM_LENGTH * kT 
    if c_salt is not None:
        if c_salt.check('[substance] [length]**-3'):
            KAPPA=1./np.sqrt(8*units.pi*BJERRUM_LENGTH*N_A*c_salt)
        elif c_salt.check('[length]**-3'):
            KAPPA=1./np.sqrt(8*units.pi*BJERRUM_LENGTH*c_salt)
        else:
            raise ValueError('Unknown units for c_salt, supported units for salt concentration are [mol / volume] or [particle / volume]', c_salt)
        
        logging.info(f"Debye kappa {KAPPA.to('nm')} = {KAPPA.to('reduced_length')}")

    if params is None:
        params = {}

    if method == 'p3m':
        logging.debug("*** Setting up Coulomb electrostatics using the P3M method ***")
        coulomb = espressomd.electrostatics.P3M(prefactor = COULOMB_PREFACTOR.m_as("reduced_length * reduced_energy"), 
                                                accuracy=accuracy,
                                                verbose=verbose,
                                                tune=tune_p3m,
                                                **params)

        if tune_p3m:
            espresso_system.time_step=0.01
            if espressomd.version.friendly() == "4.2":
                espresso_system.actors.add(coulomb)
            else:
                espresso_system.electrostatics.solver = coulomb


            # save the optimal parameters and add them by hand

            p3m_params = coulomb.get_params()
            if espressomd.version.friendly() == "4.2":
                espresso_system.actors.remove(coulomb)
            else:
                espresso_system.electrostatics.solver = None
            coulomb = espressomd.electrostatics.P3M(prefactor = COULOMB_PREFACTOR.m_as("reduced_length * reduced_energy"),
                                                    accuracy = accuracy,
                                                    mesh = p3m_params['mesh'],
                                                    alpha = p3m_params['alpha'] ,
                                                    cao = p3m_params['cao'],
                                                    r_cut = p3m_params['r_cut'],
                                                    tune = False)

    elif method == 'dh':
        logging.debug("*** Setting up Debye-Hückel electrostatics ***")
        if params:
            r_cut = params['r_cut']
        else:
            r_cut = 3*KAPPA.to('reduced_length').magnitude
            
        coulomb = espressomd.electrostatics.DH(prefactor = COULOMB_PREFACTOR.m_as("reduced_length * reduced_energy"), 
                                               kappa = (1./KAPPA).to('1/ reduced_length').magnitude, 
                                               r_cut = r_cut)
    if espressomd.version.friendly() == "4.2":
        espresso_system.actors.add(coulomb)
    else:
        espresso_system.electrostatics.solver = coulomb
    logging.debug("*** Electrostatics successfully added to the system ***")
    return

def relax_espresso_system(espresso_system, seed, gamma=1e-3, Nsteps_steepest_descent=5000, max_displacement=0.01, Nsteps_iter_relax=500):
    """
    Relaxes the energy of the given ESPResSo system by performing the following steps:
    (1) Steepest descent energy minimization, to remove large forces and relax the system to a local minimum.
    (2) A Langevin Dynamics run, to further relax the system and ensure that it is in thermal equilibrium.

    This function is useful to avoid code repetition in the sample scripts of pyMBE, but it is by no means general-purpose.
    Similarly, the default parameters are not universal and should be adapted to the specific system at hand.
    In general, system relaxation is a complex procedure and should be adapted for each particular application.
    If you experience crashes or unexpected behavior, please consider using your own relaxation procedure.

    Args:
        espresso_system (`espressomd.system.System`): system object of espressomd library.
        seed (`int`): Seed for the random number generator for the thermostat.
        gamma (`float`, optional): Starting damping constant for Langevin dynamics. Defaults to  1e-3 reduced time**-1.
        Nsteps_steepest_descent (`int`, optional): Total number of steps for steepest descent minimization. Defaults to 5000.
        max_displacement (`float`, optional): Maximum particle displacement allowed during minimization. Defaults to 0.01 reduced length.
        Nsteps_iter_relax (`int`, optional): Number of steps per iteration for Langevin dynamics relaxation. Defaults to 500.

    Return:
        (`float`): minimum distance between particles in the system after the relaxation

    Note:
        The thermostat is turned off by the end of the procedure. 
        Make sure the system is initialized properly before calling this function.
    """
    # Sanity checks
    if gamma <= 0:
        raise ValueError("The damping constant 'gamma' must be positive.")
    if Nsteps_steepest_descent <= 0 or Nsteps_iter_relax <= 0:
        raise ValueError("Step counts must be positive integers.")
    if max_displacement <= 0:
        raise ValueError("'max_displacement' must be positive.")
    logging.debug("*** Relaxing the energy of the system... ***")
    logging.debug("*** Starting steepest descent minimization ***")
    espresso_system.thermostat.turn_off()
    espresso_system.integrator.set_steepest_descent(f_max=0,
                                                    gamma=gamma, 
                                                    max_displacement=max_displacement)
    espresso_system.integrator.run(Nsteps_steepest_descent)
    logging.debug("*** Finished steepest descent minimization ***")
    logging.debug("*** Starting Langevin Dynamics relaxation ***")
    
    espresso_system.integrator.set_vv()
    espresso_system.thermostat.set_langevin(kT=1., gamma=gamma, seed=seed)
    espresso_system.integrator.run(Nsteps_iter_relax)
    espresso_system.thermostat.turn_off()

    logging.debug("*** Finished Langevin Dynamics relaxation ***")
    logging.info(f"*** Minimum particle distance after relaxation: {espresso_system.analysis.min_dist()} ***")
    logging.debug("*** Relaxation finished ***")
    return espresso_system.analysis.min_dist()

def setup_langevin_dynamics(espresso_system, kT, seed,time_step=1e-2, gamma=1, tune_skin=True, min_skin=1, max_skin=None, tolerance=1e-3, int_steps=200, adjust_max_skin=True):
    """
    Sets up Langevin Dynamics for an ESPResSo simulation system.

    Args:
        espresso_system (`espressomd.system.System`): system object of espressomd library.
        kT (`pint.Quantity`): Target temperature in reduced energy units.
        seed (`int`): Seed for the random number generator for the thermostat.
        time_step (`float`, optional): Integration time step. Defaults to 1e-2.
        gamma (`float`, optional): Damping coefficient for the Langevin thermostat. Defaults to 1.
        tune_skin (`bool`, optional): Whether to optimize the skin parameter. Defaults to True.
        min_skin (`float`, optional): Minimum skin value for optimization. Defaults to 1.
        max_skin (`float`, optional): Maximum skin value for optimization. Defaults to None, which is handled by setting its value to box length / 2.
        tolerance (`float`, optional): Tolerance for skin optimization. Defaults to 1e-3.
        int_steps (`int`, optional): Number of integration steps for tuning. Defaults to 200.
        adjust_max_skin (`bool`, optional): Whether to adjust the maximum skin value during tuning. Defaults to True.
    """        
    if not isinstance(seed, int):
        raise TypeError("seed must be an integer.")
    if not isinstance(time_step, (float, int)) or time_step <= 0:
        raise ValueError("time_step must be a positive number.")
    if not isinstance(gamma, (float, int)) or gamma <= 0:
        raise ValueError("gamma must be a positive number.")
    if max_skin is None:
        max_skin=espresso_system.box_l[0]/2
    if min_skin >= max_skin:
        raise ValueError("min_skin must be smaller than max_skin.")
    espresso_system.time_step=time_step
    espresso_system.integrator.set_vv()
    espresso_system.thermostat.set_langevin(kT= kT.to('reduced_energy').magnitude, 
                                            gamma= gamma, 
                                            seed= seed)
    # Optimize the value of skin
    if tune_skin:
        logging.debug("*** Optimizing skin ... ***")
        espresso_system.cell_system.tune_skin(min_skin=min_skin, 
                                              max_skin=max_skin, 
                                              tol=tolerance, 
                                              int_steps=int_steps, 
                                              adjust_max_skin=adjust_max_skin)
        logging.info(f"Optimized skin value: {espresso_system.cell_system.skin}")
    return

def get_number_of_particles(espresso_system, ptype):
    import espressomd.version
    if espressomd.version.friendly() == "4.2":
        args = (ptype,)
        kwargs = {}
    else:
        args = ()
        kwargs = {"type": ptype}
    return espresso_system.number_of_particles(*args, **kwargs)

def do_reaction(algorithm, steps):
    import espressomd.version
    if espressomd.version.friendly() == '4.2':
        algorithm.reaction(reaction_steps=steps)
    else:
        algorithm.reaction(steps=steps)
