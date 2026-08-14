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

import espressomd
import espressomd.electrostatics
import espressomd.version
import warnings
from  typing import List
import numpy as np
import logging
from pyMBE.simulation_builder.base_engine import SimulationEngine
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
from pyMBE.storage.templates.lj import LJInteractionTemplate
from pyMBE.storage.pint_quantity import PintQuantity

class EspressoSimulation(SimulationEngine):
    def __init__(self,box_l,db,espresso_system,units,kT,Kw,seed):
        self.db=db
        self.box_l: List[float]=box_l
        self.espresso_system=espresso_system
        self.units=units
        self.kT=kT
        self.Kw=Kw
        self.seed=seed
        pass

    def _add_angle(self,particle_id1,particle_id2,particle_id3, angle_inst):
        """ 
        Helper function to add angle instances to ESPResSo.

        Args:
            particle_id1 (int): pid of the particle instance 1 that composes the angle 
            particle_id2 (int): pid of the particle instance 2 that composes the angle 
            particle_id3 (int): pid of the particle instance 3 that composes the angle 
            angle_inst (pmb.AngleInstance):  dataclass containing information of an angle instance
        """
        angle_tpl=self.db.get_template(name=angle_inst.name, 
                                        pmb_type="angle")
        espresso_angle_inst=self._get_angle_instance(angle_template=angle_tpl)

        self.espresso_system.part.by_id(particle_id2).add_bond((espresso_angle_inst, particle_id1, particle_id3))
        self.db._update_instance(instance_id=angle_inst.angle_id,
                                    pmb_type='angle',
                                    attribute='added_to_engine',
                                    value=True)
        
    def _add_bond(self,particle_id1,particle_id2,bond_inst):
        """
        Helper function to add bond instances to ESPResSo.

        Args:
            particle_id1 (int): pid of the particle instance 1 that composes the bond 
            particle_id2 (int): pid of the particle instance 2 that composes the bond
            bond_inst (pmb.BondInstance): dataclass containing information of a bond instance
        """
        bond_tpl=self.db.get_template(name=bond_inst.name, 
                                        pmb_type="bond")
        espresso_bond_inst=self._get_bond_instance(bond_template=bond_tpl)
        self.espresso_system.part.by_id(particle_id1).add_bond((espresso_bond_inst, particle_id2))
        self.db._update_instance(instance_id=bond_inst.bond_id,
                                    pmb_type='bond',
                                    attribute='added_to_engine',
                                    value=True)
    
    def _add_particle(self,particle_id):
        """
        Helper function to add particle instances to ESPResSo

        Args:
            particle_id (int): pid of the particle instance
        """
        particle_instance=self.db.get_instance(pmb_type='particle',
                                 instance_id=particle_id)
        part_state = self.db.get_template(pmb_type="particle_state",
                                        name=particle_instance.initial_state)
        kwargs = dict(id=particle_id, 
                        pos=particle_instance.position, 
                        type=part_state.es_type, 
                        q=part_state.z,
                        fix=particle_instance.fix)        
        self.espresso_system.part.add(**kwargs)
        self.db._update_instance(instance_id=particle_id,
                                    pmb_type='particle',
                                    attribute='added_to_engine',
                                    value=True)
    
    def _check_bond_inputs(self, bond_type, bond_parameters):
        """
        Checks that the input bond parameters are valid within the current pyMBE implementation.

        Args:
            bond_type ('str'): 
                label to identify the potential to model the bond.
            
            bond_parameters ('dict'): 
                parameters of the potential of the bond.
        """
        valid_bond_types   = ["harmonic", "FENE"] 
        if bond_type not in valid_bond_types:
            raise NotImplementedError(f"Bond type '{bond_type}' currently not implemented in pyMBE, accepted types are {valid_bond_types}")
        required_parameters = {"harmonic": ["r_0","k"],
                                "FENE": ["r_0","k","d_r_max"]}
        for required_parameter in required_parameters[bond_type]:
            if required_parameter not in bond_parameters.keys():
                raise ValueError(f"Missing required parameter {required_parameter} for {bond_type} bond")
    
    def _create_bond_instance(self, bond_type, bond_parameters):
        """
        Creates an ESPResSo bond instance.

        Args:
            bond_type ('str'): 
                label to identify the potential to model the bond.

            bond_parameters ('dict'): 
                parameters of the potential of the bond.

        Notes:
            Currently, only HARMONIC and FENE bonds are supported.

            For a HARMONIC bond the dictionary must contain:
                - k ('Pint.Quantity')      : Magnitude of the bond. It should have units of energy/length**2 
                using the 'pmb.units' UnitRegistry.
                - r_0 ('Pint.Quantity')    : Equilibrium bond length. It should have units of length using 
                the 'pmb.units' UnitRegistry.
           
            For a FENE bond the dictionary must additionally contain:
                - d_r_max ('Pint.Quantity'): Maximal stretching length for FENE. It should have 
                units of length using the 'pmb.units' UnitRegistry. Default 'None'.

        Returns:
            ('espressomd.interactions'): instance of an ESPResSo bond object
        """
        
        self._check_bond_inputs(bond_parameters=bond_parameters,
                                bond_type=bond_type)
        if bond_type == 'harmonic':
            bond_instance = espressomd.interactions.HarmonicBond(k = bond_parameters["k"].m_as("reduced_energy/reduced_length**2"),
                                                                r_0 = bond_parameters["r_0"].m_as("reduced_length"))
        elif bond_type == 'FENE':
            bond_instance = espressomd.interactions.FeneBond(k = bond_parameters["k"].m_as("reduced_energy/reduced_length**2"),
                                                            r_0 = bond_parameters["r_0"].m_as("reduced_length"),
                                                            d_r_max = bond_parameters["d_r_max"].m_as("reduced_length"))    
        return bond_instance
    
    def _delete_particles(self, particle_ids):
        """
        Remove a list of particles from an ESPResSo simulation system.

        Args:
            particle_ids  ('Iterable[int]'):
                A list (or other iterable) of ESPResSo particle IDs to remove.


        Notess:
            - This method removes particles only from the ESPResSo simulation,
            **not** from the pyMBE database. Database cleanup must be handled
            separately by the caller.
            - Attempting to remove a non-existent particle ID will raise
            an ESPResSo error.
        """
        for pid in particle_ids:
            self.espresso_system.part.by_id(pid).remove()
            self.db._update_instance(instance_id=pid,
                                     pmb_type='particle',
                                     attribute='added_to_engine',
                                     value=False)

    def _get_bond_instance(self, bond_template):
        """
        Retrieve or create a bond instance in an ESPResSo system for a given pair of particle names.

        Args:
            bond_template ('BondTemplate'): 
                BondTemplate object from the pyMBE database.

        Returns:
            ('espressomd.interactions.BondedInteraction'): 
                The ESPResSo bond instance object.

        Notes:
            When a new bond instance is created, it is not added to the ESPResSo system.
        """
        if bond_template.name in self.db.espresso_bond_instances.keys():
            bond_inst = self.db.espresso_bond_instances[bond_template.name]
        else:   
            # Create an instance of the bond 
            bond_inst = self._create_bond_instance(bond_type=bond_template.bond_type,
                                                   bond_parameters=bond_template.get_parameters(self.units))
            self.db.espresso_bond_instances[bond_template.name]= bond_inst
            self.espresso_system.bonded_inter.add(bond_inst)
        return bond_inst
    
    def _get_angle_instance(self,angle_template):
        """
        Retrieve or create an angle interaction in an ESPResSo system for a given angle template.

        Args:
            angle_template ('AngleTemplate'): The angle template to use.

        Returns:
            ('espressomd.interactions.BondedInteraction'): The ESPResSo angle interaction object.
        """
        if angle_template.name in self.db.espresso_angle_instances:
                return self.db.espresso_angle_instances[angle_template.name]
    
        angle_inst = self._create_angle_instance(angle_type=angle_template.angle_type,
                                                    angle_parameters=angle_template.get_parameters(self.units))
        self.db.espresso_angle_instances[angle_template.name] = angle_inst
        self.espresso_system.bonded_inter.add(angle_inst)
    
        return angle_inst 
    
    def _create_angle_instance(self, angle_type, angle_parameters):
        """
        Creates an ESPResSo angle interaction object.

        Args:
            angle_type ('str'): Type of angle potential ("harmonic", "cosine", "harmonic_cosine").
            angle_parameters ('dict'): Parameters of the angle potential (k, phi_0).

        Returns:
            ('espressomd.interactions.BondedInteraction'): The ESPResSo angle interaction object.
        """
        from espressomd import interactions

        k = angle_parameters["k"].m_as("reduced_energy")
        phi_0 = float(angle_parameters["phi_0"].magnitude)

        if angle_type == "harmonic":
            return interactions.AngleHarmonic(bend=k, phi0=phi_0)
        elif angle_type == "cosine":
            return interactions.AngleCosine(bend=k, phi0=phi_0)
        elif angle_type == "harmonic_cosine":
            return interactions.AngleCossquare(bend=k, phi0=phi_0)

    
    def _get_particle_ids_in_espresso(self):
        """
        Gets a list of all the particles_id in espresso_system instance.
        
        Returns:
            espresso_particles_id(list): list of pids of the particles that are saved in ESPResSo
        """
        espresso_particles=self.espresso_system.part.all()
        return espresso_particles.id    
    
    
    def _get_particle_pos_espresso(self,id):
        """
        Gets the position of particle with `id` in ESPResSo.

        Args:
            id (int): pid of the particle that we want to check that exists within ESPResSo.

        Returns:
            espresso_particles_id(List[List[float,float,float]]): returns a nested list of floats containing x,y,z coordinates for each particle in ESPResSo
        """
        return self.espresso_system.part.by_id(id).pos
    
    def _update_particle_position(self, particle_id, position):
        """
        Updates the position of a particle already stored in ESPResSo.

        Args:
            particle_id ('int'):
                ESPResSo particle ID to update.

            position ('array-like'):
                New Cartesian coordinates of the particle.
        """
        self.espresso_system.part.by_id(particle_id).pos = position

    def get_box_side_length(self):
        """
        Gets the dimensions of the simulation box used in the simulation engine.
        
        Returns:
            box_l(list[float,float,float]): return a list of floats regarding the dimensions of the box
        """
        return self.box_l
    
    def calculate_net_charge(self,object_name,pmb_type,dimensionless=False):
        """
        Calculates the net charge per instance of a given pmb object type.

        Args:
            object_name (str):
                Name of the object (e.g. molecule, residue, peptide, protein).
            pmb_type (str):
                Type of object to analyze. Must be molecule-like.
            dimensionless (bool, optional):
                If True, return charge as a pure number.
                If False, return a quantity with reduced_charge units.

        Returns:
            dict:
                {"mean": mean_net_charge, "instances": {instance_id: net_charge}}
        """
        id_map = self.db.get_particle_id_map(object_name=object_name)
        label = self.db._get_label_id_map(pmb_type=pmb_type)
        instance_map = id_map[label]
        charges = {}
        for instance_id, particle_ids in instance_map.items():
            if dimensionless:
                net_charge = 0.0
            else:
                net_charge = 0 * self.units.Quantity(1, "reduced_charge")
            for pid in particle_ids:
                q = self.espresso_system.part.by_id(pid).q
                if not dimensionless:
                    q *= self.units.Quantity(1, "reduced_charge")
                net_charge += q
            charges[instance_id] = net_charge
        # Mean charge
        if dimensionless:
            mean_charge = float(np.mean(list(charges.values())))
        else:
            mean_charge = (np.mean([q.magnitude for q in charges.values()])* self.units.Quantity(1, "reduced_charge"))
        return {"mean": mean_charge, "instances": charges}
    
    def change_volume_and_rescale_particles(self, d_new, dir="xyz"):
        """
        Change the volume for a particular dimension into the ESPResSo system.

        Args:
            d_new(float): new value for the dimension
            dir(Literal[x,y,z]): coordinate in which to set the new dimension. 
        
        """
        rescale_factor=np.array([1,1,1])
        if d_new<=0:
            raise ValueError("The dimension cannot be negative, neither 0")
        if "x" in dir:
            rescale_factor[0]=d_new/self.box_l[0]
        if "y" in dir:
            rescale_factor[1]=d_new/self.box_l[1]
        if "z" in dir:
            rescale_factor[2]=d_new/self.box_l[2]

        instances=self.db._get_instances_df(pmb_type='particle')
        for pid in range(instances.index.size):
            es_pos=self.db.get_instance(instance_id=pid,
                                        pmb_type='particle').position
            rescaled_position=es_pos*rescale_factor
            self.db._update_instance(instance_id=pid,
                                     pmb_type='particle',
                                     attribute='position',
                                     value=rescaled_position)
            
        self.espresso_system.change_volume_and_rescale_particles(d_new=d_new,
                                                                 dir=dir)
    
    
    def do_reaction(self,algorithm, steps):
        """
        Executes reaction steps using an ESPResSo reaction algorithm with
        version-compatible calling semantics.

        This function wraps the `reaction` method of an ESPResSo reaction
        algorithm to account for differences in the method signature between
        ESPResSo versions.

        Args:
            algorithm ('espressomd.reaction_methods'):
                ESPResSo reaction algorithm object (e.g. constant pH,
                reaction ensemble, or similar).
            steps ('int'):
                Number of reaction steps to perform.

        Notes:
            - In ESPResSo 4.2, the `reaction` method expects the number of steps
            to be passed as the keyword argument `reaction_steps`.
            - In newer ESPResSo versions, the keyword argument is `steps`.
            - This helper function provides a stable interface across ESPResSo
            versions by dispatching to the appropriate keyword internally.
        """
        import espressomd.version
        if espressomd.version.friendly() == '4.2':
            algorithm.reaction(reaction_steps=steps)
        else:
            algorithm.reaction(steps=steps)

    def enable_motion_of_rigid_object(self, instance_id, pmb_type):
        """
        Enables translational and rotational motion of a rigid pyMBE object instance
        in an ESPResSo system.This method creates a rigid-body center particle at the center of mass of
        the specified pyMBE object and attaches all constituent particles to it
        using ESPResSo virtual sites. The resulting rigid object can translate and
        rotate as a single body.

        Args:
            instance_id ('int'):
                Instance ID of the pyMBE object whose rigid-body motion is enabled.

            pmb_type ('str'):
                pyMBE object type of the instance (e.g. '"molecule"', '"peptide"',
                '"protein"', or any assembly-like type).

        Notess:
            - This method requires ESPResSo to be compiled with the following
            features enabled:
                - '"VIRTUAL_SITES_RELATIVE"'
                - '"MASS"'
            - A new ESPResSo particle is created to represent the rigid-body center.
            - The mass of the rigid-body center is set to the number of particles
            belonging to the object.
            - The rotational inertia tensor is approximated from the squared
            distances of the particles to the center of mass.
        """
        logging.info('enable_motion_of_rigid_object requires that espressomd has the following features activated: ["VIRTUAL_SITES_RELATIVE", "MASS"]')
        inst = self.db.get_instance(pmb_type=pmb_type,
                                    instance_id=instance_id)
        label = self.db._get_label_id_map(pmb_type=pmb_type)
        particle_ids_list = self.db.get_particle_id_map(object_name=inst.name)[label][instance_id]
        center_of_mass = self.calculate_center_of_mass (instance_id=instance_id,
                                                        pmb_type=pmb_type)
        rigid_object_center = self.espresso_system.part.add(pos=center_of_mass,
                                                        rotation=[True,True,True], 
                                                        type=self.db.propose_unused_type())
        rigid_object_center.mass = len(particle_ids_list)
        momI = 0
        for pid in particle_ids_list:
            momI += np.power(np.linalg.norm(center_of_mass - self.espresso_system.part.by_id(pid).pos), 2)
        rigid_object_center.rinertia = np.ones(3) * momI        
        for particle_id in particle_ids_list:
            pid = self.espresso_system.part.by_id(particle_id)
            pid.vs_auto_relate_to(rigid_object_center.id)

    def get_number_of_particles(self, ptype):
        """
        Returns the number of particles of a given ESPResSo particle type.

        Args:
            ptype ('int'):
                ESPResSo particle type identifier.

        Returns:
            ('int'):
                Number of particles in `espresso_system` with particle type `ptype`.

        Notes:
            - In ESPResSo 4.2, `number_of_particles` expects the particle type
            as a positional argument.
            - In later ESPResSo versions, the particle type must be passed as a
            keyword argument (`type=ptype`).
            - This helper function hides these API differences and provides
            a uniform interface across ESPResSo versions.
        """
        import espressomd.version
        if espressomd.version.friendly() == "4.2":
            args = (ptype,)
            kwargs = {}
        else:
            args = ()
            kwargs = {"type": ptype}
        return self.espresso_system.number_of_particles(*args, **kwargs)
    
    def relax_espresso_system(self, seed, gamma=1e-3, Nsteps_steepest_descent=5000, max_displacement=0.01, Nsteps_iter_relax=500):
        """
        Relaxes the energy of the given ESPResSo system by performing the following steps:
        (1) Steepest descent energy minimization, to remove large forces and relax the system to a local minimum.
        (2) A Langevin Dynamics run, to further relax the system and ensure that it is in thermal equilibrium.

        This function is useful to avoid code repetition in the sample scripts of pyMBE, but it is by no means general-purpose.
        Similarly, the default parameters are not universal and should be adapted to the specific system at hand.
        In general, system relaxation is a complex procedure and should be adapted for each particular application.
        If you experience crashes or unexpected behavior, please consider using your own relaxation procedure.

        Args:

            seed (`int`): 
                Seed for the random number generator for the thermostat.

            gamma (`float`, optional): 
                Starting damping constant for Langevin dynamics. Defaults to  1e-3 reduced time**-1.

            Nsteps_steepest_descent (`int`, optional): 
                Total number of steps for steepest descent minimization. Defaults to 5000.

            max_displacement (`float`, optional): 
                Maximum particle displacement allowed during minimization. Defaults to 0.01 reduced length.

            Nsteps_iter_relax (`int`, optional): 
                Number of steps per iteration for Langevin dynamics relaxation. Defaults to 500.

        Return:
            (`float`): 
                minimum distance between particles in the system after the relaxation

        Notes:
            - The thermostat is turned off by the end of the procedure. 
            - Make sure the system is initialized properly before calling this function.
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
        self.espresso_system.thermostat.turn_off()
        self.espresso_system.integrator.set_steepest_descent(f_max=0,
                                                        gamma=gamma, 
                                                        max_displacement=max_displacement)
        self.espresso_system.integrator.run(Nsteps_steepest_descent)
        logging.debug("*** Finished steepest descent minimization ***")
        logging.debug("*** Starting Langevin Dynamics relaxation ***")
        self.espresso_system.integrator.set_vv()
        self.espresso_system.thermostat.set_langevin(kT=1., gamma=gamma, seed=seed)
        self.espresso_system.integrator.run(Nsteps_iter_relax)
        self.espresso_system.thermostat.turn_off()
        logging.debug("*** Finished Langevin Dynamics relaxation ***")
        logging.info(f"*** Minimum particle distance after relaxation: {self.espresso_system.analysis.min_dist()} ***")
        logging.debug("*** Relaxation finished ***")
        return self.espresso_system.analysis.min_dist()
    
    def setup_electrostatic_interactions(self,units, kT, c_salt=None, solvent_permittivity=78.5, method='p3m', tune_p3m=True, accuracy=1e-3, params=None, verbose=False):
        """
        Sets up electrostatic interactions in an ESPResSo system.

        Args:
            units (`pint.UnitRegistry`): 
                Unit registry for handling physical units.


            kT (`pint.Quantity`): 
                Thermal energy.

            c_salt (`pint.Quantity`): 
                Added salt concentration. If provided, the program outputs the debye screening length. It is a mandatory parameter for the Debye-Hückel method.

            solvent_permittivity (`float`): 
                Solvent relative permittivity. Defaults to 78.5, correspoding to its value in water at 298.15 K.

            method (`str`): 
                Method for computing electrostatic interactions. Defaults to "p3m". 

            tune_p3m (`bool`): 
                If True, tunes P3M parameters for efficiency. Defaults to True. 

            accuracy (`float`): 
                Desired accuracy for electrostatics. Defaults to 1e-3.

            params (`dict`): 
                Additional parameters for the electrostatic method. For P3M, it can include 'mesh', 'alpha', 'cao' and `r_cut`. For Debye-Hückel, it can include 'r_cut'.

            verbose (`bool`): 
                If True, enables verbose output for P3M tuning. Defaults to False.

        Notes:
            - `c_salt` is a mandatory argument for setting up the Debye-Hückel electrostatic potential.
            - The calculated Bjerrum length is ouput to the log. If `c_salt` is provided, the calculated Debye screening length is also output to the log.
            - Currently, the only supported electrostatic methods are P3M ("p3m") and Debye-Hückel ("dh").
        """
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
                self.espresso_system.time_step=0.01
                if espressomd.version.friendly() == "4.2":
                    self.espresso_system.actors.add(coulomb)
                else:
                    self.espresso_system.electrostatics.solver = coulomb


                # save the optimal parameters and add them by hand

                p3m_params = coulomb.get_params()
                if espressomd.version.friendly() == "4.2":
                    self.espresso_system.actors.remove(coulomb)
                else:
                    self.espresso_system.electrostatics.solver = None
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
            self.espresso_system.actors.add(coulomb)
        else:
            self.espresso_system.electrostatics.solver = coulomb
        logging.debug("*** Electrostatics successfully added to the system ***")

    def setup_cpH (self, counter_ion, constant_pH, exclusion_range=None, use_exclusion_radius_per_type = False):
        """
        Sets up the Acid/Base reactions for acidic/basic particles defined in the pyMBE database
        to be sampled in the constant pH ensemble. 

        Args:
            counter_ion ('str'): 
                'name' of the counter_ion 'particle'.

            constant_pH ('float'): 
                pH-value.

            exclusion_range ('pint.Quantity', optional): 
                Below this value, no particles will be inserted.

            use_exclusion_radius_per_type ('bool', optional): 
                Controls if one exclusion_radius for each espresso_type is used. Defaults to 'False'.

        Returns:
            ('reaction_methods.ConstantpHEnsemble'): 
                Instance of a reaction_methods.ConstantpHEnsemble object from the espressomd library.
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.db.get_radius_map().values())*2.0
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.db.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        RE = reaction_methods.ConstantpHEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                exclusion_range=exclusion_range, 
                                                seed=self.seed, 
                                                constant_pH=constant_pH,
                                                exclusion_radius_per_type = exclusion_radius_per_type)
        conterion_tpl = self.db.get_template(name=counter_ion,
                                             pmb_type="particle")
        conterion_state = self.db.get_template(name=conterion_tpl.initial_state,
                                               pmb_type="particle_state")
        for reaction in self.db.get_reactions():
            if reaction.reaction_type not in ["monoprotic_acid", "monoprotic_base"]:
                continue
            default_charges = {}
            reactant_types  = []
            product_types   = []
            for participant in reaction.participants:
                state_tpl = self.db.get_template(name=participant.state_name,
                                                 pmb_type="particle_state")
                default_charges[state_tpl.es_type] = state_tpl.z
                if participant.coefficient < 0:
                    reactant_types.append(state_tpl.es_type)
                elif participant.coefficient > 0:
                    product_types.append(state_tpl.es_type)
            # Add counterion to the products
            if conterion_state.es_type not in product_types:
                product_types.append(conterion_state.es_type)
                default_charges[conterion_state.es_type] = conterion_state.z
                reaction.add_participant(particle_name=counter_ion,
                                         state_name=conterion_tpl.initial_state,
                                         coefficient=1)
            gamma=10**-reaction.pK
            RE.add_reaction(gamma=gamma,
                            reactant_types=reactant_types,
                            product_types=product_types,
                            default_charges=default_charges)
            reaction.add_simulation_method(simulation_method="cpH")
        return RE
    
    def setup_gcmc(self, c_salt_res, salt_cation_name, salt_anion_name, activity_coefficient, exclusion_range=None, use_exclusion_radius_per_type = False):
        """
        Sets up grand-canonical coupling to a reservoir of salt.
        For reactive systems coupled to a reservoir, the grand-reaction method has to be used instead.

        Args:
            c_salt_res ('pint.Quantity'): 
                Concentration of monovalent salt (e.g. NaCl) in the reservoir.

            salt_cation_name ('str'): 
                Name of the salt cation (e.g. Na+) particle.

            salt_anion_name ('str'): 
                Name of the salt anion (e.g. Cl-) particle.

            activity_coefficient ('callable'): 
                A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.

            exclusion_range('pint.Quantity', optional): 
                For distances shorter than this value, no particles will be inserted.

            use_exclusion_radius_per_type('bool',optional): 
                Controls if one exclusion_radius for each espresso_type is used. Defaults to 'False'.

        Returns:
            ('reaction_methods.ReactionEnsemble'): 
                Instance of a reaction_methods.ReactionEnsemble object from the espressomd library.
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.db.get_radius_map().values())*2.0
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.db.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                               exclusion_range=exclusion_range, 
                                               seed=self.seed, 
                                               exclusion_radius_per_type = exclusion_radius_per_type)
        # Determine the concentrations of the various species in the reservoir and the equilibrium constants
        determined_activity_coefficient = activity_coefficient(c_salt_res)
        K_salt = (c_salt_res.to('1/(N_A * reduced_length**3)')**2) * determined_activity_coefficient
        cation_tpl = self.db.get_template(pmb_type="particle",
                                          name=salt_cation_name)
        cation_state = self.db.get_template(pmb_type="particle_state",
                                            name=cation_tpl.initial_state)
        anion_tpl = self.db.get_template(pmb_type="particle",
                                          name=salt_anion_name)
        anion_state = self.db.get_template(pmb_type="particle_state",
                                            name=anion_tpl.initial_state)
        salt_cation_es_type = cation_state.es_type
        salt_anion_es_type = anion_state.es_type     
        salt_cation_charge = cation_state.z
        salt_anion_charge = anion_state.z
        if salt_cation_charge <= 0:
            raise ValueError('ERROR salt cation charge must be positive, charge ', salt_cation_charge)
        if salt_anion_charge >= 0:
            raise ValueError('ERROR salt anion charge must be negative, charge ', salt_anion_charge)
        # Grand-canonical coupling to the reservoir
        RE.add_reaction(gamma = K_salt.magnitude,
                        reactant_types = [],
                        reactant_coefficients = [],
                        product_types = [ salt_cation_es_type, salt_anion_es_type ],
                        product_coefficients = [ 1, 1 ],
                        default_charges = {salt_cation_es_type: salt_cation_charge, 
                                           salt_anion_es_type: salt_anion_charge})
        rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=salt_cation_name,
                                                            state_name=cation_state.name,
                                                            coefficient=1),
                                        ReactionParticipant(particle_name=salt_anion_name,
                                                            state_name=anion_state.name,
                                                            coefficient=1)],
                           pK=-np.log10(K_salt.magnitude),
                           reaction_type="ion_insertion",
                           simulation_method="GCMC")
        self.db._register_reaction(rx_tpl)
        return RE
    
    def setup_grxmc_reactions(self, pH_res, c_salt_res, proton_name, hydroxide_name, salt_cation_name, salt_anion_name, activity_coefficient, exclusion_range=None, use_exclusion_radius_per_type = False):
        """
        Sets up acid/base reactions for acidic/basic monoprotic particles defined in the pyMBE database, 
        as well as a grand-canonical coupling to a reservoir of small ions. 
        
        Args:
            pH_res ('float'): 
                pH-value in the reservoir.

            c_salt_res ('pint.Quantity'): 
                Concentration of monovalent salt (e.g. NaCl) in the reservoir.

            proton_name ('str'): 
                Name of the proton (H+) particle.

            hydroxide_name ('str'): 
                Name of the hydroxide (OH-) particle.

            salt_cation_name ('str'): 
                Name of the salt cation (e.g. Na+) particle.

            salt_anion_name ('str'): 
                Name of the salt anion (e.g. Cl-) particle.

            activity_coefficient ('callable'): 
                A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.

            exclusion_range('pint.Quantity', optional): 
                For distances shorter than this value, no particles will be inserted.

            use_exclusion_radius_per_type('bool', optional): 
                Controls if one exclusion_radius for each espresso_type is used. Defaults to 'False'.

        Returns:
            'tuple(reaction_methods.ReactionEnsemble,pint.Quantity)':

                'reaction_methods.ReactionEnsemble':  
                    espressomd reaction_methods object with all reactions necesary to run the GRxMC ensamble.
                
                'pint.Quantity': 
                    Ionic strength of the reservoir (useful for calculating partition coefficients).

        Notess:
            - This implementation uses the original formulation of the grand-reaction method by Landsgesell et al. [1].

        [1] Landsgesell, J., Hebbeker, P., Rud, O., Lunkad, R., Košovan, P., & Holm, C. (2020). Grand-reaction method for simulations of ionization equilibria coupled to ion partitioning. Macromolecules, 53(8), 3007-3020.
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.db.get_radius_map().values())*2.0
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.db.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                               exclusion_range=exclusion_range, 
                                               seed=self.seed, 
                                               exclusion_radius_per_type = exclusion_radius_per_type)
        # Determine the concentrations of the various species in the reservoir and the equilibrium constants
        cH_res, cOH_res, cNa_res, cCl_res = self.determine_reservoir_concentrations(pH_res, c_salt_res, activity_coefficient)
        ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
        determined_activity_coefficient = activity_coefficient(ionic_strength_res)
        K_W = cH_res.to('1/(N_A * reduced_length**3)') * cOH_res.to('1/(N_A * reduced_length**3)') * determined_activity_coefficient
        K_NACL = cNa_res.to('1/(N_A * reduced_length**3)') * cCl_res.to('1/(N_A * reduced_length**3)') * determined_activity_coefficient
        K_HCL = cH_res.to('1/(N_A * reduced_length**3)') * cCl_res.to('1/(N_A * reduced_length**3)') * determined_activity_coefficient
        cation_tpl = self.db.get_template(pmb_type="particle",
                                          name=salt_cation_name)
        cation_state = self.db.get_template(pmb_type="particle_state",
                                            name=cation_tpl.initial_state)
        anion_tpl = self.db.get_template(pmb_type="particle",
                                          name=salt_anion_name)
        anion_state = self.db.get_template(pmb_type="particle_state",
                                            name=anion_tpl.initial_state)
        proton_tpl = self.db.get_template(pmb_type="particle",
                                          name=proton_name)
        proton_state = self.db.get_template(pmb_type="particle_state",
                                            name=proton_tpl.initial_state)
        hydroxide_tpl = self.db.get_template(pmb_type="particle",
                                             name=hydroxide_name)
        hydroxide_state = self.db.get_template(pmb_type="particle_state",
                                               name=hydroxide_tpl.initial_state)
        proton_es_type = proton_state.es_type
        hydroxide_es_type = hydroxide_state.es_type
        salt_cation_es_type = cation_state.es_type
        salt_anion_es_type = anion_state.es_type
        proton_charge = proton_state.z
        hydroxide_charge = hydroxide_state.z          
        salt_cation_charge = cation_state.z
        salt_anion_charge = anion_state.z      
        if proton_charge <= 0:
            raise ValueError('ERROR proton charge must be positive, charge ', proton_charge)
        if salt_cation_charge <= 0:
            raise ValueError('ERROR salt cation charge must be positive, charge ', salt_cation_charge)
        if hydroxide_charge >= 0:
            raise ValueError('ERROR hydroxide charge must be negative, charge ', hydroxide_charge)
        if salt_anion_charge >= 0:
            raise ValueError('ERROR salt anion charge must be negative, charge ', salt_anion_charge)
        # Grand-canonical coupling to the reservoir
        # 0 = H+ + OH-
        RE.add_reaction(gamma = K_W.magnitude,
                        reactant_types = [],
                        reactant_coefficients = [],
                        product_types = [ proton_es_type, hydroxide_es_type ],
                        product_coefficients = [ 1, 1 ],
                        default_charges = {proton_es_type: proton_charge, 
                                           hydroxide_es_type: hydroxide_charge})
        rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=proton_name,
                                                            state_name=proton_state.name,
                                                            coefficient=1),
                                        ReactionParticipant(particle_name=hydroxide_name,
                                                            state_name=hydroxide_state.name,
                                                            coefficient=1)],
                           pK=-np.log10(K_W.magnitude),
                           reaction_type="ion_insertion",
                           simulation_method="GRxMC")
        self.db._register_reaction(rx_tpl)
        # 0 = Na+ + Cl-
        RE.add_reaction(gamma = K_NACL.magnitude,
                        reactant_types = [],
                        reactant_coefficients = [],
                        product_types = [ salt_cation_es_type, salt_anion_es_type ],
                        product_coefficients = [ 1, 1 ],
                        default_charges = {salt_cation_es_type: salt_cation_charge, 
                                        salt_anion_es_type: salt_anion_charge})
        rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=salt_cation_name,
                                                            state_name=cation_state.name,
                                                            coefficient=1),
                                        ReactionParticipant(particle_name=salt_anion_name,
                                                            state_name=anion_state.name,
                                                            coefficient=1)],
                           pK=-np.log10(K_NACL.magnitude),
                           reaction_type="ion_insertion",
                           simulation_method="GRxMC")
        self.db._register_reaction(rx_tpl)
        # 0 = Na+ + OH-
        RE.add_reaction(gamma = (K_NACL * K_W / K_HCL).magnitude,
                        reactant_types = [],
                        reactant_coefficients = [],
                        product_types = [ salt_cation_es_type, hydroxide_es_type ],
                        product_coefficients = [ 1, 1 ],
                        default_charges = {salt_cation_es_type: salt_cation_charge, 
                                           hydroxide_es_type: hydroxide_charge})
        rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=salt_cation_name,
                                                            state_name=cation_state.name,
                                                            coefficient=1),
                                        ReactionParticipant(particle_name=hydroxide_name,
                                                            state_name=hydroxide_state.name,
                                                            coefficient=1)],
                           pK=-np.log10((K_NACL * K_W / K_HCL).magnitude),
                           reaction_type="ion_insertion",
                           simulation_method="GRxMC")
        self.db._register_reaction(rx_tpl)
        # 0 = H+ + Cl-
        RE.add_reaction(gamma = K_HCL.magnitude,
                        reactant_types = [],
                        reactant_coefficients = [],
                        product_types = [ proton_es_type, salt_anion_es_type ],
                        product_coefficients = [ 1, 1 ],
                        default_charges = {proton_es_type: proton_charge, 
                                           salt_anion_es_type: salt_anion_charge})
        rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=proton_name,
                                                            state_name=proton_state.name,
                                                            coefficient=1),
                                        ReactionParticipant(particle_name=salt_anion_name,
                                                            state_name=anion_state.name,
                                                            coefficient=1)],
                           pK=-np.log10(K_HCL.magnitude),
                           reaction_type="ion_insertion",
                           simulation_method="GRxMC")
        self.db._register_reaction(rx_tpl)
        # Annealing moves to ensure sufficient sampling
        # Cation annealing H+ = Na+
        RE.add_reaction(gamma = (K_NACL / K_HCL).magnitude,
                        reactant_types = [proton_es_type],
                        reactant_coefficients = [ 1 ],
                        product_types = [ salt_cation_es_type ],
                        product_coefficients = [ 1 ],
                        default_charges = {proton_es_type: proton_charge, 
                                           salt_cation_es_type: salt_cation_charge})
        rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=proton_name,
                                                            state_name=proton_state.name,
                                                            coefficient=-1),
                                        ReactionParticipant(particle_name=salt_cation_name,
                                                            state_name=cation_state.name,
                                                            coefficient=1)],
                           pK=-np.log10((K_NACL / K_HCL).magnitude),
                           reaction_type="particle replacement",
                           simulation_method="GRxMC")
        self.db._register_reaction(rx_tpl)
        # Anion annealing OH- = Cl- 
        RE.add_reaction(gamma = (K_HCL / K_W).magnitude,
                        reactant_types = [hydroxide_es_type],
                        reactant_coefficients = [ 1 ],
                        product_types = [ salt_anion_es_type ],
                        product_coefficients = [ 1 ],
            default_charges = {hydroxide_es_type: hydroxide_charge, 
                               salt_anion_es_type: salt_anion_charge})
        rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=hydroxide_name,
                                                            state_name=hydroxide_state.name,
                                                            coefficient=-1),
                                        ReactionParticipant(particle_name=salt_anion_name,
                                                            state_name=anion_state.name,
                                                            coefficient=1)],
                           pK=-np.log10((K_HCL / K_W).magnitude),
                           reaction_type="particle replacement",
                           simulation_method="GRxMC")
        self.db._register_reaction(rx_tpl)
        for reaction in self.db.get_reactions():
            if reaction.reaction_type not in ["monoprotic_acid", "monoprotic_base"]:
                continue
            default_charges = {}
            reactant_types  = []
            product_types   = []
            for participant in reaction.participants:
                state_tpl = self.db.get_template(name=participant.state_name,
                                                 pmb_type="particle_state")
                default_charges[state_tpl.es_type] = state_tpl.z
                if participant.coefficient < 0:
                    reactant_types.append(state_tpl.es_type)
                    reactant_name=state_tpl.particle_name
                    reactant_state_name=state_tpl.name
                elif participant.coefficient > 0:
                    product_types.append(state_tpl.es_type)
                    product_name=state_tpl.particle_name
                    product_state_name=state_tpl.name

            Ka = (10**-reaction.pK * self.units.mol/self.units.l).to('1/(N_A * reduced_length**3)')
            # Reaction in terms of proton: HA = A + H+
            RE.add_reaction(gamma=Ka.magnitude,
                            reactant_types=reactant_types,
                            reactant_coefficients=[1],
                            product_types=product_types+[proton_es_type],
                            product_coefficients=[1, 1],
                            default_charges= default_charges | {proton_es_type: proton_charge})
            reaction.add_participant(particle_name=proton_name,
                                     state_name=proton_state.name,
                                     coefficient=1)
            reaction.add_simulation_method("GRxMC")
            # Reaction in terms of salt cation: HA = A + Na+
            RE.add_reaction(gamma=(Ka * K_NACL / K_HCL).magnitude,
                            reactant_types=reactant_types,
                            reactant_coefficients=[1],
                            product_types=product_types+[salt_cation_es_type],
                            product_coefficients=[1, 1],
                            default_charges=default_charges | {salt_cation_es_type: salt_cation_charge})
            rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=reactant_name,
                                                                state_name=reactant_state_name,
                                                                coefficient=-1),
                                            ReactionParticipant(particle_name=product_name,
                                                                state_name=product_state_name,
                                                                coefficient=1),
                                            ReactionParticipant(particle_name=salt_cation_name,
                                                                state_name=cation_state.name,
                                                                coefficient=1),],
                              pK=-np.log10((Ka * K_NACL / K_HCL).magnitude),
                              reaction_type=reaction.reaction_type+"_salt",
                              simulation_method="GRxMC")
            self.db._register_reaction(rx_tpl)
            # Reaction in terms of hydroxide: OH- + HA = A
            RE.add_reaction(gamma=(Ka / K_W).magnitude,
                            reactant_types=reactant_types+[hydroxide_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=product_types,
                            product_coefficients=[1],
                            default_charges=default_charges | {hydroxide_es_type: hydroxide_charge})
            rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=reactant_name,
                                                                state_name=reactant_state_name,
                                                                coefficient=-1),
                                            ReactionParticipant(particle_name=product_name,
                                                                state_name=product_state_name,
                                                                coefficient=1),
                                            ReactionParticipant(particle_name=hydroxide_name,
                                                                state_name=hydroxide_state.name,
                                                                coefficient=-1),],
                              pK=-np.log10((Ka / K_W).magnitude),
                              reaction_type=reaction.reaction_type+"_conjugate",
                              simulation_method="GRxMC")
            self.db._register_reaction(rx_tpl)
            # Reaction in terms of salt anion: Cl- + HA = A
            RE.add_reaction(gamma=(Ka / K_HCL).magnitude,
                            reactant_types=reactant_types+[salt_anion_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=product_types,
                            product_coefficients=[1],
                            default_charges=default_charges | {salt_anion_es_type: salt_anion_charge})
            rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=reactant_name,
                                                                state_name=reactant_state_name,
                                                                coefficient=-1),
                                            ReactionParticipant(particle_name=product_name,
                                                                state_name=product_state_name,
                                                                coefficient=1),
                                            ReactionParticipant(particle_name=salt_anion_name,
                                                                state_name=anion_state.name,
                                                                coefficient=-1),],
                              pK=-np.log10((Ka / K_HCL).magnitude),
                              reaction_type=reaction.reaction_type+"_salt",
                              simulation_method="GRxMC")
            self.db._register_reaction(rx_tpl)
        return RE, ionic_strength_res
    
    def setup_grxmc_unified(self, pH_res, c_salt_res, cation_name, anion_name, activity_coefficient, exclusion_range=None, use_exclusion_radius_per_type = False):
        """
        Sets up acid/base reactions for acidic/basic 'particles' defined in the pyMBE database, as well as a grand-canonical coupling to a 
        reservoir of small ions using a unified formulation for small ions.

        Args:
            pH_res ('float'): 
                pH-value in the reservoir.

            c_salt_res ('pint.Quantity'): 
                Concentration of monovalent salt (e.g. NaCl) in the reservoir.

            cation_name ('str'): 
                Name of the cationic particle.

            anion_name ('str'): 
                Name of the anionic particle.

            activity_coefficient ('callable'): 
                A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.

            exclusion_range('pint.Quantity', optional): 
                Below this value, no particles will be inserted.
            
            use_exclusion_radius_per_type('bool', optional): 
                Controls if one exclusion_radius per each espresso_type. Defaults to 'False'.

        Returns:
            'tuple(reaction_methods.ReactionEnsemble,pint.Quantity)':

                'reaction_methods.ReactionEnsemble':  
                    espressomd reaction_methods object with all reactions necesary to run the GRxMC ensamble.
                
                'pint.Quantity': 
                    Ionic strength of the reservoir (useful for calculating partition coefficients).

        Notes:
            - This implementation uses the formulation of the grand-reaction method by Curk et al. [1], which relies on "unified" ion types X+ = {H+, Na+} and X- = {OH-, Cl-}. 
            - A function that implements the original version of the grand-reaction method by Landsgesell et al. [2] is also available under the name 'setup_grxmc_reactions'.

        [1] Curk, T., Yuan, J., & Luijten, E. (2022). Accelerated simulation method for charge regulation effects. The Journal of Chemical Physics, 156(4).
        [2] Landsgesell, J., Hebbeker, P., Rud, O., Lunkad, R., Košovan, P., & Holm, C. (2020). Grand-reaction method for simulations of ionization equilibria coupled to ion partitioning. Macromolecules, 53(8), 3007-3020.
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.db.get_radius_map().values())*2.0
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.db.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                               exclusion_range=exclusion_range, 
                                               seed=self.seed, 
                                               exclusion_radius_per_type = exclusion_radius_per_type)
        # Determine the concentrations of the various species in the reservoir and the equilibrium constants
        cH_res, cOH_res, cNa_res, cCl_res = self.determine_reservoir_concentrations(pH_res, c_salt_res, activity_coefficient)
        ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
        determined_activity_coefficient = activity_coefficient(ionic_strength_res)
        a_hydrogen = (10 ** (-pH_res) * self.units.mol/self.units.l).to('1/(N_A * reduced_length**3)')
        a_cation = (cH_res+cNa_res).to('1/(N_A * reduced_length**3)') * np.sqrt(determined_activity_coefficient)
        a_anion = (cH_res+cNa_res).to('1/(N_A * reduced_length**3)') * np.sqrt(determined_activity_coefficient)
        K_XX = a_cation * a_anion
        cation_tpl = self.db.get_template(pmb_type="particle",
                                          name=cation_name)
        cation_state = self.db.get_template(pmb_type="particle_state",
                                            name=cation_tpl.initial_state)
        anion_tpl = self.db.get_template(pmb_type="particle",
                                          name=anion_name)
        anion_state = self.db.get_template(pmb_type="particle_state",
                                            name=anion_tpl.initial_state)
        cation_es_type = cation_state.es_type
        anion_es_type = anion_state.es_type     
        cation_charge = cation_state.z
        anion_charge = anion_state.z
        if cation_charge <= 0:
            raise ValueError('ERROR cation charge must be positive, charge ', cation_charge)
        if anion_charge >= 0:
            raise ValueError('ERROR anion charge must be negative, charge ', anion_charge)
        # Coupling to the reservoir: 0 = X+ + X-
        RE.add_reaction(gamma = K_XX.magnitude,
                        reactant_types = [],
                        reactant_coefficients = [],
                        product_types = [ cation_es_type, anion_es_type ],
                        product_coefficients = [ 1, 1 ],
                        default_charges = {cation_es_type: cation_charge, 
                                           anion_es_type: anion_charge})
        rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=cation_name,
                                                            state_name=cation_state.name,
                                                            coefficient=1),
                                        ReactionParticipant(particle_name=anion_name,
                                                            state_name=anion_state.name,
                                                            coefficient=1)],
                           pK=-np.log10(K_XX.magnitude),
                           reaction_type="ion_insertion",
                           simulation_method="GCMC")
        self.db._register_reaction(rx_tpl)
        for reaction in self.db.get_reactions():
            if reaction.reaction_type not in ["monoprotic_acid", "monoprotic_base"]:
                continue
            default_charges = {}
            reactant_types  = []
            product_types   = []
            for participant in reaction.participants:
                state_tpl = self.db.get_template(name=participant.state_name,
                                                 pmb_type="particle_state")
                default_charges[state_tpl.es_type] = state_tpl.z
                if participant.coefficient < 0:
                    reactant_types.append(state_tpl.es_type)
                    reactant_name=state_tpl.particle_name
                    reactant_state_name=state_tpl.name
                elif participant.coefficient > 0:
                    product_types.append(state_tpl.es_type)
                    product_name=state_tpl.particle_name
                    product_state_name=state_tpl.name

            Ka = (10**-reaction.pK * self.units.mol/self.units.l).to('1/(N_A * reduced_length**3)')
            gamma_K_AX = Ka.to('1/(N_A * reduced_length**3)').magnitude * a_cation / a_hydrogen
            # Reaction in terms of small cation: HA = A + X+
            RE.add_reaction(gamma=gamma_K_AX.magnitude,
                            reactant_types=reactant_types,
                            reactant_coefficients=[1],
                            product_types=product_types+[cation_es_type],
                            product_coefficients=[1, 1],
                            default_charges=default_charges|{cation_es_type: cation_charge})
            reaction.add_participant(particle_name=cation_name,
                                     state_name=cation_state.name,
                                     coefficient=1)
            reaction.add_simulation_method("GRxMC")
            # Reaction in terms of small anion: X- + HA = A
            RE.add_reaction(gamma=gamma_K_AX.magnitude / K_XX.magnitude,
                            reactant_types=reactant_types+[anion_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=product_types,
                            product_coefficients=[1],
                            default_charges=default_charges|{anion_es_type: anion_charge})
            rx_tpl = Reaction(participants=[ReactionParticipant(particle_name=reactant_name,
                                                                state_name=reactant_state_name,
                                                                coefficient=-1),
                                            ReactionParticipant(particle_name=product_name,
                                                                state_name=product_state_name,
                                                                coefficient=1),
                                            ReactionParticipant(particle_name=anion_name,
                                                                state_name=anion_state.name,
                                                                coefficient=-1),],
                              pK=-np.log10(gamma_K_AX.magnitude / K_XX.magnitude),
                              reaction_type=reaction.reaction_type+"_conjugate",
                              simulation_method="GRxMC")
            self.db._register_reaction(rx_tpl)
        return RE, ionic_strength_res
    
    def setup_langevin_dynamics(self, kT, seed,time_step=1e-2, gamma=1, tune_skin=True, min_skin=1, max_skin=None, tolerance=1e-3, int_steps=200, adjust_max_skin=True):
        """
        Sets up Langevin Dynamics for an ESPResSo simulation system.

        Args:

            kT (`pint.Quantity`): 
                Target temperature in reduced energy units.

            seed (`int`): 
                Seed for the random number generator for the thermostat.

            time_step (`float`, optional): 
                Integration time step. Defaults to 1e-2.

            gamma (`float`, optional): 
                Damping coefficient for the Langevin thermostat. Defaults to 1.

            tune_skin (`bool`, optional): 
                Whether to optimize the skin parameter. Defaults to True.

            min_skin (`float`, optional): 
                Minimum skin value for optimization. Defaults to 1.

            max_skin (`float`, optional): 
                Maximum skin value for optimization. Defaults to None, which is handled by setting its value to box length / 2.

            tolerance (`float`, optional): 
                Tolerance for skin optimization. Defaults to 1e-3.

            int_steps (`int`, optional): 
                Number of integration steps for tuning. Defaults to 200.

            adjust_max_skin (`bool`, optional): 
                Whether to adjust the maximum skin value during tuning. Defaults to True.
        """        
        if not isinstance(seed, int):
            raise TypeError("seed must be an integer.")
        if not isinstance(time_step, (float, int)) or time_step <= 0:
            raise ValueError("time_step must be a positive number.")
        if not isinstance(gamma, (float, int)) or gamma <= 0:
            raise ValueError("gamma must be a positive number.")
        if max_skin is None:
            max_skin=self.espresso_system.box_l[0]/2
        if min_skin >= max_skin:
            raise ValueError("min_skin must be smaller than max_skin.")
        self.espresso_system.time_step=time_step
        self.espresso_system.integrator.set_vv()
        self.espresso_system.thermostat.set_langevin(kT= kT.to('reduced_energy').magnitude, 
                                                gamma= gamma, 
                                                seed= seed)
        # Optimize the value of skin
        if tune_skin:
            logging.debug("*** Optimizing skin ... ***")
            self.espresso_system.cell_system.tune_skin(min_skin=min_skin, 
                                                max_skin=max_skin, 
                                                tol=tolerance, 
                                                int_steps=int_steps, 
                                                adjust_max_skin=adjust_max_skin)
            logging.info(f"Optimized skin value: {self.espresso_system.cell_system.skin}")

    def setup_lj_interactions(self, shift_potential=True, combining_rule='Lorentz-Berthelot'):
        """
        Sets up the Lennard-Jones (LJ) potential between all pairs of particle states defined in the pyMBE database.

        Args:
            shift_potential('bool', optional): 
                If True, a shift will be automatically computed such that the potential is continuous at the cutoff radius. Otherwise, no shift will be applied. Defaults to True.

            combining_rule('string', optional): 
                combining rule used to calculate 'sigma' and 'epsilon' for the potential between a pair of particles. Defaults to 'Lorentz-Berthelot'.

            warning('bool', optional): 
                switch to activate/deactivate warning messages. Defaults to True.

        Notes:
            - Currently, the only 'combining_rule' supported is Lorentz-Berthelot.
            - Check the documentation of ESPResSo for more info about the potential https://espressomd.github.io/doc4.2.0/inter_non-bonded.html

        """
        from itertools import combinations_with_replacement
        particle_templates = self.db.get_templates("particle")
        shift = "auto" if shift_potential else 0
        if shift == "auto":
            shift_tpl = shift
        else:
            shift_tpl = PintQuantity.from_quantity(q=shift*self.units.reduced_length,
                                                   expected_dimension="length",
                                                   ureg=self.units)
        # Get all particle states registered in pyMBE
        state_entries = []
        for tpl in particle_templates.values():
            for state in self.db.get_particle_states_templates(particle_name=tpl.name).values():
                state_entries.append((tpl, state))

        # Iterate over all unique state pairs
        for (tpl1, state1), (tpl2, state2) in combinations_with_replacement(state_entries, 2):

            lj_parameters = self.db.get_lj_parameters(particle_name1=tpl1.name,
                                                   particle_name2=tpl2.name,
                                                   combining_rule=combining_rule)
            if not lj_parameters:
                continue

            self.espresso_system.non_bonded_inter[state1.es_type, state2.es_type].lennard_jones.set_params(
                epsilon=lj_parameters["epsilon"].to("reduced_energy").magnitude,
                sigma=lj_parameters["sigma"].to("reduced_length").magnitude,
                cutoff=lj_parameters["cutoff"].to("reduced_length").magnitude,
                offset=lj_parameters["offset"].to("reduced_length").magnitude,
                shift=shift)
                
            lj_template = LJInteractionTemplate(state1=state1.name,
                                                state2=state2.name,
                                                sigma=PintQuantity.from_quantity(q=lj_parameters["sigma"],
                                                                                 expected_dimension="length",
                                                                                 ureg=self.units),
                                                epsilon=PintQuantity.from_quantity(q=lj_parameters["epsilon"],
                                                                                   expected_dimension="energy",
                                                                                   ureg=self.units),
                                                cutoff=PintQuantity.from_quantity(q=lj_parameters["cutoff"],
                                                                                  expected_dimension="length",
                                                                                  ureg=self.units),
                                                offset=PintQuantity.from_quantity(q=lj_parameters["offset"],
                                                                                  expected_dimension="length",
                                                                                  ureg=self.units),
                                                shift=shift_tpl)
            self.db._register_template(lj_template)
    
    def update_particle_id(self, old_pid, new_pid):
        """
        Method to be called if particles have been previously added to a simulation engine without using pyMBE.

        Args:
            old_pid (int): old particle id to be replaced.
            new_pid (int): new particle id to be assigned to instance in the pyMBE database. 
        """
        
        if self.espresso_system == None:
            raise ValueError('No simulation engine has been set to pymbe')
        
        self.db._update_instance(instance_id=old_pid,
                                 pmb_type='particle',
                                 attribute='particle_id',
                                 value=new_pid)
    
        particle_id1_instances_ids=self.db._find_instance_ids_by_attribute(pmb_type='bond', 
                                                                           attribute='particle_id1', 
                                                                           value=old_pid)
        for bond_id in particle_id1_instances_ids:
            self.db._update_instance(instance_id=bond_id,
                                     pmb_type='bond',
                                     attribute='particle_id1',
                                     value=new_pid)
        particle_id2_instances_ids=self.db._find_instance_ids_by_attribute(pmb_type='bond', 
                                                                attribute='particle_id2', 
                                                                value=old_pid)
        for bond_id in particle_id2_instances_ids:
            self.db._update_instance(instance_id=bond_id,
                                     pmb_type='bond',
                                     attribute='particle_id2',
                                     value=new_pid)

    
    def add_instances_to_engine(self):
        """
            This method adds the set of particles instances and bond instances and angle instances
            that are not present in the pymbe data base 
        """
        
        missing_particle_ids=self.db._find_instance_ids_by_attribute(pmb_type='particle', 
                                                                     attribute='added_to_engine', 
                                                                     value=False)
        if not missing_particle_ids:
            raise RuntimeError('No particles instances in the pyMBE database set to be added to the simulation engine')
        
        added_particle_ids=self.db._find_instance_ids_by_attribute(pmb_type='particle',
                                                                   attribute='added_to_engine', 
                                                                   value=True)
        ids_in_espresso = list(self._get_particle_ids_in_espresso())
        
        overlapping_ids = set(ids_in_espresso).intersection(set(missing_particle_ids))
        for overlapping_id in overlapping_ids:
            missing_particle_ids.remove(overlapping_id)
            new_id = max(ids_in_espresso+added_particle_ids+missing_particle_ids)+1
            self.update_particle_id(old_pid=overlapping_id, 
                                    new_pid=new_id)
            missing_particle_ids.append(new_id)

        if overlapping_ids:
            warnings.warn("""Please review your setup, you have previously added a set of particles to ESPResSo. 
                          The particle ids of the pyMBE database have been updated taking into account the id of 
                          the last particle from ESPResSo. The following particle ids were updated: {}. """.format(overlapping_ids))
        
        for id in missing_particle_ids:
            self._add_particle(id)
        
        missing_bond_ids=self.db._find_instance_ids_by_attribute(pmb_type='bond', 
                                                                 attribute='added_to_engine', 
                                                                 value=False)
        
        if len(missing_bond_ids)>0:
            for id in missing_bond_ids:
                bond_instance=self.db.get_instance(pmb_type='bond',
                                                   instance_id=id)
                particle_id1=bond_instance.particle_id1
                particle_id2=bond_instance.particle_id2
                self._add_bond(particle_id1,
                               particle_id2,
                               bond_instance)
                
        missing_angle_ids=self.db._find_instance_ids_by_attribute(pmb_type='angle', 
                                                                 attribute='added_to_engine', 
                                                                 value=False)
        if len(missing_angle_ids)>0:
            for id in missing_angle_ids:
                angle_instance=self.db.get_instance(pmb_type='angle',
                                                   instance_id=id)
                particle_id1=angle_instance.particle_id1
                particle_id2=angle_instance.particle_id2
                particle_id3=angle_instance.particle_id3
                self._add_angle(particle_id1,
                               particle_id2,
                               particle_id3,
                               angle_instance)
    
