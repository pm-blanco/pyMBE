#
# Copyright (C) 2023-2026 pyMBE-dev team
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
import re
import json
import pint
import numpy as np
import pandas as pd
import scipy.constants
import scipy.optimize
import logging
import importlib.resources

# Database
from pyMBE.storage.manager import Manager
from pyMBE.storage.pint_quantity import PintQuantity
## Templates
from pyMBE.storage.templates.particle import ParticleTemplate, ParticleStateTemplate
from pyMBE.storage.templates.residue import ResidueTemplate
from pyMBE.storage.templates.molecule import MoleculeTemplate
from pyMBE.storage.templates.peptide import PeptideTemplate
from pyMBE.storage.templates.protein import ProteinTemplate
from pyMBE.storage.templates.hydrogel import HydrogelTemplate, HydrogelNode, HydrogelChain
from pyMBE.storage.templates.bond import BondTemplate
from pyMBE.storage.templates.angle import AngleTemplate
## Instances
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.instances.residue import ResidueInstance
from pyMBE.storage.instances.molecule import MoleculeInstance
from pyMBE.storage.instances.peptide import PeptideInstance
from pyMBE.storage.instances.protein import ProteinInstance
from pyMBE.storage.instances.bond import BondInstance
from pyMBE.storage.instances.angle import AngleInstance
from pyMBE.storage.instances.hydrogel import HydrogelInstance

from pyMBE.simulation_builder.espresso_engine import EspressoSimulation
from pyMBE.simulation_builder.lammps_engine import LammpsSimulation 
from pyMBE.simulation_builder.base_engine import DummyEngine
from pyMBE.simulation_builder.engine_protocol import LammpsProtocol
## Reactions
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
# Utilities
import pyMBE.lib.handy_functions as hf
import pyMBE.storage.io as io

class pymbe_library():
    """
    Core library of the Molecular Builder for ESPResSo (pyMBE).

    Attributes:
        N_A ('pint.Quantity'):
            Avogadro number.

        kB ('pint.Quantity'):
            Boltzmann constant.

        e ('pint.Quantity'):
            Elementary charge.

        kT ('pint.Quantity'):
            Thermal energy corresponding to the set temperature.

        Kw ('pint.Quantity'):
            Ionic product of water, used in G-RxMC and Donnan-related calculations.

        db ('Manager'):
            Database manager holding all pyMBE templates, instances and reactions.

        rng ('numpy.random.Generator'):
            Random number generator initialized with the provided seed.

        units ('pint.UnitRegistry'):
            Pint unit registry used for unit-aware calculations.

        lattice_builder ('pyMBE.lib.lattice.LatticeBuilder'):
            Optional lattice builder object (initialized as ''None'').
            
        root ('importlib.resources.abc.Traversable'):
            Root path to the pyMBE package resources.
    """

    def __init__(self, seed, temperature=None, unit_length=None, unit_charge=None, Kw=None):
        """
        Initializes the pyMBE library.

        Args:
            seed ('int'):
                Seed for the random number generator.

            temperature ('pint.Quantity', optional):
                Simulation temperature. If ''None'', defaults to 298.15 K.

            unit_length ('pint.Quantity', optional):
                Reference length for reduced units. If ''None'', defaults to
                0.355 nm.

            unit_charge ('pint.Quantity', optional):
                Reference charge for reduced units. If ''None'', defaults to
                one elementary charge.

            Kw ('pint.Quantity', optional):
                Ionic product of water (typically in mol²/L²). If ''None'',
                defaults to 1e-14 mol²/L².
        """
        # Seed and RNG
        self.seed=seed
        self.rng = np.random.default_rng(seed)
        self.units=pint.UnitRegistry()
        self.N_A=scipy.constants.N_A / self.units.mol
        self.kB=scipy.constants.k * self.units.J / self.units.K
        self.e=scipy.constants.e * self.units.C
        self.set_reduced_units(unit_length=unit_length, 
                               unit_charge=unit_charge,
                               temperature=temperature, 
                               Kw=Kw)
        
        self.db = Manager(units=self.units)
        self.simulation_engine = DummyEngine()
        self.lattice_builder = None
        self.root = importlib.resources.files(__package__)

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
            
    def _check_dimensionality(self, variable, expected_dimensionality):
        """
        Checks if the dimensionality of 'variable' matches 'expected_dimensionality'.

        Args:
            variable ('pint.Quantity'): 
                Quantity to be checked.

            expected_dimensionality ('str'): 
                Expected dimension of the variable.

        Returns:
            ('bool'): 
                'True' if the variable if of the expected dimensionality, 'False' otherwise.

        Notes:
            - 'expected_dimensionality' takes dimensionality following the Pint standards [docs](https://pint.readthedocs.io/en/0.10.1/wrapping.html?highlight=dimensionality#checking-dimensionality).
            - For example, to check for a variable corresponding to a velocity 'expected_dimensionality = "[length]/[time]"'    
        """
        correct_dimensionality=variable.check(f"{expected_dimensionality}")      
        if not correct_dimensionality:
            raise ValueError(f"The variable {variable} should have a dimensionality of {expected_dimensionality}, instead the variable has a dimensionality of {variable.dimensionality}")
        return correct_dimensionality   
            
    def _check_pka_set(self, pka_set):
        """
        Checks that 'pka_set' has the formatting expected by pyMBE.
       
        Args:
            pka_set ('dict'): 
                {"name" : {"pka_value": pka, "acidity": acidity}}
        """
        required_keys=['pka_value','acidity']
        for required_key in required_keys:
            for pka_name, pka_entry in pka_set.items():
                if required_key not in pka_entry:
                    raise ValueError(f'missing a required key "{required_key}" in entry "{pka_name}" of pka_set ("{pka_entry}")')

    def _create_hydrogel_chain(self, hydrogel_chain, nodes,box_l, use_default_bond=False, gen_angle=False):
        """
        Creates a chain between two nodes of a hydrogel.

        Args:
            hydrogel_chain ('HydrogelChain'): 
                template of a hydrogel chain
            nodes ('dict'): 
                {node_index: {"name": node_particle_name, "pos": node_position, "id": node_particle_instance_id}}
            box_l('list[float,float,float]'): side length of the simulation box for x,y and z coordinates.
            use_default_bond ('bool', optional): 
                If True, use a default bond template if no specific template exists. Defaults to False.

            gen_angle ('bool', optional):
                If True, generate the angle potentials internal to the created
                chain molecule. Junction angles near the hydrogel crosslinkers
                are handled separately at the hydrogel level.

        Return:
            ('int'): 
                molecule_id of the created hydrogel chian.

        Notes:
            - If the chain is defined between node_start = ''[0 0 0]'' and node_end = ''[1 1 1]'', the chain will be placed between these two nodes.
            - The chain will be placed in the direction of the vector between 'node_start' and 'node_end'. 
        """
        if self.lattice_builder is None:
            raise ValueError("LatticeBuilder is not initialized. Use 'initialize_lattice_builder' first.")
        molecule_tpl = self.db.get_template(pmb_type="molecule",
                                            name=hydrogel_chain.molecule_name)
        residue_list = molecule_tpl.residue_list
        molecule_name = molecule_tpl.name
        node_start = hydrogel_chain.node_start
        node_end = hydrogel_chain.node_end
        node_start_label = self.lattice_builder._create_node_label(node_start)
        node_end_label = self.lattice_builder._create_node_label(node_end)
        _, reverse = self.lattice_builder._get_node_vector_pair(node_start, node_end)
        if node_start == node_end and residue_list != residue_list[::-1]:
            raise ValueError(f"Aborted creation of hydrogel chain between '{node_start}' and '{node_end}' because pyMBE could not resolve a unique topology for that chain")
        if reverse:
            reverse_residue_order=True
        else:
            reverse_residue_order=False
        start_node_id = nodes[node_start_label]["id"]
        end_node_id = nodes[node_end_label]["id"]
        # Finding a backbone vector between node_start and node_end
        vec_between_nodes = np.array(nodes[node_end_label]["pos"]) - np.array(nodes[node_start_label]["pos"])
        vec_between_nodes = vec_between_nodes - self.lattice_builder.box_l * np.round(vec_between_nodes/self.lattice_builder.box_l)
        backbone_vector = vec_between_nodes / (self.lattice_builder.mpc+1)
        if reverse_residue_order:
            vec_between_nodes *= -1.0
        # Calculate the start position of the chain
        chain_residues = self.db.get_template(pmb_type="molecule",
                                              name=molecule_name).residue_list
        part_start_chain_name = self.db.get_template(pmb_type="residue",
                                                     name=chain_residues[0]).central_bead
        lj_parameters = self.get_lj_parameters(particle_name1=nodes[node_start_label]["name"],
                                               particle_name2=part_start_chain_name)
        bond_tpl = self.get_bond_template(particle_name1=nodes[node_start_label]["name"],
                                          particle_name2=part_start_chain_name,
                                          use_default_bond=use_default_bond)
        l0 = hf.calculate_initial_bond_length(lj_parameters=lj_parameters,
                                              bond_type=bond_tpl.bond_type,
                                              bond_parameters=bond_tpl.get_parameters(ureg=self.units))
        first_bead_pos = np.array((nodes[node_start_label]["pos"])) + np.array(backbone_vector)*l0
        mol_id = self.create_molecule(name=molecule_name,  # Use the name defined earlier
                                      number_of_molecules=1,  # Creating one chain
                                      box_l=box_l, ### Add lattice_builder box length size, this should be box_l=[self.lattice_builder.box_l]*3
                                      list_of_first_residue_positions=[first_bead_pos.tolist()], #Start at the first node
                                      backbone_vector=np.array(backbone_vector)/l0,
                                      use_default_bond=use_default_bond,
                                      reverse_residue_order=reverse_residue_order,
                                      gen_angle=gen_angle)[0]
        chain_pids = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                             attribute="molecule_id",
                                                             value=mol_id)
        self.create_bond(particle_id1=start_node_id,particle_id2=chain_pids[0],use_default_bond=use_default_bond)
        self.create_bond(particle_id1=chain_pids[-1],particle_id2=end_node_id,use_default_bond=use_default_bond)
        return mol_id

    def _generate_hydrogel_crosslinker_angles(self,  central_particle_ids):
        """
        Generate hydrogel angles centered on crosslinkers and adjacent terminal beads.

        If the user defines any explicit angle template for such junction
        triplets, then all required junction triplets must be defined. If none
        are defined, hydrogel construction proceeds without crosslinker-adjacent
        angles.
        """
        particle_instances = self.db.get_instances(pmb_type="particle")
        bonded_neighbors = {}
        for bond in self.db.get_instances(pmb_type="bond").values():
            bonded_neighbors.setdefault(bond.particle_id1, set()).add(bond.particle_id2)
            bonded_neighbors.setdefault(bond.particle_id2, set()).add(bond.particle_id1)

        triplets = []
        for central_particle_id in sorted(set(central_particle_ids)):
            neighbors = sorted(bonded_neighbors.get(central_particle_id, set()))
            central_name = particle_instances[central_particle_id].name
            for idx_i in range(len(neighbors)):
                for idx_k in range(idx_i + 1, len(neighbors)):
                    side_particle_id1 = neighbors[idx_i]
                    side_particle_id3 = neighbors[idx_k]
                    side_name1 = particle_instances[side_particle_id1].name
                    side_name3 = particle_instances[side_particle_id3].name
                    angle_key = AngleTemplate.make_angle_key(side1=side_name1,
                                                             central=central_name,
                                                             side2=side_name3)
                    triplets.append((side_particle_id1,
                                     central_particle_id,
                                     side_particle_id3,
                                     angle_key))

        defined_angle_templates = self.db.get_templates(pmb_type="angle")
        defined_angle_keys = {angle_key for _, _, _, angle_key in triplets if angle_key in defined_angle_templates}
        if not defined_angle_keys:
            logging.warning("No angle templates defined for hydrogel crosslinkers")
            return
        missing_angle_keys = sorted({angle_key for _, _, _, angle_key in triplets if angle_key not in defined_angle_keys})
        if missing_angle_keys:
            raise ValueError("Hydrogel crosslinker-adjacent angle templates must be defined for all required triplets. "
                             f"Missing definitions for: {missing_angle_keys}")
        for side_particle_id1, central_particle_id, side_particle_id3, _ in triplets:
            self.create_angular_potential(particle_id1=side_particle_id1,
                                          particle_id2=central_particle_id,
                                          particle_id3=side_particle_id3,
                                          use_default_angle=False)

    def _create_hydrogel_node(self, node_index, node_name,box_l):
        """
        Set a node residue type.
        
        Args:
            node_index ('str'): 
                Lattice node index in the form of a string, e.g. "[0 0 0]".

            node_name ('str'): 
                name of the node particle defined in pyMBE.
            
            box_l('list[float,float,float]'): list of floats with the dimensions of the box

        Returns:
            ('tuple(list,int)'):
                ('list'): Position of the node in the lattice.
                ('int'): Particle ID of the node.
        """
        if self.lattice_builder is None:
            raise ValueError("LatticeBuilder is not initialized. Use 'initialize_lattice_builder' first.")
        node_position = np.array(node_index)*0.25*self.lattice_builder.box_l
        p_id = self.create_particle(name = node_name,
                                    box_l=box_l,
                                    number_of_particles=1,
                                    position = [node_position])
        key = self.lattice_builder._get_node_by_label(f"[{node_index[0]} {node_index[1]} {node_index[2]}]")
        self.lattice_builder.nodes[key] = node_name
        return node_position.tolist(), p_id[0]

    def _get_residue_list_from_sequence(self, sequence):
        """
        Convenience function to get a 'residue_list' from a protein or peptide 'sequence'.

        Args:
            sequence ('lst'): 
                 Sequence of the peptide or protein.

        Returns:
            residue_list ('list' of 'str'): 
                List of the 'name's of the 'residue's  in the sequence of the 'molecule'.             
        """
        residue_list = []
        for item in sequence:
            residue_name='AA-'+item
            residue_list.append(residue_name)
        return residue_list
    
    def _get_template_type(self, name, allowed_types):
        """
        Validate that a template name resolves unambiguously to exactly one
        allowed pmb_type in the pyMBE database and return it.

        Args:
            name ('str'): 
                Name of the template to validate.

            allowed_types ('set[str]'):  
                Set of allowed pmb_type values (e.g. {"molecule", "peptide"}).

        Returns:
            ('str'): 
                Resolved pmb_type.

        Notes:
            - This method does *not* return the template itself, only the validated pmb_type. 
        """
        registered_pmb_types_with_name = self.db._find_template_types(name=name)
        filtered_types = allowed_types.intersection(registered_pmb_types_with_name)
        if len(filtered_types) > 1:
            raise ValueError(f"Ambiguous template name '{name}': found {len(filtered_types)} templates in the pyMBE database. Molecule creation aborted.")  
        if len(filtered_types) == 0:
            raise ValueError(f"No {allowed_types} template found with name '{name}'. Found templates of types: {filtered_types}.")
        return next(iter(filtered_types))

    def _delete_particles_from_engine(self, particle_ids):
        """
        Remove a list of particles from an ESPResSo simulation system.

        Args:
            particle_ids  ('Iterable[int]'):
                A list (or other iterable) of ESPResSo particle IDs to remove.

        Notes:
            - This method removes particles only from the ESPResSo simulation,
            **not** from the pyMBE database. Database cleanup must be handled
            separately by the caller.
            - Attempting to remove a non-existent particle ID will raise
            an ESPResSo error.
        """
        self.simulation_engine._delete_particles(particle_ids)

    def add_instances_to_engine(self):
        self.simulation_engine.add_instances_to_engine()

    def calculate_center_of_mass(self, instance_id, pmb_type):
        """
        Calculates the center of mass of a pyMBE object instance in an ESPResSo system.

        Args:
            instance_id ('int'):
                pyMBE instance ID of the object whose center of mass is calculated.

            pmb_type ('str'):
                Type of the pyMBE object. Must correspond to a particle-aggregating
                template type (e.g. '"molecule"', '"residue"', '"peptide"', '"protein"').

        Returns:
            ('numpy.ndarray'):
                Array of shape '(3,)' containing the Cartesian coordinates of the
                center of mass.

        Notes:
            - This method assumes equal mass for all particles.
            - Periodic boundary conditions are *not* unfolded; positions are taken
            directly from ESPResSo particle coordinates.
        """
        return self.simulation_engine.calculate_center_of_mass(instance_id=instance_id,
                                                               pmb_type=pmb_type)
    
    def calculate_HH(self, template_name, pH_list=None, pka_set=None):
        """
        Calculates the charge in the template object according to the ideal  Henderson–Hasselbalch titration curve.

        Args:
            template_name ('str'):
                Name of the template.

            pH_list ('list[float]', optional):
                pH values at which the charge is evaluated.
                Defaults to 50 values between 2 and 12.

            pka_set ('dict', optional):
                Mapping: {particle_name: {"pka_value": 'float', "acidity": "acidic"|"basic"}}

        Returns:
            'list[float]':
                Net molecular charge at each pH value.
        """
        if pH_list is None:
            pH_list = np.linspace(2, 12, 50)
        if pka_set is None:
            pka_set = self.get_pka_set()
        self._check_pka_set(pka_set=pka_set)
        particle_counts = self.db.get_particle_templates_under(template_name=template_name,
                                                               return_counts=True)
        if not particle_counts:
            return [None] * len(pH_list)
        charge_number_map = self.get_charge_number_map()
        def formal_charge(particle_name):
            tpl = self.db.get_template(name=particle_name, 
                                       pmb_type="particle")
            state = self.db.get_template(name=tpl.initial_state,
                                         pmb_type="particle_state")
            return charge_number_map[state.es_type]
        Z_HH = []
        for pH in pH_list:
            Z = 0.0
            for particle, multiplicity in particle_counts.items():
                if particle in pka_set:
                    pka = pka_set[particle]["pka_value"]
                    acidity = pka_set[particle]["acidity"]
                    if acidity == "acidic":
                        psi = -1
                    elif acidity == "basic":
                        psi = +1
                    else:
                        raise ValueError(f"Unknown acidity '{acidity}' for particle '{particle}'")
                    charge = psi / (1.0 + 10.0 ** (psi * (pH - pka)))
                    Z += multiplicity * charge
                else:
                    Z += multiplicity * formal_charge(particle)
            Z_HH.append(Z)
        return Z_HH   

    def calculate_HH_Donnan(self, c_macro, c_salt, pH_list=None, pka_set=None):
        """
        Computes macromolecular charges using the Henderson–Hasselbalch equation
        coupled to ideal Donnan partitioning.

        Args:
            c_macro ('dict'):
                Mapping of macromolecular species names to their concentrations
                in the system:
                '{molecule_name: concentration}'.
                
            c_salt ('float' or 'pint.Quantity'):
                Salt concentration in the reservoir.

            pH_list ('list[float]', optional):
                List of pH values in the reservoir at which the calculation is
                performed. If 'None', 50 equally spaced values between 2 and 12
                are used.

            pka_set ('dict', optional):
                Dictionary defining the acid–base properties of titratable particle
                types:
                '{particle_name: {"pka_value": float, "acidity": "acidic" | "basic"}}'.
                If 'None', the pKa set is taken from the pyMBE database.

        Returns:
            'dict':
                Dictionary containing:
                - '"charges_dict"' ('dict'):
                    Mapping '{molecule_name: list}' of Henderson–Hasselbalch–Donnan
                    charges evaluated at each pH value.
                - '"pH_system_list"' ('list[float]'):
                    Effective pH values inside the system phase after Donnan
                    partitioning.
                - '"partition_coefficients"' ('list[float]'):
                    Partition coefficients of monovalent cations at each pH value.

        Notes:
            - This method assumes **ideal Donnan equilibrium** and **monovalent salt**.
            - The ionic strength of the reservoir includes both salt and
            pH-dependent H⁺/OH⁻ contributions.
            - All charged macromolecular species present in the system must be
            included in 'c_macro'; missing species will lead to incorrect results.
            - The nonlinear Donnan equilibrium equation is solved using a scalar
            root finder ('brentq') in logarithmic form for numerical stability.
            - This method is intended for **two-phase systems**; for single-phase
            systems use 'calculate_HH' instead.
        """
        if pH_list is None:
            pH_list=np.linspace(2,12,50)
        if pka_set is None:
            pka_set=self.get_pka_set() 
        self._check_pka_set(pka_set=pka_set)
        partition_coefficients_list = []
        pH_system_list = []
        Z_HH_Donnan={}
        for key in c_macro:
            Z_HH_Donnan[key] = []
        def calc_charges(c_macro, pH):
            """
            Calculates the charges of the different kinds of molecules according to the Henderson-Hasselbalch equation.

            Args:
                c_macro ('dict'): 
                    {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 

                pH ('float'): 
                    pH-value that is used in the HH equation.

            Returns:
                ('dict'): 
                    {"molecule_name": charge}
            """
            charge = {}
            for name in c_macro:
                charge[name] = self.calculate_HH(name, [pH], pka_set)[0]
            return charge

        def calc_partition_coefficient(charge, c_macro):
            """
            Calculates the partition coefficients of positive ions according to the ideal Donnan theory.

            Args:
                charge ('dict'): 
                    {"molecule_name": charge}

                c_macro ('dict'): 
                    {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 
            """
            nonlocal ionic_strength_res
            charge_density = 0.0
            for key in charge:
                charge_density += charge[key] * c_macro[key]
            return (-charge_density / (2 * ionic_strength_res) + np.sqrt((charge_density / (2 * ionic_strength_res))**2 + 1)).magnitude
        for pH_value in pH_list:    
            # calculate the ionic strength of the reservoir
            if pH_value <= 7.0:
                ionic_strength_res = 10 ** (-pH_value) * self.units.mol/self.units.l + c_salt 
            elif pH_value > 7.0:
                ionic_strength_res = 10 ** (-(14-pH_value)) * self.units.mol/self.units.l + c_salt
            #Determine the partition coefficient of positive ions by solving the system of nonlinear, coupled equations
            #consisting of the partition coefficient given by the ideal Donnan theory and the Henderson-Hasselbalch equation.
            #The nonlinear equation is formulated for log(xi) since log-operations are not supported for RootResult objects.
            equation = lambda logxi: logxi - np.log10(calc_partition_coefficient(calc_charges(c_macro, pH_value - logxi), c_macro))
            logxi = scipy.optimize.root_scalar(equation, bracket=[-1e2, 1e2], method="brentq")
            partition_coefficient = 10**logxi.root
            charges_temp = calc_charges(c_macro, pH_value-np.log10(partition_coefficient))
            for key in c_macro:
                Z_HH_Donnan[key].append(charges_temp[key])
            pH_system_list.append(pH_value - np.log10(partition_coefficient))
            partition_coefficients_list.append(partition_coefficient)
        return {"charges_dict": Z_HH_Donnan, "pH_system_list": pH_system_list, "partition_coefficients": partition_coefficients_list}

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
        return self.simulation_engine.calculate_net_charge(object_name,
                                                           pmb_type,
                                                           dimensionless)
    
    def center_object_in_simulation_box(self, instance_id, box_l,pmb_type):
        """
        Centers a pyMBE object instance in the simulation box of an ESPResSo system.
        The object is translated such that its center of mass coincides with the
        geometric center of the ESPResSo simulation box.

        Args:
            instance_id ('int'):
                ID of the pyMBE object instance to be centered.
            
            box_l('list[float,float,float]'): list of floats with the dimensions of the box

            pmb_type ('str'):
                Type of the pyMBE object.

        Notes:
            - Works for both cubic and non-cubic simulation boxes.
        """
        inst = self.db.get_instance(instance_id=instance_id,
                                    pmb_type=pmb_type)
        center_of_mass = self.calculate_center_of_mass(instance_id=instance_id,
                                                       pmb_type=pmb_type)
        box_center = [box_l[0]/2.0,
                      box_l[1]/2.0,
                      box_l[2]/2.0]
        particle_id_list = self.get_particle_id_map(object_name=inst.name)["all"]
        for pid in particle_id_list:
            es_pos=self.db.get_instance(instance_id=pid,
                                    pmb_type='particle').position
            centered_position=es_pos - center_of_mass + box_center

            self.db._update_instance(instance_id=pid,
                                     pmb_type='particle',
                                     attribute='position',
                                     value=centered_position)
            if isinstance(self.simulation_engine, EspressoSimulation):
                self.simulation_engine._update_particle_position(
                    particle_id=pid, position=centered_position)

    def create_added_salt(self, box_l, cation_name, anion_name, c_salt):    
        """
        Creates a 'c_salt' concentration of 'cation_name' and 'anion_name' ions into the 'espresso_system'.

        Args:
            cation_name('str'): 'name' of a particle with a positive charge.
            anion_name('str'): 'name' of a particle with a negative charge.
            c_salt('float'): Salt concentration.
            
        Returns:
            c_salt_calculated('float'): Calculated salt concentration added to 'espresso_system'.
        """ 
        cation_tpl = self.db.get_template(pmb_type="particle",
                                          name=cation_name)
        cation_state = self.db.get_template(pmb_type="particle_state",
                                            name=cation_tpl.initial_state)
        cation_charge = cation_state.z
        anion_tpl = self.db.get_template(pmb_type="particle",
                                          name=anion_name)
        anion_state = self.db.get_template(pmb_type="particle_state",
                                            name=anion_tpl.initial_state)
        anion_charge = anion_state.z
        if cation_charge <= 0:
            raise ValueError(f'ERROR cation charge must be positive, charge {cation_charge}')
        if anion_charge >= 0:
            raise ValueError(f'ERROR anion charge must be negative, charge {anion_charge}')
        # Calculate the number of ions in the simulation box
        volume=self.units.Quantity(np.prod(box_l), 'reduced_length**3')
        if c_salt.check('[substance] [length]**-3'):
            N_ions= int((volume*c_salt.to('mol/reduced_length**3')*self.N_A).magnitude)
            c_salt_calculated=N_ions/(volume*self.N_A)
        elif c_salt.check('[length]**-3'):
            N_ions= int((volume*c_salt.to('reduced_length**-3')).magnitude)
            c_salt_calculated=N_ions/volume
        else:
            raise ValueError('Unknown units for c_salt, please provided it in [mol / volume] or [particle / volume]', c_salt)
        N_cation = N_ions*abs(anion_charge)
        N_anion = N_ions*abs(cation_charge)
        self.create_particle(box_l=box_l, 
                             name=cation_name, 
                             number_of_particles=N_cation)
        self.create_particle(box_l=box_l, 
                             name=anion_name, 
                             number_of_particles=N_anion)
        if c_salt_calculated.check('[substance] [length]**-3'):
            logging.info(f"added salt concentration of {c_salt_calculated.to('mol/L')} given by {N_cation} cations and {N_anion} anions")
        elif c_salt_calculated.check('[length]**-3'):
            logging.info(f"added salt concentration of {c_salt_calculated.to('reduced_length**-3')} given by {N_cation} cations and {N_anion} anions")
        return c_salt_calculated

    def create_bond(self, particle_id1, particle_id2, use_default_bond=False):
        """
        Creates a bond between two particle instances in an ESPResSo system and registers it in the pyMBE database.

        This method performs the following steps:
            1. Retrieves the particle instances corresponding to 'particle_id1' and 'particle_id2' from the database.
            2. Retrieves or creates the corresponding ESPResSo bond instance using the bond template.
            3. Adds the ESPResSo bond instance to the ESPResSo system if it was newly created.
            4. Adds the bond to the first particle's bond list in ESPResSo.
            5. Creates a 'BondInstance' in the database and registers it.

        Args:
            particle_id1 ('int'): 
                pyMBE and ESPResSo ID of the first particle.

            particle_id2 ('int'): 
                pyMBE and ESPResSo ID of the second particle.

            use_default_bond ('bool', optional): 
                If True, use a default bond template if no specific template exists. Defaults to False.

        Returns:
            ('int'): 
                bond_id of the bond instance created in the pyMBE database.
        """
        particle_inst_1 = self.db.get_instance(pmb_type="particle",
                                               instance_id=particle_id1)
        particle_inst_2 = self.db.get_instance(pmb_type="particle",
                                               instance_id=particle_id2)
        bond_tpl = self.get_bond_template(particle_name1=particle_inst_1.name,
                                          particle_name2=particle_inst_2.name,
                                          use_default_bond=use_default_bond)
        bond_id = self.db._propose_instance_id(pmb_type="bond")
        pmb_bond_instance = BondInstance(bond_id=bond_id,
                                         name=bond_tpl.name,
                                         particle_id1=particle_id1,
                                         particle_id2=particle_id2)
        self.db._register_instance(instance=pmb_bond_instance)

    def create_counterions(self, object_name, cation_name, anion_name, box_l):
        """
        Creates particles of 'cation_name' and 'anion_name' in 'espresso_system' to counter the net charge of 'object_name'.
        
        Args:
            object_name ('str'): 
                'name' of a pyMBE object.

            cation_name ('str'): 
                'name' of a particle with a positive charge.

            anion_name ('str'): 
                'name' of a particle with a negative charge.
            
            box_l('list[float,float,float]'): list of floats with the dimensions of the box

        Returns: 
            ('dict'): 
                {"name": number}

        Notes:
            This function currently does not support the creation of counterions for hydrogels.
        """ 
        cation_tpl = self.db.get_template(pmb_type="particle",
                                          name=cation_name)
        cation_state = self.db.get_template(pmb_type="particle_state",
                                            name=cation_tpl.initial_state)
        cation_charge = cation_state.z
        anion_tpl = self.db.get_template(pmb_type="particle",
                                          name=anion_name)
        anion_state = self.db.get_template(pmb_type="particle_state",
                                            name=anion_tpl.initial_state)
        anion_charge = anion_state.z
        object_ids = self.get_particle_id_map(object_name=object_name)["all"]
        counterion_number={}
        object_charge={}
        for name in ['positive', 'negative']:
            object_charge[name]=0
        for id in object_ids:
            object_name = self.db.get_instance(pmb_type="particle", 
                                               instance_id=id).name
            object_tpl = self.db.get_template(pmb_type="particle",
                                              name=object_name)
            object_state = self.db.get_template(pmb_type="particle_state",
                                                name=object_tpl.initial_state)
            object_z = object_state.z
            if object_z > 0:
                object_charge['positive']+=1*(np.abs(object_z ))
            elif object_z < 0:
                object_charge['negative']+=1*(np.abs(object_z ))
        if object_charge['positive'] % abs(anion_charge) == 0:
            counterion_number[anion_name]=int(object_charge['positive']/abs(anion_charge))
        else:
            raise ValueError('The number of positive charges in the pmb_object must be divisible by the  charge of the anion')
        if object_charge['negative'] % abs(cation_charge) == 0:
            counterion_number[cation_name]=int(object_charge['negative']/cation_charge)
        else:
            raise ValueError('The number of negative charges in the pmb_object must be divisible by the  charge of the cation')
        if counterion_number[cation_name] > 0: 
            self.create_particle(box_l=box_l, 
                                 name=cation_name, 
                                 number_of_particles=counterion_number[cation_name])
        else:
            counterion_number[cation_name]=0
        if counterion_number[anion_name] > 0:
            self.create_particle(box_l=box_l, 
                                 name=anion_name, 
                                 number_of_particles=counterion_number[anion_name])
        else:
            counterion_number[anion_name] = 0
        logging.info('the following counter-ions have been created: ')
        for name in counterion_number.keys():
            logging.info(f'Ion type: {name} created number: {counterion_number[name]}')
        return counterion_number


    def create_hydrogel(self, name, box_l, use_default_bond=False, gen_angle=False):
        """ 
        Creates a hydrogel in espresso_system using a pyMBE hydrogel template given by 'name'

        Args:
            box_l('list[float,float,float]'): list of floats with the dimensions of the box

            name ('str'): 
                name of the hydrogel template in the pyMBE database.

            use_default_bond ('bool', optional): 
                If True, use a default bond template if no specific template exists. Defaults to False.

            gen_angle ('bool', optional):
                If True, generate angle potentials for the internal hydrogel
                chains and, when explicitly defined, for all crosslinker-adjacent
                triplets. Defaults to False.

        Returns:
            ('int'): id of the hydrogel instance created.
        """
        if not self.db._has_template(name=name, pmb_type="hydrogel"):
            raise ValueError(f"Hydrogel template with name '{name}' is not defined in the pyMBE database.")
        hydrogel_tpl = self.db.get_template(pmb_type="hydrogel",
                                            name=name)
        assembly_id = self.db._propose_instance_id(pmb_type="hydrogel")
        # Create the nodes
        nodes = {}
        hydrogel_angle_centers = set()
        node_topology = hydrogel_tpl.node_map
        for node in node_topology:
            node_index = node.lattice_index
            node_name = node.particle_name
            node_pos, node_id = self._create_hydrogel_node(node_index=node_index,
                                                          node_name=node_name,
                                                          box_l=box_l)
            node_label = self.lattice_builder._create_node_label(node_index=node_index)
            nodes[node_label] = {"name": node_name, "id": node_id, "pos": node_pos} 
            self.db._update_instance(instance_id=node_id,
                                     pmb_type="particle",
                                     attribute="assembly_id",
                                     value=assembly_id)
        for hydrogel_chain in hydrogel_tpl.chain_map:
            molecule_id = self._create_hydrogel_chain(hydrogel_chain=hydrogel_chain,
                                                      nodes=nodes,
                                                      box_l=box_l,
                                                      use_default_bond=use_default_bond,
                                                      gen_angle=gen_angle,
                                                      )
            self.db._update_instance(instance_id=molecule_id,
                                     pmb_type="molecule",
                                     attribute="assembly_id",
                                     value=assembly_id)
            if gen_angle:
                residue_ids = self.db._find_instance_ids_by_attribute(pmb_type="residue",
                                                                      attribute="molecule_id",
                                                                      value=molecule_id)
                first_residue_id = min(residue_ids)
                last_residue_id = max(residue_ids)
                first_residue = self.db.get_instance(pmb_type="residue",
                                                     instance_id=first_residue_id)
                last_residue = self.db.get_instance(pmb_type="residue",
                                                    instance_id=last_residue_id)
                first_central_bead_name = self.db.get_template(pmb_type="residue",
                                                               name=first_residue.name).central_bead
                last_central_bead_name = self.db.get_template(pmb_type="residue",
                                                              name=last_residue.name).central_bead
                particle_instances = self.db.get_instances(pmb_type="particle")
                first_residue_particle_ids = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                                                      attribute="residue_id",
                                                                                      value=first_residue_id)
                last_residue_particle_ids = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                                                     attribute="residue_id",
                                                                                     value=last_residue_id)
                first_bead_id = None
                for particle_id in first_residue_particle_ids:
                    if particle_instances[particle_id].name == first_central_bead_name:
                        first_bead_id = particle_id
                        break

                last_bead_id = None
                for particle_id in last_residue_particle_ids:
                    if particle_instances[particle_id].name == last_central_bead_name:
                        last_bead_id = particle_id
                        break
                node_start_label = self.lattice_builder._create_node_label(hydrogel_chain.node_start)
                node_end_label = self.lattice_builder._create_node_label(hydrogel_chain.node_end)
                hydrogel_angle_centers.update({
                    nodes[node_start_label]["id"],
                    nodes[node_end_label]["id"],
                    first_bead_id,
                    last_bead_id,
                })
        self.db._propagate_id(root_type="hydrogel", 
                                root_id=assembly_id, 
                                attribute="assembly_id", 
                                value=assembly_id)
        if gen_angle:
            self._generate_hydrogel_crosslinker_angles(
                                                       central_particle_ids=hydrogel_angle_centers)
        # Register an hydrogel instance in the pyMBE databasegit 
        self.db._register_instance(HydrogelInstance(name=name,
                                                    assembly_id=assembly_id))
        return assembly_id


    def create_molecule(self, name, number_of_molecules, box_l, list_of_first_residue_positions=None, backbone_vector=None, use_default_bond=False, reverse_residue_order = False, gen_angle=False):
        """
        Creates instances of a given molecule template name into ESPResSo.

        Args:
            name ('str'): 
                Label of the molecule type to be created. 'name'.

            box_l('list[float,float,float]'): list of floats with the dimensions of the box

            number_of_molecules ('int'): 
                Number of molecules or peptides of type 'name' to be created.

            list_of_first_residue_positions ('list', optional): 
                List of coordinates where the central bead of the first_residue_position will be created, random by default.

            backbone_vector ('list' of 'float'): 
                Backbone vector of the molecule, random by default. Central beads of the residues in the 'residue_list' are placed along this vector. 

            use_default_bond('bool', optional): 
                Controls if a bond of type 'default' is used to bond particles with undefined bonds in the pyMBE database.

            reverse_residue_order('bool', optional): 
                Creates residues in reverse sequential order than the one defined in the molecule template. Defaults to False.

        Returns:
            ('list' of 'int'): 
                List with the 'molecule_id' of the pyMBE molecule instances created into 'espresso_system'.

        Notes:
            - This function can be used to create both molecules and peptides.    
        """
        pmb_type = self._get_template_type(name=name,
                                           allowed_types={"molecule", "peptide"})
        if number_of_molecules <= 0:
            return {}
        if list_of_first_residue_positions is not None:
            for item in list_of_first_residue_positions:
                if not isinstance(item, list):
                    raise ValueError("The provided input position is not a nested list. Should be a nested list with elements of 3D lists, corresponding to xyz coord.")
                elif len(item) != 3:
                    raise ValueError("The provided input position is formatted wrong. The elements in the provided list does not have 3 coordinates, corresponding to xyz coord.")

            if len(list_of_first_residue_positions) != number_of_molecules:
                raise ValueError(f"Number of positions provided in {list_of_first_residue_positions} does not match number of molecules desired, {number_of_molecules}")
        # Generate an arbitrary random unit vector
        if backbone_vector is None:
            backbone_vector = self.generate_random_points_in_a_sphere(center=[0,0,0],
                                                                      radius=1, 
                                                                      n_samples=1,
                                                                      on_surface=True)[0]
        else:
            backbone_vector = np.array(backbone_vector)
        first_residue = True
        molecule_tpl = self.db.get_template(pmb_type=pmb_type,
                                            name=name)
        if reverse_residue_order:
            residue_list = molecule_tpl.residue_list[::-1]
        else:
            residue_list = molecule_tpl.residue_list
        pos_index = 0 
        molecule_ids = []
        for n_mol in range(number_of_molecules):        
            molecule_id = self.db._propose_instance_id(pmb_type=pmb_type)
            for residue in residue_list:
                if first_residue:
                    if list_of_first_residue_positions is None:
                        central_bead_pos = None
                    else:
                        central_bead_pos = [np.array(list_of_first_residue_positions[n_mol])]
                            
                    residue_id = self.create_residue(name=residue,
                                                     box_l=box_l, 
                                                     central_bead_position=central_bead_pos,  
                                                     use_default_bond= use_default_bond, 
                                                     backbone_vector=backbone_vector)
                    
                    # Add molecule_id to the residue instance and all particles associated
                    self.db._propagate_id(root_type="residue", 
                                          root_id=residue_id,
                                          attribute="molecule_id", 
                                          value=molecule_id)
                    particle_ids_in_residue = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                                                      attribute="residue_id",
                                                                                      value=residue_id)
                    prev_central_bead_id = particle_ids_in_residue[0]
                    prev_central_bead_name = self.db.get_instance(pmb_type="particle", 
                                                                  instance_id=prev_central_bead_id).name
                    prev_central_bead_pos = self.db.get_instance(pmb_type="particle", 
                                                                  instance_id=prev_central_bead_id).position
                    # prev_central_bead_pos = espresso_system.part.by_id(prev_central_bead_id).pos
                    first_residue = False          
                else:
                    
                    # Calculate the starting position of the new residue
                    residue_tpl = self.db.get_template(pmb_type="residue",
                                                       name=residue)
                    lj_parameters = self.get_lj_parameters(particle_name1=prev_central_bead_name,
                                                           particle_name2=residue_tpl.central_bead)
                    bond_tpl = self.get_bond_template(particle_name1=prev_central_bead_name,
                                                      particle_name2=residue_tpl.central_bead,
                                                      use_default_bond=use_default_bond)
                    l0 = hf.calculate_initial_bond_length(lj_parameters=lj_parameters,
                                                          bond_type=bond_tpl.bond_type,
                                                          bond_parameters=bond_tpl.get_parameters(ureg=self.units))
                    central_bead_pos = prev_central_bead_pos+backbone_vector*l0
                    # Create the residue
                    residue_id = self.create_residue(name=residue, 
                                                     box_l=box_l, 
                                                     central_bead_position=[central_bead_pos],
                                                     use_default_bond= use_default_bond, 
                                                     backbone_vector=backbone_vector)
                    # Add molecule_id to the residue instance and all particles associated
                    self.db._propagate_id(root_type="residue", 
                                          root_id=residue_id, 
                                          attribute="molecule_id", 
                                          value=molecule_id)
                    particle_ids_in_residue = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                                                      attribute="residue_id",
                                                                                      value=residue_id)
                    central_bead_id = particle_ids_in_residue[0]

                    # Bond the central beads of the new and previous residues
                    self.create_bond(particle_id1=prev_central_bead_id,
                                     particle_id2=central_bead_id,
                                     use_default_bond=use_default_bond)
                    
                    prev_central_bead_id = central_bead_id                    
                    prev_central_bead_name = self.db.get_instance(pmb_type="particle", instance_id=central_bead_id).name
                    prev_central_bead_pos =central_bead_pos
            # Create a Peptide or Molecule instance and register it on the pyMBE database
            if pmb_type == "molecule":
                inst = MoleculeInstance(molecule_id=molecule_id,
                                        name=name)
            elif pmb_type == "peptide":
                inst = PeptideInstance(name=name,
                                       molecule_id=molecule_id)
            self.db._register_instance(inst)
            if gen_angle:
                self._generate_angles_for_entity(
                    entity_id=molecule_id,
                    entity_id_col='molecule_id')
            first_residue = True
            pos_index+=1
            molecule_ids.append(molecule_id)
        return molecule_ids
    
    def create_particle(self, name, box_l, number_of_particles, position=None, fix=False):
        """
        Creates one or more particles in an ESPResSo system based on the particle template in the pyMBE database.
        
        Args:
            name ('str'): 
                Label of the particle template in the pyMBE database. 
            
            box_l('list[float,float,float]'): list of floats with the dimensions of the box

            number_of_particles ('int'): 
                Number of particles to be created.

            position (list of ['float','float','float'], optional): 
                Initial positions of the particles. If not given, particles are created in random positions. Defaults to None.

            fix ('bool', optional): 
                Controls if the particle motion is frozen in the integrator, it is used to create rigid objects. Defaults to False.

        Returns:
            ('list' of 'int'): 
                List with the ids of the particles created into 'espresso_system'.
        """       
        if number_of_particles <=0:
            return []
        if not self.db._has_template(name=name, pmb_type="particle"):
            raise ValueError(f"Particle template with name '{name}' is not defined in the pyMBE database.")
        
        part_tpl = self.db.get_template(pmb_type="particle",
                                        name=name)
        part_state = self.db.get_template(pmb_type="particle_state",
                                         name=part_tpl.initial_state)
        name_state=part_state.name

        if fix is False:
            fix=[fix]*3

        created_pid_list=[]
        for index in range(number_of_particles):
            if position is None:
                particle_position = self.rng.random((1, 3))[0] *np.copy(box_l)
            else:
                particle_position = np.array(position[index])
            
            particle_id = self.db._propose_instance_id(pmb_type="particle")
            created_pid_list.append(particle_id)
            part_inst = ParticleInstance(name=name,
                                         particle_id=particle_id,
                                         initial_state=name_state,
                                         position=particle_position,
                                         fix=fix)
            self.db._register_instance(part_inst)
                              
        return created_pid_list

    def create_protein(self, name, number_of_proteins, box_l, topology_dict):
        """
        Creates one or more protein molecules in an ESPResSo system based on the 
        protein template in the pyMBE database and a provided topology.

        Args:
            name (str):
                Name of the protein template stored in the pyMBE database.
            
            number_of_proteins (int):
                Number of protein molecules to generate.

            box_l('list[float,float,float]'): list of floats with the dimensions of the box  
            
            topology_dict (dict):
                Dictionary defining the internal structure of the protein. Expected format:
                    {"ResidueName1": {"initial_pos": np.ndarray,
                                      "chain_id": int,
                                      "radius": float},
                     "ResidueName2": { ... },
                        ...
                    }
                The '"initial_pos"' entry is required and represents the residue’s
                reference coordinates before shifting to the protein's center-of-mass.

        Returns:
            ('list' of 'int'): 
                List of the molecule_id of the Protein instances created into ESPResSo.

        Notes:
            - Particles are created using 'create_particle()' with 'fix=True',
            meaning they are initially immobilized.
            - The function assumes all residues in 'topology_dict' correspond to
            particle templates already defined in the pyMBE database.
            - Bonds between residues are not created here; it assumes a rigid body representation of the protein.
        """
        if number_of_proteins <= 0:
            return
        if not self.db._has_template(name=name, pmb_type="protein"):
            raise ValueError(f"Protein template with name '{name}' is not defined in the pyMBE database.")
        protein_tpl = self.db.get_template(pmb_type="protein", name=name)
        box_half = box_l[0] / 2.0
        # Create protein
        mol_ids = []
        for _ in range(number_of_proteins):
            # create a molecule identifier in pyMBE
            molecule_id = self.db._propose_instance_id(pmb_type="protein")
            # place protein COM randomly
            protein_center = self.generate_coordinates_outside_sphere(radius=1,
                                                                      max_dist=box_half,
                                                                      n_samples=1,
                                                                      center=[box_half]*3)[0]
            residues = hf.get_residues_from_topology_dict(topology_dict=topology_dict,
                                                         model=protein_tpl.model)
            # CREATE RESIDUES + PARTICLES
            for _, rdata in residues.items():
                base_resname = rdata["resname"]  
                residue_name = f"AA-{base_resname}"
                # residue instance ID
                residue_id = self.db._propose_instance_id("residue")
                # register ResidueInstance
                self.db._register_instance(ResidueInstance(name=residue_name,
                                                           residue_id=residue_id,
                                                           molecule_id=molecule_id))
                # PARTICLE CREATION
                for bead_id in rdata["beads"]:
                    bead_type = re.split(r'\d+', bead_id)[0]
                    relative_pos = topology_dict[bead_id]["initial_pos"]
                    absolute_pos = relative_pos + protein_center
                    particle_id = self.create_particle(name=bead_type,
                                                       box_l=box_l,
                                                       number_of_particles=1,
                                                       position=[absolute_pos],
                                                       fix=[True,True,True])[0]
                    # update metadata
                    self.db._update_instance(instance_id=particle_id,
                                             pmb_type="particle",
                                             attribute="molecule_id",
                                             value=molecule_id)
                    self.db._update_instance(instance_id=particle_id,
                                             pmb_type="particle",
                                             attribute="residue_id",
                                             value=residue_id)
            protein_inst = ProteinInstance(name=name,
                                           molecule_id=molecule_id)
            self.db._register_instance(protein_inst)
            mol_ids.append(molecule_id)
        return mol_ids

    def create_residue(self, name, box_l, central_bead_position=None,use_default_bond=False, backbone_vector=None, gen_angle=False):
        """
        Creates a residue  into ESPResSo.

        Args:
            name ('str'): 
                Label of the residue type to be created. 

            central_bead_position ('list' of 'float'): 
                Position of the central bead.
            
            box_l('list[float,float,float]'): list of floats with the dimensions of the box

            use_default_bond ('bool'): 
                Switch to control if a bond of type 'default' is used to bond a particle whose bonds types are not defined in the pyMBE database.

            backbone_vector ('list' of 'float'): 
                Backbone vector of the molecule. All side chains are created perpendicularly to 'backbone_vector'.

        Returns:
            (int): 
                residue_id of the residue created.
        """
        if not self.db._has_template(name=name, pmb_type="residue"):
            raise ValueError(f"Residue template with name '{name}' is not defined in the pyMBE database.")
        res_tpl = self.db.get_template(pmb_type="residue",
                                       name=name)
        # Assign a residue_id
        residue_id = self.db._propose_instance_id(pmb_type="residue")
        res_inst = ResidueInstance(name=name,
                                   residue_id=residue_id)
        self.db._register_instance(res_inst)
        # create the principal bead   
        central_bead_name = res_tpl.central_bead 
        central_bead_id = self.create_particle(name=central_bead_name,
                                               box_l=box_l,
                                               position=central_bead_position,
                                               number_of_particles = 1)[0]
        
        central_bead_position = self.db.get_instance(pmb_type="particle", 
                                                                  instance_id=central_bead_id).position
        # # central_bead_position=espresso_system.part.by_id(central_bead_id).pos

        # Assigns residue_id to the central_bead particle created.
        self.db._update_instance(pmb_type="particle",
                                 instance_id=central_bead_id,
                                 attribute="residue_id",
                                 value=residue_id)
        
        # create the lateral beads  
        side_chain_list = res_tpl.side_chains
        side_chain_beads_ids = []
        for side_chain_name in side_chain_list:
            pmb_type = self._get_template_type(name=side_chain_name,
                                               allowed_types={"particle", "residue"})
            if pmb_type == 'particle':
                lj_parameters = self.get_lj_parameters(particle_name1=central_bead_name,
                                                       particle_name2=side_chain_name)
                bond_tpl = self.get_bond_template(particle_name1=central_bead_name,
                                                  particle_name2=side_chain_name,
                                                  use_default_bond=use_default_bond)
                l0 = hf.calculate_initial_bond_length(lj_parameters=lj_parameters,
                                                      bond_type=bond_tpl.bond_type,
                                                      bond_parameters=bond_tpl.get_parameters(ureg=self.units))               
                if backbone_vector is None:
                    bead_position=self.generate_random_points_in_a_sphere(center=central_bead_position, 
                                                                radius=l0, 
                                                                n_samples=1,
                                                                on_surface=True)[0]
                else:
                    bead_position=central_bead_position+self.generate_trial_perpendicular_vector(vector=np.array(backbone_vector),
                                                                                                magnitude=l0)
                    
                side_bead_id = self.create_particle(name=side_chain_name, 
                                                    box_l=box_l,
                                                    position=[bead_position], 
                                                    number_of_particles=1)[0]
                side_chain_beads_ids.append(side_bead_id)
                self.db._update_instance(pmb_type="particle",
                                         instance_id=side_bead_id,
                                         attribute="residue_id",
                                         value=residue_id)
                self.create_bond(particle_id1=central_bead_id,
                                 particle_id2=side_bead_id,
                                 use_default_bond=use_default_bond)
                
            elif pmb_type == 'residue':

                side_residue_tpl = self.db.get_template(name=side_chain_name,
                                                        pmb_type=pmb_type)
                central_bead_side_chain = side_residue_tpl.central_bead
                lj_parameters = self.get_lj_parameters(particle_name1=central_bead_name,
                                                       particle_name2=central_bead_side_chain)
                bond_tpl = self.get_bond_template(particle_name1=central_bead_name,
                                                  particle_name2=central_bead_side_chain,
                                                  use_default_bond=use_default_bond)
                l0 = hf.calculate_initial_bond_length(lj_parameters=lj_parameters,
                                                      bond_type=bond_tpl.bond_type,
                                                      bond_parameters=bond_tpl.get_parameters(ureg=self.units))
                if backbone_vector is None:
                    residue_position=self.generate_random_points_in_a_sphere(center=central_bead_position, 
                                                                radius=l0, 
                                                                n_samples=1,
                                                                on_surface=True)[0]
                else:
                    residue_position=central_bead_position+self.generate_trial_perpendicular_vector(vector=backbone_vector,
                                                                                                    magnitude=l0)
                side_residue_id = self.create_residue(name=side_chain_name, 
                                                      box_l=box_l,
                                                      central_bead_position=[residue_position],
                                                      use_default_bond=use_default_bond)
                # Find particle ids of the inner residue
                side_chain_beads_ids = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                                               attribute="residue_id",
                                                                               value=side_residue_id)
                # Change the residue_id of the residue in the side chain to the one of the outer residue
                for particle_id in side_chain_beads_ids:
                    self.db._update_instance(instance_id=particle_id,
                                             pmb_type="particle",
                                             attribute="residue_id",
                                             value=residue_id)
                # Remove the instance of the inner residue
                self.db.delete_instance(pmb_type="residue",
                                        instance_id=side_residue_id)
                self.create_bond(particle_id1=central_bead_id,
                                 particle_id2=side_chain_beads_ids[0],
                                 use_default_bond=use_default_bond)
        if gen_angle:
            self._generate_angles_for_entity(
                                             entity_id=residue_id,
                                             entity_id_col="residue_id")
        return  residue_id

    def define_bond(self, bond_type, bond_parameters, particle_pairs):
        """
        Defines bond templates for each particle pair in 'particle_pairs' in the pyMBE database.

        Args:
            bond_type ('str'): 
                label to identify the potential to model the bond.

            bond_parameters ('dict'): 
                parameters of the potential of the bond.

            particle_pairs ('lst'): 
                list of the 'names' of the 'particles' to be bonded.

        Notes:
            -Currently, only HARMONIC and FENE bonds are supported.
            - For a HARMONIC bond the dictionary must contain the following parameters:
                - k ('pint.Quantity')      : Magnitude of the bond. It should have units of energy/length**2 
                using the 'pmb.units' UnitRegistry.
                - r_0 ('pint.Quantity')    : Equilibrium bond length. It should have units of length using 
                the 'pmb.units' UnitRegistry.
           - For a FENE bond the dictionary must contain the same parameters as for a HARMONIC bond and:              
                - d_r_max ('pint.Quantity'): Maximal stretching length for FENE. It should have 
                units of length using the 'pmb.units' UnitRegistry. Default 'None'.
        """
        self._check_bond_inputs(bond_parameters=bond_parameters,
                                bond_type=bond_type)
        parameters_expected_dimensions={"r_0": "length",
                                        "k": "energy/length**2",
                                        "d_r_max": "length"}

        parameters_tpl = {}
        for key in bond_parameters.keys():
            parameters_tpl[key]= PintQuantity.from_quantity(q=bond_parameters[key],
                                                            expected_dimension=parameters_expected_dimensions[key],
                                                            ureg=self.units)

        bond_names=[]
        for particle_name1, particle_name2 in particle_pairs:
            
            tpl = BondTemplate(particle_name1=particle_name1,
                               particle_name2=particle_name2,
                               parameters=parameters_tpl,
                               bond_type=bond_type)
            tpl._make_name()
            if tpl.name in bond_names:
                raise RuntimeError(f"Bond {tpl.name} has already been defined, please check the list of particle pairs")
            bond_names.append(tpl.name)
            self.db._register_template(tpl)

    
    def define_default_bond(self, bond_type, bond_parameters):
        """
        Defines a bond template as a "default" template in the pyMBE database.
        
        Args:
            bond_type ('str'): 
                label to identify the potential to model the bond.

            bond_parameters ('dict'): 
                parameters of the potential of the bond.
            
        Notes:
            - Currently, only harmonic and FENE bonds are supported. 
        """
        self._check_bond_inputs(bond_parameters=bond_parameters,
                                bond_type=bond_type)
        parameters_expected_dimensions={"r_0": "length",
                                        "k": "energy/length**2",
                                        "d_r_max": "length"}
        parameters_tpl = {}
        for key in bond_parameters.keys():
            parameters_tpl[key]= PintQuantity.from_quantity(q=bond_parameters[key],
                                                            expected_dimension=parameters_expected_dimensions[key],
                                                            ureg=self.units)
        tpl = BondTemplate(parameters=parameters_tpl,
                               bond_type=bond_type)
        tpl.name = "default"
        self.db._register_template(tpl)

    def define_angular_potential(self, angle_type, angle_parameters, particle_triplets):
        """
        Defines angle potential templates for each particle triplet in `particle_triplets`.

        Args:
            angle_type ('str'):
                Type of angle potential. Supported: "harmonic", "cosine", "harmonic_cosine".

            angle_parameters ('dict'):
                Parameters of the angle potential. Must contain:
                    - "k" ('pint.Quantity'): Bending stiffness with dimensions of energy.
                    - "phi_0" ('float'): Equilibrium angle in radians.

            particle_triplets ('list[tuple[str,str,str]]'):
                List of (side_particle1, central_particle, side_particle2) triplets.
        """
        valid_angle_types = ["harmonic", "cosine", "harmonic_cosine"]
        if angle_type not in valid_angle_types:
            raise NotImplementedError(f"Angle potential type '{angle_type}' currently not implemented in pyMBE, accepted types are {valid_angle_types}")

        if "k" not in angle_parameters:
            raise ValueError("Magnitude of the angle potential (k) is missing")
        if "phi_0" not in angle_parameters:
            raise ValueError("Equilibrium angle (phi_0) is missing")

        parameters_tpl = {"k": PintQuantity.from_quantity(q=angle_parameters["k"],
                                                          expected_dimension="energy",
                                                          ureg=self.units),
                          "phi_0": PintQuantity.from_quantity(q=angle_parameters["phi_0"],
                                                              expected_dimension="dimensionless",
                                                              ureg=self.units),}
        angle_names = []
        for side1, central, side2 in particle_triplets:
            tpl = AngleTemplate(side_particle1=side1,
                                central_particle=central,
                                side_particle2=side2,
                                parameters=parameters_tpl,
                                angle_type=angle_type)
            tpl._make_name()
            if tpl.name in angle_names:
                raise RuntimeError(f"Angle {tpl.name} has already been defined, please check the list of particle triplets")
            angle_names.append(tpl.name)
            self.db._register_template(tpl)

    def define_default_angular_potential(self, angle_type, angle_parameters):
        """
        Defines an angle template as a "default" template in the pyMBE database.

        Args:
            angle_type ('str'):
                Type of angle potential. Supported: "harmonic", "cosine", "harmonic_cosine".

            angle_parameters ('dict'):
                Parameters of the angle potential (k, phi_0).
        """
        valid_angle_types = ["harmonic", "cosine", "harmonic_cosine"]
        if angle_type not in valid_angle_types:
            raise NotImplementedError(f"Angle potential type '{angle_type}' currently not implemented in pyMBE, accepted types are {valid_angle_types}")
        if "k" not in angle_parameters:
            raise ValueError("Magnitude of the angle potential (k) is missing")
        if "phi_0" not in angle_parameters:
            raise ValueError("Equilibrium angle (phi_0) is missing")
        parameters_tpl = {"k": PintQuantity.from_quantity(q=angle_parameters["k"],
                                                          expected_dimension="energy",
                                                          ureg=self.units),
                         "phi_0": PintQuantity.from_quantity(q=angle_parameters["phi_0"],
                                                             expected_dimension="dimensionless",
                                                             ureg=self.units),}
        tpl = AngleTemplate(parameters=parameters_tpl,
                            angle_type=angle_type)
        tpl.name = "default"
        self.db._register_template(tpl)

    def create_angular_potential(self, particle_id1, particle_id2, particle_id3, use_default_angle=False):
        """
        Creates an angle between three particle instances in an ESPResSo system
        and registers it in the pyMBE database.

        Args:
            particle_id1 ('int'): ID of the first side particle.
            particle_id2 ('int'): ID of the central particle.
            particle_id3 ('int'): ID of the second side particle.
            use_default_angle ('bool', optional): If True, use the default angle if no specific one is found.
        """
        particle_inst_1 = self.db.get_instance(pmb_type="particle", instance_id=particle_id1)
        particle_inst_2 = self.db.get_instance(pmb_type="particle", instance_id=particle_id2)
        particle_inst_3 = self.db.get_instance(pmb_type="particle", instance_id=particle_id3)

        # Verify that bonds exist between side particles and central particle
        bond_instances = self.db.get_instances(pmb_type="bond")
        bonded_pairs = set()
        for bond in bond_instances.values():
            pair = frozenset([bond.particle_id1, bond.particle_id2])
            bonded_pairs.add(pair)
        if frozenset([particle_id1, particle_id2]) not in bonded_pairs:
            raise ValueError(f"Cannot create angle: no bond exists between particle {particle_id1} and central particle {particle_id2}.")
        if frozenset([particle_id3, particle_id2]) not in bonded_pairs:
            raise ValueError(f"Cannot create angle: no bond exists between particle {particle_id3} and central particle {particle_id2}.")

        angle_tpl = self.get_angle_template(side_name1=particle_inst_1.name,
                                            central_name=particle_inst_2.name,
                                            side_name2=particle_inst_3.name,
                                            use_default_angle=use_default_angle)
        angle_id = self.db._propose_instance_id(pmb_type="angle")
        pmb_angle_instance = AngleInstance(angle_id=angle_id,
                                           name=angle_tpl.name,
                                           particle_id1=particle_id1,
                                           particle_id2=particle_id2,
                                           particle_id3=particle_id3)
        self.db._register_instance(instance=pmb_angle_instance)

    def get_angle_template(self, side_name1, central_name, side_name2, use_default_angle=False):
        """
        Retrieves an angle template connecting three particle templates.

        Args:
            side_name1 ('str'): Name of the first side particle.
            central_name ('str'): Name of the central particle.
            side_name2 ('str'): Name of the second side particle.
            use_default_angle ('bool', optional): If True, fall back to the default angle template.

        Returns:
            ('AngleTemplate'): The matching angle template.
        """
        angle_key = AngleTemplate.make_angle_key(side1=side_name1, central=central_name, side2=side_name2)
        try:
            return self.db.get_template(name=angle_key, pmb_type="angle")
        except ValueError:
            pass

        if use_default_angle:
            return self.db.get_template(name="default", pmb_type="angle")

        raise ValueError(f"No angle template found for '{side_name1}-{central_name}-{side_name2}', and default angles are deactivated.")

    def _generate_angles_for_entity(self, entity_id, entity_id_col):
        """
        Auto-generates angles from bond topology for an entity (molecule or residue).

        For each particle in the entity that has two or more bonded neighbors,
        this method finds all neighbor pairs and applies any matching angle potential.

        Args:
            entity_id ('int'): The molecule_id or residue_id to generate angles for.
            entity_id_col ('str'): Either "molecule_id" or "residue_id".
        """
        # Get all particle IDs for this entity
        particle_ids = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                               attribute=entity_id_col,
                                                               value=entity_id)
        if not particle_ids:
            return

        # Build neighbor map from bond instances
        neighbors = {pid: set() for pid in particle_ids}
        pid_set = set(particle_ids)
        bond_instances = self.db.get_instances(pmb_type="bond")
        for bond in bond_instances.values():
            i, j = bond.particle_id1, bond.particle_id2
            if i in pid_set and j in pid_set:
                neighbors[i].add(j)
                neighbors[j].add(i)

        # For each particle with 2+ neighbors, generate angles
        for j in particle_ids:
            nbs = sorted(neighbors[j])
            if len(nbs) < 2:
                continue

            for idx_i in range(len(nbs)):
                for idx_k in range(idx_i + 1, len(nbs)):
                    i = nbs[idx_i]
                    k = nbs[idx_k]
                    try:
                        self.create_angular_potential(particle_id1=i,
                                          particle_id2=j,
                                          particle_id3=k,
                                          use_default_angle=True)
                    except ValueError:
                        # No angle template defined for this triplet — skip
                        continue

    def define_hydrogel(self, name, node_map, chain_map):
        """
        Defines a hydrogel template in the pyMBE database.

        Args:
            name ('str'): 
                Unique label that identifies the 'hydrogel'.

            node_map ('list of dict'): 
                [{"particle_name": , "lattice_index": }, ... ]

            chain_map ('list of dict'): 
                [{"node_start": , "node_end": , "residue_list": , ... ]
        """
        # Sanity tests
        node_indices = {tuple(entry['lattice_index']) for entry in node_map}                
        chain_map_connectivity = set()
        for entry in chain_map:
            start = self.lattice_builder.node_labels[entry['node_start']]
            end = self.lattice_builder.node_labels[entry['node_end']]
            chain_map_connectivity.add((start,end))
        if self.lattice_builder.lattice.connectivity != chain_map_connectivity:
            raise ValueError("Incomplete hydrogel: A diamond lattice must contain correct 16 lattice index pairs")
        diamond_indices = {tuple(row) for row in self.lattice_builder.lattice.indices}
        if node_indices != diamond_indices:
            raise ValueError(f"Incomplete hydrogel: A diamond lattice must contain exactly 8 lattice indices, {diamond_indices} ")
        # Register information in the pyMBE database
        nodes=[]
        for entry in node_map:
            nodes.append(HydrogelNode(particle_name=entry["particle_name"],
                                      lattice_index=entry["lattice_index"]))
        chains=[]
        for chain in chain_map:
            chains.append(HydrogelChain(node_start=chain["node_start"],
                                        node_end=chain["node_end"],
                                        molecule_name=chain["molecule_name"]))
        tpl = HydrogelTemplate(name=name,
                               node_map=nodes,
                               chain_map=chains)
        self.db._register_template(tpl)

    def define_molecule(self, name, residue_list):
        """
        Defines a molecule template in the pyMBE database.

        Args:
            name('str'): 
                Unique label that identifies the 'molecule'.

            residue_list ('list' of 'str'): 
                List of the 'name's of the 'residue's  in the sequence of the 'molecule'.  
        """
        tpl = MoleculeTemplate(name=name,
                               residue_list=residue_list)
        self.db._register_template(tpl)

    def define_monoprototic_acidbase_reaction(self, particle_name, pka, acidity, metadata=None):
        """
        Defines an acid-base reaction for a monoprototic particle in the pyMBE database.

        Args:
            particle_name ('str'): 
                Unique label that identifies the particle template. 

            pka ('float'): 
                pka-value of the acid or base.

            acidity ('str'): 
                Identifies whether if the particle is 'acidic' or 'basic'.

            metadata ('dict', optional): 
                Additional information to be stored in the reaction. Defaults to None.
        """
        supported_acidities = ["acidic", "basic"]
        if acidity not in supported_acidities:
            raise ValueError(f"Unsupported acidity '{acidity}' for particle '{particle_name}'. Supported acidities are {supported_acidities}.")
        reaction_type = "monoprotic"
        if acidity == "basic":
            reaction_type += "_base"
        else:
            reaction_type += "_acid"
        reaction = Reaction(participants=[ReactionParticipant(particle_name=particle_name,
                                                              state_name=f"{particle_name}H", 
                                                              coefficient=-1),
                                          ReactionParticipant(particle_name=particle_name,
                                                              state_name=f"{particle_name}",
                                                              coefficient=1)],
                            reaction_type=reaction_type,
                            pK=pka,
                            metadata=metadata)
        self.db._register_reaction(reaction)

    def define_monoprototic_particle_states(self, particle_name, acidity):
        """
        Defines particle states for a monoprotonic particle template including the charges in each of its possible states. 

        Args:
            particle_name ('str'): 
                Unique label that identifies the particle template. 

            acidity ('str'): 
                Identifies whether the particle is 'acidic' or 'basic'.
        """
        acidity_valid_keys = ['acidic', 'basic']
        if not pd.isna(acidity):
            if acidity not in acidity_valid_keys:
                raise ValueError(f"Acidity {acidity} provided for particle name  {particle_name} is not supported. Valid keys are: {acidity_valid_keys}")
        if acidity == "acidic":
            states = [{"name": f"{particle_name}H", "z": 0}, 
                      {"name": f"{particle_name}",  "z": -1}]
            
        elif acidity == "basic":
            states = [{"name": f"{particle_name}H", "z": 1}, 
                      {"name": f"{particle_name}",  "z": 0}]
        self.define_particle_states(particle_name=particle_name, 
                                    states=states)

    def define_particle(self, name,  sigma, epsilon, z=0, acidity=pd.NA, pka=pd.NA, cutoff=pd.NA, offset=pd.NA):
        """
        Defines a particle template in the pyMBE database.

        Args:
            name('str'):
                 Unique label that identifies this particle type.  

            sigma('pint.Quantity'): 
                Sigma parameter used to set up Lennard-Jones interactions for this particle type. 

            epsilon('pint.Quantity'): 
                Epsilon parameter used to setup Lennard-Jones interactions for this particle tipe.

            z('int', optional): 
                Permanent charge number of this particle type. Defaults to 0.

            acidity('str', optional): 
                Identifies whether if the particle is 'acidic' or 'basic', used to setup constant pH simulations. Defaults to pd.NA.

            pka('float', optional):
                If 'particle' is an acid or a base, it defines its  pka-value. Defaults to pd.NA.

            cutoff('pint.Quantity', optional): 
                Cutoff parameter used to set up Lennard-Jones interactions for this particle type. Defaults to pd.NA.

            offset('pint.Quantity', optional): 
                Offset parameter used to set up Lennard-Jones interactions for this particle type. Defaults to pd.NA.
            
        Notes:
            - 'sigma', 'cutoff' and 'offset' must have a dimensitonality of '[length]' and should be defined using pmb.units.
            - 'epsilon' must have a dimensitonality of '[energy]' and should be defined using pmb.units.
            - 'cutoff' defaults to '2**(1./6.) reduced_length'. 
            - 'offset' defaults to 0.
            - For more information on 'sigma', 'epsilon', 'cutoff' and 'offset' check 'pmb.setup_lj_interactions()'.
        """ 
        # If 'cutoff' and 'offset' are not defined, default them to the following values
        if pd.isna(cutoff):
            cutoff=self.units.Quantity(2**(1./6.), "reduced_length")
        if pd.isna(offset):
            offset=self.units.Quantity(0, "reduced_length")
        # Define particle states
        if acidity is pd.NA:
            states = [{"name": f"{name}",  "z": z}]
            self.define_particle_states(particle_name=name, 
                                        states=states)
            initial_state = name
        else:
            self.define_monoprototic_particle_states(particle_name=name,
                                                  acidity=acidity)
            initial_state = f"{name}H"
            if pka is not pd.NA:
                self.define_monoprototic_acidbase_reaction(particle_name=name,
                                                           acidity=acidity,
                                                           pka=pka)
        tpl = ParticleTemplate(name=name, 
                               sigma=PintQuantity.from_quantity(q=sigma, expected_dimension="length", ureg=self.units), 
                               epsilon=PintQuantity.from_quantity(q=epsilon, expected_dimension="energy", ureg=self.units),
                               cutoff=PintQuantity.from_quantity(q=cutoff, expected_dimension="length", ureg=self.units), 
                               offset=PintQuantity.from_quantity(q=offset, expected_dimension="length", ureg=self.units),
                               initial_state=initial_state)
        self.db._register_template(tpl)
    
    def define_particle_states(self, particle_name, states):
        """
        Define the chemical states of an existing particle template.

        Args:
            particle_name ('str'):
                Name of a particle template. 

            states ('list' of 'dict'):
                List of dictionaries defining the particle states. Each dictionary
                must contain:
                - 'name' ('str'): Name of the particle state (e.g. '"H"', '"-"',
                '"neutral"').
                - 'z' ('int'): Charge number of the particle in this state.
                Example:
                states = [{"name": "AH", "z": 0},     # protonated
                         {"name": "A-", "z": -1}]    # deprotonated
        Notes:
            - Each state is assigned a unique Espresso 'es_type' automatically.
            - Chemical reactions (e.g. acid–base equilibria) are **not** created by
            this method and must be defined separately (e.g. via
            'set_particle_acidity()' or custom reaction definitions).
            - Particles without explicitly defined states are assumed to have a
            single, implicit state with their default charge.
        """
        for s in states:
            state = ParticleStateTemplate(particle_name=particle_name,
                                          name=s["name"],
                                          z=s["z"],
                                          es_type=self.propose_unused_type())
            self.db._register_template(state)

    def define_peptide(self, name, sequence, model):
        """
        Defines a peptide template in the pyMBE database.

        Args:
            name ('str'): 
                Unique label that identifies the peptide.

            sequence ('str'): 
                Sequence of the peptide.

            model ('str'): 
                Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.
        """
        valid_keys = ['1beadAA','2beadAA']
        if model not in valid_keys:
            raise ValueError('Invalid label for the peptide model, please choose between 1beadAA or 2beadAA')
        clean_sequence = hf.protein_sequence_parser(sequence=sequence)    
        residue_list = self._get_residue_list_from_sequence(sequence=clean_sequence)
        tpl = PeptideTemplate(name=name,
                            residue_list=residue_list,
                            model=model,
                            sequence=sequence)
        self.db._register_template(tpl)        
    
    def define_protein(self, name, sequence, model):
        """
        Defines a protein template in the pyMBE database.

        Args:
            name ('str'): 
                Unique label that identifies the protein.

            sequence ('str'): 
                Sequence of the protein.

            model ('string'): 
                Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.

        Notes:
            - Currently, only 'lj_setup_mode="wca"' is supported. This corresponds to setting up the WCA potential.
        """
        valid_model_keys = ['1beadAA','2beadAA']
        if model not in valid_model_keys:
            raise ValueError('Invalid key for the protein model, supported models are {valid_model_keys}')
        
        residue_list = self._get_residue_list_from_sequence(sequence=sequence)
        tpl = ProteinTemplate(name=name,
                              model=model,
                              residue_list=residue_list,
                              sequence=sequence)
        self.db._register_template(tpl)
    
    def define_residue(self, name, central_bead, side_chains):
        """
        Defines a residue template in the pyMBE database.

        Args:
            name ('str'): 
                Unique label that identifies the residue.

            central_bead ('str'): 
                'name' of the 'particle' to be placed as central_bead of the residue.

            side_chains('list' of 'str'): 
                List of 'name's of the pmb_objects to be placed as side_chains of the residue. Currently, only pyMBE objects of type 'particle' or 'residue' are supported.
        """
        tpl = ResidueTemplate(name=name,
                              central_bead=central_bead,
                              side_chains=side_chains)
        self.db._register_template(tpl)

    
    def delete_instances_in_system(self, instance_id, pmb_type):
        """
        Deletes the instance with instance_id from the ESPResSo system. 
        Related assembly, molecule, residue, particles and bond instances will also be deleted from the pyMBE dataframe.

        Args:
            instance_id ('int'): 
                id of the assembly to be deleted. 

            pmb_type ('str'): 
                the instance type to be deleted. 

            espresso_system ('espressomd.system.System'): 
                Instance of a system class from espressomd library.
        """
        if pmb_type == "particle":
            instance_identifier = "particle_id"
        elif pmb_type == "residue":
            instance_identifier = "residue_id"
        elif pmb_type in self.db._molecule_like_types:
            instance_identifier = "molecule_id"
        elif pmb_type in self.db._assembly_like_types:
            instance_identifier = "assembly_id"
        particle_ids = self.db._find_instance_ids_by_attribute(pmb_type="particle",
                                                               attribute=instance_identifier,
                                                               value=instance_id)
        self._delete_particles_from_engine(particle_ids=particle_ids)
        self.db.delete_instance(pmb_type=pmb_type,
                                instance_id=instance_id)

    def determine_reservoir_concentrations(self, pH_res, c_salt_res, activity_coefficient_monovalent_pair, max_number_sc_runs=200):
        """
        Determines ionic concentrations in the reservoir at fixed pH and salt concentration.

        Args:
            pH_res ('float'):
                Target pH value in the reservoir.

            c_salt_res ('pint.Quantity'):
                Concentration of monovalent salt (e.g., NaCl) in the reservoir.

            activity_coefficient_monovalent_pair ('callable'):
                Function returning the activity coefficient of a monovalent ion pair
                as a function of ionic strength:
                'gamma = activity_coefficient_monovalent_pair(I)'.

            max_number_sc_runs ('int', optional):
                Maximum number of self-consistent iterations allowed before
                convergence is enforced. Defaults to 200.

        Returns:
            tuple:
                (cH_res, cOH_res, cNa_res, cCl_res)
                - cH_res ('pint.Quantity'): Concentration of H⁺ ions.
                - cOH_res ('pint.Quantity'): Concentration of OH⁻ ions.
                - cNa_res ('pint.Quantity'): Concentration of Na⁺ ions.
                - cCl_res ('pint.Quantity'): Concentration of Cl⁻ ions.

        Notess:
            - The algorithm enforces electroneutrality in the reservoir.
            - Water autodissociation is included via the equilibrium constant 'Kw'.
            - Non-ideal effects enter through activity coefficients depending on
            ionic strength.
            - The implementation follows the self-consistent scheme described in
            Landsgesell (PhD thesis, Sec. 5.3, doi:10.18419/opus-10831), adapted
            from the original code (doi:10.18419/darus-2237).
        """
        cH_res, cOH_res, cNa_res, cCl_res = self.simulation_engine.determine_reservoir_concentrations( pH_res, c_salt_res, activity_coefficient_monovalent_pair, max_number_sc_runs)
        return cH_res, cOH_res, cNa_res, cCl_res

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
        self.simulation_engine.enable_motion_of_rigid_object(instance_id, pmb_type)

    def generate_coordinates_outside_sphere(self, center, radius, max_dist, n_samples):
        """
        Generates random coordinates outside a sphere and inside a larger bounding sphere.

        Args:
            center ('array-like'):
                Coordinates of the center of the spheres.

            radius ('float'):
                Radius of the inner exclusion sphere. Must be positive.

            max_dist ('float'):
                Radius of the outer sampling sphere. Must be larger than 'radius'.

            n_samples ('int'):
                Number of coordinates to generate.

        Returns:
            'list' of 'numpy.ndarray':
                List of coordinates lying outside the inner sphere and inside the
                outer sphere.

        Notess:
            - Points are uniformly sampled inside a sphere of radius 'max_dist' centered at 'center' 
            and only those with a distance greater than or equal to 'radius' from the center are retained.
        """
        if not radius > 0: 
            raise ValueError (f'The value of {radius} must be a positive value')
        if not radius < max_dist:
            raise ValueError(f'The min_dist ({radius} must be lower than the max_dist ({max_dist}))')
        coord_list = []
        counter = 0
        while counter<n_samples:
            coord = self.generate_random_points_in_a_sphere(center=center, 
                                            radius=max_dist,
                                            n_samples=1)[0]
            if np.linalg.norm(coord-np.asarray(center))>=radius:
                coord_list.append (coord)
                counter += 1
        return coord_list
    
    def generate_random_points_in_a_sphere(self, center, radius, n_samples, on_surface=False):
        """
        Generates uniformly distributed random points inside or on the surface of a sphere.

        Args:
            center ('array-like'):
                Coordinates of the center of the sphere.

            radius ('float'):
                Radius of the sphere.

            n_samples ('int'):
                Number of sample points to generate.

            on_surface ('bool', optional):
                If True, points are uniformly sampled on the surface of the sphere.
                If False, points are uniformly sampled within the sphere volume.
                Defaults to False.

        Returns:
            'numpy.ndarray':
                Array of shape '(n_samples, d)' containing the generated coordinates,
                where 'd' is the dimensionality of 'center'.
        Notes:
            - Points are sampled in a space whose dimensionality is inferred 
            from the length of 'center'.
        """
        # initial values
        center=np.array(center)
        d = center.shape[0]
        # sample n_samples points in d dimensions from a standard normal distribution
        samples = self.rng.normal(size=(n_samples, d))
        # make the samples lie on the surface of the unit hypersphere
        normalize_radii = np.linalg.norm(samples, axis=1)[:, np.newaxis]
        samples /= normalize_radii
        if not on_surface:
            # make the samples lie inside the hypersphere with the correct density
            uniform_points = self.rng.uniform(size=n_samples)[:, np.newaxis]
            new_radii = np.power(uniform_points, 1/d)
            samples *= new_radii
        # scale the points to have the correct radius and center
        samples = samples * radius + center
        return samples 

    def generate_trial_perpendicular_vector(self,vector,magnitude):
        """
        Generates a random vector perpendicular to a given vector.

        Args:
            vector ('array-like'):
                Reference vector to which the generated vector will be perpendicular.

            magnitude ('float'):
                Desired magnitude of the perpendicular vector.

        Returns:
            'numpy.ndarray':
                Vector orthogonal to 'vector' with norm equal to 'magnitude'.
        """ 
        np_vec = np.array(vector) 
        if np.all(np_vec == 0):
            raise ValueError('Zero vector')
        np_vec /= np.linalg.norm(np_vec) 
        # Generate a random vector 
        random_vector = self.generate_random_points_in_a_sphere(radius=1, 
                                                                center=[0,0,0],
                                                                n_samples=1, 
                                                                on_surface=True)[0]
        # Project the random vector onto the input vector and subtract the projection
        projection = np.dot(random_vector, np_vec) * np_vec
        perpendicular_vector = random_vector - projection
        # Normalize the perpendicular vector to have the same magnitude as the input vector
        perpendicular_vector /= np.linalg.norm(perpendicular_vector) 
        return perpendicular_vector*magnitude
            
    def get_bond_template(self, particle_name1, particle_name2, use_default_bond=False) :
        """
        Retrieves a bond template connecting two particle templates.

        Args:
            particle_name1 ('str'):
                Name of the first particle template.

            particle_name2 ('str'):
                Name of the second particle template.

            use_default_bond ('bool', optional):
                If True, returns the default bond template when no specific bond
                template is found. Defaults to False.

        Returns:
            'BondTemplate':
                Bond template object retrieved from the pyMBE database.
            
        Notes:
            - This method searches the pyMBE database for a bond template defined between particle templates with names 'particle_name1' and 'particle_name2'. 
            - If no specific bond template is found and 'use_default_bond' is enabled, a default bond template is returned instead.
        """
        # Try to find a specific bond template
        bond_key = BondTemplate.make_bond_key(pn1=particle_name1,
                                              pn2=particle_name2)
        try:
            return self.db.get_template(name=bond_key, 
                                        pmb_type="bond")
        except ValueError:
            pass

        #  Fallback to default bond if allowed
        if use_default_bond:
            return self.db.get_template(name="default", 
                                        pmb_type="bond")

        # No bond template found
        raise ValueError(f"No bond template found between '{particle_name1}' and '{particle_name2}', and default bonds are deactivated.")
    
    def get_charge_number_map(self):
        """
        Construct a mapping from ESPResSo particle types to their charge numbers.

        Returns:
            'dict[int, float]':
                Dictionary mapping ESPResSo particle types to charge numbers,
                ''{es_type: z}''.

        Notess:
            - The mapping is built from particle *states*, not instances.
            - If multiple templates define states with the same ''es_type'',
            the last encountered definition will overwrite previous ones.
            This behavior is intentional and assumes database consistency.
            - Neutral particles (''z = 0'') are included in the map.
        """
        charge_number_map = {}
        particle_templates = self.db.get_templates("particle")
        for tpl in particle_templates.values():
            for state in self.db.get_particle_states_templates(particle_name=tpl.name).values():
                charge_number_map[state.es_type] = state.z
        return charge_number_map

    def get_instances_df(self, pmb_type):
        """
        Returns a dataframe with all instances of type 'pmb_type' in the pyMBE database.

        Args:
            pmb_type ('str'): 
                pmb type to search instances in the pyMBE database.
        
        Returns:
            ('Pandas.Dataframe'): 
                Dataframe with all instances of type 'pmb_type'.
        """
        return self.db._get_instances_df(pmb_type=pmb_type)

    def get_lj_parameters(self, particle_name1, particle_name2, combining_rule='Lorentz-Berthelot'):
        """
        Returns the Lennard-Jones parameters for the interaction between the particle types given by 
        'particle_name1' and 'particle_name2' in the pyMBE database, calculated according to the provided combining rule.

        Args:
            particle_name1 ('str'): 
                label of the type of the first particle type

            particle_name2 ('str'): 
                label of the type of the second particle type

            combining_rule ('string', optional): 
                combining rule used to calculate 'sigma' and 'epsilon' for the potential betwen a pair of particles. Defaults to 'Lorentz-Berthelot'.

        Returns:
            ('dict'):
                {"epsilon": epsilon_value, "sigma": sigma_value, "offset": offset_value, "cutoff": cutoff_value}

        Notes:
            - Currently, the only 'combining_rule' supported is Lorentz-Berthelot.
            - If the sigma value of 'particle_name1' or 'particle_name2' is 0, the function will return an empty dictionary. No LJ interactions are set up for particles with sigma = 0.
        """
        lj_parameters=self.db.get_lj_parameters(particle_name1=particle_name1,particle_name2=particle_name2,combining_rule=combining_rule)
        return lj_parameters

    def get_particle_id_map(self, object_name):
        """
        Collect all particle IDs associated with an object of given name in the
        pyMBE database. 

        Args:
            object_name ('str'): 
                Name of the object.

        Returns:
            ('dict'): 
                {"all": [particle_ids],
                 "residue_map": {residue_id: [particle_ids]},
                 "molecule_map": {molecule_id: [particle_ids]},
                 "assembly_map": {assembly_id: [particle_ids]},}

        Notess:
            - Works for all supported pyMBE templates.
            - Relies in the internal method Manager.get_particle_id_map, see method for the detailed code.
        """
        return self.db.get_particle_id_map(object_name=object_name)

    def get_pka_set(self):
        """
        Retrieve the pKa set for all titratable particles in the pyMBE database.

        Returns:
            ('dict'): 
                Dictionary of the form:
                {"particle_name": {"pka_value": float,
                                   "acidity": "acidic" | "basic"}}
        Notes:
            - If a particle participates in multiple acid/base reactions, an error is raised.
        """
        pka_set = {}
        supported_reactions = ["monoprotic_acid",
                               "monoprotic_base"]
        for reaction in self.db._reactions.values():
            if reaction.reaction_type not in supported_reactions:
                continue
            # Identify involved particle(s)
            particle_names = {participant.particle_name for participant in reaction.participants}
            particle_name = particle_names.pop()
            if particle_name in pka_set:
                raise ValueError(f"Multiple acid/base reactions found for particle '{particle_name}'.")
            pka_set[particle_name] = {"pka_value": reaction.pK}
            if reaction.reaction_type == "monoprotic_acid":
                acidity = "acidic"
            elif reaction.reaction_type == "monoprotic_base":
                acidity = "basic"
            pka_set[particle_name]["acidity"] = acidity
        return pka_set
    
    def get_radius_map(self, dimensionless=True):
        """
        Gets the effective radius of each particle defined in the pyMBE database. 

        Args:
            dimensionless ('bool'):
                If ``True``, return magnitudes expressed in ``reduced_length``.
                If ``False``, return Pint quantities with units.
        
        Returns:
            ('dict'): 
                {espresso_type: radius}.

        Notes:
            - The radius corresponds to (sigma+offset)/2
        """
        return self.db.get_radius_map(dimensionless)

    def get_reactions_df(self):
        """
        Returns a dataframe with all reaction templates in the pyMBE database.

        Returns:
            (Pandas.Dataframe): 
                Dataframe with all  reaction templates.
        """
        return self.db._get_reactions_df()

    def get_reduced_units(self):
        """
        Returns the  current set of reduced units defined in pyMBE.

        Returns:
            reduced_units_text ('str'): 
                text with information about the current set of reduced units.

        """
        unit_length=self.units.Quantity(1,'reduced_length')
        unit_energy=self.units.Quantity(1,'reduced_energy')
        unit_charge=self.units.Quantity(1,'reduced_charge')
        reduced_units_text = "\n".join(["Current set of reduced units:",
                                       f"{unit_length.to('nm'):.5g} = {unit_length}",
                                       f"{unit_energy.to('J'):.5g} = {unit_energy}",
                                       f"{unit_charge.to('C'):.5g} = {unit_charge}",
                                       f"Temperature: {(self.kT/self.kB).to('K'):.5g}"])   
        return reduced_units_text

    def get_templates_df(self, pmb_type):
        """
        Returns a dataframe with all templates of type 'pmb_type' in the pyMBE database.

        Args:
            pmb_type ('str'): 
                pmb type to search templates in the pyMBE database.
        
        Returns:
            ('Pandas.Dataframe'): 
                Dataframe with all templates of type given by 'pmb_type'.
        """
        return self.db._get_templates_df(pmb_type=pmb_type)

    def get_type_map(self):
        """
        Return the mapping of ESPResSo types for all particle states defined in the pyMBE database.
        
        Returns:
            'dict[str, int]':
                A dictionary mapping each particle state to its corresponding ESPResSo type:
                {state_name: es_type, ...}
        """
        
        return self.db.get_es_types_map()

    def initialize_lattice_builder(self, diamond_lattice):
        """
        Initialize the lattice builder with the DiamondLattice object.

        Args:
            diamond_lattice ('DiamondLattice'): 
                DiamondLattice object from the 'lib/lattice' module to be used in the LatticeBuilder.
        """
        from .lib.lattice import LatticeBuilder, DiamondLattice
        if not isinstance(diamond_lattice, DiamondLattice):
            raise TypeError("Currently only DiamondLattice objects are supported.")
        self.lattice_builder = LatticeBuilder(lattice=diamond_lattice)
        logging.info(f"LatticeBuilder initialized with mpc={diamond_lattice.mpc} and box_l={diamond_lattice.box_l}")
        return self.lattice_builder

    def load_database(self, folder, format='csv'):
        """
        Loads a pyMBE database stored in 'folder'.

        Args:
            folder ('str' or 'Path'): 
                Path to the folder where the pyMBE database was stored.

            format ('str', optional): 
                Format of the database to be loaded. Defaults to 'csv'.

        Return:
            ('dict'): 
                metadata with additional information about the source of the information in the database.

        Notes:
            - The folder must contain the files generated by 'pmb.save_database()'.
            - Currently, only 'csv' format is supported.
        """
        supported_formats = ['csv']
        if format not in supported_formats:
            raise ValueError(f"Format {format} not supported. Supported formats are {supported_formats}")
        if format == 'csv':
            metadata =io._load_database_csv(self.db, 
                                            folder=folder)
        return metadata
        
    def load_pka_set(self, filename):
        """
        Load a pKa set and attach chemical states and acid–base reactions
        to existing particle templates.

        Args:
            filename ('str'): 
                Path to a JSON file containing the pKa set. Expected format:
                {"metadata": {...},
                  "data": {"A": {"acidity": "acidic", "pka_value": 4.5},
                           "B": {"acidity": "basic",  "pka_value": 9.8}}}

        Returns:
            ('dict'): 
                Dictionary with bibliographic metadata about the original work were the pKa set was determined.

        Notes:
            - This method is designed for monoprotic acids and bases only.
        """
        with open(filename, "r") as f:
            pka_data = json.load(f)
        pka_set = pka_data["data"]
        metadata = pka_data.get("metadata", {})
        self._check_pka_set(pka_set)
        for particle_name, entry in pka_set.items():
            acidity = entry["acidity"]
            pka = entry["pka_value"]
            self.define_monoprototic_acidbase_reaction(particle_name=particle_name,
                                                       pka=pka,
                                                       acidity=acidity,
                                                       metadata=metadata)
        return metadata
            
    def propose_unused_type(self):
        """
        Propose an unused ESPResSo particle type.

        Returns:
            ('int'): 
                The next available integer ESPResSo type. Returns ''0'' if no integer types are currently defined.
        """    
        return self.db.propose_unused_type()
       
    def read_protein_vtf(self, filename, unit_length=None):
        """
        Loads a coarse-grained protein model from a VTF file.

        Args:
            filename ('str'): 
                Path to the VTF file.
                
            unit_length ('Pint.Quantity'): 
                Unit of length for coordinates (pyMBE UnitRegistry). Defaults to Angstrom.

        Returns:
            ('tuple'):
                ('dict'): Particle topology.    
                ('str'):  One-letter amino-acid sequence (including n/c ends).
        """
        logging.info(f"Loading protein coarse-grain model file: {filename}")
        if unit_length is None:
            unit_length = 1 * self.units.angstrom
        atoms = {}        # atom_id -> atom info
        coords = []       # ordered coordinates
        residues = {}     # resid -> resname (first occurrence)
        has_n_term = False
        has_c_term = False
        aa_3to1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
                   "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G",
                   "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
                   "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                   "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
                   "n": "n", "c": "c"}
        # --- parse VTF ---
        with open(filename, "r") as f:
            for line in f:
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == "atom":
                    atom_id = int(fields[1])
                    atom_name = fields[3]
                    resname = fields[5]
                    resid = int(fields[7])
                    chain_id = fields[9]
                    radius = float(fields[11]) * unit_length
                    atoms[atom_id] = {"name": atom_name,
                                     "resname": resname,
                                     "resid": resid,
                                     "chain_id": chain_id,
                                     "radius": radius}
                    if resname == "n":
                        has_n_term = True
                    elif resname == "c":
                        has_c_term = True
                    # register residue 
                    if resid not in residues:
                        residues[resid] = resname
                elif fields[0].isnumeric():
                    xyz = [(float(x) * unit_length).to("reduced_length").magnitude
                        for x in fields[1:4]]
                    coords.append(xyz)
        sequence = ""
        # N-terminus
        if has_n_term:
            sequence += "n"
        # protein residues only
        protein_resids = sorted(resid for resid, resname in residues.items()  if resname not in ("n", "c", "Ca"))
        for resid in protein_resids:
            resname = residues[resid]
            try:
                sequence += aa_3to1[resname]
            except KeyError:
                raise ValueError(f"Unknown residue name '{resname}' in VTF file")
        # C-terminus
        if has_c_term:
            sequence += "c"
        last_resid = max(protein_resids)
        # --- build topology ---
        topology_dict = {}
        for atom_id in sorted(atoms.keys()):
            atom = atoms[atom_id]
            resname = atom["resname"]
            resid = atom["resid"]
            # apply labeling rules
            if resname == "n":
                label_resid = 0
            elif resname == "c":
                label_resid = last_resid + 1
            elif resname == "Ca":
                label_resid = last_resid + 2
            else:
                label_resid = resid  # preserve original resid 
            label = f"{atom['name']}{label_resid}"
            if label in topology_dict:
                raise ValueError(f"Duplicate particle label '{label}'. Check VTF residue definitions.")
            topology_dict[label] = {"initial_pos": coords[atom_id - 1], "chain_id": atom["chain_id"], "radius": atom["radius"],}
        return topology_dict, sequence

    
    def save_database(self, folder, format='csv'):
        """
        Saves the current pyMBE database into a file 'filename'.

        Args:
            folder ('str' or 'Path'): 
                Path to the folder where the database files will be saved.

        """
        supported_formats = ['csv']
        if format not in supported_formats:
            raise ValueError(f"Format {format} not supported. Supported formats are: {supported_formats}")
        if format == 'csv':
            io._save_database_csv(self.db, 
                                folder=folder)

    def set_particle_initial_state(self, particle_name, state_name):
        """
        Sets the default initial state of a particle template defined in the pyMBE database.

        Args:
            particle_name ('str'): 
                Unique label that identifies the particle template. 

            state_name ('str'): 
                Name of the state to be set as default initial state.
        """
        part_tpl = self.db.get_template(name=particle_name,
        
                                        pmb_type="particle")
        part_tpl.initial_state = state_name
        logging.info(f"Default initial state of particle {particle_name} set to {state_name}.")

    def set_reduced_units(self, unit_length=None, unit_charge=None, temperature=None, Kw=None):
        """
        Sets the set of reduced units used by pyMBE.units and it prints it.

        Args:
            unit_length ('pint.Quantity', optional): 
                Reduced unit of length defined using the 'pmb.units' UnitRegistry. Defaults to None. 

            unit_charge ('pint.Quantity', optional): 
                Reduced unit of charge defined using the 'pmb.units' UnitRegistry. Defaults to None. 

            temperature ('pint.Quantity', optional): 
                Temperature of the system, defined using the 'pmb.units' UnitRegistry. Defaults to None. 

            Kw ('pint.Quantity', optional): 
                Ionic product of water in mol^2/l^2. Defaults to None. 

        Notes:
            - If no 'temperature' is given, a value of 298.15 K is assumed by default.
            - If no 'unit_length' is given, a value of 0.355 nm is assumed by default.
            - If no 'unit_charge' is given, a value of 1 elementary charge is assumed by default. 
            - If no 'Kw' is given, a value of 10^(-14) * mol^2 / l^2 is assumed by default. 
        """
        if unit_length is None:
            unit_length= 0.355*self.units.nm
        if temperature is None:
            temperature = 298.15 * self.units.K
        if unit_charge is None:
            unit_charge = scipy.constants.e * self.units.C
        if Kw is None:
            Kw = 1e-14
        # Sanity check
        variables=[unit_length,temperature,unit_charge]
        dimensionalities=["[length]","[temperature]","[charge]"]
        for variable,dimensionality in zip(variables,dimensionalities):
            self._check_dimensionality(variable,dimensionality)
        self.Kw=Kw*self.units.mol**2 / (self.units.l**2)
        self.kT=temperature*self.kB
        self.units._build_cache()
        self.units.define(f'reduced_energy = {self.kT} ')
        self.units.define(f'reduced_length = {unit_length}')
        self.units.define(f'reduced_charge = {unit_charge}')
        logging.info(self.get_reduced_units())

    def set_simulation_engine(self,simulation_engine,box_l=None):
        """ 
        Sets the instance attribute simulation_engine to an instance of a class of type SimulationEngine.
            
        Args:
            simulation_engine (Any): object which contains the methods to setup molecular dynamics and montecarlo simulations
            box_l('list[float,float,float]'): list of floats with the dimensions of the box
        """

        if isinstance(simulation_engine, espressomd.System):
            self.simulation_engine=EspressoSimulation(box_l=simulation_engine.box_l,
                                                      db=self.db,
                                                      espresso_system=simulation_engine,
                                                      units=self.units,
                                                      kT=self.kT,
                                                      Kw=self.Kw,
                                                      seed=self.seed)
        elif isinstance(simulation_engine,LammpsProtocol): 
            self.simulation_engine=LammpsSimulation(box_l=box_l,
                                                      db=self.db,
                                                      lammps=simulation_engine,
                                                      units=self.units,
                                                      kT=self.kT,
                                                      Kw=self.Kw,
                                                      seed=self.seed)
        else:
            raise ValueError('The specified simulation engine is not implemented yet')

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
       
        RE = self.simulation_engine.setup_cpH(counter_ion=counter_ion, 
                                        constant_pH=constant_pH, 
                                        exclusion_range=exclusion_range, 
                                        use_exclusion_radius_per_type = use_exclusion_radius_per_type)
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
        RE = self.simulation_engine.setup_gcmc(c_salt_res=c_salt_res, 
                                                salt_anion_name=salt_anion_name, 
                                                salt_cation_name=salt_cation_name, 
                                                activity_coefficient=activity_coefficient, 
                                                exclusion_range=exclusion_range, 
                                                use_exclusion_radius_per_type = use_exclusion_radius_per_type)
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
            For the Espresso system class:
            Output ('tuple(reaction_methods.ReactionEnsemble,pint.Quantity)'):

                'reaction_methods.ReactionEnsemble':  
                    espressomd reaction_methods object with all reactions necesary to run the GRxMC ensamble.
                
                'pint.Quantity': 
                    Ionic strength of the reservoir (useful for calculating partition coefficients).

        Notess:
            - This implementation uses the original formulation of the grand-reaction method by Landsgesell et al. [1].

        [1] Landsgesell, J., Hebbeker, P., Rud, O., Lunkad, R., Košovan, P., & Holm, C. (2020). Grand-reaction method for simulations of ionization equilibria coupled to ion partitioning. Macromolecules, 53(8), 3007-3020.
        """     
        output=self.simulation_engine.setup_grxmc_reactions(pH_res=pH_res, 
                                                            c_salt_res=c_salt_res, 
                                                            proton_name=proton_name, 
                                                            hydroxide_name=hydroxide_name, 
                                                            salt_cation_name=salt_cation_name, 
                                                            salt_anion_name=salt_anion_name, 
                                                            activity_coefficient=activity_coefficient, 
                                                            exclusion_range=exclusion_range, 
                                                            use_exclusion_radius_per_type=use_exclusion_radius_per_type)
        
        return output
        
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
            For the Espresso system class:
            Output ('tuple(reaction_methods.ReactionEnsemble,pint.Quantity)'):

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
        output=self.simulation_engine.setup_grxmc_unified(pH_res=pH_res, 
                                                         c_salt_res=c_salt_res, 
                                                         cation_name=cation_name, 
                                                         anion_name=anion_name, 
                                                         activity_coefficient=activity_coefficient, 
                                                         exclusion_range=exclusion_range, 
                                                         use_exclusion_radius_per_type = use_exclusion_radius_per_type)
        return output

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
        self.simulation_engine.setup_lj_interactions(shift_potential=shift_potential, 
                                                     combining_rule=combining_rule)