#
# Copyright (C) 2023-2025 pyMBE-dev team
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

import re
import json
import pint
import numpy as np
import pandas as pd
import scipy.constants
import scipy.optimize
import logging
import importlib.resources
from pyMBE.storage.df_management import _DFManagement as _DFm

class pymbe_library():
    """
    The library for the Molecular Builder for ESPResSo (pyMBE)

    Attributes:
        N_A(`pint.Quantity`): Avogadro number.
        Kb(`pint.Quantity`): Boltzmann constant.
        e(`pint.Quantity`): Elementary charge.
        df(`Pandas.Dataframe`): Dataframe used to bookkeep all the information stored in pyMBE. Typically refered as `pmb.df`. 
        kT(`pint.Quantity`): Thermal energy.
        Kw(`pint.Quantity`): Ionic product of water. Used in the setup of the G-RxMC method.
    """

    def __init__(self, seed, temperature=None, unit_length=None, unit_charge=None, Kw=None):
        """
        Initializes the pymbe_library by setting up the reduced unit system with `temperature` and `reduced_length` 
        and sets up  the `pmb.df` for bookkeeping.

        Args:
            temperature(`pint.Quantity`,optional): Value of the temperature in the pyMBE UnitRegistry. Defaults to None.
            unit_length(`pint.Quantity`, optional): Value of the unit of length in the pyMBE UnitRegistry. Defaults to None.
            unit_charge (`pint.Quantity`,optional): Reduced unit of charge defined using the `pmb.units` UnitRegistry. Defaults to None. 
            Kw (`pint.Quantity`,optional): Ionic product of water in mol^2/l^2. Defaults to None. 
        
        Note:
            - If no `temperature` is given, a value of 298.15 K is assumed by default.
            - If no `unit_length` is given, a value of 0.355 nm is assumed by default.
            - If no `unit_charge` is given, a value of 1 elementary charge is assumed by default. 
            - If no `Kw` is given, a value of 10^(-14) * mol^2 / l^2 is assumed by default. 
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
        
        self.df = _DFm._setup_df()
        self.lattice_builder = None
        self.root = importlib.resources.files(__package__)   

    def _define_particle_entry_in_df(self,name):
        """
        Defines a particle entry in pmb.df.

        Args:
            name(`str`): Unique label that identifies this particle type.

        Returns:
            index(`int`): Index of the particle in pmb.df  
        """

        if  _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
            index = self.df[self.df['name']==name].index[0]                                   
        else:
            index = len(self.df)
            self.df.at [index, 'name'] = name
            self.df.at [index,'pmb_type'] = 'particle'
        self.df.fillna(pd.NA, inplace=True)
        return index

    def _check_supported_molecule(self, molecule_name,valid_pmb_types):
        """
        Checks if the molecule name `molecule_name` is supported by a method of pyMBE.

        Args:
            molecule_name(`str`): pmb object type to be checked.
            valid_pmb_types(`list` of `str`): List of valid pmb types supported by the method.

        Returns:
            pmb_type(`str`): pmb_type of the molecule.
        """
        pmb_type=self.df.loc[self.df['name']==molecule_name].pmb_type.values[0]
        if pmb_type not in valid_pmb_types:
            raise ValueError("The pyMBE object with name {molecule_name} has a pmb_type {pmb_type}. This function only supports pyMBE types {valid_pmb_types}")      
        return pmb_type

    def _check_if_name_has_right_type(self, name, expected_pmb_type, hard_check=True):
        """
        Checks if `name` is of the expected pmb type.

        Args:
            name(`str`): label to check if defined in `pmb.df`.
            expected_pmb_type(`str`): pmb object type corresponding to `name`.
            hard_check(`bool`, optional): If `True`, the raises a ValueError  if `name` is corresponds to an objected defined in the pyMBE DataFrame under a different object type than `expected_pmb_type`.

        Returns:
            `bool`: `True` for success, `False` otherwise.
        """
        pmb_type=self.df.loc[self.df['name']==name].pmb_type.values[0]
        if pmb_type == expected_pmb_type:
            return True
        else:
            if hard_check:
                raise ValueError(f"The name {name} has been defined in the pyMBE DataFrame with a pmb_type = {pmb_type}. This function only supports pyMBE objects with pmb_type = {expected_pmb_type}")
            return False

    def add_bonds_to_espresso(self, espresso_system) :
        """
        Adds all bonds defined in `pmb.df` to `espresso_system`.

        Args:
            espresso_system(`espressomd.system.System`): system object of espressomd library
        """

        if 'bond' in self.df["pmb_type"].values:
            bond_df = self.df.loc[self.df ['pmb_type'] == 'bond']
            bond_list = bond_df.bond_object.values.tolist()
            for bond in bond_list:
                espresso_system.bonded_inter.add(bond)
        else:
            logging.warning('there are no bonds defined in pymbe.df')
        return   
    
    def calculate_center_of_mass_of_molecule(self, molecule_id, espresso_system):
        """
        Calculates the center of the molecule with a given molecule_id.

        Args:
            molecule_id(`int`): Id of the molecule whose center of mass is to be calculated.
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
        
        Returns:
            center_of_mass(`lst`): Coordinates of the center of mass.
        """
        center_of_mass = np.zeros(3)
        axis_list = [0,1,2]
        molecule_name = self.df.loc[(self.df['molecule_id']==molecule_id) & (self.df['pmb_type'].isin(["molecule","protein"]))].name.values[0]
        particle_id_list = self.get_particle_id_map(object_name=molecule_name)["all"]
        for pid in particle_id_list:
            for axis in axis_list:
                center_of_mass [axis] += espresso_system.part.by_id(pid).pos[axis]
        center_of_mass = center_of_mass / len(particle_id_list)
        return center_of_mass

    def calculate_HH(self, molecule_name, pH_list=None, pka_set=None):
        """
        Calculates the charge per molecule according to the ideal Henderson-Hasselbalch titration curve 
        for molecules with the name `molecule_name`.

        Args:
            molecule_name(`str`): name of the molecule to calculate the ideal charge for
            pH_list(`lst`): pH-values to calculate. 
            pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}

        Returns:
            Z_HH(`lst`): Henderson-Hasselbalch prediction of the charge of `sequence` in `pH_list`

        Note:
            - This function supports objects with pmb types: "molecule", "peptide" and "protein".
            - If no `pH_list` is given, 50 equispaced pH-values ranging from 2 to 12 are calculated
            - If no `pka_set` is given, the pKa values are taken from `pmb.df`
            - This function should only be used for single-phase systems. For two-phase systems `pmb.calculate_HH_Donnan`  should be used.
        """
        _DFm._check_if_name_is_defined_in_df(name = molecule_name,
                                             df = self.df)
        self._check_supported_molecule(molecule_name = molecule_name,
                                       valid_pmb_types = ["molecule","peptide","protein"])
        if pH_list is None:
            pH_list=np.linspace(2,12,50)
        if pka_set is None:
            pka_set=self.get_pka_set()
        index = self.df.loc[self.df['name'] == molecule_name].index[0].item() 
        residue_list = self.df.at [index,('residue_list','')].copy()
        particles_in_molecule = []
        for residue in residue_list:
            list_of_particles_in_residue = self.search_particles_in_residue(residue)
            if len(list_of_particles_in_residue) == 0:
                logging.warning(f"The residue {residue} has no particles defined in the pyMBE DataFrame, it will be ignored.")
                continue
            particles_in_molecule += list_of_particles_in_residue
        if len(particles_in_molecule) == 0:
            return [None]*len(pH_list)
        self.check_pka_set(pka_set=pka_set)
        charge_number_map = self.get_charge_number_map()
        Z_HH=[]
        for pH_value in pH_list:    
            Z=0
            for particle in particles_in_molecule:
                if particle in pka_set.keys():
                    if pka_set[particle]['acidity'] == 'acidic':
                        psi=-1
                    elif pka_set[particle]['acidity']== 'basic':
                        psi=+1
                    Z+=psi/(1+10**(psi*(pH_value-pka_set[particle]['pka_value'])))                      
                else:
                    state_one_type = self.df.loc[self.df['name']==particle].state_one.es_type.values[0]
                    Z+=charge_number_map[state_one_type]
            Z_HH.append(Z)
        return Z_HH

    def calculate_HH_Donnan(self, c_macro, c_salt, pH_list=None, pka_set=None):
        """
        Calculates the charge on the different molecules according to the Henderson-Hasselbalch equation coupled to the Donnan partitioning.

        Args:
            c_macro('dict'): {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 
            c_salt('float'): Salt concentration in the reservoir.
            pH_list('lst'): List of pH-values in the reservoir. 
            pka_set('dict'): {"name": {"pka_value": pka, "acidity": acidity}}.

        Returns:
            {"charges_dict": {"name": charges}, "pH_system_list": pH_system_list, "partition_coefficients": partition_coefficients_list}
            pH_system_list ('lst'): List of pH_values in the system.
            partition_coefficients_list ('lst'): List of partition coefficients of cations.

        Note:
            - If no `pH_list` is given, 50 equispaced pH-values ranging from 2 to 12 are calculated
            - If no `pka_set` is given, the pKa values are taken from `pmb.df`
            - If `c_macro` does not contain all charged molecules in the system, this function is likely to provide the wrong result.
        """
        if pH_list is None:
            pH_list=np.linspace(2,12,50)
        if pka_set is None:
            pka_set=self.get_pka_set() 
        self.check_pka_set(pka_set=pka_set)

        partition_coefficients_list = []
        pH_system_list = []
        Z_HH_Donnan={}
        for key in c_macro:
            Z_HH_Donnan[key] = []

        def calc_charges(c_macro, pH):
            """
            Calculates the charges of the different kinds of molecules according to the Henderson-Hasselbalch equation.

            Args:
                c_macro('dic'): {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 
                pH('float'): pH-value that is used in the HH equation.

            Returns:
                charge('dict'): {"molecule_name": charge}
            """
            charge = {}
            for name in c_macro:
                charge[name] = self.calculate_HH(name, [pH], pka_set)[0]
            return charge

        def calc_partition_coefficient(charge, c_macro):
            """
            Calculates the partition coefficients of positive ions according to the ideal Donnan theory.

            Args:
                charge('dict'): {"molecule_name": charge}
                c_macro('dict'): {"name": concentration} - A dict containing the concentrations of all charged macromolecular species in the system. 
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

    def calculate_initial_bond_length(self, bond_object, bond_type, epsilon, sigma, cutoff, offset):
        """
        Calculates the initial bond length that is used when setting up molecules,
        based on the minimum of the sum of bonded and short-range (LJ) interactions.
        
        Args:
            bond_object(`espressomd.interactions.BondedInteractions`): instance of a bond object from espressomd library
            bond_type(`str`): label identifying the used bonded potential
            epsilon(`pint.Quantity`): LJ epsilon of the interaction between the particles
            sigma(`pint.Quantity`): LJ sigma of the interaction between the particles
            cutoff(`pint.Quantity`): cutoff-radius of the LJ interaction 
            offset(`pint.Quantity`): offset of the LJ interaction
        """    
        def truncated_lj_potential(x, epsilon, sigma, cutoff,offset):
            if x>cutoff:
                return 0.0
            else:
                return 4*epsilon*((sigma/(x-offset))**12-(sigma/(x-offset))**6) - 4*epsilon*((sigma/cutoff)**12-(sigma/cutoff)**6)

        epsilon_red=epsilon.to('reduced_energy').magnitude
        sigma_red=sigma.to('reduced_length').magnitude
        cutoff_red=cutoff.to('reduced_length').magnitude
        offset_red=offset.to('reduced_length').magnitude

        if bond_type == "harmonic":
            r_0 = bond_object.params.get('r_0')
            k = bond_object.params.get('k')
            l0 = scipy.optimize.minimize(lambda x: 0.5*k*(x-r_0)**2 + truncated_lj_potential(x, epsilon_red, sigma_red, cutoff_red, offset_red), x0=r_0).x
        elif bond_type == "FENE":
            r_0 = bond_object.params.get('r_0')
            k = bond_object.params.get('k')
            d_r_max = bond_object.params.get('d_r_max')
            l0 = scipy.optimize.minimize(lambda x: -0.5*k*(d_r_max**2)*np.log(1-((x-r_0)/d_r_max)**2) + truncated_lj_potential(x, epsilon_red, sigma_red, cutoff_red,offset_red), x0=1.0).x
        return l0

    def calculate_net_charge(self, espresso_system, molecule_name, dimensionless=False):
        '''
        Calculates the net charge per molecule of molecules with `name` = molecule_name. 
        Returns the net charge per molecule and a maps with the net charge per residue and molecule.

        Args:
            espresso_system(`espressomd.system.System`): system information 
            molecule_name(`str`): name of the molecule to calculate the net charge
            dimensionless(`bool'): sets if the charge is returned with a dimension or not

        Returns:
            (`dict`): {"mean": mean_net_charge, "molecules": {mol_id: net_charge_of_mol_id, }, "residues": {res_id: net_charge_of_res_id, }}

        Note:
            - The net charge of the molecule is averaged over all molecules of type `name` 
            - The net charge of each particle type is averaged over all particle of the same type in all molecules of type `name`
        '''        
        self._check_supported_molecule(molecule_name=molecule_name,
                                        valid_pmb_types=["molecule","protein","peptide"])

        id_map = self.get_particle_id_map(object_name=molecule_name)
        def create_charge_map(espresso_system,id_map,label):
            charge_number_map={}
            for super_id in id_map[label].keys():
                if dimensionless:
                    net_charge=0
                else:
                    net_charge=0 * self.units.Quantity(1,'reduced_charge')
                for pid in id_map[label][super_id]:
                    if dimensionless:
                        net_charge+=espresso_system.part.by_id(pid).q
                    else:
                        net_charge+=espresso_system.part.by_id(pid).q * self.units.Quantity(1,'reduced_charge')
                charge_number_map[super_id]=net_charge
            return charge_number_map
        net_charge_molecules=create_charge_map(label="molecule_map",
                                                espresso_system=espresso_system,
                                                id_map=id_map)
        net_charge_residues=create_charge_map(label="residue_map",
                                                espresso_system=espresso_system,
                                                id_map=id_map)       
        if dimensionless:
            mean_charge=np.mean(np.array(list(net_charge_molecules.values())))
        else:
            mean_charge=np.mean(np.array([value.magnitude for value in net_charge_molecules.values()]))*self.units.Quantity(1,'reduced_charge')
        return {"mean": mean_charge, "molecules": net_charge_molecules, "residues": net_charge_residues}

    def center_molecule_in_simulation_box(self, molecule_id, espresso_system):
        """
        Centers the pmb object matching `molecule_id` in the center of the simulation box in `espresso_md`.
        
        Args:
            molecule_id(`int`): Id of the molecule to be centered.
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
        """
        if len(self.df.loc[self.df['molecule_id']==molecule_id].pmb_type) == 0:
            raise ValueError("The provided molecule_id is not present in the pyMBE dataframe.")      
        center_of_mass = self.calculate_center_of_mass_of_molecule(molecule_id=molecule_id,espresso_system=espresso_system)
        box_center = [espresso_system.box_l[0]/2.0,
                      espresso_system.box_l[1]/2.0,
                      espresso_system.box_l[2]/2.0]
        molecule_name = self.df.loc[(self.df['molecule_id']==molecule_id) & (self.df['pmb_type'].isin(["molecule","protein"]))].name.values[0]
        particle_id_list = self.get_particle_id_map(object_name=molecule_name)["all"]
        for pid in particle_id_list:
            es_pos = espresso_system.part.by_id(pid).pos
            espresso_system.part.by_id(pid).pos = es_pos - center_of_mass + box_center
        return 

    def check_aminoacid_key(self, key):
        """
        Checks if `key` corresponds to a valid aminoacid letter code.

        Args:
            key(`str`): key to be checked.

        Returns:
            `bool`: True if `key` is a valid aminoacid letter code, False otherwise.
        """
        valid_AA_keys=['V', #'VAL'
                       'I', #'ILE'
                       'L', #'LEU'
                       'E', #'GLU'
                       'Q', #'GLN'
                       'D', #'ASP'
                       'N', #'ASN'
                       'H', #'HIS'
                       'W', #'TRP'
                       'F', #'PHE'
                       'Y', #'TYR'
                       'R', #'ARG' 
                       'K', #'LYS'
                       'S', #'SER'
                       'T', #'THR'
                       'M', #'MET'
                       'A', #'ALA'
                       'G', #'GLY'
                       'P', #'PRO'
                       'C'] #'CYS'
        if key in valid_AA_keys:
            return True
        else:
            return False

    def check_dimensionality(self, variable, expected_dimensionality):
        """
        Checks if the dimensionality of `variable` matches `expected_dimensionality`.

        Args:
            variable(`pint.Quantity`): Quantity to be checked.
            expected_dimensionality(`str`): Expected dimension of the variable.

        Returns:
            (`bool`): `True` if the variable if of the expected dimensionality, `False` otherwise.

        Note:
            - `expected_dimensionality` takes dimensionality following the Pint standards [docs](https://pint.readthedocs.io/en/0.10.1/wrapping.html?highlight=dimensionality#checking-dimensionality).
            - For example, to check for a variable corresponding to a velocity `expected_dimensionality = "[length]/[time]"`    
        """
        correct_dimensionality=variable.check(f"{expected_dimensionality}")      
        if not correct_dimensionality:
            raise ValueError(f"The variable {variable} should have a dimensionality of {expected_dimensionality}, instead the variable has a dimensionality of {variable.dimensionality}")
        return correct_dimensionality   

    def check_if_metal_ion(self,key):
        """
        Checks if `key` corresponds to a label of a supported metal ion.

        Args:
            key(`str`): key to be checked

        Returns:
            (`bool`): True if `key`  is a supported metal ion, False otherwise.
        """
        if key in self.get_metal_ions_charge_number_map().keys():
            return True
        else:
            return False

    def check_pka_set(self, pka_set):
        """
        Checks that `pka_set` has the formatting expected by the pyMBE library.
       
        Args:
            pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}
        """
        required_keys=['pka_value','acidity']
        for required_key in required_keys:
            for pka_name, pka_entry in pka_set.items():
                if required_key not in pka_entry:
                    raise ValueError(f'missing a required key "{required_key}" in entry "{pka_name}" of pka_set ("{pka_entry}")')
        return

    def create_added_salt(self, espresso_system, cation_name, anion_name, c_salt):    
        """
        Creates a `c_salt` concentration of `cation_name` and `anion_name` ions into the `espresso_system`.

        Args:
            espresso_system(`espressomd.system.System`): instance of an espresso system object.
            cation_name(`str`): `name` of a particle with a positive charge.
            anion_name(`str`): `name` of a particle with a negative charge.
            c_salt(`float`): Salt concentration.
            
        Returns:
            c_salt_calculated(`float`): Calculated salt concentration added to `espresso_system`.
        """
        for name in [cation_name, anion_name]:
            if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
                logging.warning(f"Object with name '{name}' is not defined in the DataFrame, no ions will be created.")
                return
        self._check_if_name_has_right_type(name=cation_name, 
                                           expected_pmb_type="particle") 
        self._check_if_name_has_right_type(name=anion_name,
                                           expected_pmb_type="particle") 
        cation_name_charge = self.df.loc[self.df['name']==cation_name].state_one.z.values[0]
        anion_name_charge = self.df.loc[self.df['name']==anion_name].state_one.z.values[0]     
        if cation_name_charge <= 0:
            raise ValueError('ERROR cation charge must be positive, charge ',cation_name_charge)
        if anion_name_charge >= 0:
            raise ValueError('ERROR anion charge must be negative, charge ', anion_name_charge)
        # Calculate the number of ions in the simulation box
        volume=self.units.Quantity(espresso_system.volume(), 'reduced_length**3')
        if c_salt.check('[substance] [length]**-3'):
            N_ions= int((volume*c_salt.to('mol/reduced_length**3')*self.N_A).magnitude)
            c_salt_calculated=N_ions/(volume*self.N_A)
        elif c_salt.check('[length]**-3'):
            N_ions= int((volume*c_salt.to('reduced_length**-3')).magnitude)
            c_salt_calculated=N_ions/volume
        else:
            raise ValueError('Unknown units for c_salt, please provided it in [mol / volume] or [particle / volume]', c_salt)
        N_cation = N_ions*abs(anion_name_charge)
        N_anion = N_ions*abs(cation_name_charge)
        self.create_particle(espresso_system=espresso_system, name=cation_name, number_of_particles=N_cation)
        self.create_particle(espresso_system=espresso_system, name=anion_name, number_of_particles=N_anion)
        if c_salt_calculated.check('[substance] [length]**-3'):
            logging.info(f"added salt concentration of {c_salt_calculated.to('mol/L')} given by {N_cation} cations and {N_anion} anions")
        elif c_salt_calculated.check('[length]**-3'):
            logging.info(f"added salt concentration of {c_salt_calculated.to('reduced_length**-3')} given by {N_cation} cations and {N_anion} anions")
        return c_salt_calculated

    def create_bond_in_espresso(self, bond_type, bond_parameters):
        '''
        Creates either a harmonic or a FENE bond in ESPResSo

        Args:
            bond_type(`str`): label to identify the potential to model the bond.
            bond_parameters(`dict`): parameters of the potential of the bond.

        Note:
            Currently, only HARMONIC and FENE bonds are supported.

            For a HARMONIC bond the dictionary must contain:

                - k (`obj`)      : Magnitude of the bond. It should have units of energy/length**2 
                using the `pmb.units` UnitRegistry.
                - r_0 (`obj`)    : Equilibrium bond length. It should have units of length using 
                the `pmb.units` UnitRegistry.
           
            For a FENE bond the dictionary must additionally contain:
                
                - d_r_max (`obj`): Maximal stretching length for FENE. It should have 
                units of length using the `pmb.units` UnitRegistry. Default 'None'.

        Returns:
              bond_object (`obj`): an ESPResSo bond object
        '''
        from espressomd import interactions

        valid_bond_types   = ["harmonic", "FENE"]
        
        if 'k' in bond_parameters:
            bond_magnitude     = bond_parameters['k'].to('reduced_energy / reduced_length**2')
        else:
            raise ValueError("Magnitude of the potential (k) is missing")
        
        if bond_type == 'harmonic':
            if 'r_0' in bond_parameters:
                bond_length        = bond_parameters['r_0'].to('reduced_length')
            else:
                raise ValueError("Equilibrium bond length (r_0) is missing")
            bond_object    = interactions.HarmonicBond(k   = bond_magnitude.magnitude,
                                                       r_0 = bond_length.magnitude)
        elif bond_type == 'FENE':
            if 'r_0' in bond_parameters:
                bond_length        = bond_parameters['r_0'].to('reduced_length').magnitude
            else:
                logging.warning("no value provided for r_0. Defaulting to r_0 = 0")
                bond_length=0
            if 'd_r_max' in bond_parameters:
                max_bond_stret = bond_parameters['d_r_max'].to('reduced_length')
            else:
                raise ValueError("Maximal stretching length (d_r_max) is missing")
            bond_object    = interactions.FeneBond(r_0     = bond_length, 
                                                   k       = bond_magnitude.magnitude,
                                                   d_r_max = max_bond_stret.magnitude)
        else:
            raise NotImplementedError(f"Bond type '{bond_type}' currently not implemented in pyMBE, accepted types are {valid_bond_types}")
        return bond_object


    def create_counterions(self, object_name, cation_name, anion_name, espresso_system):
        """
        Creates particles of `cation_name` and `anion_name` in `espresso_system` to counter the net charge of `pmb_object`.
        
        Args:
            object_name(`str`): `name` of a pymbe object.
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            cation_name(`str`): `name` of a particle with a positive charge.
            anion_name(`str`): `name` of a particle with a negative charge.

        Returns: 
            counterion_number(`dict`): {"name": number}

        Note:
            This function currently does not support the creation of counterions for hydrogels.
        """ 
        for name in [object_name, cation_name, anion_name]:
            if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
                logging.warning(f"Object with name '{name}' is not defined in the DataFrame, no counterions will be created.")
                return
        for name in [cation_name, anion_name]:
            self._check_if_name_has_right_type(name=name, expected_pmb_type="particle")
        self._check_supported_molecule(molecule_name=object_name,
                                        valid_pmb_types=["molecule","peptide","protein"])
        

        cation_charge = self.df.loc[self.df['name']==cation_name].state_one.z.iloc[0]
        anion_charge = self.df.loc[self.df['name']==anion_name].state_one.z.iloc[0]
        object_ids = self.get_particle_id_map(object_name=object_name)["all"]
        counterion_number={}
        object_charge={}
        for name in ['positive', 'negative']:
            object_charge[name]=0
        for id in object_ids:
            if espresso_system.part.by_id(id).q > 0:
                object_charge['positive']+=1*(np.abs(espresso_system.part.by_id(id).q ))
            elif espresso_system.part.by_id(id).q < 0:
                object_charge['negative']+=1*(np.abs(espresso_system.part.by_id(id).q ))
        if object_charge['positive'] % abs(anion_charge) == 0:
            counterion_number[anion_name]=int(object_charge['positive']/abs(anion_charge))
        else:
            raise ValueError('The number of positive charges in the pmb_object must be divisible by the  charge of the anion')
        if object_charge['negative'] % abs(cation_charge) == 0:
            counterion_number[cation_name]=int(object_charge['negative']/cation_charge)
        else:
            raise ValueError('The number of negative charges in the pmb_object must be divisible by the  charge of the cation')
        if counterion_number[cation_name] > 0: 
            self.create_particle(espresso_system=espresso_system, name=cation_name, number_of_particles=counterion_number[cation_name])
        else:
            counterion_number[cation_name]=0
        if counterion_number[anion_name] > 0:
            self.create_particle(espresso_system=espresso_system, name=anion_name, number_of_particles=counterion_number[anion_name])
        else:
            counterion_number[anion_name] = 0
        logging.info('the following counter-ions have been created: ')
        for name in counterion_number.keys():
            logging.info(f'Ion type: {name} created number: {counterion_number[name]}')
        return counterion_number

    def create_hydrogel(self, name, espresso_system):
        """ 
        creates the hydrogel `name` in espresso_system
        Args:
            name(`str`): Label of the hydrogel to be created. `name` must be defined in the `pmb.df`
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.

        Returns:
            hydrogel_info(`dict`):  {"name":hydrogel_name, "chains": {chain_id1: {residue_id1: {'central_bead_id': central_bead_id, 'side_chain_ids': [particle_id1,...]},...,"node_start":node_start,"node_end":node_end}, chain_id2: {...},...}, "nodes":{node1:[node1_id],...}}     
        """
        if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
            logging.warning(f"Hydrogel with name '{name}' is not defined in the DataFrame, no hydrogel will be created.")
            return
        self._check_if_name_has_right_type(name=name, 
                                           expected_pmb_type="hydrogel")
        hydrogel_info={"name":name, "chains":{}, "nodes":{}}
        # placing nodes
        node_positions = {}
        node_topology = self.df[self.df["name"]==name]["node_map"].iloc[0]
        for node_info in node_topology:
            node_index = node_info["lattice_index"]
            node_name = node_info["particle_name"]
            node_pos, node_id = self.create_hydrogel_node(self.format_node(node_index), node_name, espresso_system)
            hydrogel_info["nodes"][self.format_node(node_index)]=node_id
            node_positions[node_id[0]]=node_pos
        
        # Placing chains between nodes
        # Looping over all the 16 chains
        chain_topology = self.df[self.df["name"]==name]["chain_map"].iloc[0]
        for chain_info in chain_topology:
            node_s = chain_info["node_start"]
            node_e = chain_info["node_end"]
            molecule_info = self.create_hydrogel_chain(node_s, node_e, node_positions, espresso_system)
            for molecule_id in molecule_info:
                hydrogel_info["chains"][molecule_id] = molecule_info[molecule_id]
                hydrogel_info["chains"][molecule_id]["node_start"]=node_s
                hydrogel_info["chains"][molecule_id]["node_end"]=node_e
        return hydrogel_info

    def create_hydrogel_chain(self, node_start, node_end, node_positions, espresso_system):
        """
        Creates a chain between two nodes of a hydrogel.

        Args:
            node_start(`str`): name of the starting node particle at which the first residue of the chain will be attached. 
            node_end(`str`): name of the ending node particle at which the last residue of the chain will be attached.
            node_positions(`dict`): dictionary with the positions of the nodes. The keys are the node names and the values are the positions.
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.

        Note:
            For example, if the chain is defined between node_start = ``[0 0 0]`` and node_end = ``[1 1 1]``, the chain will be placed between these two nodes.
            The chain will be placed in the direction of the vector between `node_start` and `node_end`. 
        """
        if self.lattice_builder is None:
            raise ValueError("LatticeBuilder is not initialized. Use `initialize_lattice_builder` first.")

        molecule_name = "chain_"+node_start+"_"+node_end
        sequence = self.df[self.df['name']==molecule_name].residue_list.values [0]
        assert len(sequence) != 0 and not isinstance(sequence, str)
        assert len(sequence) == self.lattice_builder.mpc

        key, reverse = self.lattice_builder._get_node_vector_pair(node_start, node_end)
        assert node_start != node_end or sequence == sequence[::-1], \
            (f"chain cannot be defined between '{node_start}' and '{node_end}' since it "
            "would form a loop with a non-symmetric sequence (under-defined stereocenter)")

        if reverse:
            sequence = sequence[::-1]

        node_start_pos = np.array(list(int(x) for x in node_start.strip('[]').split()))*0.25*self.lattice_builder.box_l
        node_end_pos = np.array(list(int(x) for x in node_end.strip('[]').split()))*0.25*self.lattice_builder.box_l
        node1 = espresso_system.part.select(lambda p: (p.pos == node_start_pos).all()).id
        node2 = espresso_system.part.select(lambda p: (p.pos == node_end_pos).all()).id

        if not node1[0] in node_positions or not node2[0] in node_positions:
            raise ValueError("Set node position before placing a chain between them")
         
        # Finding a backbone vector between node_start and node_end
        vec_between_nodes = np.array(node_positions[node2[0]]) - np.array(node_positions[node1[0]])
        vec_between_nodes = vec_between_nodes - self.lattice_builder.box_l * np.round(vec_between_nodes/self.lattice_builder.box_l)
        backbone_vector = list(vec_between_nodes/(self.lattice_builder.mpc + 1))
        node_start_name = self.df[(self.df["particle_id"]==node1[0]) & (self.df["pmb_type"]=="particle")]["name"].values[0]
        first_res_name = self.df[(self.df["pmb_type"]=="residue") & (self.df["name"]==sequence[0])]["central_bead"].values[0]
        l0 = self.get_bond_length(node_start_name, first_res_name, hard_check=True)
        chain_molecule_info = self.create_molecule(
                name=molecule_name,  # Use the name defined earlier
                number_of_molecules=1,  # Creating one chain
                espresso_system=espresso_system,
                list_of_first_residue_positions=[list(np.array(node_positions[node1[0]]) + np.array(backbone_vector))],#Start at the first node
                backbone_vector=np.array(backbone_vector)/l0,
                use_default_bond=False  # Use defaut bonds between monomers
            )
        # Collecting ids of beads of the chain/molecule
        chain_ids = []
        residue_ids = []
        for molecule_id in chain_molecule_info:
            for residue_id in chain_molecule_info[molecule_id]:
                residue_ids.append(residue_id)
                bead_id = chain_molecule_info[molecule_id][residue_id]['central_bead_id']
                chain_ids.append(bead_id)

        self.lattice_builder.chains[key] = sequence
        # Search bonds between nodes and chain ends
        BeadType_near_to_node_start = self.df[(self.df["residue_id"] == residue_ids[0]) & (self.df["central_bead"].notnull())]["central_bead"].drop_duplicates().iloc[0]
        BeadType_near_to_node_end = self.df[(self.df["residue_id"] == residue_ids[-1]) & (self.df["central_bead"].notnull())]["central_bead"].drop_duplicates().iloc[0]
        bond_node1_first_monomer = self.search_bond(particle_name1 = self.lattice_builder.nodes[node_start],
                                                    particle_name2 = BeadType_near_to_node_start,
                                                    hard_check=False,
                                                    use_default_bond=False)
        bond_node2_last_monomer = self.search_bond(particle_name1 = self.lattice_builder.nodes[node_end],
                                                    particle_name2 = BeadType_near_to_node_end,
                                                    hard_check=False,
                                                    use_default_bond=False)

        espresso_system.part.by_id(node1[0]).add_bond((bond_node1_first_monomer, chain_ids[0]))
        espresso_system.part.by_id(node2[0]).add_bond((bond_node2_last_monomer, chain_ids[-1]))
        # Add bonds to data frame
        self.df, bond_index1 = _DFm._add_bond_in_df(df = self.df,
                                                    particle_id1 = node1[0],
                                                    particle_id2 = chain_ids[0],
                                                    use_default_bond = False)
        _DFm._add_value_to_df(df = self.df,
                              key = ('molecule_id',''),
                              index = int(bond_index1),
                              new_value = molecule_id,
                              overwrite = True)
        _DFm._add_value_to_df(df = self.df,
                              key = ('residue_id',''),
                              index = int(bond_index1),
                              new_value = residue_ids[0],
                              overwrite = True)
        self.df, bond_index2 = _DFm._add_bond_in_df(df = self.df,
                                                    particle_id1 = node2[0],
                                                    particle_id2 = chain_ids[-1],
                                                    use_default_bond = False)
        _DFm._add_value_to_df(df = self.df,
                              key = ('molecule_id',''),
                              index = int(bond_index2),
                              new_value = molecule_id,
                              overwrite = True)
        _DFm._add_value_to_df(df = self.df,
                              key = ('residue_id',''),
                              index = int(bond_index2),
                              new_value = residue_ids[-1],
                              overwrite = True)
        return chain_molecule_info
    
    def create_hydrogel_node(self, node_index, node_name, espresso_system):
        """
        Set a node residue type.
        
        Args:
            node_index(`str`): Lattice node index in the form of a string, e.g. "[0 0 0]".
            node_name(`str`): name of the node particle defined in pyMBE.
        Returns:
            node_position(`list`): Position of the node in the lattice.
            p_id(`int`): Particle ID of the node.
        """
        if self.lattice_builder is None:
            raise ValueError("LatticeBuilder is not initialized. Use `initialize_lattice_builder` first.")

        node_position = np.array(list(int(x) for x in node_index.strip('[]').split()))*0.25*self.lattice_builder.box_l
        p_id = self.create_particle(name = node_name,
                         espresso_system=espresso_system,
                         number_of_particles=1,
                         position = [node_position])
        key = self.lattice_builder._get_node_by_label(node_index)
        self.lattice_builder.nodes[key] = node_name

        return node_position.tolist(), p_id

    def create_molecule(self, name, number_of_molecules, espresso_system, list_of_first_residue_positions=None, backbone_vector=None, use_default_bond=False):
        """
        Creates `number_of_molecules` molecule of type `name` into `espresso_system` and bookkeeps them into `pmb.df`.

        Args:
            name(`str`): Label of the molecule type to be created. `name` must be defined in `pmb.df`
            espresso_system(`espressomd.system.System`): Instance of a system object from espressomd library.
            number_of_molecules(`int`): Number of molecules of type `name` to be created.
            list_of_first_residue_positions(`list`, optional): List of coordinates where the central bead of the first_residue_position will be created, random by default.
            backbone_vector(`list` of `float`): Backbone vector of the molecule, random by default. Central beads of the residues in the `residue_list` are placed along this vector. 
            use_default_bond(`bool`, optional): Controls if a bond of type `default` is used to bond particle with undefined bonds in `pymbe.df`

        Returns:
            molecules_info(`dict`):  {molecule_id: {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids": [particle_id1, ...]}}} 

        Note:
            Despite its name, this function can be used to create both molecules and peptides.    
        """
        if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
            logging.warning(f"Molecule with name '{name}' is not defined in the pyMBE DataFrame, no molecule will be created.")
            return {}
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
        
        # This function works for both molecules and peptides
        if not self._check_if_name_has_right_type(name=name,  expected_pmb_type="molecule", hard_check=False):
            self._check_if_name_has_right_type(name=name, expected_pmb_type="peptide")
        
        # Generate an arbitrary random unit vector
        if backbone_vector is None:
            backbone_vector = self.generate_random_points_in_a_sphere(center=[0,0,0],
                                                    radius=1, 
                                                    n_samples=1,
                                                    on_surface=True)[0]
        else:
            backbone_vector = np.array(backbone_vector)
        first_residue = True
        molecules_info = {}
        residue_list = self.df[self.df['name']==name].residue_list.values [0]
        self.df = _DFm._copy_df_entry(df = self.df,
                                      name = name,
                                      column_name = 'molecule_id',
                                      number_of_copies = number_of_molecules)

        molecules_index = np.where(self.df['name']==name)
        molecule_index_list =list(molecules_index[0])[-number_of_molecules:]
        pos_index = 0 
        for molecule_index in molecule_index_list:        
            molecule_id = _DFm._assign_molecule_id(df = self.df, 
                                                   molecule_index = molecule_index)
            molecules_info[molecule_id] = {}
            for residue in residue_list:
                if first_residue:
                    if list_of_first_residue_positions is None:
                        residue_position = None
                    else:
                        for item in list_of_first_residue_positions:
                            residue_position = [np.array(list_of_first_residue_positions[pos_index])]
                            
                    residues_info = self.create_residue(name=residue,
                                                        espresso_system=espresso_system, 
                                                        central_bead_position=residue_position,  
                                                        use_default_bond= use_default_bond, 
                                                        backbone_vector=backbone_vector)
                    residue_id = next(iter(residues_info))
                    # Add the correct molecule_id to all particles in the residue
                    for index in self.df[self.df['residue_id']==residue_id].index:
                        _DFm._add_value_to_df(df = self.df,
                                              key = ('molecule_id',''),
                                              index = int (index),
                                              new_value = molecule_id,
                                              overwrite = True)
                    central_bead_id = residues_info[residue_id]['central_bead_id']
                    previous_residue = residue
                    residue_position = espresso_system.part.by_id(central_bead_id).pos
                    previous_residue_id = central_bead_id
                    first_residue = False          
                else:                    
                    previous_central_bead_name=self.df[self.df['name']==previous_residue].central_bead.values[0]
                    new_central_bead_name=self.df[self.df['name']==residue].central_bead.values[0]
                    bond = self.search_bond(particle_name1=previous_central_bead_name, 
                                            particle_name2=new_central_bead_name, 
                                            hard_check=True, 
                                            use_default_bond=use_default_bond)                
                    l0 = self.get_bond_length(particle_name1=previous_central_bead_name, 
                                            particle_name2=new_central_bead_name, 
                                            hard_check=True, 
                                            use_default_bond=use_default_bond)                
                    
                    residue_position = residue_position+backbone_vector*l0
                    residues_info = self.create_residue(name=residue, 
                                                        espresso_system=espresso_system, 
                                                        central_bead_position=[residue_position],
                                                        use_default_bond= use_default_bond, 
                                                        backbone_vector=backbone_vector)
                    residue_id = next(iter(residues_info))      
                    for index in self.df[self.df['residue_id']==residue_id].index:
                        _DFm._add_value_to_df(df = self.df,
                                              key = ('molecule_id',''),
                                              index = int(index),
                                              new_value = molecule_id,
                                              overwrite = True)            
                    central_bead_id = residues_info[residue_id]['central_bead_id']
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, previous_residue_id))
                    self.df, bond_index = _DFm._add_bond_in_df(df = self.df,
                                                               particle_id1 = central_bead_id,
                                                               particle_id2 = previous_residue_id,
                                                               use_default_bond = use_default_bond) 
                    _DFm._add_value_to_df(df = self.df,
                                          key = ('molecule_id',''),
                                          index = int(bond_index),
                                          new_value = molecule_id,
                                          overwrite = True)           
                    previous_residue_id = central_bead_id
                    previous_residue = residue                    
                molecules_info[molecule_id][residue_id] = residues_info[residue_id]
            first_residue = True
            pos_index+=1
        
        return molecules_info
    
    def create_particle(self, name, espresso_system, number_of_particles, position=None, fix=False):
        """
        Creates `number_of_particles` particles of type `name` into `espresso_system` and bookkeeps them into `pymbe.df`.
        
        Args:
            name(`str`): Label of the particle type to be created. `name` must be a `particle` defined in `pmb_df`. 
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            number_of_particles(`int`): Number of particles to be created.
            position(list of [`float`,`float`,`float`], optional): Initial positions of the particles. If not given, particles are created in random positions. Defaults to None.
            fix(`bool`, optional): Controls if the particle motion is frozen in the integrator, it is used to create rigid objects. Defaults to False.
        Returns:
            created_pid_list(`list` of `float`): List with the ids of the particles created into `espresso_system`.
        """       
        if number_of_particles <=0:
            return []
        if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
            logging.warning(f"Particle with name '{name}' is not defined in the pyMBE DataFrame, no particle will be created.")
            return []
        self._check_if_name_has_right_type(name=name,
                                           expected_pmb_type="particle")
        # Copy the data of the particle `number_of_particles` times in the `df`
        self.df = _DFm._copy_df_entry(df = self.df,
                                      name = name,
                                      column_name = 'particle_id',
                                      number_of_copies = number_of_particles)
        # Get information from the particle type `name` from the df
        z = self.df.loc[self.df['name'] == name].state_one.z.values[0]
        z = 0. if z is None else z
        es_type = self.df.loc[self.df['name'] == name].state_one.es_type.values[0]
        # Get a list of the index in `df` corresponding to the new particles to be created
        index = np.where(self.df['name'] == name)
        index_list = list(index[0])[-number_of_particles:]
        # Create the new particles into  `espresso_system`
        created_pid_list=[]
        for index in range(number_of_particles):
            df_index = int(index_list[index])
            _DFm._clean_df_row(df = self.df, 
                               index = df_index)
            if position is None:
                particle_position = self.rng.random((1, 3))[0] *np.copy(espresso_system.box_l)
            else:
                particle_position = position[index]
            if len(espresso_system.part.all()) == 0:
                bead_id = 0
            else:
                bead_id = max (espresso_system.part.all().id) + 1
            created_pid_list.append(bead_id)
            kwargs = dict(id=bead_id, pos=particle_position, type=es_type, q=z)
            if fix:
                kwargs["fix"] = 3 * [fix]
            espresso_system.part.add(**kwargs)
            _DFm._add_value_to_df(df = self.df,
                                  key = ('particle_id',''),
                                  index = df_index,
                                  new_value = bead_id)                  
        return created_pid_list

    def create_protein(self, name, number_of_proteins, espresso_system, topology_dict):
        """
        Creates `number_of_proteins` molecules of type `name` into `espresso_system` at the coordinates in `positions`

        Args:
            name(`str`): Label of the protein to be created. 
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            number_of_proteins(`int`): Number of proteins to be created.
            positions(`dict`): {'ResidueNumber': {'initial_pos': [], 'chain_id': ''}}
        """

        if number_of_proteins <=0:
            return
        if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
            logging.warning(f"Protein with name '{name}' is not defined in the pyMBE DataFrame, no protein will be created.")
            return
        self._check_if_name_has_right_type(name=name,
                                           expected_pmb_type="protein")

        self.df = _DFm._copy_df_entry(df = self.df,
                                      name = name,
                                      column_name = 'molecule_id',
                                      number_of_copies = number_of_proteins)
        protein_index = np.where(self.df['name'] == name)
        protein_index_list = list(protein_index[0])[-number_of_proteins:]
        box_half = espresso_system.box_l[0] / 2.0
        for molecule_index in protein_index_list:     
            molecule_id = _DFm._assign_molecule_id(df = self.df,   
                                                   molecule_index = molecule_index)
            protein_center = self.generate_coordinates_outside_sphere(radius = 1, 
                                                                      max_dist = box_half, 
                                                                      n_samples = 1, 
                                                                      center = [box_half]*3)[0]
            for residue in topology_dict.keys():
                residue_name = re.split(r'\d+', residue)[0]
                residue_number = re.split(r'(\d+)', residue)[1]
                residue_position = topology_dict[residue]['initial_pos']
                position = residue_position + protein_center
                particle_id = self.create_particle(name=residue_name,
                                                            espresso_system=espresso_system,
                                                            number_of_particles=1,
                                                            position=[position], 
                                                            fix = True)
                index = self.df[self.df['particle_id']==particle_id[0]].index.values[0]
                _DFm._add_value_to_df(df = self.df,
                                      key = ('residue_id',''),
                                      index = int(index),
                                      new_value = int(residue_number),
                                      overwrite = True)
                _DFm._add_value_to_df(df = self.df,
                                      key = ('molecule_id',''),
                                      index = int(index),
                                      new_value = molecule_id,
                                      overwrite = True)

    def create_residue(self, name, espresso_system, central_bead_position=None,use_default_bond=False, backbone_vector=None):
        """
        Creates a residue of type `name` into `espresso_system` and bookkeeps them into `pmb.df`.

        Args:
            name(`str`): Label of the residue type to be created. `name` must be defined in `pmb.df`
            espresso_system(`espressomd.system.System`): Instance of a system object from espressomd library.
            central_bead_position(`list` of `float`): Position of the central bead.
            use_default_bond(`bool`): Switch to control if a bond of type `default` is used to bond a particle whose bonds types are not defined in `pmb.df`
            backbone_vector(`list` of `float`): Backbone vector of the molecule. All side chains are created perpendicularly to `backbone_vector`.

        Returns:
            residues_info(`dict`): {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids":[particle_id1, ...]}}
        """
        if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
            logging.warning(f"Residue with name '{name}' is not defined in the pyMBE DataFrame, no residue will be created.")
            return
        self._check_if_name_has_right_type(name=name,
                                           expected_pmb_type="residue")
        
        # Copy the data of a residue in the `df
        self.df = _DFm._copy_df_entry(df = self.df,
                                      name = name,
                                      column_name = 'residue_id',
                                      number_of_copies = 1)
        residues_index = np.where(self.df['name']==name)
        residue_index_list =list(residues_index[0])[-1:]
        # search for defined particle and residue names
        particle_and_residue_df = self.df.loc[(self.df['pmb_type']== "particle") | (self.df['pmb_type']== "residue")]
        particle_and_residue_names = particle_and_residue_df["name"].tolist()
        for residue_index in residue_index_list:  
            side_chain_list = self.df.loc[self.df.index[residue_index]].side_chains.values[0]
            for side_chain_element in side_chain_list:
                if side_chain_element not in particle_and_residue_names:              
                    raise ValueError (f"{side_chain_element} is not defined")
        # Internal bookkepping of the residue info (important for side-chain residues)
        # Dict structure {residue_id:{"central_bead_id":central_bead_id, "side_chain_ids":[particle_id1, ...]}}
        residues_info={}
        for residue_index in residue_index_list:     
            _DFm._clean_df_row(df = self.df,
                               index = int(residue_index))
            # Assign a residue_id
            if self.df['residue_id'].isnull().all():
                residue_id=0
            else:
                residue_id = self.df['residue_id'].max() + 1
            _DFm._add_value_to_df(df = self.df,
                                  key = ('residue_id',''),
                                  index = int(residue_index),
                                  new_value = residue_id)
            # create the principal bead   
            central_bead_name = self.df.loc[self.df['name']==name].central_bead.values[0]            
            central_bead_id = self.create_particle(name=central_bead_name,
                                                                espresso_system=espresso_system,
                                                                position=central_bead_position,
                                                                number_of_particles = 1)[0]
            central_bead_position=espresso_system.part.by_id(central_bead_id).pos
            #assigns same residue_id to the central_bead particle created.
            index = self.df[self.df['particle_id']==central_bead_id].index.values[0]
            self.df.at [index,'residue_id'] = residue_id
            # Internal bookkeeping of the central bead id
            residues_info[residue_id]={}
            residues_info[residue_id]['central_bead_id']=central_bead_id
            # create the lateral beads  
            side_chain_list = self.df.loc[self.df.index[residue_index]].side_chains.values[0]
            side_chain_beads_ids = []
            for side_chain_element in side_chain_list:  
                pmb_type = self.df[self.df['name']==side_chain_element].pmb_type.values[0] 
                if pmb_type == 'particle':
                    bond = self.search_bond(particle_name1=central_bead_name, 
                                            particle_name2=side_chain_element, 
                                            hard_check=True, 
                                            use_default_bond=use_default_bond)
                    l0 = self.get_bond_length(particle_name1=central_bead_name, 
                                              particle_name2=side_chain_element, 
                                              hard_check=True, 
                                              use_default_bond=use_default_bond)

                    if backbone_vector is None:
                        bead_position=self.generate_random_points_in_a_sphere(center=central_bead_position, 
                                                                 radius=l0, 
                                                                 n_samples=1,
                                                                 on_surface=True)[0]
                    else:
                        bead_position=central_bead_position+self.generate_trial_perpendicular_vector(vector=np.array(backbone_vector),
                                                                                                    magnitude=l0)
                     
                    side_bead_id = self.create_particle(name=side_chain_element, 
                                                                    espresso_system=espresso_system,
                                                                    position=[bead_position], 
                                                                    number_of_particles=1)[0]
                    index = self.df[self.df['particle_id']==side_bead_id].index.values[0]
                    _DFm._add_value_to_df(df = self.df,
                                          key = ('residue_id',''),
                                          index = int(index),
                                          new_value = residue_id, 
                                          overwrite = True)
                    side_chain_beads_ids.append(side_bead_id)
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, side_bead_id))
                    self.df, index = _DFm._add_bond_in_df(df = self.df,
                                                          particle_id1 = central_bead_id,
                                                          particle_id2 = side_bead_id,
                                                          use_default_bond = use_default_bond)
                    _DFm._add_value_to_df(df = self.df,
                                          key = ('residue_id',''),
                                          index = int(index),
                                          new_value = residue_id, 
                                          overwrite = True)

                elif pmb_type == 'residue':
                    central_bead_side_chain = self.df[self.df['name']==side_chain_element].central_bead.values[0]
                    bond = self.search_bond(particle_name1=central_bead_name, 
                                            particle_name2=central_bead_side_chain, 
                                            hard_check=True, 
                                            use_default_bond=use_default_bond)
                    l0 = self.get_bond_length(particle_name1=central_bead_name, 
                                              particle_name2=central_bead_side_chain, 
                                              hard_check=True, 
                                              use_default_bond=use_default_bond)
                    if backbone_vector is None:
                        residue_position=self.generate_random_points_in_a_sphere(center=central_bead_position, 
                                                                    radius=l0, 
                                                                    n_samples=1,
                                                                    on_surface=True)[0]
                    else:
                        residue_position=central_bead_position+self.generate_trial_perpendicular_vector(vector=backbone_vector,
                                                                                                        magnitude=l0)
                    lateral_residue_info = self.create_residue(name=side_chain_element, 
                                                                espresso_system=espresso_system,
                                                                central_bead_position=[residue_position],
                                                                use_default_bond=use_default_bond)
                    lateral_residue_dict=list(lateral_residue_info.values())[0]
                    central_bead_side_chain_id=lateral_residue_dict['central_bead_id']
                    lateral_beads_side_chain_ids=lateral_residue_dict['side_chain_ids']
                    residue_id_side_chain=list(lateral_residue_info.keys())[0]
                    # Change the residue_id of the residue in the side chain to the one of the bigger residue
                    index = self.df[(self.df['residue_id']==residue_id_side_chain) & (self.df['pmb_type']=='residue') ].index.values[0]
                    _DFm._add_value_to_df(df = self.df,
                                          key = ('residue_id',''),
                                          index = int(index),
                                          new_value = residue_id, 
                                          overwrite = True)
                    # Change the residue_id of the particles in the residue in the side chain
                    side_chain_beads_ids+=[central_bead_side_chain_id]+lateral_beads_side_chain_ids
                    for particle_id in side_chain_beads_ids:
                        index = self.df[(self.df['particle_id']==particle_id) & (self.df['pmb_type']=='particle')].index.values[0]
                        _DFm._add_value_to_df(df = self.df,
                                              key = ('residue_id',''),
                                              index = int (index),
                                              new_value = residue_id, 
                                              overwrite = True)
                    espresso_system.part.by_id(central_bead_id).add_bond((bond, central_bead_side_chain_id))
                    self.df, index = _DFm._add_bond_in_df(df = self.df,
                                                          particle_id1 = central_bead_id,
                                                          particle_id2 = central_bead_side_chain_id,
                                                          use_default_bond = use_default_bond)
                    _DFm._add_value_to_df(df = self.df,
                                          key = ('residue_id',''),
                                          index = int(index),
                                          new_value = residue_id, 
                                          overwrite = True)        
                    # Change the residue_id of the bonds in the residues in the side chain to the one of the bigger residue
                    for index in self.df[(self.df['residue_id']==residue_id_side_chain) & (self.df['pmb_type']=='bond') ].index:        
                        _DFm._add_value_to_df(df = self.df,
                                              key = ('residue_id',''),
                                              index = int(index),
                                              new_value = residue_id, 
                                              overwrite = True)
            # Internal bookkeeping of the side chain beads ids
            residues_info[residue_id]['side_chain_ids']=side_chain_beads_ids
        return  residues_info

    

    def define_AA_residues(self, sequence, model):
        """
        Defines in `pmb.df` all the different residues in `sequence`.

        Args:
            sequence(`lst`):  Sequence of the peptide or protein.
            model(`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.

        Returns:
            residue_list(`list` of `str`): List of the `name`s of the `residue`s  in the sequence of the `molecule`.             
        """
        residue_list = []
        for residue_name in sequence:
            if model == '1beadAA':
                central_bead = residue_name
                side_chains = []
            elif model == '2beadAA':
                if residue_name in ['c','n', 'G']: 
                    central_bead = residue_name
                    side_chains = []
                else:
                    central_bead = 'CA'              
                    side_chains = [residue_name]
            if residue_name not in residue_list:   
                self.define_residue(name = 'AA-'+residue_name, 
                                    central_bead = central_bead,
                                    side_chains = side_chains)              
            residue_list.append('AA-'+residue_name)
        return residue_list

    def define_bond(self, bond_type, bond_parameters, particle_pairs):
        
        '''
        Defines a pmb object of type `bond` in `pymbe.df`.

        Args:
            bond_type(`str`): label to identify the potential to model the bond.
            bond_parameters(`dict`): parameters of the potential of the bond.
            particle_pairs(`lst`): list of the `names` of the `particles` to be bonded.

        Note:
            Currently, only HARMONIC and FENE bonds are supported.

            For a HARMONIC bond the dictionary must contain the following parameters:

                - k (`obj`)      : Magnitude of the bond. It should have units of energy/length**2 
                using the `pmb.units` UnitRegistry.
                - r_0 (`obj`)    : Equilibrium bond length. It should have units of length using 
                the `pmb.units` UnitRegistry.
           
            For a FENE bond the dictionary must contain the same parameters as for a HARMONIC bond and:
                
                - d_r_max (`obj`): Maximal stretching length for FENE. It should have 
                units of length using the `pmb.units` UnitRegistry. Default 'None'.
        '''

        bond_object=self.create_bond_in_espresso(bond_type, bond_parameters)
        for particle_name1, particle_name2 in particle_pairs:

            lj_parameters=self.get_lj_parameters(particle_name1 = particle_name1,
                                                 particle_name2 = particle_name2,
                                                 combining_rule = 'Lorentz-Berthelot')
      
            l0 = self.calculate_initial_bond_length(bond_object = bond_object,
                                                    bond_type   = bond_type,
                                                    epsilon     = lj_parameters["epsilon"],
                                                    sigma       = lj_parameters["sigma"],
                                                    cutoff      = lj_parameters["cutoff"],
                                                    offset      = lj_parameters["offset"],)
            index = len(self.df)
            for label in [f'{particle_name1}-{particle_name2}', f'{particle_name2}-{particle_name1}']:
                _DFm._check_if_multiple_pmb_types_for_name(name=label, 
                                                           pmb_type_to_be_defined="bond",
                                                           df=self.df)
            name=f'{particle_name1}-{particle_name2}'
            _DFm._check_if_multiple_pmb_types_for_name(name=name, 
                                                       pmb_type_to_be_defined="bond", 
                                                       df=self.df)
            self.df.at [index,'name']= name
            self.df.at [index,'bond_object'] = bond_object
            self.df.at [index,'l0'] = l0
            _DFm._add_value_to_df(df = self.df,
                                  index = index,
                                  key = ('pmb_type',''),
                                  new_value = 'bond')
            _DFm._add_value_to_df(df = self.df,
                                  index = index,
                                  key = ('parameters_of_the_potential',''),
                                  new_value = bond_object.get_params(),
                                  non_standard_value = True)
        self.df.fillna(pd.NA, inplace=True)
        return

    
    def define_default_bond(self, bond_type, bond_parameters, epsilon=None, sigma=None, cutoff=None, offset=None):
        """
        Asigns `bond` in `pmb.df` as the default bond.
        The LJ parameters can be optionally provided to calculate the initial bond length

        Args:
            bond_type(`str`): label to identify the potential to model the bond.
            bond_parameters(`dict`): parameters of the potential of the bond.
            sigma(`float`, optional): LJ sigma of the interaction between the particles.
            epsilon(`float`, optional): LJ epsilon for the interaction between the particles.
            cutoff(`float`, optional): cutoff-radius of the LJ interaction.
            offset(`float`, optional): offset of the LJ interaction.

        Note:
            - Currently, only harmonic and FENE bonds are supported. 
        """
        
        bond_object=self.create_bond_in_espresso(bond_type, bond_parameters)

        if epsilon is None:
            epsilon=1*self.units('reduced_energy')
        if sigma is None:
            sigma=1*self.units('reduced_length')
        if cutoff is None:
            cutoff=2**(1.0/6.0)*self.units('reduced_length')
        if offset is None:
            offset=0*self.units('reduced_length')
        l0 = self.calculate_initial_bond_length(bond_object = bond_object, 
                                                bond_type   = bond_type, 
                                                epsilon     = epsilon,
                                                sigma       = sigma,
                                                cutoff      = cutoff,
                                                offset      = offset)

        _DFm._check_if_multiple_pmb_types_for_name(name='default',
                                                   pmb_type_to_be_defined='bond', 
                                                   df=self.df)

        index = max(self.df.index, default=-1) + 1
        self.df.at [index,'name']        = 'default'
        self.df.at [index,'bond_object'] = bond_object
        self.df.at [index,'l0']          = l0
        _DFm._add_value_to_df(df = self.df,
                              index = index,
                              key = ('pmb_type',''),
                              new_value = 'bond')
        _DFm._add_value_to_df(df = self.df,
                              index = index,
                              key = ('parameters_of_the_potential',''),
                              new_value = bond_object.get_params(),
                              non_standard_value=True)
        self.df.fillna(pd.NA, inplace=True)
        return
    
    def define_hydrogel(self, name, node_map, chain_map):
        """
        Defines a pyMBE object of type `hydrogel` in `pymbe.df`.

        Args:
            name(`str`): Unique label that identifies the `hydrogel`.
            node_map(`list of ict`): [{"particle_name": , "lattice_index": }, ... ]
            chain_map(`list of dict`): [{"node_start": , "node_end": , "residue_list": , ... ]
        """
        node_indices = {tuple(entry['lattice_index']) for entry in node_map}
        diamond_indices = {tuple(row) for row in self.lattice_builder.lattice.indices}
        if node_indices != diamond_indices:
            raise ValueError(f"Incomplete hydrogel: A diamond lattice must contain exactly 8 lattice indices, {diamond_indices} ")
        
        chain_map_connectivity = set()
        for entry in chain_map:
            start = self.lattice_builder.node_labels[entry['node_start']]
            end = self.lattice_builder.node_labels[entry['node_end']]
            chain_map_connectivity.add((start,end))

        if self.lattice_builder.lattice.connectivity != chain_map_connectivity:
            raise ValueError("Incomplete hydrogel: A diamond lattice must contain correct 16 lattice index pairs")

        _DFm._check_if_multiple_pmb_types_for_name(name=name,
                                                   pmb_type_to_be_defined='hydrogel',
                                                   df=self.df)

        index = len(self.df)
        self.df.at [index, "name"] = name
        self.df.at [index, "pmb_type"] = "hydrogel"
        _DFm._add_value_to_df(df = self.df,
                              index = index,
                              key = ('node_map',''),
                              new_value = node_map,
                              non_standard_value = True)
        _DFm._add_value_to_df(df = self.df,
                              index = index,
                              key = ('chain_map',''),
                              new_value = chain_map,
                              non_standard_value = True)
        for chain_label in chain_map:
            node_start = chain_label["node_start"]
            node_end = chain_label["node_end"]
            residue_list = chain_label['residue_list']
            # Molecule name
            molecule_name = "chain_"+node_start+"_"+node_end
            self.define_molecule(name=molecule_name, residue_list=residue_list)
        return

    def define_molecule(self, name, residue_list):
        """
        Defines a pyMBE object of type `molecule` in `pymbe.df`.

        Args:
            name(`str`): Unique label that identifies the `molecule`.
            residue_list(`list` of `str`): List of the `name`s of the `residue`s  in the sequence of the `molecule`.  
        """
        _DFm._check_if_multiple_pmb_types_for_name(name=name,
                                                   pmb_type_to_be_defined='molecule',
                                                   df=self.df)

        index = len(self.df)
        self.df.at [index,'name'] = name
        self.df.at [index,'pmb_type'] = 'molecule'
        self.df.at [index,('residue_list','')] = residue_list
        self.df.fillna(pd.NA, inplace=True)
        return   

    def define_particle(self, name, z=0, acidity=pd.NA, pka=pd.NA, sigma=pd.NA, epsilon=pd.NA, cutoff=pd.NA, offset=pd.NA,overwrite=False):
        """
        Defines the properties of a particle object.

        Args:
            name(`str`): Unique label that identifies this particle type.  
            z(`int`, optional): Permanent charge number of this particle type. Defaults to 0.
            acidity(`str`, optional): Identifies whether if the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to pd.NA.
            pka(`float`, optional): If `particle` is an acid or a base, it defines its  pka-value. Defaults to pd.NA.
            sigma(`pint.Quantity`, optional): Sigma parameter used to set up Lennard-Jones interactions for this particle type. Defaults to pd.NA.
            cutoff(`pint.Quantity`, optional): Cutoff parameter used to set up Lennard-Jones interactions for this particle type. Defaults to pd.NA.
            offset(`pint.Quantity`, optional): Offset parameter used to set up Lennard-Jones interactions for this particle type. Defaults to pd.NA.
            epsilon(`pint.Quantity`, optional): Epsilon parameter used to setup Lennard-Jones interactions for this particle tipe. Defaults to pd.NA.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False.

        Note:
            - `sigma`, `cutoff` and `offset` must have a dimensitonality of `[length]` and should be defined using pmb.units.
            - `epsilon` must have a dimensitonality of `[energy]` and should be defined using pmb.units.
            - `cutoff` defaults to `2**(1./6.) reduced_length`. 
            - `offset` defaults to 0.
            - The default setup corresponds to the Weeks−Chandler−Andersen (WCA) model, corresponding to purely steric interactions.
            - For more information on `sigma`, `epsilon`, `cutoff` and `offset` check `pmb.setup_lj_interactions()`.
        """ 
        index=self._define_particle_entry_in_df(name=name)
        _DFm._check_if_multiple_pmb_types_for_name(name=name,
                                                   pmb_type_to_be_defined='particle',
                                                   df=self.df)

        # If `cutoff` and `offset` are not defined, default them to the following values
        if pd.isna(cutoff):
            cutoff=self.units.Quantity(2**(1./6.), "reduced_length")
        if pd.isna(offset):
            offset=self.units.Quantity(0, "reduced_length")

        # Define LJ parameters
        parameters_with_dimensionality={"sigma":{"value": sigma, "dimensionality": "[length]"},
                                        "cutoff":{"value": cutoff, "dimensionality": "[length]"},
                                        "offset":{"value": offset, "dimensionality": "[length]"},
                                        "epsilon":{"value": epsilon, "dimensionality": "[energy]"},}

        for parameter_key in parameters_with_dimensionality.keys():
            if not pd.isna(parameters_with_dimensionality[parameter_key]["value"]):
                self.check_dimensionality(variable=parameters_with_dimensionality[parameter_key]["value"], 
                                          expected_dimensionality=parameters_with_dimensionality[parameter_key]["dimensionality"])
                _DFm._add_value_to_df(df = self.df,
                                      key = (parameter_key,''),
                                      index = index,
                                      new_value = parameters_with_dimensionality[parameter_key]["value"],
                                      overwrite = overwrite)

        # Define particle acid/base properties
        self.set_particle_acidity(name=name, 
                                acidity=acidity, 
                                default_charge_number=z, 
                                pka=pka,
                                overwrite=overwrite)
        self.df.fillna(pd.NA, inplace=True)
        return 
    
    def define_particles(self, parameters, overwrite=False):
        '''
        Defines a particle object in pyMBE for each particle name in `particle_names`

        Args:
            parameters(`dict`):  dictionary with the particle parameters. 
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False. 

        Note:
            - parameters = {"particle_name1: {"sigma": sigma_value, "epsilon": epsilon_value, ...}, particle_name2: {...},}
        '''
        if not parameters:
            return 0
        for particle_name in parameters.keys():
            parameters[particle_name]["overwrite"]=overwrite
            self.define_particle(**parameters[particle_name])
        return
    
    def define_peptide(self, name, sequence, model):
        """
        Defines a pyMBE object of type `peptide` in the `pymbe.df`.

        Args:
            name (`str`): Unique label that identifies the `peptide`.
            sequence (`string`): Sequence of the `peptide`.
            model (`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.
        """
        _DFm._check_if_multiple_pmb_types_for_name(name = name, 
                                                   pmb_type_to_be_defined='peptide',
                                                   df=self.df)

        valid_keys = ['1beadAA','2beadAA']
        if model not in valid_keys:
            raise ValueError('Invalid label for the peptide model, please choose between 1beadAA or 2beadAA')
        clean_sequence = self.protein_sequence_parser(sequence=sequence)    
        residue_list = self.define_AA_residues(sequence=clean_sequence,
                                                model=model)
        self.define_molecule(name = name, residue_list=residue_list)
        index = self.df.loc[self.df['name'] == name].index.item() 
        self.df.at [index,'model'] = model
        self.df.at [index,('sequence','')] = clean_sequence
        self.df.at [index,'pmb_type'] = "peptide"
        self.df.fillna(pd.NA, inplace=True)
        
    
    def define_protein(self, name,model, topology_dict, lj_setup_mode="wca", overwrite=False):
        """
        Defines a globular protein pyMBE object  in `pymbe.df`.

        Args:
            name (`str`): Unique label that identifies the protein.
            model (`string`): Model name. Currently only models with 1 bead '1beadAA' or with 2 beads '2beadAA' per amino acid are supported.
            topology_dict (`dict`): {'initial_pos': coords_list, 'chain_id': id, 'radius': radius_value}
            lj_setup_mode(`str`): Key for the setup for the LJ potential. Defaults to "wca".
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False. 

        Note:
            - Currently, only `lj_setup_mode="wca"` is supported. This corresponds to setting up the WCA potential.
        """

        # Sanity checks
        _DFm._check_if_multiple_pmb_types_for_name(name = name,
                                                   pmb_type_to_be_defined='protein',
                                                   df=self.df)
        valid_model_keys = ['1beadAA','2beadAA']
        valid_lj_setups = ["wca"]
        if model not in valid_model_keys:
            raise ValueError('Invalid key for the protein model, supported models are {valid_model_keys}')
        if lj_setup_mode not in valid_lj_setups:
            raise ValueError('Invalid key for the lj setup, supported setup modes are {valid_lj_setups}')
        if lj_setup_mode == "wca":
            sigma = 1*self.units.Quantity("reduced_length")
            epsilon = 1*self.units.Quantity("reduced_energy")
        part_dict={}
        sequence=[]
        metal_ions_charge_number_map=self.get_metal_ions_charge_number_map()
        for particle in topology_dict.keys():
            particle_name = re.split(r'\d+', particle)[0] 
            if particle_name not in part_dict.keys():
                if lj_setup_mode == "wca":
                    part_dict[particle_name]={"sigma": sigma,
                                        "offset": topology_dict[particle]['radius']*2-sigma,
                                        "epsilon": epsilon,
                                        "name": particle_name}
                if self.check_if_metal_ion(key=particle_name):
                    z=metal_ions_charge_number_map[particle_name]
                else:
                    z=0
                part_dict[particle_name]["z"]=z
            
            if self.check_aminoacid_key(key=particle_name):
                sequence.append(particle_name) 
            
        self.define_particles(parameters=part_dict,
                            overwrite=overwrite)
        residue_list = self.define_AA_residues(sequence=sequence, 
                                               model=model)
        index = len(self.df)
        self.df.at [index,'name'] = name
        self.df.at [index,'pmb_type'] = 'protein'
        self.df.at [index,'model'] = model
        self.df.at [index,('sequence','')] = sequence  
        self.df.at [index,('residue_list','')] = residue_list
        self.df.fillna(pd.NA, inplace=True)    
        return 
    
    def define_residue(self, name, central_bead, side_chains):
        """
        Defines a pyMBE object of type `residue` in `pymbe.df`.

        Args:
            name(`str`): Unique label that identifies the `residue`.
            central_bead(`str`): `name` of the `particle` to be placed as central_bead of the `residue`.
            side_chains(`list` of `str`): List of `name`s of the pmb_objects to be placed as side_chains of the `residue`. Currently, only pmb_objects of type `particle`s or `residue`s are supported.
        """
        _DFm._check_if_multiple_pmb_types_for_name(name=name,
                                                   pmb_type_to_be_defined='residue',
                                                   df=self.df)

        index = len(self.df)
        self.df.at [index, 'name'] = name
        self.df.at [index,'pmb_type'] = 'residue'
        self.df.at [index,'central_bead'] = central_bead
        self.df.at [index,('side_chains','')] = side_chains
        self.df.fillna(pd.NA, inplace=True)
        return    

    def delete_molecule_in_system(self, molecule_id, espresso_system):
        """
        Deletes the molecule with `molecule_id` from the `espresso_system`, including all particles and residues associated with that particles.
        The ids of the molecule, particle and residues deleted are also cleaned from `pmb.df`

        Args:
            molecule_id(`int`): id of the molecule to be deleted. 
            espresso_system(`espressomd.system.System`): Instance of a system class from espressomd library.

        """
        # Sanity checks 
        id_mask = (self.df['molecule_id'] == molecule_id) & (self.df['pmb_type'].isin(["molecule", "peptide"]))
        molecule_row = self.df.loc[id_mask]
        if molecule_row.empty:
            raise ValueError(f"No molecule found with molecule_id={molecule_id} in the DataFrame.")
        # Clean molecule from pmb.df
        self.df = _DFm._clean_ids_in_df_row(df  = self.df, 
                                            row = molecule_row)
        # Delete particles and residues in the molecule
        residue_mask = (self.df['molecule_id'] == molecule_id) & (self.df['pmb_type'] == "residue")
        residue_rows = self.df.loc[residue_mask]
        residue_ids = set(residue_rows["residue_id"].values)
        for residue_id in residue_ids:
            self.delete_residue_in_system(residue_id=residue_id,
                                           espresso_system=espresso_system)
        
        # Clean deleted backbone bonds from pmb.df
        bond_mask = (self.df['molecule_id'] == molecule_id) & (self.df['pmb_type'] == "bond")
        number_of_bonds = len(self.df.loc[bond_mask])
        for _ in range(number_of_bonds):
            bond_mask = (self.df['molecule_id'] == molecule_id) & (self.df['pmb_type'] == "bond")
            bond_rows = self.df.loc[bond_mask]
            row = bond_rows.loc[[bond_rows.index[0]]]
            self.df = _DFm._clean_ids_in_df_row(df = self.df, 
                                                row = row)

    def delete_particle_in_system(self, particle_id, espresso_system):
        """
        Deletes the particle with `particle_id` from the `espresso_system`.
        The particle ids of the particle and residues deleted are also cleaned from `pmb.df`

        Args:
            particle_id(`int`): id of the molecule to be deleted. 
            espresso_system(`espressomd.system.System`): Instance of a system class from espressomd library.

        """
        # Sanity check if there is a particle with the input particle id
        id_mask = (self.df['particle_id'] == particle_id) & (self.df['pmb_type'] == "particle")
        particle_row = self.df.loc[id_mask]
        if particle_row.empty:
            raise ValueError(f"No particle found with particle_id={particle_id} in the DataFrame.")
        espresso_system.part.by_id(particle_id).remove()
        self.df = _DFm._clean_ids_in_df_row(df = self.df, 
                                            row = particle_row)

    def delete_residue_in_system(self, residue_id, espresso_system):
        """
        Deletes the residue with `residue_id`, and the particles associated with it from the `espresso_system`.
        The ids of the residue and particles deleted are also cleaned from `pmb.df`

        Args:
            residue_id(`int`): id of the residue to be deleted. 
            espresso_system(`espressomd.system.System`): Instance of a system class from espressomd library.
        """
        # Sanity check if there is a residue with the input residue id
        id_mask = (self.df['residue_id'] == residue_id) & (self.df['pmb_type'] == "residue")
        residue_row = self.df.loc[id_mask]
        if residue_row.empty:
            raise ValueError(f"No residue found with residue_id={residue_id} in the DataFrame.")
        residue_map=self.get_particle_id_map(object_name=residue_row["name"].values[0])["residue_map"]
        particle_ids = residue_map[residue_id]
        # Clean residue from pmb.df
        self.df = _DFm._clean_ids_in_df_row(df = self.df, 
                                            row = residue_row)
        # Delete particles in the residue
        for particle_id in particle_ids:
            self.delete_particle_in_system(particle_id=particle_id,
                                           espresso_system=espresso_system)
        # Clean deleted bonds from pmb.df
        bond_mask = (self.df['residue_id'] == residue_id) & (self.df['pmb_type'] == "bond")
        number_of_bonds = len(self.df.loc[bond_mask])
        for _ in range(number_of_bonds):
            bond_mask = (self.df['residue_id'] == residue_id) & (self.df['pmb_type'] == "bond")
            bond_rows = self.df.loc[bond_mask]
            row = bond_rows.loc[[bond_rows.index[0]]]
            self.df = _DFm._clean_ids_in_df_row(df = self.df, 
                                                row = row)

    def determine_reservoir_concentrations(self, pH_res, c_salt_res, activity_coefficient_monovalent_pair, max_number_sc_runs=200):
        """
        Determines the concentrations of the various species in the reservoir for given values of the pH and salt concentration.
        To do this, a system of nonlinear equations involving the pH, the ionic product of water, the activity coefficient of an
        ion pair and the concentrations of the various species is solved numerically using a self-consistent approach.
        More details can be found in chapter 5.3 of Landsgesell (doi.org/10.18419/opus-10831).
        This is a modified version of the code by Landsgesell et al. (doi.org/10.18419/darus-2237).

        Args:
            pH_res('float'): pH-value in the reservoir.
            c_salt_res('pint.Quantity'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            activity_coefficient_monovalent_pair('callable', optional): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.

        Returns:
            cH_res('pint.Quantity'): Concentration of H+ ions.
            cOH_res('pint.Quantity'): Concentration of OH- ions.
            cNa_res('pint.Quantity'): Concentration of Na+ ions.
            cCl_res('pint.Quantity'): Concentration of Cl- ions.
        """

        self_consistent_run = 0
        cH_res = 10**(-pH_res) * self.units.mol/self.units.l #initial guess for the concentration of H+

        def determine_reservoir_concentrations_selfconsistently(cH_res, c_salt_res):
            #Calculate and initial guess for the concentrations of various species based on ideal gas estimate
            cOH_res = self.Kw / cH_res 
            cNa_res = None
            cCl_res = None
            if cOH_res>=cH_res:
                #adjust the concentration of sodium if there is excess OH- in the reservoir:
                cNa_res = c_salt_res + (cOH_res-cH_res)
                cCl_res = c_salt_res
            else:
                # adjust the concentration of chloride if there is excess H+ in the reservoir
                cCl_res = c_salt_res + (cH_res-cOH_res)
                cNa_res = c_salt_res
                
            def calculate_concentrations_self_consistently(cH_res, cOH_res, cNa_res, cCl_res):
                nonlocal max_number_sc_runs, self_consistent_run
                if self_consistent_run<max_number_sc_runs:
                    self_consistent_run+=1
                    ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
                    cOH_res = self.Kw / (cH_res * activity_coefficient_monovalent_pair(ionic_strength_res))
                    if cOH_res>=cH_res:
                        #adjust the concentration of sodium if there is excess OH- in the reservoir:
                        cNa_res = c_salt_res + (cOH_res-cH_res)
                        cCl_res = c_salt_res
                    else:
                        # adjust the concentration of chloride if there is excess H+ in the reservoir
                        cCl_res = c_salt_res + (cH_res-cOH_res)
                        cNa_res = c_salt_res
                    return calculate_concentrations_self_consistently(cH_res, cOH_res, cNa_res, cCl_res)
                else:
                    return cH_res, cOH_res, cNa_res, cCl_res
            return calculate_concentrations_self_consistently(cH_res, cOH_res, cNa_res, cCl_res)

        cH_res, cOH_res, cNa_res, cCl_res = determine_reservoir_concentrations_selfconsistently(cH_res, c_salt_res)
        ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
        determined_pH = -np.log10(cH_res.to('mol/L').magnitude * np.sqrt(activity_coefficient_monovalent_pair(ionic_strength_res)))

        while abs(determined_pH-pH_res)>1e-6:
            if determined_pH > pH_res:
                cH_res *= 1.005
            else:
                cH_res /= 1.003
            cH_res, cOH_res, cNa_res, cCl_res = determine_reservoir_concentrations_selfconsistently(cH_res, c_salt_res)
            ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
            determined_pH = -np.log10(cH_res.to('mol/L').magnitude * np.sqrt(activity_coefficient_monovalent_pair(ionic_strength_res)))
            self_consistent_run=0

        return cH_res, cOH_res, cNa_res, cCl_res

    def enable_motion_of_rigid_object(self, name, espresso_system):
        '''
        Enables the motion of the rigid object `name` in the `espresso_system`.

        Args:
            name(`str`): Label of the object.
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.

        Note:
            - It requires that espressomd has the following features activated: ["VIRTUAL_SITES_RELATIVE", "MASS"].
        '''
        logging.info('enable_motion_of_rigid_object requires that espressomd has the following features activated: ["VIRTUAL_SITES_RELATIVE", "MASS"]')
        self._check_supported_molecule(molecule_name=name,
                                        valid_pmb_types= ['protein'])
        molecule_ids_list = self.df.loc[self.df['name']==name].molecule_id.to_list()
        for molecule_id in molecule_ids_list:    
            particle_ids_list = self.df.loc[self.df['molecule_id']==molecule_id].particle_id.dropna().to_list()
            center_of_mass = self.calculate_center_of_mass_of_molecule ( molecule_id=molecule_id,espresso_system=espresso_system)
            rigid_object_center = espresso_system.part.add(pos=center_of_mass,
                                                           rotation=[True,True,True], 
                                                           type=self.propose_unused_type())
            
            rigid_object_center.mass = len(particle_ids_list)
            momI = 0
            for pid in particle_ids_list:
                momI += np.power(np.linalg.norm(center_of_mass - espresso_system.part.by_id(pid).pos), 2)

            rigid_object_center.rinertia = np.ones(3) * momI
            
            for particle_id in particle_ids_list:
                pid = espresso_system.part.by_id(particle_id)
                pid.vs_auto_relate_to(rigid_object_center.id)
        return

    def filter_df(self, pmb_type):
        """
        Filters `pmb.df` and returns a sub-set of it containing only rows with pmb_object_type=`pmb_type` and non-NaN columns.
        
        Args:
            pmb_type(`str`): pmb_object_type to filter in `pmb.df`.

        Returns:
            pmb_type_df(`Pandas.Dataframe`): filtered `pmb.df`.
        """
        pmb_type_df = self.df.loc[self.df['pmb_type']== pmb_type]
        pmb_type_df = pmb_type_df.dropna( axis=1, thresh=1)
        return pmb_type_df   

    def find_value_from_es_type(self, es_type, column_name):
        """
        Finds a value in `pmb.df` for a `column_name` and `es_type` pair.

        Args:
            es_type(`int`): value of the espresso type
            column_name(`str`): name of the column in `pymbe.df`

        Returns:
            Value in `pymbe.df` matching  `column_name` and `es_type`
        """
        idx = pd.IndexSlice
        for state in ['state_one', 'state_two']:            
            index = self.df.loc[self.df[(state, 'es_type')] == es_type].index
            if len(index) > 0:
                if column_name == 'label':
                    label = self.df.loc[idx[index[0]], idx[(state,column_name)]]
                    return label
                else: 
                    column_name_value = self.df.loc[idx[index[0]], idx[(column_name,'')]]
                    return column_name_value

    def format_node(self, node_list):
        return "[" + " ".join(map(str, node_list)) + "]"


    def generate_coordinates_outside_sphere(self, center, radius, max_dist, n_samples):
        """
        Generates coordinates outside a sphere centered at `center`.

        Args:
            center(`lst`): Coordinates of the center of the sphere.
            radius(`float`): Radius of the sphere.
            max_dist(`float`): Maximum distance from the center of the spahre to generate coordinates.
            n_samples(`int`): Number of sample points.

        Returns:
            coord_list(`lst`): Coordinates of the sample points.
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
        Uniformly samples points from a hypersphere. If on_surface is set to True, the points are
        uniformly sampled from the surface of the hypersphere.
        
        Args:
            center(`lst`): Array with the coordinates of the center of the spheres.
            radius(`float`): Radius of the sphere.
            n_samples(`int`): Number of sample points to generate inside the sphere.
            on_surface (`bool`, optional): If set to True, points will be uniformly sampled on the surface of the hypersphere.

        Returns:
            samples(`list`): Coordinates of the sample points inside the hypersphere.
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
        Generates an orthogonal vector to the input `vector`.

        Args:
            vector(`lst`): arbitrary vector.
            magnitude(`float`): magnitude of the orthogonal vector.
            
        Returns:
            (`lst`): Orthogonal vector with the same magnitude as the input vector.
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

    def get_bond_length(self, particle_name1, particle_name2, hard_check=False, use_default_bond=False) :
        """
        Searches for bonds between the particle types given by `particle_name1` and `particle_name2` in `pymbe.df` and returns the initial bond length.
        If `use_default_bond` is activated and a "default" bond is defined, returns the length of that default bond instead.
        If no bond is found, it prints a message and it does not return anything. If `hard_check` is activated, the code stops if no bond is found.

        Args:
            particle_name1(str): label of the type of the first particle type of the bonded particles.
            particle_name2(str): label of the type of the second particle type of the bonded particles.
            hard_check(bool, optional): If it is activated, the code stops if no bond is found. Defaults to False. 
            use_default_bond(bool, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to False. 

        Returns:
            l0(`pint.Quantity`): bond length
        
        Note:
            - If `use_default_bond`=True and no bond is defined between `particle_name1` and `particle_name2`, it returns the default bond defined in `pmb.df`.
            - If `hard_check`=`True` stops the code when no bond is found.
        """
        bond_key = _DFm._find_bond_key(df = self.df,
                                       particle_name1 = particle_name1, 
                                       particle_name2 = particle_name2, 
                                       use_default_bond = use_default_bond)
        if bond_key:
            return self.df[self.df['name'] == bond_key].l0.values[0]
        else:
            msg = f"Bond not defined between particles {particle_name1} and {particle_name2}"
            if hard_check:
                raise ValueError(msg)     
            else:
                logging.warning(msg)
                return

    def get_charge_number_map(self):
        '''
        Gets the charge number of each `espresso_type` in `pymbe.df`.
        
        Returns:
            charge_number_map(`dict`): {espresso_type: z}.
        '''
        if self.df.state_one['es_type'].isnull().values.any():         
            df_state_one = self.df.state_one.dropna()     
            df_state_two = self.df.state_two.dropna()  
        else:    
            df_state_one = self.df.state_one
            if self.df.state_two['es_type'].isnull().values.any():
                df_state_two = self.df.state_two.dropna()   
            else:
                df_state_two = self.df.state_two
        state_one = pd.Series (df_state_one.z.values,index=df_state_one.es_type.values)
        state_two = pd.Series (df_state_two.z.values,index=df_state_two.es_type.values)
        charge_number_map  = pd.concat([state_one,state_two],axis=0).to_dict()
        return charge_number_map

    def get_lj_parameters(self, particle_name1, particle_name2, combining_rule='Lorentz-Berthelot'):
        """
        Returns the Lennard-Jones parameters for the interaction between the particle types given by 
        `particle_name1` and `particle_name2` in `pymbe.df`, calculated according to the provided combining rule.

        Args:
            particle_name1 (str): label of the type of the first particle type
            particle_name2 (str): label of the type of the second particle type
            combining_rule (`string`, optional): combining rule used to calculate `sigma` and `epsilon` for the potential betwen a pair of particles. Defaults to 'Lorentz-Berthelot'.

        Returns:
            {"epsilon": epsilon_value, "sigma": sigma_value, "offset": offset_value, "cutoff": cutoff_value}

        Note:
            - Currently, the only `combining_rule` supported is Lorentz-Berthelot.
            - If the sigma value of `particle_name1` or `particle_name2` is 0, the function will return an empty dictionary. No LJ interactions are set up for particles with sigma = 0.
        """
        supported_combining_rules=["Lorentz-Berthelot"]
        lj_parameters_keys=["sigma","epsilon","offset","cutoff"]
        if combining_rule not in supported_combining_rules:
            raise ValueError(f"Combining_rule {combining_rule} currently not implemented in pyMBE, valid keys are {supported_combining_rules}")
        lj_parameters={}
        for key in lj_parameters_keys:
            lj_parameters[key]=[]
        # Search the LJ parameters of the type pair
        for name in [particle_name1,particle_name2]:
            for key in lj_parameters_keys:
                lj_parameters[key].append(self.df[self.df.name == name][key].values[0])
        # If one of the particle has sigma=0, no LJ interations are set up between that particle type and the others    
        if not all(sigma_value.magnitude for sigma_value in lj_parameters["sigma"]):
            return {}
        # Apply combining rule
        if combining_rule == 'Lorentz-Berthelot':
            lj_parameters["sigma"]=(lj_parameters["sigma"][0]+lj_parameters["sigma"][1])/2
            lj_parameters["cutoff"]=(lj_parameters["cutoff"][0]+lj_parameters["cutoff"][1])/2
            lj_parameters["offset"]=(lj_parameters["offset"][0]+lj_parameters["offset"][1])/2
            lj_parameters["epsilon"]=np.sqrt(lj_parameters["epsilon"][0]*lj_parameters["epsilon"][1])
        return lj_parameters

    def get_metal_ions_charge_number_map(self):
        """
        Gets a map with the charge numbers of all the metal ions supported.

        Returns:
            metal_charge_number_map(dict): Has the structure {"metal_name": metal_charge_number}

        """
        metal_charge_number_map = {"Ca": 2}
        return metal_charge_number_map

    def get_particle_id_map(self, object_name):
        '''
        Gets all the ids associated with the object with name `object_name` in `pmb.df`

        Args:
            object_name(`str`): name of the object
        
        Returns:
            id_map(`dict`): dict of the structure {"all": [all_ids_with_object_name], "residue_map": {res_id: [particle_ids_in_res_id]}, "molecule_map": {mol_id: [particle_ids_in_mol_id]}, }
        '''
        object_type=self._check_supported_molecule(molecule_name=object_name,
                                                   valid_pmb_types= ['particle','residue','molecule',"peptide","protein"])
        id_list = []
        mol_map = {}
        res_map = {}
        def do_res_map(res_ids):
            for res_id in res_ids:
                res_list=self.df.loc[(self.df['residue_id']== res_id) & (self.df['pmb_type']== "particle")].particle_id.dropna().tolist()
                res_map[res_id]=res_list
            return res_map
        if object_type in ['molecule', 'protein', 'peptide']:
            mol_ids = self.df.loc[self.df['name']== object_name].molecule_id.dropna().tolist()
            for mol_id in mol_ids:
                res_ids = set(self.df.loc[(self.df['molecule_id']== mol_id) & (self.df['pmb_type']== "particle") ].residue_id.dropna().tolist())
                res_map=do_res_map(res_ids=res_ids)    
                mol_list=self.df.loc[(self.df['molecule_id']== mol_id) & (self.df['pmb_type']== "particle")].particle_id.dropna().tolist()
                id_list+=mol_list
                mol_map[mol_id]=mol_list
        elif object_type == 'residue':     
            res_ids = self.df.loc[self.df['name']== object_name].residue_id.dropna().tolist()
            res_map=do_res_map(res_ids=res_ids)
            id_list=[]
            for res_id_list in res_map.values():
                id_list+=res_id_list
        elif object_type == 'particle':
            id_list = self.df.loc[self.df['name']== object_name].particle_id.dropna().tolist()
        return {"all": id_list, "molecule_map": mol_map, "residue_map": res_map}

    def get_pka_set(self):
        '''
        Gets the pka-values and acidities of the particles with acid/base properties in `pmb.df`
        
        Returns:
            pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}
        '''
        titratables_AA_df = self.df[[('name',''),('pka',''),('acidity','')]].drop_duplicates().dropna()
        pka_set = {}
        for index in titratables_AA_df.name.keys():
            name = titratables_AA_df.name[index]
            pka_value = titratables_AA_df.pka[index]
            acidity = titratables_AA_df.acidity[index]   
            pka_set[name] = {'pka_value':pka_value,'acidity':acidity}
        return pka_set 
    
    def get_radius_map(self, dimensionless=True):
        '''
        Gets the effective radius of each `espresso type` in `pmb.df`. 

        Args:
            dimensionless('bool'): controlls if the returned radii have a dimension. Defaults to False.
        
        Returns:
            radius_map(`dict`): {espresso_type: radius}.

        Note:
            The radius corresponds to (sigma+offset)/2
        '''
        df_state_one = self.df[[('sigma',''),('offset',''),('state_one','es_type')]].dropna().drop_duplicates()
        df_state_two = self.df[[('sigma',''),('offset',''),('state_two','es_type')]].dropna().drop_duplicates()
        state_one = pd.Series((df_state_one.sigma.values+df_state_one.offset.values)/2.0,index=df_state_one.state_one.es_type.values)
        state_two = pd.Series((df_state_two.sigma.values+df_state_two.offset.values)/2.0,index=df_state_two.state_two.es_type.values)
        radius_map  = pd.concat([state_one,state_two],axis=0).to_dict()  
        if dimensionless:
            for key in radius_map:
                radius_map[key] = radius_map[key].magnitude
        return radius_map

    def get_reduced_units(self):
        """
        Returns the  current set of reduced units defined in pyMBE.

        Returns:
            reduced_units_text(`str`): text with information about the current set of reduced units.

        """
        unit_length=self.units.Quantity(1,'reduced_length')
        unit_energy=self.units.Quantity(1,'reduced_energy')
        unit_charge=self.units.Quantity(1,'reduced_charge')
        reduced_units_text = "\n".join(["Current set of reduced units:",
                                       f"{unit_length.to('nm'):.5g} = {unit_length}",
                                       f"{unit_energy.to('J'):.5g} = {unit_energy}",
                                       f"{unit_charge.to('C'):.5g} = {unit_charge}",
                                       f"Temperature: {(self.kT/self.kB).to('K'):.5g}"
                                        ])   
        return reduced_units_text

    def get_type_map(self):
        """
        Gets all different espresso types assigned to particles  in `pmb.df`.
        
        Returns:
            type_map(`dict`): {"name": espresso_type}.
        """
        df_state_one = self.df.state_one.dropna(how='all')     
        df_state_two = self.df.state_two.dropna(how='all')  
        state_one = pd.Series (df_state_one.es_type.values,index=df_state_one.label)
        state_two = pd.Series (df_state_two.es_type.values,index=df_state_two.label)
        type_map  = pd.concat([state_one,state_two],axis=0).to_dict()
        return type_map

    def initialize_lattice_builder(self, diamond_lattice):
        """
        Initialize the lattice builder with the DiamondLattice object.

        Args:
            diamond_lattice(`DiamondLattice`): DiamondLattice object from the `lib/lattice` module to be used in the LatticeBuilder.
        """
        from .lib.lattice import LatticeBuilder, DiamondLattice
        if not isinstance(diamond_lattice, DiamondLattice):
            raise TypeError("Currently only DiamondLattice objects are supported.")
        self.lattice_builder = LatticeBuilder(lattice=diamond_lattice)
        logging.info(f"LatticeBuilder initialized with mpc={diamond_lattice.mpc} and box_l={diamond_lattice.box_l}")
        return self.lattice_builder

    def load_interaction_parameters(self, filename, overwrite=False):
        """
        Loads the interaction parameters stored in `filename` into `pmb.df`
        
        Args:
            filename(`str`): name of the file to be read
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False. 
        """
        without_units = ['z','es_type']
        with_units = ['sigma','epsilon','offset','cutoff']

        with open(filename, 'r') as f:
            interaction_data = json.load(f)
            interaction_parameter_set = interaction_data["data"]

        for key in interaction_parameter_set:
            param_dict=interaction_parameter_set[key]
            object_type=param_dict.pop('object_type')
            if object_type == 'particle':
                not_required_attributes={}    
                for not_required_key in without_units+with_units:
                    if not_required_key in param_dict.keys():
                        if not_required_key in with_units:
                            not_required_attributes[not_required_key] = _DFm._create_variable_with_units(variable=param_dict.pop(not_required_key), 
                                                                                                         units_registry=self.units)
                        elif not_required_key in without_units:
                            not_required_attributes[not_required_key]=param_dict.pop(not_required_key)
                    else:
                        not_required_attributes[not_required_key]=pd.NA
                self.define_particle(name=param_dict.pop('name'),
                                z=not_required_attributes.pop('z'),
                                sigma=not_required_attributes.pop('sigma'),
                                offset=not_required_attributes.pop('offset'),
                                cutoff=not_required_attributes.pop('cutoff'),
                                epsilon=not_required_attributes.pop('epsilon'),
                                overwrite=overwrite)
            elif object_type == 'residue':
                self.define_residue(**param_dict)
            elif object_type == 'molecule':
                self.define_molecule(**param_dict)
            elif object_type == 'peptide':
                self.define_peptide(**param_dict)
            elif object_type == 'bond':
                particle_pairs = param_dict.pop('particle_pairs')
                bond_parameters = param_dict.pop('bond_parameters')
                bond_type = param_dict.pop('bond_type')
                if bond_type == 'harmonic':
                    k =  _DFm._create_variable_with_units(variable=bond_parameters.pop('k'), 
                                                          units_registry=self.units)
                    r_0 = _DFm._create_variable_with_units(variable=bond_parameters.pop('r_0'), 
                                                          units_registry=self.units)
                    bond = {'r_0'    : r_0,
                            'k'      : k,
                            }

                elif bond_type == 'FENE':
                    k = _DFm._create_variable_with_units(variable=bond_parameters.pop('k'), 
                                                         units_registry=self.units)
                    r_0 = _DFm._create_variable_with_units(variable=bond_parameters.pop('r_0'), 
                                                           units_registry=self.units)
                    d_r_max = _DFm._create_variable_with_units(variable=bond_parameters.pop('d_r_max'), 
                                                               units_registry=self.units)
                    bond =  {'r_0'    : r_0,
                             'k'      : k,
                             'd_r_max': d_r_max,
                             }
                else:
                    raise ValueError("Current implementation of pyMBE only supports harmonic and FENE bonds")
                self.define_bond(bond_type=bond_type, 
                                bond_parameters=bond, 
                                particle_pairs=particle_pairs)
            else:
                raise ValueError(object_type+' is not a known pmb object type')
            
        return
    
    def load_pka_set(self, filename, overwrite=True):
        """
        Loads the pka_set stored in `filename` into `pmb.df`.
        
        Args:
            filename(`str`): name of the file with the pka set to be loaded. Expected format is {name:{"acidity": acidity, "pka_value":pka_value}}.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to True. 
        """
        with open(filename, 'r') as f:
            pka_data = json.load(f)
            pka_set = pka_data["data"]

        self.check_pka_set(pka_set=pka_set)

        for key in pka_set:
            acidity = pka_set[key]['acidity']
            pka_value = pka_set[key]['pka_value']
            self.set_particle_acidity(name=key, 
                                      acidity=acidity, 
                                      pka=pka_value, 
                                      overwrite=overwrite)
        return


    def propose_unused_type(self):
        """
        Searches in `pmb.df` all the different particle types defined and returns a new one.

        Returns:
            unused_type(`int`): unused particle type
        """
        type_map = self.get_type_map()
        if not type_map:  
            unused_type = 0
        else:
            valid_values = [v for v in type_map.values() if pd.notna(v)]  # Filter out pd.NA values
            unused_type = max(valid_values) + 1 if valid_values else 0  # Ensure max() doesn't fail if all values are NA
        return unused_type

    def protein_sequence_parser(self, sequence):
        '''
        Parses `sequence` to the one letter code for amino acids.
        
        Args:
            sequence(`str` or `lst`): Sequence of the amino acid. 

        Returns:
            clean_sequence(`lst`): `sequence` using the one letter code.
        
        Note:
            - Accepted formats for `sequence` are:
                - `lst` with one letter or three letter code of each aminoacid in each element
                - `str` with the sequence using the one letter code
                - `str` with the squence using the three letter code, each aminoacid must be separated by a hyphen "-"
        
        '''
        # Aminoacid key
        keys={"ALA": "A",
                "ARG": "R",
                "ASN": "N",
                "ASP": "D",
                "CYS": "C",
                "GLU": "E",
                "GLN": "Q",
                "GLY": "G",
                "HIS": "H",
                "ILE": "I",
                "LEU": "L",
                "LYS": "K",
                "MET": "M",
                "PHE": "F",
                "PRO": "P",
                "SER": "S",
                "THR": "T",
                "TRP": "W",
                "TYR": "Y",
                "VAL": "V",
                "PSER": "J",
                "PTHR": "U",
                "PTyr": "Z",
                "NH2": "n",
                "COOH": "c"}
        clean_sequence=[]
        if isinstance(sequence, str):
            if sequence.find("-") != -1:
                splited_sequence=sequence.split("-")
                for residue in splited_sequence:
                    if len(residue) == 1:
                        if residue in keys.values():
                            residue_ok=residue
                        else:
                            if residue.upper() in keys.values():
                                residue_ok=residue.upper()
                            else:
                                raise ValueError("Unknown one letter code for a residue given: ", residue, " please review the input sequence")
                        clean_sequence.append(residue_ok)
                    else:
                        if residue in keys.keys():
                            clean_sequence.append(keys[residue])
                        else:
                            if residue.upper() in keys.keys():
                                clean_sequence.append(keys[residue.upper()])
                            else:
                                raise ValueError("Unknown  code for a residue: ", residue, " please review the input sequence")
            else:
                for residue in sequence:
                    if residue in keys.values():
                        residue_ok=residue
                    else:
                        if residue.upper() in keys.values():
                            residue_ok=residue.upper()
                        else:
                            raise ValueError("Unknown one letter code for a residue: ", residue, " please review the input sequence")
                    clean_sequence.append(residue_ok)
        if isinstance(sequence, list):
            for residue in sequence:
                if residue in keys.values():
                    residue_ok=residue
                else:
                    if residue.upper() in keys.values():
                        residue_ok=residue.upper()
                    elif (residue.upper() in keys.keys()):
                        residue_ok= keys[residue.upper()]
                    else:
                        raise ValueError("Unknown code for a residue: ", residue, " please review the input sequence")
                clean_sequence.append(residue_ok)
        return clean_sequence
    
    def read_pmb_df (self,filename):
        """
        Reads a pyMBE's Dataframe stored in `filename` and stores the information into pyMBE.

        Args:
            filename(`str`): path to file.

        Note:
            This function only accepts files with CSV format. 
        """
        if filename.rsplit(".", 1)[1] != "csv":
            raise ValueError("Only files with CSV format are supported")
        df = pd.read_csv (filename,header=[0, 1], index_col=0)
        self.df = _DFm._setup_df()
        columns_names = pd.MultiIndex.from_frame(self.df)
        columns_names = columns_names.names
        multi_index = pd.MultiIndex.from_tuples(columns_names)
        df.columns = multi_index
        _DFm._convert_columns_to_original_format(df=df, 
                                                 units_registry=self.units)
        self.df = df            
        self.df.fillna(pd.NA, 
                       inplace=True)
        return self.df
    
    def read_protein_vtf_in_df (self,filename,unit_length=None):
        """
        Loads a coarse-grained protein model in a vtf file `filename` into the `pmb.df` and it labels it with `name`.

        Args:
            filename(`str`): Path to the vtf file with the coarse-grained model.
            unit_length(`obj`): unit of length of the the coordinates in `filename` using the pyMBE UnitRegistry. Defaults to None.

        Returns:
            topology_dict(`dict`): {'initial_pos': coords_list, 'chain_id': id, 'sigma': sigma_value}

        Note:
            - If no `unit_length` is provided, it is assumed that the coordinates are in Angstrom.
        """

        logging.info(f'Loading protein coarse grain model file: {filename}')

        coord_list = []
        particles_dict = {}

        if unit_length is None:
            unit_length = 1 * self.units.angstrom 

        with open (filename,'r') as protein_model:
            for line in protein_model :
                line_split = line.split()
                if line_split:
                    line_header = line_split[0]
                    if line_header == 'atom':
                        atom_id  = line_split[1]
                        atom_name = line_split[3]
                        atom_resname = line_split[5]
                        chain_id = line_split[9]
                        radius = float(line_split [11])*unit_length 
                        particles_dict [int(atom_id)] = [atom_name , atom_resname, chain_id, radius]
                    elif line_header.isnumeric(): 
                        atom_coord = line_split[1:] 
                        atom_coord = [(float(i)*unit_length).to('reduced_length').magnitude for i in atom_coord]
                        coord_list.append (atom_coord)

        numbered_label = []
        i = 0   
        
        for atom_id in particles_dict.keys():
    
            if atom_id == 1:
                atom_name = particles_dict[atom_id][0]
                numbered_name = [f'{atom_name}{i}',particles_dict[atom_id][2],particles_dict[atom_id][3]]
                numbered_label.append(numbered_name)

            elif atom_id != 1: 
            
                if particles_dict[atom_id-1][1] != particles_dict[atom_id][1]:
                    i += 1                    
                    count = 1
                    atom_name = particles_dict[atom_id][0]
                    numbered_name = [f'{atom_name}{i}',particles_dict[atom_id][2],particles_dict[atom_id][3]]
                    numbered_label.append(numbered_name)
                    
                elif particles_dict[atom_id-1][1] == particles_dict[atom_id][1]:
                    if count == 2 or particles_dict[atom_id][1] == 'GLY':
                        i +=1  
                        count = 0
                    atom_name = particles_dict[atom_id][0]
                    numbered_name = [f'{atom_name}{i}',particles_dict[atom_id][2],particles_dict[atom_id][3]]
                    numbered_label.append(numbered_name)
                    count +=1

        topology_dict = {}

        for i in range (0, len(numbered_label)):   
            topology_dict [numbered_label[i][0]] = {'initial_pos': coord_list[i] ,
                                                    'chain_id':numbered_label[i][1],
                                                    'radius':numbered_label[i][2] }

        return topology_dict

    def search_bond(self, particle_name1, particle_name2, hard_check=False, use_default_bond=False) :
        """
        Searches for bonds between the particle types given by `particle_name1` and `particle_name2` in `pymbe.df` and returns it.
        If `use_default_bond` is activated and a "default" bond is defined, returns that default bond instead.
        If no bond is found, it prints a message and it does not return anything. If `hard_check` is activated, the code stops if no bond is found.

        Args:
            particle_name1(`str`): label of the type of the first particle type of the bonded particles.
            particle_name2(`str`): label of the type of the second particle type of the bonded particles.
            hard_check(`bool`, optional): If it is activated, the code stops if no bond is found. Defaults to False. 
            use_default_bond(`bool`, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to False. 

        Returns:
            bond(`espressomd.interactions.BondedInteractions`): bond object from the espressomd library.
        
        Note:
            - If `use_default_bond`=True and no bond is defined between `particle_name1` and `particle_name2`, it returns the default bond defined in `pmb.df`.
            - If `hard_check`=`True` stops the code when no bond is found.
        """

        bond_key = _DFm._find_bond_key(df = self.df,
                                       particle_name1 = particle_name1,
                                       particle_name2 = particle_name2,
                                       use_default_bond = use_default_bond)
        if use_default_bond:
            if not _DFm._check_if_name_is_defined_in_df(name="default", df=self.df):
                raise ValueError(f"use_default_bond is set to {use_default_bond} but no default bond has been defined. Please define a default bond with pmb.define_default_bond")
        if bond_key:
            return self.df[self.df['name']==bond_key].bond_object.values[0]
        else:
            msg= f"Bond not defined between particles {particle_name1} and {particle_name2}"
            if hard_check:
                raise ValueError(msg)
            else:
                logging.warning(msg)
            return None
    def search_particles_in_residue(self, residue_name):
        '''
        Searches for all particles in a given residue of name `residue_name`.

        Args:
            residue_name (`str`): name of the residue to be searched

        Returns:
            list_of_particles_in_residue (`lst`): list of the names of all particles in the residue

        Note:
            - The function returns a name per particle in residue, i.e. if there are multiple particles with the same type `list_of_particles_in_residue` will have repeated items.
            - The function will return an empty list if the residue is not defined in `pmb.df`.
            - The function will return an empty list if the particles are not defined in the pyMBE DataFrame.
        '''
        if not _DFm._check_if_name_is_defined_in_df(name=residue_name, df=self.df):
            logging.warning(f"Residue {residue_name} not defined in pmb.df")
            return []
        self._check_if_name_has_right_type(name=residue_name, expected_pmb_type="residue")
        index_residue = self.df.loc[self.df['name'] == residue_name].index[0].item() 
        central_bead = self.df.at [index_residue, ('central_bead', '')]
        list_of_side_chains = self.df.at[index_residue, ('side_chains', '')]
        list_of_particles_in_residue = []
        if central_bead is not pd.NA:
            if _DFm._check_if_name_is_defined_in_df(name=central_bead, df=self.df):
                if self._check_if_name_has_right_type(name=central_bead, expected_pmb_type="particle", hard_check=False):
                    list_of_particles_in_residue.append(central_bead)
        if list_of_side_chains is not pd.NA:
            for side_chain in list_of_side_chains:
                if _DFm._check_if_name_is_defined_in_df(name=side_chain, df=self.df): 
                    object_type = self.df[self.df['name']==side_chain].pmb_type.values[0]
                else:
                    continue
                if object_type == "residue":
                    list_of_particles_in_side_chain_residue = self.search_particles_in_residue(side_chain)
                    list_of_particles_in_residue += list_of_particles_in_side_chain_residue
                elif object_type == "particle":
                    if side_chain is not pd.NA:
                        list_of_particles_in_residue.append(side_chain)
        return list_of_particles_in_residue        

    def set_particle_acidity(self, name, acidity=pd.NA, default_charge_number=0, pka=pd.NA, overwrite=True):
        """
        Sets the particle acidity including the charges in each of its possible states. 

        Args:
            name(`str`): Unique label that identifies the `particle`. 
            acidity(`str`): Identifies whether the particle is `acidic` or `basic`, used to setup constant pH simulations. Defaults to None.
            default_charge_number (`int`): Charge number of the particle. Defaults to 0.
            pka(`float`, optional):  If `particle` is an acid or a base, it defines its pka-value. Defaults to pandas.NA.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False. 
     
        Note:
            - For particles with  `acidity = acidic` or `acidity = basic`, `state_one` and `state_two` correspond to the protonated and 
        deprotonated states, respectively. 
            - For particles without an acidity `acidity = pandas.NA`, only `state_one` is defined.
            - Each state has the following properties as sub-indexes: `label`,`charge` and `es_type`. 
        """
        acidity_valid_keys = ['inert','acidic', 'basic']
        if not pd.isna(acidity):
            if acidity not in acidity_valid_keys:
                raise ValueError(f"Acidity {acidity} provided for particle name  {name} is not supproted. Valid keys are: {acidity_valid_keys}")
            if acidity in ['acidic', 'basic'] and pd.isna(pka):
                raise ValueError(f"pKa not provided for particle with name {name} with acidity {acidity}. pKa must be provided for acidic or basic particles.")   
            if acidity == "inert":
                acidity = pd.NA
                logging.warning("the keyword 'inert' for acidity has been replaced by setting acidity = pd.NA. For backwards compatibility, acidity has been set to pd.NA. Support for `acidity = 'inert'` may be deprecated in future releases of pyMBE")

        self._define_particle_entry_in_df(name=name)
        
        for index in self.df[self.df['name']==name].index:       
            if pka is not pd.NA:
                _DFm._add_value_to_df(df = self.df,
                                      key = ('pka',''),
                                      index = index,
                                      new_value = pka, 
                                      overwrite = overwrite)

            _DFm._add_value_to_df(df = self.df,
                                  key = ('acidity',''),
                                  index = index,
                                  new_value = acidity, 
                                  overwrite = overwrite) 
            if not _DFm._check_if_df_cell_has_a_value(df=self.df, index=index,key=('state_one','es_type')):
                _DFm._add_value_to_df(df = self.df,
                                      key = ('state_one','es_type'),
                                      index = index,
                                      new_value = self.propose_unused_type(),
                                      overwrite = overwrite)
            if pd.isna(self.df.loc [self.df['name']  == name].acidity.iloc[0]):
                _DFm._add_value_to_df(df = self.df,
                                      key = ('state_one','z'),
                                      index = index,
                                      new_value = default_charge_number,
                                      overwrite = overwrite)
                _DFm._add_value_to_df(df = self.df,
                                      key = ('state_one','label'),
                                      index = index,
                                      new_value = name,
                                      overwrite = overwrite)
            else:
                protonated_label = f'{name}H'
                _DFm._add_value_to_df(df = self.df,
                                      key = ('state_one','label'),
                                      index = index,
                                      new_value = protonated_label,
                                      overwrite = overwrite)
                _DFm._add_value_to_df(df = self.df,
                                      key = ('state_two','label'),
                                      index = index,
                                      new_value = name,
                                      overwrite = overwrite)
                if not _DFm._check_if_df_cell_has_a_value(df=self.df, index=index,key=('state_two','es_type')):
                    _DFm._add_value_to_df(df = self.df,
                                          key = ('state_two','es_type'),
                                          index = index,
                                          new_value = self.propose_unused_type(),
                                          overwrite = overwrite)
                if self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'acidic':
                    _DFm._add_value_to_df(df = self.df,
                                          key = ('state_one','z'),
                                          index = index,
                                          new_value = 0,
                                          overwrite = overwrite)
                    _DFm._add_value_to_df(df = self.df,
                                          key = ('state_two','z'),
                                          index = index,
                                          new_value = -1,
                                          overwrite = overwrite)
                elif self.df.loc [self.df['name']  == name].acidity.iloc[0] == 'basic':
                    _DFm._add_value_to_df(df = self.df,
                                          key = ('state_one','z'),
                                          index = index,
                                          new_value = +1,
                                          overwrite = overwrite)
                    _DFm._add_value_to_df(df = self.df,
                                          key = ('state_two','z'),
                                          index = index,
                                          new_value = 0,
                                          overwrite = overwrite)
        self.df.fillna(pd.NA, inplace=True)
        return
    
    def set_reduced_units(self, unit_length=None, unit_charge=None, temperature=None, Kw=None):
        """
        Sets the set of reduced units used by pyMBE.units and it prints it.

        Args:
            unit_length(`pint.Quantity`,optional): Reduced unit of length defined using the `pmb.units` UnitRegistry. Defaults to None. 
            unit_charge(`pint.Quantity`,optional): Reduced unit of charge defined using the `pmb.units` UnitRegistry. Defaults to None. 
            temperature(`pint.Quantity`,optional): Temperature of the system, defined using the `pmb.units` UnitRegistry. Defaults to None. 
            Kw(`pint.Quantity`,optional): Ionic product of water in mol^2/l^2. Defaults to None. 

        Note:
            - If no `temperature` is given, a value of 298.15 K is assumed by default.
            - If no `unit_length` is given, a value of 0.355 nm is assumed by default.
            - If no `unit_charge` is given, a value of 1 elementary charge is assumed by default. 
            - If no `Kw` is given, a value of 10^(-14) * mol^2 / l^2 is assumed by default. 
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
            self.check_dimensionality(variable,dimensionality)
        self.Kw=Kw*self.units.mol**2 / (self.units.l**2)
        self.kT=temperature*self.kB
        self.units._build_cache()
        self.units.define(f'reduced_energy = {self.kT} ')
        self.units.define(f'reduced_length = {unit_length}')
        self.units.define(f'reduced_charge = {unit_charge}')
        logging.info(self.get_reduced_units())
        return

    def setup_cpH (self, counter_ion, constant_pH, exclusion_range=None, pka_set=None, use_exclusion_radius_per_type = False):
        """
        Sets up the Acid/Base reactions for acidic/basic `particles` defined in `pmb.df` to be sampled in the constant pH ensemble. 

        Args:
            counter_ion(`str`): `name` of the counter_ion `particle`.
            constant_pH(`float`): pH-value.
            exclusion_range(`pint.Quantity`, optional): Below this value, no particles will be inserted.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius for each espresso_type is used. Defaults to `False`.
            pka_set(`dict`,optional): Desired pka_set, pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}. Defaults to None.

        Returns:
            RE(`reaction_methods.ConstantpHEnsemble`): Instance of a reaction_methods.ConstantpHEnsemble object from the espressomd library.
            sucessfull_reactions_labels(`lst`): Labels of the reactions set up by pyMBE.
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.get_radius_map().values())*2.0
        if pka_set is None:
            pka_set=self.get_pka_set()    
        self.check_pka_set(pka_set=pka_set)
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        
        RE = reaction_methods.ConstantpHEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range, 
                                                    seed=self.seed, 
                                                    constant_pH=constant_pH,
                                                    exclusion_radius_per_type = exclusion_radius_per_type
                                                    )
        sucessfull_reactions_labels=[]
        charge_number_map = self.get_charge_number_map()
        for name in pka_set.keys():
            if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
                logging.warning(f'The acid-base reaction of {name} has not been set up because its particle type is not defined in the pyMBE DataFrame.')
                continue
            gamma=10**-pka_set[name]['pka_value']
            state_one_type   = self.df.loc[self.df['name']==name].state_one.es_type.values[0]
            state_two_type   = self.df.loc[self.df['name']==name].state_two.es_type.values[0]
            counter_ion_type = self.df.loc[self.df['name']==counter_ion].state_one.es_type.values[0]
            RE.add_reaction(gamma=gamma,
                            reactant_types=[state_one_type],
                            product_types=[state_two_type, counter_ion_type],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            counter_ion_type: charge_number_map[counter_ion_type]})
            sucessfull_reactions_labels.append(name)
        return RE, sucessfull_reactions_labels

    def setup_gcmc(self, c_salt_res, salt_cation_name, salt_anion_name, activity_coefficient, exclusion_range=None, use_exclusion_radius_per_type = False):
        """
        Sets up grand-canonical coupling to a reservoir of salt.
        For reactive systems coupled to a reservoir, the grand-reaction method has to be used instead.

        Args:
            c_salt_res ('pint.Quantity'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            salt_cation_name ('str'): Name of the salt cation (e.g. Na+) particle.
            salt_anion_name ('str'): Name of the salt anion (e.g. Cl-) particle.
            activity_coefficient ('callable'): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
            exclusion_range(`pint.Quantity`, optional): For distances shorter than this value, no particles will be inserted.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius for each espresso_type is used. Defaults to `False`.

        Returns:
            RE (`reaction_methods.ReactionEnsemble`): Instance of a reaction_methods.ReactionEnsemble object from the espressomd library.
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.get_radius_map().values())*2.0
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range, 
                                                    seed=self.seed, 
                                                    exclusion_radius_per_type = exclusion_radius_per_type
                                                    )

        # Determine the concentrations of the various species in the reservoir and the equilibrium constants
        determined_activity_coefficient = activity_coefficient(c_salt_res)
        K_salt = (c_salt_res.to('1/(N_A * reduced_length**3)')**2) * determined_activity_coefficient

        salt_cation_es_type = self.df.loc[self.df['name']==salt_cation_name].state_one.es_type.values[0]
        salt_anion_es_type = self.df.loc[self.df['name']==salt_anion_name].state_one.es_type.values[0]     

        salt_cation_charge = self.df.loc[self.df['name']==salt_cation_name].state_one.z.values[0]
        salt_anion_charge = self.df.loc[self.df['name']==salt_anion_name].state_one.z.values[0]     

        if salt_cation_charge <= 0:
            raise ValueError('ERROR salt cation charge must be positive, charge ', salt_cation_charge)
        if salt_anion_charge >= 0:
            raise ValueError('ERROR salt anion charge must be negative, charge ', salt_anion_charge)

        # Grand-canonical coupling to the reservoir
        RE.add_reaction(
            gamma = K_salt.magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ salt_cation_es_type, salt_anion_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                salt_cation_es_type: salt_cation_charge, 
                salt_anion_es_type: salt_anion_charge, 
            }
        )

        return RE

    def setup_grxmc_reactions(self, pH_res, c_salt_res, proton_name, hydroxide_name, salt_cation_name, salt_anion_name, activity_coefficient, exclusion_range=None, pka_set=None, use_exclusion_radius_per_type = False):
        """
        Sets up Acid/Base reactions for acidic/basic 'particles' defined in 'pmb.df', as well as a grand-canonical coupling to a 
        reservoir of small ions. 
        This implementation uses the original formulation of the grand-reaction method by Landsgesell et al. [1].

        [1] Landsgesell, J., Hebbeker, P., Rud, O., Lunkad, R., Košovan, P., & Holm, C. (2020). Grand-reaction method for simulations of ionization equilibria coupled to ion partitioning. Macromolecules, 53(8), 3007-3020.

        Args:
            pH_res ('float): pH-value in the reservoir.
            c_salt_res ('pint.Quantity'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            proton_name ('str'): Name of the proton (H+) particle.
            hydroxide_name ('str'): Name of the hydroxide (OH-) particle.
            salt_cation_name ('str'): Name of the salt cation (e.g. Na+) particle.
            salt_anion_name ('str'): Name of the salt anion (e.g. Cl-) particle.
            activity_coefficient ('callable'): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
            exclusion_range(`pint.Quantity`, optional): For distances shorter than this value, no particles will be inserted.
            pka_set(`dict`,optional): Desired pka_set, pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}. Defaults to None.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius for each espresso_type is used. Defaults to `False`.

        Returns:
            RE (`obj`): Instance of a reaction_ensemble.ReactionEnsemble object from the espressomd library.
            sucessful_reactions_labels(`lst`): Labels of the reactions set up by pyMBE.
            ionic_strength_res ('pint.Quantity'): Ionic strength of the reservoir (useful for calculating partition coefficients).
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.get_radius_map().values())*2.0
        if pka_set is None:
            pka_set=self.get_pka_set()    
        self.check_pka_set(pka_set=pka_set)
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range, 
                                                    seed=self.seed, 
                                                    exclusion_radius_per_type = exclusion_radius_per_type
                                                    )

        # Determine the concentrations of the various species in the reservoir and the equilibrium constants
        cH_res, cOH_res, cNa_res, cCl_res = self.determine_reservoir_concentrations(pH_res, c_salt_res, activity_coefficient)
        ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
        determined_activity_coefficient = activity_coefficient(ionic_strength_res)
        K_W = cH_res.to('1/(N_A * reduced_length**3)') * cOH_res.to('1/(N_A * reduced_length**3)') * determined_activity_coefficient
        K_NACL = cNa_res.to('1/(N_A * reduced_length**3)') * cCl_res.to('1/(N_A * reduced_length**3)') * determined_activity_coefficient
        K_HCL = cH_res.to('1/(N_A * reduced_length**3)') * cCl_res.to('1/(N_A * reduced_length**3)') * determined_activity_coefficient

        proton_es_type = self.df.loc[self.df['name']==proton_name].state_one.es_type.values[0]
        hydroxide_es_type = self.df.loc[self.df['name']==hydroxide_name].state_one.es_type.values[0]     
        salt_cation_es_type = self.df.loc[self.df['name']==salt_cation_name].state_one.es_type.values[0]
        salt_anion_es_type = self.df.loc[self.df['name']==salt_anion_name].state_one.es_type.values[0]     

        proton_charge = self.df.loc[self.df['name']==proton_name].state_one.z.values[0]
        hydroxide_charge = self.df.loc[self.df['name']==hydroxide_name].state_one.z.values[0]     
        salt_cation_charge = self.df.loc[self.df['name']==salt_cation_name].state_one.z.values[0]
        salt_anion_charge = self.df.loc[self.df['name']==salt_anion_name].state_one.z.values[0]     

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
        RE.add_reaction(
            gamma = K_W.magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ proton_es_type, hydroxide_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                proton_es_type: proton_charge, 
                hydroxide_es_type: hydroxide_charge, 
            }
        )

        # 0 = Na+ + Cl-
        RE.add_reaction(
            gamma = K_NACL.magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ salt_cation_es_type, salt_anion_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                salt_cation_es_type: salt_cation_charge, 
                salt_anion_es_type: salt_anion_charge, 
            }
        )

        # 0 = Na+ + OH-
        RE.add_reaction(
            gamma = (K_NACL * K_W / K_HCL).magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ salt_cation_es_type, hydroxide_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                salt_cation_es_type: salt_cation_charge, 
                hydroxide_es_type: hydroxide_charge, 
            }
        )

        # 0 = H+ + Cl-
        RE.add_reaction(
            gamma = K_HCL.magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ proton_es_type, salt_anion_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                proton_es_type: proton_charge, 
                salt_anion_es_type: salt_anion_charge, 
            }
        )

        # Annealing moves to ensure sufficient sampling
        # Cation annealing H+ = Na+
        RE.add_reaction(
            gamma = (K_NACL / K_HCL).magnitude,
            reactant_types = [proton_es_type],
            reactant_coefficients = [ 1 ],
            product_types = [ salt_cation_es_type ],
            product_coefficients = [ 1 ],
            default_charges = {
                proton_es_type: proton_charge, 
                salt_cation_es_type: salt_cation_charge, 
            }
        )

        # Anion annealing OH- = Cl- 
        RE.add_reaction(
            gamma = (K_HCL / K_W).magnitude,
            reactant_types = [hydroxide_es_type],
            reactant_coefficients = [ 1 ],
            product_types = [ salt_anion_es_type ],
            product_coefficients = [ 1 ],
            default_charges = {
                hydroxide_es_type: hydroxide_charge, 
                salt_anion_es_type: salt_anion_charge, 
            }
        )

        sucessful_reactions_labels=[]
        charge_number_map = self.get_charge_number_map()
        for name in pka_set.keys():
            if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
                logging.warning(f'The acid-base reaction of {name} has not been set up because its particle type is not defined in the dataframe.')
                continue

            Ka = (10**-pka_set[name]['pka_value'] * self.units.mol/self.units.l).to('1/(N_A * reduced_length**3)')

            state_one_type = self.df.loc[self.df['name']==name].state_one.es_type.values[0]
            state_two_type = self.df.loc[self.df['name']==name].state_two.es_type.values[0]

            # Reaction in terms of proton: HA = A + H+
            RE.add_reaction(gamma=Ka.magnitude,
                            reactant_types=[state_one_type],
                            reactant_coefficients=[1],
                            product_types=[state_two_type, proton_es_type],
                            product_coefficients=[1, 1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            proton_es_type: proton_charge})

            # Reaction in terms of salt cation: HA = A + Na+
            RE.add_reaction(gamma=(Ka * K_NACL / K_HCL).magnitude,
                            reactant_types=[state_one_type],
                            reactant_coefficients=[1],
                            product_types=[state_two_type, salt_cation_es_type],
                            product_coefficients=[1, 1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            salt_cation_es_type: salt_cation_charge})

            # Reaction in terms of hydroxide: OH- + HA = A
            RE.add_reaction(gamma=(Ka / K_W).magnitude,
                            reactant_types=[state_one_type, hydroxide_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=[state_two_type],
                            product_coefficients=[1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            hydroxide_es_type: hydroxide_charge})

            # Reaction in terms of salt anion: Cl- + HA = A
            RE.add_reaction(gamma=(Ka / K_HCL).magnitude,
                            reactant_types=[state_one_type, salt_anion_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=[state_two_type],
                            product_coefficients=[1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            salt_anion_es_type: salt_anion_charge})

            sucessful_reactions_labels.append(name)
        return RE, sucessful_reactions_labels, ionic_strength_res

    def setup_grxmc_unified(self, pH_res, c_salt_res, cation_name, anion_name, activity_coefficient, exclusion_range=None, pka_set=None, use_exclusion_radius_per_type = False):
        """
        Sets up Acid/Base reactions for acidic/basic 'particles' defined in 'pmb.df', as well as a grand-canonical coupling to a 
        reservoir of small ions. 
        This implementation uses the formulation of the grand-reaction method by Curk et al. [1], which relies on "unified" ion types X+ = {H+, Na+} and X- = {OH-, Cl-}. 
        A function that implements the original version of the grand-reaction method by Landsgesell et al. [2] is also available under the name 'setup_grxmc_reactions'.

        [1] Curk, T., Yuan, J., & Luijten, E. (2022). Accelerated simulation method for charge regulation effects. The Journal of Chemical Physics, 156(4).
        [2] Landsgesell, J., Hebbeker, P., Rud, O., Lunkad, R., Košovan, P., & Holm, C. (2020). Grand-reaction method for simulations of ionization equilibria coupled to ion partitioning. Macromolecules, 53(8), 3007-3020.

        Args:
            pH_res ('float'): pH-value in the reservoir.
            c_salt_res ('pint.Quantity'): Concentration of monovalent salt (e.g. NaCl) in the reservoir.
            cation_name ('str'): Name of the cationic particle.
            anion_name ('str'): Name of the anionic particle.
            activity_coefficient ('callable'): A function that calculates the activity coefficient of an ion pair as a function of the ionic strength.
            exclusion_range(`pint.Quantity`, optional): Below this value, no particles will be inserted.
            pka_set(`dict`,optional): Desired pka_set, pka_set(`dict`): {"name" : {"pka_value": pka, "acidity": acidity}}. Defaults to None.
            use_exclusion_radius_per_type(`bool`,optional): Controls if one exclusion_radius per each espresso_type. Defaults to `False`.

        Returns:
            RE (`reaction_ensemble.ReactionEnsemble`): Instance of a reaction_ensemble.ReactionEnsemble object from the espressomd library.
            sucessful_reactions_labels(`lst`): Labels of the reactions set up by pyMBE.
            ionic_strength_res ('float'): Ionic strength of the reservoir (useful for calculating partition coefficients).
        """
        from espressomd import reaction_methods
        if exclusion_range is None:
            exclusion_range = max(self.get_radius_map().values())*2.0
        if pka_set is None:
            pka_set=self.get_pka_set()    
        self.check_pka_set(pka_set=pka_set)
        if use_exclusion_radius_per_type:
            exclusion_radius_per_type = self.get_radius_map()
        else:
            exclusion_radius_per_type = {}
        
        RE = reaction_methods.ReactionEnsemble(kT=self.kT.to('reduced_energy').magnitude,
                                                    exclusion_range=exclusion_range, 
                                                    seed=self.seed, 
                                                    exclusion_radius_per_type = exclusion_radius_per_type
                                                    )

        # Determine the concentrations of the various species in the reservoir and the equilibrium constants
        cH_res, cOH_res, cNa_res, cCl_res = self.determine_reservoir_concentrations(pH_res, c_salt_res, activity_coefficient)
        ionic_strength_res = 0.5*(cNa_res+cCl_res+cOH_res+cH_res)
        determined_activity_coefficient = activity_coefficient(ionic_strength_res)
        a_hydrogen = (10 ** (-pH_res) * self.units.mol/self.units.l).to('1/(N_A * reduced_length**3)')
        a_cation = (cH_res+cNa_res).to('1/(N_A * reduced_length**3)') * np.sqrt(determined_activity_coefficient)
        a_anion = (cH_res+cNa_res).to('1/(N_A * reduced_length**3)') * np.sqrt(determined_activity_coefficient)
        K_XX = a_cation * a_anion

        cation_es_type = self.df.loc[self.df['name']==cation_name].state_one.es_type.values[0]
        anion_es_type = self.df.loc[self.df['name']==anion_name].state_one.es_type.values[0]     
        cation_charge = self.df.loc[self.df['name']==cation_name].state_one.z.values[0]
        anion_charge = self.df.loc[self.df['name']==anion_name].state_one.z.values[0]     
        if cation_charge <= 0:
            raise ValueError('ERROR cation charge must be positive, charge ', cation_charge)
        if anion_charge >= 0:
            raise ValueError('ERROR anion charge must be negative, charge ', anion_charge)

        # Coupling to the reservoir: 0 = X+ + X-
        RE.add_reaction(
            gamma = K_XX.magnitude,
            reactant_types = [],
            reactant_coefficients = [],
            product_types = [ cation_es_type, anion_es_type ],
            product_coefficients = [ 1, 1 ],
            default_charges = {
                cation_es_type: cation_charge, 
                anion_es_type: anion_charge, 
            }
        )

        sucessful_reactions_labels=[]
        charge_number_map = self.get_charge_number_map()
        for name in pka_set.keys():
            if not _DFm._check_if_name_is_defined_in_df(name=name, df=self.df):
                logging.warning(f'The acid-base reaction of {name} has not been set up because its particle type is not defined in the dataframe.')
                continue

            Ka = 10**-pka_set[name]['pka_value'] * self.units.mol/self.units.l
            gamma_K_AX = Ka.to('1/(N_A * reduced_length**3)').magnitude * a_cation / a_hydrogen

            state_one_type = self.df.loc[self.df['name']==name].state_one.es_type.values[0]
            state_two_type = self.df.loc[self.df['name']==name].state_two.es_type.values[0]

            # Reaction in terms of small cation: HA = A + X+
            RE.add_reaction(gamma=gamma_K_AX.magnitude,
                            reactant_types=[state_one_type],
                            reactant_coefficients=[1],
                            product_types=[state_two_type, cation_es_type],
                            product_coefficients=[1, 1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            cation_es_type: cation_charge})

            # Reaction in terms of small anion: X- + HA = A
            RE.add_reaction(gamma=gamma_K_AX.magnitude / K_XX.magnitude,
                            reactant_types=[state_one_type, anion_es_type],
                            reactant_coefficients=[1, 1],
                            product_types=[state_two_type],
                            product_coefficients=[1],
                            default_charges={state_one_type: charge_number_map[state_one_type],
                            state_two_type: charge_number_map[state_two_type],
                            anion_es_type: anion_charge})

            sucessful_reactions_labels.append(name)
        return RE, sucessful_reactions_labels, ionic_strength_res

    def setup_lj_interactions(self, espresso_system, shift_potential=True, combining_rule='Lorentz-Berthelot'):
        """
        Sets up the Lennard-Jones (LJ) potential between all pairs of particle types with values for `sigma`, `offset`, and `epsilon` stored in `pymbe.df`.

        Args:
            espresso_system(`espressomd.system.System`): Instance of a system object from the espressomd library.
            shift_potential(`bool`, optional): If True, a shift will be automatically computed such that the potential is continuous at the cutoff radius. Otherwise, no shift will be applied. Defaults to True.
            combining_rule(`string`, optional): combining rule used to calculate `sigma` and `epsilon` for the potential between a pair of particles. Defaults to 'Lorentz-Berthelot'.
            warning(`bool`, optional): switch to activate/deactivate warning messages. Defaults to True.

        Note:
            - LJ interactions will only be set up between particles with defined values of `sigma` and `epsilon` in the pmb.df. 
            - Currently, the only `combining_rule` supported is Lorentz-Berthelot.
            - Check the documentation of ESPResSo for more info about the potential https://espressomd.github.io/doc4.2.0/inter_non-bonded.html

        """
        from itertools import combinations_with_replacement
        compulsory_parameters_in_df = ['sigma','epsilon']
        shift=0
        if shift_potential:
            shift="auto"
        # List which particles have sigma and epsilon values defined in pmb.df and which ones don't
        particles_types_with_LJ_parameters = []
        non_parametrized_labels= []
        for particle_type in self.get_type_map().values():
            check_list=[]
            for key in compulsory_parameters_in_df:
                value_in_df=self.find_value_from_es_type(es_type=particle_type,
                                                        column_name=key)
                check_list.append(pd.isna(value_in_df))
            if any(check_list):
                non_parametrized_labels.append(self.find_value_from_es_type(es_type=particle_type, 
                                                                            column_name='label'))
            else:
                particles_types_with_LJ_parameters.append(particle_type)
        # Set up LJ interactions between all particle types
        for type_pair in combinations_with_replacement(particles_types_with_LJ_parameters, 2): 
            particle_name1 = self.find_value_from_es_type(es_type=type_pair[0],
                                                        column_name="name")
            particle_name2 = self.find_value_from_es_type(es_type=type_pair[1],
                                                        column_name="name")
            lj_parameters= self.get_lj_parameters(particle_name1 = particle_name1,
                                                 particle_name2 = particle_name2,
                                                 combining_rule = combining_rule)
            
            # If one of the particle has sigma=0, no LJ interations are set up between that particle type and the others    
            if not lj_parameters:
                continue
            espresso_system.non_bonded_inter[type_pair[0],type_pair[1]].lennard_jones.set_params(epsilon = lj_parameters["epsilon"].to('reduced_energy').magnitude, 
                                                                                    sigma = lj_parameters["sigma"].to('reduced_length').magnitude, 
                                                                                    cutoff = lj_parameters["cutoff"].to('reduced_length').magnitude,
                                                                                    offset = lj_parameters["offset"].to("reduced_length").magnitude, 
                                                                                    shift = shift)                                                                                          
            index = len(self.df)
            label1 = self.find_value_from_es_type(es_type=type_pair[0], column_name="label")
            label2 = self.find_value_from_es_type(es_type=type_pair[1], column_name="label")
            self.df.at [index, 'name'] = f'LJ: {label1}-{label2}'
            lj_params=espresso_system.non_bonded_inter[type_pair[0], type_pair[1]].lennard_jones.get_params()

            _DFm._add_value_to_df(df = self.df,
                                  index = index,
                                  key = ('pmb_type',''),
                                  new_value = 'LennardJones')

            _DFm._add_value_to_df(df = self.df,
                                  index = index,
                                  key = ('parameters_of_the_potential',''),
                                  new_value = lj_params,
                                  non_standard_value = True)
        if non_parametrized_labels:
            logging.warning(f'The following particles do not have a defined value of sigma or epsilon in pmb.df: {non_parametrized_labels}. No LJ interaction has been added in ESPResSo for those particles.')
        return

    def write_pmb_df (self, filename):
        '''
        Writes the pyMBE dataframe into a csv file
        
        Args:
            filename(`str`): Path to the csv file 
        '''

        columns_with_list_or_dict = ['residue_list','side_chains', 'parameters_of_the_potential','sequence']
        df = self.df.copy(deep=True)
        for column_name in columns_with_list_or_dict:
            df[column_name] = df[column_name].apply(lambda x: json.dumps(x) if isinstance(x, (np.ndarray, tuple, list, dict)) or pd.notnull(x) else x)
        df['bond_object'] = df['bond_object'].apply(lambda x: f'{x.__class__.__name__}({json.dumps({**x.get_params(), "bond_id": x._bond_id})})' if pd.notnull(x) else x)
        df.fillna(pd.NA, inplace=True)
        df.to_csv(filename)
        return
