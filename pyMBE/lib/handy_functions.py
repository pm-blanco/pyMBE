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

import re
import numpy as np
import scipy

def calculate_initial_bond_length(bond_parameters, bond_type, lj_parameters):
    """
    Calculate an initial bond length for molecule setup.

    Args:
        bond_parameters ('dict'):
            Parameters defining the bonded interaction (e.g. equilibrium
            distance, force constant), as required by the selected
            ``bond_type``.

        bond_type ('str'):
            Label identifying the bonded potential used to connect the
            particles (e.g. ``"harmonic"``).

        lj_parameters ('dict'):
            Parameters of the Lennard-Jones interaction between the bonded
            particles. Expected entries include ``epsilon``, ``sigma``,
            ``cutoff``, and optionally ``offset``, typically given as
            ``pint.Quantity`` objects.

    Returns:
        ('pint.Quantity'):
            Initial bond length resulting from the minimum of the bonded
            and Lennard-Jones interactions.

    Notes:
        - This function is intended for geometry initialization and does not
          affect the interaction parameters used during the simulation.
        - The exact interpretation of ``bond_parameters`` depends on
          ``bond_type``.
    """  
    def truncated_lj_potential(x, epsilon, sigma, cutoff,offset):
        if x>cutoff:
            return 0.0
        else:
            return 4*epsilon*((sigma/(x-offset))**12-(sigma/(x-offset))**6) - 4*epsilon*((sigma/cutoff)**12-(sigma/cutoff)**6)
    epsilon=lj_parameters["epsilon"].m_as("reduced_energy")
    sigma=lj_parameters["sigma"].m_as("reduced_length")
    cutoff=lj_parameters["cutoff"].m_as("reduced_length")
    offset=lj_parameters["offset"].m_as("reduced_length")
    if bond_type == "harmonic":
        r_0 = bond_parameters['r_0'].m_as("reduced_length")
        k = bond_parameters['k'].m_as("reduced_energy/reduced_length**2")
        l0 = scipy.optimize.minimize(lambda x: 0.5*k*(x-r_0)**2 + truncated_lj_potential(x, epsilon, sigma, cutoff, offset), x0=r_0).x
    elif bond_type == "FENE":
        r_0 = bond_parameters['r_0'].m_as("reduced_length")
        k = bond_parameters['k'].m_as("reduced_energy/reduced_length**2")
        d_r_max = bond_parameters['d_r_max'].m_as("reduced_length")
        l0 = scipy.optimize.minimize(lambda x: -0.5*k*(d_r_max**2)*np.log(1-((x-r_0)/d_r_max)**2) + truncated_lj_potential(x, epsilon, sigma, cutoff,offset), x0=1.0).x
    return l0

def check_aminoacid_key(key):
    """
    Checks if `key` corresponds to a valid aminoacid letter code.

    Args:
        key (`str`): 
            key to be checked.

    Returns:
        (`bool`): 
            True if `key` is a valid aminoacid letter code, False otherwise.
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
                    'C', #'CYS'
                    "n", # n terminus
                    "c", # c terminus
                    ] 
    if key in valid_AA_keys:
        return True
    else:
        return False

def check_if_metal_ion(key):
    """
    Checks if `key` corresponds to a label of a supported metal ion.

    Args:
        key(`str`): 
            key to be checked

    Returns:
        (`bool`): 
            True if `key`  is a supported metal ion, False otherwise.
    """
    if key in get_metal_ions_charge_number_map().keys():
        return True
    else:
        return False

def define_protein_AA_particles(topology_dict, pmb, pka_set,  lj_setup_mode="wca"):
    """
    Defines particle templates in pyMBE for all unique residue/atom types appearing
    in a protein topology dictionary.

    Args:
        topology_dict ('dict'):
            Dictionary defining the structure of a protein.

        pmb ('pyMBE.pymbe_library'):
            Instance of the pyMBE library.

        pka_set ('dict'):
                Set of pka_values for the protein aminoacids and their corresponding acidities

        lj_setup_mode ('str', optional):
            Determines how Lennard-Jones parameters are assigned. Defaults to `"wca"`.           

    Notes:
        - Particle names are extracted by stripping trailing digits
          (e.g., `"ALA1"` → `"ALA"`).
        - For metal ions (identified via `check_if_metal_ion()`), the correct
          ionic charge is retrieved from the metal-ion charge map.
        - The Lennard-Jones offset is computed as:
                offset = 2 * radius - sigma
    """
    valid_lj_setups = ["wca"]
    if lj_setup_mode not in valid_lj_setups:
        raise ValueError('Invalid key for the lj setup, supported setup modes are {valid_lj_setups}')
    if lj_setup_mode == "wca":
        sigma = 1*pmb.units.Quantity("reduced_length")
        epsilon = 1*pmb.units.Quantity("reduced_energy")
    part_dict={}
    metal_ions_charge_number_map=get_metal_ions_charge_number_map()
    defined_particles=[]
    for particle in topology_dict.keys():
        particle_name = re.split(r'\d+', particle)[0] 
        if particle_name not in defined_particles:
            part_dict = {"name" : particle_name}
            if lj_setup_mode == "wca":
                part_dict["sigma"] = sigma
                part_dict["offset"]= topology_dict[particle]['radius']*2-sigma
                part_dict["epsilon"] = epsilon
            if particle_name in pka_set.keys():
                part_dict["acidity"] = pka_set[particle_name]["acidity"]
            else:
                if check_if_metal_ion(key=particle_name):
                    z=metal_ions_charge_number_map[particle_name]
                else:
                    z=0
                part_dict["z"]=z
        if particle_name not in defined_particles:
            pmb.define_particle(**part_dict)
            defined_particles.append(particle_name)

def define_protein_AA_residues(sequence, model, pmb):
    """
    Define residue templates in the pyMBE database for a protein topology dict.

    Args:
        sequence ('str'):
                Protein sequence, following the one letter amino acid convention.
               
        model ('str'):
            Coarse-grained representation to use. Supported options:
                - `"1beadAA"`
                - `"2beadAA"`

        pmb ('pyMBE.pymbe_library'):
            Instance of the pyMBE library.
    Return:
        ('list of str'): 
            List of the defined residue names

    Notes:
        - Supported models:
            - `"1beadAA"`: Each amino acid is represented by a single bead.  
                The central bead is the amino-acid name itself, and no side chains are used.
            - `"2beadAA"`: Each amino acid is represented by two beads, except for terminal or special residues:
                * `"c"`, `"n"`, and `"G"` (glycine) are treated as single-bead residues.
                * All other residues use `"CA"` (central bead) plus one side-chain bead named after the amino acid.

        - Residue names are constructed as `"AA-<residue>"`, e.g., `"AA-A"`, `"AA-L"`.
    """
    residue_list = []
    for item in sequence:
        if model == '1beadAA':
            central_bead = item
            side_chains = []
        elif model == '2beadAA':
            if item in ['c','n', 'G']: 
                central_bead = item
                side_chains = []
            else:
                central_bead = 'CA'              
                side_chains = [item]
        residue_name='AA-'+item
        if residue_name not in residue_list:   
            pmb.define_residue(name = residue_name, 
                                central_bead = central_bead,
                                side_chains = side_chains)              
        residue_list.append(residue_name)
    return residue_list

def define_peptide_AA_residues(sequence,model, pmb):
    """
    Define residue templates in the pyMBE database for a given model.

    Args:
        sequence ('list of str'):
            Ordered amino-acid sequence of the peptide or protein. Each element must
            be a residue identifier compatible with the selected model.

        model ('str'):
            Coarse-grained representation to use. Supported options:
                - `"1beadAA"`
                - `"2beadAA"`

        pmb ('pyMBE.pymbe_library'):
            Instance of the pyMBE library.

    Notes:
        - Supported models:
            - `"1beadAA"`: Each amino acid is represented by a single bead.  
                The central bead is the amino-acid name itself, and no side chains are used.
            - `"2beadAA"`: Each amino acid is represented by two beads, except for terminal or special residues:
                * `"c"`, `"n"`, and `"G"` (glycine) are treated as single-bead residues.
                * All other residues use `"CA"` (central bead) plus one side-chain bead named after the amino acid.

        - Residue names are constructed as `"AA-<residue>"`, e.g., `"AA-A"`, `"AA-L"`.
    """
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
        residue_name='AA-'+residue_name
        if "residue" in pmb.db._templates:
            if residue_name not in pmb.db._templates["residue"]:   
                pmb.define_residue(name = residue_name, 
                                    central_bead = central_bead,
                                    side_chains = side_chains)
        else:
            pmb.define_residue(name = residue_name, 
                                    central_bead = central_bead,
                                    side_chains = side_chains)
def get_residues_from_topology_dict(topology_dict, model):
    """
    Groups beads from a topology dictionary into residues and assigns residue names.

    Args:
        topology_dict ('dict'):
            Dictionary describing the molecular topology, where keys are bead
            identifiers (e.g. "CA12", "SC12") that encode both residue type and
            residue index.

        model ('str'):
            Protein model identifier. Supported values are:
            - `"1beadAA"`: single-bead-per-amino-acid model.
            - `"2beadAA"`: two-bead-per-amino-acid model, where CA beads are excluded
              from residue name assignment.

    Returns:
        ('dict'):
            Dictionary mapping residue indices (as strings) to residue data:
            {
                resid: {
                    "beads": [bead_id1, bead_id2, ...],
                    "resname": residue_name
                },
                ...
            }

    Notes:
        - Bead identifiers are parsed by separating alphabetic prefixes
          (residue or bead type) from numeric residue indices.
        - For the `"2beadAA"` model, beads named `"CA"` are excluded when
          determining the residue name.
        - Residues that only contain CA beads (i.e., no side-chain beads)
          are assigned the residue name `"G"` (glycine).
        - Residue indices are returned as strings, consistent with the parsed
          bead identifiers.
    """
    if model not in {"1beadAA", "2beadAA"}:
        raise ValueError(f"Unknown protein model '{model}'")    
    if model == "1beadAA":
        excluded_residue_names = []
    elif model == "2beadAA":
        excluded_residue_names = ["CA"]
    # GROUP BEADS BY RESIDUE
    residues = {}
    for bead_id in topology_dict.keys():
        # extract prefix and index number
        prefix = re.split(r'\d+', bead_id)[0]         
        index_match = re.findall(r'\d+', bead_id)
        if not index_match:
            raise ValueError(f"Topology key '{bead_id}' does not contain a residue index.")
        resid = index_match[0]
        if resid not in residues:
            residues[resid] = {"beads": []}
        residues[resid]["beads"].append(bead_id)
        if prefix not in excluded_residue_names:
            residues[resid]["resname"] = prefix
    
    # Assign name to glycine residues (only with CA beads)
    for bead_id in residues:
        if "resname" not in residues[bead_id]:
            residues[bead_id]["resname"] = "G"
    return residues

def get_metal_ions_charge_number_map():
    """
    Gets a map with the charge numbers of all the metal ions supported.

    Returns:
        ('dict'): 
            Has the structure {"metal_name": metal_charge_number}

    """
    metal_charge_number_map = {"Ca": 2}
    return metal_charge_number_map

def protein_sequence_parser(sequence):
    """
    Parses `sequence` to the one letter code for amino acids.
    
    Args:
        sequence(`str` or `lst`): 
            Sequence of the amino acid. 

    Returns:
        (`lst`): `
            sequence` using the one letter code.
    
    Notes:
        - Accepted formats for `sequence` are:
            - `lst` with one letter or three letter code of each aminoacid in each element
            - `str` with the sequence using the one letter code
            - `str` with the squence using the three letter code, each aminoacid must be separated by a hyphen "-"
    """
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