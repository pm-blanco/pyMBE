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


import os
import json
from pathlib import Path
from typing import Any, Dict
import pandas as pd
import logging

from pyMBE.storage.pint_quantity import PintQuantity
from pyMBE.storage.templates.particle import ParticleTemplate, ParticleStateTemplate
from pyMBE.storage.templates.residue import ResidueTemplate
from pyMBE.storage.templates.molecule import MoleculeTemplate
from pyMBE.storage.templates.bond import BondTemplate
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.instances.residue import ResidueInstance
from pyMBE.storage.instances.molecule import MoleculeInstance
from pyMBE.storage.instances.bond import BondInstance
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
from pyMBE.storage.templates.peptide import PeptideTemplate
from pyMBE.storage.instances.peptide import PeptideInstance
from pyMBE.storage.templates.protein import ProteinTemplate
from pyMBE.storage.instances.protein import ProteinInstance
from pyMBE.storage.templates.hydrogel import HydrogelTemplate, HydrogelNode, HydrogelChain
from pyMBE.storage.instances.hydrogel import HydrogelInstance
from pyMBE.storage.templates.nanoparticle import NanoparticleTemplate
from pyMBE.storage.instances.nanoparticle import NanoparticleInstance
from pyMBE.storage.templates.lj import LJInteractionTemplate

def _decode(s):
    """
    Decodes a JSON-like object or string.

    Handles various input types and converts them to a Python object.

    Args:
        s (Any): Input value to decode. Can be None, float('nan'), dict, list, number, or string.

    Returns:
        Any: 
            - None if input is None, NaN, empty string, or non-string unrecognized type.
            - Decoded Python object if input is a JSON string.
            - Original object if it is already a dict, list, int, or bool.
    """
    # None / pandas NA / nan handling
    if s is None:
        return None
    # pandas often gives float('nan') or numpy.nan
    if isinstance(s, float):
        # NaN -> None
        if pd.isna(s):
            return None
        return s
    # If already native Python container
    if isinstance(s, (dict, list, int, bool)):
        return s
    # Must be a string to parse JSON
    if not isinstance(s, str):
        return None
    s_str = s.strip()
    if s_str == "" or s_str.lower() == "nan":
        return None
    try:
        return json.loads(s_str)
    except Exception:
        # If it fails, try to interpret as plain string
        return s_str

def _encode(obj):
    """
    Encodes a Python object as a JSON string.

    Special handling for PintQuantity and Pydantic models.

    Args:
        obj (Any): Object to encode. Can be None, PintQuantity, Pydantic model, or standard Python object.

    Returns:
        str: JSON string representation of the object.
             Returns empty string for None.
    """
    if obj is None:
        return ""
    # PintQuantity dataclass (has to_dict)
    if isinstance(obj, PintQuantity):
        return json.dumps(obj.to_dict(), separators=(",", ":"), ensure_ascii=False)

    # If it's already a dict/list/scalar, json-dump it
    try:
        return json.dumps(obj, separators=(",", ":"), ensure_ascii=False)
    except TypeError:
        # Last resort: convert to string
        return json.dumps(str(obj), separators=(",", ":"), ensure_ascii=False)

def _load_database_csv(db, folder):
    """
    Loads CSV files from a folder into a database instance.

    This function populates the `templates`, `instances`, and `reactions` attributes
    of the provided database object in place. Supports various pyMBE types.

    Args:
        db (Manager): Database manager object to populate.
        folder (str or Path): Path to the folder containing CSV files.

    Return:
        (dict): metadata with additional information about the source of the information in the database.

    Notes:
        - PintQuantity objects are reconstructed from their dictionary representation.
        - Supports particle, residue, molecule, peptide, protein, bond, hydrogel,
          nanoparticle, and lj types.
    """
    folder = Path(folder)
    if not folder.exists():
        raise FileNotFoundError(f"Folder '{folder}' does not exist.")
    # target pmb types we support
    pyMBE_types = ["particle", 
                   "particle_state",
                   "residue", 
                   "molecule", 
                   "bond", 
                   "peptide",
                   "protein",
                   "hydrogel",
                   "nanoparticle",
                   "lj"]
    # TEMPLATES
    for pmb_type in pyMBE_types:
        csv_file = folder / f"templates_{pmb_type}.csv"
        if not csv_file.exists():
            continue
        df = pd.read_csv(csv_file, dtype=str).fillna("")
        templates: Dict[str, Any] = {}
        for _, row in df.iterrows():
            # row values are strings (or empty string)
            if pmb_type == "particle":
                sigma_d = _decode(row["sigma"])
                epsilon_d = _decode(row["epsilon"])
                cutoff_d = _decode(row["cutoff"])
                offset_d = _decode(row["offset"])
                sigma = PintQuantity.from_dict(sigma_d) if sigma_d is not None else None
                epsilon = PintQuantity.from_dict(epsilon_d) if epsilon_d is not None else None
                cutoff = PintQuantity.from_dict(cutoff_d) if cutoff_d is not None else None
                offset = PintQuantity.from_dict(offset_d) if offset_d is not None else None
                tpl = ParticleTemplate(name=row["name"],
                                       sigma=sigma,
                                       epsilon=epsilon,
                                       cutoff=cutoff,
                                       offset=offset,
                                       initial_state=row["initial_state"])
                templates[tpl.name] = tpl
            elif pmb_type == "particle_state":
                tpl = ParticleStateTemplate(particle_name=row["particle_name"],
                                            z=row["z"],
                                            name=row["name"],
                                            es_type=row["es_type"])
                templates[tpl.name] = tpl
            elif pmb_type == "residue":
                sc = _decode(row.get("side_chains", "")) or []
                if not isinstance(sc, list):
                    sc = list(sc)
                tpl = ResidueTemplate(name=row["name"],
                                      central_bead=row.get("central_bead", ""),
                                      side_chains=sc)
                templates[tpl.name] = tpl
            elif pmb_type == "molecule":
                rl = _decode(row.get("residue_list", "")) or []
                if not isinstance(rl, list):
                    rl = list(rl)
                tpl = MoleculeTemplate(name=row["name"],
                                       residue_list=rl)
                templates[tpl.name] = tpl
            elif pmb_type == "peptide":
                rl = _decode(row.get("residue_list", "")) or []
                if not isinstance(rl, list):
                    rl = list(rl)
                tpl = PeptideTemplate(name=row["name"],
                                      model=row.get("model", ""),
                                      residue_list=rl,
                                      sequence=row["sequence"])
                templates[tpl.name] = tpl
            elif pmb_type == "protein":
                rl = _decode(row.get("residue_list", "")) or []
                if not isinstance(rl, list):
                    rl = list(rl)
                tpl = ProteinTemplate(name=row["name"],
                                      model=row.get("model", ""),
                                      residue_list=rl,
                                      sequence=row["sequence"])
                templates[tpl.name] = tpl
            elif pmb_type == "bond":
                params_raw = _decode(row.get("parameters", "")) or {}
                parameters: Dict[str, Any] = {}
                # For the default bond case, map empty particle names to None
                particle_name1 = row.get("particle_name1", "") or ""
                particle_name2 = row.get("particle_name2", "") or ""
                for k, v in params_raw.items():
                    # if v is a dict, assume PintQuantity dict
                    if isinstance(v, dict) and {"magnitude", "units", "dimension"}.issubset(v.keys()):
                        parameters[k] = PintQuantity.from_dict(v)
                tpl = BondTemplate(name=row["name"],
                                   bond_type=row.get("bond_type", ""),
                                   particle_name1=None if particle_name1 == "" else particle_name1,
                                   particle_name2=None if particle_name2 == "" else particle_name2,
                                   parameters=parameters)
                templates[tpl.name] = tpl
            elif pmb_type == "hydrogel":
                node_map_raw = _decode(row.get("node_map", "")) or []
                chain_map_raw = _decode(row.get("chain_map", "")) or []
                
                node_map = [HydrogelNode(**n) for n in node_map_raw if isinstance(n, dict)]
                chain_map = [HydrogelChain(**c) for c in chain_map_raw if isinstance(c, dict)]
                tpl = HydrogelTemplate(name=row["name"],
                                       node_map=node_map,
                                       chain_map=chain_map)
                templates[tpl.name] = tpl
            elif pmb_type == "nanoparticle":
                surface_density_d = _decode(row.get("surface_density_of_sites", ""))
                surface_density = PintQuantity.from_dict(surface_density_d) if surface_density_d is not None else None
                secondary_site = row.get("secondary_site_particle_name", "") or None
                tpl = NanoparticleTemplate(name=row["name"],
                                           core_particle_name=row["core_particle_name"],
                                           surface_density_of_sites=surface_density,
                                           primary_site_particle_name=row["primary_site_particle_name"],
                                           fraction_primary_sites=float(row["fraction_primary_sites"]),
                                           number_of_patches_of_primary_sites=int(row["number_of_patches_of_primary_sites"]),
                                           secondary_site_particle_name=secondary_site)
                templates[tpl.name] = tpl
            elif pmb_type == "lj":
                sigma_d = _decode(row["sigma"])
                epsilon_d = _decode(row["epsilon"])
                cutoff_d = _decode(row["cutoff"])
                offset_d = _decode(row["offset"])
                state1 = row["state1"]
                state2 = row["state2"]
                shift_d = _decode(row.get("shift", ""))
                sigma = PintQuantity.from_dict(sigma_d) if sigma_d is not None else None
                epsilon = PintQuantity.from_dict(epsilon_d) if epsilon_d is not None else None
                cutoff = PintQuantity.from_dict(cutoff_d) if cutoff_d is not None else None
                offset = PintQuantity.from_dict(offset_d) if offset_d is not None else None
                if isinstance(shift_d, dict) and {"magnitude", "units", "dimension"}.issubset(shift_d.keys()):
                    shift = PintQuantity.from_dict(shift_d)
                else:
                    shift = shift_d  # could be "auto" or None
                tpl = LJInteractionTemplate(state1=state1,
                                            state2=state2,
                                            sigma=sigma,
                                            epsilon=epsilon,
                                            cutoff=cutoff,
                                            offset=offset,
                                            shift=shift)
                templates[tpl.name] = tpl
        db._templates[pmb_type] = templates

    # INSTANCES
    for pmb_type in pyMBE_types:
        csv_file = folder / f"instances_{pmb_type}.csv"
        if not csv_file.exists():
            continue
        df = pd.read_csv(csv_file, dtype=str).fillna("")

        instances: Dict[Any, Any] = {}

        for _, row in df.iterrows():
            if pmb_type == "particle":
                # some fields might be empty strings -> map to None
                residue_val = row.get("residue_id", "") or ""
                molecule_val = row.get("molecule_id", "") or ""
                assembly_val = row.get("assembly_id", "") or ""
                inst = ParticleInstance(name=row["name"],
                                        particle_id=int(row["particle_id"]),
                                        initial_state=row["initial_state"],
                                        residue_id=None if residue_val == "" else int(residue_val),
                                        molecule_id=None if molecule_val == "" else int(molecule_val),
                                        assembly_id=None if assembly_val == "" else int(assembly_val))
                instances[inst.particle_id] = inst
            elif pmb_type == "residue":
                mol_val = row.get("molecule_id", "") or ""
                assembly_val = row.get("assembly_id", "") or ""
                inst = ResidueInstance(name=row["name"],
                                       residue_id=int(row["residue_id"]),
                                       molecule_id=None if mol_val == "" else int(mol_val),
                                       assembly_id=None if assembly_val == "" else int(assembly_val))
                instances[inst.residue_id] = inst
            elif pmb_type == "molecule":
                assembly_val = row.get("assembly_id", "") or ""
                inst = MoleculeInstance(name=row["name"],
                                        molecule_id=int(row["molecule_id"]),
                                        assembly_id=None if assembly_val == "" else int(assembly_val))
                instances[inst.molecule_id] = inst
            elif pmb_type == "peptide":
                assembly_val = row.get("assembly_id", "") or ""
                inst = PeptideInstance(name=row["name"],
                                       molecule_id=int(row["molecule_id"]),
                                       assembly_id=None if assembly_val == "" else int(assembly_val))
                instances[inst.molecule_id] = inst
            elif pmb_type == "protein":
                assembly_val = row.get("assembly_id", "") or ""
                inst = ProteinInstance(name=row["name"],
                                       molecule_id=int(row["molecule_id"]),
                                       assembly_id=None if assembly_val == "" else int(assembly_val))
                instances[inst.molecule_id] = inst
            elif pmb_type == "bond":
                inst = BondInstance(name=row["name"],
                                    bond_id=int(row["bond_id"]),
                                    particle_id1=int(row["particle_id1"]),
                                    particle_id2=int(row["particle_id2"]))
                instances[inst.bond_id] = inst
            elif pmb_type == "hydrogel":
                inst = HydrogelInstance(name=row["name"],
                                        assembly_id=int(row["assembly_id"]))
                instances[inst.assembly_id] = inst
            elif pmb_type == "nanoparticle":
                molecule_val = row.get("molecule_id", "") or ""
                assembly_val = row.get("assembly_id", "") or ""
                inst = NanoparticleInstance(name=row["name"],
                                            molecule_id=int(molecule_val),
                                            assembly_id=None if assembly_val == "" else int(assembly_val))
                instances[inst.molecule_id] = inst
        db._instances[pmb_type] = instances

    # REACTIONS
    rx_file = folder / "reactions.csv"
    reactions: Dict[str, Reaction] = {}
    if rx_file.exists():
        df = pd.read_csv(rx_file, dtype=str).fillna("")
        for _, row in df.iterrows():
            participants_raw = _decode(row.get("participants", "")) or []
            participants = [ReactionParticipant(**p) for p in participants_raw]
            metadata = _decode(row.get("metadata", "")) or None
            rx = Reaction(name=row["name"],
                          participants=participants,
                          pK=float(row["pK"]) if (row.get("pK") not in (None, "", "nan")) else None,
                          reaction_type=row.get("reaction_type", None),
                          simulation_method=row.get("simulation_method", None),
                          metadata=metadata)
            db._reactions[rx.name] = rx

    # Metadata
    json_file = folder / "metadata.json"
    if json_file.exists():
        try:
            with open(json_file, "r", encoding="utf-8") as fh:
                metadata = json.load(fh)
            if not isinstance(metadata, dict):
                raise ValueError("metadata.json must contain a JSON object")
        except Exception as err:
            logging.warning(f"Failed to read metadata file '{json_file}': {err}. Metadata will be ignored.")
            metadata = {}
    else:
        metadata = {}
    return metadata

def _save_database_csv(db, folder):
    """
    Saves the database content into CSV files in a folder.

    This function serializes all templates, instances, and reactions.

    Args:
        db (Manager): Database object containing templates, instances, and reactions.
        folder (str or Path): Path to the folder where CSV files will be saved.

    """
    os.makedirs(folder, exist_ok=True)
    # TEMPLATES
    for pmb_type, tpl_dict in db._templates.items():
        rows = []
        for tpl in tpl_dict.values():
            # PARTICLE TEMPLATE: explicit custom encoding
            if pmb_type == "particle" and isinstance(tpl, ParticleTemplate):
                rows.append({"name": tpl.name,
                            "sigma": _encode(tpl.sigma),
                            "epsilon": _encode(tpl.epsilon),
                            "cutoff": _encode(tpl.cutoff),
                            "offset": _encode(tpl.offset),
                            "initial_state": tpl.initial_state})
            # RESIDUE TEMPLATE
            elif pmb_type == "residue" and isinstance(tpl, ResidueTemplate):
                rows.append({"name": tpl.name,
                            "central_bead": tpl.central_bead,
                            "side_chains": _encode(tpl.side_chains)})
            # MOLECULE TEMPLATE
            elif pmb_type == "molecule" and isinstance(tpl, MoleculeTemplate):
                rows.append({"name": tpl.name,
                             "residue_list": _encode(tpl.residue_list)})
            elif pmb_type == "peptide" and isinstance(tpl, PeptideTemplate):
                rows.append({"name": tpl.name,
                            "model": tpl.model,
                            "residue_list": _encode(tpl.residue_list),
                            "sequence": tpl.sequence})
            elif pmb_type == "protein" and isinstance(tpl, ProteinTemplate):
                rows.append({"name": tpl.name,
                            "model": tpl.model,
                            "residue_list": _encode(tpl.residue_list),
                            "sequence": tpl.sequence})
            # BOND TEMPLATE
            elif pmb_type == "bond" and isinstance(tpl, BondTemplate):
                # parameters: dict[str, scalar or PintQuantity]
                params_serial = {}
                for k, v in tpl.parameters.items():
                    if isinstance(v, PintQuantity):
                        params_serial[k] = v.to_dict()
                rows.append({"name": tpl.name,
                            "particle_name1": tpl.particle_name1,
                            "particle_name2": tpl.particle_name2,
                            "bond_type": tpl.bond_type,
                            "parameters": _encode(params_serial)})
            # HYDROGEL TEMPLATE
            elif pmb_type == "hydrogel" and isinstance(tpl, HydrogelTemplate):
                rows.append({"name": tpl.name,
                            "node_map": _encode([node.dict() for node in tpl.node_map]),
                            "chain_map": _encode([chain.dict() for chain in tpl.chain_map])})
            elif pmb_type == "nanoparticle" and isinstance(tpl, NanoparticleTemplate):
                rows.append({"name": tpl.name,
                            "core_particle_name": tpl.core_particle_name,
                            "surface_density_of_sites": _encode(tpl.surface_density_of_sites),
                            "primary_site_particle_name": tpl.primary_site_particle_name,
                            "fraction_primary_sites": tpl.fraction_primary_sites,
                            "number_of_patches_of_primary_sites": tpl.number_of_patches_of_primary_sites,
                            "secondary_site_particle_name": tpl.secondary_site_particle_name if tpl.secondary_site_particle_name is not None else ""})
            # LJ TEMPLATE
            elif pmb_type == "lj" and isinstance(tpl, LJInteractionTemplate):
                rows.append({"name":   tpl.name,
                            "state1": tpl.state1,
                            "state2": tpl.state2,
                            "sigma":  _encode(tpl.sigma),
                            "epsilon":_encode(tpl.epsilon),
                            "cutoff": _encode(tpl.cutoff),
                            "offset": _encode(tpl.offset),
                            "shift": _encode(tpl.shift)})
            else:
                # Generic fallback: try dict()
                try:
                    rows.append(tpl.dict())
                except Exception:
                    rows.append({"name": getattr(tpl, "name", None)})

        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(folder, f"templates_{pmb_type}.csv"), index=False)

    # INSTANCES
    for pmb_type, inst_dict in db._instances.items():
        rows = []
        for inst in inst_dict.values():
            if pmb_type == "particle" and isinstance(inst, ParticleInstance):
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "particle_id": int(inst.particle_id),
                            "initial_state": inst.initial_state,
                            "residue_id": int(inst.residue_id) if inst.residue_id is not None else "",
                            "molecule_id": int(inst.molecule_id) if inst.molecule_id is not None else "",
                            "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else ""})
            elif pmb_type == "residue" and isinstance(inst, ResidueInstance):
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "residue_id": int(inst.residue_id),
                            "molecule_id": int(inst.molecule_id) if inst.molecule_id is not None else "",
                            "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else ""})
            elif pmb_type == "molecule" and isinstance(inst, MoleculeInstance):
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "molecule_id": int(inst.molecule_id),
                            "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else ""})
            elif pmb_type == "peptide" and isinstance(inst, PeptideInstance):
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "molecule_id": int(inst.molecule_id),
                            "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else ""})
            elif pmb_type == "protein" and isinstance(inst, ProteinInstance):
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "molecule_id": int(inst.molecule_id),
                            "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else ""})
            elif pmb_type == "bond" and isinstance(inst, BondInstance):
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "bond_id": int(inst.bond_id),
                            "particle_id1": int(inst.particle_id1),
                            "particle_id2": int(inst.particle_id2)})
            elif pmb_type == "hydrogel" and isinstance(inst, HydrogelInstance):
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "assembly_id": int(inst.assembly_id)})
            elif pmb_type == "nanoparticle" and isinstance(inst, NanoparticleInstance):
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "molecule_id": int(inst.molecule_id),
                            "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else ""})
            else:
                # fallback to dict
                try:
                    rows.append(inst.dict())
                except Exception:
                    rows.append({"name": getattr(inst, "name", None)})

        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(folder, f"instances_{pmb_type}.csv"), index=False)

    # REACTIONS
    rows = []
    for rx in db._reactions.values():
        rows.append({"name": rx.name,
                    "participants": _encode([p.dict() for p in rx.participants]),
                    "pK": rx.pK if hasattr(rx, "pK") else None,
                    "reaction_type": rx.reaction_type,
                    "metadata": _encode(rx.metadata) if getattr(rx, "metadata", None) is not None else ""})
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(folder, "reactions.csv"), index=False)
