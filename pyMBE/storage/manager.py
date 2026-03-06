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

import pandas as pd
from collections import defaultdict
from typing import Dict, Any
from pyMBE.storage.templates.particle import ParticleTemplate
from pyMBE.storage.templates.residue import ResidueTemplate
from pyMBE.storage.templates.molecule import MoleculeTemplate
from pyMBE.storage.templates.bond import BondTemplate
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.instances.residue import ResidueInstance
from pyMBE.storage.instances.molecule import MoleculeInstance
from pyMBE.storage.instances.bond import BondInstance
from pyMBE.storage.reactions.reaction import Reaction
from pyMBE.storage.templates.peptide import PeptideTemplate
from pyMBE.storage.instances.peptide import PeptideInstance
from pyMBE.storage.templates.protein import ProteinTemplate
from pyMBE.storage.instances.protein import ProteinInstance
from pyMBE.storage.templates.hydrogel import HydrogelTemplate
from pyMBE.storage.instances.hydrogel import HydrogelInstance
from pyMBE.storage.templates.nanoparticle import NanoparticleTemplate
from pyMBE.storage.instances.nanoparticle import NanoparticleInstance
from pyMBE.storage.templates.lj import LJInteractionTemplate
from pyMBE.storage.pint_quantity import PintQuantity

TemplateType = Any  # union of template classes (ParticleTemplate, ResidueTemplate, ...)
InstanceType = Any  # union of instance classes (ParticleInstance, ResidueInstance, ...)

class Manager:
    """
    The canonical database manager for pyMBE.

    Attributes:
    
        _units ('pint.UnitRegistry'): 
            Pint unit registry used to reconstruct physical quantities from storage.

        _templates ('dict[str, dict[str, TemplateType]]'):
            Templates indexed by type and name.
        
        _instances ('dict[str, dict[int, InstanceType]]'):
            Instances indexed by type and id.

        _reactions ('dict[str, Reaction]'):
            Chemical reactions keyed by reaction name.

        _molecule_like_types ('list'):
            List of pyMBE object types that belong to the 'molecule' category in the pyMBE hierarchy.

        _assembly_like_types ('list'):
            List of pyMBE object types that belong to the 'assembly' category in the pyMBE hierarchy.

        _pmb_types ('list'):
            List of all supported pyMBE object types.

        espresso_bond_instances ('dict[int,espressomd.interactions.BondedInteraction]'):
            List of active instances of bonded interactions from ESPResSo.
    """

    def __init__(self,units):
        """
        Initialize an empty structured database.

        Args:
            units ('pint.UnitRegistry'): 
                Pint unit registry used to reconstruct physical quantities from storage.
        """
        self._units = units
        self._templates: Dict[str, Dict[str, TemplateType]] = {}
        self._instances: Dict[str, Dict[int, InstanceType]] = {}
        self._reactions: Dict[str, Reaction] = {}
        self._molecule_like_types = ["molecule",
                                     "peptide",
                                     "protein",
                                     "nanoparticle"]
        self._assembly_like_types = ["hydrogel"]
        self._pmb_types =  ["particle", "residue"] + self._molecule_like_types + self._assembly_like_types
        self.espresso_bond_instances= {}

    def _collect_particle_templates(self, name, pmb_type):
        """
        Recursively collect particle template names reachable from a given
        template in the hierarchy, accounting for their multiplicity.

        Args:
            name ('str'):
                Name of the template being processed.

            pmb_type ('str'):
                Type of the template. 

        Returns:
            ('collections.defaultdict[str, int]'):
                A mapping from particle template names to their occurrence counts
                in the hierarchy reachable from the given template.

        Notes:
            - If ``pmb_type == "particle"``, the particle itself is counted once.
            - If ``pmb_type == "particle_state"``, the state is resolved to its
            parent particle template, which is counted once.
            - Residue templates contribute their central bead and all side-chain
            particles.
            - Molecule-like templates contribute the particles from all residues
            in their ``residue_list``.
        """
        counts = defaultdict(int)
        if pmb_type == "particle":
            counts[name] += 1
            return counts
        if pmb_type == "particle_state":
            particle_name = self.get_template(name=name,pmb_type=pmb_type).particle_name
            counts[particle_name] += 1
            return counts
        if pmb_type == "residue":
            tpl = self.get_template(name=name, pmb_type="residue")
            # central bead is always a particle
            sub = self._collect_particle_templates(name=tpl.central_bead,
                                                   pmb_type="particle")
            for k, v in sub.items():
                counts[k] += v
            # side chains can be particles OR residues
            for sc_name in tpl.side_chains:
                if sc_name in self._templates.get("particle", {}):
                    sc_type = "particle"
                elif sc_name in self._templates.get("residue", {}):
                    sc_type = "residue"
                sub = self._collect_particle_templates(name=sc_name,
                                                       pmb_type=sc_type)
                for k, v in sub.items():
                    counts[k] += v
            return counts
        if pmb_type in self._molecule_like_types:
            tpl = self.get_template(name=name, 
                                    pmb_type=pmb_type)
            if pmb_type == "nanoparticle":
                counts[tpl.core_particle_name] += 1
                counts[tpl.primary_site_particle_name] += 1
                if tpl.secondary_site_particle_name is not None:
                    counts[tpl.secondary_site_particle_name] += 1
                return counts
            for res_name in tpl.residue_list:
                sub = self._collect_particle_templates(name=res_name,
                                                       pmb_type="residue")
                for k, v in sub.items():
                    counts[k] += v
            return counts
        raise NotImplementedError(f"Method not implemented for pmb_type='{pmb_type}'")

    def _delete_bonds_of_particle(self, pid):
        """
        Delete all bond instances involving a given particle instance.

        Args:
            pid ('int'): 
                The particle ID whose associated bonds should be deleted.

        Notes:
            - If no `"bond"` instances are present in the database, the method
            exits immediately.
            - This method does not raise errors if no bonds involve the particle.
            - It is intended for internal use by cascade-deletion routines.
        """
        if "bond" not in self._instances:
            return
        bonds_to_delete = [b_id for b_id, b in list(self._instances["bond"].items()) if b.particle_id1 == pid or b.particle_id2 == pid]
        for b_id in bonds_to_delete:
            del self._instances["bond"][b_id]
        if "bond" in self._instances and not self._instances["bond"]:
            del self._instances["bond"]

    def _find_instance_ids_by_attribute(self, pmb_type, attribute, value):
        """
        Return a list of instance IDs for a given pmb_type where a given attribute
        matches the requested value.

        Args:
            pmb_type ('str'): 
                The pyMBE type to search within.
                
            attribute ('str'): 
                The attribute name to match on (e.g. "residue_id", "molecule_id").

            value ('Any'): 
                The attribute value to match.

        Returns:
            ('List[int]'): 
                IDs of matching instances.
        """
        if pmb_type not in self._instances:
            return []
        results = []
        for inst_id, inst in self._instances[pmb_type].items():
            if hasattr(inst, attribute) and getattr(inst, attribute) == value:
                results.append(inst_id)
        return results

    

    def _find_template_types(self, name):
        """
        Return all pyMBE template categories that contain a template
        with a given name.

        Args:
            name ('str'):
                The template name to search for.

        Returns:
            ('list[str]'):
                A list of PMB types (e.g., ``["particle", "residue"]``) in
                which a template named ``name`` exists. The list is empty if
                no such template is found.
        """
        found = []
        for pmb_type, group in self._templates.items():
            if name in group:
                found.append(pmb_type)
        return found


    def _get_instances_df(self, pmb_type):
        """
        Returns a DataFrame containing all instance objects of a given pyMBE type.

        Args:
            pmb_type ('str'):
                The instance type to query. Must be a key in
                `self._instances`, such as `"particle"` or `"residue"`.

        Returns:
            ('pandas.DataFrame'):
                A DataFrame where each row corresponds to one registered
                instance of the specified PMB type. If no instances exist,
                an empty DataFrame is returned.

        Notes:
            - Missing integer identifiers (e.g., `residue_id`) are stored as
            `pandas.NA` to ensure proper nullable integer handling.
            - Particle and residue instances receive custom row structures;
            all other instance types use direct model dumps.
        """
        rows = []
        if pmb_type not in self._instances:
            return pd.DataFrame(rows)
        for inst in self._instances[pmb_type].values():
            if pmb_type == "particle":
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "particle_id": inst.particle_id,
                            "initial_state": inst.initial_state,
                            "residue_id": int(inst.residue_id) if inst.residue_id is not None else pd.NA,
                            "molecule_id": int(inst.molecule_id) if inst.molecule_id is not None else pd.NA,
                            "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else pd.NA })
            elif pmb_type == "residue":
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "residue_id": inst.residue_id,
                            "molecule_id": int(inst.molecule_id) if inst.molecule_id is not None else pd.NA,
                            "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else pd.NA})
            elif pmb_type in ["molecule","peptide","protein"]:
                rows.append({"pmb_type": pmb_type,
                            "name": inst.name,
                            "molecule_id": inst.molecule_id,
                            "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else pd.NA})
            else:
                # Generic representation for other types
                rows.append(inst.dict())
        return pd.DataFrame(rows)

    def _get_reactions_df(self):
        """
        Returns a DataFrame summarizing all registered chemical reactions.

        Returns:
            ('pandas.DataFrame'):
                A DataFrame where each row corresponds to one reaction.

        Notes:
            - Participant objects are expected to expose ``state_name`` and
            ``coefficient`` attributes.
            - Stoichiometry is stored as a single dictionary per row to allow
            flexible downstream manipulation.
        """  
        rows = []
        for r in self._reactions.values():
            stoich = {f"{p.state_name}": p.coefficient  for p in r.participants}
            rows.append({"reaction": r.name,
                        "stoichiometry": stoich,
                        "pK": r.pK,
                        "reaction_type": r.reaction_type,
                        "metadata": r.metadata,
                        "simulation_method": r.simulation_method})
        return pd.DataFrame(rows)

    def _get_templates_df(self, pmb_type):
        """
        Returns a DataFrame containing all template definitions of a PMB type.

        Args:
            pmb_type ('str'):
                The template type to query, e.g. `"particle"`, `"residue"`,
                `"molecule"`.

        Returns:
            ('pandas.DataFrame'):
                A DataFrame representing all templates of the given type.
                Particle templates expand to multiple rows, one per state.
                Empty DataFrame if no templates for that type exist.

        Notes:
            - Unit-bearing fields are converted to plain quantities through
            ``to_quantity(self._units)`` to maintain consistent I/O.
        """
        rows = []
        if pmb_type not in self._templates:
            return pd.DataFrame(rows)
        for tpl in self._templates[pmb_type].values():
            if pmb_type == "particle":
                rows.append({"pmb_type": tpl.pmb_type,
                             "name": tpl.name,
                             "sigma": tpl.sigma.to_quantity(self._units),
                             "epsilon": tpl.epsilon.to_quantity(self._units),
                             "cutoff": tpl.cutoff.to_quantity(self._units),
                             "offset": tpl.offset.to_quantity(self._units),
                             "initial_state": tpl.initial_state})  
            elif pmb_type == "lj":
                shift = tpl.shift
                if isinstance(shift, PintQuantity):
                    shift = tpl.shift.to_quantity(self._units)
                rows.append({"pmb_type": tpl.pmb_type,
                             "name": tpl.name,
                             "state1": tpl.state1,
                             "state2": tpl.state2,
                             "sigma": tpl.sigma.to_quantity(self._units),
                             "epsilon": tpl.epsilon.to_quantity(self._units),
                             "cutoff": tpl.cutoff.to_quantity(self._units),
                             "offset": tpl.offset.to_quantity(self._units),
                             "shift": shift})
            elif pmb_type == "bond":
                parameters = {}
                for key in tpl.parameters.keys():
                    parameters[key] = tpl.parameters[key].to_quantity(self._units)
                rows.append({"pmb_type": tpl.pmb_type,
                             "name": tpl.name,
                             "bond_type": tpl.bond_type,
                             "particle_name1": tpl.particle_name1,
                             "particle_name2": tpl.particle_name2,
                             "parameters": parameters})
            elif pmb_type == "nanoparticle":
                rows.append({"pmb_type": tpl.pmb_type,
                             "name": tpl.name,
                             "core_particle_name": tpl.core_particle_name,
                             "surface_density_of_sites": tpl.surface_density_of_sites.to_quantity(self._units),
                             "primary_site_particle_name": tpl.primary_site_particle_name,
                             "fraction_primary_sites": tpl.fraction_primary_sites,
                             "number_of_patches_of_primary_sites": tpl.number_of_patches_of_primary_sites,
                             "secondary_site_particle_name": tpl.secondary_site_particle_name})
            else:
                # Generic representation for other types
                rows.append(tpl.dict())
        return pd.DataFrame(rows)

    def _has_instance(self, pmb_type, instance_id):
        """
        Check whether an instance with a given ID exists under a specific pyMBE type.

        Args:
            pmb_type ('str'):
                The instance category to search in. 

            instance_id ('int'):
                The unique identifier of the instance.

        Returns:
            ('bool'):
                ``True`` if the instance exists in the given category,
                ``False`` otherwise.
        """
        if pmb_type not in self._instances:
            raise ValueError(f"Instance type '{pmb_type}' not found in the database.")
        return instance_id in self._instances[pmb_type]

    def _has_template(self, pmb_type, name):
        """
        Check whether a template with a given name exists within a specific pyMBE type.

        Args:
            pmb_type ('str'):
                The template category to search in (e.g. ``"particle"``,
                ``"bond"``, ``"molecule"``, ``"lj"``, etc.).

            name ('str'):
                The template name to check for.

        Returns:
            ('bool'):
                ``True`` if a template named ``name`` exists under ``pmb_type``;
                ``False`` otherwise.
        """
        if pmb_type not in self._templates:
            raise ValueError(f"Template type '{pmb_type}' not found in the database.")
        template_in_db = name in self._templates.get(pmb_type, {})
        return template_in_db

    def _register_instance(self, instance):
        """
        Register an instance of a pyMBE object.

        Args:
            instance: Any instance conforming to the pyMBE instance models.
        """
        # infer pmb_type from instance class
        if isinstance(instance, ParticleInstance):
            pmb_type = "particle"
            iid = instance.particle_id
        elif isinstance(instance, ResidueInstance):
            pmb_type = "residue"
            iid = instance.residue_id
        elif isinstance(instance, MoleculeInstance):
            pmb_type = "molecule"
            iid = instance.molecule_id
        elif isinstance(instance, PeptideInstance):
            pmb_type = "peptide"
            iid = instance.molecule_id
        elif isinstance(instance, ProteinInstance):
            pmb_type = "protein"
            iid = instance.molecule_id
        elif isinstance(instance, BondInstance):
            pmb_type = "bond"
            iid = instance.bond_id
        elif isinstance(instance, HydrogelInstance):
            pmb_type = "hydrogel"
            iid = instance.assembly_id
        elif isinstance(instance, NanoparticleInstance):
            pmb_type = "nanoparticle"
            iid = instance.molecule_id
        else:
            raise TypeError("Unsupported instance type")
        self._instances.setdefault(pmb_type, {})
        if iid in self._instances[pmb_type]:
            raise ValueError(f"Instance id {iid} already exists in type '{pmb_type}'")
        # validate template exists
        if instance.name not in self._templates.get(pmb_type, {}):
            raise ValueError(f"Template '{instance.name}' not found for type '{pmb_type}'")
        self._instances[pmb_type][iid] = instance

    def _register_reaction(self, reaction):
        """
        Register a chemical or physical reaction.

        Args:
            reaction ('Reaction'): 
                Reaction template from the pyMBE database.
        """
        if reaction.name in self._reactions:
            raise ValueError(f"Reaction '{reaction.name}' already exists.")

        self._reactions[reaction.name] = reaction

    def _register_template(self, template):
        """
        Register a template.

        Args:
            template ('TemplateType'): 
                Any template object conforming to the pyMBE template models.

        """
        pmb_type = getattr(template, "pmb_type", None)
        self._templates.setdefault(pmb_type, {})
        if template.name in self._templates[pmb_type]:
            raise ValueError(f"Template '{template.name}' exists in '{pmb_type}'")

        self._templates[pmb_type][template.name] = template

    def _update_instance(self, instance_id, pmb_type, attribute, value):
        """
        Updates a single attribute of a registered instance.
        Only a restricted set of attributes is allowed for each PMB type,
        ensuring database consistency. 

        Args:
            instance_id ('int'):
                Unique identifier of the instance to update. 

            pmb_type ('str'):
                Instance category, such as ``"particle"`` or ``"residue"``.

            attribute ('str'):
                Name of the field to update. 
                
            value ('Any'):
                New value to assign to the specified attribute.

        Notes:
            - Allowed updates:
                * ``particle``: ``initial_state``, ``residue_id``, ``molecule_id``, ``assembly_id`` 
                * ``residue``: ``molecule_id``, ``assembly_id``
                * ``molecule``: ``assembly_id``
                * All other types: no attribute updates allowed.
            - The method replaces the instance with a new Pydantic model
            using ``copy(update=...)`` to maintain immutability and
            avoid partial mutations of internal state.
        """
        if instance_id not in self._instances[pmb_type]:
            raise ValueError(f"Instance '{instance_id}' not found for type '{pmb_type}' in the pyMBE database.")                                
        if pmb_type == "particle":
            allowed = ["initial_state", "residue_id", "molecule_id", "assembly_id"]
        elif pmb_type == "residue":
            allowed = ["molecule_id", "assembly_id"]
        elif pmb_type in self._molecule_like_types:
            allowed = ["assembly_id"]
        else:
            allowed = [None]  # No attributes allowed for other types        
        if attribute not in allowed:
            raise ValueError(f"Attribute '{attribute}' not allowed for {pmb_type}. Allowed attributes: {allowed}")
        self._instances[pmb_type][instance_id] = self._instances[pmb_type][instance_id].copy(update={attribute: value})

    def _propagate_id(self, root_type, root_id, attribute, value):
        """
        Recursively updates an attribute (e.g., molecule_id, assembly_id)
        on an instance and all of its hierarchical descendants.

        Args:
            root_type ('str'):
                One of {"assembly", "molecule", "residue", "particle"}.

            root_id ('int'):
                Instance ID of the root object to update.

            attribute ('str'):
                The attribute to update (e.g., "molecule_id", "assembly_id").

            value ('Any'):
                The new value to assign.

        Returns:
            ('list[int]'):
                A flat list of all instance IDs updated (including root).

        Notes:
            - Supported relationships:
            assembly → molecules → residues → particles
            molecule → residues → particles
            residue  → particles
            particle → (nothing)
        """
        updated = []
        # Map each type to its own identity attribute
        self_id_attribute = {"hydrogel": "assembly_id",
                            "molecule": "molecule_id",
                            "peptide": "molecule_id",
                            "protein": "molecule_id",
                            "residue": "residue_id",
                            "particle": "particle_id",}
        assembly_types = ["hydrogel"]
        molecule_types = ["molecule", "peptide", "protein"]
        # 1) Update ROOT (unless attribute corresponds to its own ID)
        if attribute != self_id_attribute.get(root_type):
            self._update_instance(instance_id=root_id,
                                  pmb_type=root_type,
                                  attribute=attribute,
                                  value=value,)
            updated.append((root_type, root_id))
        # 2) Descendants: assembly → molecules
        if root_type in assembly_types:
            for mtype in molecule_types:
                molecule_ids = self._find_instance_ids_by_attribute(pmb_type=mtype,
                                                                    attribute="assembly_id",
                                                                    value=root_id)
                for mid in molecule_ids:
                    updated += self._propagate_id(root_type=mtype,   
                                                  root_id=mid,
                                                  attribute=attribute,
                                                  value=value)
        # 3) Descendants: molecule → residues
        if root_type in molecule_types:
            residue_ids = self._find_instance_ids_by_attribute(pmb_type="residue",
                                                               attribute="molecule_id",
                                                               value=root_id)
            for rid in residue_ids:
                updated += self._propagate_id(root_type="residue",
                                              root_id=rid,
                                              attribute=attribute,
                                              value=value)
        # 4) Descendants: residue → particles
        if root_type == "residue":
            particle_ids = self._find_instance_ids_by_attribute(pmb_type="particle",
                                                                attribute="residue_id",
                                                                value=root_id,)
            for pid in particle_ids:
                self._update_instance(instance_id=pid,
                                      pmb_type="particle",
                                      attribute=attribute,
                                      value=value,)
                updated.append(("particle", pid))
        return updated
    
    def _propose_instance_id(self, pmb_type):
        """
        Propose the next available id for a new TypeInstance.

        Returns:
            ('int'): 
                A non-negative integer that is not already used in the pyMBE database.

        Notes:
            - The method does not fill gaps; it always returns ``max + 1``.
        """
        if pmb_type in self._molecule_like_types:
            used_ids = []
            for t in self._molecule_like_types:
                if t in self._instances:
                    used_ids.extend(self._instances[t].keys())
            if not used_ids:
                return 0
        else:
            if pmb_type not in self._instances or len(self._instances[pmb_type]) == 0:
                return 0
            used_ids = list(self._instances[pmb_type].keys())
        return max(used_ids) + 1
    
    def delete_instance(self, pmb_type, instance_id):
        """
        Delete an instance from the pyMBE database.

        Args:
            pmb_type ('str'):
                Category of the instance (particle, residue, molecule, peptide,
                protein, hydrogel, bond).

            instance_id ('int'):
                Unique identifier of the instance.

        Notes:
            - It applies cascade deletion through the hierarchy:
            assembly → molecules → residues → particles → bonds
            molecule → residues → particles → bonds
            residue  → particles → bonds
            particle → bonds
            bond     → nothing
        """
        # ---- Basic checks ----
        if pmb_type not in self._instances:
            raise ValueError(f"Instance type '{pmb_type}' not found.")
        if instance_id not in self._instances[pmb_type]:
            raise ValueError(f"Instance ID '{instance_id}' not found in '{pmb_type}'.")
        inst = self._instances[pmb_type][instance_id]
        # ===============  CASCADE DELETION  =========================
        # --- Delete children of ASSEMBLY-like objects ---
        if pmb_type in self._assembly_like_types:
            for mtype in self._molecule_like_types:
                mids = self._find_instance_ids_by_attribute(pmb_type=mtype,
                                                            attribute="assembly_id",
                                                            value=instance_id,)
                for mid in mids:
                    self.delete_instance(pmb_type=mtype, 
                                        instance_id=mid)
            # delete particles inside the assembly *even if they have no residue/molecule* (e.g. nodes)
            pids = self._find_instance_ids_by_attribute(pmb_type="particle", 
                                                        attribute="assembly_id",
                                                        value=instance_id)
            for pid in pids:
                self.delete_instance(pmb_type="particle", 
                                    instance_id=pid)
        # --- Delete children of MOLECULE-like objects ---
        if pmb_type in self._molecule_like_types:
            residues = self._find_instance_ids_by_attribute(pmb_type="residue",
                                                            attribute="molecule_id",
                                                            value=instance_id,)
            for rid in residues:
                self.delete_instance(pmb_type="residue", 
                                    instance_id=rid)
        # --- Delete children of RESIDUE ---
        if pmb_type == "residue":
            particles = self._find_instance_ids_by_attribute(pmb_type="particle",
                                                            attribute="residue_id",
                                                            value=instance_id)
            for pid in particles:
                self.delete_instance(pmb_type="particle", 
                                    instance_id=pid)
        # --- Delete children of PARTICLE (only bonds) ---
        if pmb_type == "particle":
            self._delete_bonds_of_particle(instance_id)
        # ===============  FINAL DELETION STEP  ======================
        del self._instances[pmb_type][instance_id]
        if not self._instances[pmb_type]:
            del self._instances[pmb_type]

    def delete_instances(self, pmb_type):
        """
        Remove all instances registered for a given pyMBE type.

        Args:
            pmb_type ('str'):
                Instance category (e.g. ``"particle"``, ``"residue"``,
                ``"molecule"``, ``"protein"``, ``"hydrogel"``).

        Notes:
            - Deletion order is deterministic and safe.
            - If no instances exist for the given type, the method is a no-op.
        """
        if pmb_type not in self._instances:
            return

        # Copy IDs to avoid modifying dict during iteration
        instance_ids = list(self._instances[pmb_type].keys())

        for instance_id in instance_ids:
            self.delete_instance(pmb_type=pmb_type,
                                instance_id=instance_id)

    def delete_reaction(self, reaction_name):
        """
        Delete a reaction template from the pyMBE database.

        Args:
            reaction_name ('str'): 
                label identifying the reaction template in the database.
        """
        if reaction_name not in self._reactions:
            raise ValueError(f"Reaction '{reaction_name}' not found in the pyMBE database.")
        del self._reactions[reaction_name]

    def delete_reactions(self):
        """
        Deletes all reaction templates from the pyMBE database.
        """
        keys = list(self._reactions.keys())
        for key in keys:
            self.delete_reaction(reaction_name=key)

    def delete_template(self, pmb_type, name):
        """
        Delete a template from the pyMBE database.

        Args:
            pmb_type ('str'): 
                The template category.

            name ('str'): 
                The name of the template to delete.
        """
        # Check template exists
        if pmb_type not in self._templates:
            raise ValueError(f"Template type '{pmb_type}' not found.")
        if name not in self._templates[pmb_type]:
            raise ValueError(f"Template '{name}' not found in type '{pmb_type}'.")
        # Check if any instance depends on this template
        if pmb_type in self._instances:
            for inst in self._instances[pmb_type].values():
                if getattr(inst, "name", None) == name:
                    raise ValueError(f"Cannot delete template '{name}' from '{pmb_type}': Instance with ID {getattr(inst, pmb_type + '_id')} depends on it.")
        # Delete
        del self._templates[pmb_type][name]
        # if it is a bond template delete also stored espresso bond instances
        if pmb_type == "bond":
            if name in self.espresso_bond_instances.keys():
                del self.espresso_bond_instances[name]
        # Delete empty groups
        if not self._templates[pmb_type]:
            del self._templates[pmb_type]

    def delete_templates(self, pmb_type):
        """
        Remove all templates registered in the pyMBE database for a given pyMBE type.

        Args:
            pmb_type ('str'):
                Template category (e.g. ``"particle"``, ``"residue"``,
                ``"molecule"``, ``"hydrogel"``).

        Notes:
            - This operation is irreversible.
            - Instance data is not affected.
            - If no templates exist for the given type, the method is a no-op.
        """
        if pmb_type in self._templates:
            templates = list(self._templates[pmb_type].keys())
            for template in templates:
                self.delete_template(pmb_type=pmb_type,
                                     name=template)

    def find_instance_ids_by_name(self, pmb_type, name):
        """
        Return the IDs of all instances of a given pyMBE type that use a
        specific template name.

        Args:
            pmb_type ('str'):
                The instance category to search within.

            name ('str'):
                The template name associated with the instances of interest.

        Returns:
            ('list[int]'):
                A list of instance IDs whose underlying template name matches
                ``name``. The list is empty if no such instances exist.

        Notes:
            - Only exact name matches are considered.
            - This method does not validate whether the corresponding template
            actually exists; it only inspects registered *instances*.
        """
        if pmb_type not in self._instances:
            return []
        result = []
        for iid, inst in self._instances[pmb_type].items():
            if hasattr(inst, "name") and inst.name == name:
                result.append(iid)
        return result

    def get_instance(self, pmb_type, instance_id):
        """
        Retrieve a stored instance by type and instance_id.

        Args:
            pmb_type ('str'): 
                The instance pyMBE category.

            instance_id ('int'): 
                The unique id identifying the given instance

        Returns:
            ('InstanceType'): 
                The stored InstanceTemplate instance corresponding to the provided type and name.

        """
        if instance_id not in self._instances[pmb_type]:
            raise ValueError(f"InstanceTemplate with id = '{instance_id}' not found in type '{pmb_type}'.")
        else:
            return self._instances[pmb_type][instance_id]

    def get_instances(self, pmb_type):
        """
        Return all instances registered for a given pyMBE type.

        Args:
            pmb_type ('str'):
                The pyMBE type (e.g. 'particle', 'residue', 'molecule', 'hydrogel').

        Returns:
            ('dict'):
                Mapping {instance_id: instance_object}.
                Returns an empty dict if no instances exist for the given type.
        """
        return self._instances.get(pmb_type, {}).copy()

    def get_reaction(self,  name):
        """
        Retrieve a reaction stored in the pyMBE database by  name.

        Args:
            name ('str'): 
                The unique id of the reaction to retrieve.

        Returns:
            'Reaction': 
                The stored reaction instance corresponding to the provided name.

        """
        if name not in self._reactions:
            raise ValueError(f"Reaction '{name}' not found in the pyMBE database.")
        else:
            return self._reactions[name]
        
    def get_reactions(self):
        """
        Retrieve all reactions stored in the pyMBE database.

        Returns:
            ('list of Reaction'): 
                List with all stored reaction templates in the pyMBE database.
        """
        return list(self._reactions.values())

    def get_particle_templates_under(self, template_name, pmb_type=None, return_counts=False):
        """
        Returns the names of all particle templates associated with a given
        template by traversing the template hierarchy downward.

        Args:
            template_name ('str'):
                Name of the starting template.

            pmb_type ('str', optional):
                Type of the starting template. If not provided, the type is
                inferred from the database. In this case, the template name
                must be unique across all template types.

            return_counts ('bool', optional):
                If False (default), returns a set of unique particle template
                names. If True, returns a dictionary mapping particle template
                names to the number of times they appear in the hierarchy.

        Returns:
            ('set[str]' or 'dict[str, int]'):
                - If `return_counts=False`: unique particle template names
                - If `return_counts=True`: particle template multiplicities

        Notes:
            - Counting reflects **structural multiplicity**, not instantiated
            particle counts.
            - The returned set contains particle template names only; particle
            states are resolved to their parent particle templates.
        """
        if pmb_type is None:
            pmb_types = self._find_template_types(template_name)
            if len(pmb_types) != 1:
                raise ValueError(f"Template name '{template_name}' is ambiguous: {pmb_types}")
            pmb_type = pmb_types[0]
        counts = self._collect_particle_templates(name=template_name, pmb_type=pmb_type)
        if return_counts:
            return dict(counts)
        return set(counts.keys())


    def get_template(self, pmb_type, name):
        """
        Retrieve a stored template by type and name.

        Args:
            pmb_type ('str'): 
                The template pyMBE category.

            name ('str'): 
                The unique id of the template to retrieve.

        Returns:
            ('TemplateType'): 
                The stored template instance corresponding to the provided type and name.
        """
        if pmb_type not in self._templates:
            raise ValueError(f"There are no {pmb_type} templates defined in the database")

        if name not in self._templates[pmb_type]:
            raise ValueError(f"Template '{name}' not found in type '{pmb_type}'.")
        else:
            return self._templates[pmb_type][name]

    def get_templates(self, pmb_type):
        """
        Return all templates registered for a given pyMBE type.

        Args:
            pmb_type ('str'):
                The pyMBE type (e.g. 'particle', 'residue', 'molecule', 'hydrogel').

        Returns:
            ('dict'):
                Mapping {template_name: template_instance}.
                Returns an empty dict if no templates exist for the given type.
        """
        return self._templates.get(pmb_type, {}).copy()

    def get_es_types_map(self):
        """
        Iterates over all particle templates and extracts the ESPResSo type (`es_type`)
        defined for each state. 

        Returns:
            ('dict[str, int]'):
                A dictionary mapping each particle state to its corresponding ESPResSo type.

        """
        if "particle_state" not in self._templates:
            return {}          
        result = {}
        for _, tpl in self._templates["particle_state"].items():
            result[tpl.name] = tpl.es_type
        return result
    
    def get_particle_id_map(self, object_name):
        """
        Collect all particle IDs associated with an object of given name in the
        pyMBE database. Works for particles, residues, molecules, proteins,
        peptides, and assemblies.

        Args:
            object_name ('str'): 
                Name of the pyMBE object.

        Returns:
            ('dict'): 
                {"all": [particle_ids],
                 "residue_map": {residue_id: [particle_ids]},
                 "molecule_map": {molecule_id: [particle_ids]},
                 "assembly_map": {assembly_id: [particle_ids]},}
        """
        # --- Determine object type by searching in the DB ------------------------
        object_type = None
        object_ids = []
        for pmb_type in self._pmb_types:
            if pmb_type in self._instances:
                for inst_id, inst in self._instances[pmb_type].items():
                    if getattr(inst, "name", None) == object_name:
                        object_type = pmb_type
                        object_ids.append(inst_id)

        if object_type is None:
            raise ValueError(f"No object named '{object_name}' found in database.")
        
        # Maps to return
        id_list = []
        residue_map = {}
        molecule_map = {}
        assembly_map = {}
        # Shortcut access to all particle instances
        particles = self._instances.get("particle", {})
        # Helper: group particle IDs by attribute (molecule_id, residue_id, assembly_id)
        def add_to_map(target_map, key, pid):
            if key is None:
                return
            target_map.setdefault(key, []).append(pid)
        # Case 1: object is a molecule-like type (molecule, protein, peptide)
        if object_type in self._molecule_like_types:
            for mol_id in object_ids:
                molecule_map[mol_id] = []
                for pid, p in particles.items():
                    if p.molecule_id == mol_id:
                        id_list.append(pid)
                        molecule_map[mol_id].append(pid)
                        add_to_map(residue_map, p.residue_id, pid)
                        add_to_map(assembly_map, p.assembly_id, pid)
        # Case 2: object is a residue
        elif object_type == "residue":
            for res_id in object_ids:
                residue_map[res_id] = []
                for pid, p in particles.items():
                    if p.residue_id == res_id:
                        id_list.append(pid)
                        residue_map[res_id].append(pid)
                        add_to_map(molecule_map, p.molecule_id, pid)
                        add_to_map(assembly_map, p.assembly_id, pid)
        # Case 3: object is a particle
        elif object_type == "particle":
            id_list.extend(object_ids)
            for pid in object_ids:
                p = particles[pid]
                add_to_map(residue_map, p.residue_id, pid)
                add_to_map(molecule_map, p.molecule_id, pid)
                add_to_map(assembly_map, p.assembly_id, pid)
        # Case 4: object is an assembly
        elif object_type in self._assembly_like_types:
            for assembly_id in object_ids:
                assembly_map[assembly_id] = []
                for pid, p in particles.items():
                    if p.assembly_id == assembly_id:
                        id_list.append(pid)
                        assembly_map[assembly_id].append(pid)
                        add_to_map(molecule_map, p.molecule_id, pid)
                        add_to_map(residue_map, p.residue_id, pid)
        # Deduplicate + sort IDs
        id_list = sorted(set(id_list))
        return {"all": id_list, "molecule_map": molecule_map,  "residue_map": residue_map, "assembly_map": assembly_map,}
    
    def get_particle_states_templates(self, particle_name):
        """
        Retrieve all particle state templates associated with a given particle.

        Args:
            particle_name ('str'):  
                Name of the particle template.

        Returns:
            'Dict[str, ParticleState]':
                Dictionary mapping state names to `ParticleState` templates.
        """
        states = self._templates.get("particle_state", {})
        particle_states = {state.name: state for state in states.values()
                           if state.particle_name == particle_name}
        return particle_states
