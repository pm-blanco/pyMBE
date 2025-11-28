#
# Copyright (C) 2024 pyMBE-dev team
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

import numpy as np
import espressomd
from pyMBE.lib.handy_functions import get_number_of_particles
import logging
import io

# Create an in-memory log stream
log_stream = io.StringIO()
logging.basicConfig(level=logging.INFO, 
                    format="%(levelname)s: %(message)s",
                    handlers=[logging.StreamHandler(log_stream)] )

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(seed=42)

# Define a set of ions
pmb.define_particle(name="Na", 
                    z=1)
pmb.define_particle(name="Ca", 
                    z=2)
pmb.define_particle(name="Cl", 
                    z=-1)
pmb.define_particle(name="SO4", 
                    z=-2)

type_map=pmb.get_type_map()
# System parameters
c_salt_input = 0.01 * pmb.units.mol/ pmb.units.L
N_SALT_ION_PAIRS = 50
volume = N_SALT_ION_PAIRS/(pmb.N_A*c_salt_input)
L = volume ** (1./3.) # Side of the simulation box

# Create an instance of an espresso system
espresso_system=espressomd.System (box_l = [L.to('reduced_length').magnitude]*3)
espresso_system.setup_type_map(type_list=type_map.values())

#### Unit tests for the added salt

def check_salt_concentration(espresso_system,cation_name,anion_name,c_salt,N_SALT_ION_PAIRS):
    charge_number_map=pmb.get_charge_number_map()
    type_map=pmb.get_type_map()
    espresso_system.setup_type_map(type_list=type_map.values())
    c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,
                                            cation_name=cation_name,
                                            anion_name=anion_name,
                                            c_salt=c_salt)

    np.testing.assert_equal(get_number_of_particles(espresso_system, type_map[cation_name]),N_SALT_ION_PAIRS*abs(charge_number_map[type_map[anion_name]]))
    np.testing.assert_equal(get_number_of_particles(espresso_system, type_map[anion_name]),N_SALT_ION_PAIRS*abs(charge_number_map[type_map[cation_name]]))
    np.testing.assert_almost_equal(c_salt_calculated.m_as("mol/L"), c_salt.m_as("mol/L"))
    cation_ids = pmb.get_particle_id_map(object_name=cation_name)["all"]
    anion_ids  = pmb.get_particle_id_map(object_name=anion_name)["all"]
    for id in cation_ids+anion_ids:
        pmb.delete_particle_in_system(particle_id=id,
                                        espresso_system=espresso_system)


print("*** Unit test: test that create_added_salt works for a 1:1 salt (NaCl-like). Should print the added salt concentration and number of ions ***")
check_salt_concentration(espresso_system=espresso_system,
                        cation_name="Na",
                        anion_name="Cl",
                        c_salt=c_salt_input,
                        N_SALT_ION_PAIRS=N_SALT_ION_PAIRS)
print("*** Unit test passed***")
print("*** Unit test: test that create_added_salt works for a 2:1 salt (CaCl_2-like) ***")
check_salt_concentration(espresso_system=espresso_system,
                        cation_name="Ca",
                        anion_name="Cl",
                        c_salt=c_salt_input,
                        N_SALT_ION_PAIRS=N_SALT_ION_PAIRS)
print("*** Unit test passed***")
print("*** Unit test: test that create_added_salt works for a 1:2 salt (Na_2SO_4-like) ***")
check_salt_concentration(espresso_system=espresso_system,
                        cation_name="Na",
                        anion_name="SO4",
                        c_salt=c_salt_input,
                        N_SALT_ION_PAIRS=N_SALT_ION_PAIRS)
print("*** Unit test passed***")
print("*** Unit test: test that create_added_salt works for a 2:2 salt (CaSO_4-like) ***")
check_salt_concentration(espresso_system=espresso_system,
                        cation_name="Ca",
                        anion_name="SO4",
                        c_salt=c_salt_input,
                        N_SALT_ION_PAIRS=N_SALT_ION_PAIRS)
print("*** Unit test passed***")
print("*** Unit test: check that create_added_salt works for an input c_salt in [particle/lenght**3]. Should print the concentration and number of ions")
c_salt_part=c_salt_input*pmb.N_A
espresso_system.setup_type_map(type_list=type_map.values())
c_salt_calculated = pmb.create_added_salt(espresso_system=espresso_system,
                                            cation_name="Na",
                                            anion_name="Cl",
                                            c_salt=c_salt_part)
np.testing.assert_equal(get_number_of_particles(espresso_system, type_map["Na"]),N_SALT_ION_PAIRS)
np.testing.assert_equal(get_number_of_particles(espresso_system, type_map["Cl"]),N_SALT_ION_PAIRS)
np.testing.assert_almost_equal(c_salt_calculated.m_as("reduced_length**-3"), c_salt_part.m_as("reduced_length**-3"))
cation_ids = pmb.get_particle_id_map(object_name="Na")["all"]
anion_ids = pmb.get_particle_id_map(object_name="Cl")["all"]
for id in cation_ids+anion_ids:
    pmb.delete_particle_in_system(particle_id=id,
                                    espresso_system=espresso_system)

print("*** Unit test: check that create_added_salt raises a ValueError if one provides a cation_name of an object that has been defined with a non-positive charge ***")
input_parameters={"cation_name":"Cl",
                    "anion_name":"SO4",
                    "c_salt":c_salt_input,
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_added_salt, **input_parameters)
print("*** Unit test passed ***")

print("*** Unit test: check that create_added_salt raises a ValueError if one provides a anion_name of an object that has been defined with a non-negative charge ***")
input_parameters={"cation_name":"Na",
                    "anion_name":"Ca",
                    "c_salt":c_salt_input,
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_added_salt, **input_parameters)
print("*** Unit test passed ***")

print("*** Unit test: check that create_added_salt raises a ValueError if one provides a c_salt with the wrong dimensionality ***")
input_parameters={"cation_name":"Na",
                    "anion_name":"Cl",
                    "c_salt":1*pmb.units.nm,
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_added_salt, **input_parameters)
print("*** Unit test passed ***")

# Test that no salt ions are created if the wrong object names are provided
pmb.create_added_salt(espresso_system=espresso_system,
                            cation_name="X",
                            anion_name="Cl",
                            c_salt=c_salt_part)
log_contents = log_stream.getvalue()
assert "Object with name 'X' is not defined in the DataFrame, no ions will be created." in log_contents

# Test that no salt ions are created if the wrong object names are provided
pmb.create_added_salt(espresso_system=espresso_system,
                            cation_name="Na",
                            anion_name="X",
                            c_salt=c_salt_part)
log_contents = log_stream.getvalue()
assert "Object with name 'X' is not defined in the DataFrame, no ions will be created." in log_contents

### Unit tests for the counter ions:


pmb.define_particle(name='0P',
                        z=0)
pmb.define_particle(name='+1P',
                    z=+1)
pmb.define_particle(name='-1P',
                    z=-1)
pmb.define_residue(
    name = 'R1',
    central_bead = '0P',
    side_chains = ['+1P']
    )

pmb.define_residue(
    name = 'R2',
    central_bead = '0P',
    side_chains = ['-1P']
    )

bond_type = 'harmonic'
generic_bond_length=0.4 * pmb.units.nm
generic_harmonic_constant = 400 * pmb.units('reduced_energy / reduced_length**2')

harmonic_bond = {'r_0'    : generic_bond_length,
                'k'      : generic_harmonic_constant,
                }

pmb.define_default_bond(bond_type = bond_type, 
                        bond_parameters = harmonic_bond)
# Add all bonds to espresso system
pmb.add_bonds_to_espresso(espresso_system=espresso_system)
molecule_name = 'positive_polyampholyte'
pmb.define_molecule(name=molecule_name, 
                    residue_list = ['R1']*3+['R2']*2)

molecule_name = 'isoelectric_polyampholyte'
pmb.define_molecule(name=molecule_name, 
                    residue_list = ['R1']*3+['R2']*3)

molecule_name = 'negative_polyampholyte'
pmb.define_molecule(name=molecule_name, 
                    residue_list = ['R1']*2+['R2']*3)

def test_counterions(molecule_name, cation_name, anion_name, espresso_system, expected_numbers):
    pmb.create_molecule(name=molecule_name,
                        number_of_molecules= 2,
                        espresso_system=espresso_system,
                        use_default_bond=True)
    pmb.create_counterions(object_name=molecule_name,
                            cation_name=cation_name,
                            anion_name=anion_name,
                            espresso_system=espresso_system)
    espresso_system.setup_type_map(type_list=type_map.values())
    np.testing.assert_equal(get_number_of_particles(espresso_system, type_map[cation_name]),expected_numbers[cation_name])
    np.testing.assert_equal(get_number_of_particles(espresso_system, type_map[anion_name]),expected_numbers[anion_name])
    molecule_ids = list(pmb.get_particle_id_map(object_name=molecule_name)["molecule_map"].keys())
    for mol_id in molecule_ids:
        pmb.delete_molecule_in_system(molecule_id=mol_id,
                                        espresso_system=espresso_system)
    cation_ids = pmb.get_particle_id_map(object_name=cation_name)["all"]
    anion_ids  = pmb.get_particle_id_map(object_name=anion_name)["all"]
    for id in cation_ids+anion_ids:
        pmb.delete_particle_in_system(particle_id=id,
                                        espresso_system=espresso_system)

print("*** Unit test: check that create_counterions creates the right number of monovalent counter ions for a polyampholyte with positive net charge. Should print the number of ions. ***")

test_counterions(molecule_name='positive_polyampholyte', 
                cation_name="Na", 
                anion_name="Cl", 
                espresso_system=espresso_system, 
                expected_numbers={"Na":4,
                                  "Cl":6})

print("*** Unit test passed ***")

print("*** Unit test: check that create_counterions creates the right number of monovalent counter ions for a polyampholyte at its isoelectric point ***")

test_counterions(molecule_name='isoelectric_polyampholyte', 
                cation_name="Na", 
                anion_name="Cl", 
                espresso_system=espresso_system, 
                expected_numbers={"Na":6,
                                  "Cl":6})

print("*** Unit test passed ***")

print("*** Unit test: check that create_counterions creates the right number of monovalent counter ions for a polyampholyte with a negative net charge ***")

test_counterions(molecule_name='negative_polyampholyte', 
                cation_name="Na", 
                anion_name="Cl", 
                espresso_system=espresso_system, 
                expected_numbers={"Na":6,
                                  "Cl":4})

print("*** Unit test passed ***")

print("*** Unit test: check that create_counterions creates the right number of multivalent counter ions for a polyampholyte ***")

test_counterions(molecule_name='negative_polyampholyte', 
                cation_name="Ca", 
                anion_name="Cl", 
                espresso_system=espresso_system, 
                expected_numbers={"Ca":3,
                                  "Cl":4})

print("*** Unit test passed ***")

print("*** Unit test: check that create_counterions raises a ValueError if the charge number of the cation is not divisible by the negative charge of the polyampholyte ***")
pmb.create_molecule(name='isoelectric_polyampholyte',
                        number_of_molecules= 1,
                        espresso_system=espresso_system,
                        use_default_bond=True)
input_parameters={"cation_name":"Ca",
                    "anion_name":"Cl",
                    "object_name":'isoelectric_polyampholyte',
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_counterions, **input_parameters)
print("*** Unit test passed ***")
print("*** Unit test: check that create_counterions raises a ValueError if the charge number of the anion is not divisible by the positive charge of the polyampholyte ***")
input_parameters={"cation_name":"Na",
                    "anion_name":"SO4",
                    "object_name":'isoelectric_polyampholyte',
                   "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_counterions, **input_parameters)
pmb.delete_molecule_in_system(espresso_system=espresso_system,
                              molecule_id=0)

print("*** Unit test passed ***")
print("*** Unit test: check that no create_counterions does not create counterions for molecules with no charge")
pmb.define_particle(name='0P',
                        z=0)
pmb.define_residue(
    name = 'R0',
    central_bead = '0P',
    side_chains = []
    )

pmb.define_molecule(name='neutral_molecule', 
                    residue_list = ['R0'])

pmb.create_molecule(name='neutral_molecule',
                    number_of_molecules= 1,
                    espresso_system=espresso_system)
pmb.create_counterions(object_name='neutral_molecule',
                            cation_name="Na",
                            anion_name="Cl",
                            espresso_system=espresso_system)


espresso_system.setup_type_map(type_list=type_map.values())
np.testing.assert_equal(get_number_of_particles(espresso_system, type_map["Na"]),0)
np.testing.assert_equal(get_number_of_particles(espresso_system, type_map["Cl"]),0)

# Assert that no counterions are created if the wrong object names are provided
pmb.create_counterions(object_name='test',
                            cation_name="Na",
                            anion_name="Cl",
                            espresso_system=espresso_system)

log_contents = log_stream.getvalue()
assert "Object with name 'test' is not defined in the DataFrame, no counterions will be created." in log_contents


pmb.create_counterions(object_name='isoelectric_polyampholyte',
                            cation_name="Z",
                            anion_name="Cl",
                            espresso_system=espresso_system)

log_contents = log_stream.getvalue()
assert "Object with name 'Z' is not defined in the DataFrame, no counterions will be created." in log_contents

pmb.create_counterions(object_name='isoelectric_polyampholyte',
                            cation_name="Na",
                            anion_name="X",
                            espresso_system=espresso_system)
log_contents = log_stream.getvalue()
assert "Object with name 'X' is not defined in the DataFrame, no counterions will be created." in log_contents

input_parameters={"object_name":'isoelectric_polyampholyte',
                "cation_name":"isoelectric_polyampholyte",
                "anion_name":"Cl",
                "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_counterions, **input_parameters)
input_parameters={"object_name":'isoelectric_polyampholyte',
                "cation_name":"Na",
                "anion_name":'isoelectric_polyampholyte',
                "espresso_system":espresso_system}
np.testing.assert_raises(ValueError, pmb.create_counterions, **input_parameters)

print("*** Unit test passed ***")
