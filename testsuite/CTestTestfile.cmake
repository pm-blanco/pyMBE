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
#

cmake_minimum_required(VERSION 3.22.1)
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter NumPy)

file(REAL_PATH "CTestTestfile.cmake" CMAKE_CURRENT_SOURCE_FILE)
cmake_path(GET CMAKE_CURRENT_SOURCE_FILE PARENT_PATH CMAKE_CURRENT_SOURCE_DIR)
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH CMAKE_SOURCE_DIR)

function(pymbe_add_test)
  cmake_parse_arguments(TEST "" "PATH;SUFFIX;THREADS" "LABELS;ARGUMENTS" ${ARGN})
  cmake_path(GET TEST_PATH STEM TEST_NAME)
  if(TEST_SUFFIX)
    set(TEST_NAME "${TEST_NAME}_${TEST_SUFFIX}")
  endif()
  if(DEFINED ENV{COVERAGE} AND "$ENV{COVERAGE}" STREQUAL "1")
    list(APPEND PYTHON_ARGUMENTS "-m" "coverage" "run" "--parallel-mode" "--source=${CMAKE_SOURCE_DIR}")
  endif()
  add_test(${TEST_NAME} "${Python_EXECUTABLE}" ${PYTHON_ARGUMENTS} "${TEST_PATH}" ${TEST_ARGUMENTS})
  set_tests_properties(${TEST_NAME} PROPERTIES SKIP_RETURN_CODE 5)
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS ${TEST_LABELS})
  if(DEFINED TEST_THREADS)
    set_tests_properties(${TEST_NAME} PROPERTIES PROCESSORS ${TEST_THREADS})
  endif()
endfunction()

# functional tests, e.g. long simulations and ensemble averages
pymbe_add_test(PATH hydrogel_tests.py LABELS long THREADS 2 SUFFIX pressure  ARGUMENTS HydrogelTest.test_pressure)
pymbe_add_test(PATH hydrogel_tests.py LABELS long THREADS 2 SUFFIX titration ARGUMENTS HydrogelTest.test_titration)
pymbe_add_test(PATH cph_ideal_tests.py LABELS long THREADS 2)
pymbe_add_test(PATH grxmc_ideal_tests.py LABELS long THREADS 2)
pymbe_add_test(PATH gcmc_tests.py LABELS long)
pymbe_add_test(PATH globular_protein_functional_tests.py LABELS long beyer2024 THREADS 2)
pymbe_add_test(PATH peptide_tests.py LABELS long beyer2024 THREADS 2)
pymbe_add_test(PATH weak_polyelectrolyte_dialysis_test.py LABELS long beyer2024)

# unit tests
pymbe_add_test(PATH test_in_out_pmb_df.py)
pymbe_add_test(PATH serialization_test.py)
pymbe_add_test(PATH test_global_variables.py)
pymbe_add_test(PATH lj_tests.py)
pymbe_add_test(PATH set_particle_acidity_test.py)
pymbe_add_test(PATH bond_tests.py)
pymbe_add_test(PATH generate_perpendicular_vectors_test.py)
pymbe_add_test(PATH define_and_create_molecules_unit_tests.py)
pymbe_add_test(PATH create_molecule_position_test.py)
pymbe_add_test(PATH seed_test.py)
pymbe_add_test(PATH read-write-df_test.py)
pymbe_add_test(PATH parameter_test.py)
pymbe_add_test(PATH henderson_hasselbalch_tests.py)
pymbe_add_test(PATH calculate_net_charge_unit_test.py)
pymbe_add_test(PATH setup_salt_ions_unit_tests.py)
pymbe_add_test(PATH analysis_tests.py)
pymbe_add_test(PATH test_handy_functions.py)
pymbe_add_test(PATH charge_number_map_tests.py)
pymbe_add_test(PATH generate_coordinates_tests.py)
pymbe_add_test(PATH reaction_methods_unit_tests.py)
pymbe_add_test(PATH determine_reservoir_concentrations_unit_test.py)
pymbe_add_test(PATH globular_protein_unit_tests.py)
pymbe_add_test(PATH lattice_builder.py)
pymbe_add_test(PATH hydrogel_builder.py)
