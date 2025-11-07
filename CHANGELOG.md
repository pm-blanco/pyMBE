# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- Methods that interact directly with the pyMBE dataframe are now private and stored in a dedicated module in `storage/df_management`. These methods also have been refactored to be stateless methods, i.e. making it impossible for them to change behavior during the pyMBE object lifetime or for the user to change the pyMBE dataframe unless explicitely calling them. This includes the methods: `add_bond_in_df`, `add_value_to_df`, `assign_molecule_id`, `check_if_df_cell_has_a_value`, `check_if_name_is_defined_in_df`, `check_if_multiple_pmb_types_for_name`, `clean_df_row`, `clean_ids_in_df_row`, `copy_df_entry`, `create_variable_with_units`, `convert_columns_to_original_format`, `convert_str_to_bond_object`, `delete_entries_in_df`, `find_bond_key`, `setup_df`. (#142)
- `define_particle_entry_in_df` is now a private method in pyMBE, as it is a convenience method for internal use. (#142)
- The custom `NumpyEncoder` is nos a private class in the private module `storage/df_management` because it is  only internally used in pyMBE for serialization/deserialization. (#142)

## [1.0.0] - 2025-10-08

### Changed
- Parameter sets (pKa values and force field-related parameters) from previous work are now belong to the pyMBE package and are directly accesible from the root path of the package. For example, the path to the Lunkad2021 data set can be accessed as `pmb.root / "parameters" / "peptides" / "Lunkad2021.json"`. (#132)
- Peptide molecules have now their own pyMBE type pmb_type = "peptide" instead that sharing type with custom molecules pmb_type = "molecule". (#126)
- Unified exception handling, soft errors and internal information messages are now handled by loggins and hard errors by raise statements. (#126)
- `lib.handy_functions.setup_langevin_dynamics` now takes `seed` as input argument instead than `SEED` for coherence with the rest of the library. (#118)
- `pd.NA` is now enforced as value for empty cells in `pmb.df`, this prevents transformation of variable types from `int` to `float` which was breaking the code when reading `pmb.df` from file (See #102). (#116)
- The sample script `plot_HH.py` has been replaced for specific examples on how to plot data post-processed with pyMBE: `plot_branched_polyampholyte.py`, `plot_peptide.py`, and `plot_peptide_mixture_grxmc_ideal.py`. (#95)
- Sample scripts now take the pH as an argparse input instead of looping over several pH values. This enables paralization of the sample scripts and avoids conflicts with the current post-processing pipeline. (#95)
- Switched from `os.makedirs` to `Path().mkdir()` to prevent ocasional failure of the scripts when running them in paralel. (#91)
- `pmb.set_reduced_units()` now redefines the reduced units instead of creating a new instance of `pint.UnitRegistry`. Therefore, the user can do operations between objects defining before and after changing the set of reduced units without getting a `ValueError` (#89)
- `pmb.set_reduced_units()` now checks that the arguments provided by the user have the right dimensionality. (#89)
- The constants stored as attributes in `pyMBE.pymbe_library` are now using their values established in the 2019 SI. Their value are taken directly from `scipy.constants` instead of being hard-coded constants. (#86)
- Switched to CTest for testing, allowing to run the tests on paralel (#87)

### Added
- Case-specific methods to delete particles, residues and molecules in pyMBE: `delete_particle_in_system`, `delete_residue_in_system`, `delete_molecule_in_system`. (#137)
- Support for conda and miniconda virtual environments. (#134)
- Private methods for sanity checks, used in various methods to ensure that the inputs are pyMBE objects of the expected type. (#126)
- New benchmark for hydrogels, including scripts to reproduce the data `samples/Landsgesell2022/run_simulations.py` and `samples/Landsgesell2022/plot_pH_vs_alpha.py` and `samples/Landsgesell2022/plot_P_vs_V.py` (#103)
- New sample scripts for hydrogels `samples/build_hydrogel.py` and  `samples/weak_polyacid_hydrogel_grxmc.py` (#103)
- New methods to support building hydrogels with pyMBE `pmb.define_hydrogel`, `pmb.create_hydrogel`, `pmb.initialize_lattice_builder`, `pmb.create_hydrogel_chain`, `pmb.create_hydrogel_node`. (#103)
- CI testing for functions in `lib.handy_functions`. (#118)
- sanity tests for `lib.handy_functions`. (#118)
- Use of `logging`  in `lib.handy_functions` to handle output and error logs. (#118)
- new function to relax the espreso_system `lib.handy_functions.relax_espresso_system`
- Helper method to delete entries from `pmb.df` (#112)
- Code of conduct of our community `CODE_OF_CONDUCT.md`, adhering to the Contributor Covenant v2.1 (#104) 
- New optional argument `backbone_vector` enabling to build molecules along an input vector using `pmb.create_molecule` and `pmb.create_pmb_object` (#99)
- New boolean flag `--ideal` as argparse argument of `samples/globular_protein.py` enabling to run the script without setting up interactions.
- Unit tests for `pmb.create_protein`, `pmb.enable_motion_of_rigid_object`, `pmb.protein_sequence_parser`, `pmb.define_protein`, `pmb.read_protein_vtf_in_df` (#101)
- Library `lattice.py`, a general builder for crystalline lattices. This library is part of on-going project to support hydrogels in pyMBE. (#93)
- New sample script showing how to use the analysis tools in pyMBE for post-processing time series from the sample scripts `analyze_time_series.py` (#95) 
- A new optional argument `ignore_files`  for `lib.analysis.analyze_time_series`, enabling to provide a list of files to be ignored for post-processing of time series. (#95)
- Functional testing for all sample scripts. (#95)
- Unit testing for reaction methods, bonds, object serialization. (#86, #113)

### Fixed
- Occassional crashes when the user provided `backbone_vector` as list instead than a `ndarray` in `pmb.create_molecule()`  (#120)
- Docs in `lib.handy_functions`. (#118)
- Writing and reading `pmb.df` from file does no longer change the variable type from `int` to `float` when there are empty cells in the column. (#116)
- Espresso bond objects stored in `pmb.df` now retain the same value for the  `._bond_id` attribute as the original Espresso objects. (#116)
- Warning handling and coverage in `setup_lj_interaction` (#112)
- Wrong parsing in `pmb.protein_sequence_parser` of input sequences provided as a list of aminoacids using the three letter code. (#101)
- Wrong setup of the rigid object in `pmb.enable_motion_of_rigid_object`, leading to crashes in `samples/globular_protein.py` when enabling the protein motion. (#101)
- The argparse argument `--move_protein` of `samples/globular_protein.py` is now a boolean flag instead of taking arbitrary float values. (#101)
- `lib.analysis.get_dt` now raises a ValueError if the two first two rows of the dataframe have the same values for the time, which break the subsequent code. (#95)
- Removed global state variables, instead they are now created by the constructor of `pyMBE.pymbe_library`. This prevents two instances of the pyMBE library to share the same memory address for their attributes. (#89)
- Required Python dependency versions compatible with ESPResSo 4.2 (#84)
- NumPy 2, Pandas 2 and the development version of ESPResSo are now fully supported. (#106)
- Fixed several deprecated paths and function names in `tutorials/pyMBE_tutorial.ipynb`. (#77, #78, #79, #80, #81)
- Fixed error handling mechanism in the bond search function when the searched bond doesn't exist. (#113)
- Unsupported ESPResSo features now raise a `NotImplementedError` instead of a `ValueError`. (#113)

### Removed
- `pmb.destroy_pmb_object` has been deprecated in favor of case-specific deletion methods. (#137)
- `pmb.create_pmb_object` has been deprecated in favor of case-specific creation methods. (#137)
- `pmb.get_resource()` is no longer needed because pyMBE has been restructured as a package and therefore all of its resources are internally accessible. (#135)
- `lib/create_cg_from_pdb.py` because it was not covered by CI testing and it was imposing a too restrictive coarse-graning of globular proteins to be used for the general public. Future development plans include substituting this functunality for a more general parser. (#135)
- `tutorials/solution_tutorial.ipynb` to avoid code repetition. Instead, the exercise in `pyMBE_tutorial.ipynb` has been adapted to match one of the example in `samples/branched_polyampholyte.py`. (#130)
- `verbose` optional argument has been deprecated in most of pyMBE methods, except those where is needed to silence verbose from ESPResSo. Now the module `loggins` is used instead for handling pyMBE's logs. (#119)
- `lib.handy_functions.minimize_espresso_system_energy` because its name was confusing (it was not only minimizing the system energy) and was changing the parameters of the integrator under the hood. (#118)
- print statements and most of `verbose` in `lib.handy_functions`. (#118)
- `pmb.parse_sequence_from_file` has been removed since it is no longer necesary to parse the sequence from pmb.df (#110)
- `handy_functions.create_random_seed` no longer needed because now instances of pyMBE take the random seed as input (#111)
- `handy_functions.visualize_espresso_system` because it was not used anywhere in the library (#111)
- `handy_functions.do_snapshot_espresso_system` moved to the tutorial because it was a function specific for it (#111)

## [0.8.0] - 2024-06-18

### Added

* Initial release of pyMBE.
