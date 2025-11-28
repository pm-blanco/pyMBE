<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/pyMBE-dev/pyMBE/blob/logos/logo_banner_dark_mode.png">
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/pyMBE-dev/pyMBE/blob/logos/logo_banner.png">
  <img alt="pyMBE logo" src="https://github.com/pyMBE-dev/pyMBE/blob/logos/logo_banner.png">
</picture>

# pyMBE: the Python-based Molecule Builder for ESPResSo

![GitHub Actions](https://github.com/pyMBE-dev/pyMBE/actions/workflows/testsuite.yml/badge.svg)
[![codecov](https://codecov.io/gh/pyMBE-dev/pyMBE/branch/main/graph/badge.svg)](https://codecov.io/gh/pyMBE-dev/pyMBE)
[![docs](https://badgen.net/static/docs/up-to-date/blue/)](pymbe-dev.github.io/pyMBE/pyMBE.html) 
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md) 

[comment]: # (introduction start)

pyMBE provides tools to facilitate building up molecules with complex architectures in the Molecular Dynamics software [ESPResSo](https://espressomd.org/wordpress/). Some examples of molecules that can be set up with pyMBE are polyelectrolytes, hydrogels, peptides and globular proteins. pyMBE bookkeeps all the information about the molecule topology, permitting to link each particle to its corresponding residue and molecule. pyMBE uses the [Pint](https://pint.readthedocs.io/en/stable/) library to enable input parameters in any arbitrary unit system, which is later transformed in the reduced unit system used in ESPResSo.

[comment]: # (introduction stop)

An up-to-date documentation of all methods of the library can be found [here](pymbe-dev.github.io/pyMBE/pyMBE.html) and in the source code.

## Dependencies

- [ESPResSo](https://espressomd.org/wordpress/)
- [Pint](https://pint.readthedocs.io/en/stable/)
- [Pandas](https://pandas.pydata.org/)
- [Pint-Pandas](https://pypi.org/project/Pint-Pandas/)
- [Numpy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [pdoc](https://pdoc.dev/) (for building the docs)
- [CMake](https://cmake.org/) (for running the testsuite)
- any virtual environment manager: [venv](https://docs.python.org/3/library/venv.html),
  [virtualenv](https://virtualenv.pypa.io/en/latest/),
  [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-python),
  [miniconda](https://docs.anaconda.com/free/miniconda),
  [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html),
  etc.

Version requirements are documented in [`requirements.txt`](requirements.txt).

## Contents

- `pyMBE/`: source code of pyMBE
- `pyMBE/parameters/`: folder with various sets of parameters from previous works.
- `figs/`: folder with various images used in the tutorials of pyMBE.
- `maintainer/`: folder with various scripts used by the maintainers.
- `samples/`: folder with various sample scripts showcasing how to use pyMBE to setup different systems.
- `testsuite/`: folder with several test scripts and data for continous integration of the library.
- `tutorials/`: folder with the available tutorials on pyMBE.
- `visualization/`: folder with helper scripts to aid the visualization of vtf trajectories from constant pH and Grand reaction simulations with [VMD](https://www.ks.uiuc.edu/Research/vmd/).
- `AUTHORS.md`: list of authors and contributors of pyMBE.
- `CONTRIBUTING.md`: rules on how to contribute to pyMBE.
- `LICENSE.txt`: license of pyMBE.
- `requirements.txt`: list of required libraries to use pyMBE.

[comment]: # (usage start)

## Usage

### Set up the pyMBE virtual environment

To use pyMBE in your simulations, first clone this repository locally:

```sh
git clone git@github.com:pyMBE-dev/pyMBE.git
```

Please, be aware that pyMBE is intended to be a supporting tool to setup simulations with ESPResSo.
Thus, for most of its functionalities ESPResSo must also be available.
pyMBE supports ESPResSo 4.2 and ESPResSo 4.3-dev.
Following the NEP29 guidelines, we recommend using Python3.10+.
Both NumPy 1 and NumPy 2 are supported.

The pyMBE module needs a Python virtual environment to avoid compatibility issues with its dependencies.
Any virtual environment manager should work, but this readme will use `venv`, which can be installed on Ubuntu as follows:

```sh
sudo apt install python3-venv
```

To set up pyMBE, users need to configure a virtual environment. This is achieved by installing the Python dependencies and setting the path to the ESPResSo build folder, as follows:

```sh
python3 -m venv pymbe  # create a local folder named pymbe containing the environment files
source pymbe/bin/activate  # activate the virtual environment
python3 -m pip install -r requirements.txt "numpy<2.0" "pandas<2.0"
python3 maintainer/configure_venv.py --espresso_path=/home/user/espresso/build # please adapt the espresso path accordingly
python3 simulation_script.py # run a simulation script
deactivate  # deactivate the virtual environment
```

NumPy 2 users should adapt the pip command as follows:

```sh
python3 -m pip install -r requirements.txt "numpy>=2.1" "pandas>=2.0"
```

We highlight that the path `/home/user/espresso/build` is just an example of a possible path to the ESPResSo build folder. 
The user should change this path to match the local absolute path where ESPResSo was built.
Also, ESPResSo must be built with the same NumPy version as the one installed in the environment to avoid API version mismatch.
For more details on how to install ESPResSo, please consult the [ESPResSo installation guide](https://espressomd.github.io/doc4.2.2/installation.html).

The pyMBE virtual environment can be deactivated at any moment as follows:

```sh
deactivate
```

Cluster users who rely on module files to load dependencies should opt for the
following alternative:

```sh
module load ESPResSo/4.2.2-foss-2023a # adapt release if needed
python3 -m venv --system-site-packages pymbe
source pymbe/bin/activate
python3 maintainer/configure_venv.py
python3 -m pip install -r requirements.txt
deactivate
module purge
```

Please note the module files need to be loaded before every activation
of the virtual environment.

Now you can use pyMBE and ESPResSo by activating the virtual environment:

```sh
$ source pymbe/bin/activate
(pymbe) $ python3 -c "import espressomd.version; print(espressomd.version.friendly())"
4.2
(pymbe) $ python3 -c "import pyMBE; print(pyMBE.__file__)"
/home/user/Documents/pyMBE/pyMBE/__init__.py
$ deactivate
```

To use pyMBE in JupyterLab, install extra dependencies and register the virtual environment in a new kernel:

```sh
source pymbe/bin/activate
python3 -m pip install ipykernel "jupyterlab>=4.0.8" "PyOpenGL>=3.1.5" "ipympl>=0.9.3"
python3 -m ipykernel install --user --name=pyMBE
deactivate
```

Please be aware the pyMBE kernel will be registered outside the environment,
typically in your home folder. You can later inspect the list of registered
kernels and delete unwanted ones with the following commands:

```sh
jupyter kernelspec list
jupyter kernelspec uninstall pymbe
```

The JupyterLab main menu will now show a new Python kernel called "pyMBE"
that uses the virtual environment.

### Run simulation scripts

You can run the branched polyampholyte sample with the following commands:

```sh
source pymbe/bin/activate
python3 samples/branched_polyampholyte.py --pH 6
python3 samples/analyze_time_series.py --data_folder samples/time_series/branched_polyampholyte
python3 samples/plot_branched_polyampholyte.py
deactivate
```

### Run tutorials

You can run the interactive tutorials with the following commands:

```sh
source pymbe/bin/activate
jupyter-lab
deactivate
```

In the Jupyter interface, open the `tutorials` folder and then the `pyMBE_tutorial` file.
It will guide you through the creation of polyelectrolytes with pyMBE.
Be sure to use the pyMBE kernel instead of the default Python3 kernel.
The currently active kernel is usually displayed in the top right corner of the notebook.

### Run the testsuite

To make sure your code is valid, please run the testsuite before submitting your contribution:

```sh
source pymbe/bin/activate
make tests -j4
deactivate
```

Here, `-j4` instructs CTest to run the test cases in parallel using 4 CPU cores.
This number can be adjusted depending on your hardware specifications.
You can use `make unit_tests -j4` to run the subset of fast tests, but keep in mind those
won't be able to detect more serious bugs that only manifest themselves in long simulations.
You can also run individual test cases directly, for example with `python3 testsuite/parameter_test.py`.

When contributing new features, consider adding a unit test in the `testsuite/`
folder and a corresponding line in the `testsuite/CTestTestfile.cmake` file.

Every contribution is automatically tested in CI using EESSI (https://www.eessi.io)
and the [EESSI GitHub Action](https://github.com/marketplace/actions/eessi).

[comment]: # (usage stop)

[comment]: # (references start)

## References

Check out the corresponding [paper](https://doi.org/10.1063/5.0216389) to learn more about pyMBE.
If you use pyMBE in your research, please cite our paper:

```bibtex
@article{beyer2024pymbe,
  author = {Beyer, David and Torres, Paola B. and Pineda, Sebastian P. and
            Narambuena, Claudio F. and Grad, Jean-No{\"e}l and Ko{\v{s}}ovan,
            Peter and Blanco, Pablo M.},
  title = {{pyMBE}: The {P}ython-based molecule builder for {ESPResSo}},
  journal = {The Journal of Chemical Physics},
  volume = {161},
  number = {2},
  pages = {022502},
  year = {2024},
  month = jul,
  issn = {0021-9606},
  doi = {10.1063/5.0216389},
}
```

When using a released version of pyMBE, we recommend citing the corresponding
[Zenodo record](https://doi.org/10.5281/zenodo.12102634) in addition to the pyMBE paper,
for example: "We set up our coarse-grained models using pyMBE v0.8.0
\cite{beyer2024pymbe,zenodo2024pymbe}".

Please also make sure to properly cite the original authors if you use the resources provided in the `pyMBE/parameters/` folder.
The relevant references are provided as metadata in the corresponding files.

[comment]: # (references stop)

### Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) to discover the different ways you can
contribute to pyMBE and find our authorship policy and code of conduct.

## License

Copyright (C) 2023-2025 pyMBE-dev team

pyMBE is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

You should have received a [copy](LICENSE.txt) of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
