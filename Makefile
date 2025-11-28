#!/bin/bash
#
# Copyright (C) 2023-2024 pyMBE-dev team
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

.PHONY: sample
.PHONY: visual 
.PHONY: clean

# whether to run unit tests with code coverage
COVERAGE = 0

# output directory for the code coverage HTML files
COVERAGE_HTML = coverage

# Python executable or launcher, possibly with command line arguments
PYTHON = python3

# number of threads
THREADS = $(shell echo $(MAKEFLAGS) | grep -oP "\\-j *\\d+")

docs:
	mkdir -p ./documentation
	PYTHONWARNINGS=error PDOC_ALLOW_EXEC=0 ${PYTHON} -m pdoc ./pyMBE -o ./documentation --docformat google

unit_tests:
	COVERAGE=$(COVERAGE) ctest --output-on-failure $(THREADS) --test-dir testsuite -LE long --timeout 300

functional_tests:
	COVERAGE=$(COVERAGE) ctest --output-on-failure $(THREADS) --test-dir testsuite -L long

tests:
	COVERAGE=$(COVERAGE) ctest --output-on-failure $(THREADS) --test-dir testsuite

coverage_xml:
	${PYTHON} -m coverage combine testsuite
	${PYTHON} -m coverage report
	${PYTHON} -m coverage xml

coverage_html:
	${PYTHON} -m coverage combine testsuite
	${PYTHON} -m coverage html --directory="${COVERAGE_HTML}"

sample:
	${PYTHON} samples/peptide.py

visual:
	${PYTHON} visualization/vmd-traj.py
	vmd -e visualization.tcl

tutorial:
	jupyter-lab tutorials/pyMBE_tutorial.ipynb

pylint:
	${PYTHON} -m pylint pyMBE/ testsuite/ samples/ maintainer/ visualization/
