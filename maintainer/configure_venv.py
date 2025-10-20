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
#

import os
import sys
import enum
import argparse
import sysconfig
try:
    import espressomd  # pylint: disable=unused-import
    espressomd_found = True
except ModuleNotFoundError:
    espressomd_found = False


parser = argparse.ArgumentParser(description="Configure pyMBE and ESPResSo module paths")
parser.add_argument("--espresso_path", type=str, required=not espressomd_found,
                    help="Path to the ESPResSo build folder")
args = parser.parse_args()


def make_pth(name, path):
    if not os.path.isdir(path):
        raise ValueError(f"Folder '{path}' doesn't exist")
    site_packages = sysconfig.get_path("platlib")
    with open(os.path.join(site_packages, f"{name}.pth"), "w") as f:
        f.write(os.path.realpath(path))


if not hasattr(enum, "StrEnum") and sys.version_info < (3, 11):
    # backport StrEnum
    class StrEnum(str, enum.Enum):
        @staticmethod
        def _generate_next_value_(name, *args, **kwargs):  # pylint: disable=unused-argument
            return name.lower()
    enum.StrEnum = StrEnum


class PyVirtualEnv(enum.StrEnum):
    venv = enum.auto()
    conda = enum.auto()


def detect_py_virtual_environment():
    virtual_env_family_tried = set()
    virtual_env_family_found = set()
    # detect venv
    key = PyVirtualEnv.venv
    virtual_env_family_tried.add(key)
    if os.environ.get("VIRTUAL_ENV"):
        virtual_env_family_found.add(key)
    # detect conda/miniconda
    key = PyVirtualEnv.conda
    virtual_env_family_tried.add(key)
    if int(os.environ.get("CONDA_SHLVL", "0")) >= 1:
        virtual_env_family_found.add(key)
    # process results
    if len(virtual_env_family_found) > 1:
        raise RuntimeError(f"This script should be run in a Python virtual environment, but more than one was detected (found {', '.join(sorted(virtual_env_family_found))})")
    if len(virtual_env_family_found) == 0:
        raise RuntimeError(f"This script should be run in a Python virtual environment, but none was detected (tried {', '.join(sorted(virtual_env_family_tried))})")
    return list(virtual_env_family_found)[0]


virtual_env_family = detect_py_virtual_environment()

if not espressomd_found:
    make_pth("espresso", os.path.join(args.espresso_path, "src", "python"))
make_pth("pymbe", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
