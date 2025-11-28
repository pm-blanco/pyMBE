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

import argparse
import pathlib
from pyMBE.lib import analysis

parser = argparse.ArgumentParser(description='Sample script analyze time series from the other sample scripts using the binning method.')
parser.add_argument('--data_folder',
                    type=pathlib.Path,
                    required=True,
                    help='path to the data folder with the time series')
args = parser.parse_args()

# Read and analyze time series
analyzed_data=analysis.analyze_time_series(path_to_datafolder=args.data_folder,
                                            ignore_files=["analyzed_data.csv","df.csv"])
analyzed_data.to_csv(args.data_folder / "analyzed_data.csv",
                        index=False)
