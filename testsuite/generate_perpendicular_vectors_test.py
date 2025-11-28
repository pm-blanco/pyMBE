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
import pyMBE
from itertools import combinations
pmb = pyMBE.pymbe_library(seed=42)


def check_if_different_perpendicular_vectors_are_generated(vector,magnitude,n=50):
    """
    Checks if pmb.generate_trial_perpendicular_vector generates perpendicular vectors to `vector`.

    Args:
        vector(`array`,`float`): vector to which new perpendicular vectors will be generated.
        magnitude(`float`): magnitude of the perpendicular vectors.
        n(`int`): number of perpendicular vectors to be generated.
    """
    perpendicular_vectors = []
    for _ in range(n):
        perpendicular_vector = pmb.generate_trial_perpendicular_vector(vector=vector,
                                                                        magnitude=magnitude)
        perpendicular_vectors.append(perpendicular_vector)
    # Check that the function creates {n} perpendicular vectors to an input vector
    for pvector in perpendicular_vectors:
        np.testing.assert_almost_equal(actual = np.dot(pvector, vector),
                                desired = 0,
                                decimal = 5,
                                verbose = True)
    # Check that the {n} perpendicular vectors are different
    for vector_pair in combinations(perpendicular_vectors, 2):
        assert not np.array_equal(vector_pair[0], vector_pair[1]), \
            f"Error: pmb.generate_trial_perpendicular_vector generated two equal perpendicular vectors v1   = {vector_pair[0]} v2 = {vector_pair[1]}"
        
    # Check that the perpendicular vectors have the same magnitude as the input magnitude
    for pvector in perpendicular_vectors:
        np.testing.assert_almost_equal(actual = np.linalg.norm(pvector), 
                                        desired = magnitude, 
                                        decimal = 5, 
                                        verbose = True)


print("*** generate_trial_perpendicular_vector unit tests ***")
print("*** Unit test: Check that the function raises a ValueError when provided a zero vector  ***")
input_parameters={"vector": [0,0,0],
                   "magnitude":1.0}
np.testing.assert_raises(ValueError, pmb.generate_trial_perpendicular_vector, **input_parameters)
print("*** Unit test passed ***")
print("*** Unit test: Check that the function creates perpendicular vectors to an arbitrary vector of the same magnitude  ***")
vector = pmb.generate_random_points_in_a_sphere(center=[0,0,0],
                                                radius=1, 
                                                n_samples=1, 
                                                on_surface=True)[0]
check_if_different_perpendicular_vectors_are_generated(vector=vector,
                                                    magnitude=1)
print("*** Unit test passed ***")
print("*** Unit test: Check that the function creates perpendicular vectors to a vector with an arbitrary origin ***")
vector = pmb.generate_random_points_in_a_sphere(center=[1,2,3],
                                                radius=1, 
                                                n_samples=1, 
                                                on_surface=True)[0]
check_if_different_perpendicular_vectors_are_generated(vector=vector,
                                                    magnitude=1)
print("*** Unit test passed ***")
print("*** Unit test: Check that the function creates perpendicular vectors with a different magnitude  ***")
vector = pmb.generate_random_points_in_a_sphere(center=[1,2,3],
                                                radius=2, 
                                                n_samples=1, 
                                                on_surface=True)[0]
check_if_different_perpendicular_vectors_are_generated(vector=vector,
                                                        magnitude=3)
print("*** Unit test passed ***")
print("*** All unit tests passed ***")
