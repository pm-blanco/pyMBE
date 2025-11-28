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
pmb = pyMBE.pymbe_library(seed=42)

def test_arrays_less_equal(arr1, arr2, rtol=1e-7, atol=1e-7):
    """
    Test if all elements in arr1 are less than or equal to the corresponding elements in arr2.

    Parameters:
    arr1 (numpy.ndarray): First array.
    arr2 (numpy.ndarray): Second array.
    rtol (float): Relative tolerance.
    atol (float): Absolute tolerance.
    """
    # Ensure arr1 is less than or equal to arr2 within the given tolerances
    np.testing.assert_array_less(arr1, arr2 + atol + np.abs(arr2) * rtol)
    # For equality, use np.minimum to ensure each element in arr1 is not greater than the corresponding element in arr2
    np.testing.assert_allclose(arr1, np.minimum(arr1, arr2), rtol=rtol, atol=atol)



print("*** generate_random_points_in_a_sphere unit tests ***")
print("*** Unit test: Check that the function creates samples inside a sphere centered around the origin ***")
center=[0,0,0]
n_samples=1000
radius=5.5
samples = pmb.generate_random_points_in_a_sphere(center=center,
                                                 radius=radius,
                                                 n_samples=n_samples)
test_arrays_less_equal(np.linalg.norm(samples, axis=1), radius)
print("*** Unit test passed ***")
print("*** Unit test: Check that the function creates samples inside a sphere centered around a random point ***")
center=[-5.6,1.3444223,np.exp(1)]
n_samples=1000
radius=2.3
samples = pmb.generate_random_points_in_a_sphere(center=center,
                                                 radius=radius,
                                                 n_samples=n_samples)
test_arrays_less_equal(np.linalg.norm(samples-center, axis=1), radius)
print("*** Unit test passed ***")
print("*** Unit test: Check that the function creates samples on the surface of a sphere centered around the origin ***")
center=[0,0,0]
n_samples=1000
radius=4.89
samples = pmb.generate_random_points_in_a_sphere(center = center,
                                                 radius = radius,
                                                 n_samples = n_samples,
                                                 on_surface = True)
np.testing.assert_almost_equal(actual = np.linalg.norm(samples, axis=1),
                               desired = radius,
                               decimal = 5,
                               verbose = True)
print("*** Unit test passed ***")
print("*** Unit test: Check that the function creates samples on the surface of a sphere centered around a random point ***")
center=[1,np.pi,-7.234]
n_samples=1000
radius=3.7
samples = pmb.generate_random_points_in_a_sphere(center = center,
                                                 radius = radius,
                                                 n_samples = n_samples,
                                                 on_surface = True)
np.testing.assert_almost_equal(actual = np.linalg.norm(samples-center, axis=1),
                               desired = radius,
                               decimal = 5,
                               verbose = True)
print("*** Unit test passed ***")
print("*** All unit tests passed ***\n")

print("*** generate_coordinates_outside_sphere unit tests ***")
print("*** Unit test: Check that the function creates samples between two spheres centered around the origin ***")
center=[0,0,0]
n_samples=1000
inner_radius=4.3
outer_radius=5.5
samples = pmb.generate_coordinates_outside_sphere(center=center,
                                                 radius=inner_radius,
                                                 max_dist=outer_radius,
                                                 n_samples=n_samples)
test_arrays_less_equal(np.linalg.norm(samples, axis=1), outer_radius)
test_arrays_less_equal(inner_radius, np.linalg.norm(samples, axis=1))
print("*** Unit test passed ***")
print("*** Unit test: Check that the function creates samples between two spheres centered around a random point ***")
center=[-1,3.5,7.222]
n_samples=1000
inner_radius=3.4
outer_radius=7.2
samples = pmb.generate_coordinates_outside_sphere(center=center,
                                                 radius=inner_radius,
                                                 max_dist=outer_radius,
                                                 n_samples=n_samples)
test_arrays_less_equal(np.linalg.norm(samples-np.asarray(center), axis=1), outer_radius)
test_arrays_less_equal(inner_radius, np.linalg.norm(samples-np.asarray(center), axis=1))
print("*** Unit test passed ***")
print("*** Check that the function raises a ValueError if the user provides a negative radius ***")
np.testing.assert_raises(ValueError, 
                         pmb.generate_coordinates_outside_sphere,
                         center=center,
                         radius=-1.0,
                         max_dist=5.0,
                         n_samples=1)
print("*** Unit test passed ***")
print("*** Check that the function raises a ValueError if the user provides a larger inner radius than outer radius ***")
np.testing.assert_raises(ValueError, 
                         pmb.generate_coordinates_outside_sphere,
                         center=center,
                         radius=10.0,
                         max_dist=5.0,
                         n_samples=1)
print("*** Unit test passed ***")

print("*** All unit tests passed ***\n")
