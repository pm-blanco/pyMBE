
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
import random
from scipy.spatial import cKDTree

# Distributing points evenly on the surface of a sphere

def uniform_distribution_sites_on_sphere(number_of_edges=2, tolerance=1e-6):
    """
    Distribute points approximately uniformly on the surface of a unit sphere.

    The algorithm uses iterative force-based relaxation, conceptually similar
    to the Thomson problem.

    Args:
        number_of_edges ('int', optional):
            Number of points to distribute on the sphere surface.

        tolerance ('float', optional):
            Convergence tolerance for the iterative relaxation.

    Returns:
        ('list[tuple[float, float, float]]'):
            Point coordinates on the unit sphere centered at ``[0, 0, 0]``.

    Notes:
        – Thomson problem — Wikipedia (overview of the physics problem), Wikipedia.
        – Simple schemes for uniform point distribution. Cheng Guan Koay, J Comput Sci. 2011 Dec;2(4):377–381. doi: 10.1016/j.jocs.2011.06.007.
    """

    # Generates initial configuration

    random.seed(42)
    edges = []
    for i in range(number_of_edges):
        
        theta = random.random() * 2*np.pi     
        phi   = np.arcsin(random.random() * 2 - 1)
        edges.append((np.cos(theta)*np.cos(phi), np.sin(theta)*np.cos(phi), np.sin(phi)))

    # Iterates until fullfil the tolerance

    while 1:
        # Determine the total force acting on each point.
        forces = []
        for i in range(len(edges)):
            p = edges[i]
            f = (0,0,0)
            ftotal = 0
            for j in range(len(edges)):
                if j == i: continue
                q = edges[j]

                # Find the distance vector, and its length.
                dv = (p[0]-q[0], p[1]-q[1], p[2]-q[2])
                dl = np.sqrt(dv[0]**2 + dv[1]**2 + dv[2]**2)

                # The force vector is dv divided by dl^3. (We divide by dl once to make dv a unit vector, then by dl^2 to make its length correspond to the force.)
                dl3 = dl ** 3
                fv = (dv[0]/dl3, dv[1]/dl3, dv[2]/dl3)

                # Add to the total force on the point p.
                f = (f[0]+fv[0], f[1]+fv[1], f[2]+fv[2])

            # Add in the forces array.
            forces.append(f)

            # Add to the running sum of the total forces/distances.
            ftotal = ftotal + np.sqrt(f[0]**2 + f[1]**2 + f[2]**2)

        # Scale the forces to ensure the points do not move too far in one go. Otherwise there will be chaotic jumping around and never any convergence.
        
        if ftotal > 0.25:
            fscale = 0.25 / ftotal
        else:
            fscale = 1

        # Move each point, and normalise. While we do this, also track the distance each point ends up moving.
        
        dist = 0
        for i in range(len(edges)):
            p = edges[i]
            f = forces[i]
            p2 = (p[0] + f[0]*fscale, p[1] + f[1]*fscale, p[2] + f[2]*fscale)
            pl = np.sqrt(p2[0]**2 + p2[1]**2 + p2[2]**2)
            p2 = (p2[0] / pl, p2[1] / pl, p2[2] / pl)
            dv = (p[0]-p2[0], p[1]-p2[1], p[2]-p2[2])
            dl = np.sqrt(dv[0]**2 + dv[1]**2 + dv[2]**2)
            dist = dist + dl
            edges[i] = p2
      
        # Check for convergence and finish.
        
        if dist < tolerance:
            break

    return edges

# Auxiliary functions to calculate the patchy distribution 

def calculate_distance_vector_point(A,p):
    """
    Compute Euclidean distances between a set of 3D points and one 3D point.

    Args:
        A ('iterable'):
            Iterable of 3D points.

        p ('iterable'):
            Reference 3D point.

    Returns:
        ('list[float]'):
            Euclidean distance from each point in ``A`` to ``p``.
    """
    C = []
    for a in A:
        C.append(((a[0] - p[0])**2 + (a[1] - p[1])**2 + (a[2] - p[2])**2)**(1/2))
    return C

def define_patch(points,central_point,patch_size):
    """
    Select the nearest ``patch_size`` points around a central point.

    Args:
        points ('iterable'):
            Iterable of candidate 3D points.

        central_point ('iterable'):
            3D point used as the patch center.

        patch_size ('int'):
            Number of points to include in the patch.

    Returns:
        ('tuple[list[float], list[tuple[float, float, float]]]'):
            Pair with:
            - Distances from all input points to ``central_point``.
            - Coordinates of the selected patch points.
    """
    site_positions            = []
    distance_to_central_point = calculate_distance_vector_point(points,central_point)
    points_index              = sorted(range(len(distance_to_central_point)), key=lambda sub: distance_to_central_point[sub])[:patch_size]
    for index in points_index:
        site_positions.append((points[index][0],points[index][1],points[index][2]))
    return distance_to_central_point, site_positions

def check_patch_overlaps(sites_positions,number_patches):
    """
    Check for overlapping site coordinates between patches.

    Args:
        sites_positions ('list'):
            List containing one list of coordinates per patch.

        number_patches ('int'):
            Number of patches to compare.

    Returns:
        ('int'):
            Returns ``0`` when no overlap is detected.

    Raises:
        ValueError:
            If overlapping coordinates are found between any pair of patches.
    """
    overlapped_sites = []
    for i in range(number_patches-1):
        for j in range(i + 1, number_patches):
            overlapped_sites.append(set(sites_positions[i]) & set(sites_positions[j]))
    for overlap in overlapped_sites:
        if len(overlap) != 0:
            raise ValueError("The patchies are overlapping in {} sites. Please adjust the angle between them.\n".format(len(overlap)));
        else:
            0
    return 0 

def calculate_distance_between_points_on_sphere(points):
    """
    Compute nearest-neighbor distance statistics for points on a sphere.

    Args:
        points ('iterable'):
            Nested iterable with 3D point coordinates.

    Returns:
        ('tuple[float, float, float]'):
            Tuple ``(mean, std, stderr)`` of nearest-neighbor distances.
    """
    points = np.vstack(points)
    tree = cKDTree(points)
    distances, _ = tree.query(points, k=7)  # Assuming 6 neighbors plus the point itself k = 7
    nearest_neighbors_dist = distances[:, 1]      # Removing self (distance = 0)
    avg_dis = np.mean(nearest_neighbors_dist)
    dev_dis = np.std(nearest_neighbors_dist)
    err_dis = dev_dis / np.sqrt(len(nearest_neighbors_dist))
    return avg_dis, dev_dis, err_dis

def calculate_dipole_moment(charges, positions):
    """
    Compute dipole moment for a set of point charges.

    Args:
        charges ('iterable'):
            Charge values.

        positions ('iterable'):
            3D coordinates matching ``charges``.

    Returns:
        ('tuple[numpy.ndarray, float]'):
            Dipole vector and dipole magnitude.
    """
    dipole_moment = np.sum(np.array(charges)[:, None] * np.array(positions), axis=0)
    dipole_magnitude = np.linalg.norm(dipole_moment)
    return dipole_moment, dipole_magnitude

def get_nanoparticle_properties(pmb, name):
    """
    Return the computed properties of a nanoparticle template.

    Args:
        pmb ('pyMBE.pymbe_library'):
            Active pyMBE object with a populated template database.

        name ('str'):
            Name of the nanoparticle template.

    Returns:
        ('dict'): Dictionary with geometric, compositional, and electrostatic
            properties as returned by
            :meth:`~pyMBE.storage.templates.nanoparticle.NanoparticleTemplate.calculate_nanoparticle_properties`.

    Raises:
        ValueError:
            If no nanoparticle template with the given name exists in the database.
    """
    tpl = pmb.db.get_template(name=name, pmb_type="nanoparticle")
    return tpl.calculate_nanoparticle_properties(pmb)


def print_nanoparticle_properties(properties, name=""):
    """
    Print nanoparticle properties in a human-readable format.

    Args:
        properties ('dict'):
            Dictionary of nanoparticle properties as returned by
            :func:`get_nanoparticle_properties`.

        name ('str', optional):
            Nanoparticle name shown in the header. Defaults to empty string.
    """
    header = f"Properties of nanoparticle '{name}'" if name else "Nanoparticle properties"
    print(header)
    for key, value in properties.items():
        print(f"  {key}: {value}")


def calculate_quadrupole_moment(charges, positions):
    """
    Compute quadrupole moment tensor for a set of point charges.

    Args:
        charges ('iterable'):
            Charge values.

        positions ('iterable'):
            3D coordinates matching ``charges``.

    Returns:
        ('tuple[numpy.ndarray, float, numpy.ndarray]'):
            Quadrupole tensor (3x3), Frobenius norm, and eigenvalues.
    """
    Q = np.zeros((3, 3))
    positions = np.array(positions)
    for q, r in zip(charges, positions):
        r_outer = np.outer(r, r)
        Q += q * (3 * r_outer - np.eye(3) * np.dot(r, r))
    quadrupole_magnitude = np.linalg.norm(Q)
    eigenvalues, _ = np.linalg.eigh(Q)
    return Q, quadrupole_magnitude, eigenvalues


def relax_nanoparticle_overlaps(espresso_system, pmb, nanoparticle_name, seed,
                                delta_offset=0.5, verbose=False):
    """
    Gradually grows the LJ offset of the nanoparticle core particle against all
    other non-site particles from zero to its target value, relaxing the system
    after each increment.

    This avoids hard overlaps between the nanoparticle core and salt ions that
    arise at high salt concentrations when the full offset is applied immediately.
    Must be called after :func:`pyMBE.pymbe_library.setup_lj_interactions`, which
    registers the target LJ parameters in ESPResSo.

    Args:
        espresso_system ('espressomd.system.System'):
            ESPResSo system object.

        pmb ('pyMBE.pymbe_library'):
            Active pyMBE object with a populated database.

        nanoparticle_name ('str'):
            Name of the nanoparticle template in the pyMBE database.

        seed ('int'):
            Random seed forwarded to :func:`~pyMBE.lib.handy_functions.relax_espresso_system`.

        delta_offset ('float', optional):
            Offset increment per iteration in reduced length units. Defaults to 0.5.

        verbose ('bool', optional):
            If True, prints progress after each relaxation step. Defaults to False.

    Notes:
        - The Langevin thermostat is left OFF on exit, consistent with
          :func:`~pyMBE.lib.handy_functions.relax_espresso_system`.
        - Re-initialize the thermostat after this call before starting production MD.
    """
    from pyMBE.lib.handy_functions import relax_espresso_system

    # Identify core particle and its ESPResSo type
    np_tpl = pmb.db.get_template(name=nanoparticle_name, pmb_type="nanoparticle")
    core_name = np_tpl.core_particle_name
    type_map = pmb.get_type_map()
    core_type = type_map[core_name]

    # Collect target LJ parameters for core vs every other particle type.
    # Pairs involving particles with sigma=0 (sites) return {} and are skipped.
    particle_templates = pmb.db.get_templates("particle")
    target_lj = {}  # {other_es_type: lj_params_dict}
    for name in particle_templates:
        lj_params = pmb.get_lj_parameters(core_name, name)
        if not lj_params:
            continue
        for state in pmb.db.get_particle_states_templates(particle_name=name).values():
            target_lj[state.es_type] = lj_params

    if not target_lj:
        return

    # Maximum target offset drives the loop (core-core is the largest pair)
    max_target = max(p["offset"].to("reduced_length").magnitude for p in target_lj.values())

    # Reset all core-particle offsets to zero before growing
    for other_type, lj in target_lj.items():
        espresso_system.non_bonded_inter[core_type, other_type].lennard_jones.set_params(
            epsilon=lj["epsilon"].to("reduced_energy").magnitude,
            sigma=lj["sigma"].to("reduced_length").magnitude,
            cutoff=lj["cutoff"].to("reduced_length").magnitude,
            offset=0.0,
            shift="auto")

    current_offset = 0.0
    while current_offset < max_target:
        for other_type, lj in target_lj.items():
            per_pair_target = lj["offset"].to("reduced_length").magnitude
            setup_offset = min(current_offset, per_pair_target)
            espresso_system.non_bonded_inter[core_type, other_type].lennard_jones.set_params(
                epsilon=lj["epsilon"].to("reduced_energy").magnitude,
                sigma=lj["sigma"].to("reduced_length").magnitude,
                cutoff=lj["cutoff"].to("reduced_length").magnitude,
                offset=setup_offset,
                shift="auto")
        relax_espresso_system(espresso_system=espresso_system, seed=seed)
        if verbose:
            print(f"Growing NP core offset: {current_offset:.2f} / {max_target:.2f} (reduced length)")
        current_offset += delta_offset

    # Final pass: restore exact target offsets and do one last relaxation
    for other_type, lj in target_lj.items():
        espresso_system.non_bonded_inter[core_type, other_type].lennard_jones.set_params(
            epsilon=lj["epsilon"].to("reduced_energy").magnitude,
            sigma=lj["sigma"].to("reduced_length").magnitude,
            cutoff=lj["cutoff"].to("reduced_length").magnitude,
            offset=lj["offset"].to("reduced_length").magnitude,
            shift="auto")
    relax_espresso_system(espresso_system=espresso_system, seed=seed)
