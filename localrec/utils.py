# *****************************************************************************
# *
# * Authors:  Serban Ilca
# *           Juha T. Huiskonen (juha@strubi.ox.ac.uk)
# *           J.M. de la Rosa Trevin
# *           Vahid Abrishami (vahid.abrishami@helsinki.fi)
# *
# * Laboratory of Structural Biology,
# * Helsinki Institute of Life Science HiLIFE
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# *****************************************************************************

import math
import random
import numpy as np
from numpy.linalg import inv
from itertools import izip
import xml.etree.ElementTree

from pyworkflow.em.convert.transformations import (vector_norm, unit_vector,
                                           euler_matrix, euler_from_matrix)
from pyworkflow.em.constants import (SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL,
                                     SYM_TETRAHEDRAL, SYM_I222, SYM_I222r,
                                     SYM_In25, SYM_In25r)
from pyworkflow.em.convert.symmetry import getSymmetryMatrices
from pyworkflow.em.data import Coordinate
import pyworkflow.em as em


class Vector3:
    def __init__(self):
        self.vector = np.empty((3,), dtype=float)
        self.length = 0
        self.matrix = np.identity(4, dtype=float)

    def set_vector(self, v):
        self.vector = np.array(v)

    def get_length(self):
        return self.length

    def get_matrix(self):
        return self.matrix[0:3, 0:3]

    def set_length(self, d):
        self.length = float(d)

    def compute_length(self):
        self.set_length(vector_norm(self.vector))

    def compute_unit_vector(self):
        self.set_vector(unit_vector(self.vector))

    def compute_matrix(self):
        """ Compute rotation matrix to align Z axis to this vector. """

        if abs(self.vector[0]) < 0.00001 and abs(self.vector[1]) < 0.00001:
            rot = math.radians(0.00)
            tilt = math.radians(0.00)
        else:
            rot = math.atan2(self.vector[1], self.vector[0])
            tilt = math.acos(self.vector[2])

        psi = 0
        self.matrix = euler_matrix(-rot, -tilt, -psi)

    def print_vector(self):
        print(self.vector[2])
        print("[%.3f,%.3f,%.3f]" % (self.vector[0], self.vector[1], self.vector[2]))


def getSymMatricesXmipp(symmetryGroup):
    """ Return symmetry matrices related to a point group"""

    symmetryGroupLetter = symmetryGroup[0].upper()

    if symmetryGroupLetter == 'I':
        symmetryGroupNum = int(symmetryGroup[1])
        if symmetryGroupNum == 1:
            return getSymmetryMatrices(SYM_I222)
        if symmetryGroupNum == 2:
            return getSymmetryMatrices(SYM_I222r)
        if symmetryGroupNum == 3:
            return getSymmetryMatrices(SYM_In25)
        if symmetryGroupNum == 4:
            return getSymmetryMatrices(SYM_In25r)

    if symmetryGroupLetter == 'C':
        symmetryGroupNum = int(symmetryGroup[1])
        return getSymmetryMatrices(SYM_CYCLIC, n=symmetryGroupNum)

    if symmetryGroupLetter == 'D':
        symmetryGroupNum = int(symmetryGroup[1])
        return getSymmetryMatrices(SYM_DIHEDRAL, n=symmetryGroupNum)

    if symmetryGroupLetter == 'T':
        return getSymmetryMatrices(SYM_TETRAHEDRAL, n=symmetryGroupNum)

    if symmetryGroupLetter == 'O':
        return getSymmetryMatrices(SYM_TETRAHEDRAL, n=SYM_OCTAHEDRAL)


def geometryFromMatrix(matrix):
    from pyworkflow.em.convert.transformations import translation_from_matrix, euler_from_matrix

    shifts = translation_from_matrix(matrix)
    angles = -1.0 * np.ones(3) * euler_from_matrix(matrix, axes='szyz')
    return shifts, angles


def matrixFromGeometry(shifts, angles, inverseTransform):
    """ Create the transformation matrix from a given
    2D shifts in X and Y...and the 3 euler angles.
    """

    from pyworkflow.em.convert.transformations import euler_matrix
    # angles list is in radians, but sign changed
    radAngles = -angles

    M = euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    if inverseTransform:
        from numpy.linalg import inv
        M[:3, 3] = -shifts[:3]
        M = inv(M)
    else:
        M[:3, 3] = shifts[:3]

    return M


def load_vectors(cmm_file, vectors_str, distances_str, angpix):
    """ Load subparticle vectors either from Chimera CMM file or from
    a vectors string. Distances can also be specified for each vector
    in the distances_str. """

    if cmm_file:
        subparticle_vector_list = vectors_from_cmm(cmm_file, angpix)
    else:
        subparticle_vector_list = vectors_from_string(vectors_str)

    if distances_str:
        # Change distances from A to pixel units
        subparticle_distances = [float(x) / angpix for x in
                                 distances_str.split(',')]

        if len(subparticle_distances) != len(subparticle_vector_list):
            raise Exception("Error: The number of distances does not match "
                            "the number of vectors!")

        for vector, distance in izip(subparticle_vector_list,
                                     subparticle_distances):
            if distance > 0:
                vector.set_length(distance)
            else:
                vector.compute_length()
    else:
        for vector in subparticle_vector_list:
            vector.compute_length()

    print("Using vectors:")

    for subparticle_vector in subparticle_vector_list:
        print("Vector: ")
        subparticle_vector.compute_unit_vector()
        subparticle_vector.compute_matrix()
        subparticle_vector.print_vector()
        print("")
        print("Length: %.2f pixels" % subparticle_vector.get_length())
    print("")

    return subparticle_vector_list


def vectors_from_cmm(input_cmm, angpix):
    """function that obtains the input vector from a cmm file"""

    # coordinates in the CMM file need to be in Angstrom
    vector_list = []
    e = xml.etree.ElementTree.parse(input_cmm).getroot()

    for marker in e.findall('marker'):
        x = float(marker.get('x')) / angpix
        y = float(marker.get('y')) / angpix
        z = float(marker.get('z')) / angpix
        id = int(marker.get('id'))
        if id != 1:
            vector = Vector3()
            x = x - x0
            y = y - y0
            z = z - z0
            vector.set_vector([x, y, z])
            vector_list.append(vector)
        else:
            x0 = x
            y0 = y
            z0 = z

    return vector_list


def vectors_from_string(input_str):
    """ Function to parse vectors from an string.
    Our (arbitrary) convention is:
    x1,y1,z1; x2,y2,z2 ... etc
    """
    vectors = []

    for vectorStr in input_str.split(';'):
        v = Vector3()
        v.set_vector([float(x) for x in vectorStr.split(',')])
        vectors.append(v)

    return vectors


def within_mindist(p1, p2, mindist):
    """ Returns True if two particles are closer to each other
    than the given distance in the projection. """

    coordinate_p1 = p1.getCoordinate()
    coordinate_p2 = p2.getCoordinate()
    x1 = coordinate_p1.getX()
    y1 = coordinate_p1.getY()
    x2 = coordinate_p2.getX()
    y2 = coordinate_p2.getY()
    distance_sqr = (x1 - x2) ** 2 + (y1 - y2) ** 2
    mindist_sqr = mindist ** 2

    return distance_sqr < mindist_sqr


def vector_from_two_eulers(rot, tilt):
    """function that obtains a vector from the first two Euler angles"""

    x = math.sin(tilt) * math.cos(rot)
    y = math.sin(tilt) * math.sin(rot)
    z = math.cos(tilt)

    return [x, y, z]


def within_unique(p1, p2, unique):
    """ Returns True if two particles are closer to each other
    than the given angular distance. """

    v1 = vector_from_two_eulers(p1._angles[2], p1._angles[1])
    v2 = vector_from_two_eulers(p2._angles[2], p2._angles[1])
    dp = np.inner(v1, v2) / (vector_norm(v1)) * (vector_norm(v2))

    if dp < -1:
        dp = -1.000

    if dp > 1:
        dp = 1.000

    angle = math.acos(dp)
    # print(np.degrees(angle))

    return angle <= math.radians(unique)


def filter_unique(subparticles, subpart, unique):
    """ Return True if subpart is not close to any other subparticle
        by unique (angular distance).
        For this function we assume that subpart is not contained
        inside."""

    for sp in subparticles:
        if within_unique(sp, subpart, unique):
            return False

    return True


def filter_mindist(subparticles, subpart, mindist):
    """ Return True if subpart is not close to any other subparticle
    by mindist. """

    for sp in subparticles:
        if (sp._id != subpart._id and
                within_mindist(sp, subpart, mindist)):
            return False

    return True


def filter_side(subpart, side):
    return (abs(abs(subpart._angles[1]) - math.radians(90))) < math.radians(side)


def filter_top(subpart, top):
    return (abs(abs(subpart._angles[1]) - math.radians(180))) < math.radians(top)


def filter_subparticles(subparticles, filters):
    return [sp for sp in subparticles
            if all(f(subparticles, sp) for f in filters)]


def create_subparticles(particle, symmetry_matrices, subparticle_vector_list,
                        part_image_size, randomize, subparticles_total,
                        align_subparticles):
    """ Obtain all subparticles from a given particle and set
    the properties of each such subparticle. """

    # Euler angles that take particle to the orientation of the model
    matrix_particle = inv(particle.getTransform().getMatrix())
    shifts, angles = geometryFromMatrix(matrix_particle)

    subparticles = []
    subparticles_total += 1
    symmetry_matrix_ids = range(1, len(symmetry_matrices) + 1)

    if randomize:
        # randomize the order of symmetry matrices, prevents preferred views
        random.shuffle(symmetry_matrix_ids)

    for subparticle_vector in subparticle_vector_list:
        matrix_from_subparticle_vector = subparticle_vector.get_matrix()

        for symmetry_matrix_id in symmetry_matrix_ids:
            # symmetry_matrix_id can be later written out to find out
            # which symmetry matrix created this subparticle
            symmetry_matrix = np.array(symmetry_matrices[symmetry_matrix_id - 1][0:3, 0:3])

            subpart = particle.clone()
            m = np.matmul(matrix_particle[0:3, 0:3], (np.matmul(symmetry_matrix.transpose(),
                                                                matrix_from_subparticle_vector.transpose())))
            angles_org = -1.0 * np.ones(3) * euler_from_matrix(m, 'szyz')
            if align_subparticles:
                  angles = -1.0 * np.ones(3) * euler_from_matrix(m, 'szyz')
            else:
                m2 = np.matmul(matrix_particle[0:3, 0:3], symmetry_matrix.transpose())
                angles = -1.0 * np.ones(3) * euler_from_matrix(m2, 'szyz')

            # subparticle origin
            d = subparticle_vector.get_length()
            x = -m[0, 2] * d + shifts[0]
            y = -m[1, 2] * d + shifts[1]
            z = -m[2, 2] * d

            # save the subparticle coordinates (integer part) relative to the
            # user given image size and as a small shift in the origin (decimal part)
            x_d, x_i = math.modf(x)
            y_d, y_i = math.modf(y)
            alignment = em.Transform()
            alignmentOrg = em.Transform()
            M = matrixFromGeometry(np.array([-x_d, -y_d, 0]), angles, True)
            MOrg = matrixFromGeometry(np.array([-x_d, -y_d, 0]), angles_org, True)
            alignment.setMatrix(M)
            alignmentOrg.setMatrix(MOrg)
            subpart._transorg = alignmentOrg.clone()
            subpart.setTransform(alignment)
            coord = Coordinate()
            coord.setObjId(None)
            coord.setX(int(part_image_size / 2) - x_i)
            coord.setY(int(part_image_size / 2) - y_i)
            coord.setMicId(particle.getObjId())

            if subpart.hasCTF():
                ctf = subpart.getCTF()
                ctf.setDefocusU(subpart.getDefocusU() + z)
                ctf.setDefocusV(subpart.getDefocusV() + z)

            subpart.setCoordinate(coord)
            coord._subparticle = subpart.clone()
            subparticles.append(subpart)

    return subparticles
