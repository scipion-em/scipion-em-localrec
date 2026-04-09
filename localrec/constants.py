# **************************************************************************
# *
# * Authors:    Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# **************************************************************************
from pwem.constants import (SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL,
                            SYM_TETRAHEDRAL, SYM_I222, SYM_I222r,
                            SYM_In25, SYM_In25r, SCIPION_SYM_NAME,
                            SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r)

# we declarate global constants to multiple usage
LOCALREC_HOME = 'LOCALREC_HOME'

# Supported versions
V2_0 = '2.0'

# Import vector is from Chimera or string
CMM = 0
HAND = 1

# Intepolation method for stitch sub-particles
LINEAR = 0
BSPLINE = 1

# Symmetry dictionary (id: key)
symDict = {0: 'C', 1: 'D', 2: 'T', 3: 'O',
           4: 'I1', 5: 'I2', 6: 'I3', 7: 'I4',
           8: 'I5', 9: 'I6', 10: 'I7', 11: 'I8'}
# inverse symmetry dictionary (key: id)
symDictFromKey = {'C': 0, 'D': 1, 'T': 2, 'O': 3,
                  'I1': 4, 'I2': 5, 'I3': 6, 'I4': 7,
                  'I5': 8, 'I6': 9, 'I7': 10, 'I8': 11}
# Symmetry dictionary to scipion (id: scipion_id)
symDictScipion = {0: SYM_CYCLIC,
                  1: SYM_DIHEDRAL,
                  2: SYM_TETRAHEDRAL,
                  3: SYM_OCTAHEDRAL,
                  4: SYM_I222,
                  5: SYM_I222r,
                  6: SYM_In25,
                  7: SYM_In25r,
                  8: SYM_I2n3,
                  9: SYM_I2n3r,
                  10: SYM_I2n5,
                  11: SYM_I2n5r}
# Symmetry dictionary to scipion from key (key: scipion_id)
symDictScipionFromKey = {'C': SYM_CYCLIC,
                         'D': SYM_DIHEDRAL,
                         'T': SYM_TETRAHEDRAL,
                         'O': SYM_OCTAHEDRAL,
                         'I1': SYM_I222,
                         'I2': SYM_I222r,
                         'I3': SYM_In25,
                         'I4': SYM_In25r,
                         'I5': SYM_I2n3,
                         'I6': SYM_I2n3r,
                         'I7': SYM_I2n5,
                         'I8': SYM_I2n5r}
# valid scipion keys on April-2026 
# SYM_CYCLIC = 0
# SYM_DIHEDRAL = 1
# SYM_DIHEDRAL_X = SYM_DIHEDRAL
# SYM_DIHEDRAL_Y = 2
# SYM_TETRAHEDRAL = 3
# SYM_TETRAHEDRAL_222 = SYM_TETRAHEDRAL
# SYM_TETRAHEDRAL_Z3 = 4
# SYM_TETRAHEDRAL222R  is a degenerated symmetry, that is
# SYM_TETRAHEDRAL222R = 14
# SYM_TETRAHEDRAL_Z3R = 5
# SYM_OCTAHEDRAL = 6  # O
# SYM_I222 = 7
# SYM_I222r = 8
# SYM_In25 = 9  # '2-fold symmetry along y and 5-fold along z
# SYM_In25r = 10  # 'n25' with 180 degree rotation about x
# SYM_I2n3 = 11  # Two-fold symmetry along X and 3-fold along Z
# SYM_I2n3r = 12  # Idem but rotated 180 degree about Y
# SYM_I2n5 = 13  # Two-fold symmetry along Y and 5-fold along Z
# SYM_I2n5r = 14  # Idem but rotated 180 degree about X
# SYM_HELICAL = 15
