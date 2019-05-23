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

# we declarate global constants to multiple usage
LOCALREC_HOME = 'LOCALREC_HOME'

# Supported versions
V1_2_0 = '1.2.0'

# Import vector is from Chimera or string
CMM = 0
HAND = 1

from pyworkflow.em.constants import (
    SYM_CYCLIC, SYM_DIHEDRAL, SYM_TETRAHEDRAL, SYM_OCTAHEDRAL, SYM_I222,
    SYM_I222r, SYM_In25, SYM_In25r)

LOCALREC_SYM_NAME = {}
LOCALREC_SYM_NAME[SYM_CYCLIC] = 'Cn'
LOCALREC_SYM_NAME[SYM_DIHEDRAL] = 'Dn'
LOCALREC_SYM_NAME[SYM_TETRAHEDRAL] = 'T'
LOCALREC_SYM_NAME[SYM_OCTAHEDRAL] = 'O'
LOCALREC_SYM_NAME[SYM_I222] = 'I1'
LOCALREC_SYM_NAME[SYM_I222r] = 'I2'
LOCALREC_SYM_NAME[SYM_In25] = 'I3'
LOCALREC_SYM_NAME[SYM_In25r] = 'I4'