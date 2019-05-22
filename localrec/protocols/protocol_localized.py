# *****************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *              Vahid Abrishami (vahid.abrishami@helsinki.fi)
# *              roberto Marabini roberto.marabini@uam.es
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
# *****************************************************************************

from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import (PointerParam, BooleanParam, StringParam,
                                        IntParam,
                                        EnumParam, NumericRangeParam,
                                        PathParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol import ProtParticles, ProtParticlePicking
from localrec.utils import (load_vectors, create_subparticles)
from localrec.constants import CMM, HAND

from pyworkflow.em.constants import (SYM_I222, SYM_I222r, SYM_In25, SYM_In25r,
                                     SYM_CYCLIC, SYM_DIHEDRAL, SYM_TETRAHEDRAL,
                                     SYM_OCTAHEDRAL, SCIPION_SYM_NAME)
from xmipp3.constants import XMIPP_SYM_NAME
from pyworkflow.em.convert.symmetry import getSymmetryMatrices


class ProtLocalizedRecons(ProtParticlePicking, ProtParticles):
    """ This class contains a re-implementation to a method for the
    localized three-dimensional reconstruction of such subunits.
    After determining the particle orientations, local areas
    corresponding to the subunits can be extracted and treated as
    single particles.
    """
    _label = 'localized subparticles'
    _lastUpdateVersion = VERSION_1_1

    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)
        ProtParticles.__init__(self, **args)
        self.allowMpi = False
        self.allowThreads = False

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label="Input particles",
                      help='Select the input images from the project.')

        group = form.addGroup('Symmetry')
        group.addHidden('symmetryGroup', StringParam, default='c1')
        group.addParam('symmetryGroup2', EnumParam,
                      choices=[XMIPP_SYM_NAME[SYM_CYCLIC] +
                               " (" + SCIPION_SYM_NAME[SYM_CYCLIC] + ")",
                               XMIPP_SYM_NAME[SYM_DIHEDRAL] +
                               " (" + SCIPION_SYM_NAME[SYM_DIHEDRAL] + ")",
                               XMIPP_SYM_NAME[SYM_TETRAHEDRAL] +
                               " (" + SCIPION_SYM_NAME[SYM_TETRAHEDRAL] + ")",
                               XMIPP_SYM_NAME[SYM_OCTAHEDRAL] +
                               " (" + SCIPION_SYM_NAME[SYM_OCTAHEDRAL] + ")",
                               XMIPP_SYM_NAME[SYM_I222] +
                               " (" + SCIPION_SYM_NAME[SYM_I222] + ")",
                               XMIPP_SYM_NAME[SYM_I222r] +
                               " (" + SCIPION_SYM_NAME[SYM_I222r] + ")",
                               XMIPP_SYM_NAME[SYM_In25] +
                               " (" + SCIPION_SYM_NAME[SYM_In25] + ")",
                               XMIPP_SYM_NAME[SYM_In25r] +
                               " (" + SCIPION_SYM_NAME[SYM_In25r] + ")"],
                      default=SYM_I222r,
                      label="Symmetry",
                      help="See https://github.com/I2PC/xmipp-portal/wiki/Symmetry"
                           "Symmetry for a description of the symmetry groups "
                           "format in Xmipp.\n"
                           "If no symmetry is present, use _c1_."
                      )
        group.addParam('symmetryOrder', IntParam, default=1,
                      condition='symmetryGroup2<=%d' % SYM_DIHEDRAL,
                      label='Symmetry Order',
                      help='Order of cyclic/dihedral symmetry.')

        group.addParam('randomize', BooleanParam, default=False,
                       label='Randomize the order of the symmetry matrices?',
                       help='Useful for preventing preferred orientations.')
        group.addParam('relaxSym', BooleanParam, default=False,
                       expertLevel=LEVEL_ADVANCED,
                       label='Relax symmetry?',
                       help='Create one random subparticle for each particle ')

        group = form.addGroup('Vectors')
        group.addParam('defineVector', EnumParam, default=CMM,
                       label='Is vector defined by?',
                       choices=['cmm file', 'string'],
                       display=EnumParam.DISPLAY_HLIST)
        group.addParam('vector', NumericRangeParam, default='0,0,1',
                       label='Location vectors', condition="defineVector==1",
                       help='Vector defining the location of the '
                            'subparticles. The vector is defined by 3 '
                            'values x,y,z separated by comma. \n'
                            'More than one vector can be specified separated by'
                            'semicolon. For example: \n'
                            '0,0,1            # Defines only one vector.\n'
                            '0,0,1; 1,0,0;    # Defines two vectors.'
                       )
        group.addParam('vectorFile', PathParam, default='',
                       condition="defineVector==0",
                       label='file obtained by Chimera: ',
                       help='CMM file defining the location(s) of the '
                            'sub-particle(s). Use instead of vector. ')
        group.addParam('length', StringParam, default=-1,
                       label='Alternative length of the vector (A)',
                       help='Use to adjust the sub-particle center. If it '
                            'is <= 0, the length of the given vector. '
                            'Different values must be separated by commas.')

        group = form.addGroup('Sub-particles')
        group.addParam('alignSubparticles', BooleanParam, default=True,
                      label='Align the subparticles?',
                      help='Align sub-particles to the standard orientation. ')


        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -----------------------------
    def createOutputStep(self):
        # return
        inputSet = self._getInputParticles()
        outputSet = self._createSetOfCoordinates(inputSet)
        params = {"symmetryGroup": self.symmetryGroup2.get(),
                  "symmetryOrder": self.symmetryOrder.get(),
                  "vector": self.vector.get(),
                  "vectorFile": self.vectorFile.get(),
                  "length": self.length.get(),
                  "pxSize": self.inputParticles.get().getSamplingRate(),
                  "dim": self.inputParticles.get().getXDim()
                  }

        symMatrices = getSymmetryMatrices(sym = self.symmetryGroup2.get(),
                                          n = self.symmetryOrder.get())
        #symMatrices = getSymMatricesXmipp(self.symmetryGroup.get())

        if self.defineVector == CMM:
            cmmFn = params["vectorFile"]
            vector = " "
        else:
            cmmFn = ""
            vector = params["vector"]

        subpartVectorList = load_vectors(cmmFn, vector,
                                         params["length"],
                                         params["pxSize"])

        for part in inputSet:

            subparticles = create_subparticles(part,
                                               symMatrices,
                                               subpartVectorList,
                                               params["dim"],
                                               self.randomize,
                                               0,
                                               self.alignSubparticles)

            for subpart in subparticles:
                coord = subpart.getCoordinate()
                outputSet.append(coord)
                if part.hasAttribute('_rlnRandomSubset'):
                    coord._subparticle.copyAttributes(part, '_rlnRandomSubset')

        self._defineOutputs(outputCoordinates=outputSet)
        self._defineSourceRelation(self.inputParticles, outputSet)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        pass

    def _citations(self):
        return ['Ilca2015']

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        return []

    # -------------------------- UTILS functions ------------------------------
    def _getInputParticles(self):
        return self.getInputParticlesPointer().get()

    def getInputParticlesPointer(self):
        return self.inputParticles
