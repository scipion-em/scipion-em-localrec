# *****************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *              Vahid Abrishami (vahid.abrishami@helsinki.fi)
# *              Roberto Marabini (roberto@cnb.csic.es)
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
from pyworkflow.protocol.params import (PointerParam, BooleanParam, StringParam,
                                        EnumParam, NumericRangeParam, IntParam,
                                        PathParam)
from pwem.constants import (SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL,
                            SYM_TETRAHEDRAL, SYM_I222, SYM_I222r,
                            SYM_In25, SYM_In25r, SCIPION_SYM_NAME,
                            SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r)
from pwem.convert.symmetry import getSymmetryMatrices
from pwem.protocols import ProtParticles, ProtParticlePicking
import pwem.emlib.metadata as md

from localrec.utils import load_vectors, create_subparticles
from localrec.constants import CMM, HAND
from pyworkflow.utils import ProgressBar
from pwem.objects import SetOfCoordinates


class ProtExtractCoordSubparticles(ProtParticlePicking, ProtParticles):
    """ This protocol extract a set of coordinates from a set of subparticles.
    The output set of coordinates can be represented on the particles. Note
    that the output coordinates are the position of the subparticles.
    """
    _label = 'extract subparticles coordinates'
    OUTPUTCOORDINATESNAME = "outputCoordinates"
    _possibleOutputs = {OUTPUTCOORDINATESNAME: SetOfCoordinates}

    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)
        ProtParticles.__init__(self, **args)
        self.allowMpi = False
        self.allowThreads = False

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubParticles', PointerParam,
                      pointerCondition='getIsSubparticles',
                      pointerClass='SetOfParticles',
                      important=True,
                      label="Subparticles",
                      help='Select the input set of subparticles from which the coordinates will be extracted')

        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      important=True,
                      label="Particles (as micrographs)",
                      help='Select the set of particles that will be used as micrographs. The subparticles are'
                           'cropped images from this set')
        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions -----------------------------
    def createOutputStep(self):
        inputSubParticlesSet = self.inputSubParticles.get()
        inputParticlesSet = self.inputParticles.get()
        outputSet = self._createSetOfCoordinates(inputParticlesSet)
        # we set the boxsize the coordinates as 1/0 of the boxsize of the subparticles.
        boxsize = 0.1 * inputParticlesSet.getXDim()
        outputSet.setBoxSize(boxsize)
        idx = 1  # counter to make coordinates id unique

        # the problem that we try to solve here is
        # that subparticle coordinates id  are not unique
        # if we join all the sets of subparticle particle coordinates
        for part in inputSubParticlesSet:
            coor = part.getCoordinate()
            # micid is used to relate particle and subparticle coordinates
            # we just clean the mic name since it is standard practice
            # in other localrec protocols.
            coor.setMicName(None)
            coor.setObjId(idx)
            outputSet.append(coor)
            idx += 1

        self._defineOutputs(**{self.OUTPUTCOORDINATESNAME: outputSet})
        self._defineSourceRelation(inputSubParticlesSet, outputSet)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        subparticles = self.inputSubParticles.get()
        particles = self.inputParticles.get()
        firstSRate = subparticles.getSamplingRate()
        secondSRate = particles.getSamplingRate()
        sRateTol = 1e-4
        if abs(firstSRate - secondSRate) >= sRateTol:
            validateMsgs.append('Particles and subparticles do not have the same sampling rate.'
                                'The sampling rate of both sets of coordinates is expected to be the same within '
                                'tolerance %.4f: %.4f != %.4f' % (sRateTol, firstSRate, secondSRate))
        if not subparticles.getIsSubparticles():
            validateMsgs.append('Set assigned as subparticles is not a subparticles set. '
                                'This error message may appear for old versions of localrec prior to 3.1.0. '
                                'If this is the case please re-extract the subparticles')

        return validateMsgs

    def _summary(self):
        summary = []
        return summary

