# *****************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
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
from pwem.objects import SetOfCoordinates, Class2D, SetOfParticles
import numpy as np
from scipy.spatial.distance import cdist


class ProtLocalizedSubparticleDistribution(ProtParticlePicking, ProtParticles):
    """ This protocol determines the distibution of number of subparticles.
    Consider a set of proteins and several subparticles in each of them. The number of
    subparticles per particle (protein) is variable. This protocol estimates the histogram
    of the number of subparticles per particle.
    """
    _label = 'subparticles distribution'

    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)
        ProtParticles.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfParticles',
                      important=True,
                      label="Subparticles",
                      help='Select the input set of coordinates of subparticles to be classified.')
        form.addParam('maxNumSubpart', IntParam,
                      important=True,
                      default=60,
                      label="Maximum number of subparticles per particle",
                      help='Maximum number of subparticles per particle.')
        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.groupByCardinalOfClasses)

    # -------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        inputSubParticles = self.inputSet.get()

        # The set of coordinates are read and stored as a dictonary. All subparticles with the same micId
        # are stored in the same dictionary key
        self.subPartsDict = {}

        for part in inputSubParticles.iterItems():
            particleKey = str(part.getMicId())
            if particleKey not in self.subPartsDict:
                self.subPartsDict[particleKey] = []
            self.subPartsDict[particleKey].append(part.clone())

    def groupByCardinalOfClasses(self):

        maxNumber = self.maxNumSubpart.get()
        self.subparticlesclasses = []
        for i in range(0, maxNumber): self.subparticlesclasses.append([])
        # for key in self.coordsDict:
        #    subParts_in_Part = self.coordsDict[key]
        #    self.subparticlesclasses[len(subParts_in_Part)-1].append(subParts_in_Part)
        for key in self.subPartsDict:
            subParts_in_Part = self.subPartsDict[key]
            self.subparticlesclasses[len(subParts_in_Part) - 1].append(subParts_in_Part)

        fnHistogram = self._getExtraPath('histogram.txt')
        f = open(fnHistogram, "a")

        for i in range(0, maxNumber):
            line = ' %i,%i \n' % (i+1, len(self.subparticlesclasses[i]))
            f.writelines(line)
        f.close()

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        pass

    def _citations(self):
        return ['Ilca2015', 'Abrishami2020']

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        return []
