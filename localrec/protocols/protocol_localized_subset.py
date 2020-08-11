# **************************************************************************
# *
# * Authors:   Vahid Abrishami (vahid.abrishami@helsinki.fi)
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
# **************************************************************************

from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import PointerParam, FloatParam
from pwem.protocols import ProtParticles

class ProtParticleSubset(ProtParticles):
    """ This protocol make a subset of particles for which there is at least
    one sub-particle.
    """
    _label = 'particles subset by subparticles'
    _lastUpdateVersion = VERSION_1_1

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label='Input particles',
                      help='Set of particles from which the subset will be taken')
        form.addParam('inputSubParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label='Input sub-particles',
                      help='Only the particles in this set of sub-particles will be output')

        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):

        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):

        inputParticles = self.inputParticles.get()
        inputSubParticles = self.inputSubParticles.get()

        outputSet = self._createSetOfParticles()
        outputSet.copyInfo(inputParticles)

        subParticleIds = []

        for subParticle in inputSubParticles:
            coordinate = subParticle.getCoordinate()
            particleId = int(coordinate._micId)
            if not particleId in subParticleIds:
                subParticleIds.append(particleId)

        for particle in inputParticles:
            if particle._objId in subParticleIds:
                outputSet.append(particle)

        self._defineOutputs(outputParticles=outputSet)
        self._defineTransformRelation(inputParticles, outputSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        inputParticles = self.inputSubParticles.get()
        firstSubParticle = inputParticles.getFirstItem()
        if not firstSubParticle.hasAttribute('_transorg'):
            errors.append('The selected input particles does not are the '
                          'output from a localized-extraction protocol.')
        return errors

    def _citations(self):
        return ['Ilca2015', 'Abrishami2020']

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        return []
