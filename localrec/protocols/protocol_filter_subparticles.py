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

import numpy as np
import math

from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import PointerParam, FloatParam
from pyworkflow.em.protocol import ProtParticles
from pyworkflow.em.data import SetOfParticles

from localrec.utils import *


class ProtFilterSubParts(ProtParticles):
    """ This protocol mainly filters output particles from two protocols:
    extend symmetry and localized subparticles. It can filter the particles
    (sub-particles) according to spatial distance, view, and angular distance.
    """
    _label = 'filter_subunits'
    _lastUpdateVersion = VERSION_1_1

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfCoordinates, SetOfParticles',
                      label='Input type')
        form.addSection('Sub-particles')
        form.addParam('unique', FloatParam, default=-1,
                      label='Angle to keep unique sub-particles (deg)',
                      help='Keep only unique subparticles within angular '
                           'distance. It is useful to remove overlapping '
                           'sub-particles on symmetry axis.')
        form.addParam('mindist', FloatParam, default=-1,
                      label='Minimum distance between sub-particles (px)',
                      help='In pixels. Minimum distance between the '
                           'subparticles in the image. All overlapping ones '
                           'will be discarded.')
        form.addParam('side', FloatParam, default=-1,
                      label='Angle to keep sub-particles from side views (deg)',
                      help='Keep only particles within specified angular '
                           'distance from side views. All others will be '
                           'discarded. ')
        form.addParam('top', FloatParam, default=-1,
                      label='Angle to keep sub-particles from top views (deg)',
                      help='Keep only particles within specified angular '
                           'distance from top views. All others will be '
                           'discarded. ')

        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):

        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):

        inputSet = self.inputSet.get()
        params = {"unique": self.unique.get(),
                  "mindist": self.mindist.get(),
                  "side": self.side.get(),
                  "top": self.top.get()
                  }
        if isinstance(inputSet, SetOfParticles):
            self._particlesOutputStep(inputSet, params)
        else:
            self._coordinateOutputStep(inputSet, params)

    def _particlesOutputStep(self, inputSet, params):

        outputParts = self._createSetOfParticles()
        outputParts.copyInfo(inputSet)

        lastPartId = None
        particlesList = []

        for particle in inputSet.iterItems(orderBy=['_index']):

            partId = int(particle._index)
            if partId != lastPartId:
                lastPartId = partId
                particlesList = []

            subpart = particle.clone()
            _, cAngles = geometryFromMatrix(inv(particle.getTransform().getMatrix()))
            subpart._angles = cAngles

            if not self._filterParticles(params, particlesList, subpart):
                continue

            particlesList.append(subpart)
            outputParts.append(subpart)

        self._defineOutputs(outputParticles=outputParts)
        self._defineSourceRelation(self.inputSet, outputParts)

    def _coordinateOutputStep(self, inputSet, params):

        inputMics = inputSet.getMicrographs()
        outputSet = self._createSetOfCoordinates(inputMics)
        outputSet.copyInfo(inputSet)

        lastPartId = None
        subParticles = []
        coordArr = []
        subParticleId = 0

        for coord in inputSet.iterItems(orderBy=['_subparticle._micId',
                                                 '_micId', 'id']):
            # The original particle id is stored in the sub-particle as micId
            partId = coord._micId.get()

            # Load the particle if it has changed from the last sub-particle
            if partId != lastPartId:
                self._genOutputCoordinates(subParticles, coordArr, outputSet, params["mindist"])
                subParticleId = 0
                coordArr = []
                subParticles = []
                lastPartId = partId

            subParticle = coord._subparticle
            subpart = subParticle.clone()
            _, cAngles = geometryFromMatrix(inv(subParticle._transorg.getMatrix()))
            subpart._angles = cAngles
            subParticleId += 1
            subpart._id = subParticleId
            if not self._filterParticles(params, subParticles, subpart):
                continue
            subParticles.append(subpart)
            coordArr.append(coord.clone())

        self._genOutputCoordinates(subParticles, coordArr, outputSet, params["mindist"])
        self._defineOutputs(outputCoordinates=outputSet)
        self._defineTransformRelation(self.inputSet, self.outputCoordinates)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        inputSet = self.inputSet.get()
        if isinstance(inputSet, SetOfParticles):
            if not inputSet.hasAlignmentProj():
                errors.append('The selected input particles do not have '
                              'alignment information.')
        else:
            firstCoord = inputSet.getFirstItem()
            if not firstCoord.hasAttribute('_subparticle'):
                errors.append('The selected input coordinates are not the '
                              'output from a localized-subparticles protocol.')

        return errors

    def _citations(self):
        return ['Ilca2015']

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        return []

    # -------------------------- UTILS functions ------------------------------
    def _genOutputCoordinates(self, subParticles, coordArr, outputSet, minDist):

        for index, coordinate in enumerate(coordArr):
            if minDist:
                subpart = subParticles[index]
                if filter_mindist(subParticles, subpart, minDist):
                    outputSet.append(coordinate.clone())
            else:
                outputSet.append(coordinate.clone())

    def _filterParticles(self, params, subParticles, subPart):

        if params["side"] > 0:
            # print("Side Filter")
            if not filter_side(subPart, params["side"]):
                return False
        if params["top"] > 0:
            # print("Top Filter", params["top"])
            if not filter_top(subPart, params["top"]):
                return False
        if params["unique"] >= 0:
            # print(np.degrees(subPart._angles))
            if not filter_unique(subParticles, subPart, params["unique"]):
                return False

        return True
