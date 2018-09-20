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
from pyworkflow.em import ImageHandler
from pyworkflow.protocol.params import (PointerParam, BooleanParam, StringParam,
                                        EnumParam, NumericRangeParam,
                                        PathParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol import ProtParticles, ProtParticlePicking

from localrec.utils import *


class ProtFilterSubParts(ProtParticlePicking, ProtParticles):
    """ Extract computed sub-particles from a SetOfParticles. """
    _label = 'filter_subunits'
    _lastUpdateVersion = VERSION_1_1

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label='Input coordinates')
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

        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep',
                                 self.inputCoordinates.get().getObjId())

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self, coordsId):
        """ Create the input file in STAR format as expected by Relion.
        Params:
            particlesId: use this parameters just to force redo of convert if
                the input particles are changed.
        """

        inputCoords = self.inputCoordinates.get()
        inputMics = inputCoords.getMicrographs()
        outputSet = self._createSetOfCoordinates(inputMics)
        outputSet.copyInfo(inputCoords)
        params = {"unique": self.unique.get(),
                  "mindist": self.mindist.get(),
                  "side": self.side.get(),
                  "top": self.top.get()
                  }

        lastPartId = None
        subParticles = []

        for coord in inputCoords.iterItems(orderBy=['_subparticle._micId',
                                                    '_micId', 'id']):
            # The original particle id is stored in the sub-particle as micId
            partId = coord._micId.get()
            print("Mic num", partId)

            # Load the particle if it has changed from the last sub-particle
            if partId != lastPartId:
                subParticles = []
                lastPartId = partId

            subParticle = coord._subparticle
            subpart = subParticle.clone()
            _, cAngles = geometryFromMatrix(inv(subParticle.getTransform().getMatrix()))
            subpart._angles = cAngles
            if params["side"] > 0:
                print("Side Filter")
                if not filter_side(subpart, params["side"]):
                    continue
            if params["top"] > 0:
                print("Top Filter", params["top"])
                if not filter_top(subpart, params["top"]):
                    print(subpart._angles)
                    continue
            if params["unique"] >= 0:
                print("number_of subparticles", len(subParticles))
                if not filter_unique(subParticles, subpart, params["unique"]):
                    continue
            if params["mindist"] > 0:
                print("Mindist Filter", params["mindist"])
                flag, index = filter_mindist(subParticles, subpart, params["mindist"])
                if not flag:
                    print index
                    continue;
            subParticles.append(subpart)
            outputSet.append(coord)


        self._defineOutputs(outputCoordinates=outputSet)
        self._defineTransformRelation(self.inputCoordinates, self.outputCoordinates)


    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        inputCoords = self.inputCoordinates.get()
        firstCoord = inputCoords.getFirstItem()

        if not firstCoord.hasAttribute('_subparticle'):
            errors.append('The selected input coordinates does not are the '
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
    def _getInputParticles(self):
        return self.inputParticles.get()
