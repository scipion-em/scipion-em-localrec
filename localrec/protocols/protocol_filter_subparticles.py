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

import sys
import pprint

from pyworkflow import VERSION_3_0
from pyworkflow.protocol.params import PointerParam, FloatParam, BooleanParam
from pwem.protocols import ProtParticles
from pwem.objects.data import SetOfParticles

from localrec.utils import *
# eventually progressbar will be move to scipion core
from pyworkflow.utils import ProgressBar


class ProtFilterSubParts(ProtParticles):
    """ This protocol mainly filters output particles from two protocols:
    extend symmetry and localized subparticles. It can filter the particles
    (sub-particles) according to spatial distance, view, and angular distance.
    """
    _label = 'filter subparticles'
    _lastUpdateVersion = VERSION_3_0

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
        form.addParam('keepRedundant', BooleanParam, default=False,
                      condition='mindist>0',
                      label='keep overlapping particles in simmetry axis',
                      help="In order to break symmetry constraints, sometimes you want"
                           " all the repetitions of your particle related by symmetry."
                           " but not particles that overlap"
                      )
        form.addParam('distorigin', FloatParam, default=-1,
                      label='Minimum distance to origin (px)',
                      help='In pixels. Minimum distance from subparticle to origin'
                           ' If positive it will drop the subparticles closer'
                           ' to the origin')
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
                  "keepRedundant": self.keepRedundant.get(),
                  "distorigin": self.distorigin.get(),
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
        isSubPart = False;
        firstParticle = inputSet.getFirstItem()
        if firstParticle.hasAttribute('_transorg'):
            isSubPart = True;

        print("Processing particles:")
        progress = ProgressBar(total=len(inputSet), fmt=ProgressBar.NOBAR)
        progress.start()

        sys.stdout.flush()
        step = max(100, len(inputSet) // 100)
        for i, particle in enumerate(inputSet.iterItems(orderBy=['_index'])):

            if i % step == 0:
                progress.update(i+1)

            if (isSubPart):
                partId = int(particle._coordinate._micId)
            else:
                partId = int(particle._index)
            if partId != lastPartId:
                lastPartId = partId
                particlesList = []

            subpart = particle.clone()
            if (isSubPart):
                _, cAngles = geometryFromMatrix(inv(particle._transorg.getMatrix()))
            else:
                _, cAngles = geometryFromMatrix(inv(particle.getTransform().getMatrix()))
            subpart._angles = cAngles

            if not self._filterParticles(params, particlesList, subpart):
                continue

            particlesList.append(subpart)
            outputParts.append(subpart)

        progress.finish()

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

        progress = ProgressBar(len(inputSet), fmt=ProgressBar.NOBAR)
        print("Processing coordinates:")
        progress.start()
        step = max(100, len(inputSet) // 100)
        for i, coord in enumerate(inputSet.iterItems(orderBy=['_subparticle._micId',
                                                     '_micId', 'id'])):
            if i % step == 0:
                progress.update(i+1)

            # The original particle id is stored in the sub-particle as micId
            partId = coord._micId.get()

            # Load the particle if it has changed from the last sub-particle
            if partId != lastPartId:
                self._genOutputCoordinates(subParticles, coordArr, outputSet,
                                           params["mindist"], params["keepRedundant"],
                                           params["distorigin"])
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
        progress.finish()
        self._genOutputCoordinates(subParticles, coordArr, outputSet,
                                   params["mindist"], params["keepRedundant"],
                                   params["distorigin"])
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
        return ['Ilca2015', 'Abrishami2020']

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        return []

    # -------------------------- UTILS functions ------------------------------

    def _genOutputCoordinates(self, subParticles, coordArr,
                              outputSet, minDist, keepRedundant,
                              distorigin):
        for index, coordinate in enumerate(coordArr):
            if minDist > 0 or distorigin > 0:
                subpart = subParticles[index]
                if not filter_mindist(subParticles, subpart, minDist, keepRedundant):
                    continue
                if not filter_distorigin(subParticles, subpart, distorigin):
                    continue
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
