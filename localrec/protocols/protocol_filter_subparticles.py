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
from pwem.objects import SetOfCoordinates


class ProtFilterSubParts(ProtParticles):
    """
    Filters sub-particles or sub-particle coordinates according to spatial,
    angular, and geometrical criteria in order to retain biologically meaningful
    localized regions while removing redundant or overlapping information.

    AI Generated:

    Filter Subparticles (ProtFilterSubParts) - User Manual
        Overview

        The Filter Subparticles protocol is designed to refine and curate sets
        of localized sub-particles generated from symmetry expansion or
        localized reconstruction workflows. Its main purpose is to remove
        redundant, overlapping, or biologically irrelevant sub-particles before
        downstream processing steps such as classification, refinement, focused
        reconstruction, or structural interpretation.

        In cryo-EM workflows involving symmetry expansion or localized
        reconstruction, a single particle may generate many related
        sub-particles corresponding to different orientations or symmetry
        operators. While this greatly increases the amount of available data,
        it can also introduce strong redundancy and geometrical overlap. This
        protocol helps reduce those effects by selecting only the most relevant
        sub-particles according to angular orientation, proximity, and spatial
        distribution.

        Inputs and Biological Context

        The protocol accepts either a set of particles or a set of coordinates
        associated with localized sub-particles. These inputs are typically
        produced after workflows that isolate flexible domains, asymmetric
        regions, membrane components, or localized structural features from
        larger complexes.

        From a biological perspective, filtering becomes especially important
        when studying assemblies with high symmetry or repeated structural
        motifs. Without careful filtering, downstream analyses may become
        biased by duplicated information originating from equivalent symmetry
        positions rather than independent observations.

        The protocol preserves the relationship between the filtered
        sub-particles and their original parent particles, allowing users to
        continue processing while maintaining biologically interpretable
        provenance.

        Angular Redundancy and Unique Sub-particles

        One of the most important filtering strategies provided by the protocol
        is the removal of angularly redundant sub-particles. In symmetry-based
        workflows, several sub-particles may correspond to nearly identical
        orientations and therefore contribute little additional information.

        The angular uniqueness criterion allows the user to retain only one
        representative sub-particle within a specified angular distance. This
        is particularly useful near symmetry axes, where multiple equivalent
        orientations may otherwise dominate the dataset.

        Biologically, this filtering helps prevent overrepresentation of
        structurally repeated views and improves the statistical balance of
        downstream analyses. It is especially valuable when studying flexible
        domains or subtle conformational differences that could otherwise be
        masked by symmetry-related redundancy.

        Spatial Distance Filtering

        The protocol can also remove sub-particles that are spatially too close
        to each other within the image plane. This option is intended to reduce
        overlap between neighboring extracted regions that may not represent
        independent signal contributions.

        In practical cryo-EM processing, overlapping sub-particles can produce
        correlations that artificially inflate classification confidence or
        refinement stability. By enforcing a minimum spatial separation, the
        protocol promotes a more independent and biologically meaningful
        dataset.

        Users may optionally preserve overlapping particles associated with
        symmetry-related positions when symmetry-breaking analyses are desired.
        This can be useful when investigating asymmetric features embedded in
        otherwise symmetric assemblies.

        Distance from the Particle Origin

        The protocol also allows filtering according to the distance between a
        sub-particle and the origin of the parent particle. This criterion is
        useful when excluding regions too close to the center of the complex or
        when focusing specifically on peripheral domains.

        In biological applications, this option is often employed to isolate
        membrane-associated regions, distal flexible domains, or external
        structural elements while excluding central density that may not be
        relevant for the intended analysis.

        Orientation-Based Filtering

        Orientation filtering enables selection of sub-particles corresponding
        preferentially to side views or top views. This capability becomes
        particularly important when studying directional flexibility,
        anisotropic sampling, or preferred orientation effects.

        Side-view filtering is often useful for elongated or membrane protein
        complexes where lateral projections provide the most informative
        structural features. Top-view filtering, in contrast, may help isolate
        projections that better reveal rotational symmetry or central channel
        architecture.

        From a biological standpoint, orientation selection can substantially
        improve focused classification and targeted refinement by emphasizing
        the most informative viewing directions.

        Outputs and Their Interpretation

        Depending on the input type, the protocol produces either a filtered
        set of particles or a filtered set of coordinates. The resulting
        dataset preserves the geometrical and biological relationships required
        for subsequent localized reconstruction workflows.

        The filtered outputs are intended to represent a cleaner and less
        redundant subset of the original localized regions. In practice, this
        often leads to improved refinement stability, more interpretable
        classifications, and reduced bias introduced by excessive symmetry
        duplication.

        Practical Recommendations

        In most biological workflows, it is advisable to begin with moderate
        filtering thresholds and visually inspect the resulting particle
        distribution before applying stricter criteria. Excessively aggressive
        filtering may remove meaningful conformational variability, whereas
        insufficient filtering may leave strong redundancy in the dataset.

        Angular uniqueness filtering is typically one of the most effective
        starting points for symmetry-expanded datasets. Minimum distance
        filtering becomes particularly important when extracted regions are
        densely packed or partially overlapping.

        Orientation filtering should be applied carefully and only when there
        is a clear biological or methodological reason to prioritize specific
        projection directions. Over-filtering orientations may unintentionally
        reduce angular coverage and compromise reconstruction quality.

        Final Perspective

        For many cryo-EM studies involving localized reconstruction, symmetry
        expansion, or focused analysis of flexible domains, sub-particle
        filtering is an essential preparation step rather than a simple cleanup
        operation. Careful control of redundancy, overlap, and geometrical
        distribution can strongly influence the interpretability and biological
        reliability of downstream structural analyses.
    """
    _label = 'filter subparticles'
    _lastUpdateVersion = VERSION_3_0
    OUTPUTCOORDINATESNAME = "outputCoordinates"
    _possibleOutputs = {OUTPUTCOORDINATESNAME: SetOfCoordinates}

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
        missingProvenanceWarned = False
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
            if (isSubPart and not missingProvenanceWarned and
                    not has_subparticle_provenance(subpart)):
                self.warning("Subparticle provenance fields "
                             "(_symmetryGroup/_symmetryOperatorId) are "
                             "missing for some items. Continuing without "
                             "failing.")
                missingProvenanceWarned = True
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
        missingProvenanceWarned = False

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
            if (not missingProvenanceWarned and
                    not has_subparticle_provenance(subParticle)):
                self.warning("Subparticle provenance fields "
                             "(_symmetryGroup/_symmetryOperatorId) are "
                             "missing for some items. Continuing without "
                             "failing.")
                missingProvenanceWarned = True
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
        self._defineOutputs(**{self.OUTPUTCOORDINATESNAME: outputSet})
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
