# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

from pyworkflow import VERSION_1_2
from pwem.emlib.image import ImageHandler
from pyworkflow.protocol.params import PointerParam
from pwem.protocols import ProtParticles
from pyworkflow.protocol.params import IntParam

# eventually progressbar will be move to scipion core
from pyworkflow.utils import ProgressBar
from pwem.objects import SetOfParticles
from localrec.utils import has_subparticle_provenance


class ProtLocalizedExtraction(ProtParticles):
    """
    Extracts localized subparticles from a larger particle dataset using
    previously computed subparticle coordinates and orientation
    information.

    AI Generated:

    Extract Subparticles (ProtLocalizedExtraction) - User Manual
        Overview

        The Extract Subparticles protocol generates individual localized
        particle images from a larger set of aligned particles. Its main
        purpose is to isolate biologically relevant regions from complex
        macromolecular assemblies so they can be processed independently
        in downstream cryo-EM workflows.

        In many structural biology projects, the full particle contains
        flexible domains, repeated subunits, or localized structural
        features that are difficult to resolve through conventional
        whole-particle refinement. By extracting smaller focused regions,
        researchers can improve local alignment accuracy, enhance
        structural interpretability, and study conformational variability
        with greater precision.

        Inputs and Biological Context

        The protocol requires a set of particles with known projection
        orientations together with a set of coordinates defining the
        locations of the desired subparticles. These coordinates are
        typically generated through localized reconstruction workflows
        that identify specific structural regions inside larger
        assemblies.

        From a biological perspective, this approach is especially useful
        for viruses, ribosomes, membrane complexes, symmetric particles,
        and assemblies containing flexible or partially occupied domains.
        Local extraction allows researchers to focus computational effort
        on the region of biological interest rather than on the entire
        particle.

        Subparticle Box Size

        One of the most important parameters in localized extraction is
        the subparticle box size. This parameter defines the spatial
        region surrounding the selected feature that will be included in
        the extracted image.

        A box size that is too small may truncate important structural
        information or remove neighboring density required for accurate
        alignment and refinement. Conversely, an excessively large box
        may introduce unnecessary noise, increase computational cost, and
        reduce the benefit of focused analysis.

        In practical cryo-EM workflows, the optimal box size usually
        includes the target domain together with a modest surrounding
        margin. The selected dimensions should reflect both the size of
        the biological feature and the expected flexibility of the
        surrounding regions.

        Coordinate Consistency and Data Integrity

        The protocol assumes that the provided coordinates correspond
        correctly to the original particle images and their associated
        orientations. Maintaining this relationship is biologically
        important because the extracted subparticles inherit their
        contextual orientation and structural meaning from the parent
        particles.

        During extraction, regions that would extend beyond the image
        boundaries are excluded to avoid generating incomplete or
        misleading subparticle images. This is particularly relevant when
        localization vectors approach the edges of the particle box or
        when the selected subunits exhibit large positional variability.

        In many biological datasets, a small fraction of subparticles may
        be discarded for this reason. Such exclusions are generally
        expected and usually do not compromise the overall quality of the
        localized reconstruction workflow.

        Preservation of Subparticle Identity

        The extracted outputs preserve the identity and provenance of the
        original localized regions. This continuity is essential for
        downstream analyses such as focused classification, local
        refinement, variability analysis, or symmetry-related studies.

        From a biological standpoint, preserving the relationship between
        parent particles and extracted subparticles allows researchers to
        trace structural variability back to the original assembly and to
        compare equivalent regions across different conformational states
        or experimental conditions.

        Outputs and Their Interpretation

        The protocol produces a new set of particle images containing the
        extracted subparticles. These outputs can be processed using
        standard cryo-EM refinement and classification methods as if they
        were conventional particle datasets.

        Biologically, the extracted dataset represents a focused view of
        specific structural regions within the larger complex. This often
        enables improved local resolution, clearer interpretation of
        flexible domains, and more accurate characterization of local
        conformational changes.

        Practical Recommendations

        In most applications, users should verify the quality of the
        input particle alignments before performing localized extraction.
        Errors in orientation assignment or coordinate placement can
        propagate directly into the extracted dataset and negatively
        affect downstream reconstruction quality.

        Careful selection of the subparticle box size is also critical.
        Pilot extractions followed by visual inspection are often useful
        for determining whether the selected dimensions appropriately
        capture the biological feature of interest.

        When working with highly flexible regions, additional
        classification after extraction is frequently necessary to
        separate distinct conformational states and maximize structural
        interpretability.

        Final Perspective

        Localized extraction is a key step in focused cryo-EM analysis
        workflows because it transforms structural regions embedded
        within large assemblies into independent particle datasets.
        By enabling targeted refinement of specific domains or subunits,
        the protocol supports detailed biological interpretation of
        structural heterogeneity that would otherwise remain obscured in
        conventional whole-particle reconstructions.
    """
    _label = 'extract subparticles'
    _lastUpdateVersion = VERSION_1_2
    OUTPUTPARTICLESNAME = "outputParticles"
    _possibleOutputs = {OUTPUTPARTICLESNAME: SetOfParticles}

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label="Input particles",
                      help='Select the input images from the project.')

        form.addParam('inputCoordinates', PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label='Input coordinates')

        form.addParam('boxSize', IntParam,
                      label='Subparticle box size (px)',
                      help='Select the amount of pixels to extract the '
                           'sub-particles.')

        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):

        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):

        ih = ImageHandler()
        outputStack = self._getPath('particles.mrcs')
        outputImg = ih.createImage()

        inputParticles = self.inputParticles.get()
        inputCoords = self.inputCoordinates.get()
        outputSet = self._createSetOfParticles()
        outputSet.copyInfo(inputParticles)

        boxSize = self.boxSize.get()
        b2 = int(round(boxSize / 2))
        center = np.zeros((boxSize, boxSize))

        ih = ImageHandler()

        i = 0
        outliers = 0
        partIdExcluded = []
        lastPartId = None
        missingProvenanceWarned = False

        progress = ProgressBar(len(inputCoords), fmt=ProgressBar.NOBAR)
        progress.start()
        step = max(100, len(inputCoords) // 100)
        for i, coord in enumerate(inputCoords.iterItems(orderBy=['_subparticle._micId',
                                                    '_micId', 'id'])):
            if i % step == 0:
                progress.update(i+1)

            # The original particle id is stored in the sub-particle as micId
            partId = coord._micId.get()

            # Load the particle if it has changed from the last sub-particle
            if partId != lastPartId:
                particle = inputParticles[partId]

                if particle is None:
                    partIdExcluded.append(partId)
                    self.info("WARNING: Missing particle with id %s from "
                              "input particles set" % partId)
                else:
                    # Now load the particle image to extract later sub-particles
                    img = ih.read(particle)
                    x, y, _, _ = img.getDimensions()
                    data = img.getData()

                lastPartId = partId

            # If particle is not in inputParticles, subparticles will not be
            # generated. Now, subtract from a subset of original particles is
            # supported.
            if partId not in partIdExcluded:
                xpos = coord.getX()
                ypos = coord.getY()

                # Check that the sub-particle will not lay out of the particle
                if (ypos - b2 < 0 or ypos + b2 > y or
                        xpos - b2 < 0 or xpos + b2 > x):
                    outliers += 1
                    continue

                # Crop the sub-particle data from the whole particle image
                center[:, :] = data[ypos - b2:ypos + b2, xpos - b2:xpos + b2]
                outputImg.setData(center)
                i += 1
                outputImg.write((i, outputStack))
                subpart = coord._subparticle
                if (not missingProvenanceWarned and
                        not has_subparticle_provenance(subpart)):
                    self.warning("Subparticle provenance fields "
                                 "(_symmetryGroup/_symmetryOperatorId) are "
                                 "missing for some items. Continuing without "
                                 "failing.")
                    missingProvenanceWarned = True
                subpart.setLocation(
                    (i, outputStack))  # Change path to new stack
                subpart.setObjId(i)  # Ids will be always the same no mater the number of outliers 
                outputSet.append(subpart)

        progress.finish()
        if outliers:
            self.info("WARNING: Discarded %s particles because laid out of the "
                      "particle (for a box size of %d" % (outliers, boxSize))
        outputSet.setIsSubparticles(True)
        self._defineOutputs(**{self.OUTPUTPARTICLESNAME: outputSet})
        self._defineSourceRelation(self.inputParticles, outputSet)

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
        return ['Serban2015', 'Abrishami2020']

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        return []
