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
# from pwem.constants import (SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL,
                            # SYM_TETRAHEDRAL, SYM_I222, SYM_I222r,
                            # SYM_In25, SYM_In25r, SCIPION_SYM_NAME,
                            # SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r)
# from pwem.convert.symmetry import getSymmetryMatrices
from pwem.protocols import ProtParticles, ProtParticlePicking
# import pwem.emlib.metadata as md

# from localrec.utils import load_vectors, create_subparticles
# from localrec.constants import CMM, HAND
# from pyworkflow.utils import ProgressBar
from pwem.objects import SetOfCoordinates


class ProtExtractCoordSubparticles(ProtParticlePicking, ProtParticles):
    """
    Extracts coordinates corresponding to subparticles from an input set of
    subparticles and maps them back onto a parent particle set, where the
    parent particles are treated as micrograph-like references.

    AI Generated:

    Extract Subparticles Coordinates (ProtExtractCoordSubparticles) — User Manual
        Overview

        This protocol generates a set of coordinates representing the spatial
        positions of subparticles relative to a parent particle set. It is
        especially useful in workflows where subparticles have been extracted
        from larger particles and one wants to recover their original spatial
        positions for visualization, analysis, or downstream local
        reconstruction tasks.

        In practical cryo-EM analysis, subparticles are often created by
        cropping localized regions from larger particles. These cropped images
        preserve structural detail but lose the direct visual context of their
        original parent particle. This protocol restores that spatial link by
        converting subparticle metadata into coordinate objects.

        Inputs and General Workflow

        The protocol requires two input datasets.

        The first input is a set of subparticles. These are particles that
        must explicitly be marked as subparticles in their metadata. Each
        subparticle contains coordinate information describing its position
        within the original parent particle.

        The second input is a set of particles. This set acts as the spatial
        reference system and is treated internally in a way similar to a set
        of micrographs. The extracted coordinates will be defined relative to
        these particles.

        During execution, the protocol iterates over all subparticles,
        retrieves the coordinate associated with each one, clones it, and
        appends it to a new output coordinate set.

        Coordinate Extraction Logic

        A practical issue addressed by this protocol is that coordinate
        identifiers in independent subparticle sets may not be globally
        unique. When multiple subparticle groups are merged, duplicate
        coordinate identifiers can appear.

        To avoid this ambiguity, the protocol assigns a new unique identifier
        to every extracted coordinate. This guarantees that each coordinate in
        the output set is uniquely traceable.

        In addition, the micrograph name field is cleared. This is intentional
        and follows the conventions used by other local reconstruction
        protocols, where the parent particle—not an external micrograph—is the
        relevant spatial reference.

        The protocol also stores a clone of the originating subparticle inside
        each coordinate object. This preserves the direct association between
        the coordinate and the original subparticle image.

        Box Size Handling

        The output coordinate set receives a box size derived from the parent
        particle dimensions.

        Specifically, the box size is set to one tenth of the X dimension of
        the parent particles. This is not intended to redefine the physical
        size of the subparticle itself, but rather to provide a small and
        practical coordinate marker size for display and visualization.

        Biological Interpretation

        From a biological perspective, this protocol does not modify particle
        images, orientations, or alignment parameters.

        Its sole purpose is spatial bookkeeping.

        This becomes particularly useful in local reconstruction workflows
        where one wants to inspect where subparticles originate within larger
        macromolecular assemblies. Typical applications include:

        - mapping local domains within large complexes,
        - visualizing symmetry-related subparticle distributions,
        - validating subparticle extraction strategies,
        - relating local refinements back to the parent particle geometry.

        Validation and Consistency Checks

        Before execution, the protocol performs two important checks.

        Sampling Rate Consistency

        The parent particle set and the subparticle set must share the same
        sampling rate. Because coordinates are interpreted geometrically,
        mismatched sampling rates would make the extracted positions
        biologically incorrect.

        The protocol compares both sampling rates using a small tolerance.
        If they differ beyond that threshold, execution is blocked.

        Subparticle Metadata Verification

        The input subparticle set must explicitly carry the metadata flag that
        identifies it as a true subparticle set.

        Older datasets generated by very old versions of local reconstruction
        software may lack this flag. In those cases, the protocol reports the
        issue and recommends re-extracting the subparticles.

        Outputs and Their Interpretation

        The protocol produces a single output:

        outputCoordinates

        This output is a set of coordinates corresponding to the spatial
        locations of all input subparticles.

        Each coordinate:

        - refers to the parent particle reference system,
        - has a unique identifier,
        - preserves a link to the original subparticle.

        These coordinates can be displayed over the parent particles to
        visually inspect where each subparticle was extracted.

        Practical Recommendations

        In biological workflows, this protocol is most useful after
        subparticle extraction and before downstream interpretation.

        It is particularly valuable when:

        - validating whether subparticle extraction targeted the intended
          structural region,
        - inspecting symmetry expansion strategies,
        - checking spatial coverage of local refinement regions.

        For reliable results, always ensure that the parent particles and
        subparticles originate from the same preprocessing pipeline and share
        identical sampling rates.

        Final Perspective

        For most cryo-EM users, this protocol is a coordinate recovery tool
        rather than an image-processing operation.

        Its importance lies in reconnecting local structural information to
        the global geometry of the original particles. This makes it a useful
        interpretative step in subparticle-based analysis, especially when
        studying large assemblies, local flexibility, or region-specific
        reconstructions.
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
            coor = part.getCoordinate().clone()
            # micid is used to relate particle and subparticle coordinates
            # we just clean the mic name since it is standard practice
            # in other localrec protocols.
            coor.setMicName(None)
            coor.setObjId(idx)
            coor._subparticle = part.clone()
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
