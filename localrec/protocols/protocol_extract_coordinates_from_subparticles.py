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
    Extracts coordinate information from a set of subparticles and maps those
    coordinates onto a parent particle dataset. The protocol is intended for
    workflows in which locally refined or symmetry-expanded subparticles need
    to be represented spatially relative to their original particle context.

    AI Generated:

    Extract Subparticles Coordinates (ProtExtractCoordSubparticles) — User Manual
        Overview

        The Extract Subparticles Coordinates protocol generates a coordinate
        dataset from previously created subparticles. Its main purpose is to
        recover the spatial positions of local regions extracted during
        subparticle analysis so they can be visualized, interpreted, or reused
        in downstream workflows. In cryo-EM studies involving local refinement,
        symmetry expansion, or focused reconstruction, this protocol provides a
        bridge between the refined local regions and the original particle
        organization.

        For biological users, this protocol is particularly useful when studying
        flexible domains, symmetric assemblies, membrane complexes, or particles
        containing regions that benefit from localized processing. Subparticles
        often represent biologically meaningful structural regions extracted
        from larger particles, and recovering their coordinates allows users to
        interpret where these refined regions are located within the parent
        particle population.

        Inputs and Biological Context

        The protocol requires two related datasets. The first is a set of
        subparticles generated from a previous local reconstruction or
        subparticle extraction workflow. The second is the parent particle set
        from which those subparticles originated. The relationship between both
        datasets is biologically important because the extracted coordinates are
        interpreted relative to the original particle framework.

        In practical cryo-EM workflows, subparticles may correspond to flexible
        domains, repeating symmetric units, receptor subcomplexes, or localized
        membrane regions. Recovering their coordinates enables users to inspect
        spatial distributions, evaluate local occupancy, or integrate the
        information into additional refinement or classification procedures.

        Coordinate Representation

        The generated coordinates correspond to the positions of the
        subparticles rather than the centers of the original particles. This
        distinction is important biologically because the resulting dataset
        represents localized structural regions that may occupy different
        spatial positions within the same molecular assembly.

        The output coordinates can be visualized in relation to the parent
        particles, allowing users to understand how local structural regions
        are distributed across the dataset. This is particularly valuable in
        studies of symmetry-expanded particles, conformational heterogeneity,
        or assemblies containing multiple repeating domains.

        Sampling Consistency and Data Integrity

        A critical requirement of the protocol is that the parent particles and
        the subparticles share the same sampling rate. Maintaining consistent
        voxel scaling is essential to preserve biologically meaningful spatial
        relationships. If the sampling rates differ, the extracted coordinates
        may no longer correctly represent the original geometry of the system.

        The protocol is therefore most reliable when subparticles originate
        directly from the provided particle dataset without intermediate
        rescaling or modification. In standard local reconstruction workflows,
        preserving this consistency ensures that the resulting coordinates can
        be interpreted accurately in downstream visualization and analysis.

        Outputs and Their Interpretation

        After execution, the protocol produces a set of coordinates associated
        with the parent particle dataset. These coordinates represent the
        spatial positions of the extracted subparticles and can be used in
        visualization tools, particle inspection workflows, or additional local
        processing pipelines.

        Biologically, the resulting coordinate distribution may reveal patterns
        of structural organization, preferred orientations of flexible domains,
        or spatial occupancy of localized regions. In symmetry-related systems,
        the coordinates can also help users evaluate whether equivalent regions
        behave consistently across the dataset.

        Practical Recommendations

        In routine workflows, users should verify that the input subparticles
        truly originate from the provided parent particle set and that both
        datasets maintain identical sampling conditions. This is especially
        important when combining results from different processing sessions or
        software environments.

        The protocol is most informative when integrated with local refinement
        or focused classification analyses. Visual inspection of the extracted
        coordinates in the context of the parent particles often provides
        valuable biological insight into domain flexibility, symmetry behavior,
        or localized structural variability.

        Final Perspective

        For many cryo-EM studies, subparticle analysis enables structural
        interpretation beyond what can be achieved with whole-particle
        refinement alone. Extracting and recovering the coordinates of these
        localized regions is therefore an important step for connecting refined
        local information back to the biological architecture of the original
        molecular assembly. The protocol provides a practical way to preserve
        and interpret those spatial relationships in downstream structural
        analysis.
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
