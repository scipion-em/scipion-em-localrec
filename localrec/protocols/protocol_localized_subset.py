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
from pwem.objects import SetOfParticles

class ProtParticleSubset(ProtParticles):
    """
    Creates a subset of particles associated with a corresponding set of
    localized sub-particles, allowing users to recover only the original
    particles that contributed to a localized reconstruction workflow.

    AI Generated:

    Particle Subset by Subparticles (ProtParticleSubset) - User Manual
        Overview

        The Particle Subset protocol selects a reduced collection of
        particles from a larger particle dataset based on their
        relationship to a set of localized sub-particles. Its main
        purpose is to identify which original particles contributed to
        localized reconstruction or sub-particle extraction workflows.

        In cryo-EM studies of highly symmetric or structurally complex
        assemblies, localized reconstruction approaches are frequently
        used to improve the resolution of specific regions such as
        capsomers, flexible domains, membrane regions, or repeating
        structural motifs. During these workflows, only a subset of the
        original particles may successfully contribute to the final
        localized analysis. This protocol allows researchers to recover
        precisely those parent particles.

        Biological Context and Motivation

        In practical cryo-EM workflows, localized sub-particles are
        generated from larger particles using symmetry expansion,
        geometric extraction, or focused refinement strategies. However,
        not all original particles necessarily contribute equally to the
        resulting localized dataset. Some particles may fail quality
        filters, alignment requirements, or geometric constraints during
        extraction and refinement.

        The protocol addresses this problem by identifying the original
        particles associated with successful sub-particle generation.
        This produces a biologically meaningful subset of the initial
        dataset that can be reused in downstream analyses.

        For example, a researcher studying an icosahedral virus may
        perform localized refinement on capsid regions and later wish to
        analyze only the original particles that contributed to the
        refined localized map. Similarly, in studies of flexible
        macromolecular assemblies, the protocol can help isolate the
        subset of particles associated with a particular structural
        region or conformational state.

        Inputs and Dataset Relationships

        The protocol requires two related datasets: the original set of
        particles and a set of localized sub-particles derived from
        them. The relationship between these datasets is essential
        because the protocol uses the provenance of the localized
        sub-particles to determine which original particles should be
        retained.

        Biologically, this relationship preserves the connection between
        local structural refinement and the intact particle population.
        This is particularly important when researchers need to revisit
        the original particle images for additional classification,
        validation, focused refinement, or metadata analysis.

        The protocol assumes that the sub-particles originate from a
        compatible localized extraction workflow. Maintaining proper
        provenance throughout the cryo-EM processing pipeline is
        therefore critical for obtaining reliable results.

        Outputs and Interpretation

        The output is a reduced particle dataset containing only the
        original particles linked to the selected localized
        sub-particles. The resulting subset preserves the original
        particle identities and metadata while restricting the dataset
        to biologically relevant contributors.

        This filtered particle collection can be used in several
        downstream applications. Researchers may employ it for focused
        classification, particle polishing, reconstruction comparison,
        heterogeneity analysis, or additional localized refinement
        workflows.

        From a biological perspective, the resulting subset often
        represents particles that contain better-preserved structural
        features in the region of interest. Consequently, the subset may
        provide improved signal quality for analyses centered on that
        localized region.

        Practical Recommendations

        In routine cryo-EM practice, it is advisable to verify that the
        localized sub-particle dataset originates from a validated
        localized extraction workflow before applying this protocol.
        Broken provenance relationships or incompatible datasets may
        lead to incomplete or biologically inconsistent particle
        selection.

        Researchers should also remember that the resulting subset is
        not necessarily representative of the full structural diversity
        of the original dataset. Because the selection is driven by
        successful localized extraction, the output may be enriched for
        particles with specific orientations, conformations, or local
        structural quality.

        When studying structural heterogeneity, users are encouraged to
        compare the selected subset against the original particle set to
        assess whether any biologically meaningful populations were
        preferentially retained or excluded.

        Final Perspective

        The Particle Subset protocol provides a practical mechanism for
        reconnecting localized cryo-EM analyses with the original
        particle population from which they were derived. By isolating
        the particles that contributed to successful localized
        reconstruction, the protocol supports more targeted biological
        interpretation and enables cleaner downstream processing
        workflows focused on regions of specific structural interest.
    """
    _label = 'particles subset by subparticles'
    _lastUpdateVersion = VERSION_1_1
    OUTPUTPARTICLESNAME = "outputParticles"
    _possibleOutputs = {OUTPUTPARTICLESNAME: SetOfParticles}
    
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

        self._defineOutputs(**{self.OUTPUTPARTICLESNAME: outputSet})
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
