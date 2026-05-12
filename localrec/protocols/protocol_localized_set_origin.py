# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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
# **************************************************************************


import os

import pyworkflow.protocol.params as params
from pyworkflow import VERSION_3_0
from pwem.protocols import EMProtocol
import pwem.objects as emobj
import pyworkflow.utils as pwutils
from os.path import basename

from localrec.utils import vectors_from_string
from localrec.constants import CMM
from localrec.utils import load_vectors
from pwem.objects import Volume


class ProtLocalOrigSampling(EMProtocol):
    """
    Sets the origin and sampling information of a reconstructed subvolume
    so that it remains spatially consistent with the original full map.
    This protocol is especially useful in localized reconstruction workflows
    where extracted regions or subparticles must preserve their correct
    position within the parent structure.

    AI Generated:

    Set Origin To Subvolume (ProtLocalOrigSampling) - User Manual
        Overview

        The Set Origin To Subvolume protocol restores or adjusts the spatial
        reference system of a reconstructed subvolume so that it matches the
        coordinate framework of the original macromolecular assembly. In
        localized reconstruction workflows, subparticles are often reconstructed
        independently from the larger complex. Although these focused maps may
        contain high-resolution information for specific regions, their spatial
        origin may no longer correspond to the original biological structure.

        This protocol addresses that problem by assigning a biologically
        meaningful origin to the reconstructed subvolume using the localization
        vectors defined during subparticle generation. The resulting map can
        then be visualized, interpreted, and compared in the correct structural
        context.

        Biological Context

        In cryo-EM studies of symmetric or modular assemblies, localized
        reconstruction is frequently used to improve the resolution of flexible
        domains, membrane regions, ligand-binding sites, or symmetry-related
        subunits. After reconstruction, these focused maps often need to be
        repositioned relative to the original assembly for interpretation,
        docking, segmentation, or figure preparation.

        By restoring the proper origin, the protocol allows the focused map to
        appear in its biologically correct location when combined with the
        consensus reconstruction or atomic models. This becomes especially
        important when comparing multiple subunit reconstructions or when
        integrating local maps into structural biology visualization software.

        Inputs and Workflow

        The protocol requires an input volume corresponding to the localized
        reconstruction and a localized reconstruction protocol containing the
        vectors that define the subparticle positions. These vectors represent
        the spatial relationship between the extracted subvolume and the
        original particle.

        During execution, the protocol computes the proper origin coordinates
        using the localization information and assigns them to the output map.
        The resulting volume preserves the structural placement of the
        subparticle relative to the original assembly.

        Optionally, the protocol can also assign a new sampling rate to the
        output volume. This is useful when the reconstruction has been processed
        externally or when metadata corrections are needed for downstream
        visualization or interpretation workflows.

        Sampling Considerations

        Correct sampling information is essential for accurate biological
        interpretation. If the voxel size is inconsistent with the true
        experimental sampling, measurements such as distances, fitting quality,
        or structural comparisons may become unreliable.

        In most standard workflows, the original sampling inherited from the
        reconstruction is appropriate. However, when maps have been rescaled,
        imported from external software, or reconstructed using modified
        parameters, manually setting the sampling rate may be necessary.

        Origin Assignment and Structural Interpretation

        The most important role of this protocol is the reassignment of the
        spatial origin. Without a correct origin, focused reconstructions may
        appear displaced relative to the parent structure, making biological
        interpretation difficult or misleading.

        Correct origin placement allows subvolumes to be visualized directly
        together with the consensus map, fitted atomic structures, or symmetry-
        related regions. This facilitates interpretation of conformational
        variability, ligand occupancy, flexible domains, or local structural
        rearrangements.

        File Management and Output

        The protocol can either create a physical copy of the input volume or
        generate a filesystem link to the original map. Linking is efficient
        for routine workflows and reduces storage requirements, while copying
        may be preferable for archival purposes or independent data management.

        The output consists of a volume with updated spatial metadata,
        including the corrected origin and optionally updated sampling
        information. The protocol also allows association of an atomic
        structure for integrated visualization purposes.

        Practical Recommendations

        In most localized reconstruction workflows, this protocol should be
        applied before detailed visualization, figure generation, or atomic
        fitting. Ensuring that subvolumes occupy their correct biological
        position greatly improves interpretability and reduces the risk of
        structural misinterpretation.

        Users should verify that the localization vectors correspond to the
        intended subparticle region and that the sampling information remains
        consistent across all related maps. When multiple focused maps are
        compared together, maintaining a common coordinate framework becomes
        especially important.

        Final Perspective

        Localized reconstruction provides a powerful strategy for resolving
        structural heterogeneity and improving local resolution in cryo-EM.
        However, the biological interpretation of focused maps depends strongly
        on preserving their spatial relationship to the original assembly. By
        restoring accurate origin and sampling information, this protocol helps
        ensure that localized reconstructions remain meaningful within the
        broader structural context.
    """
    _label = 'Set origin to subvolume'
    _program = ""
    _lastUpdateVersion = VERSION_3_0
    OUTPUTVOLUMENAME = "outputVolume"
    _posibleOutputs = {OUTPUTVOLUMENAME: Volume}

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inVolume', params.PointerParam, label="Input Volume", important=True,
                      pointerClass='Volume',
                      help='input volume')
        form.addParam('inputProtocol', params.PointerParam,
                      label="Input protocol", important=True,
                      pointerClass='ProtLocalizedRecons',
                      help="protocol that contains the vector defining "
                           "the location of the subparticles")
        form.addParam('copyFiles', params.BooleanParam,
                      label="copy (true)/link(false) volume:",
                      help="Option YES:\nA new volume file will be copied"
                      "otherwise a link to the input volume is maded\n"
                      "default = false", expertLevel=params.LEVEL_ADVANCED,
                      default=False)
        form.addParam('setSampling', params.BooleanParam,
                      label="Set SamplingRate",
                      help="Option YES:\nA new volume object will be created with "
                           "the given SamplinRate."
                           "This SamplinRate will NOT be set in the map file header.\n\n",
                      default=False)
        form.addParam('samplingRate', params.FloatParam,
                      condition='setSampling',
                      label=pwutils.Message.LABEL_SAMP_RATE)
        form.addParam('atomStruct', params.PointerParam,
                      pointerClass="AtomStruct", allowsNull=True,
                      important=True,
                      label='Atomic structure',
                      help="Thia Atomic Structure will be shown when visualizing"
                           " the results")
    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('assignStep')

    # --------------------------- STEPS functions -----------------------------

    def assignStep(self):
        # Create a Volume object
        self.inVol = self.inVolume.get()
        self.outVol = emobj.Volume()
        self.outVol.copy(self.inVol)

        # Set sampling Rate (or copy from input)
        if self.setSampling.get():
            samplingRate = self.samplingRate
        else:
            samplingRate = self.inVol.getSamplingRate()

        self.outVol.setSamplingRate(samplingRate)

        # set Origin
        # get location of the subparticles (px)
        # this is a string with 3 coor
        inputProtocol = self.inputProtocol.get()
        if inputProtocol.defineVector == CMM:
            cmmFn = inputProtocol.vectorFile.get()
            vector = " "
        else:
            cmmFn = ""
            vector = inputProtocol.vector.get()

        length = inputProtocol.length.get()
        pxSize = inputProtocol.inputParticles.get().getSamplingRate()

        subpartVectorList = load_vectors(cmmFn, vector,
                                         length,
                                         pxSize)

        vector = subpartVectorList[0]
        # coordinates with origin as defined in localrec
        # convert from string to array
        length = vector.get_length()
        x = vector.vector[0] * length
        y = vector.vector[1] * length
        z = vector.vector[2] * length
        # 0.0        107.3775 281.1066 pixel
        # get subvolume
        vol = self.inVolume.get()
        shifts = vol.getOrigin(force=True).getShifts()
        sampling = vol.getSamplingRate()

        # 100 100 100
        x += shifts[0]/sampling
        y += shifts[1]/sampling
        z += shifts[2]/sampling

        self.origin = emobj.Transform()
        samplingRate = samplingRate.get()
        self.origin.setShifts(x * samplingRate,
                         y * samplingRate,
                         z * samplingRate)
        self.outVol.setOrigin(self.origin)

        # copy or link
        copyOrLink = self.getCopyOrLink()
        fileName = self.inVol.getFileName()
        fileName = fileName.replace(':mrc', '')
        imgBase = basename(fileName)
        imgDst = self._getExtraPath(imgBase)
        copyOrLink(os.path.abspath(fileName), imgDst)
        self.outVol.setLocation(imgDst)

        # save
        self._defineOutputs(**{self.OUTPUTVOLUMENAME: self.outVol})
        self._defineSourceRelation(self.inVol, self.outVol)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        message = []
        # message.append("You must set a path to export.")
        # Nothing to append
        return message

    def _summary(self):
        # volume copy or linked
        message = []
        if self.setSampling.get():
            message.append("New Sampling: %f\n"%
                                   self.samplingRate)
        return message

    def _methods(self):
        return []

# --------------------------- UTILS functions ---------------------------------

    def getCopyOrLink(self):
        # Set a function to copyFile or createLink
        # depending in the user selected option
        if self.copyFiles:
            return pwutils.copyFile
        else:
            return pwutils.createAbsLink

    def getFnPath(self, label='volume'):
        return os.path.join(self.filesPath.get(),
                            self._getFileName(label))

