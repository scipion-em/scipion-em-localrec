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
    Assigns a corrected origin and optionally a new sampling rate to a
    reconstructed subvolume so that it is positioned consistently within
    the coordinate system of the original larger volume.

    This protocol is intended for localized reconstruction workflows in
    which a subvolume has been reconstructed from extracted subparticles.
    Because localized reconstructions are generated around a shifted
    subparticle center rather than the center of the original particle,
    the resulting map must be repositioned in the global reference frame
    of the parent volume before biological interpretation, fitting, or
    visualization.

    The protocol requires an input volume together with a
    ProtLocalizedRecons protocol that defines the subparticle vectors
    used during subparticle generation. These vectors determine where
    the localized reconstruction lies relative to the original particle.
    The first available vector is used as the reference position for
    computing the new origin.

    During execution, the input volume metadata is first copied into a
    new output volume object. The sampling rate is then assigned. By
    default, the protocol preserves the original sampling of the input
    volume, but the user may optionally define a new sampling rate.
    This new value is stored in the output volume object only and is
    not written into the map file header.

    To compute the new origin, the protocol reads the subparticle
    vector either from a Chimera CMM file or from manually defined
    vector coordinates. The vector length may optionally be adjusted.
    After conversion into pixel units using the sampling rate of the
    original particle images, the vector defines the displacement of
    the localized reconstruction relative to the parent particle
    center.

    The protocol then combines this subparticle displacement with the
    origin shifts already present in the input volume. These shifts are
    converted into pixel coordinates, added to the subparticle vector,
    and finally converted back into physical units using the output
    sampling rate. The resulting transformation is assigned as the new
    origin of the output map.

    The map file itself is not recomputed. Instead, the original volume
    file is either copied or linked into the protocol working
    directory, depending on the selected option. The output volume then
    points to this copied or linked file while carrying the corrected
    origin and sampling metadata.

    An optional atomic structure may also be associated with the
    protocol so that the corrected subvolume can be visualized together
    with structural models in downstream analysis.

    The output is therefore a volume whose voxel data remain unchanged,
    but whose metadata correctly place the localized reconstruction in
    the spatial context of the original full particle. This is
    particularly useful for interpretation of localized densities,
    fitting atomic models, and comparing subvolumes against the parent
    reconstruction.
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

