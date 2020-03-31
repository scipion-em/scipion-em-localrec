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


class ProtLocalOrigSampling(EMProtocol):
    """ Set the origin and sampling values assigned to a 3D map
        so that the subvolume fits the original, larger volume
    """
    _label = 'Set origin to subvolume'
    _program = ""
    _lastUpdateVersion = VERSION_3_0

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
        self._defineOutputs(outputVolume=self.outVol)
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

