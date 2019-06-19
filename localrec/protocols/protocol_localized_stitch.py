# **************************************************************************
# *
# * Authors:         Vahid Abrishami (vahid.abrishami@helsinki.fi)
# *                  Juha Huiskonen (juha.huiskonen@helsinki.fi)
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
from pyworkflow.protocol.params import (EnumParam, IntParam, StringParam, BooleanParam,
                                        NumericRangeParam, PathParam, Positive, MultiPointerParam)
from pyworkflow.em.convert.transformations import euler_from_matrix
from pyworkflow.em.protocol import ProtPreprocessVolumes
from pyworkflow.em.data import *
from pyworkflow.protocol.constants import STEPS_PARALLEL

import xmippLib

from localrec.constants import CMM, LINEAR, symDict
from localrec.utils import load_vectors


class ProtLocalizedStich(ProtPreprocessVolumes):
    """
        Generate a full volume from a sub-volume applying a
        point group symmetry operation.

        An example of usage is to generate the adenovirus Capsid
        from its assymetric unit.
        """

    _label = 'stitch subparticles'

    LINEAR = 0
    BSPLINE = 1

    def __init__(self, *args,**kwargs):
        ProtPreprocessVolumes.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubVolumesHalf1', MultiPointerParam,
                      pointerClass='Volume',
                      label="Input sub-volume for half 1", important=True,
                      help='Select the input sub-volume for half 1')
        form.addParam('inputSubVolumesHalf2', MultiPointerParam,
                      pointerClass='Volume',
                      label="Input sub_volumes for half 2", important=True,
                      help='Select the input sub-volume for half 2')
        form.addParam('symMasks', MultiPointerParam, pointerClass='VolumeMask',
                      label='Masks', allowsNull=True,
                      help='Mask to normalize final volume. By default a sphere'
                           ' is used as the mask')
        form.addParam('interpMethod', EnumParam, choices=['linear', 'bspline'],
                      default=LINEAR,
                      label='Method for interpolation')
        form.addParam('outDim', IntParam,
                      label='Output volume size',
                      validators=[Positive],
                      help='This is size of the output volume after symmetrization')
        form.addParam('usePreRun', BooleanParam,
                      label="Use previous localrec run(s)", default=False)
        form.addParam('preRuns', MultiPointerParam, pointerClass='ProtLocalizedRecons',
                      label='Localrec previous runs', allowsNull=True,
                      condition="usePreRun",
                      help="This is the symmetirzed masked of the assymetric unit"
                           " which is used for overlapping normalization")
        form.addParam('symmetryGroup', StringParam, default='I1',
                       label="Symmetry", condition="not usePreRun",
                       help='There are multiple possibilities for '
                            'icosahedral symmetry: \n'
                            '* I1: No-Crowther 222 (standard in Heymann, '
                            'Chagoyen & Belnap, JSB, 151 (2005) 196-207)\n'
                            '* I2: Crowther 222 \n'
                            '* I3: 52-setting (as used in SPIDER?) \n'
                            '* I4: A different 52 setting \n')
        form.addParam('alignSubParticles', BooleanParam,
                      label="Sub-volumes are aligned?",condition="not usePreRun",
                      default=False,
                      help='If you aligned the sub-partilces with z-axis '
                           'to apply symmetry')

        group = form.addGroup('Vectors', condition="not usePreRun")
        group.addParam('defineVector', EnumParam, default=CMM,
                       label='Is vector defined by?',
                       choices=['cmm file', 'string'],
                       display=EnumParam.DISPLAY_HLIST)
        group.addParam('vector', NumericRangeParam, default='0,0,1',
                       label='Location vectors', condition="defineVector==1",
                       help='Vector defining the location of the '
                            'subparticles. The vector is defined by 3 '
                            'values x,y,z separated by comma. \n'
                            'More than one vector can be specified separated by'
                            'semicolon. For example: \n'
                            '0,0,1            # Defines only one vector.\n'
                            '0,0,1; 1,0,0;    # Defines two vectors.'
                       )
        group.addParam('vectorFile', PathParam, default='',
                       condition="defineVector==0",
                       label='file obtained by Chimera: ',
                       help='CMM file defining the location(s) of the '
                            'sub-particle(s). Use instead of vector. ')
        group.addParam('length', StringParam, default=-1,
                       label='Alternative length of the vector (A)',
                       help='Use to adjust the sub-particle center. If it '
                            'is <= 0, the length of the given vector. '
                            'Different values must be separated by commas.')

        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):


        inputStepId = self._insertFunctionStep('convertInputStep')
        depsSymVolHalf1 = []
        depsSymVolHalf2 = []
        depsSymMask = []
        doAlign = False

        if self.usePreRun:
            localRecProt = self.preRuns[0].get()
            localRecSymGrp = localRecProt.symGrp.get()
            localRecSymOrd = localRecProt.symmetryOrder.get()
            if localRecSymGrp == 0 or localRecSymGrp == 1:
                localRecSym = "%s%d" % (symDict[localRecSymGrp], localRecSymOrd)
            else:
                localRecSym = symDict[localRecSymGrp]
            doAlign = localRecProt.alignSubParticles
        else:
            localRecSym = self.symmetryGroup.get()
            doAlign = self.alignSubParticles

        for i, (vol, vol2) in enumerate(zip(self.inputSubVolumesHalf1, self.inputSubVolumesHalf2)):
            # Check if we have a mask for this volume
            maskFn = None
            volFn = vol.get().getFileName()
            volFn2 = vol2.get().getFileName()
            if len(self.symMasks) >= i+1:
                maskFn = self.symMasks[i].get().getFileName()

            maskVolumeHalf1Id = self._insertFunctionStep('maskVolume', volFn, maskFn,
                                                         i, 'half1',
                                                         prerequisites=[inputStepId])
            maskVolumeHalf2Id = self._insertFunctionStep('maskVolume', volFn2, maskFn,
                                                         i, 'half2',
                                                         prerequisites=[inputStepId])
            maskPreprationId = self._insertFunctionStep('prepareMask', i, doAlign,
                                                        prerequisites=[maskVolumeHalf1Id ,maskVolumeHalf2Id])
            volPreparationHalf1Id = self._insertFunctionStep('prepareVol',
                                                             i, 'half1', doAlign,
                                                             prerequisites=[maskVolumeHalf1Id])
            volPreparationHalf2Id = self._insertFunctionStep('prepareVol',
                                                         i, 'half2', doAlign,
                                                         prerequisites=[maskVolumeHalf2Id])
            symMaskStepId = self._insertFunctionStep('symmetrizeMask', i, localRecSym,
                                                     prerequisites=[maskPreprationId])

            symVolHalf1StepId = self._insertFunctionStep('symmetrizeVolume', i,
                                                         localRecSym, 'half1',
                                                         prerequisites=[maskPreprationId,
                                                                        volPreparationHalf1Id])
            symVolHalf2StepId = self._insertFunctionStep('symmetrizeVolume', i,
                                                         localRecSym, 'half2',
                                                         prerequisites=[maskPreprationId,
                                                                        volPreparationHalf2Id])
            depsSymVolHalf1.append(symVolHalf1StepId)
            depsSymVolHalf2.append(symVolHalf2StepId)
            depsSymMask.append(symMaskStepId)

        stitchStepHalf1Id = self._insertFunctionStep('stitchParticles', 'half1',
                                                     prerequisites=depsSymMask + depsSymVolHalf1)
        stitchStepHalf2Id = self._insertFunctionStep('stitchParticles', 'half2',
                                                     prerequisites=depsSymMask + depsSymVolHalf2)
        self._insertFunctionStep('createOutputStep', prerequisites=[stitchStepHalf1Id,
                                                                    stitchStepHalf2Id])

    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self):

        # Read voxel size
        self.pxSize = self.inputSubVolumesHalf1[0].get().getSamplingRate()
        interpMethod = self.interpMethod
        self.interpString = 'linear' if interpMethod.get() == LINEAR else 'spline'
        self.interpNum = 1 if interpMethod.get() == LINEAR else 3
        self.subVolCenterVec = []

        """ Compute the center of the sub-particle (the same as the vector
         for localized coordinate extraction)"""
        if self.usePreRun:
            for i, run in enumerate(self.preRuns):
                self.subVolCenterVec += self.createVector(self.preRuns[i].get())
        else:
            self.subVolCenterVec = self.createVector()

    def stitchParticles(self, halfString):

        # define the required file name to make a proper mask
        binarizedMaskFn = self._getTmpPath('binarized_mask_%s.vol' % halfString)
        erodedMaskFn = self._getTmpPath('eroded_mask_%s.vol' % halfString)
        softMaskFn = self._getTmpPath('soft_edge_mask_%s.vol' % halfString)
        volWithouMask = self._getTmpPath('vol_without_mask._%s.vol' % halfString)

        ih = ImageHandler()
        sumImg = ih.createImage()
        img = ih.createImage()
        sumMask = ih.createImage()
        imgMask = ih.createImage()
        finalMask = ih.createImage()
        volSym = self._getSymFn('volume', 1, halfString)
        maskSym = self._getSymFn('mask', 1)
        sumImg.read(volSym)
        sumMask.read(maskSym)
        sumImg.convert2DataType(ih.DT_DOUBLE)
        sumMask.convert2DataType(ih.DT_DOUBLE)

        # Loop over the halfX subvolumes and sum them up (and mask)
        for i, vol in enumerate(self.inputSubVolumesHalf1):
             if i==0:
                 continue
             volSym = self._getSymFn('volume', i+1, halfString)
             maskSym = self._getSymFn('mask', i+1)
             img.read(volSym)
             imgMask.read(maskSym)
             img.convert2DataType(ih.DT_DOUBLE)
             imgMask.convert2DataType(ih.DT_DOUBLE)
             sumImg.inplaceAdd(img)
             sumMask.inplaceAdd(imgMask)
        sumMaskData = sumMask.getData()
        sumMaskDataTmp = sumMask.getData()
        sumMaskData[sumMaskData<1.0] = 1.0
        sumMask.setData(sumMaskData)
        sumImg.inplaceDivide(sumMask)
        # This is to generate a smooth mask
        sumMaskDataTmp[sumMaskDataTmp<1.0] = 0.0
        sumMaskDataTmp[sumMaskDataTmp>=1.0] = 1.0
        finalMask.setData(sumMaskDataTmp)
        sumImg.write(volWithouMask)
        finalMask.write(binarizedMaskFn)
        # Apply morphology to mask to then apply softedge
        program = 'xmipp_transform_morphology'
        args = '-i %s --size 2.0 -o %s --binaryOperation erosion' % (binarizedMaskFn, erodedMaskFn)
        self.runJob(program,args)
        # Apply soft edge to the mask
        program = 'xmipp_transform_filter'
        args = '-i %s --fourier real_gaussian 2.0 -o %s' % (erodedMaskFn, softMaskFn)
        self.runJob(program,args)
        program = 'xmipp_image_operate'
        args = '-i %s --mult %s -o %s' % (volWithouMask, softMaskFn, self._getOutputVol(halfString))
        self.runJob(program,args)

    def maskVolume(self, volFn, maskFn, index, halfString):
        ih = ImageHandler()
        # Check if there is a mask if not then use a spherical mask

        if maskFn is None:
            # Compute sphere radius and make a spherical mask
            maskFn = self._getMaskFn(index+1)
            volSize = ih.getDimensions(volFn)[0]
            radius = volSize/2 - 1
            xmippLib.createEmptyFile(maskFn, volSize, volSize, volSize)
            program =  "xmipp_transform_mask"
            args = "-i %s --mask circular %d --create_mask %s " % (maskFn, -1 * radius, maskFn)
            self.runJob(program,args)
        # Apply mask to the volume
        if volFn.endswith(".mrc"):
            volFn = volFn + ":mrc"
        volMasked =  self._getMaskedVolFn(index+1, halfString)
        program = 'xmipp_image_operate'
        args = '-i %s --mult %s -o %s' % (volFn, maskFn, volMasked)
        self.runJob(program,args)

    def prepareMask(self, index, doAlign):

        shiftX, shiftY, shiftZ, rotMatrix = self.readVector(index)
        rot, tilt, psi = np.rad2deg(euler_from_matrix(rotMatrix.transpose(), 'szyz'))

        # Window the mask to the output volume
        maskWin = self._getWinFn('mask', index+1)
        maskFn = self._getMaskFn(index+1)
        program = 'xmipp_transform_window'
        args = '-i %s --size %d -o %s' % (maskFn, self.outDim, maskWin)
        self.runJob(program,args)

        #If sub-particles are aligned along z
        if doAlign:
            program = 'xmipp_transform_geometry'
            args = '-i %s --rotate_volume euler %f %f %f -o %s'\
                   % (maskWin, rot, tilt, psi, maskWin)
            self.runJob(program,args)

        # Shift the mask to the center of sub-volume
        maskShifted = self._getShiftedFn('mask', index+1)
        program = 'xmipp_transform_geometry'
        args = ('-i %s --shift %f %f %f -o %s --dont_wrap --interp %s'
               % (maskWin, shiftX, shiftY, shiftZ,
                  maskShifted, self.interpString))
        self.runJob(program,args)

    def prepareVol(self, index, halfString, doAlign):

        shiftX, shiftY, shiftZ, rotMatrix = self.readVector(index)
        rot, tilt, psi = np.rad2deg(euler_from_matrix(rotMatrix.transpose(), 'szyz'))

        # Window the sub-volume to the output volume size
        volWin = self._getWinFn('volume', index+1, halfString)
        volMasked =  self._getMaskedVolFn(index+1, halfString)

        program = 'xmipp_transform_window'

        args = '-i %s --size %d -o %s' % (volMasked, self.outDim, volWin)
        self.runJob(program,args)

        # If sub-particles are aligned along z
        if doAlign:
            program = 'xmipp_transform_geometry'
            args = '-i %s --rotate_volume euler %f %f %f -o %s'\
                   % (volWin, rot, tilt, psi, volWin)
            self.runJob(program,args)

        # Shift the sub-volume to its center in the volume
        volShifted = self._getShiftedFn('volume', index+1, halfString)
        program = 'xmipp_transform_geometry'
        args = ('-i %s --shift %f %f %f -o %s --dont_wrap --interp %s'
               % (volWin, shiftX, shiftY, shiftZ,
                  volShifted, self.interpString))
        self.runJob(program,args)

    def symmetrizeMask(self, index, localRecSym):

        volShifted = self._getShiftedFn('mask', index+1)
        # Apply symmetry operation to the mask
        maskSym = self._getSymFn('mask', index+1)
        program = 'xmipp_transform_symmetrize'
        args = '-i %s --sym %s -o %s --dont_wrap --sum --spline %d' % (volShifted, localRecSym, maskSym,  self.interpNum)
        self.runJob(program,args)

    def symmetrizeVolume(self, index, localRecSym, halfString):

        # Apply symmetry operation to the volume
        volShifted = self._getShiftedFn('volume', index+1, halfString)
        volSym = self._getSymFn('volume', index+1, halfString)
        program = 'xmipp_transform_symmetrize'
        args = '-i %s --sym %s -o %s --dont_wrap --sum --spline %d' % (volShifted, localRecSym, volSym, self.interpNum)
        self.runJob(program,args)

    def readVector(self, index):

        length = self.subVolCenterVec[index].get_length()
        shiftX = self.subVolCenterVec[index].vector[0] * length
        shiftY = self.subVolCenterVec[index].vector[1] * length
        shiftZ = self.subVolCenterVec[index].vector[2] * length
        rotMatrix = self.subVolCenterVec[index].get_matrix()
        return shiftX, shiftY, shiftZ, rotMatrix

    def createVector(self, protocolSitich=None):

        vector = ""
        cmmFn = ""
        if protocolSitich == None:
            protocolSitich = self
        if  protocolSitich.defineVector== CMM:
            cmmFn = protocolSitich.vectorFile.get()
        else:
            vector = protocolSitich.vector.get()
        return load_vectors(cmmFn, vector, protocolSitich.length.get(), self.pxSize)

    def createOutputStep(self):
        vol = self.inputSubVolumesHalf1[0]
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(self.pxSize)
        for halfString in ['half1', 'half2']:
            outVol = Volume()
            outVol.setFileName(self._getOutputVol(halfString))
            volumes.append(outVol)
        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(vol, volumes)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []

        if len(self.inputSubVolumesHalf1) != len(self.inputSubVolumesHalf1):
            validateMsgs.append("The number of sub-volumes for"
                                "half1 and half2 must be equal")
        if len(self.inputSubVolumesHalf1)>1:
            pxSize = self.inputSubVolumesHalf1[0].get().getSamplingRate()
            for i, vol in enumerate(self.inputSubVolumesHalf1):
                if i==0:
                    continue
                if pxSize != vol.get().getSamplingRate():
                    validateMsgs.append("The sampling rate of Volumes"
                                        " *MUST BE EQUAL*")

        if self.usePreRun:
            if len(self.inputSubVolumesHalf1) != len(self.preRuns):
                validateMsgs.append("You must assign each sub-volume"
                                    " a previous run of localrec")
        return validateMsgs

    def _summary(self):
        summary = ["Stitch %d sub-volumes to make a full volume of size %d"
                   % (len(self.inputSubVolumesHalf1), self.outDim)]
        return summary

    def _methods(self):
        messages = []
        return messages

    #--------------------------- UTILS functions -------------------------------
    def _getMaskFn(self, index):
        return self._getTmpPath('output_mask_%d.vol' % index)

    def _getOutputVol(self, halfString):
        return self._getExtraPath('output_volume_%s.vol' % halfString)

    def _getMaskedVolFn(self, index, halfString):
        return self._getTmpPath('output_volume_masked_%d_%s.vol' % (index, halfString))

    def _getWinFn(self, typeStr, index, halfString=''):
        if halfString == '':
            return self._getTmpPath('output_%s_windowed_%d.vol' % (typeStr, index))
        else:
            return self._getTmpPath('output_%s_windowed_%d_%s.vol' % (typeStr, index, halfString))

    def _getShiftedFn(self, typeStr, index, halfString=''):
        if halfString == '':
            return self._getTmpPath('output_%s_shifted_%d.vol' % (typeStr, index))
        else:
            return self._getTmpPath('output_%s_shifted_%d_%s.vol' % (typeStr, index, halfString))

    def _getSymFn(self, typeStr, index, halfString=''):
        if halfString == '':
            return self._getTmpPath('output_%s_symmetrized_%d.vol' % (typeStr, index))
        else:
            return self._getTmpPath('output_%s_symmetrized_%d_%s.vol' % (typeStr, index, halfString))
