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
from pwem.emlib import DT_DOUBLE
from pyworkflow.protocol.params import (EnumParam, IntParam, StringParam, BooleanParam,
                                        NumericRangeParam, PathParam, Positive, MultiPointerParam)
from pwem.convert.transformations import euler_from_matrix
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtPreprocessVolumes
from pwem.objects.data import *
from pyworkflow.protocol.constants import STEPS_PARALLEL
import pyworkflow.utils.path as putils

from pwem import emlib

from localrec.constants import CMM, LINEAR, symDict
from localrec.utils import load_vectors


class ProtLocalizedStich(ProtPreprocessVolumes):
    """
        Generate a full volume from a sub-volume applying a
        point group symmetry operation.

        An example of usage is to generate the adenovirus capsid
        from its asymmetric unit.
        """

    _label = 'stitch subvolumes'

    LINEAR = 0
    BSPLINE = 1

    def __init__(self, **kwargs):
        ProtPreprocessVolumes.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('useHalMaps', BooleanParam,
                      label="Use two half maps?",
                      default=False,
                      help='Stitch half-maps separately')
        form.addParam('inputSubVolumes', MultiPointerParam,
                      pointerClass='Volume', condition="not useHalMaps",
                      label="Input sub-volumes ", allowsNull=True,
                      help='Select input sub-volume for stitching')
        form.addParam('inputSubVolumesHalf1', MultiPointerParam,
                      pointerClass='Volume', condition="useHalMaps",
                      label="Input sub-volume for half-map 1", allowsNull=True,
                      help='Select the input sub-volume for half-map 1')
        form.addParam('inputSubVolumesHalf2', MultiPointerParam,
                      pointerClass='Volume', condition="useHalMaps",
                      label="Input sub_volumes for half-map 2", allowsNull=True,
                      help='Select the input sub-volume for half-map 2')
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
                      help="Previous Localrec runs used to extract the parameters")
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
                      label="Sub-volumes are aligned?", condition="not usePreRun",
                      default=False,
                      help='Set to Yes if you aligned the sub-particles with the z-axis '
                           'earlier. Note that the you can mix sub-particles with and '
                           'without this additional alignment. ')

        group = form.addGroup('Vectors', condition="not usePreRun")
        group.addParam('defineVector', EnumParam, default=CMM,
                       label='Is vector defined by?',
                       choices=['cmm file', 'string'],
                       display=EnumParam.DISPLAY_HLIST)
        group.addParam('vector', NumericRangeParam, default='0,0,1',
                       label='Location vectors', condition="defineVector==1",
                       help='Vector defining the location of the '
                            'sub-particles. The vector is defined by 3 '
                            'values x,y,z separated by comma. \n'
                            'More than one vector can be specified separated by a '
                            'semicolon. For example: \n'
                            '0,0,1            # Defines only one vector.\n'
                            '0,0,1; 1,0,0;    # Defines two vectors.'
                       )
        group.addParam('vectorFile', PathParam, default='',
                       condition="defineVector==0",
                       label='file obtained by Chimera: ',
                       help='CMM file defining the location(s) of the '
                            'sub-particle(s). Use instead of a vector. ')
        group.addParam('length', StringParam, default=-1,
                       label='Alternative length of the vector (A)',
                       help='Use to adjust the sub-particle center. If it '
                            'is <= 0, the length of the given vector is used. '
                            'Multiple values must be separated by commas.')

        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):

        inputStepId = self._insertFunctionStep('convertInputStep')
        depsSymVolHalf1 = []
        depsSymVolHalf2 = []
        depsSymVol = []
        depsSymMask = []

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

        if self.useHalMaps:
            for i, (vol, vol2) in enumerate(zip(self.inputSubVolumesHalf1, self.inputSubVolumesHalf2)):
                # Check if we have a mask for this volume
                maskFn = None
                volFn = vol.get().getFileName()
                volFn2 = vol2.get().getFileName()
                if len(self.symMasks) >= i+1:
                    maskFn = self.symMasks[i].get().getFileName()

                maskVolumeHalf1Id = self._insertFunctionStep('maskVolume', volFn,
                                                        maskFn,i, 'half1',
                                                        prerequisites=[inputStepId])
                maskVolumeHalf2Id = self._insertFunctionStep('maskVolume', volFn2,
                                                             maskFn, i, 'half2',
                                                             prerequisites=[inputStepId])
                maskPreprationId = self._insertFunctionStep('prepareMask', i, doAlign,
                                                            prerequisites=[maskVolumeHalf1Id, maskVolumeHalf2Id])
                volPreparationHalf1Id = self._insertFunctionStep('prepareVol',
                                                                 i, 'half1', doAlign,
                                                                 prerequisites=[maskVolumeHalf1Id])
                volPreparationHalf2Id = self._insertFunctionStep('prepareObj',
                                                                 i, 'half2', doAlign, 'volume',
                                                                 prerequisites=[maskVolumeHalf2Id])

                depsSymVolHalf1.append(volPreparationHalf1Id)
                depsSymVolHalf2.append(volPreparationHalf2Id)
                depsSymMask.append(maskPreprationId)

            genAsymUnitHalf1Id = self._insertFunctionStep('genAsymUnit','half1',
                                                          prerequisites=depsSymVolHalf1 + depsSymMask)
            genAsymUnitHalf2Id = self._insertFunctionStep('genAsymUnit','half2',
                                                          prerequisites=depsSymVolHalf2 + depsSymMask)

            symMaskStepId = self._insertFunctionStep('symmetrizeObj',
                                                     localRecSym, '', 'mask',
                                                     prerequisites=[genAsymUnitHalf1Id, genAsymUnitHalf2Id])
            symVolHalf1StepId = self._insertFunctionStep('symmetrizeObj',
                                                         localRecSym, 'half1', 'volume',
                                                         prerequisites=[genAsymUnitHalf1Id])
            symVolHalf2StepId = self._insertFunctionStep('symmetrizeObj',
                                                         localRecSym, 'half2', 'volume',
                                                         prerequisites=[genAsymUnitHalf2Id])

            stitchStepHalf1Id = self._insertFunctionStep('stitchParticles', 'half1',
                                                        prerequisites=[symMaskStepId, symVolHalf1StepId])
            stitchStepHalf2Id = self._insertFunctionStep('stitchParticles', 'half2',
                                                        prerequisites=[symMaskStepId, symVolHalf2StepId])
            self._insertFunctionStep('createOutputStep', prerequisites=[stitchStepHalf1Id,
                                                                        stitchStepHalf2Id])
        else:
            for i, vol in enumerate(self.inputSubVolumes):
                maskFn = None
                volFn = vol.get().getFileName()
                if len(self.symMasks) >= i+1:
                    maskFn = self.symMasks[i].get().getFileName()

                #Generate mask and apply it to the volume
                maskVolumeId = self._insertFunctionStep('maskVolume', volFn,
                                                        maskFn, i, '',
                                                        prerequisites=[inputStepId])
                maskPreprationId = self._insertFunctionStep('prepareObj',
                                                            i, '', doAlign, 'mask',
                                                            prerequisites=[maskVolumeId])
                volPreparationId = self._insertFunctionStep('prepareObj',
                                                            i, '', doAlign, 'volume',
                                                            prerequisites=[maskVolumeId])

                depsSymVol.append(volPreparationId)
                depsSymMask.append(maskPreprationId)

            genAsymUnitId = self._insertFunctionStep('genAsymUnit','',
                                                      prerequisites=depsSymVol + depsSymMask)
            symMaskStepId = self._insertFunctionStep('symmetrizeObj',
                                                     localRecSym, '', 'mask',
                                                     prerequisites=[genAsymUnitId])
            symVolStepId = self._insertFunctionStep('symmetrizeObj',
                                                    localRecSym, '', 'volume',
                                                    prerequisites=[genAsymUnitId])
            stitchStepId = self._insertFunctionStep('stitchParticles', '',
                                                    prerequisites=[symMaskStepId, symVolStepId])
            self._insertFunctionStep('createOutputStep', prerequisites=[stitchStepId])

    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self):

        # Read voxel size
        if self.useHalMaps:
            self.pxSize = self.inputSubVolumesHalf1[0].get().getSamplingRate()
        else:
            self.pxSize = self.inputSubVolumes[0].get().getSamplingRate()
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

    def genAsymUnit(self, halfString):

        # Create objects to handle images
        ih = ImageHandler()
        sumImg = ih.createImage()
        img = ih.createImage()
        sumMask = ih.createImage()
        imgMask = ih.createImage()
        # Read first mask and volume and add it to sum
        maskShifted = self._getFileName('mask', 'shifted', 1)
        volShifted = self._getFileName('volume', 'shifted', 1, halfString)
        sumImg.read(volShifted)
        sumMask.read(maskShifted)
        sumImg.convert2DataType(DT_DOUBLE)
        sumMask.convert2DataType(DT_DOUBLE)

         # Loop over the halfX subvolumes and sum them up (and mask)
        if not self.useHalMaps:
            listObj = self.inputSubVolumes
        else:
            listObj = self.inputSubVolumesHalf1
        for i, vol in enumerate(listObj):
             if i==0:
                 continue
             # Read current volume and mask
             volShifted = self._getFileName('volume', 'shifted', i+1, halfString)
             maskShifted = self._getFileName('mask', 'shifted', i+1)
             img.read(volShifted)
             imgMask.read(maskShifted)
             img.convert2DataType(DT_DOUBLE)
             imgMask.convert2DataType(DT_DOUBLE)
             # Add vol and mask to related sums
             sumImg.inplaceAdd(img)
             sumMask.inplaceAdd(imgMask)
        # Write mask and volume after sum on disk
        sumMask.write(self._getFileName('mask', 'sum'))
        sumImg.write(self._getFileName('volume', 'sum', -1, halfString))

    def stitchParticles(self, halfString):

        # define the required file name to make a proper mask
        binarizedMaskFn = self._getFileName('mask', 'binarized', -1, halfString)
        erodedMaskFn = self._getFileName('mask', 'eroded', -1, halfString)
        softMaskFn = self._getFileName('mask', 'soft_edge', -1, halfString)
        volWithouMask = self._getFileName('volume', 'without_mask', -1, halfString)
        outputVol = self._getOutputFileName(halfString)

        # Read symmetrized volume and mask
        ih = ImageHandler()
        maskSymImg = ih.createImage()
        volSymImg = ih.createImage()
        finalMask = ih.createImage()
        volSymFn = self._getFileName('volume', 'symmetrized', -1, halfString)
        maskSymFn = self._getFileName('mask', 'symmetrized')
        volSymImg.read(volSymFn)
        maskSymImg.read(maskSymFn)
        volSymImg.convert2DataType(DT_DOUBLE)
        maskSymImg.convert2DataType(DT_DOUBLE)

        # Here we divide volume by mask
        maskSymData = maskSymImg.getData()
        sumMaskDataTmp = maskSymImg.getData()
        maskSymData[maskSymData<1.0] = 1.0
        maskSymImg.setData(maskSymData)
        volSymImg.inplaceDivide(maskSymImg)
        volSymImg.write(volWithouMask)

        # This is to generate a smooth mask
        sumMaskDataTmp[sumMaskDataTmp < 1.0] = 0.0
        sumMaskDataTmp[sumMaskDataTmp >= 1.0] = 1.0
        finalMask.setData(sumMaskDataTmp)
        finalMask.write(binarizedMaskFn)

        # Apply morphology to mask to then apply softedge
        program = 'xmipp_transform_morphology'
        args = '-i %s --size 2.0 -o %s --binaryOperation erosion' % (binarizedMaskFn, erodedMaskFn)
        self.runJob(program, args)
        # Apply soft edge to the mask
        program = 'xmipp_transform_filter'
        args = '-i %s --fourier real_gaussian 2.0 -o %s' % (erodedMaskFn, softMaskFn)
        self.runJob(program, args)
        program = 'xmipp_image_operate'
        args = '-i %s --mult %s -o %s' % (volWithouMask, softMaskFn, outputVol)
        self.runJob(program,args)

    def maskVolume(self, volFn, maskFn, index, halfString):
        ih = ImageHandler()

        # Check if there is a mask if not then use a spherical mask
        if maskFn is None:
            # Compute sphere radius and make a spherical mask
            maskFn = self._getFileName('mask', index+1)
            volSize = ih.getDimensions(volFn)[0]
            radius = volSize/2 - 1
            emlib.createEmptyFile(maskFn, volSize, volSize, volSize)
            program = "xmipp_transform_mask"
            args = "-i %s --mask circular %d --create_mask %s " % (maskFn, -1 * radius, maskFn)
            self.runJob(program,args)
        else:
            putils.copyFile(maskFn, self._getFileName('mask', index+1))
            maskFn = self._getFileName('mask', index+1)

        # Read the volume if it is provided
        if volFn.endswith(".mrc"):
            volFn = volFn + ":mrc"
        # Apply mask to the volume
        volMasked = self._getFileName('volume', 'masked', index+1, halfString);
        program = 'xmipp_image_operate'
        args = '-i %s --mult %s -o %s' % (volFn, maskFn, volMasked)
        self.runJob(program, args)

    def prepareObj(self, index, halfString, doAlign, objType):

        shiftX, shiftY, shiftZ, rotMatrix = self.readVector(index)
        rot, tilt, psi = np.rad2deg(euler_from_matrix(rotMatrix.transpose(), 'szyz'))

        # Window the sub-volume to the output volume size
        objWin = self._getFileName(objType, 'windowed', index+1, halfString)
        if objType == 'mask':
            objMasked = self._getFileName('mask', index+1)
        else:
            objMasked = self._getFileName(objType, 'masked', index+1, halfString)

        program = 'xmipp_transform_window'
        args = '-i %s --size %d -o %s' % (objMasked, self.outDim, objWin)
        self.runJob(program,args)

        # If sub-particles are aligned along z
        if doAlign:
            program = 'xmipp_transform_geometry'
            args = '-i %s --rotate_volume euler %f %f %f -o %s' \
                   % (objWin, -rot, -tilt, -psi, objWin)
            self.runJob(program,args)

        # Shift the sub-volume to its center in the volume
        objShifted = self._getFileName(objType, 'shifted', index+1, halfString)
        program = 'xmipp_transform_geometry'
        args = ('-i %s --shift %f %f %f -o %s --dont_wrap --interp %s'
                % (objWin, shiftX, shiftY, shiftZ,
                   objShifted, self.interpString))
        self.runJob(program,args)

    def symmetrizeObj(self, localRecSym, halfString, objType):

        # Apply symmetry operation to the volume
        objSum = self._getFileName(objType, 'sum', -1, halfString)
        objSym = self._getFileName(objType, 'symmetrized', -1, halfString)
        program = 'xmipp_transform_symmetrize'
        args = '-i %s --sym %s -o %s --dont_wrap --sum --spline %d' % (objSum, localRecSym, objSym, self.interpNum)
        self.runJob(program,args)

    def readVector(self, index):

        length = self.subVolCenterVec[index].get_length()
        [shiftX, shiftY, shiftZ] = [x * length for x in self.subVolCenterVec[index].vector]
        rotMatrix = self.subVolCenterVec[index].get_matrix()
        return shiftX, shiftY, shiftZ, rotMatrix

    def createVector(self, protocolSitich=None):

        vector = ""
        cmmFn = ""
        if protocolSitich is None:
            protocolSitich = self
        if protocolSitich.defineVector == CMM:
            cmmFn = protocolSitich.vectorFile.get()
        else:
            vector = protocolSitich.vector.get()
        return load_vectors(cmmFn, vector, protocolSitich.length.get(),
                            self.pxSize)

    def createOutputStep(self):
        if self.useHalMaps:
            vol = self.inputSubVolumesHalf1[0]
            volumes = self._createSetOfVolumes()
            volumes.setSamplingRate(self.pxSize)
            for halfString in ['half1', 'half2']:
                outVol = Volume()
                outputVolFn = self._getOutputFileName(halfString)
                outVol.setFileName(outputVolFn)
                volumes.append(outVol)
            self._defineOutputs(outputVolumes=volumes)
            self._defineSourceRelation(vol, volumes)
        else:
            vol = self.inputSubVolumes[0]
            outVol = Volume()
            outputVolFn = self._getOutputFileName()
            outVol.setSamplingRate(self.pxSize)
            outVol.setFileName(outputVolFn)
            self._defineOutputs(outputVolume=outVol)
            self._defineSourceRelation(vol,outVol)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):

        validateMsgs = []
        if self.useHalMaps:
            if len(self.inputSubVolumesHalf1) != len(self.inputSubVolumesHalf2):
                validateMsgs.append("The number of sub-volumes for"
                                    "half1 and half2 must be equal")
            listObj = self.inputSubVolumesHalf1
        else:
            listObj = self.inputSubVolumes
        if len(listObj)>1:
            pxSize = listObj[0].get().getSamplingRate()
            for i, vol in enumerate(listObj):
                if i==0:
                    continue
                if pxSize != vol.get().getSamplingRate():
                    validateMsgs.append("The sampling rate of Volumes"
                                        " *MUST BE EQUAL*")
        if self.usePreRun:
            if len(listObj) != len(self.preRuns):
                validateMsgs.append("You must assign each sub-volume"
                                    " a previous run of localrec")
        return validateMsgs

    def _citations(self):
        return ['Ilca2015', 'Abrishami2020']

    def _summary(self):
        listObj = self.inputSubVolumesHalf1
        if not self.useHalMaps:
            listObj = self.inputSubVolumes
        summary = ["Stitch %d sub-volumes to make a full volume of size %d"
                   % (len(listObj), self.outDim)]
        return summary

    def _methods(self):
        messages = []
        return messages

    #--------------------------- UTILS functions -------------------------------
    def _getFileName(self, imgType, desc='', index=-1, halfString=''):
        auxString = '' if halfString == '' else '_{}'.format(halfString)
        auxString2 = '' if index == -1 else '_{}'.format(index)
        auxString3 = '' if desc == '' else '_{}'.format(desc)
        return self._getTmpPath('output_%s%s%s%s.vol'
                                  % (imgType, auxString3,
                                     auxString2, auxString))
    def _getOutputFileName(self, halfString=''):
        auxString = '' if halfString == '' else '_{}'.format(halfString)
        return self._getExtraPath('output_volume%s.vol'
                                  % auxString)
