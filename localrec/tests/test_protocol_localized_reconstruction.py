# **************************************************************************
# *
# * Authors:   Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *            Vahid Abrishami (vahid.abrishami@helsinki.fi)
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
from numpy import inexact
from pyworkflow.em.protocol import ProtImportParticles, ProtImportVolumes
from pyworkflow.tests import *
from localrec.utils import *
from localrec.protocols import *

# Some utility functions to import micrographs that are used
# in several tests.
class TestLocalizedReconsBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_programs'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('input/Protocol_Projection_Match'
                                              'ing/ProjMatch/goldStandard/Iter'
                                              '_01/current_angles.xmd')
        cls.chimeraFile = cls.dataset.getFile('input/Protocol_Projection_Match'
                                              'ing/ico.cmm')

    @classmethod
    def runImportParticles(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportParticles,
                                     objLabel='from Xmipp ProjMatch',
                                     importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
                                     mdFile=pattern,
                                     magnification=65000,
                                     samplingRate=samplingRate,
                                     haveDataBeenPhaseFlipped=True
                                     )
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception(
                'Import of images: %s, failed. outputParticles is None.' % pattern)
        return protImport

#    @classmethod
#    def runImportVolumes(cls, pattern, samplingRate):
#        """ Run an Import particles protocol. """
#        protImport = cls.newProtocol(ProtImportVolumes,
#                                     filesPath=pattern,
#                                     samplingRate=samplingRate)
#        cls.launchProtocol(protImport)
#        return protImport


class TestLocalizedRecons(TestLocalizedReconsBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestLocalizedReconsBase.setData('xmipp_programs')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 1)

    #         cls.protImportVol = cls.runImportVolumes(cls.vol, 1)

    def _runSubparticles(self, checkSize, angles, defVector=0, **kwargs):
        label = 'localized subpartices ('
        for t in kwargs.iteritems():
            label += '%s=%s' % t
        label += ')'

        prot = self.newProtocol(ProtLocalizedRecons,
                                objLabel=label,
                                symmetryGroup2=SYM_In25,
                                defineVector=defVector,
                                **kwargs)
        prot.inputParticles.set(self.protImport.outputParticles)

        if defVector == 1:
            prot.vector.set('-0.89,0.016,0.455')
            prot.length.set(51.06)
        else:
            prot.vectorFile.set(self.chimeraFile)

        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputCoordinates,
                             "There was a problem with localized "
                             "subparticles protocol")
        self.assertEqual(checkSize, prot.outputCoordinates.getSize())

        coord = prot.outputCoordinates[10]
        cShifts, cAngles = geometryFromMatrix(inv((coord._subparticle.getTransform().getMatrix())))
        cAngles = [math.degrees(cAngles[j]) for j in range(len(cAngles))]

        self.assertAlmostEqual(first=cAngles[0],
                               second=angles[0], delta=0.1,
                               msg="Rot angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[0], angles[0]))

        self.assertAlmostEqual(first=cAngles[1],
                               second=angles[1], delta=0.1,
                               msg="Tilt angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[1], angles[1]))

        self.assertAlmostEqual(first=cAngles[2],
                               second=angles[2], delta=0.1,
                               msg="Psi angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[2], angles[2]))

        return prot

    def _runFilterSubParticles(self, checkSize, angles, subParticles, **kwargs):
        label = 'filter subpartices ('
        for t in kwargs.iteritems():
            label += '%s=%s' % t
        label += ')'

        prot = self.newProtocol(ProtFilterSubParts,
                                objLabel=label,
                                unique=5, **kwargs)
        prot.inputSet.set(subParticles.outputCoordinates)

        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputCoordinates,
                             "There was a problem with localized "
                             "subparticles protocol")
        self.assertEqual(checkSize, prot.outputCoordinates.getSize())
        coord = prot.outputCoordinates.getFirstItem()

        cShifts, cAngles = geometryFromMatrix(inv((coord._subparticle.getTransform().getMatrix())))
        cAngles = [math.degrees(cAngles[j]) for j in range(len(cAngles))]

        self.assertAlmostEqual(first=cAngles[0],
                               second=angles[0], delta=0.1,
                               msg="Rot angle is %0.1f, but should be %0.1f "
                                   "for the first subparticle."
                                   % (
                                       cAngles[0], angles[0]))

        self.assertAlmostEqual(first=cAngles[1],
                               second=angles[1], delta=0.1,
                               msg="Tilt angle is %0.1f, but should be %0.1f "
                                   "for the first subparticle."
                                   % (
                                       cAngles[1], angles[1]))

        self.assertAlmostEqual(first=cAngles[2],
                               second=angles[2], delta=0.1,
                               msg="Psi angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[2], angles[2]))

        return prot

    def testProtLocalizedReconstruction(self):
        print("Run ProtLocalized Reconstruction")

        # Test for filter sub-particles which are aligned in the z
        localSubparticles_aligned = self._runSubparticles(600, [-172.0, 5.5, -4.1], alignSubparticles=True)
        localSubparticles = self._runSubparticles(600, [-1.2, 111.6, -177.4], alignSubparticles=False)

        # Test for filter sub-particles which are aligned in the z
        self._runFilterSubParticles(90, [150.5, 64.6, 72.7], localSubparticles_aligned, mindist=10)
        self._runFilterSubParticles(50, [107.5, 112.2, -177.8], localSubparticles_aligned, side=25)
        self._runFilterSubParticles(21, [52.4, 175.0, 162.5], localSubparticles_aligned, top=50)

        # Test for filter sub-particles which are not aligned
        self._runFilterSubParticles(90, [31.6, 60.155541, -138.80597], localSubparticles, mindist=10)
        self._runFilterSubParticles(50, [103.2, 66.1, -66.9], localSubparticles, side=25)
        self._runFilterSubParticles(21, [175.2, 66.1, -66.9], localSubparticles, top=50)

        # Test extract particles
        localExtraction = self.newProtocol(ProtLocalizedExtraction, boxSize=26)
        localExtraction.inputParticles.set(self.protImport.outputParticles)
        localExtraction.inputCoordinates.set(localSubparticles.outputCoordinates)
        self.launchProtocol(localExtraction)
        self.assertIsNotNone(localExtraction.outputParticles,
                             "There was a problem with localized "
                             "extraction protocol")

# This test creates de volume instead of importing it
import pyworkflow.utils as pwutils
from pyworkflow.em.convert.symmetry import Icosahedron
from pyworkflow.em.constants import (SYM_I222r, SYM_In25, SCIPION_SYM_NAME)
import os

class TestLocalizedRecons2(TestLocalizedReconsBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        # create volume
        #    feat
        cls.volMapName = '/tmp/TestLocalizedRecons2.map'
        volFeatName = '/tmp/TestLocalizedRecons2.feat'
        cls.projName = "/tmp/proj.stk"
        cls.sym = SYM_I222r
        cls.createFeatVolume(volFeatName, cls.volMapName, sym=cls.sym)
        # create projections
        cls.createProj(cls.volMapName, cls.projName)
        # distance penton origin
        cls.radius = 55.
        # direction pentom
        # cls.pentonDir = None
        # import projections
        cls.protImport = cls.runImportParticles()

    @classmethod
    def createFeatVolume(cls, volFeatName, volMapName, sym=SYM_I222r):
        f = open(volFeatName, "w")
        f.write("""# Phantom description file, (generated with phantom help)
# General Volume Parameters:
#      Xdim      Ydim      Zdim   Background_Density Scale
       180 180 180 0 1.0
# Feature Parameters:
#Type  +/=  Density X_Center Y_Center Z_Center
""")
        icosahedron = Icosahedron(orientation=SCIPION_SYM_NAME[sym][1:])
        x = 0.;
        y = 0.;
        z = 0.
        f.write("# large sphere at the center\n")
        f.write("sph  + 1. %.3f %.3f %.3f 36.\n" % (x, y, z))
        f.write("# 5-fold\n")

        for i, vertice in enumerate(icosahedron.getVertices()):
            vertice = 55.0 * vertice
            f.write("sph  + 3 %.3f %.3f %.3f 8.25\n" %
                    (vertice[0], vertice[1], vertice[2]))
            if i==0:
                cls.pentonDir =  "%.3f, %.3f, %.3f"%(vertice[0], vertice[1], vertice[2])

        # print 3fold points
        f.write("# 3-fold\n")

        for _3fold in icosahedron.get3foldAxis():
            x, y, z = _3fold
            f.write("sph  + 0.8 %.3f %.3f %.3f 6.0\n" % (55.0 * x, 55.0 * y, 55.0 * z))

        # print 2fold points
        f.write("# 2-fold\n")
        for _2fold in icosahedron.get2foldAxis():
            x, y, z = _2fold
            f.write("sph  + 0.7 %.3f %.3f %.3f 3.0\n" %
                    (55.0 * x, 55.0 * y, 55.0 * z))
        f.close()
        #    map
        program = "xmipp_phantom_create"
        args= '-i {featFile} -o {mapFile}'.format(
            featFile=volFeatName, mapFile=volMapName)
        cls.__runXmippProgram(program, args)

    @classmethod
    def createProj(cls, volMapName, projStack='proj.stk'):
        #projection directions
        projDirections = [
            [],
            [0, 0, 0],
            [72.8262802662, 35.467688713, 115.701733501], # keep this angles under 180
                                                          # or the function that conver from matrix to angles
                                                          # will report a different pair
            [298.153046614, 53.460868127, 100.044536396],
            [55.9545190927, 2.63201459146, 236.784245424],
            [359.783825617, 347.692361404, 4.866600415],
            [318.280094133, 78.7595410907, 343.995821614],
            [325.505603341, 85.297516634, 288.615732761],
            [137.11561079, 3.14867893078, 17.8840386505],
            [5.09868873239, 11.9211673352, 202.7125587],
            [210.62035136, 157.802079063, 129.961496962]
        ]
        ## If proj file exists, delete projections
        if os.path.isfile(projStack):
            os.remove(projStack)

        f = open(projStack.replace(".stk", ".xmd"), "w")
        f.write("""# XMIPP_STAR_1 *
data_noname
loop_
 _image
 _enabled
 _angleRot
 _angleTilt
 _anglePsi
 _shiftX
 _shiftY
 _micrograph
 _micrographId
""")

        command = "xmipp_phantom_project"
        for index in range(1,11):
            phi, theta, psi = projDirections[index]
            f.write("{id} 1 {phi} {theta} {psi} 0. 0. {micrograph} {micrographId}\n".format(id="%04d@proj.stk" % index,
                                                                                            phi=phi, theta=theta,
                                                                                            psi=psi,
                                                                                            micrograph="mic%04d.mrc" % (
                                                                                                        index / 5),
                                                                                            micrographId=index / 5))
            args = "-i {volMapName} -o {index}@{projStack} --angles {phi} {theta} {psi}".\
                format(volMapName=volMapName, projStack=projStack, index=index,
                       phi=phi, theta=theta, psi=psi)
            cls.__runXmippProgram(command, args)

    @classmethod
    def __runXmippProgram(cls, program, args):
        """ Internal shortcut function to launch a Xmipp program.
        If xmipp not available o fails return False, else Tru"""
        try:
            xmipp3 = pwutils.importFromPlugin('xmipp3', doRaise=True)
            xmipp3.Plugin.runXmippProgram(program, args)
        except:
            return False
        return True

    @classmethod
    def runImportParticles(cls):
        pattern=cls.projName.replace(".stk", ".xmd")
        samplingRate=1
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportParticles,
                                     objLabel='import particles',
                                     importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
                                     mdFile=pattern,
                                     magnification=65000,
                                     samplingRate=samplingRate,
                                     haveDataBeenPhaseFlipped=True
                                     )
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception(
                'Import of images: %s, failed. outputParticles is None.' % pattern)
        return protImport

    def _runSubparticles(self, checkSize, angles, defVector=1, **kwargs):
        label = 'localized subpartices ('
        for t in kwargs.iteritems():
            label += '%s=%s' % t
        label += ')'

        prot = self.newProtocol(ProtLocalizedRecons,
                                objLabel=label,
                                symmetryGroup2=self.sym,
                                defineVector=defVector,
                                **kwargs
                                )
        prot.inputParticles.set(self.protImport.outputParticles)

        if defVector == 1:
            prot.vector.set(self.pentonDir)
            prot.length.set(self.radius)
        else:
            prot.vectorFile.set(self.chimeraFile)

        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputCoordinates,
                             "There was a problem with localized "
                             "subparticles protocol")
        self.assertEqual(checkSize, prot.outputCoordinates.getSize())

        coord = prot.outputCoordinates[61]
        print "coord", coord
        cShifts, cAngles = geometryFromMatrix(inv((coord._subparticle.getTransform().getMatrix())))
        cAngles = [math.degrees(cAngles[j]) for j in range(len(cAngles))]
        print "cAngles", cAngles

        self.assertAlmostEqual(first=cAngles[0],
                               second=angles[0], delta=0.1,
                               msg="Rot angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[0], angles[0]))

        self.assertAlmostEqual(first=cAngles[1],
                               second=angles[1], delta=0.1,
                               msg="Tilt angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[1], angles[1]))

        self.assertAlmostEqual(first=cAngles[2],
                               second=angles[2], delta=0.1,
                               msg="Psi angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[2], angles[2]))

        return prot

    def testProtA(self):
        # import projection
        # angles = [210.62035136, 157.802079063, 129.961496962]
        angles = [72.8262802662, 35.467688713, 115.701733501]
        self._runSubparticles(600, angles , defVector=1, alignSubparticles=False)


