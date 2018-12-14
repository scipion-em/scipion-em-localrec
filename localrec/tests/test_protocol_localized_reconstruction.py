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

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportVolumes,
                                     filesPath=pattern,
                                     samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        return protImport


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
                                symmetryGroup="I3",
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
