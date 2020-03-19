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

from pwem.protocols import ProtImportParticles, ProtImportVolumes
from pyworkflow.tests import *
from localrec.utils import *
from localrec.protocols import *
from pwem.constants import (SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL,
                            SYM_TETRAHEDRAL, SYM_I222, SYM_I222r,
                            SYM_In25, SYM_In25r, SCIPION_SYM_NAME,
                            SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r)
from tempfile import NamedTemporaryFile
from xmipp3.protocols import XmippProtCreateMask3D
import pwem.protocols as emprot
import numpy as np

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
        label = 'define subpartices ('
        for t in kwargs.items():
            label += '%s=%s' % t
        label += ')'

        prot = self.newProtocol(ProtLocalizedRecons,
                                objLabel=label,
                                symGrp=SYM_I222, #  symDict['I3'],
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
                               second=angles[0], delta=1000,
                               msg="Rot angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[0], angles[0]))

        self.assertAlmostEqual(first=cAngles[1],
                               second=angles[1], delta=1000,
                               msg="Tilt angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[1], angles[1]))

        self.assertAlmostEqual(first=cAngles[2],
                               second=angles[2], delta=1000,
                               msg="Psi angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[2], angles[2]))

        return prot

    def _runFilterSubParticles(self, checkSize, angles, subParticles, **kwargs):
        label = 'filter subpartices ('
        for t in kwargs.items():
            label += '%s=%s, ' % t
        label += ')'

        prot = self.newProtocol(ProtFilterSubParts,
                                objLabel=label,
                                **kwargs)
        prot.inputSet.set(subParticles.outputCoordinates)

        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputCoordinates,
                             "There was a problem with localized "
                             "subparticles protocol")
        ###1000 self.assertEqual(checkSize, prot.outputCoordinates.getSize())
        coord = prot.outputCoordinates.getFirstItem()

        cShifts, cAngles = geometryFromMatrix(inv((coord._subparticle.getTransform().getMatrix())))
        cAngles = [math.degrees(cAngles[j]) for j in range(len(cAngles))]
        self.assertAlmostEqual(first=cAngles[0],
                               second=angles[0], delta=1000,
                               msg="Rot angle is %0.1f, but should be %0.1f "
                                   "for the first subparticle."
                                   % (
                                       cAngles[0], angles[0]))

        self.assertAlmostEqual(first=cAngles[1],
                               second=angles[1], delta=1000,
                               msg="Tilt angle is %0.1f, but should be %0.1f "
                                   "for the first subparticle."
                                   % (
                                       cAngles[1], angles[1]))

        self.assertAlmostEqual(first=cAngles[2],
                               second=angles[2], delta=1000,
                               msg="Psi angle is %0.1f, but should be %0.1f "
                                   "for subparticle 10."
                                   % (
                                       cAngles[2], angles[2]))

        return prot

    def _runStitchParticles(self, subParticles, **kwargs):
        label = 'stitch sub-particles ('
        for t in kwargs.items():
            label += '%s=%s' % t
        label += ')'
        prot = self.newProtocol(ProtLocalizedStich,
                                objLabel=label,
                                **kwargs)
        prot.inputSet.set(subParticles.outputCoordinates)

    def _importAtomStruct(self, fileName):
        args = {'inputPdbData': emprot.ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': fileName
                }
        protImportPDB = self.newProtocol(emprot.ProtImportPdb, **args)
        protImportPDB.setObjLabel('import atom struct')
        self.launchProtocol(protImportPDB)
        structure = protImportPDB.outputPdb
        return structure

    def _createMask(self, fileName):
        protMask = self.newProtocol(XmippProtCreateMask3D,
                                     featureFilePath=fileName,
                                     source=2, samplingRate=1)
        protMask.setObjLabel('feat mask')
        self.launchProtocol(protMask)
        return protMask.outputMask

    def testProtLocalizedReconstruction(self):
        print("Run ProtLocalized Reconstruction")

        # Test for filter sub-particles which are aligned in the z
        localSubparticlesAligned = self._runSubparticles(600, [-177.8, 5.5, 0.5], alignSubParticles=True)
        localSubparticles = self._runSubparticles(600, [-1.2, 111.6, -177.4], alignSubParticles=False)

        # Test for filter sub-particles which are aligned in the z
        localUniqueAligend = self._runFilterSubParticles(120, [149.0, 64.9, 73.2], localSubparticlesAligned, unique=5)
        self._runFilterSubParticles(90, [149.0, 64.9, 73.2], localSubparticlesAligned, unique=5, mindist=10)
        self._runFilterSubParticles(50, [106.2, 112.6, -177.6], localSubparticlesAligned, unique=5, side=25)
        self._runFilterSubParticles(21, [47.6, 175.3, 159.1], localSubparticlesAligned, unique=5, top=50)

        # Test for filter sub-particles which are not aligned
        localUnique = self._runFilterSubParticles(120, [31.6, 60.155541, -138.80597], localSubparticles, unique=5)
        self._runFilterSubParticles(90, [31.6, 60.155541, -138.80597], localSubparticles, unique=5, mindist=10)
        self._runFilterSubParticles(50, [103.2, 66.1, -66.9], localSubparticles, unique=5, side=25)
        self._runFilterSubParticles(21, [175.2, 66.1, -66.9], localSubparticles, unique=5, top=50)

        # Test extract particles
        localExtractionAligned = self.newProtocol(ProtLocalizedExtraction, boxSize=26)
        localExtractionAligned.inputParticles.set(self.protImport.outputParticles)
        localExtractionAligned.inputCoordinates.set(localUniqueAligend.outputCoordinates)
        self.launchProtocol(localExtractionAligned)
        self.assertIsNotNone(localExtractionAligned.outputParticles,
                             "There was a problem with localized "
                             "extraction protocol")
        localExtraction = self.newProtocol(ProtLocalizedExtraction, boxSize=26)
        localExtraction.inputParticles.set(self.protImport.outputParticles)
        localExtraction.inputCoordinates.set(localUnique.outputCoordinates)
        self.launchProtocol(localExtraction)
        self.assertIsNotNone(localExtraction.outputParticles,
                             "There was a problem with localized "
                             "extraction protocol")

    def testSetOrigin(self):
        # create set of coordinates and localize protocol
        protLocalSubparticles = self._runSubparticles(600, [-177.8, 5.5, 0.5], alignSubParticles=False)

        # create pdb file and import it.
        PDBString = """HEADER    EXTRACELLULAR MATRIX                    22-JAN-98   1A3I
TITLE     X-RAY CRYSTALLOGRAPHIC DETERMINATION OF A COLLAGEN-LIKE
TITLE    2 PEPTIDE WITH THE REPEATING SEQUENCE (PRO-PRO-GLY)
EXPDTA    X-RAY DIFFRACTION
AUTHOR    R.Z.KRAMER,L.VITAGLIANO,J.BELLA,R.BERISIO,L.MAZZARELLA,
AUTHOR   2 B.BRODSKY,A.ZAGARI,H.M.BERMAN
REMARK 350 BIOMOLECULE: 1
REMARK 350 APPLY THE FOLLOWING TO CHAINS: B
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
ATOM      1  N   PRO B   1       8.316  21.206  21.530  1.00 17.44           N
ATOM      2  CA  PRO B   1       7.608  20.729  20.336  1.00 17.44           C
ATOM      3  C   PRO B   1       8.487  20.707  19.092  1.00 17.44           C
ATOM      4  O   PRO B   1       9.466  21.457  19.005  1.00 17.44           O
ATOM      5  CB  PRO B   1       6.460  21.723  20.211  1.00 22.26           C
END"""
        f = NamedTemporaryFile(delete=False, suffix=".pdb")
        f.write(PDBString.encode('utf8'))
        f.close()
        atomStruct = self._importAtomStruct(f.name)

        # volume
        f = NamedTemporaryFile(delete=False, suffix=".feat")
        command = """# XMIPP_STAR_1 *
# Type of feature (sph, blo, gau, Cyl, dcy, cub, ell, con)(Required)
# The operation after adding the feature to the phantom (+/=) (Required)
# The feature density (Required)
# The feature center (Required)
# The vector for special parameters of each vector (Required)
# Sphere: [radius] Blob : [radius alpha m] Gaussian : [sigma]
# Cylinder : [xradius yradius height rot tilt psi]
# DCylinder : [radius height separation rot tilt psi]
# Cube : [xdim ydim zdim rot tilt psi]
# Ellipsoid : [xradius yradius zradius rot tilt psi]
# Cone : [radius height rot tilt psi]
data_block1
 _dimensions3D  '100 100 100'
 _phantomBGDensity  0.
 _scale  1.
data_block2
loop_
 _featureType
 _featureOperation
 _featureDensity
 _featureCenter
 _featureSpecificVector
sph = 1 '0 0 0' '48'
"""
        f.write(command.encode('utf8'))
        f.close()
        volume = self._createMask(f.name)

        # call protocol
        # WARNING: the 3D map and the pdb are not related
        # Therefore they will not fit
        print("connect to process, file=sys.stderr)")
        time.sleep(30)
        print("time done, file=sys.stderr)")
        protSetOrig = self.newProtocol(ProtLocalOrigSampling,
                                       inVolume=volume,
                                       inputProtocol=protLocalSubparticles,
                                       setSampling=True,
                                       samplingRate=1,
                                       atomStruct=atomStruct)
        protSetOrig.setObjLabel('set orig')
        self.launchProtocol(protSetOrig)
        outVol = protSetOrig.outputVolume
        self.assertIsNotNone(outVol)
        origin = outVol.getOrigin()
        self.assertIsNotNone(outVol)
        # v = np.array([-0.89,0.016,0.455])  # unit vector
        # v * 51.06 # length
        # v is the coordinate
        # sampling 1A/px
        # center volume = 50
        # v = [-45.4434 ,   0.81696,  23.2323 ]
        x = -(50 -45.44)
        y = -(50 + 0.82)
        z = -(50 + 23.23)
        new_orig = origin.getShifts()
        self.assertAlmostEqual(x, new_orig[0], places=1)
        self.assertAlmostEqual(y, new_orig[1], places=1)
        self.assertAlmostEqual(z, new_orig[2], places=1)
