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

from pwem.protocols import (ProtImportParticles, ProtImportVolumes,
                            ProtImportPdb)
from pyworkflow.tests import (BaseTest, DataSet, setupTestProject)
from localrec.utils import (geometryFromMatrix, inv)
from localrec.protocols import (ProtLocalizedRecons, ProtFilterSubParts,
                                ProtLocalizedStich, ProtLocalizedExtraction,
                                ProtLocalOrigSampling)
from pwem.constants import (SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL,
                            SYM_TETRAHEDRAL, SYM_I222, SYM_I222r,
                            SYM_In25, SYM_In25r, SCIPION_SYM_NAME,
                            SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r,
                            SYM_HELICAL)
from tempfile import NamedTemporaryFile
from xmipp3.protocols import XmippProtCreateMask3D
import math
import os
from tempfile import mkstemp
from pwem.convert.symmetry import Icosahedron
import time

import numpy as np
# from Bio.PDB import PDBParser
from pwem.constants import (SCIPION_SYM_NAME)
from xmipp3.constants import (XMIPP_SYM_NAME, XMIPP_TO_SCIPION, XMIPP_CYCLIC,
                              XMIPP_DIHEDRAL_X, XMIPP_TETRAHEDRAL,
                              XMIPP_OCTAHEDRAL,
                              XMIPP_I222, XMIPP_I222r, XMIPP_In25, XMIPP_In25r)
from xmipp3 import Plugin
from pyworkflow.utils import runJob
# Some utility functions to import micrographs that are used
# in several tests.


# Global cache dictionary
_cache = {}


class TestLocalizedReconsBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_programs'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('input/Protocol_Projection_Match'
                                              'ing/ProjMatch/goldStandard/Iter'
                                              '_01/current_angles.xmd')
        cls.chimeraFile = cls.dataset.getFile('input/Protocol_Projection_Match'
                                              'ing/ico.cmm')
        from pwem import Domain
        try:
            cls.xmipp3 = \
                Domain.importFromPlugin('xmipp3', doRaise=True)
            cls.xmippAvailable = True
        except ImportError:
            cls.xmippAvailable = False

    @classmethod
    def runImportParticles(cls, pattern, samplingRate,
                           objLabel='from Xmipp project'):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(
            ProtImportParticles,
            importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
            mdFile=pattern,
            magnification=65000,
            samplingRate=samplingRate,
            haveDataBeenPhaseFlipped=True
            )
        protImport.setObjLabel(objLabel)
        cls.launchProtocol(protImport)
        # check that input images have been imported (a better way to do this?)
        if protImport.outputParticles is None:
            raise Exception(
                'Import of images: %s, failed. outputParticles is None.' % pattern)
        return protImport

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate,
                         objLabel='from Xmipp createvol'):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(
            ProtImportVolumes,
            filesPath=pattern,
            samplingRate=samplingRate,
            )
        protImport.setObjLabel(objLabel)
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def chimeraOperate(cls, inputVolume, inputAtomStruct,
                       objLabel='Chimera operate'):
        try:
            from chimera.protocols import ChimeraProtOperate
            # do not open chimera GUI
            extraCommands = "scipionwrite #2 " \
                            "prefix DONOTSAVESESSION\n"
            extraCommands += "scipionwrite #3 " \
                             "prefix DONOTSAVESESSION\n"
            extraCommands += "exit\n"

            args = {'extraCommands': extraCommands,
                    'inputVolume': inputVolume,
                    'pdbFileToBeRefined': inputAtomStruct
                    }
            protChimera = cls.newProtocol(ChimeraProtOperate,
                                          **args)
            protChimera.setObjLabel(objLabel)
            cls.launchProtocol(protChimera)
            return protChimera
        except ImportError:
            cls.fail("Could not import ChimeraProtOperate")

    @classmethod
    def runReconstructRelion(cls, inputParticles, symmetryGroup='C1',
                             numberOfMpis=4,
                             objLabel='Fourier reconstruction'):
        """ Run an Import particles protocol. """
        try:
            from relion.protocols import ProtRelionReconstruct
            prot = cls.newProtocol(
                ProtRelionReconstruct,
                symmetryGroup=symmetryGroup,
                numberOfMpis=numberOfMpis,
                objLabel=objLabel,
                inputParticles=inputParticles)
            _ = cls.launchProtocol(prot)
            return prot
        except ImportError:
            cls.fail("Could not import ProtRelionReconstruct")

    @classmethod
    def createFeatFile(cls, featDict, dim, bgDensity=0., scale=1.):
        f = NamedTemporaryFile(delete=False, suffix=".feat")
        command = f"""# XMIPP_STAR_1 *
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
 _dimensions3D '{dim} {dim} {dim}'
 _phantomBGDensity {bgDensity}
 _scale {scale}
data_block2
loop_
 _featureType
 _featureOperation
 _featureDensity
 _featureCenter
 _featureSpecificVector
"""
        f.write(command.encode('utf8'))
        for feat in featDict:
            f.write(f"{feat['type']} {feat['operation']} "
                    f"{feat['density']} '{feat['center']}' "
                    f"'{feat['specific']}'\n".encode('utf8'))
        f.close()
        return f.name


class TestLocalizedRecons(TestLocalizedReconsBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestLocalizedReconsBase.setData('xmipp_programs')
        cls.protImport = cls.runImportParticles(
            cls.particlesFn, 1, objLabel='Import particles from Xmipp project')

    #         cls.protImportVol = cls.runImportVolumes(cls.vol, 1)

    def _runSubparticles(self,
                         checkSize,
                         angles,
                         defVector=1,
                         vector='-0.89,0.016,0.455',
                         length=51.06,
                         alignSubParticles=False,
                         symGroup=6,
                         symmetryOrder=1,
                         label='define subparticles',
                         inputParticles=None,
                         **kwargs):
        if inputParticles is None:
            inputParticles = self.protImport.outputParticles
        prot = self.newProtocol(ProtLocalizedRecons,
                                objLabel=label,
                                symGrp=symGroup,
                                symmetryOrder=symmetryOrder,
                                alignSubParticles=alignSubParticles,
                                defineVector=defVector,
                                **kwargs)
        prot.setObjLabel(label)
        prot.inputParticles.set(inputParticles)

        if defVector == 1:
            prot.vector.set(vector)
            prot.length.set(length)
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

    def _runFilterSubParticles(self, checkSize, angles, subParticles,
                               label="filter particles", **kwargs):
        label = 'filter subpartices ('
        for t in kwargs.items():
            label += '%s=%s, ' % t
        label += ')'

        prot = self.newProtocol(ProtFilterSubParts,
                                **kwargs)
        prot.setObjLabel(label)
        prot.inputSet.set(subParticles.outputCoordinates)

        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputCoordinates,
                             "There was a problem with localized "
                             "subparticles protocol")
        # 1000 self.assertEqual(checkSize, prot.outputCoordinates.getSize())
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
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': fileName
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import atom struct')
        self.launchProtocol(protImportPDB)
        # structure = protImportPDB.outputPdb
        return protImportPDB

    def _createMask(self, fileName):
        protMask = self.newProtocol(XmippProtCreateMask3D,
                                    featureFilePath=fileName,
                                    source=2, samplingRate=1)
        protMask.setObjLabel('feat mask')
        self.launchProtocol(protMask)
        return protMask

    def testProtLocalizedReconstruction(self):
        print("Run ProtLocalized Reconstruction")

        # Test for filter sub-particles which are aligned in the z
        localSubparticlesAligned = self._runSubparticles(
            600,
            [-177.8, 5.5, 0.5],
            symGroup=6,
            alignSubParticles=True,
            label='define subparticles aligned in z')
        localSubparticles = self._runSubparticles(
            600,
            [-1.2, 111.6, -177.4],
            symGroup=6,
            alignSubParticles=False,
            label='define subparticles not aligned')

        # Test for filter sub-particles which are aligned in the z
        localUniqueAligend = self._runFilterSubParticles(
            120, [149.0, 64.9, 73.2],
            localSubparticlesAligned, unique=5, label='filter aligned unique=5')
        self._runFilterSubParticles(
            90, [149.0, 64.9, 73.2],
            localSubparticlesAligned, unique=5, mindist=6,
            label='filter aligned unique=5 mindist=6')
        self._runFilterSubParticles(
            50, [106.2, 112.6, -177.6],
            localSubparticlesAligned, unique=5, side=25,
            label='filter aligned unique=5 side=25')
        self._runFilterSubParticles(
            21, [47.6, 175.3, 159.1],
            localSubparticlesAligned, unique=5, top=50,
            label='filter aligned unique=5 top=50')

        # Test for filter sub-particles which are not aligned
        localUnique = self._runFilterSubParticles(
            120, [31.6, 60.155541, -138.80597],
            localSubparticles, unique=5, label='filter not aligned unique=5')
        self._runFilterSubParticles(
            90, [31.6, 60.155541, -138.80597],
            localSubparticles, unique=5, mindist=6,
            label='filter not aligned unique=5 mindist=6')
        self._runFilterSubParticles(
            50, [103.2, 66.1, -66.9],
            localSubparticles, unique=5, side=25,
            label='filter not aligned unique=5 side=25')
        self._runFilterSubParticles(
            21, [175.2, 66.1, -66.9],
            localSubparticles, unique=5, top=50,
            label='filter not aligned unique=5 top=50')

        # Test extract particles
        localExtractionAligned = self.newProtocol(
            ProtLocalizedExtraction, boxSize=26)
        localExtractionAligned.inputParticles.set(
            self.protImport.outputParticles)
        localExtractionAligned.inputCoordinates.set(
            localUniqueAligend.outputCoordinates)
        self.launchProtocol(localExtractionAligned)
        self.assertIsNotNone(localExtractionAligned.outputParticles,
                             "There was a problem with localized "
                             "extraction protocol")
        # reconstruct 3d map so we can visually check that everything is ok
        # do not check here since relion may be not installed
        self.runReconstructRelion(
            inputParticles=localExtractionAligned.outputParticles,
            symmetryGroup='C1',
            numberOfMpis=4,
            objLabel='Fourier reconstruction')

        localExtraction = self.newProtocol(ProtLocalizedExtraction, boxSize=26)
        localExtraction.inputParticles.set(self.protImport.outputParticles)
        localExtraction.inputCoordinates.set(localUnique.outputCoordinates)
        self.launchProtocol(localExtraction)
        self.assertIsNotNone(localExtraction.outputParticles,
                             "There was a problem with localized "
                             "extraction protocol")
        # reconstruct 3d map so we can visually check that everything is ok
        # do not check here since relion may be not installed
        self.runReconstructRelion(
            inputParticles=localExtraction.outputParticles,
            symmetryGroup='C1',
            numberOfMpis=4,
            objLabel='Fourier reconstruction')

    def testSetOrigin(self):
        # create set of coordinates and localize protocol
        # C8
        symOrder = 8
        # generate projections
        protImportVol, self.protImport = self.generate(
            SCIPION_SYM_NAME[XMIPP_TO_SCIPION[XMIPP_CYCLIC]][:1]+str(symOrder),
            XMIPP_SYM_NAME[XMIPP_CYCLIC][:1]+str(symOrder))
        localSubparticles = self._runSubparticles(
            92*symOrder,  # #particles*simmetry order
            [-1.2, 111.6, -177.4],
            symGroup=0,  # cyclic
            symmetryOrder=symOrder,
            vector='0.,1.,0.',
            length=60,
            alignSubParticles=False,
            label='define subparticles not aligned')
        localExtraction = self.newProtocol(
            ProtLocalizedExtraction, boxSize=26,
            inputParticles=self.protImport.outputParticles,
            inputCoordinates=localSubparticles.outputCoordinates)
        _ = self.launchProtocol(localExtraction)
        relionRecons = self.runReconstructRelion(
            inputParticles=localExtraction.outputParticles,
            symmetryGroup=f'C{symOrder}',
            numberOfMpis=4,
            objLabel='Fourier reconstruction')
        # create pdb file and import it.
        # it should be inside one of the large spheres
        # both in the large 3D map and in the small after
        # executing LocalOrigSampling
        PDBString = """ATOM      1  C   PRO B   1       0.000  60.000   0.000  0.00  0.00           C
TER       1      PRO B   1
END
"""
        f = NamedTemporaryFile(delete=False, suffix=".pdb")
        f.write(PDBString.encode('utf8'))
        f.close()
        protImportPDB = self._importAtomStruct(f.name)
        # # volume
        # featDict = self.generate_cyclic(order=8, offset=0.)
        # filename = self.createFeatFile(
        #     featDict, dim=3, bgDensity=0., scale=60.)
        # volume = self._createMask(filename)
        # self.chimeraOperate(volume.outputMask, protImportPDB.outputPdb,
        #                     objLabel='Chimera show pdb in volume')

        # call protocol
        # WARNING: the 3D map and the pdb are not related
        # Therefore they will not fit
        print("connect to process, file=sys.stderr)")
        # time.sleep(30)
        print("time done, file=sys.stderr)")
        protSetOrig = self.newProtocol(ProtLocalOrigSampling,
                                       inVolume=relionRecons.outputVolume,
                                       inputProtocol=localSubparticles,
                                       setSampling=True,
                                       samplingRate=1,
                                       atomStruct=protImportPDB.outputPdb)
        protSetOrig.setObjLabel('set orig')
        self.launchProtocol(protSetOrig)
        outVol = protSetOrig.outputVolume
        self.assertIsNotNone(outVol)
        origin = outVol.getOrigin()
        self.assertIsNotNone(outVol)
        # v = np.array([0,1,0])  # unit vector
        # v * 60 # length
        # v is the coordinate
        # sampling 1A/px
        # center volume = 13
        # v = [-45.4434 ,   0.81696,  23.2323 ]
        # v = [-95.4434 , -49.18304,  26.7677 ]  # in volume coordinates
        x = -13
        y = 60 - 13
        z = -13
        new_orig = origin.getShifts()
        self.assertAlmostEqual(x, new_orig[0], places=1)
        self.assertAlmostEqual(y, new_orig[1], places=1)
        self.assertAlmostEqual(z, new_orig[2], places=1)

    def createProjection(self, volFilename, min_tilt=0, max_tilt=0, label=""):
        # create projection
        # unfortunately there is no available protocol
        # so we will run xmipp_angular_project_library
        # -i kk.mrc  -o kk.mrcs  --sampling_rate 15 --method real_space
        # remove extension from volFilename
        volFilenameNoExt = volFilename[:-4]  # remove .mrc
        args = " -i %s" % volFilenameNoExt + ".mrc"
        args += " -o %s" % volFilenameNoExt + ".mrcs"
        args += " --sampling_rate 20"
        if min_tilt != max_tilt:
            args += " --min_tilt_angle  %d" % min_tilt
            args += " --max_tilt_angle %d" % max_tilt
        args += " --method real_space"
        progname = "xmipp_angular_project_library"
        self.xmipp3.Plugin.runXmippProgram(progname, args)

        # rename doc file to xmp since import assumes the
        # xmipp metadata files extension is xmd.
        os.rename(volFilenameNoExt + ".doc", volFilenameNoExt + ".xmd")
        # import projection
        protImportPart = self.runImportParticles(
            volFilenameNoExt + ".xmd",
            samplingRate=1.0,
            objLabel=f"import geometrical phantom projections {label}")
        protImportVol = self.runImportVolumes(
            volFilenameNoExt + ".mrc",
            1.0,
            objLabel=f"import volume from geometrical phantom {label}")
        return protImportVol, protImportPart

    def generate_ico(self, sym):
        print("SYM=", sym)
        icosahedron = Icosahedron(orientation=sym)  # 'i222r'
        featList = []
        #  sphere in virus center, helps centering
        feat = {'type': 'sph', 'operation': '+', 'density': 1.,
                'center': '0 0 0', 'specific': '15.'}
        featList.append(feat.copy())
        # 5fold points
        for _5fold in icosahedron.getVertices():
            _5fold = _5fold * 60  # scale to radius 60
            feat['center'] = '%.3f %.3f %.3f' % (_5fold[0],
                                                 _5fold[1],
                                                 _5fold[2])
            feat['specific'] = '6.'
            featList.append(feat.copy())
            for _5fold_copy in icosahedron.getVertices():
                # dot product of vectors to find adjacent vertices
                _5fold_copy = _5fold_copy * 60  # scale to radius 60
                dot_product = np.dot(_5fold, _5fold_copy)
                norm_5fold = np.linalg.norm(_5fold)
                norm_5fold_copy = np.linalg.norm(_5fold_copy)
                
                cos_angle = dot_product / (norm_5fold * norm_5fold_copy)
                if cos_angle > 1.0:  # avoid runding errors
                    cos_angle = 1.0
                elif cos_angle < -1.0:
                    cos_angle = -1.0
                angle_between = math.degrees(math.acos(cos_angle))
                if angle_between < 1e-3:
                    continue  # skip the same vertex
                # angle between vectors from the center to two adjacent
                #  vertices is about 63.43 degrees, but to be safe we use 65
                if angle_between > 65.0:
                    continue  # skip non-adjacent vertices
                # place a smaller sphere 15A away from the 5fold axis towards
                # each adjacent 5fold axis
                direction = (_5fold_copy - _5fold)
                direction = np.array(direction / np.linalg.norm(direction))
                x = _5fold[0] - direction[0] * 11
                y = _5fold[1] - direction[1] * 11
                z = _5fold[2] - direction[2] * 11
                feat['center'] = '%.3f %.3f %.3f' % (x, y, z)
                feat['specific'] = '3.'
                featList.append(feat.copy())
        # # 3fold points
        # for _3fold in icosahedron.get3foldAxis():
        #     _3fold = np.array(_3fold) * 60  # scale to radius 60
        #     feat['center'] = '%.3f %.3f %.3f' % (_3fold[0],
        #                                          _3fold[1],
        #                                          _3fold[2])
        #     feat['specific'] = '10.'
        #     featList.append(feat.copy())
        # # 2fold points
        # for _2fold in icosahedron.get2foldAxis():
        #     _2fold = np.array(_2fold) * 60  # scale to radius 60
        #     feat['center'] = '%.3f %.3f %.3f' % (_2fold[0],
        #                                          _2fold[1],
        #                                          _2fold[2])
        #     feat['specific'] = '7.'
        #     featList.append(feat.copy())
        return featList

    def generate_cyclic(self, order, offset):
        featList = []
        feat = {'type': 'cyl', 'operation': '+', 'density': 1.,
                'center': '0 0 +27', 'specific': '60.0 60.0 12 0 0 0'}
        featList.append(feat.copy())
        feat['density'] = -1.
        feat['specific'] = '48 48 12 0 0 0'
        featList.append(feat.copy())

        feat = {'type': 'sph', 'operation': '+', 'density': 1.,
                'center': None, 'specific': None}
        z_value = [-27., 0]
        for z in z_value:
            x = 0.
            y = 0.
            feat['center'] = f"{x:.3f} {y:.3f} {z:.3f}"
            feat['specific'] = '9.0'
            featList.append(feat.copy())
            for point in range(order):
                x = 60 * math.cos(2*point*math.pi/order + offset)
                y = 60 * math.sin(2*point*math.pi/order + offset)
                feat['center'] = f"{x:.3f} {y:.3f} {z:.3f}"
                feat['specific'] = '9.0'
                featList.append(feat.copy())
            for point in range(order):
                x = 30 * math.cos(2*point*math.pi/order + offset)
                y = 30 * math.sin(2*point*math.pi/order + offset)
                feat['center'] = f"{x:.3f} {y:.3f} {z:.3f}"
                feat['specific'] = '6.'
                featList.append(feat.copy())
            for point in range(order):
                x = 15 * math.cos(2*point*math.pi/order + offset)
                y = 15 * math.sin(2*point*math.pi/order + offset)
                feat['center'] = f"{x:.3f} {y:.3f} {z:.3f}"
                feat['specific'] = '3.0'
                featList.append(feat.copy())
        return featList

    def generate_helical(self, sym, order, riseValue, twist,
                         nelements, Radius, radius):
        """ Generate helical feature file
        4 parameters:
        - rise: in A
        - twist: in degrees
        - nturns: number of unitcell repetitions
        - n_strands: number of parallel strands (i.e. 2 in DNA)
        - Radius, helix radius
        - sphere radius (unit cells)
        """
        featList = []
        feat = {'type': 'sph', 'operation': '+', 'density': 1.,
                'center': None, 'specific': None}
        for strand in range(order):
            for element in range(-nelements//2, nelements//2+1):
                theta = math.radians(element * twist)
                x = math.cos(theta + 2 * math.pi * strand/order) * Radius
                y = math.sin(theta + 2 * math.pi * strand/order) * Radius
                z = element * riseValue
                feat['center'] = f"{x:.3f} {y:.3f} {z:.3f}"
                feat['specific'] = f'{radius:.3f}'
                featList.append(feat.copy())
        return featList

    def generate(self, sym='I2n5', label=None,
                 offset=0., riseValue=0, twist=0, extraSym='c1'):

        global _cache
        key_sym = sym+extraSym
        if key_sym in _cache:
            return _cache[key_sym]
        max_tilt = 0  # for helicalsymmetry
        min_tilt = 0
        offset = math.radians(offset)
        symPreffix = sym[:1]
        symSuffix = sym[1:]
        if label is None:
            label = f'generate sym {sym}'
        if symPreffix == 'I':
            featDict = self.generate_ico(symSuffix)
            filename = self.createFeatFile(
                featDict, dim=180, bgDensity=0., scale=1.)
        elif symPreffix == 'C':
            featDict = self.generate_cyclic(int(symSuffix), offset)
            filename = self.createFeatFile(
                featDict, dim=180, bgDensity=0., scale=1.)
            min_tilt = 0
            max_tilt = 0
        # elif symPreffix == 'D':
        #     pass  # generate_dihedral(int(symSuffix), offset,  mode, f)
        # elif symPreffix == 'T':
        #     pass  # generate_tetrahedral(mode, f)
        # elif symPreffix == 'O':
        #     pass  # generate_octahedral(mode, f)
        elif symPreffix == 'H':
            featDict = self.generate_helical(
                sym=extraSym[0],
                order=int(extraSym[1:]),
                riseValue=riseValue,
                twist=twist,
                nelements=40,
                Radius=100,
                radius=5)
            filename = self.createFeatFile(
                featDict, dim=392, bgDensity=0., scale=1.)
        # f.close()
        command = "xmipp_phantom_create "
        args = " -i %s" % filename
        _, outputFile1 = mkstemp(suffix=".mrc")
        args += " -o %s" % outputFile1
        runJob(None, command, args, env=Plugin.getEnviron())
        protImportVol, protImportPart = self.createProjection(
            outputFile1, min_tilt=min_tilt, max_tilt=max_tilt, label=label)
        _cache[key_sym] = (protImportVol, protImportPart)
        return protImportVol, protImportPart

    def testExtractCoordinatesIco(self):
        # I2
        # generate projections
        vector = '-42.544,0.000,51.039;\
-37.044,-8.899,54.438;\
-28.145,-5.500,59.938;\
-28.145,5.500,59.938;\
-37.044,8.899,54.438'

        protImportVol, protImport = self.generate(
            SCIPION_SYM_NAME[XMIPP_TO_SCIPION[XMIPP_I222r]],
            XMIPP_SYM_NAME[XMIPP_I222r])
        localSubparticles = self._runSubparticles(
            27600,  # 92 * 8 * 2
            [0, 0, 0],
            defVector=1,
            vector=vector,
            length=-1,
            alignSubParticles=False,
            symGroup=ProtLocalizedRecons.map_sym[SYM_I222r],
            label='define ico particles 5 vectors, 5 vectors',
            inputParticles=protImport.outputParticles

            )
        # extract subparticles
        localExtraction = self.newProtocol(
            ProtLocalizedExtraction, boxSize=26,
            inputParticles=protImport.outputParticles,
            inputCoordinates=localSubparticles.outputCoordinates)
        _ = self.launchProtocol(localExtraction)
        relionRecons = self.runReconstructRelion(
            inputParticles=localExtraction.outputParticles,
            symmetryGroup='I2',
            numberOfMpis=4,
            objLabel='Fourier reconstruction')

    def testExtractCoordinatesHelicalC1(self):
        # Radius = 100
        # radius = 5
        riseValue = 5  # raise is a reserver word in python
        twist = 18

        # H1
        symOrder = 1
        # generate projections
        protImportVol, protImport = self.generate(
            SCIPION_SYM_NAME[SYM_HELICAL], "helical-C1", offset=0.,
            riseValue=riseValue, twist=twist)
        # execute define subparticles
        localSubparticles = self._runSubparticles(
            3588,  # 1 + 92 * dim // (riseValue*2)
            [0, 0, 0],
            defVector=1,
            vector='1.000, 0.000, 0.000',
            length=100,
            alignSubParticles=False,
            symGroup=ProtLocalizedRecons.map_sym[SYM_HELICAL],
            symmetryOrder=symOrder,
            label='define helical subparticles C1',
            percentage=50,
            riseValue=riseValue,
            twist=twist,
            inputParticles=protImport.outputParticles)
        # extract subparticles
        localExtraction = self.newProtocol(
            ProtLocalizedExtraction, boxSize=26,
            inputParticles=protImport.outputParticles,
            inputCoordinates=localSubparticles.outputCoordinates)
        _ = self.launchProtocol(localExtraction)
        relionRecons = self.runReconstructRelion(
            inputParticles=localExtraction.outputParticles,
            symmetryGroup=f'C{symOrder}',
            numberOfMpis=4,
            objLabel='Fourier reconstruction')

    def testExtractCoordinatesHelicalC4(self):
        # Radius = 100
        # radius = 5
        riseValue = 5  # raise is a reserver word in python
        twist = 18

        # H4
        symOrder = 4
        # generate projections
        protImportVol, protImport = self.generate(
            SCIPION_SYM_NAME[SYM_HELICAL], "helical-C4", offset=0.,
            riseValue=riseValue, twist=twist, extraSym='C4')
        # execute define subparticles
        localSubparticles = self._runSubparticles(
            14352,  # symOrder * (1 + 92 * dim // (riseValue*2))
            [0, 0, 0],
            defVector=1,
            vector='1.000, 0.000, 0.000',
            length=100,
            alignSubParticles=False,
            symGroup=ProtLocalizedRecons.map_sym[SYM_HELICAL],
            symmetryOrder=symOrder,
            extraSym=ProtLocalizedRecons.map_sym[SYM_CYCLIC],
            label='define helical subparticles C4',
            percentage=50,
            riseValue=riseValue,
            twist=twist,
            inputParticles=protImport.outputParticles
            )
        # extract subparticles
        localExtraction = self.newProtocol(
            ProtLocalizedExtraction, boxSize=26,
            inputParticles=protImport.outputParticles,
            inputCoordinates=localSubparticles.outputCoordinates)
        _ = self.launchProtocol(localExtraction)
        relionRecons = self.runReconstructRelion(
            inputParticles=localExtraction.outputParticles,
            symmetryGroup=f'C{symOrder}',
            numberOfMpis=4,
            objLabel='Fourier reconstruction')
