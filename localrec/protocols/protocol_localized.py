# *****************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *              Vahid Abrishami (vahid.abrishami@helsinki.fi)
# *              Roberto Marabini (roberto@cnb.csic.es)
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
# *****************************************************************************
from pyworkflow.protocol.params import (
    PointerParam, BooleanParam, StringParam,
    EnumParam, IntParam,
    PathParam, FloatParam)
from pwem.constants import (
    SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL,
    SYM_TETRAHEDRAL, SYM_I222, SYM_I222r,
    SYM_In25, SYM_In25r, SCIPION_SYM_NAME,
    SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r, SYM_HELICAL)
from pwem.convert.symmetry import getSymmetryMatrices
from pwem.protocols import ProtParticles, ProtParticlePicking
import pwem.emlib.metadata as md

from localrec.utils import load_vectors, create_subparticles
from localrec.constants import CMM
from pyworkflow.utils import ProgressBar
from pwem.objects import SetOfCoordinates
# import numpy as np


class ProtLocalizedRecons(ProtParticlePicking, ProtParticles):
    """ This class contains a re-implementation to a method for the
    localized three-dimensional reconstruction of such subunits.
    After determining the particle orientations, local areas
    corresponding to the subunits can be extracted and treated as
    single particles.
    """
    _label = 'define subparticles'
    OUTPUTCOORDINATESNAME = "outputCoordinates"
    _possibleOutputs = {OUTPUTCOORDINATESNAME: SetOfCoordinates}
    # dictionary that convert locally used simmetries to scipion standard
    sym_map = {
        0: SYM_CYCLIC,
        1: SYM_DIHEDRAL,
        2: SYM_TETRAHEDRAL,
        3: SYM_OCTAHEDRAL,
        4: SYM_I222,
        5: SYM_I222r,
        6: SYM_In25,
        7: SYM_In25r,
        8: SYM_I2n3,
        9: SYM_I2n3r,
        10: SYM_I2n5,
        11: SYM_I2n5r,
        12: SYM_HELICAL
    }
    # dictionary to convert scipion standard symmetries to local ones
    map_sym = {v: k for k, v in sym_map.items()}

    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)
        ProtParticles.__init__(self, **args)
        self.allowMpi = False
        self.allowThreads = False

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label="Input particles",
                      help='Select the input images from the project.')

        group = form.addGroup('Symmetry')
        group.addParam(
            'symGrp',
            EnumParam,
            choices=["Cn" + " (" + SCIPION_SYM_NAME[SYM_CYCLIC] + ")",   # 0
                     "Dn" + " (" + SCIPION_SYM_NAME[SYM_DIHEDRAL] + ")",  # 1
                     "T" + " (" + SCIPION_SYM_NAME[SYM_TETRAHEDRAL] + ")",  # 2
                     "O" + " (" + SCIPION_SYM_NAME[SYM_OCTAHEDRAL] + ")",  # 3
                     "I1" + " (" + SCIPION_SYM_NAME[SYM_I222] + ")",   # 4
                     "I2" + " (" + SCIPION_SYM_NAME[SYM_I222r] + ")",  # 5
                     "I3" + " (" + SCIPION_SYM_NAME[SYM_In25] + ")",  # 6
                     "I4" + " (" + SCIPION_SYM_NAME[SYM_In25r] + ")",  # 7
                     "I5" + " (" + SCIPION_SYM_NAME[SYM_I2n3] + ")",  # 8
                     "I6" + " (" + SCIPION_SYM_NAME[SYM_I2n3r] + ")",  # 9
                     "I7" + " (" + SCIPION_SYM_NAME[SYM_I2n5] + ")",  # 10
                     "I8" + " (" + SCIPION_SYM_NAME[SYM_I2n5r] + ")",  # 11
                     "H" + " (" + SCIPION_SYM_NAME[SYM_HELICAL] + ")"  # 12
                     ],
            default=0,
            label="Symmetry",
            help="See https://scipion-em.github.io/docs/docs/developer/symmetries"
                  "Symmetry for a description of the symmetry groups "
                  "format in Xmipp.\n"
                  "If no symmetry is present, use _c1_."
        )
        group.addParam('rise', FloatParam, default=0.0,
                       condition='symGrp==%d' % self.map_sym[SYM_HELICAL],
                       label='Raise (A)',
                       help='The raise is the linear distance along the helix'
                            ' axis between consecutive repeating units.')
        group.addParam('twist', FloatParam, default=0.0,
                       condition='symGrp==%d' % self.map_sym[SYM_HELICAL],
                       label='Twist (degrees)',
                       help='The twist is the angular rotation between'
                            ' consecutive units around the helix axis.')
        group.addParam('percentage', FloatParam, default=100,
                       condition='symGrp==%d' % self.map_sym[SYM_HELICAL],
                       label='Percentage (%)',
                       help='replicate coordinates only for subparticles with'
                            ' z coordinate greater than -dim*percentage/2.'
                            ' and smaller than dim*percentage/2. '
                            ' where dim is the projection image dimensions'
                            ' If percentage = 100 dim/rise coordinates are'
                            'created.'
                       )
        group.addParam(
            'symGrpExtra',
            EnumParam,
            choices=["Cn" + " (" + SCIPION_SYM_NAME[SYM_CYCLIC] + ")",   # 0
                     "Dn" + " (" + SCIPION_SYM_NAME[SYM_DIHEDRAL] + ")"  # 1
                     ],
            condition='symGrp==%d' % self.map_sym[SYM_HELICAL],
            default=0,
            label="Extra Symmetry (helical only)",
            help="See https://scipion-em.github.io/docs/docs/developer/symmetries"
                  "Symmetry for a description of the symmetry groups "
                  "format in Xmipp.\n"
                  "If no symmetry is present, use _c1_."
        )
        group.addParam('symmetryOrder', IntParam, default=1,
                       label='Symmetry Order',
                       help='Order of cyclic symmetry.')
        group.addParam('randomize', BooleanParam, default=False,
                       label='Randomize the order of the symmetry matrices?',
                       help='Useful for preventing preferred orientations.')

        group = form.addGroup('Vectors')
        group.addParam('defineVector', EnumParam, default=CMM,
                       label='Is vector defined by?',
                       choices=['cmm file', 'string'],
                       display=EnumParam.DISPLAY_HLIST)
        group.addParam('vector', StringParam, default='0,0,1',
                       label='Location vectors (px)', condition="defineVector==1",
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

        group = form.addGroup('Sub-particles')
        group.addParam('alignSubParticles', BooleanParam, default=False,
                       label='Align the subparticles?',
                       help='Align sub-particles to the standard orientation. ')
        group.addParam('handness', BooleanParam, default=False,
                       label='Consider alternative handedness?',
                       help='If set to yes, the alternative hand is assumed'
                            ' to be correct and defocus gradient correction'
                            ' will be applied in the opposite direction . ')

        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -----------------------------
    def createOutputStep(self):
        # return
        inputSet = self.inputParticles.get()
        outputSet = self._createSetOfCoordinates(inputSet)
        params = {"vector": self.vector.get(),
                  "vectorFile": self.vectorFile.get(),
                  "length": self.length.get(),
                  "pxSize": self.inputParticles.get().getSamplingRate(),
                  "dim": self.inputParticles.get().getXDim()
                  }
        # convert symmetry to scipion
        sym = self.symGrp.get()
        sym = self.sym_map.get(sym, sym)  # Keeps original value if not found

        if sym == SYM_HELICAL:
            rise = self.rise.get()
            twist = self.twist.get()
            percentage = self.percentage.get() / 100.0
            dim = params["dim"]
            n = int((dim * percentage) / rise / 2)  # number of subparticles/2
            # ROB DEBUG
            symMatrices = getSymmetryMatrices(sym=sym,
                                              n=n,
                                              rise=rise,
                                              angle=twist
                                              )
            extraSym = self.symGrpExtra.get()
            symmetryOrder = self.symmetryOrder.get()
            extraSymmetryMatrices = getSymmetryMatrices(sym=extraSym,
                                                        n=symmetryOrder)
            newSymMatrices = []
            for mat1 in symMatrices:
                for mat2 in extraSymmetryMatrices:
                    newMat = mat2 @ mat1
                    # ROB DEBUG
                    print("type(mat1)=%s, type(mat2)=%s" % (type(mat1), type(mat2)))
                    print("Combining matrices:\n%s\nAND\n%s" % (mat2, mat1))
                    print("New symmetry matrix:\n%s" % newMat)
                    newSymMatrices.append(newMat)
            symMatrices = newSymMatrices
            # ROB DEBUG
            for i, mat in enumerate(symMatrices):
                print("Symmetry matrix %d:\n%s" % (i, mat))
        else:
            n = self.symmetryOrder.get()
            symMatrices = getSymmetryMatrices(sym=sym, n=n)

        if self.defineVector == CMM:
            cmmFn = params["vectorFile"]
            vector = " "
        else:
            cmmFn = ""
            vector = params["vector"]

        subpartVectorList = load_vectors(cmmFn, vector,
                                         params["length"],
                                         params["pxSize"])
        vectorsMd = md.MetaData()
        # Save the vectors in a metadata
        for vector in subpartVectorList:
            objId = vectorsMd.addObject()
            vectorsMd.setValue(md.MDL_SHIFT_X, vector.vector[0], objId)
            vectorsMd.setValue(md.MDL_SHIFT_Y, vector.vector[1], objId)
            vectorsMd.setValue(md.MDL_SHIFT_Z, vector.vector[2], objId)
        vectorsMd.write(self._getOutpuVecMetadata())

        # for part in inputSet:
        print("Processing coordinates:")
        progress = ProgressBar(total=len(inputSet), fmt=ProgressBar.NOBAR)
        progress.start()

        step = max(100, len(inputSet) // 100)

        for i, part in enumerate(inputSet):
            if i % step == 0:
                progress.update(i+1)

            subparticles = create_subparticles(part, symMatrices,
                                               subpartVectorList,
                                               params["dim"],
                                               self.randomize, 0,
                                               self.alignSubParticles,
                                               self.handness,
                                               params["pxSize"])

            for subpart in subparticles:
                coord = subpart.getCoordinate()
                outputSet.append(coord)
                if part.hasAttribute('_rlnRandomSubset'):
                    coord._subparticle.copyAttributes(part, '_rlnRandomSubset')
        progress.finish()

        self._defineOutputs(**{self.OUTPUTCOORDINATESNAME: outputSet})
        self._defineSourceRelation(self.inputParticles, outputSet)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        pass

    def _citations(self):
        return ['Ilca2015', 'Abrishami2020']

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        return []

    # -------------------------- UTILS functions ------------------------------
    def _getOutpuVecMetadata(self):
        return self._getExtraPath('output_vectors.xmd')
