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
from pyworkflow.protocol.params import (PointerParam, BooleanParam, StringParam,
                                        EnumParam, NumericRangeParam, IntParam,
                                        PathParam)
from pwem.constants import (SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL,
                            SYM_TETRAHEDRAL, SYM_I222, SYM_I222r,
                            SYM_In25, SYM_In25r, SCIPION_SYM_NAME,
                            SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r)
from pwem.convert.symmetry import getSymmetryMatrices
from pwem.protocols import ProtParticles, ProtParticlePicking
import pwem.emlib.metadata as md

from localrec.utils import load_vectors, create_subparticles
from localrec.constants import CMM, HAND
from pyworkflow.utils import ProgressBar


class ProtLocalizedRecons(ProtParticlePicking, ProtParticles):
    """ This class contains a re-implementation to a method for the
    localized three-dimensional reconstruction of such subunits.
    After determining the particle orientations, local areas
    corresponding to the subunits can be extracted and treated as
    single particles.
    """
    _label = 'define subparticles'

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
        group.addParam('symGrp', EnumParam,
                      choices=["Cn" + " (" + SCIPION_SYM_NAME[SYM_CYCLIC] + ")",
                               "Dn" + " (" + SCIPION_SYM_NAME[SYM_DIHEDRAL] + ")",
                               "T" + " (" + SCIPION_SYM_NAME[SYM_TETRAHEDRAL] + ")",
                               "O" + " (" + SCIPION_SYM_NAME[SYM_OCTAHEDRAL] + ")",
                               "I1" + " (" + SCIPION_SYM_NAME[SYM_I222] + ")",
                               "I2" + " (" + SCIPION_SYM_NAME[SYM_I222r] + ")",
                               "I3" + " (" + SCIPION_SYM_NAME[SYM_In25] + ")",
                               "I4" + " (" + SCIPION_SYM_NAME[SYM_In25r] + ")",
                               "I5" + " (" + SCIPION_SYM_NAME[SYM_I2n3] + ")",
                               "I6" + " (" + SCIPION_SYM_NAME[SYM_I2n3r] + ")",
                               "I7" + " (" + SCIPION_SYM_NAME[SYM_I2n5] + ")",
                               "I8" + " (" + SCIPION_SYM_NAME[SYM_I2n5r] + ")"],
                      default=0,
                      label="Symmetry",
                      help="See https://scipion-em.github.io/docs/docs/developer/symmetries"
                           "Symmetry for a description of the symmetry groups "
                           "format in Xmipp.\n"
                           "If no symmetry is present, use _c1_."
                       )
        group.addParam('symmetryOrder', IntParam, default=1,
                      condition='symGrp<=%d' % SYM_DIHEDRAL,
                      label='Symmetry Order',
                      help='Order of cyclic symmetry.')
        group.addParam('randomize', BooleanParam, default=False,
                       label='Randomize the order of the symmetry matrices?',
                       help='Useful for preventing preferred orientations.')
# ROB this parameter is never used
#        group.addParam('relaxSym', BooleanParam, default=False,
#                       expertLevel=LEVEL_ADVANCED,
#                       label='Relax symmetry?',
#                       help='Create one random subparticle for each particle ')

        group = form.addGroup('Vectors')
        group.addParam('defineVector', EnumParam, default=CMM,
                       label='Is vector defined by?',
                       choices=['cmm file', 'string'],
                       display=EnumParam.DISPLAY_HLIST)
        group.addParam('vector', NumericRangeParam, default='0,0,1',
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
        # symDict = {0: 'C', 1: 'D', 2: 'T', 3: 'O',
        # 4: 'I1', 5: 'I2', 6: 'I3', 7: 'I4'}
        if sym == 0: sym = SYM_CYCLIC
        elif sym == 1: sym = SYM_DIHEDRAL
        elif sym == 2: sym = SYM_TETRAHEDRAL
        elif sym == 3: sym = SYM_OCTAHEDRAL
        elif sym == 4: sym = SYM_I222
        elif sym == 5: sym = SYM_I222r
        elif sym == 6: sym = SYM_In25
        elif sym == 7: sym = SYM_In25r
        elif sym == 8: sym = SYM_I2n3
        elif sym == 9: sym = SYM_I2n3r
        elif sym == 10: sym = SYM_I2n5
        elif sym == 11: sym = SYM_I2n5r
        symMatrices = getSymmetryMatrices(sym=sym,
                                          n=self.symmetryOrder.get())
###ROB                                          n = self.symmetryOrder.get())
#        for mat in symMatrices:
#            print (mat)
#        exit(0)
####

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

        #for part in inputSet:
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

        self._defineOutputs(outputCoordinates=outputSet)
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

