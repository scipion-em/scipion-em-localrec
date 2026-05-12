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
from localrec.constants import (CMM, HAND, symDictScipion,
                                symDictScipionFromKey, symDictFromKey)
from pyworkflow.utils import ProgressBar
from pwem.objects import SetOfCoordinates


class ProtLocalizedRecons(ProtParticlePicking, ProtParticles):
    """
    Defines localized subparticles from previously aligned particles in
    order to enable focused reconstruction of specific structural regions
    within a larger macromolecular assembly.

    AI Generated:

    Define Subparticles (ProtLocalizedRecons) - User Manual
        Overview

        The Define Subparticles protocol generates localized subparticles
        from a set of particles with known projection orientations. Its
        main purpose is to isolate biologically relevant regions of a
        macromolecular complex so they can be analyzed independently from
        the full particle. This strategy is especially important in cryo-EM
        studies involving large symmetric assemblies, flexible domains,
        membrane complexes, or structures containing repeated subunits.

        In many biological systems, local regions may exhibit structural
        variability or conformational heterogeneity that becomes blurred
        during conventional whole-particle refinement. By extracting
        localized regions and treating them as independent particles,
        researchers can improve the interpretability and resolution of
        flexible or partially occupied components.

        Inputs and Biological Context

        The protocol requires a set of particles that already contain
        projection alignment information. These orientations define the
        spatial relationship between each experimental image and the
        reconstructed structure, allowing the protocol to determine where
        the desired subunits are located within every particle.

        From a biological perspective, this workflow is commonly applied
        to viral capsids, ribonucleoprotein assemblies, membrane protein
        complexes, and highly symmetric particles where repeated domains
        or subunits can be analyzed separately. It is particularly useful
        when the region of interest represents only a fraction of the
        total particle volume.

        Symmetry Definition

        Symmetry plays a central role in localized reconstruction because
        it determines how equivalent subunits are distributed throughout
        the structure. The protocol supports cyclic, dihedral,
        tetrahedral, octahedral, and multiple icosahedral symmetry
        conventions commonly used in cryo-EM processing.

        Correct symmetry selection is biologically essential. Using the
        proper symmetry ensures that all equivalent subunits are sampled
        consistently and that localized regions are extracted in the
        correct spatial orientations. Incorrect symmetry assignment may
        generate biologically meaningless subparticles or introduce
        orientation inconsistencies.

        The protocol also allows the selection of specific symmetry
        operators instead of using the complete symmetry set. This can be
        valuable when the biological question focuses on only a subset of
        equivalent positions or when users wish to reduce redundancy in
        downstream analyses.

        Definition of Localization Vectors

        Subparticle positions are defined through vectors describing the
        location of the region of interest relative to the particle
        center. These vectors may be provided manually or imported from
        external visualization tools such as Chimera marker files.

        In practical cryo-EM workflows, the vectors usually correspond to
        biologically meaningful structural features such as protruding
        domains, symmetry-related subunits, membrane anchors, flexible
        arms, or ligand-binding regions. Accurate vector placement is one
        of the most important factors for successful localized
        reconstruction.

        The protocol also allows adjustment of vector length. This is
        particularly useful when refining the exact center of the
        subparticle or compensating for uncertainties in the initial
        structural interpretation. Small changes in vector placement can
        substantially affect the quality of the extracted signal.

        Subparticle Alignment and Orientation

        An optional alignment mode can orient subparticles into a common
        reference frame. This facilitates downstream classification,
        averaging, and local refinement by reducing orientation
        variability between equivalent subunits.

        From a biological standpoint, aligned subparticles are often
        easier to interpret visually and can improve the detection of
        subtle conformational changes. However, users should remember
        that aggressive alignment strategies may also reduce observable
        structural diversity if the target region is highly flexible.

        The protocol additionally supports alternative handedness
        interpretation. This option is mainly relevant in workflows where
        handedness ambiguity exists or when defocus-gradient behavior
        requires inversion relative to the original orientation
        convention.

        Randomization of Symmetry Operators

        The protocol can randomize the order of symmetry operators during
        subparticle generation. This feature is particularly important in
        cases where preferred orientations or uneven sampling may bias
        downstream classification and refinement.

        In biological datasets with strong orientation preferences,
        randomization can improve angular diversity and reduce systematic
        overrepresentation of particular views. This often leads to more
        balanced local reconstructions and improved interpretability of
        flexible regions.

        Outputs and Their Interpretation

        The protocol produces a set of coordinates representing the
        localized subparticles derived from the original particle set.
        Each coordinate preserves the relationship to its parent particle
        and contains the information necessary for subsequent extraction
        and local reconstruction workflows.

        Biologically, these outputs define a new experimental dataset
        focused specifically on the selected structural regions. The
        resulting subparticles can then be subjected to classification,
        refinement, variability analysis, or focused reconstruction using
        standard cryo-EM processing tools.

        Practical Recommendations

        In most biological applications, it is advisable to begin with
        carefully validated particle orientations and an accurately
        interpreted structural model before defining localization
        vectors. Errors in orientation assignment or vector placement are
        among the most common causes of poor localized reconstruction
        results.

        For symmetric assemblies, users should verify that the selected
        symmetry convention matches the biological architecture of the
        specimen. When uncertainty exists, exploratory tests with limited
        symmetry operators can help evaluate the consistency of the
        extracted subparticles.

        Vector placement should generally target structurally stable and
        biologically meaningful regions. Extremely flexible domains may
        still benefit from localized reconstruction, but additional
        classification steps are often required to separate distinct
        conformational states.

        Final Perspective

        Localized reconstruction is one of the most powerful approaches
        for studying structural heterogeneity in large cryo-EM complexes.
        By isolating individual subunits or regions from a larger
        particle, the protocol enables focused structural analysis that
        would otherwise be obscured by symmetry averaging or conformational
        variability. Careful definition of symmetry, vector placement,
        and subparticle orientation is essential for obtaining reliable
        and biologically meaningful results.
    """
    _label = 'define subparticles'
    OUTPUTCOORDINATESNAME = "outputCoordinates"
    _possibleOutputs = {OUTPUTCOORDINATESNAME: SetOfCoordinates}
    sym_choices = [f"{k}, {SCIPION_SYM_NAME[v]}" 
                   for k, v in symDictScipionFromKey.items()]

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
                       choices=self.sym_choices,
                       #   choices=["Cn" + " (" + SCIPION_SYM_NAME[SYM_CYCLIC] + ")",
                       #    "Dn" + " (" + SCIPION_SYM_NAME[SYM_DIHEDRAL] + ")",
                       #    "T" + " (" + SCIPION_SYM_NAME[SYM_TETRAHEDRAL] + ")",
                       #    "O" + " (" + SCIPION_SYM_NAME[SYM_OCTAHEDRAL] + ")",
                       #    "I1" + " (" + SCIPION_SYM_NAME[SYM_I222] + ")",
                       #    "I2" + " (" + SCIPION_SYM_NAME[SYM_I222r] + ")",
                       #    "I3" + " (" + SCIPION_SYM_NAME[SYM_In25] + ")",
                       #    "I4" + " (" + SCIPION_SYM_NAME[SYM_In25r] + ")",
                       #    "I5" + " (" + SCIPION_SYM_NAME[SYM_I2n3] + ")",
                       #    "I6" + " (" + SCIPION_SYM_NAME[SYM_I2n3r] + ")",
                       #    "I7" + " (" + SCIPION_SYM_NAME[SYM_I2n5] + ")",
                       #    "I8" + " (" + SCIPION_SYM_NAME[SYM_I2n5r] + ")"],
                       default=0,
                       label="Symmetry",
                       help="See https://scipion-em.github.io/docs/docs/developer/symmetries"
                            "Symmetry for a description of the symmetry groups "
                            "format in Xmipp.\n"
                            "If no symmetry is present, use _c1_.")

        group.addParam('symmetryOrder', IntParam, default=1,
                       condition='symGrp<=%d' % symDictFromKey['D'],
                       label='Symmetry Order',
                       help='Order of cyclic symmetry.')

        group.addParam('symmetryMatrixIds', StringParam, default='',
                       label='Symmetry matrix IDs',
                       help='Optional comma-separated list of symmetry matrix '
                            'IDs to use (for example: 1,3,5). IDs are '
                            '1-based and must be between 1 and N, where N is '
                            'the total number of matrices for the selected '
                            'symmetry. If empty, all matrices are used.')

        group.addParam('randomize', BooleanParam, default=False,
                       label='Randomize the order of the symmetry matrices?',
                       help='Useful for preventing preferred orientations.')

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
        # symGrpParam = sym
        sym = symDictScipion[sym]  # now sym follows scipion convention
        # symDict = {0: 'C', 1: 'D', 2: 'T', 3: 'O',
        # 4: 'I1', 5: 'I2', 6: 'I3', 7: 'I4'}
        # if sym == 0: sym = SYM_CYCLIC
        # elif sym == 1: sym = SYM_DIHEDRAL
        # elif sym == 2: sym = SYM_TETRAHEDRAL
        # elif sym == 3: sym = SYM_OCTAHEDRAL
        # elif sym == 4: sym = SYM_I222
        # elif sym == 5: sym = SYM_I222r
        # elif sym == 6: sym = SYM_In25
        # elif sym == 7: sym = SYM_In25r
        # elif sym == 8: sym = SYM_I2n3
        # elif sym == 9: sym = SYM_I2n3r
        # elif sym == 10: sym = SYM_I2n5
        # elif sym == 11: sym = SYM_I2n5r
        allSymMatrices = getSymmetryMatrices(sym=sym,  
                                             n=self.symmetryOrder.get())
        # NOTE: selectedSymmetryMatrixIds are always 1-based IDs that refer to
        # the full symmetry set (allSymMatrices), never to a filtered position.
        selectedSymmetryMatrixIds = self._parseSymmetryMatrixIds(
            self.symmetryMatrixIds.get(), len(allSymMatrices))
        if selectedSymmetryMatrixIds:
            symMatrices = [allSymMatrices[matrixId - 1]
                           for matrixId in selectedSymmetryMatrixIds]
            operatorIds = selectedSymmetryMatrixIds
        else:
            symMatrices = allSymMatrices
            operatorIds = list(range(1, len(symMatrices) + 1))
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
                                               params["pxSize"],
                                               symmetry_operator_ids=operatorIds,
                                               symmetry_group=sym)

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
        validateMsgs = []
        sym = self.symGrp.get()
        sym = symDictScipion[sym]
        # if sym == 0:
            # sym = SYM_CYCLIC
        # elif sym == 1:
            # sym = SYM_DIHEDRAL
        # elif sym == 2:
            # sym = SYM_TETRAHEDRAL
        # elif sym == 3:
            # sym = SYM_OCTAHEDRAL
        # elif sym == 4:
            # sym = SYM_I222
        # elif sym == 5:
            # sym = SYM_I222r
        # elif sym == 6:
            # sym = SYM_In25
        # elif sym == 7:
            # sym = SYM_In25r
        # elif sym == 8:
            # sym = SYM_I2n3
        # elif sym == 9:
            # sym = SYM_I2n3r
        # elif sym == 10:
            # sym = SYM_I2n5
        # elif sym == 11:
            # sym = SYM_I2n5r
        symMatrices = getSymmetryMatrices(sym=sym, n=self.symmetryOrder.get())

        try:
            self._parseSymmetryMatrixIds(self.symmetryMatrixIds.get(),
                                         len(symMatrices))
        except ValueError as e:
            validateMsgs.append(str(e))

        return validateMsgs

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

    @staticmethod
    def _parseSymmetryMatrixIds(matrixIdsText, nSymmetryMatrices):
        matrixIdsText = matrixIdsText.strip()
        if not matrixIdsText:
            return []

        parsedIds = []
        for token in matrixIdsText.split(','):
            token = token.strip()
            if not token:
                raise ValueError('Symmetry matrix IDs contain empty values. '
                                 'Please provide a comma-separated list of '
                                 'integers.')
            if not token.isdigit():
                raise ValueError('Symmetry matrix IDs must be integers. '
                                 'Invalid value: %s' % token)

            matrixId = int(token)
            if matrixId < 1 or matrixId > nSymmetryMatrices:
                raise ValueError('Symmetry matrix IDs must be between 1 and '
                                 '%d. Invalid value: %d'
                                 % (nSymmetryMatrices, matrixId))
            if matrixId in parsedIds:
                raise ValueError('Symmetry matrix IDs contain duplicates. '
                                 'Repeated value: %d' % matrixId)

            parsedIds.append(matrixId)

        return parsedIds
