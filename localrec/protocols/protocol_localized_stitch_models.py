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
from turtle import distance
from pwem.emlib import DT_DOUBLE
from pyworkflow.protocol.params import (EnumParam, StringParam, BooleanParam,
                                        NumericRangeParam, PathParam, MultiPointerParam)


from pwem.convert.atom_struct import AtomicStructHandler
from pwem.convert.transformations import euler_from_matrix, euler_matrix
from pwem.protocols import EMProtocol
from pwem.objects.data import *
from pyworkflow.protocol.constants import STEPS_PARALLEL
from Bio import PDB
from pwem.convert.symmetry import getSymmetryMatrices, _applyMatrix
from glob import glob
from Bio.PDB.mmcifio import MMCIFIO
from pwem.constants import (SYM_CYCLIC, SYM_DIHEDRAL, SYM_OCTAHEDRAL,
                            SYM_TETRAHEDRAL, SYM_I222, SYM_I222r,
                            SYM_In25, SYM_In25r, SCIPION_SYM_NAME, SYM_TETRAHEDRAL_Z3,
                            SYM_I2n3, SYM_I2n3r, SYM_I2n5, SYM_I2n5r, SYM_DIHEDRAL_X, SYM_DIHEDRAL_Y)

from localrec.constants import CMM
from localrec.utils import euler2matrix,distances_from_string, generate_chain_id, load_vectors, pdbIds_from_string, vectors_from_cmm, vectors_from_string
I1 = 1

class ProtLocalizedStitchModels(EMProtocol):
    """
        Generate a full volume from a sub-model applying a
        point group symmetry operation.

        An example of usage is to generate the adenovirus capsid
        from its asymmetric unit.
    """

    _label = 'stitch models'

    LINEAR = 0
    BSPLINE = 1

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input Structure')
        group.addParam('definePdb', EnumParam, default=CMM,
                       label='Input Structure Came From?',
                       choices=['File', 'PDB ID'],
                       display=EnumParam.DISPLAY_HLIST)
        group.addParam('pdbId', StringParam, default='',
                       label='PDB ID', condition="definePdb==1",
                       help='see https://www.rcsb.org/'
                       )
        group.addParam('pdbFile', PathParam, default='',
                       condition="definePdb==0",
                       label='PDB File: ',
                       help="Directory with the files you want to import.\n"
                            "The path can also contain wildcards to selectfrom several folders.\n" 
                            "Examples:\n"
                            "\t~/project/data/day??_files/\n"
                            "\tEach '?' represents one unknown character\n"
                            "\t~/project/data/day*_files/\n"
                            "\t'*' represents any number of unknown characters\n"
                            "\t~/project/data/day##_files/\n"
                            "\t'##' represents two digits that will be used as file ID\n"
                            "NOTE: wildcard characters ('*', '?', '#') cannot appear in the actual path.\)\n"
                            "See more: https://docs.python.org/library/glob.html\n"
                        )
        
        form.addParam('usePreRun', BooleanParam,
                      label="Use previous localrec run(s)", default=False)

        form.addParam('preRuns', MultiPointerParam, pointerClass='ProtLocalizedRecons',
                      label='Localrec previous runs', allowsNull=True,
                      condition="usePreRun",
                      help="Previous Localrec runs used to extract the parameters")
        group = form.addGroup('Symmetry')
        group.addParam('symmetryName', EnumParam, 
                      default=0,
                      choices=[ SCIPION_SYM_NAME[SYM_CYCLIC],           # 'Cn'
                                SCIPION_SYM_NAME[SYM_DIHEDRAL],         # 'Dn'
                                SCIPION_SYM_NAME[SYM_TETRAHEDRAL],      # 'T222'
                                SCIPION_SYM_NAME[SYM_TETRAHEDRAL_Z3],   # 'Tz3'
                                SCIPION_SYM_NAME[SYM_OCTAHEDRAL],       #  'O'
                                 SCIPION_SYM_NAME[SYM_I222],            # 'I222'
                                SCIPION_SYM_NAME[SYM_I222r],            # 'I222r'
                                SCIPION_SYM_NAME[SYM_In25],             # 'In25'
                                SCIPION_SYM_NAME[SYM_In25r],            # 'In25r'
                                ],
                      label="Symmetry",
                      help="See https://scipion-em.github.io/docs/docs/developer/symmetries"
                           "Symmetry for a description of the symmetry groups "
                           "format in Xmipp.\n"
                           "If no symmetry is present, use _c1_."
                       )
        group.addParam('alignSubParticles', BooleanParam,
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
        group = form.addGroup('Output')   
        group.addParam('outputOnlyMatrices', BooleanParam,
                      label="Output symmetry matrices only",
                      default=False,
                      help='Set to Yes if you wat to get your output file with Biological Assembly.')
                       
        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):

        
        inputStepId = self._insertFunctionStep('convertInputStep')
        
        tfCoordId = self._insertFunctionStep('transformCoordinatesStep',
                                                    prerequisites=[inputStepId])

        # If user wants to get biological assembly in their outputfile, we continue with biological assembly step
        if self.outputOnlyMatrices.get():
            bioAssmblyId = self._insertFunctionStep('biologicalAssemblyStep', prerequisites=[tfCoordId])
            self._insertFunctionStep('createOutputStep', prerequisites=[bioAssmblyId])
        else:
            applySymId = self._insertFunctionStep('applySymmetryStep', prerequisites=[tfCoordId])
            self._insertFunctionStep('createOutputStep', prerequisites=[applySymId])
        
        

    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        """
         Convert the input to use later on in this protocol. 
        """
        self.symGroup = self.symmetryName.get()
        self.doAlign = self.alignSubParticles
        # As we don't use Dyn symmetry and we get symmetry matrices by their enum correspondings 
        # we need to add 1 if choosen symmetry comes after Dyn
        self.symGroup = self.symGroup + 1 if self.symGroup > 2 else self.symGroup
        # Sort input file names alphabetically to match them to the input vectors and distances
        ah = AtomicStructHandler()
        if self.definePdb == 1:
            ids = pdbIds_from_string(self.pdbId.get())
            self.inputFiles = [ah.readFromPDBDatabase(id, self._getTmpPath()) for id in ids]
        
        else:
            self.inputFiles =  glob(self.pdbFile.get())
        self.inputFiles.sort()


        # if self.defineVector == CMM:
        #     self.vectorsFile = self.vectorFile.get()
        #     #self.vectors = [vectors_from_cmm(self.vectorsFile, 1)]
        # else:
        #     self.vectors = vectors_from_string(self.vector.get())
        # self.distances = distances_from_string(self.length.get()) if self.length !='-1' else [-1]
        
        vector = ""
        cmmFn = ""
        if self.defineVector == CMM:
            cmmFn = self.vectorFile.get()
        else:
            vector = self.vector.get()

        self.subVolCenterVec = load_vectors(cmmFn, vector, self.length.get(), 1)


        # for vector, distance in zip(self.vectors, self.distances):
        #     if distance > 0:
        #         vector.set_length(distance)
        #     else:
        #         vector.compute_length()
        



    def transformCoordinatesStep(self):
        """
            Transform input structures and bind them to a master structure to make them ready for symmetrisation
        """
        
        listOfAtomicStructObjects = []
        # vectors = self.subparticle_vector_list
        # distances = self.distances
        for file in self.inputFiles:
            ah = AtomicStructHandler()
            print(file)
            ah.read(file)
            listOfAtomicStructObjects.append(ah.structure.copy())
        boxSize = 800
        samplingRate = 1.34
        # valueToShift = (boxSize/2)*samplingRate

        valueToShift = (((boxSize/2)/2)/2)*(samplingRate+0.04*3) # 100*1.46 = 146
        print(valueToShift)
        for struct in listOfAtomicStructObjects:
            rotMatrix = np.identity(3)
            struct.transform(rotMatrix, np.array([-valueToShift, -valueToShift, -valueToShift]))
        #vectorsForShift = [np.array([628.89361348, 390.23970294, 732.22162194]), np.array([459.99159856, 399.05464447, 733.41199339])]
        for i, struct in enumerate(listOfAtomicStructObjects):
            
            shiftX, shiftY, shiftZ, rotMatrixFromVector = self.readVector(i)
            rotMatrix = np.identity(3)

            if self.doAlign:
                rotMatrix = rotMatrixFromVector

            # rotation_matrix_transposed = np.transpose(rotMatrix)   
            vectorForShift = np.array([shiftX, shiftY, shiftZ])
            vectorForShift = np.array([0, 0, 0])
            struct.transform(rotMatrix, vectorForShift)

        for struct in listOfAtomicStructObjects:
            rotMatrix = np.identity(3)
            struct.transform(rotMatrix, np.array([valueToShift, valueToShift, valueToShift]))

        for i, struct in enumerate(listOfAtomicStructObjects):
            
            shiftX, shiftY, shiftZ, rotMatrixFromVector = self.readVector(i)
            rotMatrix = np.identity(3)
            vectorForShift = np.array([shiftX, shiftY, shiftZ])
            struct.transform(rotMatrix, vectorForShift)
            
        masterStructure = PDB.Structure.Structure("master")
        
        i = 0
        for structure in listOfAtomicStructObjects:
            for model in structure.get_models():
                new_model=model.copy()
                new_model.id=i
                new_model.serial_num=i+1
                i = i+1
                masterStructure.add(new_model)
        self.outputStructure = masterStructure
        

    def applySymmetryStep(self):
        
        structure = self.outputStructure
        symMatrices = getSymmetryMatrices(sym=self.symGroup)
        structureList = [structure.copy() for i in range(len(symMatrices))]

        for i in range((len(symMatrices))):
            atoms = structureList[i].get_atoms()
            matrix = symMatrices[i]
            for atom in atoms:
                coord = atom.get_coord()
                
                r = _applyMatrix(matrix, coord)
                coord = r[:3]
                atom.set_coord(coord)

        outputStructure = PDB.Structure.Structure("output")
        i = 0
        for structure in structureList:
            for model in structure.get_models():
                new_model=model.copy()
                new_model.id=i
                new_model.serial_num=i+1
                i = i+1
                outputStructure.add(new_model)
        self.outputStructure = outputStructure

    def biologicalAssemblyStep(self):
        ah = AtomicStructHandler()
        io=MMCIFIO()
        structure = self.outputStructure
        chains = [i.id for i in list(structure.get_chains())]

        old_chains = list(structure.get_chains())
        new_chains_id = generate_chain_id(len(list(old_chains)))
        index = 0
        for chain in old_chains:
            print("new chain:", new_chains_id[index])
            chain.id = new_chains_id[index]
            index+=1
                

        symMatrices = getSymmetryMatrices(sym=self.symGroup)
        matricesToWrite = symMatrices[:,:3, :3]
        # print(symMatrices[:,:3, :3])
        outputModelFn = self._getTmpPath("tempretureFileBeforeBioAssembly.cif")
        io.set_structure(structure)
        io.save(outputModelFn)
        # io.save("test")
        mmDict = ah.readLowLevel(outputModelFn)
        mmDict['_pdbx_struct_assembly.id'] = '1'               
        mmDict['_pdbx_struct_assembly.details'] = 'Scipion Defined Assembly'
        mmDict['_pdbx_struct_assembly.method_details'] = SCIPION_SYM_NAME[self.symGroup]
        mmDict['_pdbx_struct_assembly_gen.assembly_id'] = '1'
        mmDict['_pdbx_struct_assembly_gen.asym_id_list'] = ','.join(chains)
        mmDict['_pdbx_struct_assembly_gen.oper_expression'] = '(1-'+str(len(symMatrices))+')'

        mmDict['_pdbx_struct_oper_list.id'] = [str(i+1) for i in range(len(symMatrices))]
        mmDict['_pdbx_struct_oper_list.type'] = ['matrix'+str(i+1) for i in range(len(symMatrices))]
        mmDict['_pdbx_struct_oper_list.name'] = ['?' for i in range(len(symMatrices))]
        mmDict['_pdbx_struct_oper_list.matrix[1][1]'] = ['%.5f'%i for i in matricesToWrite[:,0, 0]]
        mmDict['_pdbx_struct_oper_list.matrix[1][2]'] = ['%.5f'%i for i in matricesToWrite[:,0, 1]]
        mmDict['_pdbx_struct_oper_list.matrix[1][3]'] = ['%.5f'%i for i in matricesToWrite[:,0, 2]]
        mmDict['_pdbx_struct_oper_list.vector[1]'] = ['0.0000000000' for i in range(len(symMatrices))]
        mmDict['_pdbx_struct_oper_list.matrix[2][1]'] = ['%.5f'%i for i in matricesToWrite[:,1, 0]]
        mmDict['_pdbx_struct_oper_list.matrix[2][2]'] = ['%.5f'%i for i in matricesToWrite[:,1, 1]]
        mmDict['_pdbx_struct_oper_list.matrix[2][3]'] = ['%.5f'%i for i in matricesToWrite[:,1, 2]]
        mmDict['_pdbx_struct_oper_list.vector[2]'] = ['0.0000000000' for i in range(len(symMatrices))]
        mmDict['_pdbx_struct_oper_list.matrix[3][1]'] = ['%.5f'%i for i in matricesToWrite[:,2, 0]]
        mmDict['_pdbx_struct_oper_list.matrix[3][2]'] = ['%.5f'%i for i in matricesToWrite[:,2, 1]]
        mmDict['_pdbx_struct_oper_list.matrix[3][3]'] = ['%.5f'%i for i in matricesToWrite[:,2, 2]]
        mmDict['_pdbx_struct_oper_list.vector[3]'] = ['0.0000000000' for i in range(len(symMatrices))]

        self.ouputDict = mmDict

        # print(type(mmDict['_atom_site.id']))
        
        # print(mmDict['_atom_site.id'])

    def createOutputStep(self):
        if self.outputOnlyMatrices:
            ah = AtomicStructHandler()
            outputModelFn = self._getOutputFileName()
            ah._writeLowLevel(outputModelFn, self.ouputDict)

        else:
            io=MMCIFIO()
            model = self.outputStructure
            outputModelFn = self._getOutputFileName()

            io.set_structure(model)
            io.save(outputModelFn)
            self._defineOutputs(outputModel = model)
        raise Exception('exceptionnnn')
        
    def readVector(self, index):
        
        length = self.subVolCenterVec[index].get_length()            
        [shiftX, shiftY, shiftZ] = [x * length for x in self.subVolCenterVec[index].vector]
        rotMatrix = self.subVolCenterVec[index].get_matrix()
        return shiftX, shiftY, shiftZ, rotMatrix


    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        fileNum = len(glob(self.pdbFile.get()))

        if self.defineVector == CMM:
            self.vectorsFile = self.vectorFile.get()
            vectorNum = len(vectors_from_cmm(self.vectorsFile, 1))
        else:
            vectorNum = len(vectors_from_string(self.vector.get()))

        distanceNum = len(distances_from_string(str(self.length.get())))
        if not (fileNum == vectorNum or vectorNum == 1) and (distanceNum ==0 or distanceNum ==1 or distanceNum == vectorNum):
            errMsg = "The number of vectors, distances and files you enter must be one or equal to the number of files entered(Number of Distances can be 0).\n Number of vectors: "+ str(vectorNum)+"\nNumber of files: "+ str(fileNum)+"\nNumber of Distances: "+ str(distanceNum)
            validateMsgs.append(errMsg)
        return validateMsgs

    def _warnings(self):

        warningsMsgs = []

        return warningsMsgs

    def _citations(self):
        return ['Ilca2015', 'Abrishami2020']

    def _summary(self):
    
        listObj = self.inputAtomicStructs
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
        return self._getTmpPath('output_%s%s%s%s.cif'
                                  % (imgType, auxString3,
                                     auxString2, auxString))
    def _getOutputFileName(self):
        return self._getExtraPath('output_model.cif')
