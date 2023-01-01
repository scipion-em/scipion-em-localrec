
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
import Bio.PDB

I1 = 1

class ProtSuperImpose(EMProtocol):
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
        group.addParam('pdbFile1', PathParam, default='',
                       label='PDB File1: ',
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

        group.addParam('pdbFile2', PathParam, default='',
                       label='PDB File2: ',
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
        
        
                       
        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):

        
        inputStepId = self._insertFunctionStep('convertInputStep')
        
        tfCoordId = self._insertFunctionStep('superImposeStep',
                                                    prerequisites=[inputStepId])

        # If user wants to get biological assembly in their outputfile, we continue with biological assembly step
        self._insertFunctionStep('createOutputStep', prerequisites=[tfCoordId])
        
        

    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        """
         Convert the input to use later on in this protocol. 
        """
        # Sort input file names alphabetically to match them to the input vectors and distances
        ah = AtomicStructHandler()
        self.inputFile1 =  glob(self.pdbFile1.get())
        self.inputFile2 =  glob(self.pdbFile2.get())
        self.inputFile1.sort()
        self.inputFile2.sort()
        self.inputFiles = self.inputFile1
        self.inputSuperImposes= self.inputFile2



    def superImposeStep(self):
        """
            Transform input structures and bind them to a master structure to make them ready for symmetrisation
        """
        
        listOfAtomicStructObjects = []
        # i = 0
        for file in self.inputFiles + self.inputSuperImposes:
            ah = AtomicStructHandler()
            print(file)
            ah.read(file)
            # matrix = ah.getTransformMatrix(self.inputSuperImposes[i])
            # i +=1
            # print(matrix)
            listOfAtomicStructObjects.append(ah.structure.copy())
        
        # for i, struct in enumerate(listOfAtomicStructObjects):
        #     c = struct.get_chains()
        #     print([i for i in c])

        start_id = 1
        end_id   = 100000
        atoms_to_be_aligned = range(start_id, end_id + 1)

        # pdb_parser = Bio.PDB.PDBParser(QUIET = True)
        # print(self.inputFiles[0])
        # print(self.inputSuperImposes[0])
        # ref_structure = pdb_parser.get_structure("reference", file = self.inputFiles[0])
        # sample_structure = pdb_parser.get_structure("sample", file = self.inputSuperImposes[0])

        ref_structure = listOfAtomicStructObjects[0]
        sample_structure = listOfAtomicStructObjects[1]
        print(ref_structure)
        print(ref_structure)
        ref_model    = ref_structure[0]
        sample_model = sample_structure[0]

        ref_atoms = []
        sample_atoms = []

        # Iterate of all chains in the model in order to find all residues
        for ref_chain in ref_model:
        # Iterate of all residues in each model in order to find proper atoms
            for ref_res in ref_chain:
                # Check if residue number ( .get_id() ) is in the list
                # if ref_res.get_id()[1] in atoms_to_be_aligned:
                # Append CA atom to list3
                for i in ref_res:
                    ref_atoms.append(i)

        # Do the same for the sample structure
        for sample_chain in sample_model:
            for sample_res in sample_chain:
                # if sample_res.get_id()[1] in atoms_to_be_aligned:
                for i in sample_res:
                    sample_atoms.append(i)
                # Now we initiate the superimposer:
        super_imposer = Bio.PDB.Superimposer()

        super_imposer.set_atoms(ref_atoms, sample_atoms)
        (rot, trans) = super_imposer.rotran
        super_imposer.apply(sample_model.get_atoms())
        # sample_model.transform(rot, trans)
        ah = AtomicStructHandler()
        outputModelFn = self._getOutputFileName()
        ah.structure = sample_model
        ah.write(outputModelFn)
        print(rot)
        print(trans)

    def createOutputStep(self):
        pass
        
    def readVector(self, index):
        
        length = self.subVolCenterVec[index].get_length()            
        [shiftX, shiftY, shiftZ] = [x * length for x in self.subVolCenterVec[index].vector]
        rotMatrix = self.subVolCenterVec[index].get_matrix()
        return shiftX, shiftY, shiftZ, rotMatrix


    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
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
