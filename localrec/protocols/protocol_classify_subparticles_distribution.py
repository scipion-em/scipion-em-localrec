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
from pwem.objects import SetOfCoordinates, Class2D, SetOfParticles
import numpy as np
from scipy.spatial.distance import cdist


class ProtLocalizedClassifyDistribution(ProtParticlePicking, ProtParticles):
    """ This protocol carries out a classification the distibution of subparticles.
    Consider a protein and several subparticles in it. This algorithm searches in the
    set of subparticles different geometrical distribution of the subparticles in the
    protein.
    """
    _label = 'classify subparticle distribution'
    OUTPUTCLASSESNAME = "outputClasses"

    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)
        ProtParticles.__init__(self, **args)
        self.allowMpi = False
        self.allowThreads = False

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfParticles',
                      important=True,
                      label="Subparticles",
                      help='Select the input set of coordinates of subparticles to be classified.')
        form.addParam('maxNumSubpart', IntParam,
                      important=True,
                      default=60,
                      label="Maximum number of subparticles per particle",
                      help='Maximum number of subparticles per particle.')
        form.addParam('classfyByDistance', BooleanParam,
                      important=True,
                      default=False,
                      label="Classify by distance",
                      help='Maximum number of subparticles per particle.')
        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)

        if self.classfyByDistance.get():
            self._insertFunctionStep(self.classifyByDistanceStep)
        else:
            self._insertFunctionStep(self.groupByCardinalOfClasses)
            self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        inputSubParticles = self.inputSet.get()

        # The set of coordinates are read and stored as a dictonary. All subparticles with the same micId
        # are stored in the same dictionary key
        # self.coordsDict = {}
        self.subPartsDict = {}

        for part in inputSubParticles.iterItems():
            particleKey = str(part.getMicId())
            coord = part.getCoordinate()
            #print(coord.getPosition())
            # if particleKey not in self.coordsDict:
            #    self.coordsDict[particleKey] = []
            if particleKey not in self.subPartsDict:
                self.subPartsDict[particleKey] = []
            # self.coordsDict[particleKey].append(coord.getPosition())
            self.subPartsDict[particleKey].append(part.clone())

    def groupByCardinalOfClasses(self):

        self.subparticlesclasses = []
        for i in range(0, self.maxNumSubpart.get()): self.subparticlesclasses.append([])
        # for key in self.coordsDict:
        #    subParts_in_Part = self.coordsDict[key]
        #    self.subparticlesclasses[len(subParts_in_Part)-1].append(subParts_in_Part)
        for key in self.subPartsDict:
            subParts_in_Part = self.subPartsDict[key]
            self.subparticlesclasses[len(subParts_in_Part) - 1].append(subParts_in_Part)

    def createOutputStep(self):
        inputSubParticles = self.inputSet.get()

        #print(self.subparticlesclasses)
        outputParticles = self._createSetOfParticles()
        #outputParticles = SetOfParticles()
        idx = 1
        for i in range(0, self.maxNumSubpart.get()):
            for subpartList in self.subparticlesclasses[i]:
                for subpart in subpartList:
                    #print(subpart)
                    subpart.setObjId(idx)
                    subpart.setClassId(i+1)
                    idx += 1
                    outputParticles.append(subpart.clone())
            print('---------------------------------------------------')
        outputParticles.write()
        #classesSet = self._createSetOfClasses2D(outputParticles)

        #classesSet.setImages(outputParticles)
        #self._fillClassesFromLevel(classesSet)
        #classesSet.classifyItems()


        #self._defineOutputs(**{self.OUTPUTCLASSESNAME: classesSet})
        #self._defineSourceRelation(inputSubParticles, classesSet)

        self._defineOutputs(**{self.OUTPUTCLASSESNAME: outputParticles})
        self._defineSourceRelation(inputSubParticles, outputParticles)

    def classifyByDistanceStep(self):

        sumDistance = {}

        for particleKey in self.subPartsDict:
            listCoords = []
            for subpart in self.subPartsDict[particleKey]:
                listCoords.append(subpart.getCoordinate().getPosition())
            print(listCoords)
            print('........................................................................')
            #coord = self.subPartsDict[particleKey].getCoordinate()
            dist = 0
            # We compute the distances between all possible pairs of coordinates of subparticles in the same Particles
            distances = cdist(listCoords, listCoords, metric='euclidean')
            Npairs = 0
            for i in range(0, len(listCoords)):
                for j in range(i + 1, len(listCoords)):
                    dist += distances[i, j]
                    Npairs += 1

            # The distance is sotred in sumDistance. And particleKey is the particleID or micId
            sumDistance[particleKey] = dist / Npairs

        # The dictionary is sorted
        sorted_dict = dict(sorted(sumDistance.items(), key=lambda item: item[1]))
        print(sorted_dict)

        lastValue = 0
        listOfClasses = []
        newclass = []
        threshold = 1

        for key, value in sorted_dict.items():
            if np.abs(value - lastValue) <= threshold:
                newclass.append(key)
            else:
                listOfClasses.append(newclass)
                newclass = []
                newclass.append(key)
            lastValue = value

        listOfClasses.append(newclass)

        print(listOfClasses)
        # Print the sorted dictionary


        # Creating the output

        inputSubParticles = self.inputSet.get()

        outputParticles = self._createSetOfParticles()

        idx = 1
        for i in range(0, len(listOfClasses)):
            cl = listOfClasses[i]
            for key in cl:
                print(key)
                subpart = self.subPartsDict[key]
                for elem in subpart:
                    print(elem)
                    elem.setObjId(idx)
                    elem.setClassId(i)
                    outputParticles.append(elem)
                    idx += 1

            print('---------------------------------------------------')
        outputParticles.write()
        #classesSet = self._createSetOfClasses2D(outputParticles)

        #classesSet.setImages(outputParticles)
        #self._fillClassesFromLevel(classesSet)
        #classesSet.classifyItems()


        #self._defineOutputs(**{self.OUTPUTCLASSESNAME: classesSet})
        #self._defineSourceRelation(inputSubParticles, classesSet)

        self._defineOutputs(**{self.OUTPUTCLASSESNAME: outputParticles})
        self._defineSourceRelation(inputSubParticles, outputParticles)




    def clasifyStep(self):
        inputSubParticles = self.inputSet.get()

        # The set of coordinates are read and stored as a dictonary. All subparticles with the same micId
        # are stored in the same dictionary key
        coordsDict = {}
        for part in inputSubParticles:
            particleKey = str(part.getMicId())
            coord = part.getCoordinate()
            if particleKey not in coordsDict:
                coordsDict[particleKey] = []
            coordsDict[particleKey].append(coord.getPosition())

        sumDistance = {}

        for particleKey in coordsDict:
            dist = 0
            # We compute the distances between all possible pairs of coordinates of subparticles in the same Particles
            distances = cdist(coordsDict[particleKey], coordsDict[particleKey], metric='euclidean')
            Npairs = 0
            for i in range(0, len(coordsDict[particleKey])):
                for j in range(i + 1, len(coordsDict[particleKey])):
                    dist += distances[i, j]
                    Npairs += 1

            # The distance is sotred in sumDistance. And particleKey is the particleID or micId
            sumDistance[particleKey] = dist / Npairs

        # The dictionary is sorted
        sorted_dict = dict(sorted(sumDistance.items(), key=lambda item: item[1]))
        print(sorted_dict)

        lastValue = 0
        listOfClasses = []
        newclass = []
        threshold = 1

        for key, value in sorted_dict.items():
            if np.abs(value - lastValue) <= threshold:
                newclass.append(key)
            else:
                listOfClasses.append(newclass)
                newclass = []
                newclass.append(key)
            lastValue = value

        listOfClasses.append(newclass)

        print(listOfClasses)
        # Print the sorted dictionary

        '''
        for i in range(0, len(keyCandidate)):
            if np.abs(sumDistance[particleKey]-sumDistance[keyCandidate]) < 1.0:
                listOfClasses.append(sumDistance[particleKey])
        keyCandidate.append(particleKey)
        '''

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
