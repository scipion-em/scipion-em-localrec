# **************************************************************************
# *
# * Authors:     Jose Luis Vilas
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
import matplotlib.pyplot as plt
import numpy as np
from pyworkflow.viewer import (ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO)
from ..protocols.protocol_subparticles_distribution import ProtLocalizedSubparticleDistribution2


class ProtLocalizedSubparticleDistributionViewer(ProtocolViewer):
    """ Visualization of the subparticles distribution
    """
    _label = 'viewer subparticle distribution'
    _targets = [ProtLocalizedSubparticleDistribution2]

    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, *args, **kwargs):
        protocol = kwargs['protocol']
        fnHistogram = protocol._getExtraPath('histogram.txt')
        file1 = open(fnHistogram, 'r')
        count = 0

        idxList = []
        wList = []
        while True:
            count += 1

            # Get next line from file
            line = file1.readline()
            if not line:
                break
            line = line.strip()
            parsed = line.split(",")
            idxList.append(int(parsed[0]))
            wList.append(int(parsed[1]))

        file1.close()

        print(idxList)
        plt.bar(idxList, wList, align='center')
        plt.title('Number of subparticles distribution')
        plt.xlabel('Number of subparticles')
        plt.ylabel('Number of particles')
        plt.show()

