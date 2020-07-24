# **************************************************************************
# *
# * Authors:  Roberto Marabini (roberto@cnb.csic.es)
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
from localrec.protocols.protocol_localized_set_origin \
    import ProtLocalOrigSampling
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer
from pwem.viewers.viewer_chimera import (Chimera)
from pwem.emlib.image import ImageHandler
from ..protocols.protocol_localized_set_origin import ProtLocalOrigSampling

import os

class ProtLocalOrigSamplingViewer(Viewer):
    """Visualize the output of ProtLocalOrigSampling"""
    _label = 'viewer localized set origin'
    _targets = [ProtLocalOrigSampling]
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **args):
        # Create minimalistic coordinate axis
        dim = 150.
        sampling = 1.
        bildFileName = self.protocol._getTmpPath("axis_output.bild")
        Chimera.createCoordinateAxisFile(dim,
                                         bildFileName=bildFileName,
                                         sampling=sampling)
        fnCmd = self.protocol._getTmpPath("chimera_output.cxc")
        f = open(fnCmd, 'w')
        # change to workingDir
        # If we do not use cd and the project name has an space
        # the protocol fails even if we pass absolute paths
        f.write('cd %s\n' % os.getcwd())
        f.write("open %s\n" % bildFileName)
        f.write("cofr 0,0,0\n")  # set center of coordinates

        # show vol
        _showVol = self.protocol.outputVolume #.get()
        showVolFileName = ImageHandler.removeFileType(_showVol.getFileName())
        f.write("open %s\n" % showVolFileName)
        if _showVol.hasOrigin():
            x, y, z = _showVol.getOrigin().getShifts()
        else:
            x, y, z = _showVol.getOrigin(force=True).getShifts()

        f.write("volume #2 style surface voxelSize %f\n"
                "volume #2 origin %0.2f,%0.2f,%0.2f\n"
                % (_showVol.getSamplingRate(), x, y, z))
        # show atom struct)
        atomstruct = self.protocol.atomStruct.get()
        if atomstruct is not None:
            f.write("open %s\n" % atomstruct.getFileName())
        # run in the background
        Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")
        return []

