# **************************************************************************
# *
# * Authors:    Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
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

"""
This sub-package contains data and protocol classes
wrapping Localized recontruction of subunits.
"""
import os
import pyworkflow.em
from pyworkflow.utils import Environ

from bibtex import _bibtex # Load bibtex dict with references
from localrec.constants import *
from localrec.convert import *


# _logo = "opic_logo.png"
_references = ['Ilca2015']

class Plugin(pyworkflow.em.Plugin):
    _homeVar = LOCALREC_HOME
    _pathVars = [LOCALREC_HOME]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(LOCALREC_HOME, 'localrec-1.2.0')

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch localrec. """
        environ = Environ(os.environ)
        print("getEnvirion(): %s"%os.environ.get(cls.getHome()))
        if ('%s' % cls.getHome()) in environ:
            environ.update({
                'PATH': cls.getHome(),
                'LD_LIBRARY_PATH': os.environ[cls.getHome()],
            }, position=Environ.BEGIN)
        else:
            # TODO: Find a generic way to warn of this situation
            print("%s variable not set on environment." % cls.getHome())
        return environ

    @classmethod
    def validateInstallation(cls):
        """ This function will be used to check if RELION package is properly
            installed."""

        # missingPaths = ["%s: %s" % (var, os.environ[var])
        #                 for var in [cls._homeVar]
        #                 if not os.path.exists(os.environ[var])]
        #
        # if missingPaths:
        #     return ["Missing variables:"] + missingPaths
        # else:
        #     print("returning no errors")
        return []  # No errors

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(V1_2_0)

    @classmethod
    def defineBinaries(cls, env):

         # Add localrec
         env.addPackage('localrec', version='1.2.0',
                        tar='localrec-1.2.0.tgz',
                        default=True)


pyworkflow.em.Domain.registerPlugin(__name__)
