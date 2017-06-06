# pylint: disable=W0612
"""
``test_bdf`` runs multiple checks on a BDF in order to make sure that:
  - no data is lost on IO
  - card field types are correct (e.g. node_ids are integers)
  - various card methods (e.g. Area) work correctly

As such, ``test_abaqus`` is very useful for debugging models.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
from six import iteritems
import numpy as np
#warnings.simplefilter('always')

np.seterr(all='raise')

from pyNastran.utils import print_bad_path, integer_types
from pyNastran.converters.abaqus.abaqus import read_abaqus

import pyNastran.bdf.test
test_path = pyNastran.bdf.test.__path__[0]


def run_abaqus(abaqus_filename):
    """
    Runs a single abaqus deck

    Parameters
    ----------
    abaqus_filename : str
       the abaqus filename to read
    """

    debug = False
    fem1 = read_abaqus(abaqus_filename, debug=debug, log=None)


def main():
    """
    The main function for the command line ``test_bdf`` script.
    """
    encoding = sys.getdefaultencoding()
    from docopt import docopt
    msg = "Usage:\n"
    msg += "  test_abaqus ABAQUS_FILENAME\n"
    msg += '  test_abaqus -h | --help\n'
    msg += '  test_abaqus -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  ABAQUS_FILENAME   path to Abaqus INP file\n"
    msg += '\n'

    #msg += 'Options:\n'
    #msg += "\n"
    msg += "Info:\n"
    msg += '  -h, --help     show this help message and exit\n'
    msg += "  -v, --version  show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    type_defaults = {
        '--nerrors' : [int, 100],
    }
    data = docopt(msg, version=ver)

    #print(data)
    import time
    t0 = time.time()
    run_abaqus(
        data['ABAQUS_FILENAME'],
    )
    print("total time:  %.2f sec" % (time.time() - t0))


if __name__ == '__main__':  # pragma: no cover
    main()
