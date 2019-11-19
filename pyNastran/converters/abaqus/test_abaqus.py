"""
``test_abaqus`` runs multiple checks on a ABAQUS in order to make sure that:
  - materials exist

As such, ``test_abaqus`` is very useful for debugging models.

"""
import os
import sys
import numpy as np
#warnings.simplefilter('always')

import pyNastran
#from pyNastran.utils.numpy_utils import integer_types
from pyNastran.converters.abaqus.abaqus import read_abaqus

np.seterr(all='raise')


def run_abaqus(abaqus_filename, write_abaqus=True, debug=False):
    """
    Runs a single abaqus deck

    Parameters
    ----------
    abaqus_filename : str
       the abaqus filename to read

    """
    fem1 = read_abaqus(abaqus_filename, debug=debug, log=None)
    if write_abaqus:
        base, ext = os.path.splitext(abaqus_filename)
        abqaqus_filename_out = '%s.test_abqaus%s' % (base, ext)
        fem1.write(abqaqus_filename_out)



def main():
    """
    The main function for the command line ``test_abaqus`` script.
    """
    #encoding = sys.getdefaultencoding()
    from docopt import docopt
    msg = "Usage:\n"
    msg += "  test_abaqus ABAQUS_FILENAME [-d] [-w]\n"
    msg += '  test_abaqus -h | --help\n'
    msg += '  test_abaqus -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  ABAQUS_FILENAME   path to Abaqus INP file\n"
    msg += '\n'
    msg += "Options:\n"
    msg += "  -d, --debug  debug mode\n"
    msg += "  -w, --write  write test.test_abaqus.inp\n"
    msg += '\n'

    #msg += 'Options:\n'
    #msg += "\n"
    msg += "Info:\n"
    msg += '  -h, --help     show this help message and exit\n'
    msg += "  -v, --version  show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
        #'--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)

    #print(data)
    import time
    time0 = time.time()
    run_abaqus(
        data['ABAQUS_FILENAME'],
        data['--write'],
        debug=data['--debug'],
    )
    print("total time:  %.2f sec" % (time.time() - time0))


if __name__ == '__main__':  # pragma: no cover
    main()
