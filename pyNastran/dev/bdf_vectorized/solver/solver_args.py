import os
import sys

import pyNastran
from docopt import docopt

def run_arg_parse(mode=''):
    msg = "Usage:\n"
    msg += "  pyNastran%s BDFNAME [--old=<OLD>] [--out=<OUT>]\n" % mode
    msg += '              [--k=<K_MATRIX>] [--m=<MASS_MATRIX>] [--f=<F_MATRIX>] [-d]\n'
    msg += '  pyNastran%s -h | --help\n' % mode
    msg += '  pyNastran%s -v | --version\n' % mode
    msg += "\n"
    msg += "Options:\n"
    msg += "  --k=<K_MATRIX>     Divide the debug Stiffness matrix by K_MATRIX, [default: 1.0] \n"
    msg += "  --m=<MASS_MATRIX>  Divide the debug Mass matrix by MASS_MATRIX, [default: 1.0] \n"
    msg += "  --f=<F_MATRIX>     Divide the Load matrix by F_MATRIX, [default: 1.0] \n"
    msg += "  --old=<OLD>        Save the old data, [default: no] \n"
    msg += "  --out=<OUT>        Creates out.f06, OUT.op2\n"
    msg += "  -d, --debug        Turns on debugging\n"
    msg += '\n'
    msg += 'Info:\n'
    msg += "  -h, --help        Show this help message and exit\n"

    msg += "  -q, --quiet       Prints debug messages (default: True)\n"
    msg += "  -v, --version     Show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    data = docopt(msg, version=ver)
    print(data)


    bdf_filename = data['BDFNAME']
    bdf_filename = os.path.abspath(bdf_filename)
    if data['--out']:
        bdf_base = data['--out']
    else:
        bdf_base = os.path.splitext(bdf_filename)[0]
    data['BDFBASE'] = bdf_base


    old = data['--old'].lower()
    if old not in ['no', 'yes']:
        raise RuntimeError('Use "old=no" or "old=yes".')

    try:
        m = float(data['--m'])
    except ValueError:
        raise ValueError('--m must be a float; m=%r' % data['--m'])
    if m < 1.0:
        raise ValueError('--m must be >= 1.0; m=%r' % m)

    try:
        f = float(data['--f'])
    except ValueError:
        raise ValueError('-f must be a float; f=%r' % data['--f'])
    if f < 1.0:
        raise ValueError('-f must be >= 1.0; f=%r' % f)

    try:
        k = float(data['--k'])
    except ValueError:
        raise ValueError('--k must be a float; k=%r' % data['--k'])
    if k < 1.0:
        raise ValueError('--k must be >= 1.0; k=%r' % k)
    data['--m'] = m
    data['--f'] = f
    data['--k'] = k

    #old = data['--old']
    #out = data['--out']
    #debug = not(data['--quiet'])
    #return (bdf_filename, old, out, debug)
    return data
