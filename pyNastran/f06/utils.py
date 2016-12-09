from __future__ import print_function
import os
import numpy as np
from pyNastran.f06.parse_flutter import plot_flutter_f06


def cmd_line_plot_flutter():  # pragma: no cover
    import sys
    from docopt import docopt
    import pyNastran
    msg = "Usage:\n"
    msg += "  f06 plot_145 F06_FILENAME [--noline] [--modes MODES] [--subcases SUB] [--xlim FREQ] [--ylim DAMP]\n"

    msg += '  f06 plot_145 -h | --help\n'
    msg += '  f06 plot_145 -v | --version\n'
    msg += '\n'

    msg += "Positional Arguments:\n"
    msg += "  F06_FILENAME    path to input F06 files\n"

    msg += 'Options:\n'
    msg += "  --modes MODES   the modes to plot (e.g. 1:10,20:22); unused\n"
    msg += "  --subcases SUB  the subcases to plot (e.g. 1,3); unused\n"
    msg += "  --xlim FREQ     the frequency limits (unused)"
    msg += "  --ylim DAMP     the damping limits (unused)"
    msg += '\n'

    msg += 'Info:\n'
    msg += '  -h, --help      show this help message and exit\n'
    msg += "  -v, --version   show program's version number and exit\n"

    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    print(data)
    f06_filename = data['F06_FILENAME']
    if not f06_filename.lower().endswith('.f06'):
        base = os.path.splitext(f06_filename)[0]
        f06_filename = base + '.f06'

    modes = data['--modes']
    modes2 = []
    if modes is not None:
        smodes = modes.strip().split(',')
        for mode in smodes:
            mode = mode.strip()
            if ':' in mode:
                smode = mode.split(':')
                if len(smode) == 2:
                    istart = int(smode[0])
                    if smode[1] == '':
                        iend = None
                        modes2 = slice(istart, None)
                        assert len(smodes) == 1, smodes
                    else:
                        iend = int(smode[1])
                        assert iend > istart, 'smode=%s; istart=%s iend=%s' % (smode, istart, iend)
                        modes2 += list(range(istart, iend + 1))
                elif len(smode) == 3:
                    istart = int(smode[0])
                    iend = int(smode[1])
                    assert iend > istart, 'smode=%s; istart=%s iend=%s' % (smode, istart, iend)
                    istep = int(smode[2])
                    modes2 += list(range(istart, iend + 1, istep))
                else:
                    raise NotImplementedError('smode=%r; len=%s' % (smode, len(smode)))
            else:
                imode = int(mode)
                modes2.append(imode)
        #modes = np.array(modes2, dtype='int32') - 1
        modes = modes2
    print('modes = %s' % modes)
    plot_flutter_f06(f06_filename, modes=modes,
                     plot_root_locus=True, plot_vg_vf=True, plot_vg=False)
    #plot_flutter_f06(f06_filename, plot_root_locus=False, plot_vg_vf=True)

def cmd_line():  # pragma: no cover
    import sys
    dev = True
    msg = 'Usage:\n'
    msg += '  f06 plot_145 F06_FILENAME [--noline] [--modes MODES] [--subcases SUB] [--xlim FREQ] [--ylim DAMP]\n'
    msg += '\n'
    msg += '  f06 plot_145 -h | --help\n'
    msg += '  f06 -v | --version\n'
    msg += '\n'

    if len(sys.argv) == 1:
        sys.exit(msg)

    #assert sys.argv[0] != 'bdf', msg

    if sys.argv[1] == 'plot_145':
        cmd_line_plot_flutter()
    else:
        sys.exit(msg)
        #raise NotImplementedError('arg1=%r' % sys.argv[1])

if __name__ == '__main__':  # pragma: no cover
    main()

