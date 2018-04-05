from __future__ import print_function
import os
import numpy as np
from pyNastran.f06.parse_flutter import plot_flutter_f06


def cmd_line_plot_flutter():  # pragma: no cover
    import sys
    from docopt import docopt
    import pyNastran
    msg = (
        'Usage:\n'
        '  f06 plot_145 F06_FILENAME [--noline] [--modes MODES] [--subcases SUB] [--xlim FREQ] [--ylimdamp DAMP] '
        '[--eas|--tas] [--kfreq] [--rootlocus] [--in_units IN][--out_units OUT] [--nopoints]\n'
        '  f06 plot_145 -h | --help\n'
        '  f06 plot_145 -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  F06_FILENAME    path to input F06 files\n'

        'Plot Types for V-g/V-f:\n'
        #'  --vgvf           plots a V-g/V-f plot\n'
        '  --rootlocus      plots a root locus\n'
        '  --kfreq          plots a kfreq-g/kfreq-f plot\n'
        '\n'
        'Plot Types for V-g/V-f:\n'
        '  --tas            plot true airspeed (default)\n'
        '  --eas            plot eqivalent airspeed\n'
        '\n'
        'Units:\n'
        '  --in_units IN    Selects the input unit system\n'
        '                   si (kg, m, s) -> m/s\n'
        '                   english_ft (slug/ft^3, ft, s) -> ft/s\n'
        '                   english_in (slinch/in^3, in, s) -> in/s (default)\n'

        '  --out_units OUT  Selects the ouptut unit system\n'
        '                   si (kg, m, s) -> m/s\n'
        '                   english_ft (slug/ft^3, ft, s) -> ft/s\n'
        '                   english_in (slinch/in^3, in, s) -> in/s (default)\n'
        '                   english_kt (slinch/in^3, nm, s) -> knots\n'
        '\n'
        'Options:\n'
        '  --modes MODES    the modes to plot (e.g. 1:10,20:22)\n'
        '  --subcases SUB   the subcases to plot (e.g. 1,3); unused\n'
        '  --xlim FREQ      the frequency limits (unused)\n'
        '  --ylimdamp DAMP  the damping limits (default=-0.3:0.3)\n'
        "  --nopoints       don't plot the points\n"
        '\n'
        'Info:\n'
        '  -h, --help      show this help message and exit\n'''
        "  -v, --version   show program's version number and exit\n"
    )
    if len(sys.argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver)
    f06_filename = data['F06_FILENAME']
    if not f06_filename.lower().endswith('.f06'):
        base = os.path.splitext(f06_filename)[0]
        f06_filename = base + '.f06'

    modes = _split_modes(data['--modes'])
    ylim_damping = [-.3, .3]
    if data['--ylimdamp']:
        ylim_damping = split_float_colons(data['--ylimdamp'])

    in_units = 'si'
    if data['--in_units']:
        in_units = data['--in_units'].lower()
    assert in_units in ['si', 'english_in', 'english_ft', 'english_kt'], 'in_units=%r' % in_units

    out_units = 'si'
    if data['--out_units']:
        out_units = data['--out_units'].lower()
    assert out_units in ['si', 'english_in', 'english_ft', 'english_kt'], 'out_units=%r' % out_units

    plot_type = 'tas'
    if data['--eas']:
        plot_type = 'eas'

    plot_kfreq_damping = data['--kfreq']
    plot_root_locus = data['--rootlocus']

    nopoints = data['--nopoints']
    print('modes = %s' % modes)
    plot_flutter_f06(f06_filename, modes=modes,
                     plot_type=plot_type,
                     f06_units=in_units,
                     out_units=out_units,
                     plot_root_locus=plot_root_locus, plot_vg_vf=True, plot_vg=False,
                     plot_kfreq_damping=plot_kfreq_damping,
                     ylim_damping=ylim_damping, nopoints=nopoints)

def split_float_colons(string_values):
    """
    Parse the following:
       1:10
       1:
       :100
       1.1:200
       1:300.
    """
    if string_values:
        assert ',' not in string_values, string_values
        sline = string_values.split(':')
        assert len(sline) <= 2, sline

        values = []
        for svalue in sline:
            if svalue == '':
                val = None
            else:
                val = float(svalue)
            values.append(val)
    else:
        assert string_values is None, None
        values = None # all values
    return values

def _split_modes(modes):
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
    return modes

def cmd_line():  # pragma: no cover
    import sys
    msg = (
        'Usage:\n'
        '  f06 plot_145 F06_FILENAME [--noline] [--modes MODES] [--subcases SUB] [--xlim FREQ] [--ylimdamp DAMP] '
        '[--eas|--tas] [--kfreq] [--rootlocus]\n'
        '\n'
        '  f06 plot_145 -h | --help\n'
        '  f06 -v | --version\n'
        '\n'
    )

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

