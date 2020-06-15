"""
defines:
 - cmd_line_plot_flutter()

"""
#import matplotlib
#matplotlib.use('Qt5Agg')
#from pyNastran.gui.matplotlib_backend import  matplotlib_backend
#print(matplotlib_backend)
#matplotlib.use(matplotlib_backend)
#from pyNastran.gui.qt_version import qt_version

PLOT_TYPES = '[--eas|--tas]'
USAGE = (
    'Usage:\n'
    '  f06 plot_145 F06_FILENAME [--noline] [--modes MODES] [--subcases SUB] [--xlim XLIM] [--ylimdamp DAMP] [--ylimfreq FREQ]'
    f'{PLOT_TYPES} [--kfreq] [--rootlocus] [--in_units IN] [--out_units OUT] [--nopoints] [--export] [--f06]\n'
)
def cmd_line_plot_flutter(argv=None, plot=True, show=True, log=None):
    """the interface to ``f06 plot_145`` on the command line"""
    import sys
    import os
    from docopt import docopt
    import pyNastran
    from pyNastran.f06.parse_flutter import plot_flutter_f06
    if argv is None:
        argv = sys.argv
    msg = (
        USAGE +
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
        '  --xlim XLIM      the velocity limits (default=no limit)\n'
        '  --ylimfreq FREQ  the damping limits (default=no limit)\n'
        '  --ylimdamp DAMP  the damping limits (default=-0.3:0.3)\n'
        "  --nopoints       don't plot the points\n"
        "  --noline         don't plot the lines\n"
        "  --export         export a zona file\n"
        "  --f06            export an F06 file (temporary)\n"
        '\n'
        'Info:\n'
        '  -h, --help      show this help message and exit\n'''
        "  -v, --version   show program's version number and exit\n"
    )
    if len(argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver, argv=argv[1:])
    f06_filename = data['F06_FILENAME']
    if not f06_filename.lower().endswith('.f06'):
        base = os.path.splitext(f06_filename)[0]
        f06_filename = base + '.f06'

    modes = split_int_colon(data['--modes'], start_value=1)

    xlim = [None, None]
    if data['--xlim']:
        xlim = split_float_colons(data['--xlim'])

    ylim_damping = [-.3, .3]
    if data['--ylimdamp']:
        ylim_damping = split_float_colons(data['--ylimdamp'])

    ylim_freq = [None, None]
    if data['--ylimfreq']:
        ylim_freq = split_float_colons(data['--ylimfreq'])

    in_units = 'si'
    if data['--in_units']:
        in_units = data['IN'].lower()
    assert in_units in ['si', 'english_in', 'english_ft', 'english_kt'], 'in_units=%r' % in_units

    out_units = 'si'
    if data['--out_units']:
        out_units = data['OUT'].lower()
    assert out_units in ['si', 'english_in', 'english_ft', 'english_kt'], 'out_units=%r' % out_units

    plot_type = 'tas'
    if data['--eas']:
        plot_type = 'eas'
    elif data['--tas']:
        plot_type = 'tas'
    else:
        sys.stderr.write('plot_type assumed to be --tas\n')

    plot_kfreq_damping = data['--kfreq']
    plot_root_locus = data['--rootlocus']

    nopoints = data['--nopoints']
    noline = data['--noline']

    export_f06 = data['--f06']
    export_zona = data['--export']
    export_f06_filename = None if export_f06 is False else 'nastran.f06'
    export_zona_filename = None if export_zona is False else 'nastran.zona'
    export_veas_filename = None if export_zona is False else 'nastran.veas'

    # TODO: need a new parameter
    vg_filename = None if export_zona is  None else 'vg_subcase_%i.png'
    vg_vf_filename = None if export_zona is  None else 'vg_vf_subcase_%i.png'
    kfreq_damping_filename = None if export_zona is  None else 'kfreq_damping_subcase_%i.png'
    root_locus_filename = None if export_zona is  None else 'root_locus_subcase_%i.png'
    if not plot:
        return
    plot_flutter_f06(f06_filename, modes=modes,
                     plot_type=plot_type,
                     f06_units=in_units,
                     out_units=out_units,
                     plot_root_locus=plot_root_locus, plot_vg_vf=True, plot_vg=False,
                     plot_kfreq_damping=plot_kfreq_damping,
                     xlim=xlim,
                     ylim_damping=ylim_damping, ylim_freq=ylim_freq,
                     nopoints=nopoints,
                     noline=noline,
                     export_veas_filename=export_veas_filename,
                     export_zona_filename=export_zona_filename,
                     export_f06_filename=export_f06_filename,
                     vg_filename=vg_filename,
                     vg_vf_filename=vg_vf_filename,
                     root_locus_filename=root_locus_filename,
                     kfreq_damping_filename=kfreq_damping_filename, show=show, log=log)

def split_float_colons(string_values):
    """
    Uses numpy-ish syntax to parse a set of two floats.  Blanks are interpreted as None.

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

def split_int_colon(modes, nmax=1000, start_value=0):
    """
    Uses numpy-ish syntax to parse a set of integers.  Values are inclusive.
    Blanks are interpreted as 0 unless start_value is specified.

    Parse the following:
       1:10
       1:
       :100
       1:5:2
       :5,11:15:2,10:20
       1,3:

    """
    modes2 = []
    if modes is not None:
        smodes = modes.strip().split(',')
        for mode in smodes:
            mode = mode.strip()
            if ':' in mode:
                smode = mode.split(':')
                if len(smode) == 2:
                    #start_str, stop_str = smode
                    if smode[0] == '':
                        smode[0] = start_value
                    istart = int(smode[0])

                    if smode[1] == '':
                        iend = None
                        modes2 = slice(istart, nmax)
                        assert len(smode) == 2, smode
                    else:
                        iend = int(smode[1])
                        assert iend > istart, 'smode=%s; istart=%s iend=%s' % (smode, istart, iend)
                        modes2 += list(range(istart, iend + 1))
                elif len(smode) == 3:
                    #start_str, stop_str, step_str = smode
                    if smode[0] == '':
                        smode[0] = start_value
                    istart = int(smode[0])

                    if smode[2] == '':
                        smode[2] = 1
                    istep = int(smode[2])

                    if smode[1] == '':
                        #iend = None
                        modes2 += list(range(istart, nmax, istep))
                    else:
                        iend = int(smode[1]) + 1
                        assert iend > istart, 'smode=%s; istart=%s iend=%s' % (smode, istart, iend)
                        modes2 += list(range(istart, iend, istep))

                else:
                    raise NotImplementedError('mode=%r; len=%s' % (mode, len(smode)))
            else:
                imode = int(mode)
                modes2.append(imode)
        #modes = np.array(modes2, dtype='int32') - 1
        modes = modes2

    #print('modes =', list(modes))
    #if None not in modes:
        #modes.sort()
    try:
        modes.sort()
    except AttributeError:
        pass
    return modes

def cmd_line(argv=None, plot=True, show=True, log=None):
    """the interface to ``f06`` on the command line"""
    import sys
    if argv is None:
        argv = sys.argv

    msg = (
        USAGE +
        '\n'
        '  f06 plot_145 -h | --help\n'
        '  f06 -v | --version\n'
        '\n'
    )
    if len(argv) == 1:
        sys.exit(msg)

    #assert sys.argv[0] != 'bdf', msg

    if argv[1] == 'plot_145':
        cmd_line_plot_flutter(argv=argv, plot=plot, show=show, log=log)
    else:
        sys.exit(msg)
        #raise NotImplementedError('arg1=%r' % argv[1])

if __name__ == '__main__':  # pragma: no cover
    cmd_line()
