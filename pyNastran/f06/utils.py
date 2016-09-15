from parse_flutter import plot_flutter_f06


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
    msg += "  --modes MODES   the modes to plot (e.g. 1:10,20:22)\n"
    msg += "  --subcases SUB  the subcases to plot (e.g. 1,3)\n"
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

    plot_flutter_f06(f06_filename, plot=True)

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

