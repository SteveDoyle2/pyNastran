from typing import Optional
from numpy import vstack
from pyNastran.converters.tecplot.tecplot import Tecplot


def merge_tecplot_files(tecplot_filenames: list[str],
                        tecplot_filename_out: Optional[str]=None,
                        log=None) -> Tecplot:
    """merges one or more tecplot files"""
    assert isinstance(tecplot_filenames, (list, tuple)), type(tecplot_filenames)
    assert len(tecplot_filenames) > 0, tecplot_filenames

    model = Tecplot(log=log)
    model.log.info('reading %s' % tecplot_filenames[0])
    model.read_tecplot(tecplot_filenames[0])
    model.log.info(f'  nzones0 = {model.nzones}')
    zones = model.zones

    if len(tecplot_filenames) == 1 and tecplot_filename_out is not None:
        model.write_tecplot(tecplot_filename_out)
        return model

    for tecplot_filename in tecplot_filenames[1:]:
        model = Tecplot(log=log)
        model.log.info('reading %s' % tecplot_filename)
        model.read_tecplot(tecplot_filename)
        model.log.info(f'  nzones = {model.nzones}')
        zones += model.zones
    model.zones = zones

    if tecplot_filename_out is not None:
        model.write_tecplot(tecplot_filename_out)
    return model

def merge_tecplot(argv: Optional[list[str]] = None) -> Tecplot:
    import sys
    if argv is None:
        argv = sys.argv[1:]

    msg = (
        'Usage:\n'
        '  merge_tecplot <TECPLOT_FILE>... [--output <OUT>] [-b]\n'
        '  merge_tecplot -h | --help\n'
        '  merge_tecplot -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  TECPLOT_FILE   path to input Tecplot file\n"
        '\n'

        'Options:\n'
        '  -o OUT, --out OUT, --output OUT  path to output Tecplot file (default=out.plt)\n'
        '  -b, --binary                     write to a binary file\n\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    import pyNastran
    ver = str(pyNastran.__version__)
    if '-v' in argv or '--version' in argv:
        sys.exit(ver)
    #if len(argv) < 2:
        #sys.exit(msg)

    import docopt
    args = docopt.docopt(docstring=msg, argv=argv, default_help=True, version=ver,
                         options_first=False, more_magic=False)

    print('args =', args)
    if args['--output'] is None:
        args['--output'] = 'out.plt'

    is_binary = args['--binary']
    tecplot_filename_out = args['--output']
    tecplot_filenames = args['<TECPLOT_FILE>']
    if len(tecplot_filenames) == 0:
        sys.exit(msg)

    print(f'tecplot_filenames={tecplot_filenames} '
          f'tecplot_filename_out={tecplot_filename_out}')
    model = merge_tecplot_files(
        tecplot_filenames,
        tecplot_filename_out=None, log=None)

    model.log.info(f'nzones_total = {model.nzones}')
    if is_binary:
        model.write_tecplot_binary(
            tecplot_filename_out, version='102')
    else:
        model.write_tecplot_ascii(
            tecplot_filename_out, res_types=None, adjust_nids=True)
    model.log.info(f'finished creating {tecplot_filename_out}')
    return model

if __name__ == '__main__':  # pragma: no cover
    merge_tecplot()
