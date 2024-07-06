import sys

def cmd_line_merge(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf_merge"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    from docopt import docopt
    import pyNastran
    msg = (
        "Usage:\n"
        '  bdf merge (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME] [--debug]\n'
        '  bdf merge -h | --help\n'
        '  bdf merge -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAMES   path to input BDF/DAT/NAS files\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    if len(argv) == 1:
        sys.exit(msg)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver, argv=argv[1:])
    if not quiet:  # pragma: no cover
        print(data)
    size = 16
    bdf_filenames = data['IN_BDF_FILENAMES']
    bdf_filename_out = data['--output']
    debug = data['--debug']
    assert debug in {True, False}, debug
    if bdf_filename_out is None:
        bdf_filename_out = 'merged.bdf'

    #cards_to_skip = [
        #'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
        #'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    cards_to_skip = []

    from cpylog import SimpleLogger
    from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
    log = SimpleLogger(level='debug')

    bdf_merge(bdf_filenames, bdf_filename_out, renumber=True,
              encoding=None, size=size, is_double=False, cards_to_skip=cards_to_skip,
              log=log)
