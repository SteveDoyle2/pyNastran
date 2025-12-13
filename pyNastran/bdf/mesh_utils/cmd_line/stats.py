import sys
from docopt import docopt

import pyNastran

from .utils import (
    get_bdf_filename_punch_log, filter_no_args,
)


def cmd_line_stats(argv=None, quiet: bool = False) -> None:
    """list the cards"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    msg = (
        'Usage:\n'
        '  bdf stats IN_BDF_FILENAME [--punch]\n'
        '  bdf stats -h | --help\n'
        '  bdf stats -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME  path to input BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        '  --punch          flag to identify a *.pch/*.inc file\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    filter_no_args(msg, argv, quiet=quiet)

    ver = str(pyNastran.__version__)
    # type_defaults = {
    #    '--nerrors' : [int, 100],
    # }
    data = docopt(msg, version=ver, argv=argv[1:])

    bdf_filename, punch, log = get_bdf_filename_punch_log(data, quiet)
    if not quiet:  # pragma: no cover
        print(data)

    from pyNastran.bdf.bdf import read_bdf, BDF
    model: BDF = read_bdf(bdf_filename, validate=True, xref=True, punch=punch,
                          encoding=None, log=log, debug=True, mode='msc')
    msg = model.get_bdf_stats()
    if not quiet:  # pragma: no cover
        print(msg)
