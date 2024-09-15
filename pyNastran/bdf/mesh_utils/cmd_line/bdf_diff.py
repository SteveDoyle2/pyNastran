from __future__ import annotations
import sys
#from typing import TYPE_CHECKING
#if TYPE_CHECKING:
#    from pyNastran.bdf.bdf import BDF


def cmd_line_diff(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf_diff"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import os
    from docopt import docopt
    import pyNastran
    msg = (
        "Usage:\n"
        '  bdf diff IN_BDF_FILENAME1 IN_BDF_FILENAME2 [--punch] [--debug]\n'
        '  bdf diff -h | --help\n'
        '  bdf diff -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME1   path to input BDF/DAT/NAS files\n'
        '  IN_BDF_FILENAME2   path to input BDF/DAT/NAS files\n'
        '\n'

        'Options:\n'
        '  --punch      uses a punch file\n'
        '  --debug      increase debug logging\n\n'

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
    bdf_filename1 = data['IN_BDF_FILENAME1']
    bdf_filename2 = data['IN_BDF_FILENAME2']

    #bdf_filename_out = data['--output']
    debug = data['--debug']
    assert debug in {True, False}, debug
    # if bdf_filename_out is None:
    #     bdf_filename_out = 'merged.bdf'

    #cards_to_skip = [
        #'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
        #'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    #cards_to_skip = []

    from pyNastran.bdf.mesh_utils.bdf_diff import get_diff_bdfs
    added_cards, removed_cards, added_model, removed_model = get_diff_bdfs(bdf_filename1, bdf_filename2)
