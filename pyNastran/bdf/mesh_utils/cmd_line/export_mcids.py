from __future__ import annotations
import os
import sys
from cpylog import SimpleLogger
from docopt import docopt

import pyNastran
from .utils import filter_no_args


def cmd_line_export_mcids(argv=None, quiet: bool=False) -> None:
    """command line interface to export_mcids"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    msg = (
        'Usage:\n'
        '  bdf export_mcids IN_BDF_FILENAME [-o OUT_CSV_FILENAME] [--iplies PLIES] [--no_x | --no_y]\n'
        '  bdf export_mcids -h | --help\n'
        '  bdf export_mcids -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        '  -o OUT, --output  OUT_CSV_FILENAME  path to output CSV file\n'
        '  --iplies PLIES                      the plies indices to export; comma separated (default=0)\n'
        '\n'

        'Data Suppression:\n'
        "  --no_x,  don't write the x axis\n"
        "  --no_y,  don't write the y axis\n"
        '\n'

        'Info:\n'
        '  -h, --help      show this help message and exit\n'
        "  -v, --version   show program's version number and exit\n"
    )
    filter_no_args(msg, argv, quiet=quiet)

    ver = str(pyNastran.__version__)
    #type_defaults = {
    #    '--nerrors' : [int, 100],
    #}
    data = docopt(msg, version=ver, argv=argv[1:])
    if not quiet:  # pragma: no cover
        print(data)
    #size = 16
    bdf_filename = data['IN_BDF_FILENAME']
    csv_filename_in = data['--output']
    if csv_filename_in is None:
        csv_filename_in = 'mcids.csv'

    export_xaxis = True
    export_yaxis = True
    if data['--no_x']:
        export_xaxis = False
    if data['--no_y']:
        export_yaxis = False
    csv_filename_base = os.path.splitext(csv_filename_in)[0]
    iplies = [0]
    if data['--iplies']:
        iplies = data['--iplies'].split(',')
        iplies = [int(iply) for iply in iplies]
        if not quiet:  # pragma: no cover
            print('iplies = %s' % iplies)

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')

    from pyNastran.bdf.bdf import read_bdf
    from pyNastran.bdf.mesh_utils.export_mcids import export_mcids

    model = read_bdf(bdf_filename, log=log, xref=False)
    model.safe_cross_reference()

    for iply in iplies:
        csv_filename = csv_filename_base + f'_ply={iply:d}.csv'
        export_mcids(model, csv_filename,
                     export_xaxis=export_xaxis, export_yaxis=export_yaxis, iply=iply)
        model.log.info(f'wrote {csv_filename}')
