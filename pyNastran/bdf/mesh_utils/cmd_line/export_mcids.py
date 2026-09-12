from __future__ import annotations
import os
import sys
import argparse
from cpylog import SimpleLogger

import pyNastran
from .utils import filter_no_args


def cmd_line_export_mcids(argv=None, quiet: bool=False) -> None:
    """command line interface to export_mcids"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    filter_no_args("bdf export_mcids: use 'bdf export_mcids -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf export_mcids',
        description='Export material coordinate system IDs to CSV',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default='mcids.csv',
                        metavar='OUT_CSV_FILENAME',
                        help='path to output CSV file (default=mcids.csv)')
    parser.add_argument('--iplies', default=None, metavar='PLIES',
                        help='the plies indices to export; comma separated (default=0)')

    axis_group = parser.add_mutually_exclusive_group()
    axis_group.add_argument('--no_x', action='store_true',
                            help="don't write the x axis")
    axis_group.add_argument('--no_y', action='store_true',
                            help="don't write the y axis")

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename = args.IN_BDF_FILENAME
    csv_filename_in = args.output

    export_xaxis = not args.no_x
    export_yaxis = not args.no_y
    csv_filename_base = os.path.splitext(csv_filename_in)[0]
    iplies = [0]
    if args.iplies is not None:
        iplies = [int(iply) for iply in args.iplies.split(',')]
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
