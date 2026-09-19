import sys
import argparse

import pyNastran
from cpylog import SimpleLogger

from .utils import filter_no_args


def cmd_line_stats(argv=None, quiet: bool = False) -> None:
    """list the cards"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    filter_no_args("bdf stats: use 'bdf stats -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf stats',
        description='List the cards in a BDF model',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')

    args = parser.parse_args(argv[2:])

    bdf_filename = args.IN_BDF_FILENAME
    punch = args.punch
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    if not quiet:  # pragma: no cover
        print(vars(args))

    from pyNastran.bdf.bdf import read_bdf, BDF
    model: BDF = read_bdf(bdf_filename, validate=True, xref=True, punch=punch,
                          encoding=None, log=log, debug=True, mode='msc')
    msg = model.get_bdf_stats()
    if not quiet:  # pragma: no cover
        print(msg)
