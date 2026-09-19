from __future__ import annotations
import sys
import argparse
from cpylog import SimpleLogger

import pyNastran
from .utils import filter_no_args


def cmd_line_split_cbars_by_pin_flag(argv=None, quiet: bool=False) -> None:
    """command line interface to split_cbars_by_pin_flag"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    filter_no_args("bdf split_cbars_by_pin_flags: use 'bdf split_cbars_by_pin_flags -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf split_cbars_by_pin_flags',
        description='Split CBAR elements by pin flag',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default='model_new.bdf',
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF file (default=model_new.bdf)')
    parser.add_argument('-p', '--pin', default='pin_flags.csv',
                        metavar='PIN_FLAGS_CSV_FILENAME',
                        help='path to pin_flags_csv file (default=pin_flags.csv)')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename_in = args.IN_BDF_FILENAME
    punch = args.punch
    bdf_filename_out = args.output
    pin_flags_filename = args.pin

    from pyNastran.bdf.mesh_utils.split_cbars_by_pin_flag import split_cbars_by_pin_flag
    split_cbars_by_pin_flag(bdf_filename_in, pin_flags_filename=pin_flags_filename,
                            bdf_filename_out=bdf_filename_out,
                            punch=punch)
