import os
import sys
import argparse

import pyNastran
from cpylog import SimpleLogger
from .utils import filter_no_args


def cmd_line_collapse_quads(argv=None, quiet: bool=False) -> None:
    """command line interface to ``delete_bad_shells``"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    filter_no_args("bdf collapse_quads: use 'bdf collapse_quads -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf collapse_quads',
        description='Collapse bad quads to tris',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default=None,
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF/DAT/NAS file')
    parser.add_argument('--size', type=int, default=8,
                        help='size of the output (default=8)')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename = args.IN_BDF_FILENAME
    punch = args.punch
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')

    base, ext = os.path.splitext(bdf_filename)
    bdf_filename_out = base + '_collapsed' + ext
    if args.output is not None:
        bdf_filename_out = args.output

    size = args.size

    from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris
    from pyNastran.bdf.bdf import read_bdf, BDF
    model: BDF = read_bdf(bdf_filename, xref=False,
                          validate=False, punch=punch, log=log)
    convert_bad_quads_to_tris(model)
    model.write_bdf(bdf_filename_out, size=size)
