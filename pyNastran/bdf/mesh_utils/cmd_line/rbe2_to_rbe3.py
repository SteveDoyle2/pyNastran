from __future__ import annotations
import os
import sys
from typing import TYPE_CHECKING
from cpylog import SimpleLogger
import numpy as np

import pyNastran
from .utils import add_argparse_arguments
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def cmd_line_rbe2_to_rbe3(argv=None, quiet: bool=False) -> None:
    """
    rbe3_to_rbe2 filename.bdf
    """
    import argparse
    FILE = os.path.abspath(__file__)
    if argv is None:
        argv = sys.argv[1:]  # ['run_jobs'] + sys.argv[2:]
    else:
        argv = [FILE] + argv[2:]  # ['run_jobs'] + sys.argv[2:]
    if not quiet:
        print(f'argv = {argv}')
    # print(f'argv = {argv}')

    parser = argparse.ArgumentParser(prog='rbe2_to_rbe3')
    parser.add_argument('bdf_filename', help='path to Nastran filename')
    parser.add_argument('-o', '--out', help='path to output Nastran filename')

    file_group = parser.add_mutually_exclusive_group(required=False)
    file_group.add_argument('--infile', help='defines list of RBE3s to update')
    # file_group.add_argument('--outfile', help='skip run the jobs')

    add_argparse_arguments(parser, ['--punch', '--lax'])
    # parser.add_argument('--test', action='store_false', help='skip run the jobs')
    # parser.add_argument('--debug', action='store_true', help='more debugging')
    # parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parser.add_argument('-v', '--version', action='version',
                        version=pyNastran.__version__)

    args = parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    bdf_filename = args.bdf_filename
    bdf_filename_out = args.out
    infile = args.infile
    print(f'infile = {infile!r}')

    from pyNastran.utils import print_bad_path
    log = SimpleLogger(level='debug')

    base, ext = os.path.splitext(bdf_filename)
    if bdf_filename_out is None:
        bdf_filename_out = f'{base}.out{ext}'
    assert args.punch is True, args
    # debug = args.debug

    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.mesh_utils.rbe_tools import rbe2_to_rbe3
    model = BDF(log=log)
    if args.lax:
        log.warning('using lax card parser')
        model.is_strict_card_parser = False
    model.read_bdf(bdf_filename, punch=args.punch, xref=False)

    if infile is None:
        eids_to_fix = list(model.rigid_elements)
    else:
        assert os.path.exists(infile), print_bad_path(infile)
        eids_to_fix = np.loadtxt(infile, dtype='int32').flatten().tolist()
    log.info(f'eids_to_fix = {eids_to_fix}')
    rbe2_to_rbe3(model, eids_to_fix)
    model.write_bdf(bdf_filename_out, write_header=False)
