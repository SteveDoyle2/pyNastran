import os
import sys
from cpylog import SimpleLogger

import pyNastran
from .utils import add_argparse_arguments, get_bdf_outfilename


def cmd_line_rbe3_to_rbe2(argv=None, quiet: bool=False) -> None:
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

    parser = argparse.ArgumentParser(prog='rbe3_to_rbe2')
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

    log = SimpleLogger(level='debug')
    bdf_filename_out = get_bdf_outfilename(
        bdf_filename, bdf_filename_out,
        tag='out')

    #assert args.punch is True, args
    # debug = args.debug

    from pyNastran.bdf.mesh_utils.rbe_tools import rbe3_to_rbe2
    from .utils_bdf import read_lax_bdf
    from .rbe_utils import load_ints_from_defaults

    model = read_lax_bdf(
        bdf_filename, punch=args.punch, xref=False,
        is_strict_card_parser=not args.lax,
        log=log)

    eids_to_fix = load_ints_from_defaults(
        model.rigid_elements, infile)
    log.info(f'eids_to_fix = {eids_to_fix}')

    rbe3_to_rbe2(model, eids_to_fix)
    model.write_bdf(bdf_filename_out, write_header=False)
