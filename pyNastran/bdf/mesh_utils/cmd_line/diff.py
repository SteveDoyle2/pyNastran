from __future__ import annotations
import sys


def cmd_line_apply_diff(argv=None, quiet: bool=False) -> None:
    if argv is None:  # pragma: no cover
        argv = sys.argv[1:]
    print(f'argv = {argv}')

    import pyNastran
    import argparse
    from pyNastran.bdf.mesh_utils.bdf_diff_remove import apply_diff

    parser = argparse.ArgumentParser()
    parser.add_argument("diff", type=str)
    parser.add_argument("bdf_filename", help='path to main BDF/DAT/NAS file', type=str)
    parser.add_argument("diff_filename", help='path to diff BDF/DAT/NAS file', type=str)
    args = parser.parse_args()
    print(f'args = {args}')
    bdf_filename = args.bdf_filename
    diff_filename = args.diff_filename
    apply_diff(bdf_filename, diff_filename)


def cmd_line_diff(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf_diff"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import argparse
    import pyNastran

    if len(argv) == 1:
        sys.exit("bdf diff: use 'bdf diff -h' for help")

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf diff',
        description='Diff two BDF files',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME1',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('IN_BDF_FILENAME2',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('--punch', action='store_true',
                        help='uses a punch file')
    parser.add_argument('--skip_cards', default=None, metavar='CARDS',
                        help='comma-separated list of cards to skip')

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    bdf_filename1 = args.IN_BDF_FILENAME1
    bdf_filename2 = args.IN_BDF_FILENAME2
    skip_cards = []
    if args.skip_cards:
        skip_cards = args.skip_cards.split(',')

    debug = None if quiet else True

    from cpylog import SimpleLogger
    from pyNastran.bdf.mesh_utils.bdf_diff import get_diff_bdfs
    level = 'warning' if debug is None else 'debug' if debug else 'info'
    log = SimpleLogger(level=level)
    print(log)
    added_cards, removed_cards, added_model, removed_model = get_diff_bdfs(
        bdf_filename1, bdf_filename2, skip_cards=skip_cards, log=log)
