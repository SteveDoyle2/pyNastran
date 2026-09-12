import sys
import argparse


def cmd_line_merge(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf_merge"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import pyNastran

    if len(argv) == 1:
        sys.exit("bdf merge: use 'bdf merge -h' for help")

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf merge',
        description='Merge multiple BDF files into one',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAMES', nargs='+',
                        help='path to input BDF/DAT/NAS files')
    parser.add_argument('-o', '--output', default='merged.bdf',
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF/DAT/NAS file (default=merged.bdf)')
    parser.add_argument('--debug', action='store_true',
                        help='enable debug output')

    args = parser.parse_args(argv[2:])
    if not quiet:  # pragma: no cover
        print(vars(args))

    size = 16
    bdf_filenames = args.IN_BDF_FILENAMES
    bdf_filename_out = args.output
    debug = args.debug

    #cards_to_skip = [
        #'AEFACT', 'CAERO1', 'CAERO2', 'SPLINE1', 'SPLINE2',
        #'AERO', 'AEROS', 'PAERO1', 'PAERO2', 'MKAERO1']
    cards_to_skip = []

    from cpylog import SimpleLogger
    from pyNastran.bdf.mesh_utils.bdf_merge import bdf_merge
    log = SimpleLogger(level='debug')

    bdf_merge(bdf_filenames, bdf_filename_out, renumber=True,
              encoding=None, size=size, is_double=False, cards_to_skip=cards_to_skip,
              log=log)
