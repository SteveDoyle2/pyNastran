import sys
import argparse
from cpylog import SimpleLogger
import pyNastran
from .utils import filter_no_args


def cmd_line_flip_shell_normals(argv=None, quiet: bool=False) -> None:
    """command line interface to flip_shell_normals"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    filter_no_args("bdf flip_shell_normals: use 'bdf flip_shell_normals -h' for help",
                   argv, quiet=quiet)

    ver = str(pyNastran.__version__)

    parser = argparse.ArgumentParser(
        prog='bdf flip_shell_normals',
        description='Flip shell element normals',
    )
    parser.add_argument('-v', '--version', action='version', version=ver)
    parser.add_argument('IN_BDF_FILENAME',
                        help='path to input BDF/DAT/NAS file')
    parser.add_argument('-o', '--output', default='flipped_shell_normals.bdf',
                        metavar='OUT_BDF_FILENAME',
                        help='path to output BDF/DAT/NAS file (default=flipped_shell_normals.bdf)')
    parser.add_argument('--punch', action='store_true',
                        help='flag to identify a *.pch/*.inc file')
    parser.add_argument('--zero_zoffset', action='store_true',
                        help='zero out the z-offset')

    args = parser.parse_args(argv[2:])

    if not quiet:  # pragma: no cover
        print(vars(args))
    size = 16

    bdf_filename = args.IN_BDF_FILENAME
    punch = args.punch
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    zero_zoffset = args.zero_zoffset
    bdf_filename_out = args.output

    #from io import StringIO
    from pyNastran.bdf.bdf import read_bdf, BDF
    from pyNastran.bdf.mesh_utils.flip_shell_normals import flip_shell_normals

    model = BDF(log=log)
    model.set_error_storage(nparse_errors=100, stop_on_parsing_error=True,
                            nxref_errors=100, stop_on_xref_error=False)
    model = read_bdf(bdf_filename, punch=punch, log=log, xref=False)
    flip_shell_normals(model, zero_zoffset)
    model.write_bdf(bdf_filename_out, encoding=None,
                    size=size, nodes_size=16, elements_size=8, loads_size=8,
                    is_double=False, interspersed=False, enddata=None, write_header=True, close=True)
