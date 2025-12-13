import sys
from docopt import docopt
import pyNastran
from .utils import filter_no_args, get_bdf_filename_punch_log


def cmd_line_flip_shell_normals(argv=None, quiet: bool=False) -> None:
    """command line interface to flip_shell_normals"""
    if argv is None:  # pragma: no cover
        argv = sys.argv
    msg = (
        "Usage:\n"
        "  bdf flip_shell_normals IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--zero_zoffset]\n"
        "  bdf flip_shell_normals IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--zero_zoffset]\n"
        '  bdf flip_shell_normals -h | --help\n'
        '  bdf flip_shell_normals -v | --version\n'
        '\n'

        "Positional Arguments:\n"
        "  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n"
        #"  OUT_BDF_FILENAME   path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        "  -o OUT, --output  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
        '  --punch                             flag to identify a *.pch/*.inc file\n'
        "\n"  # (default=0.000001)

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
    size = 16

    bdf_filename, punch, log = get_bdf_filename_punch_log(data, quiet)
    zero_zoffset = data['--zero_zoffset']
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'flipped_shell_normals.bdf'

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
