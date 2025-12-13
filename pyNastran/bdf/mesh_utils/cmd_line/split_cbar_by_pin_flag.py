from __future__ import annotations
import sys
from docopt import docopt

import pyNastran
from .utils import get_bdf_filename_punch_log, filter_no_args


def cmd_line_split_cbars_by_pin_flag(argv=None, quiet: bool=False) -> None:
    """command line interface to split_cbars_by_pin_flag"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    msg = (
        'Usage:\n'
        '  bdf split_cbars_by_pin_flags  IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [-p PIN_FLAGS_CSV_FILENAME]\n'
        '  bdf split_cbars_by_pin_flags -h | --help\n'
        '  bdf split_cbars_by_pin_flags -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME    path to input BDF/DAT/NAS file\n'
        '\n'

        'Options:\n'
        ' -o OUT, --output  OUT_BDF_FILENAME        path to output BDF file\n'
        ' -p PIN, --pin     PIN_FLAGS_CSV_FILENAME  path to pin_flags_csv file\n'
        '  --punch                                  flag to identify a *.pch/*.inc file\n'
        '\n'

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
    #size = 16
    bdf_filename_in, punch, log = get_bdf_filename_punch_log(data, quiet)
    bdf_filename_out = data['--output']
    if bdf_filename_out is None:
        bdf_filename_out = 'model_new.bdf'

    pin_flags_filename = data['--pin']
    if pin_flags_filename is None:
        pin_flags_filename = 'pin_flags.csv'

    from pyNastran.bdf.mesh_utils.split_cbars_by_pin_flag import split_cbars_by_pin_flag
    split_cbars_by_pin_flag(bdf_filename_in, pin_flags_filename=pin_flags_filename,
                            bdf_filename_out=bdf_filename_out,
                            punch=punch)
