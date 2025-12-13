import os
import sys
from docopt import docopt

import pyNastran
from .utils import (
    get_bdf_filename_punch_log, filter_no_args,
)


def cmd_line_collapse_quads(argv=None, quiet: bool=False) -> None:
    """command line interface to ``delete_bad_shells``"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    msg = (
        'Usage:\n'
        '  bdf collapse_quads IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--punch] [--size SIZE]\n'

        '  bdf collapse_quads -h | --help\n'
        '  bdf collapse_quads -v | --version\n'
        '\n'

        'Positional Arguments:\n'
        '  IN_BDF_FILENAME        path to input BDF/DAT/NAS file\n'
        #"  OUT_BDF_FILENAME  path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        '  -o     OUT   path to output BDF/DAT/NAS file\n'
        '  --size SIZE  size of the output\n'
        '  --punch      flag to identify a *.pch/*.inc file\n'

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

    bdf_filename, punch, log = get_bdf_filename_punch_log(data, quiet)

    base, ext = os.path.splitext(bdf_filename)
    bdf_filename_out = base + '_collapsed' + ext
    if data['OUT_BDF_FILENAME']:
        bdf_filename_out = data['OUT_BDF_FILENAME']

    size = 8
    if data['--size']:
        size_str = data['--size']
        size = int(size_str)

    from pyNastran.bdf.mesh_utils.collapse_bad_quads import convert_bad_quads_to_tris
    from pyNastran.bdf.bdf import read_bdf, BDF
    model: BDF = read_bdf(bdf_filename, xref=False,
                          validate=False, punch=punch, log=log)
    convert_bad_quads_to_tris(model)
    model.write_bdf(bdf_filename_out, size=size)
