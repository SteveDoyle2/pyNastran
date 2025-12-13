from __future__ import annotations
import sys
from cpylog import SimpleLogger

import pyNastran
from .utils import filter_no_args, get_is_double_large


def cmd_line_free_faces(argv=None, quiet: bool=False) -> None:
    """command line interface to bdf free_faces"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    encoding = sys.getdefaultencoding()
    usage = (
        'Usage:\n'
        '  bdf free_faces BDF_FILENAME SKIN_FILENAME [-d] [-l] [-f] [--encoding ENCODE]\n'
        '  bdf free_faces -h | --help\n'
        '  bdf free_faces -v | --version\n'
        '\n'
    )
    arg_msg = (
        "Positional Arguments:\n"
        "  BDF_FILENAME    path to input BDF/DAT/NAS file\n"
        "  SKIN_FILENAME   path to output BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        '  -l, --large        writes the BDF in large field, single precision format (default=False)\n'
        '  -d, --double       writes the BDF in large field, double precision format (default=False)\n'
        f'  --encoding ENCODE  the encoding method (default=None -> {encoding!r})\n'
        '\n'
        'Developer:\n'
        '  -f, --profile    Profiles the code (default=False)\n'
        '\n'
        "Info:\n"
        '  -h, --help     show this help message and exit\n'
        "  -v, --version  show program's version number and exit\n"
    )
    filter_no_args(arg_msg, argv, quiet=quiet)

    arg_msg += '\n'

    examples = (
        'Examples\n'
        '--------\n'
        '  bdf free_faces solid.bdf skin.bdf\n'
        '  bdf free_faces solid.bdf skin.bdf --large\n'
    )
    import argparse
    parent_parser = argparse.ArgumentParser()
    # positional arguments
    parent_parser.add_argument('BDF_FILENAME', help='path to input BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('SKIN_FILENAME', help='path to output BDF/DAT/NAS file', type=str)

    size_group = parent_parser.add_mutually_exclusive_group()
    size_group.add_argument('-d', '--double', help='writes the BDF in large field, single precision format', action='store_true')
    size_group.add_argument('-l', '--large', help='writes the BDF in large field, double precision format', action='store_true')
    size_group.add_argument('--encoding', help=f'the encoding method (default=None -> {repr(encoding)})', type=str)
    parent_parser.add_argument('--profile', help='Profiles the code', action='store_true')
    parent_parser.add_argument('-v', '--version', action='version', version=pyNastran.__version__)

    from pyNastran.utils.arg_handling import argparse_to_dict, update_message

    update_message(parent_parser, usage, arg_msg, examples)
    if not quiet:  # pragma: no cover
        print(argv)
    args = parent_parser.parse_args(args=argv[2:])
    data = argparse_to_dict(args)

    if not quiet:  # pragma: no cover
        for key, value in sorted(data.items()):
            print("%-12s = %r" % (key.strip('--'), value))

    import time
    time0 = time.time()

    size, is_double = get_is_double_large(data)
    bdf_filename = data['BDF_FILENAME']
    skin_filename = data['SKIN_FILENAME']

    from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes
    from pyNastran.bdf.mesh_utils.free_faces import write_skin_solid_faces
    tol = 1e-005
    bdf_filename_merged = 'merged.bdf'
    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    bdf_equivalence_nodes(bdf_filename, bdf_filename_merged, tol,
                          renumber_nodes=False, neq_max=10, xref=True,
                          node_set=None,
                          size=8, is_double=is_double,
                          remove_collapsed_elements=False,
                          avoid_collapsed_elements=False,
                          crash_on_collapse=False, log=log, debug=True)
    if not quiet:  # pragma: no cover
        print('done with equivalencing')
    write_skin_solid_faces(
        bdf_filename_merged, skin_filename,
        write_solids=False, write_shells=True,
        size=size, is_double=is_double, encoding=None, log=log,
    )
    if not quiet:  # pragma: no cover
        print('total time:  %.2f sec' % (time.time() - time0))
