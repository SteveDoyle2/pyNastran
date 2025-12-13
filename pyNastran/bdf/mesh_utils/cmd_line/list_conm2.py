from __future__ import annotations
import sys
from cpylog import SimpleLogger

import pyNastran
from .utils import filter_no_args, add_argparse_arguments


def cmd_line_list_conm2(argv=None, quiet=False) -> None:
    """command line interface to bdf list_conm2"""
    if argv is None:  # pragma: no cover
        argv = sys.argv
    encoding = sys.getdefaultencoding()
    usage = (
        'Usage:\n'
        '  bdf list_conm2 BDF_FILENAME [--scale SCALE] [--encoding ENCODE]\n'
        '  bdf list_conm2 -h | --help\n'
        '  bdf list_conm2 -v | --version\n'
        '\n'
    )
    arg_msg = (
        "Positional Arguments:\n"
        "  BDF_FILENAME    path to input BDF/DAT/NAS file\n"
        '\n'

        'Options:\n'
        f'  --encoding ENCODE  the encoding method (default=None -> {encoding!r})\n'
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
        '  bdf list_conm2 fem.bdf --scale 386.1\n'
    )
    import argparse
    parent_parser = argparse.ArgumentParser()
    # positional arguments
    parent_parser.add_argument('BDF_FILENAME', help='path to input BDF/DAT/NAS file', type=str)
    add_argparse_arguments(parent_parser, ['--lax'])

    size_group = parent_parser.add_mutually_exclusive_group()
    size_group.add_argument('--scale', help='scales the mass')
    #size_group.add_argument('--encoding', help=f'the encoding method (default=None -> {repr(encoding)})', type=str)
    parent_parser.add_argument('-v', '--version', action='version', version=pyNastran.__version__)

    from pyNastran.utils.arg_handling import argparse_to_dict, update_message

    update_message(parent_parser, usage, arg_msg, examples)
    if not quiet:
        print(argv)
    args = parent_parser.parse_args(args=argv[2:])
    data = argparse_to_dict(args)

    if not quiet:  # pragma: no cover
        for key, value in sorted(data.items()):
            print("%-12s = %r" % (key.strip('--'), value))

    # import time
    # time0 = time.time()

    #size, is_double = get_is_double_large(data)
    print(data)
    bdf_filename = data['BDF_FILENAME']
    mass_scale = 1.
    if 'scale' in data and data['scale'] is not None:
        mass_scale = float(data['scale'])
    print(f'mass_scale = {mass_scale}')

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')

    from pyNastran.bdf.mesh_utils.list_conm2 import list_conm2
    list_conm2(bdf_filename, mass_scale, log=log)
