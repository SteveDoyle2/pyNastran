import os
import sys
import argparse
from typing import Optional, Any
from cpylog import SimpleLogger


def get_bdf_filename_punch_log(data: dict[str, Any],
                               quiet: bool) -> tuple[str, bool, SimpleLogger]:
    """gets IN_BDF_FILENAME and --punch flag"""
    try:
        bdf_filename = data['IN_BDF_FILENAME']
    except:
        if not quiet:  # pragma: no cover
            print(data)
        raise
    punch = data['--punch'] if '--punch' in data else None

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    return bdf_filename, punch, log


def add_argparse_arguments(parser: argparse.ArgumentParser,
                           args: list[str]) -> None:
    # parent_parser.add_argument('--lax', action='store_false', help='lax card parser')
    for arg in args:
        if arg == '--punch':
            parser.add_argument('--punch', action='store_true', help='assume a punch file')
        elif arg == '--lax':
            parser.add_argument('--lax', action='store_true', help='lax card parser')
        elif arg == '--nosort':
            parser.add_argument(
                '--nosort', action='store_false',
                help='Dont sort the nodes, elements, ... (default=False -> sort)')
        elif arg == '--allow_dup':
            parser.add_argument('--allow_dup', help='allow duplicate cards -> "GRID,CONM2"')
        else:
            raise RuntimeError(arg)


def filter_no_args(msg: str, argv: list[str], quiet: bool=False):
    if len(argv) == 1:
        if quiet:
            sys.exit()
        sys.exit(msg)
    for i in range(len(argv)):
        arg = argv[i]
        if not isinstance(arg, str):
            argv[i] = str(arg)


def get_is_double_large(data: dict[str, Any]) -> tuple[int, bool]:
    is_double = False
    if data['double']:
        size = 16
        is_double = True
    elif data['large']:
        size = 16
    else:
        size = 8
    return size, is_double


def get_bdf_outfilename(bdf_filename: str,
                        bdf_filename_out: Optional[str]=None,
                        tag: str='out') -> str:
    base, ext = os.path.splitext(bdf_filename)
    if bdf_filename_out is None:
        bdf_filename_out = f'{base}.{tag}{ext}'
    return bdf_filename_out
