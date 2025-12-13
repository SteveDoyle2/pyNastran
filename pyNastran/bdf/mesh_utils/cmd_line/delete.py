from __future__ import annotations
import os
import sys
from typing import TYPE_CHECKING
from cpylog import SimpleLogger

from .utils import add_argparse_arguments
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def cmd_line_delete(argv=None, quiet: bool=False) -> None:
    """command line interface to ``bdf delete``"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import argparse
    parent_parser = argparse.ArgumentParser(
        #prog = 'pyNastranGUI',
        #usage = usage,
        #description='A foo that bars',
        #epilog="And that's how you'd foo a bar",
        #formatter_class=argparse.RawDescriptionHelpFormatter,
        #description=textwrap.dedent(text),
        #version=pyNastran.__version__,
        #add_help=False,
    )
    # positional arguments
    parent_parser.add_argument('delete', type=str)
    parent_parser.add_argument('INPUT', help='path to output BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('OUTPUT', nargs='?', help='path to output file', type=str)
    add_argparse_arguments(parent_parser, ['--punch', '--lax', '--allow_dup'])
    args = parent_parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    bdf_filename = args.INPUT
    bdf_filename_out = args.OUTPUT
    is_strict_card_parser = not args.lax
    punch = args.punch
    duplicate_cards = args.allow_dup.split(',') if args.allow_dup else []

    bdf_filename_out = ''
    if bdf_filename_out is None:
        bdf_filename_base, ext = os.path.splitext(bdf_filename)
        bdf_filename_out = f'{bdf_filename_base}.delete{ext}'

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')

    model = BDF(log=log)
    bdf_filename = 'fem.bdf'
    #encoding = None
    read_includes = False
    cards_list = model.read_cards(
        bdf_filename=bdf_filename,
        punch=punch, read_includes=read_includes,
        save_file_structure=False,
        #encoding=encoding,
    )

    is_list = False
    has_none = True
    for card in cards_list:
        card_name, comment, card_lines, (ifile, unused_iline) = card
        card_name = card_name.upper()
        card_obj, unused_card = model.create_card_object(
            card_lines, card_name,
            is_list=is_list, has_none=has_none)
        # regex to remove elements, properties, ...
