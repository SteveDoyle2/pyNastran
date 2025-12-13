from __future__ import annotations
import os
import sys
from cpylog import SimpleLogger
from .utils import add_argparse_arguments

# import pyNastran


def cmd_line_inclzip(argv=None, quiet: bool=False) -> None:
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
    parent_parser.add_argument('inclzip', type=str)
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
        bdf_filename_out = f'{bdf_filename_base}.zip{ext}'

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')

    from .utils_bdf import read_lax_bdf
    model = read_lax_bdf(
        bdf_filename, punch=args.punch, xref=False,
        validate=False,
        is_strict_card_parser=not args.lax,
        duplicate_cards=duplicate_cards,
        log=log)
    if bdf_filename_out:
        model.write_bdf(bdf_filename_out)
