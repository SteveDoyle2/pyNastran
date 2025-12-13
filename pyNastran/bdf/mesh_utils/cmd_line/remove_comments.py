from __future__ import annotations
import os
import sys
from typing import TYPE_CHECKING
from cpylog import SimpleLogger

import pyNastran
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def cmd_line_remove_comments(argv=None, quiet: bool=False) -> BDF:
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
    parent_parser.add_argument('remove_comments', type=str)
    parent_parser.add_argument('INPUT', help='path to output BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('OUTPUT', nargs='?', help='path to output file', type=str)
    #parent_parser.add_argument('--noxref', action='store_true', help='skips cross-referencing (default=True)')

    parent_parser.add_argument('-q', '--quiet', action='store_true', help='prints debug messages (default=True)')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parent_parser.add_argument('-v', '--version', action='version',
                               version=pyNastran.__version__)

    args = parent_parser.parse_args(args=argv[1:])
    if args.quiet:
        quiet = args.quiet
    if not quiet:  # pragma: no cover
        print(args)

    #xref = not args.noxref
    bdf_filename = args.INPUT
    bdf_filename_out = args.OUTPUT
    if bdf_filename_out is None:
        bdf_filename_base, ext = os.path.splitext(bdf_filename)
        bdf_filename_out = '%s.nocomments%s' % (bdf_filename_base, ext)

    from pyNastran.bdf.mesh_utils.bdf_remove_comments import bdf_remove_comments

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    model = bdf_remove_comments(
        bdf_filename, bdf_filename_out=bdf_filename_out,
        #xref=xref,
        log=log)
    return model
