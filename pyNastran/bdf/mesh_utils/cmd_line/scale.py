from __future__ import annotations
import os
import sys
from cpylog import SimpleLogger

import pyNastran


def cmd_line_scale(argv=None, quiet: bool=False) -> None:
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import argparse
    #import textwrap
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
    parent_parser.add_argument('scale', type=str)
    parent_parser.add_argument('INPUT', help='path to output BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('OUTPUT', nargs='?', help='path to output file', type=str)

    #'  --l  LENGTH_SF                    length scale factor\n'
    #'  --m  MASS_SF                      mass scale factor\n'
    #'  --f  FORCE_SF                     force scale factor\n'
    #'  --p  PRESSURE_SF                  pressure scale factor\n'
    #'  --t  TIME_SF                      time scale factor\n'
    #'  --v  VEL_SF                       velocity scale factor\n'

    parent_parser.add_argument('-l', '--length', help='length scale factor')
    parent_parser.add_argument('-m', '--mass', help='mass scale factor')
    parent_parser.add_argument('-f', '--force', help='force scale factor')
    parent_parser.add_argument('-p', '--pressure', help='pressure scale factor')
    parent_parser.add_argument('-t', '--time', help='time scale factor')
    parent_parser.add_argument('-V', '--velocity', help='velocity scale factor')
    #parent_parser.add_argument('--user_geom', type=str, help='log msg')

    #parent_parser.add_argument('-q', '--quiet', help='prints debug messages (default=True)', action='store_true')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parent_parser.add_argument('-v', '--version', action='version',
                               version=pyNastran.__version__)

    args = parent_parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    scales: list[float] = []
    terms: list[str] = []
    bdf_filename = args.INPUT
    bdf_filename_out = args.OUTPUT
    if bdf_filename_out is None:
        bdf_filename_base, ext = os.path.splitext(bdf_filename)
        bdf_filename_out = f'{bdf_filename_base}.scaled{ext}'

    #assert bdf_filename_out is not None
    if args.mass:
        scale = float(args.mass)
        scales.append(scale)
        terms.append('M')
    if args.length:
        scale = float(args.length)
        scales.append(scale)
        terms.append('L')
    if args.time:
        scale = float(args.time)
        scales.append(scale)
        terms.append('T')
    if args.force:
        scale = float(args.force)
        scales.append(scale)
        terms.append('F')
    if args.pressure:
        scale = float(args.pressure)
        scales.append(scale)
        terms.append('P')
    if args.velocity:
        scale = float(args.velocity)
        scales.append(scale)
        terms.append('V')

    from pyNastran.bdf.mesh_utils.convert import scale_by_terms

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    scale_by_terms(bdf_filename, terms, scales, bdf_filename_out=bdf_filename_out, log=log)
