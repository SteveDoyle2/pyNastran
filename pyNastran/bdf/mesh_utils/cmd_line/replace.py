import sys
import os
from itertools import count

from cpylog import SimpleLogger
from .utils import add_argparse_arguments
import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.bdf.mesh_utils.bdf_renumber import _get_bdf_model


def cmd_line_replace(argv=None, quiet: bool=False):
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
    parent_parser.add_argument('replace', type=str)
    parent_parser.add_argument('bdf_main', help='path to output BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('blk_files', nargs='+', help='path to output BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('-o', '--out', nargs='?', help='path to output file', type=str, default='')
    parent_parser.add_argument('--blk_punch', nargs='?', help='path to output file', type=str, default='')

    parent_parser.add_argument('-q', '--quiet', action='store_true', help='prints debug messages (default=True)')
    #parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parent_parser.add_argument('-v', '--version', action='version',
                               version=pyNastran.__version__)
    add_argparse_arguments(parent_parser, ['--obj', '--punch', '--lax', '--allow_dup'])

    args = parent_parser.parse_args(args=argv[1:])
    if args.quiet:
        quiet = args.quiet
    if not quiet:  # pragma: no cover
        print(args)

    #xref = not args.noxref
    punch = args.punch
    blk_punch = args.blk_punch.strip()
    bdf_filename = args.bdf_main
    replace_filenames = args.blk_files
    assert isinstance(replace_filenames, list), replace_filenames
    assert len(replace_filenames) > 0, replace_filenames

    assert isinstance(blk_punch, str), blk_punch
    if len(blk_punch) == 0:
        punch_replace = [True] * len(replace_filenames)
    else:
        punch_replace = [punchi in 'tT' for punchi in blk_punch]
        assert len(blk_punch) == len(replace_filenames), replace_filenames

    input_filenames = [bdf_filename] + replace_filenames
    for input_filename in input_filenames:
        assert os.path.exists(input_filename), print_bad_path(input_filename)

    bdf_filename_out = args.out
    if bdf_filename_out == '':
        bdf_filename_base, ext = os.path.splitext(bdf_filename)
        bdf_filename_out = '%s.replace%s' % (bdf_filename_base, ext)
    # print(f'bdf_filename_out = {bdf_filename_out}')

    level = 'debug' if not quiet else 'warning'
    log = SimpleLogger(level=level, encoding='utf-8')
    model = bdf_replace(
        bdf_filename, replace_filenames, bdf_filename_out=bdf_filename_out,
        #xref=xref,
        punch=punch,
        punch_replace=punch_replace,
        log=log)
    return model


def bdf_replace(bdf_filename, replace_filenames: list[str],
                bdf_filename_out=None,
                punch: bool=True,
                punch_replace: bool=True,
                log=None, cards_to_skip=None):
    if isinstance(punch_replace, bool):
        punch_replace = [punch_replace] * len(replace_filenames)

    model = _get_bdf_model(bdf_filename, punch=punch,
                           xref=False,
                           cards_to_skip=cards_to_skip,
                           log=log, debug=True)
    log = model.log
    for ifile, replace_filename, punchi in zip(count(), replace_filenames, punch_replace):
        modeli = _get_bdf_model(replace_filename, punch=punchi,
                                xref=False,
                                cards_to_skip=cards_to_skip,
                                log=log, debug=True)
        log.info(''.join(modeli.get_bdf_stats()))
        model.replace_cards(modeli)

    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out)
    return model
