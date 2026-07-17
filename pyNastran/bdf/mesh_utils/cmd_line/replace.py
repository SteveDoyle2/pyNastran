from __future__ import annotations
import sys
import os
from itertools import count
from typing import TYPE_CHECKING

from cpylog import SimpleLogger
from .utils import add_argparse_arguments
import pyNastran
from pyNastran.utils import print_bad_path
from pyNastran.bdf.mesh_utils.bdf_renumber import _get_bdf_model
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF


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
    parent_parser.add_argument('-o', '--out', nargs='?', help='path to output file',
                               type=str, default='')
    parent_parser.add_argument('--blk_punch', nargs='?', help='make blk_files punch files (TF for true, false for 2 files)',
                               type=str, default='')

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
                punch_replace: bool | list[bool]=True,
                log=None, cards_to_skip=None,
                write_log: bool=True):
    if isinstance(punch_replace, bool):
        punch_replace = [punch_replace] * len(replace_filenames)

    log.info(f'bdf_filename: {bdf_filename}')
    model = _get_bdf_model(bdf_filename, punch=punch,
                           xref=False,
                           cards_to_skip=cards_to_skip,
                           log=log, debug=True)
    log = model.log
    for ifile, replace_filename, punchi in zip(count(), replace_filenames, punch_replace):
        log.info(f'replace_filename: {replace_filename}')
        modeli = _get_bdf_model(replace_filename, punch=punchi,
                                xref=False,
                                cards_to_skip=cards_to_skip,
                                log=log, debug=True)
        log.info(''.join(modeli.get_bdf_stats()))
        model.replace_cards(modeli, write_log=write_log)

    if bdf_filename_out is not None:
        model.write_bdf(bdf_filename_out)
    return model


DICT_NAMES = [
    'nodes', 'elements', 'properties', 'masses', 'rigid_elements',
    'materials', 'thermal_materials',
    'tables', 'tables_d', 'tables_m', 'tables_sdamping', 'random_tables',
    'trims', 'csschds', 'flfacts', 'flutters', 'gusts',
    # 'transfer_functions'
    'desvars', 'dvprels', 'dvmrels', 'dvgrids',
    'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9', 'MATT11',
    'MATS1', 'MATS8',
    'params', 'methods', 'tsteps', 'tstepnls',
]
DICT_LIST_NAMES = [
    'spcadds', 'mpcadds', 'dconadds',
    'spcs', 'loads', 'dconstrs',
]

def replace_cards(self: BDF,
                  replace_model: BDF,
                  write_log: bool = True) -> None:
    """
    Replaces the common cards from the current (self) model from the
    ones in the new replace_model.  The intention is that you're
    going to replace things like PSHELLs and DESVARs from a pch file
    in order to update your BDF with the optimized geometry. You can
    also just add cards with this and if the ids exist, it'll overwrite.

    .. todo:: only does a subset of cards.

    Notes
    -----
    loads/spcs (not supported) are tricky because you can't replace
    cards one-to-one...not sure what to do

    """
    log = self.log
    log.info('replacing cards')
    log.info(replace_model.get_bdf_stats())
    # singeltons = ['aero', 'aeros']
    # assert len(self.nodes) > 0
    # assert len(replace_model.nodes) > 0
    ids_changed_dict = {}
    if write_log:
        for name in DICT_NAMES:
            # we're basically just doing this
            # for nid, node in replace_model.nodes.items():
            #     self.nodes[nid] = node

            my_objs = getattr(self, name)
            replace_objs = getattr(replace_model, name)
            assert isinstance(replace_objs, dict), (name, type(replace_objs))
            ids_changed = []
            for idi, replace_obj in replace_objs.items():
                assert not isinstance(replace_obj, list), replace_obj
                if idi not in my_objs:
                    log.debug(f'adding:\n{str(replace_obj)}\n'
                              '----------------------------')
                    my_objs[idi] = replace_obj
                    ids_changed.append(idi)
                    continue

                my_obj = my_objs[idi]
                if my_obj != replace_obj:
                    log.debug(f'replacing:\n{str(my_obj)}\nwith:\n{str(replace_obj)}\n'
                              '----------------------------')
                    my_objs[idi] = replace_obj
                    ids_changed.append(idi)
                # else:
                #     print(f'{name} = {idi}')

            if len(ids_changed):
                ids_changed_dict[name] = ids_changed

        for name in DICT_LIST_NAMES:
            my_objs = getattr(self, name)
            replace_objs_dict = getattr(replace_model, name)
            assert isinstance(replace_objs_dict, dict), (name, type(replace_objs_dict))
            for idi, replace_objs in replace_objs_dict.items():
                assert isinstance(replace_objs, list), (name, replace_objs)
                ## TODO: how sould I replace loads or spcs?
                ##       just whole sale replace
                log.debug(f'replacing {name}={idi}:')
                for obj in replace_objs:
                    log.debug(str(obj))
                log.debug('----------')
                log.debug(f'with {name}={idi}:')
                for obj in replace_objs:
                    log.debug(str(obj))
                log.debug('----------------------------')
                my_objs[idi] = replace_objs

    else:
        for name in DICT_NAMES:
            my_objs = getattr(self, name)
            dict_objs = getattr(replace_model, name)
            for idi, obj in dict_objs.items():
                my_objs[idi] = obj

        for name in DICT_LIST_NAMES:
            my_objs = getattr(self, name)
            dict_objs = getattr(replace_model, name)
            for idi, objs in dict_objs.items():
                assert isinstance(objs, list), (name, objs)
                ## TODO: how sould I replace loads or spcs?
                ##       just whole sale replace
                my_objs[idi] = objs

    for name, ids_changed in ids_changed_dict.items():
        log.info(f'{name}_ids_changed = {ids_changed}')

    for cid, coord in replace_model.coords.items():
        if cid == 0:
            continue
        self.coords[cid] = coord
