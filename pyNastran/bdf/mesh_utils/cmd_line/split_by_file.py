"""
Given a vehicle with many properties, it's likely that
some are going be duplicated. However, these property regions
do not correspond to logical design regions (e.g., airplane
wing vs. fuselage).

This takes the vehicle and the sub-regions (e.g., wing,
fuselage) and renumbers the properties.

bdf split_by_file bwb_saero.bdf eids.csv --obj --skip_props PSOLID
"""
from __future__ import annotations
import os
import sys
import copy
from collections import defaultdict
from typing import Optional, TYPE_CHECKING
import numpy as np

from cpylog import SimpleLogger
import pyNastran
from .utils import add_argparse_arguments
from pyNastran.utils import PathLike, print_bad_path

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def cmd_line_split_by_file(argv=None, quiet: bool=False,
                           ) -> BDF:
    """command line interface to split_by_file"""
    if argv is None:  # pragma: no cover
        argv = sys.argv

    import argparse
    from pyNastran.utils.arg_handling import argparse_to_dict  # update_message
    parent_parser = argparse.ArgumentParser()
    # positional arguments
    parent_parser.add_argument('IN_BDF_FILENAME', help='path to input BDF/DAT/NAS file', type=str)
    parent_parser.add_argument('EID_CSV_FILES', help='path to input CSV files', type=str, nargs='+')
    parent_parser.add_argument('--skip_props', help='CSV list of properties that shouldnt be considered', type=str)
    parent_parser.add_argument('-o', '--output', help='path to output BDF/DAT/NAS file', type=str)
    add_argparse_arguments(parent_parser, ['--obj', '--punch', '--lax', '--allow_dup'])
    parent_parser.add_argument('--nocheck', help='read the bdf', action='store_true')
    parent_parser.add_argument('-v', '--version', action='version', version=pyNastran.__version__)
    args = parent_parser.parse_args(args=argv[2:])
    data = argparse_to_dict(args)

    if not quiet:  # pragma: no cover
        print(data)
    #size = 16

    punch = args.punch
    is_obj = args.obj
    is_strict_card_parser = not args.lax
    level = 'debug' if not quiet else 'warning'

    bdf_filename = data['IN_BDF_FILENAME']
    eid_filenames = data['EID_CSV_FILES']
    bdf_filename_out = data['output']
    properties_to_skip = []
    if data['skip_props']:
        properties_to_skip = [
            prop.strip() for prop in
            data['skip_props'].strip(', ').upper().split(',')]

    base, ext = os.path.splitext(bdf_filename)
    obj_filename = base + '.obj'
    if bdf_filename_out is None:
        bdf_filename_out = base + f'.split_file{ext}'

    assert len(eid_filenames) > 0, eid_filenames
    log = SimpleLogger(level=level, encoding='utf-8')
    _check_missing_files(eid_filenames, log)

    from .utils_bdf import read_lax_obj
    model = read_lax_obj(bdf_filename,
        obj_filename, is_obj,
        is_strict_card_parser=is_strict_card_parser,
        xref=False, punch=punch,
        duplicate_cards=None,
        log=log)

    # {'PSOLID', 'PBUSH'}
    split_by_file(model, eid_filenames,
                  properties_to_skip=properties_to_skip,
                  bdf_filename_out=bdf_filename_out)
    log.info('done')

    if not args.nocheck:
        from pyNastran.bdf.bdf import read_bdf
        read_bdf(bdf_filename_out, log=log)
    return model


def _check_missing_files(eid_filenames: list[str],
                         log: SimpleLogger) -> None:
    is_missing = False
    for eid_filename in eid_filenames:
        if not os.path.exists(eid_filename):
            is_missing = True
            log.error(print_bad_path(eid_filename))
    if is_missing:
        raise RuntimeError('There are missing files')


def split_by_file(model: BDF,
                  eid_filenames: list[PathLike | np.ndarray],
                  properties_to_skip: Optional[list[str] | set[str]]=None,
                  bdf_filename_out: PathLike='') -> BDF:
    """
    splits a model
    assumes that model is not cross-referenced
    doesn't handle optimization
    """
    if properties_to_skip is None:
        properties_to_skip = set([])
    elif isinstance(properties_to_skip, (list, tuple)):
        properties_to_skip = set(properties_to_skip)
    assert isinstance(properties_to_skip, set), properties_to_skip

    log: SimpleLogger = model.log
    assert len(eid_filenames) > 0, eid_filenames

    # find the elements in each rgion to renumber
    eid_sets = load_eid_sets(eid_filenames, log)
    assert len(eid_sets) == len(eid_filenames)

    eid_sets = get_unique_eid_sets(eid_filenames, eid_sets, log)
    assert len(eid_sets) == len(eid_filenames)

    eids_all = list(model.elements)
    eids_set_all = np.hstack(eid_sets)
    eids_diff = np.setdiff1d(eids_set_all, eids_all)
    if len(eids_diff):
        model.log.warning(f'eids_diff = {eids_diff}')
    eids_intersect = np.intersect1d(eids_all, eids_set_all)

    eid_sets_intersect = get_intersected_sets(
        eids_intersect, eid_sets)

    pid_to_pidnew = _update_elements(
        model, properties_to_skip, eid_sets_intersect, log)

    log.debug(f'pid_to_pidnew = {pid_to_pidnew}')

    if bdf_filename_out:
        model.write_bdf(bdf_filename_out)
    return model


def _update_elements(model: BDF,
                     properties_to_skip: set[str],
                     eid_sets_intersect: dict[int, np.ndarray],
                     log: SimpleLogger) -> dict[tuple[int, int], int]:
    """
    Loop over each unique element set.
    Get the list of all elements with each property id.
    Then update the property id and remap those element pids.
    """
    pid_new = max(model.properties) + 1
    log.debug(f'pid_new0 = {pid_new}')

    pid_to_pidnew = {}
    for ifile, eids_intersect in eid_sets_intersect.items():
        log.debug(f'ifile={ifile} eids={eids_intersect} pid_new={pid_new}')
        # get pids
        pid_to_elems = defaultdict(list)
        for eid in eids_intersect:
            elem = model.elements[eid]
            pid = elem.pid
            prop = model.properties[pid]
            if prop.type in properties_to_skip:
                log.info(f'skipping eid={eid} pid={pid}; found prop_type={prop.type}')
                continue
            pid_to_elems[pid].append(elem)

        # update pid to pid_new; update element
        for pid, elems in pid_to_elems.items():
            pid_to_pidnew[(ifile, pid)] = pid_new
            prop = model.properties[pid]
            prop2 = copy.deepcopy(prop)
            prop2.pid = pid_new
            model.properties[pid_new] = prop2
            for elem in elems:
                elem.pid = pid_new
    return pid_to_pidnew


def get_intersected_sets(eids_all: np.ndarray,
                         eid_sets: list[np.ndarray]) -> dict[int, np.ndarray]:
    """make sure that all elements exist"""
    eid_sets_intersect = {}
    for i, eids in enumerate(eid_sets):
        eids_intersect = np.intersect1d(eids_all, eids)
        eids_all = np.setdiff1d(eids_all, eids)
        if len(eids_intersect):
            eid_sets_intersect[i] = eids_intersect
    return eid_sets_intersect


def load_eid_sets(eid_filenames: list[PathLike | np.ndarray],
                  log: SimpleLogger) -> list[np.ndarray]:
    """read in the data"""
    eid_sets = []
    for eid_filename in eid_filenames:
        if isinstance(eid_filename, PathLike):
            log.debug(f'loading {eid_filename}')
            assert os.path.exists(eid_filename), print_bad_path(eid_filename)
            eids = np.loadtxt(eid_filename).astype('int32')
        else:
            # assert isinstance(eid_filename, np.ndarray), eid_filename
            eids = np.asarray(eid_filename, dtype='int32')
        eid_sets.append(eids)
    return eid_sets


def get_unique_eid_sets(
        eid_filenames: list[PathLike | np.ndarray],
        eids_sets: list[np.ndarray],
        log: SimpleLogger) -> list[np.ndarray]:
    """elements get put in the first set that defines them"""
    eid_sets2 = [eids_sets[0]]
    eids_all = eids_sets[0]
    for i, eids in enumerate(eids_sets[1:]):
        diff = np.setdiff1d(eids_all, eids)
        intersect = np.intersect1d(eids_all, eids)
        if len(intersect):
            eid_filename = eid_filenames[i+1]
            if isinstance(eid_filename, PathLike):
                eid_filename = eid_filenames[i+1]
            else:
                eid_filename = f'eids={i+1}'
            log.warning(f'removing duplicate elements in {eid_filename}; eids={intersect}')
        eids_all = np.union1d(eids_all, eids)
        eid_sets2.append(diff)
    return eid_sets2
