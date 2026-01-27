from __future__ import annotations
import os
import sys
import argparse
from collections import defaultdict

import numpy as np
import pyNastran
from .utils import add_argparse_arguments

from pyNastran.bdf.bdf import BDF
from cpylog import SimpleLogger

SHELLS = {'CTRIA3', 'CTRIAR', 'CTRIA6',
          'CQUAD4', 'CQUADR', 'CQUAD', 'CQUAD8'}
LINES = {'CBAR', 'CBEAM'}
ETYPES_NOT_CONSIDERED = {
    'CBUSH', 'CHEXA', 'CTETRA', 'CPENTA',
    'CROD', 'CONROD', 'CTUBE', 'CBEAM3',
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
    'CBEND', 'CBUSH1D',
}

def float_eval(value_str: str):
    """casts a float"""
    try:
        value = float(value_str)
    except ValueError:
        value = eval(value_str)
    return float(value)

def cmd_line_nsm_split(argv=None, quiet: bool=False) -> None:
    """
    nsm_split filename.bdf [nsm_id nsm_value]...

    nsm_split model.bdf 100 '10/(12*32.174)' 20  '20/(12*32.174)'
    """
    FILE = os.path.abspath(__file__)
    if argv is None:
        argv = sys.argv[1:]  # ['run_jobs'] + sys.argv[2:]
    else:
        argv = [FILE] + argv[2:]  # ['run_jobs'] + sys.argv[2:]
    if not quiet:
        print(f'argv = {argv}')
    # print(f'argv = {argv}')

    parser = argparse.ArgumentParser(prog='nsm_split')
    parser.add_argument('bdf_filename', help='path to Nastran filename')
    # parser.add_argument('nsm_id', help='nsm id to split')
    # parser.add_argument('nsm_value', help='nsm value to apply')

    parser.add_argument(
        'nsm_value',
        metavar='nsm_id nsm_value',
        nargs='*',
        help='nsm id to split',)
    parser.add_argument('-o', '--out', help='path to output Nastran filename')

    # file_group = parser.add_mutually_exclusive_group(required=False)
    # file_group.add_argument('--infile', help='defines list of RBE2s to update')
    # file_group.add_argument('--outfile', help='skip run the jobs')

    add_argparse_arguments(parser, ['--punch', '--lax'])
    # parser.add_argument('--test', action='store_false', help='skip run the jobs')
    # parser.add_argument('--debug', action='store_true', help='more debugging')
    # parent_parser.add_argument('-h', '--help', help='show this help message and exits', action='store_true')
    parser.add_argument('-v', '--version', action='version',
                        version=pyNastran.__version__)

    args = parser.parse_args(args=argv[1:])
    if not quiet:  # pragma: no cover
        print(args)

    bdf_filename = args.bdf_filename

    assert len(args.nsm_value) > 1, args.nsm_value
    assert len(args.nsm_value) % 2 == 0, args.nsm_value
    nsm_ids = []
    nsm_values = []
    nsm_strings = []
    for i, value_str in enumerate(args.nsm_value):
        if i % 2 == 0:
            nsm_id = int(value_str)
            nsm_ids.append(nsm_id)
        else:
            nsm_strings.append(value_str)
            nsm_value = float_eval(value_str)
            nsm_values.append(nsm_value)

    log = SimpleLogger(level='debug')
    log.info(f'nsm_ids = {nsm_ids}')
    log.info(f'nsm_values = {nsm_values}')
    bdf_filename_out = args.out
    if bdf_filename_out is None:
        bdf_filename_out = 'nsm.out.bdf'
    log.info(f'bdf_filename_out = {bdf_filename_out}')

    # infile = args.infile
    # print(f'infile = {infile!r}')

    # from pyNastran.utils import print_bad_path
    comments = [f'nsm_id={nsm_id:d} nsm_value={value_str!r} -> {value}\n'
                for nsm_id, value_str, value in zip(nsm_ids, nsm_strings, nsm_values)]

    # bdf_filename_out = get_bdf_outfilename(
    #     bdf_filename, bdf_filename_out=None,
    #     tag='out')
    # debug = args.debug

    # from pyNastran.bdf.mesh_utils.rbe_tools import split_nsm
    from .utils_bdf import read_lax_bdf
    model = read_lax_bdf(
        bdf_filename, punch=args.punch, xref=True,
        is_strict_card_parser=not args.lax,
        log=log)
    # eids_to_fix = load_ints_from_defaults(
    #     model.rigid_elements, infile)
    # log.info(f'eids_to_fix = {eids_to_fix}')

    model2 = split_nsm(model, nsm_ids, nsm_values, comments=comments)
    model2.write_bdf(bdf_filename_out, write_header=False)


def split_nsm(model: BDF,
              nsm_ids: list[int],
              nsm_values: list[float],
              comments: list[str]=None) -> BDF:
    """
    Split a multi-typed card and distribute mass based on area/length

    Parameters
    ----------
    model : BDF()
        the model geometry
    nsm_ids : list[int]
        the ids to update
    nsm_values : list[float]
        the nonstructural mass to apply
    comments : list[str] or None; default=None -> [''] *len(nsm_ids)
        comments for each card

    Instead of:
    NSML1, 1, PCOMP,  1.2, 10, 11
    NSML1, 1, PSHELL, 1.2, 12

    We want:
    NSML1, 1, PROP,   1.2, 10, 11, 12

    split_nsm(model, [1], [1.2], comment=['mass=1.2 slinch, weight=463 lbm'])
    We need to calculate the area/length by property and distribute it appropriately
    We'll ignore the mass because the goal was always to apply 1.2 in this case.
    Femap wrote it differently than anticipated.
    """
    if comments is None:
        comments = [''] * len(nsm_ids)

    log = model.log
    model2 = BDF(log=log)
    assert len(nsm_ids) == len(nsm_values), (nsm_ids, nsm_values)
    any_area, upids_area, any_length, upids_length = _check_nsm(model, nsm_ids)

    # grab the element data
    area, length, all_area_pid, all_length_pid = get_area_length(
        model, log,
        any_area, upids_area,
        any_length, upids_length,
    )

    # apply each nsm
    for nsm_id, nsm_value, comment  in zip(nsm_ids, nsm_values, comments):
        log.info(comment)
        log.info(f'nsm_id={nsm_id}, nsm_value={nsm_value}')
        nsms = model.nsms[nsm_id]

        area_length_flags = set()
        pids = []
        for nsm in nsms:
            # print(nsm)
            nsm_type = nsm.nsm_type
            if nsm_type in {'PCOMP', 'PCOMPG', 'PSHELL'}:
                area_length_flags.add('area')
            elif nsm_type in {'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBEND'}:
                area_length_flags.add('length')
            else:
                area_length_flags.add('')
            pids.extend(nsm.ids)
            comment += nsm.comment
        assert len(pids) > 0, (nsm_id, nsm_value, pids)
        upids = np.unique(pids).tolist()

        log.debug(comment)
        assert len(area_length_flags) == 1, area_length_flags
        is_area = 'area' in area_length_flags
        is_length = 'length' in area_length_flags

        if is_area:
            uall_pids = np.unique(all_area_pid)
            log.debug(f'all_area_pid = {all_area_pid}')
        else:
            uall_pids = np.unique(all_length_pid)
        log.debug(f'all_length_pid = {all_length_pid}')
        log.debug(f'upids = {upids}')
        log.debug(f'uall_pids = {uall_pids}')

        value_by_type = defaultdict(float)
        pids_by_type = defaultdict(list)
        if is_area:
            for pid in upids:
                prop = model.properties[pid]
                if pid not in uall_pids:
                    log.warning(f'skipping; missing elements with pid={pid}')
                    print(pid, prop)
                    continue
                ipid = np.where(pid == all_area_pid)[0]
                areai = area[ipid].sum()
                value_by_type[prop.type] += areai
                pids_by_type[prop.type].append(pid)
                comment += f'pid={pid} area={areai}\n'
            assert len(pids_by_type), 'area'
            assert len(value_by_type), 'area'

        elif is_length:
            for pid in upids:
                prop = model.properties[pid]
                if pid not in uall_pids:
                    log.warning(f'skipping; missing elements with pid={pid}')
                    print(pid, prop)
                    continue
                ipid = np.where(pid == all_length_pid)[0]
                lengthi = length[ipid].sum()
                comment += f'pid={pid} length={lengthi}\n'
                value_by_type[prop.type] += lengthi
                pids_by_type[prop.type].append(pid)
            assert len(pids_by_type), 'length'
            assert len(value_by_type), 'length'
        else:
            raise RuntimeError((is_area, is_length))

        total = sum(value_by_type.values())

        # nsm_type = 'ELEMENT'
        log.info(f'total = {total}')
        log.info(f'value_by_type = {str(dict(value_by_type))}')
        assert len(value_by_type)
        for prop_type, valuei in value_by_type.items():
            ids = pids_by_type[prop_type]
            value = (valuei / total) * nsm_value

            out = model2.add_nsml1(nsm_id, prop_type, value, ids,
                             comment=comment)
            print(out)
            comment = ''

    # model2.add_nsml1(nsm_id, nsm_type, nsm_value, ids,
    #               comment=comment)
    return model2


def _check_nsm(model: BDF, nsm_ids: list[int]) -> tuple[bool, np.ndarray,
                                                        bool, np.ndarray,]:
    """
    verify all nsms will be valid before doing a bunch of work
    also identify all area/length elements

    Parameters
    ----------
    model : BDF()
        the nastran geometry
    nsm_ids : list[int]
        the nonstructural mass ids

    """
    any_area = False
    any_length = False
    pids_area = set()
    pids_length = set()
    assert len(nsm_ids) > 0, nsm_ids
    for nsm_id in nsm_ids:
        nsms = model.nsms[nsm_id]
        area_length_flags = set()
        for nsm in nsms:
            nsm_type = nsm.nsm_type
            if nsm_type in {'PCOMP', 'PCOMPG', 'PSHELL'}:
                any_area = True
                pids_area.update(nsm.ids)
                area_length_flags.add('area')
            elif nsm_type in {'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBEND'}:
                any_length = True
                pids_length.update(nsm.ids)
                area_length_flags.add('length')
            else:
                area_length_flags.add('')
            assert any_area ^ any_length, f'exclusive or failed; any_area={any_area}, any_length={any_length}'
        assert len(area_length_flags) == 1, area_length_flags

    upids_area = np.unique(list(pids_area))
    upids_length = np.unique(list(pids_length))
    return any_area, upids_area, any_length, upids_length

def get_area_length(model: BDF, log: SimpleLogger,
                    any_area: bool, upids_area: np.ndarray,
                    any_length: bool, upids_length: np.ndarray,
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """get the area/lengths for each element by property id"""
    area_ids = []
    length_ids = []
    areas = []
    lengths = []
    all_area_pids = []
    all_length_pids = []
    for eid, elem in model.elements.items():
        if any_area and elem.type in SHELLS:
            pid = elem.pid
            if pid in upids_area:
                areai = elem.Area()
                all_area_pids.append(pid)
                area_ids.append(eid)
                areas.append(areai)
        elif any_length and elem.type in LINES:
            pid = elem.pid
            if pid in upids_length:
                all_length_pids.append(pid)
                length_ids.append(eid)
                lengthi = elem.Length()
                lengths.append(lengthi)
        elif elem.type in LINES or elem.type in SHELLS:
            continue
        elif elem.type in ETYPES_NOT_CONSIDERED:
            continue
        else: #elif hasattr(elem, 'pid_ref'):
            #log.warning(elem.pid_ref)
            log.warning(elem)
    area = np.array(areas, dtype='float64')
    length = np.array(lengths, dtype='float64')
    all_area_pid = np.array(all_area_pids, dtype='int32')
    all_length_pid = np.array(all_length_pids, dtype='int32')
    return area, length, all_area_pid, all_length_pid
