from __future__ import annotations
from itertools import chain
from io import StringIO, IOBase
from typing import Optional, TYPE_CHECKING

import numpy as np

from pyNastran.bdf.bdf import BDF
from pyNastran.utils import PathLike
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.utils.mathematics import roundup
from pyNastran.bdf.mesh_utils.bdf_renumber import _get_bdf_model, _write_bdf

if TYPE_CHECKING:
    from cpylog import SimpleLogger


def bdf_renumber_exclude(bdf_filename: PathLike | BDF | StringIO,
                         bdf_filename_out: str,
                         size: int=8, is_double: bool=False,
                         range_id_dict: Optional[dict[str, list]]=None,
                         #round_ids: bool=False,
                         cards_to_skip: Optional[list[str]]=None,
                         punch: bool=False,
                         log: Optional[SimpleLogger]=None,
                         debug: bool=False) -> BDF:
    model = _get_bdf_model(bdf_filename, cards_to_skip=cards_to_skip,
                           punch=punch, log=log, debug=debug)

    for key, start_stop_range  in range_id_dict.items():
        assert key in {'nid', 'eid', 'pid', 'mid'}, key
        # 'set_id': max(model.sets) + 1 if model.sets else 1,
        # 'spline_id': max(model.splines) + 1 if model.splines else 1,
        # 'caero_id': max(caero.box_ids[-1, -1]
        #'cid'

    nid_map, unused_reverse_nid_map = _create_nid_maps(model, range_id_dict)
    eid_map, unused_reverse_eid_map = _create_eid_maps(model, range_id_dict)
    pid_map, unused_reverse_pid_map = _create_pid_maps(model, range_id_dict)
    #mid_map, unused_reverse_mid_map = _create_mid_maps(model, range_id_dict)

    _update_nodes(model, nid_map)
    _update_elements(model, eid_map)
    _update_properties(model, pid_map)
    #_update_materials(model, mid_map)
    _write_bdf(model, bdf_filename_out, size=size, is_double=is_double)
    return model


def _update_nodes(model: BDF, nid_map: dict[int, int]) -> None:
    """updates the nodes"""
    if len(nid_map) == 0:
        return

    for nid, node in sorted(model.nodes.items()):
        nid_new = nid_map[nid]
        # print('nid=%s -> %s' % (nid, nid_new))
        node.nid = nid_new


def _update_elements(model: BDF, eid_map: dict[int, int]) -> None:
    """updates the elements"""
    if len(eid_map) == 0:
        return

    for eid, elem in sorted(model.elements.items()):
        eid_new = eid_map[eid]
        elem.eid = eid_new
    for eid, elem in sorted(model.rigid_elements.items()):
        eid_new = eid_map[eid]
        elem.eid = eid_new
    for eid, elem in sorted(model.masses.items()):
        eid_new = eid_map[eid]
        elem.eid = eid_new


def _update_properties(model: BDF, pid_map: dict[int, int]) -> None:
    """updates the properties"""
    if len(pid_map) == 0:
        return

    for pid, prop in sorted(model.properties.items()):
        pid_new = pid_map[pid]
        prop.pid = pid_new
    for pid, prop in sorted(model.properties_mass.items()):
        pid_new = pid_map[pid]
        prop.pid = pid_new


def _create_nid_maps(model, range_id_dict: dict[str, tuple[int, int, list[np.ndarray]]],
                     ) -> tuple[dict[int, int], dict[int, int]]:
    """

    Parameters
    ----------
    model : BDF
    range_id_dict : dict[str, tuple[int, int, list[np.ndarray]]]'
        key : str
        value = (start, stop, exclude_list)

    Returns
    -------
    nid_map : dict[int, int]
        nid_new = nid_map[nid]
    reverse_nid_map : dict[int, int]
        nid = reverse_nid_map[nid_new]

    """
    key = 'nid'
    if key not in range_id_dict:
        return {}, {}

    datai = range_id_dict[key]
    start, stop, exclude_id_set = _get_start_stop_exclude_set(*datai)
    # ids_actual = set(list(model.nodes))
    # ids_output = []

    range_id_dict[key]
    # ----------------------
    spoints = list(model.spoints.keys())
    epoints = list(model.epoints.keys())
    nids = model.nodes.keys()
    nids_spoints_epoints = sorted(chain(nids, spoints, epoints))

    nid_map, reverse_nid_map = _get_map_reverse_map(
        key, start, stop, nids_spoints_epoints, exclude_id_set)
    return nid_map, reverse_nid_map


def _create_eid_maps(model, range_id_dict: dict[str, tuple[int, int, list[np.ndarray]]],
                     ) -> tuple[dict[int, int], dict[int, int]]:
    """

    Parameters
    ----------
    model : BDF
    range_id_dict : dict[str, tuple[int, int, list[np.ndarray]]]'
        key : str
        value = (start, stop, exclude_list)

    Returns
    -------
    nid_map : dict[int, int]
        nid_new = nid_map[nid]
    reverse_nid_map : dict[int, int]
        nid = reverse_nid_map[nid_new]

    """
    key = 'eid'
    if key not in range_id_dict:
        return {}, {}

    datai = range_id_dict[key]
    start, stop, exclude_id_set = _get_start_stop_exclude_set(*datai)

    range_id_dict[key]
    # ----------------------
    all_eids = []
    if len(model.elements):
        all_eids.extend(list(model.elements))
    if len(model.rigid_elements):
        all_eids.extend(list(model.rigid_elements))
    if len(model.masses):
        all_eids.extend(list(model.masses))
    all_eids.sort()

    eid_map, reverse_eid_map = _get_map_reverse_map(
        key, start, stop, all_eids, exclude_id_set)
    return eid_map, reverse_eid_map


def _create_pid_maps(model, range_id_dict: dict[str, tuple[int, int, list[np.ndarray]]],
                     ) -> tuple[dict[int, int], dict[int, int]]:
    """

    Parameters
    ----------
    model : BDF
    range_id_dict : dict[str, tuple[int, int, list[np.ndarray]]]'
        key : str
        value = (start, stop, exclude_list)

    Returns
    -------
    nid_map : dict[int, int]
        nid_new = nid_map[nid]
    reverse_nid_map : dict[int, int]
        nid = reverse_nid_map[nid_new]

    """
    key = 'pid'
    if key not in range_id_dict:
        return {}, {}

    datai = range_id_dict[key]
    start, stop, exclude_id_set = _get_start_stop_exclude_set(*datai)
    range_id_dict[key]
    # ----------------------
    all_pids = []
    if len(model.properties):
        all_pids.extend(list(model.properties))
    if len(model.properties_mass):
        all_pids.extend(list(model.properties_mass))
    all_pids.sort()

    pid_map, reverse_pid_map = _get_map_reverse_map(
        key, start, stop, all_pids, exclude_id_set)
    return pid_map, reverse_pid_map


def _get_map_reverse_map(key: str,
                         id_start: int, id_stop: int,
                         ids: list[int],
                         exclude_id_set: set[int]):
    """
    Ideally this would check to see if values are in range and
    not renumber them if they're fine
    """
    nid_map = {}
    reverse_nid_map = {}

    set_ids = set(ids)
    ids_to_update = list(set_ids - exclude_id_set)
    ids_to_update.sort()

    ids_to_lock_set = set_ids.intersection(exclude_id_set)
    ids_to_lock = list(ids_to_lock_set)
    ids_to_lock.sort()
    print(f'{key}_ids_excluded = {ids_to_lock}')
    for nid in ids_to_lock:
        nid_map[nid] = nid
        reverse_nid_map[nid] = nid

    i = id_start
    for nid in ids_to_update:
        if nid in ids_to_lock_set:
            # already handled
            continue

        # get a valid next value
        while i in exclude_id_set:
            i += 1
        nid_map[nid] = i
        reverse_nid_map[i] = nid
        i += 1
    assert i < id_stop, (key, id_stop)
    return nid_map, reverse_nid_map


def _get_start_stop_exclude_set(start: int, stop: int,
                                *exclude_list: list[np.ndarray]) -> tuple[int, int, set[int]]:
    #start, stop, excludes_list = range_id_dict[key]
    if len(exclude_list) == 1:
        exclude_ids = np.hstack(exclude_list)
    else:
        exclude_ids = exclude_list[0]
    exclude_id_set = set(exclude_ids.tolist())
    return start, stop, exclude_id_set
