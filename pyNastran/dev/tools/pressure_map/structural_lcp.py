"""
Maps structural loads (PLOAD2, PLOAD4, FORCE, MOMENT) onto load collection
points (LCPs).

This is the structural counterpart of pressure_map(method='LCP'), which
maps aero panels.  The lumping is the same (lump_loads_to_lcp in
pressure_map.py); only the load source is different.

Math
----
1. The load case (load_id, with LOAD combinations and scale factors) is
   reduced to nodal loads in the basic frame:

   - FORCE / MOMENT: the scaled vector at its grid (rectangular CID only).
   - PLOAD2 / PLOAD4 (shells): the consistent nodal loads of the pressure
     field, F_a = integral(N_a * p * n dA), so a PLOAD2 on a CQUAD4 goes
     to its 4 nodes.

       tri  : F_a = (A*n)/12 * (2*p_a + p_b + p_c)          (exact)
       quad : 2x2 Gauss on the bilinear quad (exact for a flat quad)

     PLOAD4 corner pressures P1-P4 vary linearly/bilinearly; with
     CID/N1-N3 the load acts along N (magnitude p*dA) instead of the
     element normal.  Positive pressure acts along the element normal
     (G1-G2-G3 right-hand rule), same as PLOAD2.

   The nodal loads reproduce the pressure field's total F and M exactly
   (flat elements).

2. Per group, each nodal load at x_i goes to the nearest LCP k of that
   group::

       F_k = sum F_i                         -> FORCE  (force_sid)
       M_k = sum (x_i - x_k) x F_i + M_i     -> MOMENT (moment_sid)

   so FORCE + MOMENT reproduce the group's F and M about any point.

Groups (load_control_nids)
--------------------------
Same shape as pressure_map's load_control_nids, but structural regions
are picked by property id and/or element id::

    load_control_nids = {
        'fuse': {'property_ids': [10, 11], 'nodes': fuse_lcp_nids},
        'wing': {'property_ids': [20], 'element_ids': [5001, 5002],
                 'nodes': wing_lcp_nids},
    }

 - A group's elements are the union of its property_ids and element_ids.
 - PLOAD2/PLOAD4 on a group's elements go to that group's LCPs.
 - FORCE/MOMENT go to the group whose elements use that grid.  A grid
   used by more than one group (e.g. a wing/fuselage junction) goes to
   the first group in dict order (warning).  A grid in no group is not
   mapped (warning).
 - {} or None -> one 'total' group: every load to every LCP.
 - Missing ids are dropped with a warning; an empty group raises.
   Loaded elements in no group are not mapped; elements in more than one
   group are mapped once per group (warnings).

Checks (in the log)
-------------------
 - per group: "Force/Moment Balance: Structure vs LCP (group=NAME)" with
   PASS/FAIL.  This must always pass; a FAIL is a bug.
 - total: the whole load case vs the sum of all groups.  FAIL here means
   loads were dropped or double counted by the grouping (see the
   warnings above it).

Not supported (warned and skipped): PLOAD4 on solids, PLOAD4 LINE loads,
and other load cards (PLOAD, PLOAD1, FORCE1/2, MOMENT1/2, GRAV, ...).
CTRIA6/CQUAD8 are loaded through their corner nodes only (total F/M is
still exact).

Example::

    model_out = structural_lcp_map(
        'fem.bdf', 'lcp_nodes.bdf', load_id=1,
        load_control_nids=load_control_nids,
        reference_point=np.array([807., 0., 0.]),
        lcp_load_filename='lcp_loads.bdf',
        xyz_units='in', force_units='lbf')

writes lcp_loads.bdf (FORCE/MOMENT), main_lcp_loads.bdf (SOL 101 deck
that INCLUDEs lcp_nodes.bdf and lcp_loads.bdf), and
lcp_loads.force_moment.csv.
"""
from __future__ import annotations
from typing import Optional

import numpy as np
from cpylog import SimpleLogger

from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.dev.tools.pressure_map.setup_structure import get_structure_xyz
from pyNastran.dev.tools.pressure_map.pressure_map import (
    lump_loads_to_lcp, add_lcp_force_moment_cards, _write_lcp_pressure_file)

TRI_TYPES = {'CTRIA3', 'CTRIAR', 'CTRIA6'}
QUAD_TYPES = {'CQUAD4', 'CQUADR', 'CSHEAR', 'CQUAD8'}
MIDSIDE_TYPES = {'CTRIA6', 'CQUAD8'}
GROUP_KEYS = {'property_ids', 'element_ids', 'nodes'}

# 2x2 Gauss points for the bilinear quad
_G = 1. / np.sqrt(3.)
_QUAD_XI = np.array([-1., 1., 1., -1.])
_QUAD_ETA = np.array([-1., -1., 1., 1.])
_GAUSS_2X2 = [(-_G, -_G), (_G, -_G), (_G, _G), (-_G, _G)]


def structural_lcp_map(load_filename: PathLike | BDF,
                       lcp_filename: PathLike | BDF,
                       load_id: int,
                       load_control_nids: Optional[dict[str, dict[str, list[int]]]]=None,
                       lcp_load_filename: PathLike='',
                       reference_point: Optional[np.ndarray]=None,
                       force_sid: int=2,
                       moment_sid: int=3,
                       xyz_units: str='',
                       force_units: str='',
                       pressure_str: str='',
                       log: Optional[SimpleLogger]=None) -> BDF:
    """
    Reads the models, maps load_id onto the LCPs, and writes the files.

    Parameters
    ----------
    load_filename : PathLike | BDF
        the model with the loads and the loaded elements
    lcp_filename : PathLike | BDF
        the model with the LCP GRIDs (may be the same model)
    load_id : int
        the LOAD/FORCE/PLOAD2/... set id (LOAD combinations are expanded)
    load_control_nids : dict; default=None
        load_control_nids = {
            'fuse': {'property_ids': [10, 11], 'nodes': fuse_lcp_nids},
            'wing': {'property_ids': [20], 'element_ids': [5001, 5002],
                     'nodes': wing_lcp_nids},
        }
        element ids are the union of the elements in property_ids and element_ids
    lcp_load_filename : PathLike; default=''
        '' -> no files are written
    reference_point : (3,) float ndarray; default=None -> origin
        moment reference for the logged balance only
    force_sid / moment_sid : int
        sids of the output FORCE / MOMENT cards
    xyz_units / force_units : str
        labels for the log only

    Returns
    -------
    bdf_model_out : BDF
        FORCE/MOMENT cards (basic frame) on the LCPs
    """
    load_model = _read(load_filename, log)
    lcp_model = load_model if lcp_filename is load_filename else _read(lcp_filename, load_model.log)
    bdf_model_out = structural_loads_to_lcp_model(
        load_model, lcp_model, load_id, load_control_nids=load_control_nids,
        reference_point=reference_point, force_sid=force_sid, moment_sid=moment_sid,
        xyz_units=xyz_units, force_units=force_units)

    if lcp_load_filename:
        _write_lcp_pressure_file(
            bdf_model_out, str(lcp_load_filename), pressure_str,
            lcp_filename if isinstance(lcp_filename, BDF) else str(lcp_filename),
            map_type='force_moment', method='LCP',
            force_sid=force_sid, moment_sid=moment_sid)
    return bdf_model_out


def _read(bdf_filename: PathLike | BDF, log: Optional[SimpleLogger]) -> BDF:
    if isinstance(bdf_filename, BDF):
        return bdf_filename
    return read_bdf(bdf_filename, xref=True, log=log)


def structural_loads_to_lcp_model(load_model: BDF,
                                  lcp_model: BDF,
                                  load_id: int,
                                  load_control_nids: Optional[dict[str, dict[str, list[int]]]]=None,
                                  reference_point: Optional[np.ndarray]=None,
                                  force_sid: int=2,
                                  moment_sid: int=3,
                                  xyz_units: str='',
                                  force_units: str='') -> BDF:
    """
    Maps a structural load case onto LCPs (see the module docstring).

    Parameters
    ----------
    load_model : BDF
        cross-referenced model with the loads and elements
    lcp_model : BDF
        model with the LCP GRIDs (only GRIDs are used); may be load_model

    Returns
    -------
    bdf_model_out : BDF
        one FORCE (force_sid) + MOMENT (moment_sid) per loaded LCP per group
    """
    log = load_model.log
    if reference_point is None:
        reference_point = np.zeros(3)
    reference_point = np.asarray(reference_point, dtype='float64')
    if load_control_nids is None:
        load_control_nids = {}

    load_nids, load_xyz_all = _sorted_xyz(load_model, xyz_units)
    nodal = get_nodal_loads(load_model, load_id, load_nids, load_xyz_all)

    lcp_nids_all, lcp_xyz_all = _sorted_xyz(lcp_model, xyz_units)
    groups = _simplify_structural_groups(load_model, load_control_nids, lcp_nids_all, log)

    # pressure rows (one per element node) -> groups by element id
    p_eid = nodal['pressure_eid']
    p_xyz = nodal['pressure_xyz']
    p_force = nodal['pressure_force']
    ingroup_count = np.zeros(len(p_eid), dtype='int32')

    # point loads -> groups by the grids of each group's elements
    pt_nid = nodal['point_nid']
    pt_group = _assign_point_loads(load_model, groups, pt_nid, log)

    bdf_model_out = BDF(log=log)
    total_out_f = np.zeros(3)
    total_out_m = np.zeros(3)
    for igroup, (group_name, group) in enumerate(groups.items()):
        ipressure = np.isin(p_eid, group['element_ids'])
        ingroup_count += ipressure
        ipoint = (pt_group == igroup)

        xyz = np.vstack([p_xyz[ipressure], nodal['point_xyz'][ipoint]])
        force = np.vstack([p_force[ipressure], nodal['point_force'][ipoint]])
        moment = np.vstack([np.zeros((ipressure.sum(), 3)), nodal['point_moment'][ipoint]])

        ilcp = np.searchsorted(lcp_nids_all, group['nodes'])
        nids, lcp_xyz, lcp_force, lcp_moment = lump_loads_to_lcp(
            xyz, force, moment, group['nodes'], lcp_xyz_all[ilcp, :])

        comment = f'load_id={load_id}; group={group_name} (structural_loads_to_lcp_model)'
        add_lcp_force_moment_cards(
            bdf_model_out, nids, lcp_force, lcp_moment, force_sid, moment_sid,
            fcomment='Force; ' + comment, mcomment='Moment; ' + comment)

        in_f, in_m = _sum_about(xyz, force, moment, reference_point)
        out_f, out_m = _sum_about(lcp_xyz, lcp_force, lcp_moment, reference_point)
        total_out_f += out_f
        total_out_m += out_m
        _log_structural_balance(
            log, f'Force/Moment Balance: Structure vs LCP (group={group_name})',
            in_f, in_m, out_f, out_m, reference_point, xyz_units, force_units)

    unassigned = np.unique(p_eid[ingroup_count == 0])
    if len(unassigned):
        log.warning(f'LCP: loaded element ids={unassigned} are not in any group; their load is not mapped')
    multiple = np.unique(p_eid[ingroup_count > 1])
    if len(multiple):
        log.warning(f'LCP: loaded element ids={multiple} are in more than one group; their load is mapped more than once')

    case_f, case_m = _sum_about(
        np.vstack([p_xyz, nodal['point_xyz']]),
        np.vstack([p_force, nodal['point_force']]),
        np.vstack([np.zeros(p_xyz.shape), nodal['point_moment']]),
        reference_point)
    _log_structural_balance(
        log, f'Force/Moment Balance: load_id={load_id} vs all LCP groups',
        case_f, case_m, total_out_f, total_out_m, reference_point, xyz_units, force_units,
        fail_msg='loads were dropped or double counted by the grouping (see the warnings above)')
    return bdf_model_out


def _sorted_xyz(model: BDF, xyz_units: str) -> tuple[np.ndarray, np.ndarray]:
    """node ids (sorted) and basic xyz"""
    nids, xyz = get_structure_xyz(model, xyz_units=xyz_units)
    isort = np.argsort(nids)
    return nids[isort], xyz[isort, :]


def get_nodal_loads(model: BDF,
                    load_id: int,
                    nids: np.ndarray,
                    xyz_cid0: np.ndarray) -> dict[str, np.ndarray]:
    """
    Reduces a load case to nodal loads in the basic frame.

    Parameters
    ----------
    model : BDF
        cross-referenced
    load_id : int
        LOAD combinations and scale factors are applied
    nids, xyz_cid0 : (n,) int / (n, 3) float ndarray
        sorted node ids and basic xyz (get_structure_xyz)

    Returns
    -------
    nodal : dict
        pressure_eid/nid  : (r,) int    one row per loaded element node
        pressure_xyz      : (r, 3)
        pressure_force    : (r, 3)
        point_nid         : (s,) int    one row per FORCE/MOMENT
        point_xyz/force/moment : (s, 3)
    """
    log = model.log
    loads, scale_factors, unused_is_grav = model.get_reduced_loads(load_id)

    tri_eid, tri_nids, tri_p, tri_dir = [], [], [], []
    quad_eid, quad_nids, quad_p, quad_dir = [], [], [], []
    pt_nid, pt_force, pt_moment = [], [], []
    skipped: dict[str, int] = {}
    midside_types = set()
    no_dir = np.full(3, np.nan)

    def _add_shell(eid: int, pressures: np.ndarray, direction: np.ndarray) -> None:
        elem = model.elements[eid]
        if elem.type in TRI_TYPES:
            tri_eid.append(eid)
            tri_nids.append(elem.node_ids[:3])
            tri_p.append(pressures[:3])
            tri_dir.append(direction)
        elif elem.type in QUAD_TYPES:
            quad_eid.append(eid)
            quad_nids.append(elem.node_ids[:4])
            quad_p.append(pressures[:4])
            quad_dir.append(direction)
        else:
            key = f'{elem.type} (pressure)'
            skipped[key] = skipped.get(key, 0) + 1
            return
        if elem.type in MIDSIDE_TYPES:
            midside_types.add(elem.type)

    for load, scale in zip(loads, scale_factors):
        if load.type in {'FORCE', 'MOMENT'}:
            _check_rectangular(model, load.Cid(), load)
            vector = scale * np.asarray(load.to_global(), dtype='float64')
            pt_nid.append(load.node_id)
            if load.type == 'FORCE':
                pt_force.append(vector)
                pt_moment.append(np.zeros(3))
            else:
                pt_force.append(np.zeros(3))
                pt_moment.append(vector)
        elif load.type == 'PLOAD2':
            pressures = np.full(4, scale * load.pressure)
            for eid in load.element_ids:
                _add_shell(eid, pressures, no_dir)
        elif load.type == 'PLOAD4':
            if load.surf_or_line != 'SURF':
                skipped['PLOAD4 (LINE)'] = skipped.get('PLOAD4 (LINE)', 0) + 1
                continue
            if load.g1 is not None or model.elements[load.element_ids[0]].type not in (TRI_TYPES | QUAD_TYPES):
                skipped['PLOAD4 (solid)'] = skipped.get('PLOAD4 (solid)', 0) + 1
                continue
            direction = no_dir
            nvector = np.asarray(load.nvector, dtype='float64')
            if np.linalg.norm(nvector) > 0.:
                cid = load.Cid()
                _check_rectangular(model, cid, load)
                direction = np.asarray(model.Coord(cid).transform_vector_to_global(nvector),
                                       dtype='float64').ravel()
                direction /= np.linalg.norm(direction)
            pressures = scale * np.asarray(load.pressures, dtype='float64')
            for eid in load.element_ids:
                _add_shell(eid, pressures, direction)
        else:
            skipped[load.type] = skipped.get(load.type, 0) + 1

    if skipped:
        log.warning(f'LCP: unsupported loads were skipped (card: count); {skipped}')
    if midside_types:
        log.warning(f'LCP: {sorted(midside_types)} pressure goes to the corner nodes only '
                    '(total F/M is still exact)')

    def _xyz(node_ids: np.ndarray) -> np.ndarray:
        inode = np.searchsorted(nids, node_ids)
        assert np.array_equal(nids[inode], node_ids), 'missing nodes'
        return xyz_cid0[inode]

    eid_list, nid_list, force_list = [], [], []
    for eids, node_ids, pressures, directions, nnode in [
            (tri_eid, tri_nids, tri_p, tri_dir, 3),
            (quad_eid, quad_nids, quad_p, quad_dir, 4)]:
        if not eids:
            continue
        node_ids = np.array(node_ids, dtype='int64')
        xyz = _xyz(node_ids.ravel()).reshape(-1, nnode, 3)
        nodal_force = consistent_pressure_nodal_forces(
            xyz, np.array(pressures), np.array(directions))
        eid_list.append(np.repeat(np.array(eids, dtype='int64'), nnode))
        nid_list.append(node_ids.ravel())
        force_list.append(nodal_force.reshape(-1, 3))

    if eid_list:
        p_eid = np.hstack(eid_list)
        p_nid = np.hstack(nid_list)
        p_force = np.vstack(force_list)
    else:
        p_eid = np.zeros(0, dtype='int64')
        p_nid = np.zeros(0, dtype='int64')
        p_force = np.zeros((0, 3))

    pt_nid_array = np.array(pt_nid, dtype='int64')
    nodal = {
        'pressure_eid': p_eid,
        'pressure_nid': p_nid,
        'pressure_xyz': _xyz(p_nid) if len(p_nid) else np.zeros((0, 3)),
        'pressure_force': p_force,
        'point_nid': pt_nid_array,
        'point_xyz': _xyz(pt_nid_array) if len(pt_nid_array) else np.zeros((0, 3)),
        'point_force': np.array(pt_force).reshape(-1, 3),
        'point_moment': np.array(pt_moment).reshape(-1, 3),
    }
    log.info(f'LCP: load_id={load_id}: {len(set(p_eid.tolist()))} loaded shells, '
             f'{len(pt_nid_array)} FORCE/MOMENT cards')
    return nodal


def consistent_pressure_nodal_forces(xyz: np.ndarray,
                                     pressures: np.ndarray,
                                     directions: np.ndarray) -> np.ndarray:
    """
    Consistent nodal forces F_a = integral(N_a * p * v dA) for tris/quads.

    Parameters
    ----------
    xyz : (ne, nnode, 3) float ndarray
        corner coordinates (nnode=3: tri, nnode=4: quad), element order
    pressures : (ne, nnode) float ndarray
        corner pressures (linear/bilinear field)
    directions : (ne, 3) float ndarray
        unit load direction (PLOAD4 N); a nan row -> element normal

    Returns
    -------
    force : (ne, nnode, 3) float ndarray
    """
    xyz = np.asarray(xyz, dtype='float64')
    pressures = np.asarray(pressures, dtype='float64')
    directions = np.asarray(directions, dtype='float64').reshape(-1, 3)
    ne, nnode = xyz.shape[:2]
    use_dir = ~np.isnan(directions[:, 0])

    def _area_vector(dA_vec: np.ndarray) -> np.ndarray:
        """n*dA, or |dA|*N when the load direction is given"""
        out = dA_vec.copy()
        if use_dir.any():
            mag = np.linalg.norm(dA_vec[use_dir], axis=1)
            out[use_dir] = mag[:, np.newaxis] * directions[use_dir]
        return out

    if nnode == 3:
        area_vec = 0.5 * np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
        vec = _area_vector(area_vec)
        weight = (pressures + pressures.sum(axis=1, keepdims=True)) / 12.  # (2p_a + p_b + p_c)/12
        return weight[:, :, np.newaxis] * vec[:, np.newaxis, :]

    assert nnode == 4, nnode
    force = np.zeros((ne, 4, 3))
    for xi, eta in _GAUSS_2X2:
        shape = 0.25 * (1. + xi * _QUAD_XI) * (1. + eta * _QUAD_ETA)       # (4,)
        dshape_dxi = 0.25 * _QUAD_XI * (1. + eta * _QUAD_ETA)
        dshape_deta = 0.25 * _QUAD_ETA * (1. + xi * _QUAD_XI)
        x_xi = np.einsum('a,eaj->ej', dshape_dxi, xyz)
        x_eta = np.einsum('a,eaj->ej', dshape_deta, xyz)
        vec = _area_vector(np.cross(x_xi, x_eta))                         # n dA/(dxi deta)
        p_gauss = pressures @ shape                                       # (ne,)
        force += shape[np.newaxis, :, np.newaxis] * (p_gauss[:, np.newaxis] * vec)[:, np.newaxis, :]
    return force


def _check_rectangular(model: BDF, cid: int, load) -> None:
    if cid == 0:
        return
    coord = model.Coord(cid)
    if coord.Type != 'R':
        raise NotImplementedError(f'only rectangular CIDs are supported; cid={cid} '
                                  f'({coord.type}) on\n{load}')


def _simplify_structural_groups(model: BDF,
                                load_control_nids: dict[str, dict[str, list[int]]],
                                lcp_nids: np.ndarray,
                                log: SimpleLogger) -> dict[str, dict[str, np.ndarray]]:
    """
    Returns
    -------
    {group_name: {'element_ids': (n,) int ndarray, 'nodes': (m,) int ndarray}}
        sorted, unique ids
    """
    all_eids = np.array(sorted(model.elements), dtype='int64')
    if len(load_control_nids) == 0:
        return {'total': {'element_ids': all_eids, 'nodes': lcp_nids}}

    all_pids = np.array([_pid(model.elements[eid]) for eid in all_eids], dtype='int64')
    groups = {}
    for group_name, group in load_control_nids.items():
        bad_keys = set(group) - GROUP_KEYS
        if bad_keys:
            raise ValueError(f'LCP group={group_name!r} has invalid keys={sorted(bad_keys)}; '
                             f'allowed={sorted(GROUP_KEYS)}')
        pids = np.unique(np.asarray(group.get('property_ids', []), dtype='int64'))
        eids = np.unique(np.asarray(group.get('element_ids', []), dtype='int64'))
        lcps = np.unique(np.asarray(group.get('nodes', []), dtype='int64'))
        if len(pids) == 0 and len(eids) == 0:
            raise ValueError(f'LCP group={group_name!r} needs property_ids and/or element_ids')

        missing = np.setdiff1d(pids, all_pids)
        if len(missing):
            log.warning(f'LCP group={group_name!r} has {len(missing)} missing property_ids; {missing}')
        missing = np.setdiff1d(eids, all_eids)
        if len(missing):
            log.warning(f'LCP group={group_name!r} has {len(missing)} missing element_ids; {missing}')
        missing = np.setdiff1d(lcps, lcp_nids)
        if len(missing):
            log.warning(f'LCP group={group_name!r} has {len(missing)} missing nodes; {missing}')

        group_eids = np.union1d(all_eids[np.isin(all_pids, pids)],
                                np.intersect1d(eids, all_eids))
        lcps = np.intersect1d(lcps, lcp_nids)
        if len(group_eids) == 0:
            raise ValueError(f'LCP group={group_name!r} has no elements in the load model')
        if len(lcps) == 0:
            raise ValueError(f'LCP group={group_name!r} has no nodes in the LCP model')
        groups[group_name] = {'element_ids': group_eids, 'nodes': lcps}
    return groups


def _pid(elem) -> int:
    """property id, or -1 for elements without one (e.g., CONROD, CONM2-like)"""
    try:
        pid = elem.Pid()
    except AttributeError:
        return -1
    return -1 if pid is None else pid


def _assign_point_loads(model: BDF,
                        groups: dict[str, dict[str, np.ndarray]],
                        point_nids: np.ndarray,
                        log: SimpleLogger) -> np.ndarray:
    """
    Returns the group index for each FORCE/MOMENT row (-1: not mapped).
    A grid goes to the first group whose elements use it.
    """
    igroup = np.full(len(point_nids), -1, dtype='int32')
    if len(point_nids) == 0:
        return igroup
    if list(groups) == ['total'] and len(groups['total']['element_ids']) == len(model.elements):
        igroup[:] = 0
        return igroup

    counts = np.zeros(len(point_nids), dtype='int32')
    for i, group in enumerate(groups.values()):
        group_nids = set()
        for eid in group['element_ids']:
            group_nids.update(nid for nid in model.elements[eid].node_ids if nid is not None)
        in_group = np.isin(point_nids, np.array(sorted(group_nids), dtype='int64'))
        counts += in_group
        igroup[in_group & (igroup == -1)] = i

    names = list(groups)
    shared = np.unique(point_nids[counts > 1])
    if len(shared):
        log.warning(f'LCP: FORCE/MOMENT grids={shared} are used by more than one group; '
                    f'each goes to the first one (dict order: {names})')
    unassigned = np.unique(point_nids[igroup == -1])
    if len(unassigned):
        log.warning(f'LCP: FORCE/MOMENT grids={unassigned} are not on any group element; '
                    'their load is not mapped')
    return igroup


def _sum_about(xyz: np.ndarray, force: np.ndarray, moment: np.ndarray,
               reference_point: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    xyz = np.asarray(xyz).reshape(-1, 3)
    force = np.asarray(force).reshape(-1, 3)
    moment = np.asarray(moment).reshape(-1, 3)
    total_force = force.sum(axis=0)
    total_moment = (np.cross(xyz - reference_point, force) + moment).sum(axis=0)
    return total_force, total_moment


def _log_structural_balance(log: SimpleLogger,
                            title: str,
                            in_force: np.ndarray, in_moment: np.ndarray,
                            out_force: np.ndarray, out_moment: np.ndarray,
                            reference_point: np.ndarray,
                            xyz_units: str='', force_units: str='',
                            fail_msg: str='this is a bug',
                            rtol: float=1e-8) -> bool:
    """Logs input vs LCP F/M about reference_point with a PASS/FAIL line."""
    def _f(label, vals):
        return f'  {label:22s} {vals[0]:14.6g} {vals[1]:14.6g} {vals[2]:14.6g}'

    fscale = max(np.abs(in_force).max(), np.abs(out_force).max(), 1e-30)
    mscale = max(np.abs(in_moment).max(), np.abs(out_moment).max(), 1e-30)
    is_pass = (np.abs(out_force - in_force).max() <= rtol * fscale and
               np.abs(out_moment - in_moment).max() <= rtol * mscale)

    log.info('================================================================')
    log.info(f'  {title}')
    log.info(f'  reference_point ({xyz_units}) = {reference_point}; force units={force_units!r}')
    log.info(f'  {"":22s} {"Fx":>14s} {"Fy":>14s} {"Fz":>14s}')
    log.info(_f('Structure Force', in_force))
    log.info(_f('LCP       Force', out_force))
    log.info(_f('Delta     Force', out_force - in_force))
    log.info(f'  {"":22s} {"Mx":>14s} {"My":>14s} {"Mz":>14s}')
    log.info(_f('Structure Moment', in_moment))
    log.info(_f('LCP       Moment', out_moment))
    log.info(_f('Delta     Moment', out_moment - in_moment))
    if is_pass:
        log.info(f'  PASS: LCP F/M match the structure loads (rtol={rtol})')
    else:
        log.warning(f'  FAIL: LCP F/M do not match the structure loads ({title}); {fail_msg}')
    log.info('================================================================')
    return is_pass
