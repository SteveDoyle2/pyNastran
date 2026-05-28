from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def get_bad_shells(model: BDF,
                   min_theta: float = 0.1,
                   max_theta: float = 175.,
                   max_skew: float = 70.,
                   max_aspect_ratio: float = 100.,
                   max_taper_ratio: float = 4.0,
                   max_warp: float = 90.) -> np.ndarray:
    """
    Gets element IDs of shells that exceed quality thresholds.

    Parameters
    ----------
    model : BDF
        the model (must have been setup)
    min_theta : float; default=0.1
        minimum interior angle (degrees)
    max_theta : float; default=175.
        maximum interior angle (degrees)
    max_skew : float; default=70.
        maximum skew angle (degrees)
    max_aspect_ratio : float; default=100.
        maximum aspect ratio
    max_taper_ratio : float; default=4.0
        maximum taper ratio (quads only)
    max_warp : float; default=90.
        maximum warp angle (degrees)

    Returns
    -------
    bad_eids : (n,) int ndarray
        element IDs that fail one or more quality criteria

    """
    SHELL_TYPES = {
        'CTRIA3', 'CTRIA6', 'CTRIAR',
        'CQUAD4', 'CQUAD8', 'CQUADR',
    }
    (element_ids, taper_ratio, area_ratio, max_skew_arr, aspect_ratio,
     min_theta_arr, max_theta_arr, dideal_theta, min_edge_length, max_warp_arr) = model.quality()

    shell_mask = np.zeros(len(element_ids), dtype=bool)
    for card in model.element_cards:
        if card.n > 0 and card.type in SHELL_TYPES:
            shell_mask |= np.isin(element_ids, card.element_id)

    bad = np.zeros(len(element_ids), dtype=bool)
    with np.errstate(invalid='ignore'):
        bad |= (min_theta_arr < min_theta)
        bad |= (max_theta_arr > max_theta)
        bad |= (max_skew_arr > max_skew)
        bad |= (aspect_ratio > max_aspect_ratio)
        bad |= (taper_ratio > max_taper_ratio)
        bad |= (max_warp_arr > max_warp)

    bad &= shell_mask
    bad_eids = element_ids[bad]
    return bad_eids


def delete_bad_shells(model: BDF,
                      min_theta: float = 0.1,
                      max_theta: float = 175.,
                      max_skew: float = 70.,
                      max_aspect_ratio: float = 100.,
                      max_taper_ratio: float = 4.0,
                      max_warp: float = 90.) -> np.ndarray:
    """
    Removes shell elements that exceed quality thresholds.

    Parameters
    ----------
    model : BDF
        the model (must have been setup)
    min_theta / max_theta / max_skew / max_aspect_ratio /
    max_taper_ratio / max_warp : float
        quality thresholds (degrees for angles)

    Returns
    -------
    bad_eids : (n,) int ndarray
        element IDs that were removed

    """
    bad_eids = get_bad_shells(
        model,
        min_theta=min_theta, max_theta=max_theta,
        max_skew=max_skew, max_aspect_ratio=max_aspect_ratio,
        max_taper_ratio=max_taper_ratio, max_warp=max_warp)

    if len(bad_eids) == 0:
        return bad_eids

    bad_set = set(bad_eids.tolist())
    SHELL_TYPES = {
        'CTRIA3', 'CTRIA6', 'CTRIAR',
        'CQUAD4', 'CQUAD8', 'CQUADR',
    }
    for card in model.element_cards:
        if card.n > 0 and card.type in SHELL_TYPES:
            keep = ~np.isin(card.element_id, bad_eids)
            if not np.all(keep):
                i_keep = np.where(keep)[0]
                card.__apply_slice__(card, i_keep)

    model.log.info(f'deleted {len(bad_eids)} bad shells')
    return bad_eids


def convert_bad_quads_to_tris(model: BDF,
                              min_edge_length: float = 0.0) -> int:
    """
    Converts degenerate CQUAD4 elements (with collapsed edges) to CTRIA3.

    Parameters
    ----------
    model : BDF
        the model (must have been setup)
    min_edge_length : float; default=0.0
        minimum edge length; elements with shorter edges are collapsed

    Returns
    -------
    nconverted : int
        number of quads converted to tris

    """
    cquad4 = model.cquad4
    if cquad4.n == 0:
        return 0

    nodes = cquad4.nodes
    n1 = nodes[:, 0]
    n2 = nodes[:, 1]
    n3 = nodes[:, 2]
    n4 = nodes[:, 3]

    collapsed_12 = (n1 == n2)
    collapsed_23 = (n2 == n3)
    collapsed_34 = (n3 == n4)
    collapsed_41 = (n4 == n1)
    is_collapsed = collapsed_12 | collapsed_23 | collapsed_34 | collapsed_41

    if min_edge_length > 0.0:
        xyz = model.grid.xyz_cid0()
        nid = model.grid.node_id
        inode = np.searchsorted(nid, nodes)
        p1 = xyz[inode[:, 0]]
        p2 = xyz[inode[:, 1]]
        p3 = xyz[inode[:, 2]]
        p4 = xyz[inode[:, 3]]
        len12 = np.linalg.norm(p2 - p1, axis=1)
        len23 = np.linalg.norm(p3 - p2, axis=1)
        len34 = np.linalg.norm(p4 - p3, axis=1)
        len41 = np.linalg.norm(p1 - p4, axis=1)
        short_edge = ((len12 < min_edge_length) | (len23 < min_edge_length) |
                      (len34 < min_edge_length) | (len41 < min_edge_length))
        is_collapsed |= short_edge

    if not np.any(is_collapsed):
        return 0

    i_bad = np.where(is_collapsed)[0]
    nconverted = len(i_bad)

    for idx in i_bad:
        eid = cquad4.element_id[idx]
        pid = cquad4.property_id[idx]
        ns = nodes[idx]
        unique_nodes = []
        seen = set()
        for n in ns:
            if n not in seen:
                unique_nodes.append(n)
                seen.add(n)
        if len(unique_nodes) >= 3:
            model.ctria3.add(eid, pid, unique_nodes[:3])

    model.ctria3.parse_cards()

    i_keep = np.where(~is_collapsed)[0]
    cquad4.__apply_slice__(cquad4, i_keep)

    model.log.info(f'converted {nconverted} CQUAD4s to CTRIA3s')
    return nconverted
