"""
Torsion constant for a thin-walled section defined by a soup of wall segments.

The section that comes out of a cutting plane is an unordered set of straight
segments (one per cut shell element), each with a thickness and a shear
modulus.  The torsional stiffness of such a section is *not* ``G*(Ix+Iz)``;
the polar moment is only correct for a circular section.  For a thin wall the
correct answers are

- closed (one or more cells): Bredt-Batho, ``GJ = T`` for a unit twist rate,
  found by solving for the constant shear flow in each cell
- open (no closed cell): ``GJ = sum(G_i*s_i*t_i**3)/3``

This module recovers the cell topology from the segment connectivity and
solves the multi-cell problem.

"""
from __future__ import annotations
from collections import defaultdict
from typing import Any, Optional

import numpy as np


def _weld(p1: np.ndarray, p2: np.ndarray,
          tol: float) -> tuple[np.ndarray, np.ndarray, dict[int, np.ndarray]]:
    """Merge coincident wall ends into shared node ids.

    The cutter emits an independent pair of node ids per cut element, so the
    section topology has to be recovered geometrically.  Points are binned on
    a grid of size ``tol`` and matched against the 3x3 neighborhood so that
    a pair straddling a bin boundary still welds.
    """
    pts = np.vstack([p1, p2])
    keys = np.round(pts / tol).astype('int64')
    lookup: dict[tuple[int, int], int] = {}
    xy: dict[int, np.ndarray] = {}
    ids = np.empty(len(pts), dtype='int64')
    for i, (kx, kz) in enumerate(keys):
        found = -1
        for dx in (0, -1, 1):
            for dz in (0, -1, 1):
                j = lookup.get((kx + dx, kz + dz))
                if j is not None and np.linalg.norm(pts[i] - xy[j]) <= 2. * tol:
                    found = j
                    break
            if found >= 0:
                break
        if found < 0:
            found = len(xy)
            xy[found] = pts[i]
            lookup[(kx, kz)] = found
        ids[i] = found
    n = len(p1)
    return ids[:n], ids[n:], xy


def _build_graph(nid1: np.ndarray,
                 nid2: np.ndarray,
                 xy: dict[int, np.ndarray],
                 delta: np.ndarray) -> tuple[dict[int, set[int]],
                                             dict[tuple[int, int], float]]:
    """Undirected graph keyed on node id, with parallel walls merged.

    Two walls between the same pair of nodes act as springs in parallel for
    torsion, so their compliances ``delta = s/(G*t)`` combine reciprocally.
    """
    adjacency: dict[int, set[int]] = defaultdict(set)
    compliance: dict[tuple[int, int], float] = {}
    for n1, n2, di in zip(nid1, nid2, delta):
        if n1 == n2 or not np.isfinite(di) or di <= 0.:
            continue
        key = (n1, n2) if n1 < n2 else (n2, n1)
        if key in compliance:
            # parallel walls: 1/d = 1/d1 + 1/d2
            d0 = compliance[key]
            compliance[key] = d0 * di / (d0 + di)
        else:
            compliance[key] = di
        adjacency[n1].add(n2)
        adjacency[n2].add(n1)
    return adjacency, compliance


def _prune_dangling(adjacency: dict[int, set[int]],
                    compliance: dict[tuple[int, int], float]) -> None:
    """Strip degree-1 chains; open flanges carry no closed-cell shear flow."""
    stack = [n for n, nbrs in adjacency.items() if len(nbrs) <= 1]
    while stack:
        n = stack.pop()
        nbrs = adjacency.get(n)
        if nbrs is None or len(nbrs) > 1:
            continue
        for m in list(nbrs):
            adjacency[m].discard(n)
            compliance.pop((n, m) if n < m else (m, n), None)
            if len(adjacency[m]) <= 1:
                stack.append(m)
        adjacency.pop(n, None)


def _find_faces(adjacency: dict[int, set[int]],
                xy: dict[int, np.ndarray]) -> list[list[int]]:
    """Enumerate the faces of the planar graph by half-edge traversal.

    At a node ``v`` reached from ``u``, the next edge of the face is the
    neighbor immediately clockwise from ``u`` in angular order.  This walks
    interior faces counter-clockwise (positive signed area) and the single
    outer face clockwise (negative signed area).
    """
    # angular order of the neighbors of each node
    order: dict[int, list[int]] = {}
    index: dict[int, dict[int, int]] = {}
    for n, nbrs in adjacency.items():
        p = xy[n]
        nbrs_sorted = sorted(
            nbrs, key=lambda m: np.arctan2(xy[m][1] - p[1], xy[m][0] - p[0]))
        order[n] = nbrs_sorted
        index[n] = {m: i for i, m in enumerate(nbrs_sorted)}

    faces: list[list[int]] = []
    visited: set[tuple[int, int]] = set()
    for u0, nbrs in adjacency.items():
        for v0 in nbrs:
            if (u0, v0) in visited:
                continue
            face = []
            u, v = u0, v0
            while (u, v) not in visited:
                visited.add((u, v))
                face.append(u)
                nbrs_v = order[v]
                # step clockwise from the edge v->u
                i = index[v][u]
                w = nbrs_v[(i - 1) % len(nbrs_v)]
                u, v = v, w
                if len(face) > len(visited) + 4:  # pragma: no cover
                    break
            if len(face) >= 3:
                faces.append(face)
    return faces


def _signed_area(face: list[int], xy: dict[int, np.ndarray]) -> float:
    pts = np.array([xy[n] for n in face], dtype='float64')
    x, z = pts[:, 0], pts[:, 1]
    return 0.5 * float(np.dot(x, np.roll(z, -1)) - np.dot(z, np.roll(x, -1)))


def bredt_batho_gj(xyz1: np.ndarray,
                   xyz2: np.ndarray,
                   length: np.ndarray,
                   thickness: np.ndarray,
                   gxy: np.ndarray,
                   iaxes: tuple[int, int]=(0, 2),
                   weld_tol: float=1e-4,
                   log: Optional[Any]=None) -> tuple[float, str, int]:
    """
    Torsional stiffness GJ of a thin-walled section.

    Parameters
    ----------
    xyz1, xyz2 : (nwall, 3) float ndarray
        end coordinates of each wall in the cut plane's local frame; the
        section topology is recovered by welding coincident ends
    length, thickness, gxy : (nwall,) float ndarray
        wall length, wall thickness and wall shear modulus
    iaxes : (int, int); default=(0, 2)
        the two in-plane columns of xyz (the cut plane is local x-z)
    weld_tol : float; default=1e-4
        node welding tolerance as a fraction of the median wall length
    log : logger; optional

    Returns
    -------
    GJ : float
        the torsional stiffness (already multiplied by G)
    method : str
        'closed', 'open' or 'none'
    ncells : int
        number of closed cells found

    """
    ix, iz = iaxes
    keep = ((thickness > 0.) & (gxy > 0.) & (length > 0.) &
            np.isfinite(thickness) & np.isfinite(gxy) & np.isfinite(length))
    if not keep.any():
        return 0., 'none', 0

    length = length[keep]
    thickness = thickness[keep]
    gxy = gxy[keep]
    p1 = xyz1[keep][:, [ix, iz]]
    p2 = xyz2[keep][:, [ix, iz]]

    # open-section value; also the fallback when nothing closes
    gj_open = float((gxy * length * thickness ** 3).sum() / 3.)

    tol = weld_tol * float(np.median(length))
    if tol <= 0.:  # pragma: no cover
        return gj_open, 'open', 0
    nid1, nid2, xy = _weld(p1, p2, tol)

    delta = length / (gxy * thickness)  # wall compliance, s/(G*t)
    adjacency, compliance = _build_graph(nid1, nid2, xy, delta)
    _prune_dangling(adjacency, compliance)
    if not compliance:
        return gj_open, 'open', 0

    faces = _find_faces(adjacency, xy)
    cells = [f for f in faces if _signed_area(f, xy) > 0.]
    if not cells:
        return gj_open, 'open', 0

    ncell = len(cells)
    cell_area = np.array([_signed_area(f, xy) for f in cells])

    # which cells touch each wall
    edge_cells: dict[tuple[int, int], list[int]] = defaultdict(list)
    cell_edges: list[list[tuple[int, int]]] = []
    for icell, face in enumerate(cells):
        edges = []
        for a, b in zip(face, face[1:] + face[:1]):
            key = (a, b) if a < b else (b, a)
            if key not in compliance:
                continue
            edges.append(key)
            edge_cells[key].append(icell)
        cell_edges.append(edges)

    # unit twist rate: (1/2A_i)*sum_walls (q_i - q_j)*delta = 1
    M = np.zeros((ncell, ncell), dtype='float64')
    rhs = np.ones(ncell, dtype='float64')
    for icell, edges in enumerate(cell_edges):
        scale = 1. / (2. * cell_area[icell])
        for key in edges:
            di = compliance[key]
            M[icell, icell] += scale * di
            for jcell in edge_cells[key]:
                if jcell != icell:
                    M[icell, jcell] -= scale * di

    try:
        q = np.linalg.solve(M, rhs)
    except np.linalg.LinAlgError:  # pragma: no cover
        if log is not None:
            log.warning('singular multi-cell torsion system; using open-section GJ')
        return gj_open, 'open', 0

    gj_closed = float(2. * (q * cell_area).sum())
    if not np.isfinite(gj_closed) or gj_closed <= 0.:  # pragma: no cover
        return gj_open, 'open', 0

    # the open-section term is a (small) additive contribution
    return gj_closed + gj_open, 'closed', ncell
