"""
Shear center of a thin-walled section defined by a soup of wall segments.

The shear center is the point through which a transverse shear can be applied
without twisting the section.  It is what MSC's CBEAM element axis is supposed
to pass through: ``GA + WA`` lies on the shear-center axis, and ``N1/N2`` then
locate the neutral axis relative to it (MSC Figure 16-142).  For a doubly
symmetric section the shear center, the neutral axis and the area centroid all
coincide and none of this matters; for an airfoil box they are several inches
apart and putting the beam axis in the wrong place converts bending into
torsion and back.

The section that comes out of a cutting plane is an unordered set of straight
segments, one per cut shell element, each with a thickness, a Young's modulus
and a shear modulus.  ``torsion.py`` already recovers the topology of that soup
to solve for the *torsional* shear flow; this module reuses the same welding
and graph machinery to solve for the *transverse* shear flow, which is a
different problem:

- torsion: constant flow around each cell, unknown per cell, driven by a unit
  twist rate
- transverse shear: the flow varies along each wall, driven by the axial stress
  gradient that accompanies a shear force, and is fixed by wall equilibrium
  plus the requirement that the section does not twist

Method
------
Work in the cut plane with in-plane coordinates ``(u, v)`` measured from the
modulus-weighted centroid, so that ``sum(E*u dA) = sum(E*v dA) = 0``.  For a
shear force ``(Vu, Vv)`` the axial stress gradient is::

    dsigma/dx = -E*(a*u + b*v)      [Euu Euv; Euv Evv] [a; b] = [Vu; Vv]

with ``Euu = int(E*u^2 dA)`` etc.  Wall equilibrium ``dq/ds + t*dsigma/dx = 0``
then gives the flow gradient along each wall::

    dq/ds = -E*t*(a*u + b*v) =: g(s)

The minus sign is the one that makes the flow come out in the same direction as
the applied shear.  Sanity check on a solid rectangle in pure vertical shear:
the flow has to vanish at the top and bottom fibers and peak at the neutral
axis, so ``dq/ds`` must be positive below the neutral axis and negative above
it, i.e. proportional to ``-v``.  Getting this backwards is invisible in the
shear center itself -- negating the whole flow field negates the moment too --
which is exactly why ``force_error`` compares against the *signed* unit shear.

``g`` is linear in ``s`` because ``u`` and ``v`` are, so the flow entering a
wall (``Q``) determines everything about it in closed form::

    q(s)         = Q + g0*s + (g1 - g0)*s^2/(2L)
    q(L)         = Q + G,          G  = L*(g0 + g1)/2
    int(q ds)    = Q*L + L^2*(2*g0 + g1)/6

The unknowns are the ``nwall`` entry flows ``Q``.  They are fixed by

1. flow conservation at every welded junction (this is what forces ``q = 0`` at
   a free edge, with no special-casing of open sections), and
2. zero twist around each independent cycle, ``sum(q/(G*t) ds) = 0`` -- which is
   the *definition* of the shear center.

Those two sets together are exactly determined: ``nwall`` unknowns against
``nnode - ncomponent`` independent junction equations plus
``nwall - nnode + ncomponent`` cycles.  Open, branched, single-cell and
multi-cell sections are all handled by the same assembly; the only difference
is how many cycles there are.

Finally, the resultant of the computed flow is a force ``(Vu, Vv)`` acting on
some line; its moment about the centroid places that line.  Two load cases
(``Vu = 1`` and ``Vv = 1``) give the two shear-center coordinates.  The force
resultant is *not* imposed anywhere in the solve, so comparing it against the
applied unit shear is a free end-to-end check that the topology, the welding
and the linear algebra all came out right; ``force_error`` reports it and the
caller should not trust a result with a large one.

"""
from __future__ import annotations
from collections import defaultdict
from typing import Any, NamedTuple, Optional

import numpy as np

from pyNastran.bdf.mesh_utils.cut.torsion import _weld


class ShearCenterResult(NamedTuple):
    """
    Parameters
    ----------
    xy_shear_center : (2,) float ndarray
        the shear center, in the cut plane's in-plane coordinates (the two
        ``iaxes`` columns of the local frame).  NaN if it could not be found.
    xy_neutral_axis : (2,) float ndarray
        the modulus-weighted centroid, in the same coordinates.  This is the
        neutral axis of the section, which is NOT the area centroid unless the
        section is homogeneous.
    method : str
        'closed', 'open' or 'none'
    ncells : int
        number of independent cycles in the wall graph
    force_error : float
        how far the resultant of the computed shear flow is from the applied
        unit shear, as a fraction of it.  A few times 1e-12 is normal; anything
        approaching 1 means the answer is meaningless.

    """
    xy_shear_center: np.ndarray
    xy_neutral_axis: np.ndarray
    method: str
    ncells: int
    force_error: float


def _fundamental_cycles(edges: list[tuple[int, int]],
                        nodes: list[int],
                        ) -> tuple[list[list[tuple[int, int]]], int]:
    """
    A cycle basis of the wall graph, and the number of connected components.

    Builds a spanning forest; every edge left out of it closes exactly one
    independent cycle, found by walking both endpoints up to their common
    ancestor.  Returns the cycles as lists of ``(iedge, sign)``, where the sign
    is +1 if the cycle traverses that wall from node 1 to node 2.
    """
    adjacency: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for iedge, (n1, n2) in enumerate(edges):
        adjacency[n1].append((n2, iedge))
        adjacency[n2].append((n1, iedge))

    parent: dict[int, tuple[int, int]] = {}   # node -> (parent node, iedge)
    depth: dict[int, int] = {}
    tree_edges: set[int] = set()
    ncomponent = 0
    for root in nodes:
        if root in depth:
            continue
        ncomponent += 1
        depth[root] = 0
        stack = [root]
        while stack:
            n = stack.pop()
            for m, iedge in adjacency[n]:
                if m in depth:
                    continue
                depth[m] = depth[n] + 1
                parent[m] = (n, iedge)
                tree_edges.add(iedge)
                stack.append(m)

    def _step_up(n: int) -> tuple[int, tuple[int, int]]:
        p, iedge = parent[n]
        # +1 means "traversed n1 -> n2"; going from n up to p uses the wall
        # backwards when the wall was stored as p -> n
        sign = 1 if edges[iedge][0] == n else -1
        return p, (iedge, sign)

    cycles: list[list[tuple[int, int]]] = []
    for iedge, (n1, n2) in enumerate(edges):
        if iedge in tree_edges:
            continue
        # walk both ends up to their common ancestor
        up1: list[tuple[int, int]] = []
        up2: list[tuple[int, int]] = []
        a, b = n1, n2
        while depth[a] > depth[b]:
            a, step = _step_up(a)
            up1.append(step)
        while depth[b] > depth[a]:
            b, step = _step_up(b)
            up2.append(step)
        while a != b:
            a, step = _step_up(a)
            up1.append(step)
            b, step = _step_up(b)
            up2.append(step)
        # the cycle runs n2 -> ancestor (up2 as walked), ancestor -> n1 (up1
        # walked backwards, so each step is taken against its stored sense),
        # then n1 -> n2 along the wall that closed it
        cycle = list(up2)
        cycle += [(ie, -sign) for ie, sign in reversed(up1)]
        cycle.append((iedge, 1))
        cycles.append(cycle)
    return cycles, ncomponent


def shear_center(xyz1: np.ndarray,
                 xyz2: np.ndarray,
                 length: np.ndarray,
                 thickness: np.ndarray,
                 ex: np.ndarray,
                 gxy: np.ndarray,
                 iaxes: tuple[int, int]=(0, 2),
                 weld_tol: float=1e-4,
                 log: Optional[Any]=None) -> ShearCenterResult:
    """
    Shear center and neutral axis of a thin-walled section.

    Parameters
    ----------
    xyz1, xyz2 : (nwall, 3) float ndarray
        end coordinates of each wall in the cut plane's local frame; the
        section topology is recovered by welding coincident ends
    length, thickness, ex, gxy : (nwall,) float ndarray
        wall length, wall thickness, wall Young's modulus and wall shear
        modulus.  ``ex`` weights the bending problem and ``gxy`` the twist
        constraint; a homogeneous section is the special case where both are
        constant and neither matters.
    iaxes : (int, int); default=(0, 2)
        the two in-plane columns of xyz (the cut plane is local x-z)
    weld_tol : float; default=1e-4
        node welding tolerance as a fraction of the median wall length; must
        match ``bredt_batho_gj`` so the two see the same topology
    log : logger; optional

    Returns
    -------
    result : ShearCenterResult

    Notes
    -----
    The wall's own second moment about its midline is neglected (the
    ``t^3/12`` term), which is the usual thin-wall assumption and is the same
    approximation ``bredt_batho_gj`` makes.  The variation *along* each wall is
    not neglected -- the strip integrals are exact for a straight wall.

    """
    ix, iz = iaxes
    nan2 = np.full(2, np.nan, dtype='float64')
    keep = ((thickness > 0.) & (gxy > 0.) & (length > 0.) & (ex > 0.) &
            np.isfinite(thickness) & np.isfinite(gxy) &
            np.isfinite(length) & np.isfinite(ex))
    if not keep.any():
        return ShearCenterResult(nan2, nan2, 'none', 0, np.nan)

    length = length[keep]
    thickness = thickness[keep]
    ex = ex[keep]
    gxy = gxy[keep]
    p1 = xyz1[keep][:, [ix, iz]]
    p2 = xyz2[keep][:, [ix, iz]]

    # E-weighted wall "area".  The modulus-weighted centroid is the neutral
    # axis; the area centroid is not, unless E is uniform.
    ea = ex * length * thickness
    total_ea = ea.sum()
    if total_ea <= 0.:  # pragma: no cover
        return ShearCenterResult(nan2, nan2, 'none', 0, np.nan)
    xy_na = (ea[:, np.newaxis] * 0.5 * (p1 + p2)).sum(axis=0) / total_ea

    u1, v1 = (p1 - xy_na).T
    u2, v2 = (p2 - xy_na).T

    # exact strip integrals for a straight wall of uniform E*t
    euu = float((ea * (u1 * u1 + u1 * u2 + u2 * u2) / 3.).sum())
    evv = float((ea * (v1 * v1 + v1 * v2 + v2 * v2) / 3.).sum())
    euv = float((ea * (2. * u1 * v1 + u1 * v2 + u2 * v1 + 2. * u2 * v2) / 6.).sum())
    det = euu * evv - euv * euv
    if not np.isfinite(det) or det <= 0.:
        # every wall on one line: the section has no bending stiffness about
        # one axis and the shear center is not defined
        if log is not None:
            log.warning('degenerate section (walls are collinear); '
                        'no shear center')
        return ShearCenterResult(nan2, xy_na, 'none', 0, np.nan)

    tol = weld_tol * float(np.median(length))
    if tol <= 0.:  # pragma: no cover
        return ShearCenterResult(nan2, xy_na, 'none', 0, np.nan)
    nid1, nid2, unused_xy = _weld(p1, p2, tol)

    # a wall whose ends welded together carries no flow gradient anywhere
    live = nid1 != nid2
    if not live.any():  # pragma: no cover
        return ShearCenterResult(nan2, xy_na, 'none', 0, np.nan)
    nid1 = nid1[live]
    nid2 = nid2[live]
    length = length[live]
    thickness = thickness[live]
    ex = ex[live]
    gxy = gxy[live]
    ea = ea[live]
    u1, v1, u2, v2 = u1[live], v1[live], u2[live], v2[live]

    edges = list(zip(nid1.tolist(), nid2.tolist()))
    nodes = sorted(set(nid1.tolist()) | set(nid2.tolist()))
    cycles, ncomponent = _fundamental_cycles(edges, nodes)
    ncell = len(cycles)
    if ncomponent > 1 and log is not None:
        log.warning(
            f'the section has {ncomponent:d} disconnected pieces; the shear '
            'center of the assembly is not the shear center of any one of '
            'them and this answer is probably not meaningful')

    nwall = len(edges)
    inode = {n: i for i, n in enumerate(nodes)}
    nnode = len(nodes)

    # --- assembly -------------------------------------------------------
    # rows 0..nnode-1        junction flow conservation
    # rows nnode..nnode+ncell-1   zero twist around each independent cycle
    nrow = nnode + ncell
    matrix = np.zeros((nrow, nwall), dtype='float64')
    for iedge, (n1, n2) in enumerate(edges):
        matrix[inode[n1], iedge] -= 1.   # leaves n1 carrying Q
        matrix[inode[n2], iedge] += 1.   # arrives at n2 carrying Q (+G, in rhs)

    # int(q ds) over a wall = Q*L + L^2*(2*g0 + g1)/6, so the twist row picks
    # up L/(G*t) against the unknown and the rest goes to the rhs
    gt = gxy * thickness
    for icycle, cycle in enumerate(cycles):
        row = nnode + icycle
        for iedge, sign in cycle:
            matrix[row, iedge] += sign * length[iedge] / gt[iedge]

    # --- two unit load cases --------------------------------------------
    rhs = np.zeros((nrow, 2), dtype='float64')
    g0s = np.zeros((nwall, 2), dtype='float64')
    g1s = np.zeros((nwall, 2), dtype='float64')
    et = ex * thickness
    for icase, (vu, vv) in enumerate(((1., 0.), (0., 1.))):
        a = (evv * vu - euv * vv) / det
        b = (euu * vv - euv * vu) / det
        # dq/ds = -E*t*(a*u + b*v); see the module docstring for the sign
        g0 = -et * (a * u1 + b * v1)
        g1 = -et * (a * u2 + b * v2)
        g0s[:, icase] = g0
        g1s[:, icase] = g1
        gain = length * (g0 + g1) / 2.            # G, the flow picked up
        for iedge, (unused_n1, n2) in enumerate(edges):
            rhs[inode[n2], icase] -= gain[iedge]
        moment_rest = length ** 2 * (2. * g0 + g1) / 6.
        for icycle, cycle in enumerate(cycles):
            row = nnode + icycle
            for iedge, sign in cycle:
                rhs[row, icase] -= sign * moment_rest[iedge] / gt[iedge]

    # one junction equation per component is redundant (the E-weighted first
    # moment vanishes at the centroid), so this is least-squares by shape but
    # exactly determined by rank
    qentry, unused_res, rank, unused_sv = np.linalg.lstsq(matrix, rhs, rcond=None)
    if rank < nwall:
        if log is not None:
            log.warning(
                f'shear-flow system is rank {rank:d} of {nwall:d}; the wall '
                'graph is probably broken (check the weld tolerance)')
        return ShearCenterResult(nan2, xy_na, 'none', ncell, np.nan)

    # --- resultants ------------------------------------------------------
    # a straight wall has a constant moment arm about the origin, because
    # r(s) x t_hat = r0 x t_hat for r(s) = r0 + s*t_hat
    arm = (u1 * v2 - u2 * v1) / length
    tu = (u2 - u1) / length
    tv = (v2 - v1) / length

    xy_sc = np.empty(2, dtype='float64')
    force_error = 0.
    for icase, (vu, vv) in enumerate(((1., 0.), (0., 1.))):
        qint = (qentry[:, icase] * length +
                length ** 2 * (2. * g0s[:, icase] + g1s[:, icase]) / 6.)
        fu = float((tu * qint).sum())
        fv = float((tv * qint).sum())
        moment = float((arm * qint).sum())
        force_error = max(force_error, abs(fu - vu), abs(fv - vv))
        if icase == 0:
            # V = (1, 0) applied at (su, sv) has moment -sv about the origin
            xy_sc[1] = -moment
        else:
            # V = (0, 1) applied at (su, sv) has moment +su
            xy_sc[0] = moment

    if not np.isfinite(xy_sc).all() or force_error > 1e-6:
        if log is not None:
            log.warning(
                f'shear-flow resultant is off by {force_error:g} of the '
                'applied unit shear; discarding the shear center')
        return ShearCenterResult(nan2, xy_na, 'none', ncell, force_error)

    method = 'closed' if ncell else 'open'
    return ShearCenterResult(xy_na + xy_sc, xy_na, method, ncell, force_error)
