"""
Inertia relief on a ``pyNastran.dev.bdf_vectorized3`` model.

This is a thin adapter: it pulls the per-grid lumped 6x6 mass blocks out of
the basic-frame lumped mass matrix (``build_mgg_lumped``) and hands them to the
model-agnostic ``inertia_relief2.inertia_relief_6x6``.  All of the physics
lives in ``inertia_relief2.py``, which has no ``pyNastran.bdf`` imports; a
different model type only needs to supply its own ``(ngrid,6,6)`` blocks.

Assumptions
-----------
- Same as ``inertia_relief2.inertia_relief``: rigid-body, quasi-static, free-free.
- Lumped mass (PARAM,COUPMASS,-1): structural element mass is split equally to
  its nodes, CONM2 offsets and inertias are rotated to basic.

Limitations
-----------
- GRIDs only; SPOINT/EPOINT mass and loads are ignored.
- ``inertia_relief_from_model`` builds the loads with the vectorized solver's
  ``build_Fb_from_loadid``, so it supports the same load cards (FORCE, MOMENT,
  PLOAD1, PLOAD4, GRAV, RFORCE, ...).  SPCD (enforced motion) is not a load
  and is ignored with a warning.
- Mass coupling between grids (e.g. CMASS2 between two grids) is ignored and
  logged as a warning; see ``inertia_relief_6x6``.

"""

from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse

from pyNastran.dev.bdf_vectorized3.solver.loads.build_fb import build_Fb_from_loadid
from pyNastran.dev.bdf_vectorized3.solver.matrices.mass_matrix import build_mgg_lumped
from pyNastran.dev.bdf_vectorized3.solver.utils import get_param
from pyNastran.dev.tools.inrel.inertia_relief2 import (
    get_mass_blocks,
    inertia_relief_6x6,
    inertia_relief_grav_rforce,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def inertia_relief_vectorized(
    model: BDF,
    force: np.ndarray,
    moment: np.ndarray,
    wtmass: float = 1.0,
    debug: bool = False,
    mbb: sparse.spmatrix | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the d'Alembert inertia relief loads for a bdf_vectorized3 model.

    Parameters
    ----------
    model : pyNastran.dev.bdf_vectorized3.bdf.BDF
        the cross-referenced model
    force : (ngrid,3) float ndarray
        applied forces in the basic frame, ordered by ``model.grid.node_id``
    moment : (ngrid,3) float ndarray
        applied moments in the basic frame, ordered by ``model.grid.node_id``
    wtmass : float; default=1.0
        PARAM,WTMASS; ignored if ``mbb`` is given
    debug : bool; default=False
        passed to ``inertia_relief``
    mbb : (6*ngrid, 6*ngrid) sparse matrix; default=None
        the basic-frame lumped mass matrix (already scaled by wtmass);
        built with ``build_mgg_lumped`` if None

    Returns
    -------
    grid_id : (ngrid,) int ndarray
        the grid ids defining the row order
    dforce : (ngrid,3) float ndarray
        the inertial force increment, basic frame
    dmoment : (ngrid,3) float ndarray
        the inertial moment increment, basic frame

    The balanced loads are ``force + dforce`` and ``moment + dmoment``.

    """
    grid_id = model.grid.node_id
    if mbb is None:
        mbb, grid_id, unused_total_mass = build_mgg_lumped(
            model=model, wtmass=wtmass, output_frame="basic"
        )
    assert np.array_equal(grid_id, model.grid.node_id), "grid order mismatch"
    ngrid = len(grid_id)
    assert force.shape == (ngrid, 3), (force.shape, ngrid)
    assert moment.shape == (ngrid, 3), (moment.shape, ngrid)

    mass_blocks, coupling_ratio = get_mass_blocks(mbb)
    if coupling_ratio > 1e-12:
        model.log.warning(
            f"inertia relief ignores mass coupling between grids "
            f"(off-block fraction={coupling_ratio:.3e})"
        )
    xyz = model.grid.xyz_cid0()
    dforce, dmoment = inertia_relief_6x6(mass_blocks, xyz, force, moment, debug=debug)
    return grid_id, dforce, dmoment


def get_grid_loads(
    model: BDF,
    load_id: int,
    mbb: sparse.spmatrix | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Assembles the applied grid forces and moments for a static load set id.

    Parameters
    ----------
    model : pyNastran.dev.bdf_vectorized3.bdf.BDF
        the cross-referenced model
    load_id : int
        the LOAD/FORCE/MOMENT/PLOAD4/GRAV/... set id (Case Control LOAD=...)
    mbb : (6*ngrid, 6*ngrid) sparse matrix; default=None
        basic-frame mass matrix used for GRAV/RFORCE; built if None

    Returns
    -------
    force : (ngrid,3) float ndarray
        applied forces in the basic frame, ordered by ``model.grid.node_id``
    moment : (ngrid,3) float ndarray
        applied moments in the basic frame, ordered by ``model.grid.node_id``

    """
    assert load_id > 0, load_id
    grid_id = model.grid.node_id
    ngrid = len(grid_id)
    ndof = 6 * ngrid
    dof_map = {
        (nid, component): 6 * igrid + component - 1
        for igrid, nid in enumerate(grid_id)
        for component in range(1, 7)
    }
    if mbb is None:
        wtmass = get_param(model, "WTMASS", 1.0)
        mbb, unused_grid_id, unused_total_mass = build_mgg_lumped(
            model=model, wtmass=wtmass, output_frame="basic"
        )
    xg = np.zeros(ndof, dtype="float64")
    sset_b = np.zeros(ndof, dtype="bool")
    fb, unused_mbb, unused_xyz = build_Fb_from_loadid(
        model, 0, dof_map, ndof, xg, sset_b, load_id=load_id, fdtype="float64", Mbb=mbb
    )
    if sset_b.any():
        # SPCD flags its DOFs with nan; enforced motion is not an applied load
        model.log.warning(f"inertia relief ignores SPCD in load_id={load_id}")
        fb[sset_b] = 0.0
    fb = fb.reshape(ngrid, 6)
    return fb[:, :3].copy(), fb[:, 3:].copy()


def inertia_relief_from_model(
    model: BDF,
    load_id: int,
    wtmass: float | None = None,
    debug: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Inertia relieves a static load set on a bdf_vectorized3 model.

    Parameters
    ----------
    model : pyNastran.dev.bdf_vectorized3.bdf.BDF
        the cross-referenced model
    load_id : int
        the static load set id (e.g., Case Control LOAD=load_id)
    wtmass : float; default=None
        PARAM,WTMASS; None reads it from the model (default 1.0)
    debug : bool; default=False
        passed to ``inertia_relief``

    Returns
    -------
    grid_id : (ngrid,) int ndarray
        the grid ids defining the row order
    force : (ngrid,3) float ndarray
        the applied forces, basic frame
    moment : (ngrid,3) float ndarray
        the applied moments, basic frame
    dforce : (ngrid,3) float ndarray
        the inertial force increment, basic frame
    dmoment : (ngrid,3) float ndarray
        the inertial moment increment, basic frame

    The balanced loads are ``force + dforce`` and ``moment + dmoment``.

    """
    if wtmass is None:
        wtmass = get_param(model, "WTMASS", 1.0)
    mbb, unused_grid_id, unused_total_mass = build_mgg_lumped(
        model=model, wtmass=wtmass, output_frame="basic"
    )
    force, moment = get_grid_loads(model, load_id, mbb=mbb)
    grid_id, dforce, dmoment = inertia_relief_vectorized(model, force, moment, debug=debug, mbb=mbb)
    return grid_id, force, moment, dforce, dmoment


def add_inertia_relief_cards(
    model: BDF,
    load_id: int,
    balanced_load_id: int,
    grav_id: int,
    rforce_id: int,
    wtmass: float | None = None,
    rforce_method: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Balances a static load set with one GRAV, one RFORCE, and a LOAD card.

    Instead of writing a FORCE/MOMENT pair per grid, the rigid-body relief is
    expressed as a GRAV acceleration plus an RFORCE angular acceleration
    (A=0, about the basic origin), combined with the applied load::

        LOAD, balanced_load_id, 1.0, 1.0, load_id, 1.0, grav_id, 1.0, rforce_id

    A zero GRAV or RFORCE (e.g., a load with no net force) is omitted.

    See ``inertia_relief2.inertia_relief_grav_rforce`` for how each card's
    load is modeled.  The cards are added to ``model`` and ``model.setup()``
    is re-run.

    Parameters
    ----------
    model : pyNastran.dev.bdf_vectorized3.bdf.BDF
        the cross-referenced model
    load_id : int
        the applied static load set id
    balanced_load_id : int
        the new LOAD id (select with Case Control LOAD=balanced_load_id)
    grav_id : int
        the new GRAV id
    rforce_id : int
        the new RFORCE id
    wtmass : float; default=None
        PARAM,WTMASS; None reads it from the model (default 1.0).  Nastran
        scales the mass for GRAV and RFORCE by WTMASS, so the accelerations
        are in true units.
    rforce_method : int; default=1
        the RFORCE METHOD field.  NX rejects METHOD=2 when any CONM2, shell,
        or beam has an offset (USER FATAL 6584); METHOD=1 is the QRG choice
        for lumped mass with offsets (RFORCE remark 7).

    Returns
    -------
    grav_accel : (3,) float ndarray
        the GRAV acceleration vector, basic frame
    alpha : (3,) float ndarray
        the RFORCE angular acceleration (rad/time^2), basic frame
    residual : (6,) float ndarray
        the unbalanced [force, moment] resultant about the origin (~0)

    Notes
    -----
    - GRAV can't share a SID with FORCE/MOMENT cards (QRG GRAV remark 3),
      which is why the combination goes through a LOAD card.
    - Nastran does not allow a LOAD card to reference another LOAD card, so
      if ``load_id`` is a LOAD card, its components are copied into the new
      LOAD with their scale factors (``S * Si``).

    """
    if wtmass is None:
        wtmass = get_param(model, "WTMASS", 1.0)
    mbb, grid_id, unused_total_mass = build_mgg_lumped(
        model=model, wtmass=wtmass, output_frame="basic"
    )
    assert np.array_equal(grid_id, model.grid.node_id), "grid order mismatch"
    force, moment = get_grid_loads(model, load_id, mbb=mbb)
    xyz = model.grid.xyz_cid0()
    grav_accel, alpha, residual = inertia_relief_grav_rforce(
        mbb, xyz, force, moment, rforce_method=rforce_method
    )

    load = model.load
    if load.n and load_id in load.load_id:
        # LOAD can't reference a LOAD; copy its components
        iload = np.where(load.load_id == load_id)[0][0]
        i0, i1 = load.iload[iload]
        scale_factors = list(load.scale[iload] * load.scale_factors[i0:i1])
        load_ids = list(load.load_ids[i0:i1])
    else:
        scale_factors = [1.0]
        load_ids = [load_id]

    accel_mag = np.linalg.norm(grav_accel)
    if accel_mag > 0.0:
        model.add_grav(grav_id, accel_mag, grav_accel / accel_mag)
        load_ids.append(grav_id)
        scale_factors.append(1.0)
    alpha_mag = np.linalg.norm(alpha)
    if alpha_mag > 0.0:
        # RACC is in rev/time^2 and R is not normalized (QRG RFORCE remark 2)
        model.add_rforce(
            rforce_id,
            0,
            0.0,
            alpha / alpha_mag,
            cid=0,
            method=rforce_method,
            racc=alpha_mag / (2 * np.pi),
        )
        load_ids.append(rforce_id)
        scale_factors.append(1.0)
    model.add_load(balanced_load_id, 1.0, scale_factors, load_ids)
    model.setup()
    return grav_accel, alpha, residual
