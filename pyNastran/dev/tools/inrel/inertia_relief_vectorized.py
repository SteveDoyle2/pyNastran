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
- Loads are passed in as ``(ngrid,3)`` basic-frame arrays.  Assembling them
  from a load case (FORCE/MOMENT/PLOAD4/GRAV...) is the caller's job.
  TODO: accept a load id and build the loads from the model.
- Mass coupling between grids (e.g. CMASS2 between two grids) is ignored and
  logged as a warning; see ``inertia_relief_6x6``.

"""

from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse

from pyNastran.dev.bdf_vectorized3.solver.matrices.mass_matrix import build_mgg_lumped
from pyNastran.dev.tools.inrel.inertia_relief2 import inertia_relief_6x6

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def get_mass_blocks(mbb: sparse.spmatrix) -> tuple[np.ndarray, float]:
    """
    Extracts the diagonal 6x6 block of each grid from a g-set mass matrix.

    Parameters
    ----------
    mbb : (6*ngrid, 6*ngrid) sparse matrix
        the mass matrix with 6 DOFs per grid in [T1, T2, T3, R1, R2, R3] order

    Returns
    -------
    mass_blocks : (ngrid,6,6) float ndarray
        the per-grid diagonal blocks
    coupling_ratio : float
        ``sum(|off-block terms|) / sum(|all terms|)``; 0.0 for a purely
        block-diagonal (lumped) matrix

    """
    ndof = mbb.shape[0]
    assert ndof % 6 == 0, ndof
    ngrid = ndof // 6
    coo = sparse.coo_matrix(mbb)
    igrid = coo.row // 6
    in_block = igrid == coo.col // 6

    mass_blocks = np.zeros((ngrid, 6, 6), dtype="float64")
    np.add.at(
        mass_blocks,
        (igrid[in_block], coo.row[in_block] % 6, coo.col[in_block] % 6),
        coo.data[in_block],
    )
    abs_total = np.abs(coo.data).sum()
    coupling_ratio = np.abs(coo.data[~in_block]).sum() / abs_total if abs_total > 0.0 else 0.0
    return mass_blocks, float(coupling_ratio)


def inertia_relief_vectorized(
    model: BDF,
    force: np.ndarray,
    moment: np.ndarray,
    wtmass: float = 1.0,
    debug: bool = False,
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
        PARAM,WTMASS
    debug : bool; default=False
        passed to ``inertia_relief``

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
