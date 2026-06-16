"""Inertia relief for free-free structures (vectorized).

Computes self-equilibrating loads for an unrestrained structure by removing
the rigid body force components via d'Alembert's principle. The structure
accelerates as a rigid body such that the net applied + inertial forces
produce zero resultant force and moment.

Theory
------
For a free-free structure with applied loads F and mass matrix M:

1. Rigid body mode matrix D maps 6 rigid body DOFs (3 translations +
   3 rotations about a reference point) to all g-set DOFs.

2. Rigid body acceleration that equilibrates the applied loads:
       a_r = M_rr^{-1} @ D^T @ F
   where M_rr = D^T @ M @ D is the 6x6 rigid body mass matrix.

3. Inertial loads (d'Alembert forces):
       F_inertia = -M @ D @ a_r

4. Self-equilibrated (net) loads:
       F_net = F + F_inertia

   These satisfy D^T @ F_net = 0 (zero net force and moment), so the
   structure can be solved with a constrained (SUPORT) stiffness matrix:
       K_ll @ u_l = F_net_l

The elastic deformation u represents the deformation relative to the
rigid body motion. The total displacement is: u_total = u_elastic + D @ x_r
where x_r is the rigid body displacement (arbitrary for inertia relief).

Rank-Deficient Geometry
~~~~~~~~~~~~~~~~~~~~~~~
For degenerate geometries (e.g. all grids collinear), M_rr will be
rank-deficient (< 6). The module handles this via least-squares
(minimum-norm pseudo-inverse) rather than raising an error, since the
self-equilibration D^T @ F_net = 0 still holds for the non-degenerate
DOFs. Only a truly zero mass matrix (rank 0) is an error.

Usage with Nastran SUPORT
-------------------------
In SOL 101 with INREL (inertia relief):
1. Build M (lumped or consistent) via build_mgg_lumped().
2. Compute D from grid geometry and SUPORT reference point.
3. Apply inertia relief to get self-equilibrated loads.
4. Partition to l-set (remove SUPORT DOFs) and solve.

The reference point for D determines how the rigid body acceleration
is expressed (but F_net is invariant to this choice).

Example
-------
>>> from pyNastran.dev.bdf_vectorized3.bdf import BDF
>>> from pyNastran.dev.bdf_vectorized3.solver.matrices.mass_matrix import build_mgg_lumped
>>> from pyNastran.dev.bdf_vectorized3.mesh_utils.inertia_relief import (
...     build_rigid_body_modes, compute_inertia_relief)
>>> import numpy as np
>>>
>>> # Load model and build mass matrix
>>> model = BDF()
>>> model.read_bdf('my_model.bdf')
>>> M, grid_ids, total_mass = build_mgg_lumped(model=model)
>>>
>>> # Build rigid body modes (ref at origin)
>>> node_xyz = model.grid.xyz_cid0()
>>> D = build_rigid_body_modes(grid_ids, node_xyz)
>>>
>>> # Apply gravity in -Z to all grids
>>> n_dof = 6 * len(grid_ids)
>>> F = np.zeros(n_dof)
>>> g = 9.81
>>> diag = M.diagonal()
>>> for i in range(len(grid_ids)):
...     F[6*i + 2] = -diag[6*i] * g  # Fz = -m_i * g
>>>
>>> # Compute self-equilibrated loads
>>> F_net, a_rigid, F_inertia = compute_inertia_relief(M, D, F)
>>> # F_net can now be solved with K_ll (SUPORT-constrained stiffness)

Or use the convenience wrapper:

>>> from pyNastran.dev.bdf_vectorized3.mesh_utils.inertia_relief import (
...     compute_inertia_relief_from_model)
>>> F_net, a_rigid, F_inertia, info = compute_inertia_relief_from_model(
...     model, F, ref_point=np.array([0., 0., 0.]))
>>> print(f"Total mass: {info['total_mass']:.1f} kg")
>>> print(f"Rigid body accel: {a_rigid}")
>>> print(f"Self-eq residual: {np.max(np.abs(info['residual_force'])):.2e}")

See Also
--------
pyNastran.dev.bdf_vectorized3.solver.matrices.mass_matrix : MGG assembly module.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray
from scipy import sparse
from cpylog import SimpleLogger

from pyNastran.dev.bdf_vectorized3.solver.matrices.mass_matrix import build_mgg_lumped
if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

__all__ = [
    'build_rigid_body_modes',
    'compute_inertia_relief',
    'compute_inertia_relief_from_model',
]

log = SimpleLogger(level='warning')


def build_rigid_body_modes(
    grid_ids: NDArray[np.integer],
    node_xyz: NDArray[np.floating],
    ref_point: NDArray[np.floating] | None = None,
) -> NDArray[np.floating]:
    """Build the rigid body mode matrix D (vectorized, no loops).

    D maps 6 rigid body DOFs (Tx, Ty, Tz, Rx, Ry, Rz at the reference
    point) to all g-set DOFs. Each column is one rigid body mode.

    Parameters
    ----------
    grid_ids : (n_grids,) int array
        Sorted grid IDs defining DOF ordering.
    node_xyz : (n_grids, 3) float array
        Grid positions in the basic coordinate system.
    ref_point : (3,) array or None
        Reference point for rigid body modes. If None, uses [0, 0, 0].
        The choice of reference point affects where the rigid body
        acceleration is defined but not the resulting F_net.

    Returns
    -------
    D : (n_dof, 6) ndarray
        Rigid body mode matrix. n_dof = 6 * n_grids.
        Columns: [Tx, Ty, Tz, Rx, Ry, Rz].
    """
    if ref_point is None:
        ref_point = np.zeros(3)
    ref_point = np.asarray(ref_point, dtype=np.float64)

    n_grids = len(grid_ids)
    n_dof = 6 * n_grids

    arms = node_xyz - ref_point[np.newaxis, :]
    dx = arms[:, 0]
    dy = arms[:, 1]
    dz = arms[:, 2]

    D = np.zeros((n_dof, 6), dtype=np.float64)
    base = np.arange(n_grids, dtype=np.int64) * 6

    # Translation modes (columns 0, 1, 2)
    D[base, 0] = 1.0
    D[base + 1, 1] = 1.0
    D[base + 2, 2] = 1.0

    # Rotation modes (columns 3, 4, 5): velocity = omega x r
    # For unit rotation omega_k, translation = omega_k_hat x arm_vector.
    # Rx (omega=[1,0,0]): v = [0, -dz, +dy]
    D[base + 1, 3] = -dz
    D[base + 2, 3] = dy
    # Ry (omega=[0,1,0]): v = [+dz, 0, -dx]
    D[base, 4] = dz
    D[base + 2, 4] = -dx
    # Rz (omega=[0,0,1]): v = [-dy, +dx, 0]
    D[base, 5] = -dy
    D[base + 1, 5] = dx

    # Rotational DOFs: identity for rotation modes
    D[base + 3, 3] = 1.0
    D[base + 4, 4] = 1.0
    D[base + 5, 5] = 1.0

    return D


def compute_inertia_relief(
    M: sparse.spmatrix | NDArray,
    D: NDArray[np.floating],
    F_applied: NDArray[np.floating],
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """Compute inertia relief loads for a free-free structure.

    Removes the rigid body force component from the applied loads by
    adding d'Alembert inertial forces, producing self-equilibrated loads
    that can be solved with a constrained stiffness matrix.

    Parameters
    ----------
    M : sparse or dense matrix, shape (n_dof, n_dof)
        Mass matrix (Mbb or Mgg, must match DOF ordering of D and F).
    D : (n_dof, 6) ndarray
        Rigid body mode matrix from build_rigid_body_modes().
    F_applied : (n_dof,) or (n_dof, n_load) ndarray
        Applied load vector(s). Multiple load cases supported as columns.

    Returns
    -------
    F_net : same shape as F_applied
        Self-equilibrated loads (F_applied + F_inertia).
        Satisfies D^T @ F_net = 0 (zero net force/moment).
    a_rigid : (6,) or (6, n_load) ndarray
        Rigid body acceleration [ax, ay, az, alpha_x, alpha_y, alpha_z]
        at the reference point used to build D.
    F_inertia : same shape as F_applied
        Inertial (d'Alembert) forces: F_inertia = -M @ D @ a_rigid.

    Raises
    ------
    ValueError
        If the rigid body mass matrix is singular (e.g., zero mass or
        degenerate geometry).

    Notes
    -----
    The computation is:
        M_rr = D^T @ M @ D           (6x6 rigid body mass matrix)
        a_r = M_rr^{-1} @ D^T @ F    (rigid body acceleration)
        F_inertia = -M @ D @ a_r      (d'Alembert forces)
        F_net = F + F_inertia          (self-equilibrated)

    For multiple load cases, each column of F_applied is processed
    independently (same M_rr factorization reused).
    """
    F_applied = np.asarray(F_applied)
    single_load = F_applied.ndim == 1
    if single_load:
        F_applied = F_applied[:, np.newaxis]

    n_dof = F_applied.shape[0]
    if D.shape[0] != n_dof:
        raise ValueError(
            f"D has {D.shape[0]} rows but F_applied has {n_dof} entries")

    # M_rr = D^T @ M @ D  (6x6)
    if sparse.issparse(M):
        MD = M.dot(D)
    else:
        MD = M @ D
    M_rr = D.T @ MD  # (6, 6)

    # Check for zero mass (no inertia at all)
    rank = np.linalg.matrix_rank(M_rr, tol=1e-10)
    if rank == 0:
        raise ValueError(
            "Rigid body mass matrix is zero. "
            "Check that the model has mass assigned to elements/CONM2.")

    # D^T @ F  (6, n_load)
    DtF = D.T @ F_applied

    # a_r = M_rr^{-1} @ D^T @ F (lstsq handles rank-deficient cases,
    # e.g. collinear grids with zero inertia about their axis)
    if rank < 6:
        log.warning(
            f"M_rr is rank {rank} (< 6). Degenerate geometry (e.g. collinear "
            "grids). Using minimum-norm solution for rigid body acceleration.")
        a_rigid, _, _, _ = np.linalg.lstsq(M_rr, DtF, rcond=None)
    else:
        a_rigid = np.linalg.solve(M_rr, DtF)  # (6, n_load)

    # F_inertia = -M @ D @ a_r
    D_a = D @ a_rigid  # (n_dof, n_load)
    if sparse.issparse(M):
        F_inertia = -M.dot(D_a)
    else:
        F_inertia = -M @ D_a

    # F_net = F + F_inertia
    F_net = F_applied + F_inertia

    if single_load:
        return F_net.ravel(), a_rigid.ravel(), F_inertia.ravel()
    return F_net, a_rigid, F_inertia


def compute_inertia_relief_from_model(
    model: 'BDF',
    F_applied: NDArray[np.floating],
    ref_point: NDArray[np.floating] | None = None,
    wtmass: float = 1.0,
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating], dict]:
    """Compute inertia relief directly from a BDF model.

    Convenience function that builds M and D internally, then applies
    inertia relief.

    Parameters
    ----------
    model : BDF
        Pre-loaded bdf_vectorized3 BDF object (xref=True).
    F_applied : (n_dof,) or (n_dof, n_load) ndarray
        Applied load vector(s) in g-set DOF ordering (sorted grid IDs,
        6 DOFs per grid).
    ref_point : (3,) array or None
        Reference point for rigid body modes. If None, uses [0, 0, 0].
    wtmass : float
        PARAM,WTMASS value. Default 1.0.

    Returns
    -------
    F_net : same shape as F_applied
        Self-equilibrated loads.
    a_rigid : (6,) or (6, n_load) ndarray
        Rigid body acceleration at the reference point.
    F_inertia : same shape as F_applied
        d'Alembert inertial forces.
    info : dict
        Diagnostic information:
        - 'M_rr': (6, 6) rigid body mass matrix
        - 'total_mass': float
        - 'grid_ids': sorted grid IDs
        - 'D': rigid body mode matrix
        - 'residual_force': D^T @ F_net (should be ~0)
    """
    M, grid_ids, total_mass = build_mgg_lumped(model=model, wtmass=wtmass)
    node_xyz = model.grid.xyz_cid0()
    D = build_rigid_body_modes(grid_ids, node_xyz, ref_point)

    F_net, a_rigid, F_inertia = compute_inertia_relief(M, D, F_applied)

    # Diagnostics
    if sparse.issparse(M):
        MD = M.dot(D)
    else:
        MD = M @ D
    M_rr = D.T @ MD

    F_net_2d = F_net if F_net.ndim == 2 else F_net[:, np.newaxis]
    residual = D.T @ F_net_2d

    info = {
        'M_rr': M_rr,
        'total_mass': total_mass,
        'grid_ids': grid_ids,
        'D': D,
        'residual_force': residual.ravel() if F_net.ndim == 1 else residual,
    }

    return F_net, a_rigid, F_inertia, info
