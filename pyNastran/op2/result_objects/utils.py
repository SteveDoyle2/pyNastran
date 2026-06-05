from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse

if TYPE_CHECKING:
    from pyNastran.op2.result_objects.table_object import RealTableArray


def modal_kinetic_energy_fraction(
    eigenvectors: RealTableArray,
    mgg: sparse.spmatrix | np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute fractional modal kinetic energy per grid point.

    For each mode, computes what fraction of the total modal kinetic energy
    is contributed by each grid point. This identifies which grids dominate
    each mode's motion.

    Parameters
    ----------
    eigenvectors : RealTableArray (e.g., op2.eigenvectors[subcase])
        Mode shape data. Shape of .data is (nmodes, nnodes, 6).
        Assumes mass-normalized modes (phi^T @ M @ phi = I).
    mgg : sparse matrix or ndarray, shape (ndof, ndof)
        Global mass matrix (MGG). ndof = nnodes * 6.

    Returns
    -------
    mke_fraction : ndarray, shape (nmodes, nnodes)
        Fractional kinetic energy per grid per mode. Sums to ~1.0 per mode.
    node_ids : ndarray, shape (nnodes,)
        Grid point IDs corresponding to columns of mke_fraction.
    modes : ndarray, shape (nmodes,)
        Mode numbers.

    Notes
    -----
    For mode i, the kinetic energy fraction at grid g is::

        MKE_g = phi_g^T @ M_gg_block @ phi_g / (phi^T @ M @ phi)

    where phi_g is the 6-DOF eigenvector at grid g, and M_gg_block is
    the 6x6 diagonal block of MGG for that grid.

    For mass-normalized modes, phi^T @ M @ phi = 1, so::

        MKE_g = phi_g^T @ M_gg_block @ phi_g

    However this function does NOT assume mass normalization — it
    normalizes by the total phi^T @ M @ phi for each mode.

    """
    phi_data = eigenvectors.data  # (nmodes, nnodes, 6)
    nmodes, nnodes, ndof_per_node = phi_data.shape
    assert ndof_per_node == 6, f'expected 6 DOF/node, got {ndof_per_node}'

    ndof = nnodes * 6
    if sparse.issparse(mgg):
        assert mgg.shape == (ndof, ndof), f'MGG shape {mgg.shape} != ({ndof}, {ndof})'
    else:
        assert mgg.shape == (ndof, ndof), f'MGG shape {mgg.shape} != ({ndof}, {ndof})'

    mke_fraction = np.zeros((nmodes, nnodes), dtype=np.float64)

    for i in range(nmodes):
        # Flatten mode shape to (ndof,) vector: [node1_T1..T6, node2_T1..T6, ...]
        phi = phi_data[i].ravel()  # (ndof,)

        # Total modal mass: phi^T @ M @ phi
        if sparse.issparse(mgg):
            m_phi = mgg.dot(phi)
        else:
            m_phi = mgg @ phi
        total_modal_mass = phi @ m_phi

        if total_modal_mass == 0.0:
            continue

        # Per-grid contribution: phi_g^T @ M_gg @ phi (includes coupling)
        # Use the full M @ phi product, then dot with phi per grid
        for g in range(nnodes):
            dof_start = g * 6
            dof_end = dof_start + 6
            phi_g = phi[dof_start:dof_end]
            m_phi_g = m_phi[dof_start:dof_end]
            mke_fraction[i, g] = phi_g @ m_phi_g / total_modal_mass

    node_ids = eigenvectors.node_gridtype[:, 0]
    modes = np.array(eigenvectors.modes)
    return mke_fraction, node_ids, modes


def real_modes_to_omega_freq(eigns: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    omega is also called radians
    cycles is also called frequency
    """
    omega_radians = np.sqrt(np.abs(eigns))
    abs_freqs = omega_radians / (2 * np.pi)
    return omega_radians, abs_freqs


def complex_damping_frequency(eigr: np.ndarray,
                              eigi: np.ndarray,
                              cref: float=1.0, velocity: float=1.0) -> tuple[np.ndarray, np.ndarray]:
    """eigenvalue = eigr + eigi*1j"""
    assert isinstance(eigr, np.ndarray), eigr
    eigr[eigr == -0.] = 0.
    eigi[eigi == -0.] = 0.
    damping = np.zeros(len(eigr), dtype=eigr.dtype)
    if 0:  # pragma: no cover
        denom = np.sqrt(eigr ** 2 + eigi ** 2)
        inonzero = np.where(denom != 0)[0]
        if len(inonzero):
            damping[inonzero] = -eigr[inonzero] / denom[inonzero]

        # not sure
        abs_freqs = np.sqrt(np.abs(eigi)) / (2 * np.pi)
    else:
        # flutter
        # eig = omega*gamma + omega*1j = eigr + eigi*1j
        # freq = eigi/(2*pi)
        # g = 2*eigr/eigi
        # g = eigr/ln(2) * L/V when eigi=0
        abs_freqs = abs(eigi) / (2 * np.pi)
        izero = np.where(eigi == 0)[0]
        inonzero = np.where(eigi != 0)[0]
        if len(izero):
            ln2 = np.log(2)
            damping[izero] = eigr[izero] / ln2 * cref / velocity
        if len(inonzero):
            damping[inonzero] = 2 * eigr[inonzero] / eigi[inonzero]
    return damping, abs_freqs
