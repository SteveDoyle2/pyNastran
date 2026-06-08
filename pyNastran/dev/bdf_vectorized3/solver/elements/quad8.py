"""CQUAD8 element stiffness matrix (8-node serendipity shell).

8-node quadrilateral Mindlin-Reissner shell element with:
- Serendipity (quadratic) shape functions
- 3x3 Gauss quadrature for membrane and bending
- 2x2 Gauss quadrature for transverse shear (reduced integration)

DOF ordering per node: [u, v, w, theta_x, theta_y] (5 DOF)
Total: 40 DOF per element.

Notes
-----
This is a prototype implementation. The element passes symmetry and
positive semi-definiteness checks. The 2x2 reduced integration for
transverse shear introduces 4 spurious zero-energy modes beyond the
6 rigid body modes (rank = 30 instead of 34). A stabilization or
assumed-strain shear treatment (analogous to MITC8) is needed to
eliminate these. Not yet validated against NX Nastran KELM.

References
----------
Zienkiewicz, O.C. & Taylor, R.L. (2000). "The Finite Element Method",
Vol. 2, 5th edition. Butterworth-Heinemann, Section 4.10.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# 8-node serendipity shape functions
# ---------------------------------------------------------------------------
def _shape_quad8(xi: float, eta: float) -> np.ndarray:
    """Serendipity shape functions for 8-node quad.

    Node ordering:
        4---7---3
        |       |
        8       6
        |       |
        1---5---2

    Returns (8,) array.
    """
    # Corner nodes (1-4)
    N = np.zeros(8)
    N[0] = 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1)
    N[1] = 0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1)
    N[2] = 0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1)
    N[3] = 0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1)
    # Midside nodes (5-8)
    N[4] = 0.5 * (1 - xi**2) * (1 - eta)
    N[5] = 0.5 * (1 + xi) * (1 - eta**2)
    N[6] = 0.5 * (1 - xi**2) * (1 + eta)
    N[7] = 0.5 * (1 - xi) * (1 - eta**2)
    return N


def _dshape_quad8(xi: float, eta: float) -> np.ndarray:
    """Shape function derivatives for 8-node quad.

    Returns (2, 8) array: [dN/dxi; dN/deta].
    """
    dN = np.zeros((2, 8))

    # dN/dxi
    dN[0, 0] = 0.25 * (1 - eta) * (2 * xi + eta)
    dN[0, 1] = 0.25 * (1 - eta) * (2 * xi - eta)
    dN[0, 2] = 0.25 * (1 + eta) * (2 * xi + eta)
    dN[0, 3] = 0.25 * (1 + eta) * (2 * xi - eta)
    dN[0, 4] = -xi * (1 - eta)
    dN[0, 5] = 0.5 * (1 - eta**2)
    dN[0, 6] = -xi * (1 + eta)
    dN[0, 7] = -0.5 * (1 - eta**2)

    # dN/deta
    dN[1, 0] = 0.25 * (1 - xi) * (xi + 2 * eta)
    dN[1, 1] = 0.25 * (1 + xi) * (-xi + 2 * eta)
    dN[1, 2] = 0.25 * (1 + xi) * (xi + 2 * eta)
    dN[1, 3] = 0.25 * (1 - xi) * (-xi + 2 * eta)
    dN[1, 4] = -0.5 * (1 - xi**2)
    dN[1, 5] = -(1 + xi) * eta
    dN[1, 6] = 0.5 * (1 - xi**2)
    dN[1, 7] = -(1 - xi) * eta

    return dN


def _jacobian_q8(dN: np.ndarray, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    """Jacobian (2x2) and its determinant for 8-node quad."""
    J = np.array(
        [
            [dN[0] @ x, dN[0] @ y],
            [dN[1] @ x, dN[1] @ y],
        ]
    )
    det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    return J, det_J


# ---------------------------------------------------------------------------
# Strain-displacement matrices
# ---------------------------------------------------------------------------
def _membrane_B_q8(dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Membrane B-matrix (3, 16) for 8-node quad.

    DOF order: [u1, v1, u2, v2, ..., u8, v8].
    """
    Bm = np.zeros((3, 16))
    Bm[0, 0::2] = dN_dx
    Bm[1, 1::2] = dN_dy
    Bm[2, 0::2] = dN_dy
    Bm[2, 1::2] = dN_dx
    return Bm


def _bending_B_q8(dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Bending curvature-displacement matrix (3, 24) for 8-node quad.

    DOF order per node: [w, theta_x, theta_y].
    """
    Bb = np.zeros((3, 24))
    Bb[0, 2::3] = dN_dx  # kappa_x = d(theta_y)/dx
    Bb[1, 1::3] = -dN_dy  # kappa_y = -d(theta_x)/dy
    Bb[2, 1::3] = -dN_dx  # kappa_xy
    Bb[2, 2::3] = dN_dy
    return Bb


def _shear_B_q8(N: np.ndarray, dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Transverse shear B-matrix (2, 24) for 8-node quad.

    DOF order per node: [w, theta_x, theta_y].
    gamma_xz = dw/dx + theta_y
    gamma_yz = dw/dy - theta_x
    """
    Bs = np.zeros((2, 24))
    Bs[0, 0::3] = dN_dx
    Bs[0, 2::3] = N
    Bs[1, 0::3] = dN_dy
    Bs[1, 1::3] = -N
    return Bs


# ---------------------------------------------------------------------------
# Gauss quadrature points
# ---------------------------------------------------------------------------
def _gauss_3x3() -> tuple[np.ndarray, np.ndarray]:
    """3x3 Gauss quadrature points and weights."""
    g = np.sqrt(3.0 / 5.0)
    pts_1d = np.array([-g, 0.0, g])
    wts_1d = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])

    pts = np.array([[xi, eta] for eta in pts_1d for xi in pts_1d])
    wts = np.array([wx * wy for wy in wts_1d for wx in wts_1d])
    return pts, wts


def _gauss_2x2() -> tuple[np.ndarray, np.ndarray]:
    """2x2 Gauss quadrature points and weights."""
    g = 1.0 / np.sqrt(3.0)
    pts = np.array([[-g, -g], [g, -g], [g, g], [-g, g]])
    wts = np.array([1.0, 1.0, 1.0, 1.0])
    return pts, wts


# ---------------------------------------------------------------------------
# Element stiffness
# ---------------------------------------------------------------------------
def quad8_stiffness(
    x: np.ndarray,
    y: np.ndarray,
    A_mat: np.ndarray,
    D_mat: np.ndarray,
    Ds_mat: np.ndarray,
    B_mat: np.ndarray | None = None,
) -> np.ndarray:
    """Compute the 40x40 CQUAD8 shell element stiffness matrix.

    Parameters
    ----------
    x : (8,) nodal x-coordinates in element local frame
    y : (8,) nodal y-coordinates in element local frame
    A_mat : (3, 3) membrane constitutive
    D_mat : (3, 3) bending constitutive
    Ds_mat : (2, 2) transverse shear constitutive
    B_mat : (3, 3) membrane-bending coupling, optional

    Returns
    -------
    Ke : (40, 40) element stiffness matrix
        DOF ordering per node: [u, v, w, theta_x, theta_y]
        Nodes ordered 1-8.
    """
    # Index maps into 40-DOF vector
    m_idx = np.zeros(16, dtype=int)
    b_idx = np.zeros(24, dtype=int)
    for i in range(8):
        m_idx[2 * i] = 5 * i
        m_idx[2 * i + 1] = 5 * i + 1
        b_idx[3 * i] = 5 * i + 2
        b_idx[3 * i + 1] = 5 * i + 3
        b_idx[3 * i + 2] = 5 * i + 4

    K_mm = np.zeros((16, 16))
    K_bb = np.zeros((24, 24))
    K_ss = np.zeros((24, 24))
    K_mb = np.zeros((16, 24))

    # Membrane and bending: 3x3 Gauss
    pts_3x3, wts_3x3 = _gauss_3x3()
    for gpt, wt in zip(pts_3x3, wts_3x3):
        xi, eta = gpt
        dN_dnat = _dshape_quad8(xi, eta)
        J, det_J = _jacobian_q8(dN_dnat, x, y)
        J_inv = (
            np.array(
                [
                    [J[1, 1], -J[0, 1]],
                    [-J[1, 0], J[0, 0]],
                ]
            )
            / det_J
        )
        dN_dxy = J_inv @ dN_dnat
        dN_dx = dN_dxy[0]
        dN_dy = dN_dxy[1]

        Bm = _membrane_B_q8(dN_dx, dN_dy)
        K_mm += (Bm.T @ A_mat @ Bm) * det_J * wt

        Bb = _bending_B_q8(dN_dx, dN_dy)
        K_bb += (Bb.T @ D_mat @ Bb) * det_J * wt

        if B_mat is not None:
            K_mb += (Bm.T @ B_mat @ Bb) * det_J * wt

    # Transverse shear: 2x2 reduced integration (avoids shear locking)
    pts_2x2, wts_2x2 = _gauss_2x2()
    for gpt, wt in zip(pts_2x2, wts_2x2):
        xi, eta = gpt
        N = _shape_quad8(xi, eta)
        dN_dnat = _dshape_quad8(xi, eta)
        J, det_J = _jacobian_q8(dN_dnat, x, y)
        J_inv = (
            np.array(
                [
                    [J[1, 1], -J[0, 1]],
                    [-J[1, 0], J[0, 0]],
                ]
            )
            / det_J
        )
        dN_dxy = J_inv @ dN_dnat
        dN_dx = dN_dxy[0]
        dN_dy = dN_dxy[1]

        Bs = _shear_B_q8(N, dN_dx, dN_dy)
        K_ss += (Bs.T @ Ds_mat @ Bs) * det_J * wt

    # Assemble into 40x40
    Ke = np.zeros((40, 40))
    Ke[np.ix_(m_idx, m_idx)] = K_mm
    Ke[np.ix_(b_idx, b_idx)] = K_bb + K_ss
    if B_mat is not None:
        Ke[np.ix_(m_idx, b_idx)] = K_mb
        Ke[np.ix_(b_idx, m_idx)] = K_mb.T

    return Ke


def quad8_mass(
    x: np.ndarray,
    y: np.ndarray,
    mass_per_area: float,
) -> np.ndarray:
    """Compute the 40x40 CQUAD8 consistent mass matrix (translational only).

    Parameters
    ----------
    x : (8,) nodal x-coordinates in element local frame
    y : (8,) nodal y-coordinates in element local frame
    mass_per_area : float
        rho * thickness

    Returns
    -------
    Me : (40, 40) consistent mass matrix
    """
    # Only w DOF (translation normal to shell), indices 2, 7, 12, ...
    # For simplicity, apply to all translational DOFs (u, v, w)
    Me = np.zeros((40, 40))

    pts_3x3, wts_3x3 = _gauss_3x3()
    for gpt, wt in zip(pts_3x3, wts_3x3):
        xi, eta = gpt
        N = _shape_quad8(xi, eta)
        dN_dnat = _dshape_quad8(xi, eta)
        _, det_J = _jacobian_q8(dN_dnat, x, y)

        # Consistent mass: Me_ij = rho*t * integral(N_i * N_j) dA
        # For translational DOFs u, v, w:
        NNT = np.outer(N, N) * mass_per_area * det_J * wt
        for a in range(8):
            for b in range(8):
                val = NNT[a, b]
                if abs(val) > 0.0:
                    # u DOF
                    Me[5 * a, 5 * b] += val
                    # v DOF
                    Me[5 * a + 1, 5 * b + 1] += val
                    # w DOF
                    Me[5 * a + 2, 5 * b + 2] += val

    return Me
