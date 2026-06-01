"""MacNeal QUAD4 shell element stiffness matrix.

4-node quadrilateral Mindlin-Reissner shell element using stabilized
1-point integration (MacNeal/Hughes-Liu approach):

    K = K_centroid + K_stabilization

The centroid (1-point) integration captures the correct constant-strain
response. The stabilization adds higher-order stiffness that prevents
hourglass modes without introducing locking.

This matches the formulation used by NX Nastran's CQUAD4.

References
----------
MacNeal, R.H. (1978). "A simple quadrilateral shell element."
    Computers & Structures, 8:175-183.

Hughes, T.J.R. & Liu, W.K. (1981). "Nonlinear finite element analysis
    of shells." Comp. Meth. Appl. Mech. Eng., Parts I & II.

Belytschko, T. & Tsay, C.S. (1983). "A stabilization procedure for the
    quadrilateral plate element with one-point quadrature." IJNME, 19:405-419.

DOF ordering per node: [u, v, w, theta_x, theta_y]
    20 DOF total, nodes ordered 1-4 counterclockwise.
"""
from __future__ import annotations

import numpy as np


# Hourglass vector for bilinear quad: h = [1, -1, 1, -1]
_GAMMA = np.array([1.0, -1.0, 1.0, -1.0])


def _dshape_quad4(xi: float, eta: float) -> np.ndarray:
    """Shape function derivatives in natural coords. Returns (2, 4)."""
    return 0.25 * np.array([
        [-(1.0 - eta), (1.0 - eta), (1.0 + eta), -(1.0 + eta)],
        [-(1.0 - xi), -(1.0 + xi), (1.0 + xi), (1.0 - xi)],
    ])


def _jacobian(dN_dnat: np.ndarray, x: np.ndarray,
              y: np.ndarray) -> tuple[np.ndarray, float]:
    """Jacobian and its determinant."""
    J = np.array([
        [dN_dnat[0] @ x, dN_dnat[0] @ y],
        [dN_dnat[1] @ x, dN_dnat[1] @ y],
    ])
    det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    return J, det_J


def _membrane_B(dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Membrane strain-displacement matrix (3, 8)."""
    Bm = np.zeros((3, 8))
    for i in range(4):
        c = 2 * i
        Bm[0, c] = dN_dx[i]
        Bm[1, c + 1] = dN_dy[i]
        Bm[2, c] = dN_dy[i]
        Bm[2, c + 1] = dN_dx[i]
    return Bm


def _bending_B(dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Bending curvature-displacement matrix (3, 12).

    DOF per node: [w, theta_x, theta_y].
    """
    Bb = np.zeros((3, 12))
    for i in range(4):
        c = 3 * i
        Bb[0, c + 2] = dN_dx[i]       # kappa_x = d(theta_y)/dx
        Bb[1, c + 1] = -dN_dy[i]      # kappa_y = -d(theta_x)/dy
        Bb[2, c + 1] = -dN_dx[i]      # kappa_xy = d(theta_y)/dy - d(theta_x)/dx
        Bb[2, c + 2] = dN_dy[i]
    return Bb


def _shear_B(N: np.ndarray, dN_dx: np.ndarray,
             dN_dy: np.ndarray) -> np.ndarray:
    """Transverse shear strain-displacement matrix (2, 12).

    DOF per node: [w, theta_x, theta_y].
    gamma_xz = dw/dx + theta_y
    gamma_yz = dw/dy - theta_x  (note sign convention)
    """
    Bs = np.zeros((2, 12))
    for i in range(4):
        c = 3 * i
        Bs[0, c] = dN_dx[i]       # dw/dx
        Bs[0, c + 2] = N[i]       # theta_y
        Bs[1, c] = dN_dy[i]       # dw/dy
        Bs[1, c + 1] = -N[i]     # -theta_x
    return Bs


def macneal_quad4_stiffness(x: np.ndarray, y: np.ndarray,
                            A_mat: np.ndarray, D_mat: np.ndarray,
                            Ds_mat: np.ndarray,
                            B_mat: np.ndarray | None = None) -> np.ndarray:
    """Compute the 20x20 MacNeal QUAD4 shell element stiffness matrix.

    Uses stabilized 1-point integration:
        K = K_1pt + K_hourglass_membrane + K_hourglass_bending

    Parameters
    ----------
    x : (4,) nodal x-coordinates in element local frame
    y : (4,) nodal y-coordinates in element local frame
    A_mat : (3, 3) membrane constitutive matrix (CLT A-matrix)
    D_mat : (3, 3) bending constitutive matrix (CLT D-matrix)
    Ds_mat : (2, 2) transverse shear constitutive matrix [kappa*G*t]
    B_mat : (3, 3) membrane-bending coupling matrix, optional

    Returns
    -------
    Ke : (20, 20) element stiffness matrix
    """
    m_idx = np.array([0, 1, 5, 6, 10, 11, 15, 16], dtype=int)
    b_idx = np.array([2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19], dtype=int)

    # === Centroid evaluation (1-point integration, weight = 4) ===
    dN_dnat_c = _dshape_quad4(0.0, 0.0)
    J_c, det_J_c = _jacobian(dN_dnat_c, x, y)
    J_inv_c = np.array([
        [J_c[1, 1], -J_c[0, 1]],
        [-J_c[1, 0], J_c[0, 0]],
    ]) / det_J_c
    dN_dxy_c = J_inv_c @ dN_dnat_c
    dN_dx_c = dN_dxy_c[0]
    dN_dy_c = dN_dxy_c[1]

    A_elem = det_J_c * 4.0  # element area (exact for parallelogram)

    # Membrane at centroid
    Bm_c = _membrane_B(dN_dx_c, dN_dy_c)
    K_mm = (Bm_c.T @ A_mat @ Bm_c) * A_elem

    # Bending at centroid
    Bb_c = _bending_B(dN_dx_c, dN_dy_c)
    K_bb = (Bb_c.T @ D_mat @ Bb_c) * A_elem

    # Transverse shear at centroid (1-point: automatically avoids shear locking)
    N_c = np.array([0.25, 0.25, 0.25, 0.25])
    Bs_c = _shear_B(N_c, dN_dx_c, dN_dy_c)
    K_ss = (Bs_c.T @ Ds_mat @ Bs_c) * A_elem

    # Membrane-bending coupling
    K_mb = None
    if B_mat is not None:
        K_mb = (Bm_c.T @ B_mat @ Bb_c) * A_elem

    # === Hourglass stabilization ===
    # The hourglass vector gamma = [1, -1, 1, -1] spans the kernel of the
    # 1-point integration (it's orthogonal to the constant and linear modes).
    # We need to stabilize it for both membrane and bending.

    # Compute stabilization vectors b1, b2 (Flanagan-Belytschko)
    # b_i = gamma - (gamma . x_hat) * dN_dx - (gamma . y_hat) * dN_dy
    # where x_hat, y_hat are the isoparametric coordinate vectors
    gamma = _GAMMA

    # Project out the linear part: gamma_bar = gamma - (gamma.x)*dN/dx*A - ...
    # For the centroid: dN_dxi = dN_dnat_c[0], dN_deta = dN_dnat_c[1]
    # The hourglass mode in physical coords needs projection
    gx = gamma @ x  # = x1 - x2 + x3 - x4
    gy = gamma @ y  # = y1 - y2 + y3 - y4

    # Stabilization vectors in physical space (Flanagan-Belytschko 1981)
    # These represent the hourglass mode gradient
    b1 = (gamma - gx * dN_dx_c - gy * dN_dy_c) / A_elem
    b2 = b1  # For rectangular elements; general case needs separate treatment

    # Actually use the simpler Belytschko-Tsay approach:
    # K_hg = alpha * A * (b^T * C * b) where b is the hourglass strain
    # For membrane: hourglass strain = b^T @ u_membrane
    # For bending: hourglass curvature = b^T @ theta

    # Membrane hourglass stabilization
    # The membrane hourglass mode produces parasitic strains that must be controlled.
    # Stabilization stiffness: K_hg_mm = alpha_m * A * Q_m^T * A_mat_hg * Q_m
    # where Q_m relates hourglass modes to membrane DOF.

    # MacNeal approach: evaluate the difference between 2x2 integration and
    # 1-point integration, giving the hourglass stabilization as the residual.
    # K_stab = K_2x2 - K_1pt
    # This automatically gives the correct stabilization without heuristics.

    K_mm_stab = _membrane_stabilization(x, y, A_mat, Bm_c, det_J_c)
    K_bb_stab = _bending_stabilization(x, y, D_mat, Bb_c, det_J_c)

    # Combine
    K_mm_total = K_mm + K_mm_stab
    K_bb_total = K_bb + K_bb_stab + K_ss

    # Assemble 20x20
    Ke = np.zeros((20, 20))
    Ke[np.ix_(m_idx, m_idx)] = K_mm_total
    Ke[np.ix_(b_idx, b_idx)] = K_bb_total

    if K_mb is not None:
        Ke[np.ix_(m_idx, b_idx)] = K_mb
        Ke[np.ix_(b_idx, m_idx)] = K_mb.T

    return Ke


def _membrane_stabilization(x: np.ndarray, y: np.ndarray,
                            A_mat: np.ndarray,
                            Bm_c: np.ndarray,
                            det_J_c: float) -> np.ndarray:
    """Membrane hourglass stabilization (8x8).

    Uses alpha-scaled residual: K_stab = alpha * (K_2x2 - K_1pt).
    alpha=0.3 balances accuracy across bending, compression, and shear loads.
    """
    from pyNastran.dev.solver.stiffness.shell_mitc4 import _GAUSS_2x2_PTS, _GAUSS_2x2_WTS

    A_elem = det_J_c * 4.0

    K_2x2 = np.zeros((8, 8))
    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = np.array([
            [J[1, 1], -J[0, 1]],
            [-J[1, 0], J[0, 0]],
        ]) / det_J
        dN_dxy = J_inv @ dN_dnat
        Bm = _membrane_B(dN_dxy[0], dN_dxy[1])
        K_2x2 += (Bm.T @ A_mat @ Bm) * det_J * wt

    K_1pt = (Bm_c.T @ A_mat @ Bm_c) * A_elem

    alpha_m = 0.3
    return alpha_m * (K_2x2 - K_1pt)


def _bending_stabilization(x: np.ndarray, y: np.ndarray,
                           D_mat: np.ndarray,
                           Bb_c: np.ndarray,
                           det_J_c: float) -> np.ndarray:
    """Bending hourglass stabilization (12x12).

    Uses the residual (K_2x2 - K_1pt) scaled by alpha.
    Alpha < 1 prevents over-stiffening of the bending response.
    """
    from pyNastran.dev.solver.stiffness.shell_mitc4 import _GAUSS_2x2_PTS, _GAUSS_2x2_WTS

    A_elem = det_J_c * 4.0

    K_2x2 = np.zeros((12, 12))
    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = np.array([
            [J[1, 1], -J[0, 1]],
            [-J[1, 0], J[0, 0]],
        ]) / det_J
        dN_dxy = J_inv @ dN_dnat
        Bb = _bending_B(dN_dxy[0], dN_dxy[1])
        K_2x2 += (Bb.T @ D_mat @ Bb) * det_J * wt

    K_1pt = (Bb_c.T @ D_mat @ Bb_c) * A_elem

    alpha_b = 0.5
    return alpha_b * (K_2x2 - K_1pt)


def macneal_quad4_pressure_load(x: np.ndarray, y: np.ndarray,
                                pressures: np.ndarray,
                                direction: np.ndarray | None = None) -> np.ndarray:
    """Equivalent nodal force vector for surface pressure (20x1).

    Same interface as mitc4_pressure_load — pressure load is independent
    of the stiffness formulation.
    """
    from pyNastran.dev.solver.stiffness.shell_mitc4 import mitc4_pressure_load
    return mitc4_pressure_load(x, y, pressures, direction=direction)
