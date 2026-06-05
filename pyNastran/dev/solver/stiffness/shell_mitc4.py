"""MITC4 shell element stiffness matrix.

4-node quadrilateral Mindlin-Reissner shell element with:
- Bilinear membrane (2x2 Gauss)
- Mindlin-Reissner bending (2x2 Gauss)
- MITC4 assumed natural strain for transverse shear (locking-free)
- Optional membrane-bending coupling (B-matrix from CLT)

References
----------
Dvorkin, E.N. & Bathe, K.J. (1984). "A continuum mechanics based four-node
shell element for general nonlinear analysis." Eng. Computations, 1:77-88.

Bathe, K.J. & Dvorkin, E.N. (1986). "A formulation of general shell elements --
the use of mixed interpolation of tensorial components." IJNME, 22:697-722.
"""
from __future__ import annotations

from itertools import count
from typing import TYPE_CHECKING

import numpy as np

from pyNastran.bdf.cards.elements.shell import transform_shell_material_coordinate_system
from ..utils import DOF_MAP

if TYPE_CHECKING:
    from pyNastran.nptyping_interface import NDArrayN3float
    from pyNastran.bdf.bdf import BDF, CQUAD4

# 2x2 Gauss quadrature points and weights
_GAUSS_2x2_PTS = np.array([
    [-1.0 / np.sqrt(3.0), -1.0 / np.sqrt(3.0)],
    [+1.0 / np.sqrt(3.0), -1.0 / np.sqrt(3.0)],
    [+1.0 / np.sqrt(3.0), +1.0 / np.sqrt(3.0)],
    [-1.0 / np.sqrt(3.0), +1.0 / np.sqrt(3.0)],
])
_GAUSS_2x2_WTS = np.array([1.0, 1.0, 1.0, 1.0])

# MITC4 tying points for transverse shear
# e_xi3: sampled at midpoints of eta=const edges (bottom and top)
_TYING_XI = np.array([
    [0.0, -1.0],  # A: midpoint of edge 1-2 (bottom)
    [0.0, +1.0],  # B: midpoint of edge 4-3 (top)
])
# e_eta3: sampled at midpoints of xi=const edges (left and right)
_TYING_ETA = np.array([
    [-1.0, 0.0],  # D: midpoint of edge 1-4 (left)
    [+1.0, 0.0],  # C: midpoint of edge 2-3 (right)
])


def _shape_quad4(xi: float, eta: float) -> np.ndarray:
    """Bilinear shape functions for 4-node quad.

    Returns (4,) array [N1, N2, N3, N4].
    """
    return 0.25 * np.array([
        (1.0 - xi) * (1.0 - eta),
        (1.0 + xi) * (1.0 - eta),
        (1.0 + xi) * (1.0 + eta),
        (1.0 - xi) * (1.0 + eta),
    ])


def _dshape_quad4(xi: float, eta: float) -> np.ndarray:
    """Derivatives of shape functions w.r.t. natural coordinates.

    Returns (2, 4) array:
        row 0: dN/dxi  = [dN1/dxi, dN2/dxi, dN3/dxi, dN4/dxi]
        row 1: dN/deta = [dN1/deta, dN2/deta, dN3/deta, dN4/deta]
    """
    return 0.25 * np.array([
        [-(1.0 - eta), (1.0 - eta), (1.0 + eta), -(1.0 + eta)],
        [-(1.0 - xi), -(1.0 + xi), (1.0 + xi), (1.0 - xi)],
    ])


def _jacobian(dN_dnat: np.ndarray, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    """Compute Jacobian matrix and its determinant.

    Parameters
    ----------
    dN_dnat : (2, 4) shape function derivatives in natural coords
    x, y : (4,) nodal coordinates in element local frame

    Returns
    -------
    J : (2, 2) Jacobian matrix [[dx/dxi, dy/dxi], [dx/deta, dy/deta]]
    det_J : determinant of J
    """
    J = np.array([
        [dN_dnat[0] @ x, dN_dnat[0] @ y],
        [dN_dnat[1] @ x, dN_dnat[1] @ y],
    ])
    det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    return J, det_J


def _split_membrane_constitutive(
    A: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Split membrane constitutive into volumetric and deviatoric parts.

    For selective-reduced integration: volumetric (dilatational) is integrated
    at 1 point, deviatoric at 2x2.

    The volumetric projector for plane stress:
        P_vol = [[1,1,0],[1,1,0],[0,0,0]] / 2

    A_vol = P_vol^T @ A @ P_vol (projected volumetric stiffness)
    A_dev = A - A_vol

    Returns (A_vol, A_dev).
    """
    # Volumetric projection: eps_vol = (eps_xx + eps_yy)/2 * [1, 1, 0]^T
    # The volumetric part of A relates tr(eps) to tr(sigma):
    # A_vol_ij = sum_k sum_l P_ki * A_kl * P_lj where P = [[1,1,0],[1,1,0],[0,0,0]]/2
    P = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 0]], dtype=float) / 2.0
    A_vol = P.T @ A @ P
    A_dev = A - A_vol
    return A_vol, A_dev


def _membrane_B(dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Membrane strain-displacement matrix.

    Relates {eps_x, eps_y, gamma_xy} to {u1, v1, u2, v2, u3, v3, u4, v4}.
    Returns (3, 8) matrix.
    """
    Bm = np.zeros((3, 8))
    for i in range(4):
        c = 2 * i
        Bm[0, c] = dN_dx[i]       # eps_x = du/dx
        Bm[1, c + 1] = dN_dy[i]   # eps_y = dv/dy
        Bm[2, c] = dN_dy[i]       # gamma_xy = du/dy + dv/dx
        Bm[2, c + 1] = dN_dx[i]
    return Bm


def _bending_B(dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Bending curvature-displacement matrix.

    Uses sign convention:
        u = z * theta_y,  v = -z * theta_x
        kappa_x  =  d(theta_y)/dx
        kappa_y  = -d(theta_x)/dy
        kappa_xy =  d(theta_y)/dy - d(theta_x)/dx

    DOF per node: [w, theta_x, theta_y]
    Relates {kappa_x, kappa_y, kappa_xy} to {w1,tx1,ty1, ..., w4,tx4,ty4}.
    Returns (3, 12) matrix.
    """
    Bb = np.zeros((3, 12))
    for i in range(4):
        c = 3 * i
        # kappa_x = d(theta_y)/dx
        Bb[0, c + 2] = dN_dx[i]
        # kappa_y = -d(theta_x)/dy
        Bb[1, c + 1] = -dN_dy[i]
        # kappa_xy = d(theta_y)/dy - d(theta_x)/dx
        Bb[2, c + 1] = -dN_dx[i]
        Bb[2, c + 2] = dN_dy[i]
    return Bb


def _shear_B_standard(N: np.ndarray, dN_dx: np.ndarray, dN_dy: np.ndarray) -> np.ndarray:
    """Standard transverse shear strain-displacement matrix (displacement-based).

    gamma_xz = dw/dx + theta_y
    gamma_yz = dw/dy - theta_x

    DOF per node: [w, theta_x, theta_y]
    Relates {gamma_xz, gamma_yz} to {w1,tx1,ty1, ..., w4,tx4,ty4}.
    Returns (2, 12) matrix.
    """
    Bs = np.zeros((2, 12))
    for i in range(4):
        c = 3 * i
        Bs[0, c] = dN_dx[i]     # gamma_xz: dw/dx
        Bs[0, c + 2] = N[i]     # gamma_xz: +theta_y
        Bs[1, c] = dN_dy[i]     # gamma_yz: dw/dy
        Bs[1, c + 1] = -N[i]    # gamma_yz: -theta_x
    return Bs


def _covariant_shear_row(N: np.ndarray, dN_dnat_row: np.ndarray,
                         J_row: np.ndarray) -> np.ndarray:
    """Compute one covariant transverse shear strain row.

    e_alpha3 = dw/d(alpha) + J_alpha1 * theta_y - J_alpha2 * theta_x

    where alpha is either xi or eta, and:
        dw/d(alpha) = sum_i dN_i/d(alpha) * w_i
        J_alpha1 = dx/d(alpha), J_alpha2 = dy/d(alpha)

    DOF per node: [w, theta_x, theta_y]
    Returns (12,) vector.
    """
    J_a1 = J_row[0]  # dx/d(alpha)
    J_a2 = J_row[1]  # dy/d(alpha)
    row = np.zeros(12)
    for i in range(4):
        c = 3 * i
        row[c] = dN_dnat_row[i]     # dw/d(alpha)
        row[c + 1] = -J_a2 * N[i]   # -dy/d(alpha) * theta_x
        row[c + 2] = J_a1 * N[i]    # +dx/d(alpha) * theta_y
    return row


def _mitc4_shear_B(xi: float, eta: float,
                   x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """MITC4 assumed-strain transverse shear B-matrix at point (xi, eta).

    Evaluates the covariant shear strains at tying points, interpolates,
    then transforms to Cartesian components.

    Returns (2, 12) matrix relating {gamma_xz, gamma_yz} to bending DOF.
    """
    # --- Evaluate e_xi3 at tying points A=(0,-1) and B=(0,+1) ---
    e_xi3_rows = np.zeros((2, 12))
    for itp, (xi_tp, eta_tp) in enumerate(_TYING_XI):
        N_tp = _shape_quad4(xi_tp, eta_tp)
        dN_tp = _dshape_quad4(xi_tp, eta_tp)
        J_tp, _ = _jacobian(dN_tp, x, y)
        # e_xi3 uses the xi-row of J (row 0)
        e_xi3_rows[itp] = _covariant_shear_row(N_tp, dN_tp[0], J_tp[0])

    # --- Evaluate e_eta3 at tying points D=(-1,0) and C=(+1,0) ---
    e_eta3_rows = np.zeros((2, 12))
    for itp, (xi_tp, eta_tp) in enumerate(_TYING_ETA):
        N_tp = _shape_quad4(xi_tp, eta_tp)
        dN_tp = _dshape_quad4(xi_tp, eta_tp)
        J_tp, _ = _jacobian(dN_tp, x, y)
        # e_eta3 uses the eta-row of J (row 1)
        e_eta3_rows[itp] = _covariant_shear_row(N_tp, dN_tp[1], J_tp[1])

    # --- Interpolate to current point (xi, eta) ---
    # e_xi3^MITC  = 0.5*(1-eta)*e_xi3(A) + 0.5*(1+eta)*e_xi3(B)
    e_xi3_interp = 0.5 * (1.0 - eta) * e_xi3_rows[0] + 0.5 * (1.0 + eta) * e_xi3_rows[1]
    # e_eta3^MITC = 0.5*(1-xi)*e_eta3(D) + 0.5*(1+xi)*e_eta3(C)
    e_eta3_interp = 0.5 * (1.0 - xi) * e_eta3_rows[0] + 0.5 * (1.0 + xi) * e_eta3_rows[1]

    # --- Transform covariant to Cartesian: gamma = J^{-1} * e_nat ---
    dN_here = _dshape_quad4(xi, eta)
    J_here, det_J = _jacobian(dN_here, x, y)
    J_inv = np.array([
        [J_here[1, 1], -J_here[0, 1]],
        [-J_here[1, 0], J_here[0, 0]],
    ]) / det_J

    # Bs_mitc[0,:] = J_inv[0,0]*e_xi3 + J_inv[0,1]*e_eta3  (gamma_xz)
    # Bs_mitc[1,:] = J_inv[1,0]*e_xi3 + J_inv[1,1]*e_eta3  (gamma_yz)
    Bs_mitc = np.zeros((2, 12))
    Bs_mitc[0] = J_inv[0, 0] * e_xi3_interp + J_inv[0, 1] * e_eta3_interp
    Bs_mitc[1] = J_inv[1, 0] * e_xi3_interp + J_inv[1, 1] * e_eta3_interp
    return Bs_mitc


def _eas_membrane_4mode(x: np.ndarray, y: np.ndarray,
                        A_mat: np.ndarray) -> np.ndarray:
    """Compute EAS-4 membrane correction via static condensation.

    Adds 4 incompatible strain modes (xi, eta variations in eps_xx, eps_yy)
    that are statically condensed to produce a correction to K_mm.

    Returns K_eas_correction (8, 8) to be SUBTRACTED from K_mm.
    """
    # Centroid Jacobian (reference for EAS mapping)
    dN_c = _dshape_quad4(0.0, 0.0)
    J0, det_J0 = _jacobian(dN_c, x, y)
    J0_invT = np.linalg.inv(J0).T

    # EAS modes: 4 internal DOFs (alpha)
    # Enhanced strain interpolation in natural coords:
    #   M_nat = [[xi, 0], [0, eta], [0, 0], [0, 0]]  (3x4 in Voigt: exx, eyy, gxy)
    # Mapped to physical: M_phys = det(J0)/det(J) * T @ M_nat
    # where T transforms strain from natural to physical coords

    K_aa = np.zeros((4, 4))   # internal-internal
    K_au = np.zeros((4, 8))   # internal-displacement coupling

    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt
        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)
        J_inv = np.array([
            [J[1, 1], -J[0, 1]],
            [-J[1, 0], J[0, 0]],
        ]) / det_J
        dN_dxy = J_inv @ dN_dnat

        # Standard B-matrix
        Bm = _membrane_B(dN_dxy[0], dN_dxy[1])

        # EAS strain interpolation matrix (3x4) at this Gauss point
        # Natural EAS modes: [xi, 0, eta, 0; 0, xi, 0, eta; 0, 0, 0, 0]
        # maps 4 alpha DOFs to 3 strain components
        # Transform: M_phys = (det_J0/det_J) * T_strain @ M_nat
        # T_strain transforms contravariant strains: uses J0_invT tensor product
        # For Voigt notation: T_e[3x3] from J0_invT
        j = J0_invT
        T_e = np.array([
            [j[0, 0]**2, j[0, 1]**2, 2 * j[0, 0] * j[0, 1]],
            [j[1, 0]**2, j[1, 1]**2, 2 * j[1, 0] * j[1, 1]],
            [j[0, 0] * j[1, 0], j[0, 1] * j[1, 1],
             j[0, 0] * j[1, 1] + j[0, 1] * j[1, 0]],
        ])

        # Natural EAS modes (3x4): exx~xi, exx~eta, eyy~xi, eyy~eta
        M_nat = np.array([
            [xi, eta, 0, 0],
            [0, 0, xi, eta],
            [0, 0, 0, 0],
        ])

        M_phys = (det_J0 / det_J) * (T_e @ M_nat)  # (3, 4)

        # Accumulate
        K_aa += (M_phys.T @ A_mat @ M_phys) * det_J * wt
        K_au += (M_phys.T @ A_mat @ Bm) * det_J * wt

    # Static condensation: K_correction = K_ua^T @ K_aa^{-1} @ K_au
    K_aa_inv = np.linalg.inv(K_aa)
    K_correction = K_au.T @ K_aa_inv @ K_au

    return K_correction


def mitc4_stiffness(x: np.ndarray, y: np.ndarray,
                    A_mat: np.ndarray, D_mat: np.ndarray,
                    Ds_mat: np.ndarray,
                    B_mat: np.ndarray | None = None,
                    membrane: str = 'full',
                    eas: int = 0) -> np.ndarray:
    """Compute the 20x20 MITC4 shell element stiffness matrix.

    Parameters
    ----------
    x : (4,) nodal x-coordinates in element local frame
    y : (4,) nodal y-coordinates in element local frame
    A_mat : (3, 3) membrane constitutive matrix (CLT A-matrix)
    D_mat : (3, 3) bending constitutive matrix (CLT D-matrix)
    Ds_mat : (2, 2) transverse shear constitutive matrix [kappa*G*t]
    B_mat : (3, 3) membrane-bending coupling matrix (CLT B-matrix), optional
    membrane : str
        Membrane integration scheme:
        - 'full': 2x2 Gauss (default, original MITC4)
        - 'selective': selective-reduced integration (1-pt volumetric,
          2x2 deviatoric). Matches Nastran CQUAD4 membrane behavior.
    eas : int
        Number of EAS (Enhanced Assumed Strain) modes for membrane:
        - 0: disabled (default)
        - 4: incompatible modes for eps_xx and eps_yy (xi, eta variations).
          Eliminates parasitic membrane stiffness on coarse meshes.

    Returns
    -------
    Ke : (20, 20) element stiffness matrix

    DOF ordering per node: [u, v, w, theta_x, theta_y]
    Nodes ordered 1-4 counterclockwise.

    Notes
    -----
    Bending uses 2x2 Gauss quadrature.
    Transverse shear uses MITC4 assumed strain field (tying at mid-edge points).
    Membrane uses the specified integration scheme.

    With membrane='selective', the volumetric (pressure) part of the membrane
    constitutive is integrated at 1 point (centroid) while the deviatoric (shear)
    part uses 2x2. This prevents membrane locking on distorted elements and
    matches Nastran's selective-RI approach.
    """
    # Index maps into 20-DOF vector
    m_idx = np.array([0, 1, 5, 6, 10, 11, 15, 16], dtype=int)
    b_idx = np.array([2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19], dtype=int)

    K_mm = np.zeros((8, 8))
    K_bb = np.zeros((12, 12))
    K_ss = np.zeros((12, 12))
    K_mb = np.zeros((8, 12))

    # Selective-RI: split A into volumetric and deviatoric parts
    if membrane == 'selective':
        A_vol, A_dev = _split_membrane_constitutive(A_mat)
    elif membrane != 'full':
        raise ValueError(f"membrane must be 'full' or 'selective', got '{membrane}'")

    for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
        xi, eta = gpt

        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)

        J_inv = np.array([
            [J[1, 1], -J[0, 1]],
            [-J[1, 0], J[0, 0]],
        ]) / det_J
        dN_dxy = J_inv @ dN_dnat
        dN_dx = dN_dxy[0]
        dN_dy = dN_dxy[1]

        Bm = _membrane_B(dN_dx, dN_dy)

        if membrane == 'full':
            K_mm += (Bm.T @ A_mat @ Bm) * det_J * wt
        else:
            # Deviatoric part at 2x2 Gauss points
            K_mm += (Bm.T @ A_dev @ Bm) * det_J * wt

        Bb = _bending_B(dN_dx, dN_dy)
        K_bb += (Bb.T @ D_mat @ Bb) * det_J * wt

        Bs = _mitc4_shear_B(xi, eta, x, y)
        K_ss += (Bs.T @ Ds_mat @ Bs) * det_J * wt

        if B_mat is not None:
            K_mb += (Bm.T @ B_mat @ Bb) * det_J * wt

    # Selective-RI: volumetric part at centroid (1 point, weight=4)
    if membrane == 'selective':
        dN_dnat_c = _dshape_quad4(0.0, 0.0)
        J_c, det_J_c = _jacobian(dN_dnat_c, x, y)
        J_inv_c = np.array([
            [J_c[1, 1], -J_c[0, 1]],
            [-J_c[1, 0], J_c[0, 0]],
        ]) / det_J_c
        dN_dxy_c = J_inv_c @ dN_dnat_c
        Bm_c = _membrane_B(dN_dxy_c[0], dN_dxy_c[1])
        K_mm += (Bm_c.T @ A_vol @ Bm_c) * det_J_c * 4.0

    # EAS correction: softens membrane via statically condensed internal modes
    if eas == 4:
        K_mm -= _eas_membrane_4mode(x, y, A_mat)
    elif eas != 0:
        raise ValueError(f"eas must be 0 or 4, got {eas}")

    # Assemble into 20x20
    Ke = np.zeros((20, 20))
    Ke[np.ix_(m_idx, m_idx)] += K_mm
    Ke[np.ix_(b_idx, b_idx)] += K_bb + K_ss

    if B_mat is not None:
        Ke[np.ix_(m_idx, b_idx)] += K_mb
        Ke[np.ix_(b_idx, m_idx)] += K_mb.T

    return Ke


def mitc4_geometric_stiffness(x: np.ndarray, y: np.ndarray,
                              stress_resultants: np.ndarray,
                              moment_resultants: np.ndarray | None = None,
                              thickness: float = 0.0) -> np.ndarray:
    """Compute the 20x20 MITC4 geometric (differential) stiffness matrix.

    Parameters
    ----------
    x : (4,) nodal x-coordinates in element local frame
    y : (4,) nodal y-coordinates in element local frame
    stress_resultants : (3,) or (4, 3) membrane stress resultants [Nxx, Nyy, Nxy]
        If (3,): constant stress field across element (from element centroid)
        If (4, 3): stress at each Gauss point
    moment_resultants : (3,) or (4, 3) or None
        Bending moment resultants [Mxx, Myy, Mxy].
        When provided, amplifies the w-w geometric stiffness to account for
        through-thickness stress variation (offset shell P-delta effect).
    thickness : float
        Shell thickness. Used for through-thickness scaling.

    Returns
    -------
    Kg : (20, 20) geometric stiffness matrix

    DOF ordering per node: [u, v, w, theta_x, theta_y]
    Nodes ordered 1-4 counterclockwise.

    Notes
    -----
    For shells with offset (ZOFFS), the prebuckling bending moment M = N*z0
    amplifies the w-w geometric stiffness by factor (1 + 12*M^2/(N*t^2*N)).
    Equivalently, Kg_ww uses an effective stress resultant:
        N_eff = N + 12*M^2/(N*t^2)

    This captures the P-delta destabilization where eccentric compression
    creates bending that amplifies lateral instability.

    Reference: Bathe, K.J., "Finite Element Procedures" (2014), Ch. 6.
    """
    stress_resultants = np.asarray(stress_resultants)
    if stress_resultants.ndim == 1:
        stress_at_gp = np.tile(stress_resultants, (4, 1))
    else:
        stress_at_gp = stress_resultants

    has_moments = moment_resultants is not None
    if has_moments:
        moment_resultants = np.asarray(moment_resultants)
        if moment_resultants.ndim == 1:
            moment_at_gp = np.tile(moment_resultants, (4, 1))
        else:
            moment_at_gp = moment_resultants

    t2_12 = thickness**2 / 12.0 if thickness > 0.0 else 0.0

    Kg = np.zeros((20, 20))

    u_idx = np.array([0, 5, 10, 15], dtype=int)
    v_idx = np.array([1, 6, 11, 16], dtype=int)
    w_idx = np.array([2, 7, 12, 17], dtype=int)
    tx_idx = np.array([3, 8, 13, 18], dtype=int)
    ty_idx = np.array([4, 9, 14, 19], dtype=int)

    for igp, (gpt, wt) in enumerate(zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS)):
        xi, eta = gpt

        dN_dnat = _dshape_quad4(xi, eta)
        J, det_J = _jacobian(dN_dnat, x, y)

        J_inv = np.array([
            [J[1, 1], -J[0, 1]],
            [-J[1, 0], J[0, 0]],
        ]) / det_J
        dN_dxy = J_inv @ dN_dnat  # (2, 4)

        Nxx, Nyy, Nxy = stress_at_gp[igp]
        S_N = np.array([[Nxx, Nxy],
                        [Nxy, Nyy]])

        G = dN_dxy  # (2, 4)
        kg_N = (G.T @ S_N @ G) * det_J * wt  # (4, 4)

        # N terms: u-u, v-v, w-w
        Kg[np.ix_(u_idx, u_idx)] += kg_N
        Kg[np.ix_(v_idx, v_idx)] += kg_N
        Kg[np.ix_(w_idx, w_idx)] += kg_N

        # t^2/12 * N terms for theta-theta
        if t2_12 > 0.0:
            Kg[np.ix_(tx_idx, tx_idx)] += kg_N * t2_12
            Kg[np.ix_(ty_idx, ty_idx)] += kg_N * t2_12

        # Moment contribution: amplifies w-w by 12*M^2/(N*t^2)
        # Equivalent to Kg_ww *= (1 + 12*z0^2/t^2) for M = N*z0
        if has_moments:
            Mxx, Myy, Mxy = moment_at_gp[igp]
            # Effective stress for w-w amplification from M
            # ∫_{-t/2}^{t/2} (12*M*z/t^3)^2 / (N/t) dz = 12*M^2/(N*t^2)
            # (when N != 0)
            if abs(Nxx) > 1e-20:
                S_M_amp = np.array([
                    [12 * Mxx**2 / (Nxx * thickness**2), 0.0],
                    [0.0, 0.0],
                ])
                if abs(Nyy) > 1e-20:
                    S_M_amp[1, 1] = 12 * Myy**2 / (Nyy * thickness**2)
                kg_M_amp = (G.T @ S_M_amp @ G) * det_J * wt
                Kg[np.ix_(w_idx, w_idx)] += kg_M_amp

    return Kg


def mitc4_pressure_load(x: np.ndarray, y: np.ndarray,
                        pressures: np.ndarray,
                        direction: np.ndarray | None = None) -> np.ndarray:
    """Compute the 20x1 equivalent nodal force vector for surface pressure.

    Implements PLOAD4-style loading: pressure interpolated bilinearly over
    the element surface, acting in either the element normal direction or
    a user-specified direction (CID field on PLOAD4).

    f_i = integral[ N_i * p(xi,eta) * |J| ] dxi deta * direction_local

    Parameters
    ----------
    x, y : (4,) ndarray
        Nodal coordinates in element local frame.
    pressures : (4,) ndarray or float
        Pressure at each node [p1, p2, p3, p4].
        If scalar, uniform pressure on all nodes.
    direction : (3,) ndarray or None
        Unit load direction in element local frame [dx, dy, dz].
        If None, pressure acts in +z (element normal) only.
        When CID is used, the caller transforms the CID direction vector
        into the element local frame before passing it here.

    Returns
    -------
    Fe : (20,) ndarray
        Equivalent nodal force vector. DOF order per node: [u, v, w, tx, ty].
    """
    pressures = np.atleast_1d(np.asarray(pressures, dtype=float))
    if pressures.size == 1:
        pressures = np.full(4, pressures[0])

    Fe = np.zeros(20)

    if direction is None:
        # Normal pressure: only w DOF
        for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
            xi, eta = gpt
            N = _shape_quad4(xi, eta)
            dN_dnat = _dshape_quad4(xi, eta)
            _, det_J = _jacobian(dN_dnat, x, y)
            p_gp = N @ pressures
            for i in range(4):
                Fe[5 * i + 2] += N[i] * p_gp * det_J * wt
    else:
        # Directional pressure: force = p * direction at each node
        dx, dy, dz = direction[0], direction[1], direction[2]
        for gpt, wt in zip(_GAUSS_2x2_PTS, _GAUSS_2x2_WTS):
            xi, eta = gpt
            N = _shape_quad4(xi, eta)
            dN_dnat = _dshape_quad4(xi, eta)
            _, det_J = _jacobian(dN_dnat, x, y)
            p_gp = N @ pressures
            for i in range(4):
                contrib = N[i] * p_gp * det_J * wt
                Fe[5 * i + 0] += contrib * dx  # u
                Fe[5 * i + 1] += contrib * dy  # v
                Fe[5 * i + 2] += contrib * dz  # w

    return Fe


def _get_shear_constitutive(pid_ref, thickness: float) -> np.ndarray:
    """Get the 2x2 transverse shear constitutive matrix from PSHELL properties.

    Ds = kappa * G * ts * [[1, 0], [0, 1]]

    where kappa = 5/6, ts = TS/T * thickness.
    """
    kappa = 5.0 / 6.0
    ptype = pid_ref.type
    if ptype == 'PSHELL':
        # MID3 controls transverse shear
        if pid_ref.mid3 is None or pid_ref.mid3 == 0:
            # No transverse shear flexibility: return large stiffness (Kirchhoff limit)
            # In practice, Nastran still uses a finite shear stiffness
            mid_ref = pid_ref.mid1_ref if pid_ref.mid1_ref is not None else pid_ref.mid2_ref
            G = mid_ref.G()
            ts_ratio = 1.0
        else:
            mid3_ref = pid_ref.mid3_ref
            G = mid3_ref.G()
            ts_ratio = pid_ref.tst if pid_ref.tst is not None else 0.8333333
        ts = ts_ratio * thickness
        Ds = kappa * G * ts * np.eye(2)
    elif ptype == 'PCOMP':
        # For composites, use weighted average shear modulus through thickness
        # Simplified: use the first ply's G as representative
        mid_ref = pid_ref.mids_ref[0]
        G = mid_ref.G()
        ts = kappa * thickness
        Ds = G * ts * np.eye(2)
    else:
        raise NotImplementedError(f'shear constitutive for {ptype}')
    return Ds


def _apply_shell_offset(Ke: np.ndarray, z0: float) -> np.ndarray:
    """Apply shell offset eccentricity transformation to element stiffness.

    When the shell midsurface is offset from the grid point plane by z0 in the
    element normal (+z local) direction, grid-point DOFs relate to midsurface
    DOFs via:
        u_mid = u_grid + z0 * theta_y
        v_mid = v_grid - z0 * theta_x
        w_mid = w_grid
        theta_x_mid = theta_x_grid
        theta_y_mid = theta_y_grid

    The eccentricity matrix E (per node, 5x5) transforms:
        q_mid = E @ q_grid

    and the offset stiffness at the grid level is:
        K_grid = E^T @ K_mid @ E

    Parameters
    ----------
    Ke : (20, 20) ndarray
        Element stiffness at midsurface level.
    z0 : float
        Offset distance from grid plane to midsurface (positive = +z local).

    Returns
    -------
    Ke_offset : (20, 20) ndarray
        Element stiffness at grid point level.
    """
    E_full = np.eye(20)
    for i in range(4):
        r = 5 * i
        # u_mid = u_grid + z0 * theta_y  -> E[r+0, r+4] = z0
        E_full[r + 0, r + 4] = z0
        # v_mid = v_grid - z0 * theta_x  -> E[r+1, r+3] = -z0
        E_full[r + 1, r + 3] = -z0
    return E_full.T @ Ke @ E_full


def build_kbb_cquad4_mitc4(model: BDF,
                           Kbb,
                           dof_map: DOF_MAP,
                           all_nids,
                           xyz_cid0: NDArrayN3float,
                           idtype: str = 'int32',
                           fdtype: str = 'float64') -> int:
    """Assemble CQUAD4 element stiffness using MITC4 formulation.

    Parameters
    ----------
    model : BDF
        The BDF model with cross-referenced properties.
    Kbb : sparse or dense global stiffness matrix
        Modified in place.
    dof_map : dict
        Maps (nid, component) to global DOF index.
    all_nids : array
        Sorted array of all node IDs.
    xyz_cid0 : (nnodes, 3) array
        Nodal coordinates in basic frame.

    Returns
    -------
    nelements : int
        Number of CQUAD4 elements processed.
    """
    eids = np.array(model._type_to_id_map['CQUAD4'], dtype=idtype)
    nelements = len(eids)
    if nelements == 0:
        return nelements

    eids.sort()

    # Gather coordinate system data for material orientation
    ncoords = len(model.coords)
    cids = np.zeros(ncoords, dtype='int32')
    iaxes = np.zeros((ncoords, 3), dtype='float64')
    for icid, (cid, coord) in zip(count(), sorted(model.coords.items())):
        cids[icid] = cid
        iaxes[icid, :] = coord.i

    # Gather element data
    theta_mcid = []
    nids = np.zeros((nelements, 4), dtype='int32')
    for i, eid in enumerate(eids):
        elem: CQUAD4 = model.elements[eid]
        nids[i, :] = elem.nodes
        theta_mcid.append(elem.theta_mcid)

    # Nodal positions
    inids = np.searchsorted(all_nids, nids.ravel()).reshape(nelements, 4)
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]
    p4 = xyz_cid0[inids[:, 3], :]

    # Element normals
    v13 = p1 - p3
    v24 = p2 - p4
    normal = np.cross(v13, v24)
    ni = np.linalg.norm(normal, axis=1)
    if np.any(ni == 0.0):
        ibad = np.where(ni == 0.0)[0]
        raise RuntimeError(f'elements={eids[ibad].tolist()} have zero-area normals')
    normal /= ni[:, np.newaxis]

    # Element x-axis (node 1 to node 2)
    iaxes_elem = p2 - p1

    # Material coordinate transformation (3x3 rotation matrices)
    T = transform_shell_material_coordinate_system(
        cids, iaxes_elem, theta_mcid, normal, p1, p2,
        idtype=idtype, fdtype=fdtype)

    # Process each element
    for i_elem, (eid, Ti) in enumerate(zip(eids, T)):
        elem: CQUAD4 = model.elements[eid]
        pid_ref = elem.pid_ref

        # Transform node positions to element local frame
        xy1 = Ti @ p1[i_elem]
        xy2 = Ti @ p2[i_elem]
        xy3 = Ti @ p3[i_elem]
        xy4 = Ti @ p4[i_elem]

        x_local = np.array([xy1[0], xy2[0], xy3[0], xy4[0]])
        y_local = np.array([xy1[1], xy2[1], xy3[1], xy4[1]])

        # Get constitutive matrices
        ptype = pid_ref.type
        if ptype == 'PSHELL':
            A_mat, B_mat, D_mat = pid_ref.get_individual_ABD_matrices()
            thickness = pid_ref.t
        elif ptype == 'PCOMP':
            A_mat, B_mat, D_mat = pid_ref.get_individual_ABD_matrices()
            thickness = pid_ref.Thickness()
        else:
            raise NotImplementedError(f'CQUAD4 property type {ptype}')

        Ds_mat = _get_shear_constitutive(pid_ref, thickness)

        # Check if coupling exists (from MID4 or PCOMP asymmetry)
        B_coupling = B_mat if np.any(np.abs(B_mat) > 1e-12) else None

        # Compute element stiffness in local frame (at midsurface)
        Ke_local = mitc4_stiffness(x_local, y_local, A_mat, D_mat, Ds_mat, B_coupling)

        # Shell offset (ZOFFS field on CQUAD4)
        z0 = elem.zoffset
        if abs(z0) > 0.0:
            Ke_local = _apply_shell_offset(Ke_local, z0)

        # Correct rotation transform: local [tx, ty] from global [rx, ry, rz]
        T_node_rot = Ti[0:2, :]  # (2, 3): rows 0,1 of Ti applied to rotation vector

        T_node = np.zeros((5, 6))
        T_node[0:3, 0:3] = Ti          # u,v,w from ux,uy,uz
        T_node[3:5, 3:6] = T_node_rot  # tx,ty from rx,ry,rz

        # Block-diagonal 20x24 transformation
        T_full = np.zeros((20, 24))
        for inode in range(4):
            r0 = 5 * inode
            c0 = 6 * inode
            T_full[r0:r0 + 5, c0:c0 + 6] = T_node

        # Transform to global: K_global_24x24 = T_full^T @ Ke_local @ T_full
        Ke_global = T_full.T @ Ke_local @ T_full

        # Assemble into global Kbb
        nid1, nid2, nid3, nid4 = elem.nodes
        node_ids = [nid1, nid2, nid3, nid4]
        global_dofs = []
        for nid in node_ids:
            for comp in range(1, 7):  # DOF 1-6
                global_dofs.append(dof_map[(nid, comp)])

        for ii in range(24):
            gi = global_dofs[ii]
            for jj in range(24):
                gj = global_dofs[jj]
                val = Ke_global[ii, jj]
                if abs(val) > 0.0:
                    Kbb[gi, gj] += val

    return nelements
