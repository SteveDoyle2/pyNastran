"""CHEXA element stiffness, mass, and geometric stiffness matrices.

8-node (trilinear) and 20-node (serendipity quadratic) isoparametric hexahedral
solid elements with full Gauss quadrature.

DOF ordering: 3 translational DOF per node (u, v, w).

References
----------
Bathe, K.J. (2014). "Finite Element Procedures", 2nd ed., Chapter 5.
Zienkiewicz, O.C. & Taylor, R.L. (2000). "The Finite Element Method", Vol. 1.
MYSTRAN Source: EMG5/HEXA.f90, EMG7/SHP3DH.f90
"""

from __future__ import annotations

import numpy as np


def _hexa8_shape_functions(
    xi: float, eta: float, mu: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute shape functions and derivatives for 8-node hex.

    Parameters
    ----------
    xi, eta, mu : float
        Natural coordinates in [-1, 1].

    Returns
    -------
    N : (8,) shape functions
    dNdxi : (8,) derivatives w.r.t. xi
    dNdeta : (8,) derivatives w.r.t. eta
    dNdmu : (8,) derivatives w.r.t. mu
    """
    xi_n = np.array([-1, 1, 1, -1, -1, 1, 1, -1], dtype=float)
    eta_n = np.array([-1, -1, 1, 1, -1, -1, 1, 1], dtype=float)
    mu_n = np.array([-1, -1, -1, -1, 1, 1, 1, 1], dtype=float)

    N = 0.125 * (1 + xi_n * xi) * (1 + eta_n * eta) * (1 + mu_n * mu)
    dNdxi = 0.125 * xi_n * (1 + eta_n * eta) * (1 + mu_n * mu)
    dNdeta = 0.125 * (1 + xi_n * xi) * eta_n * (1 + mu_n * mu)
    dNdmu = 0.125 * (1 + xi_n * xi) * (1 + eta_n * eta) * mu_n

    return N, dNdxi, dNdeta, dNdmu


def _isotropic_constitutive(E: float, nu: float) -> np.ndarray:
    """Build 6x6 isotropic elasticity matrix (Voigt notation).

    Ordering: [sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz].

    Parameters
    ----------
    E : float
        Young's modulus.
    nu : float
        Poisson's ratio.

    Returns
    -------
    D : (6, 6) constitutive matrix
    """
    coeff = E / ((1 + nu) * (1 - 2 * nu))
    d11 = 1 - nu
    d12 = nu
    d44 = (1 - 2 * nu) / 2
    D = coeff * np.array(
        [
            [d11, d12, d12, 0, 0, 0],
            [d12, d11, d12, 0, 0, 0],
            [d12, d12, d11, 0, 0, 0],
            [0, 0, 0, d44, 0, 0],
            [0, 0, 0, 0, d44, 0],
            [0, 0, 0, 0, 0, d44],
        ]
    )
    return D


def hexa8_stiffness(coords: np.ndarray, D: np.ndarray) -> np.ndarray:
    """Compute the 24x24 CHEXA8 element stiffness matrix.

    Parameters
    ----------
    coords : (8, 3)
        Nodal coordinates [x, y, z] for each of the 8 nodes.
    D : (6, 6)
        Material constitutive matrix (Voigt notation).

    Returns
    -------
    Ke : (24, 24) element stiffness matrix
    """
    k = 1.0 / np.sqrt(3.0)
    gauss_pts = np.array(
        [
            [-k, -k, -k],
            [k, -k, -k],
            [k, k, -k],
            [-k, k, -k],
            [-k, -k, k],
            [k, -k, k],
            [k, k, k],
            [-k, k, k],
        ]
    )

    Ke = np.zeros((24, 24))
    ix = np.arange(0, 24, 3)

    for gpt in gauss_pts:
        xi, eta, mu = gpt
        N, dNdxi, dNdeta, dNdmu = _hexa8_shape_functions(xi, eta, mu)

        # Jacobian
        J = np.array([dNdxi, dNdeta, dNdmu]) @ coords
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)

        # Shape function derivatives in physical coordinates
        dN_dxyz = Jinv @ np.array([dNdxi, dNdeta, dNdmu])
        Nx = dN_dxyz[0]
        Ny = dN_dxyz[1]
        Nz = dN_dxyz[2]

        # B-matrix (6x24)
        B = np.zeros((6, 24))
        B[0, ix] = Nx
        B[1, ix + 1] = Ny
        B[2, ix + 2] = Nz
        B[3, ix] = Ny
        B[3, ix + 1] = Nx
        B[4, ix + 1] = Nz
        B[4, ix + 2] = Ny
        B[5, ix] = Nz
        B[5, ix + 2] = Nx

        Ke += (B.T @ D @ B) * detJ

    return Ke


def hexa8_mass(coords: np.ndarray, rho: float) -> np.ndarray:
    """Compute the 24x24 CHEXA8 consistent mass matrix.

    Parameters
    ----------
    coords : (8, 3)
        Nodal coordinates.
    rho : float
        Material density.

    Returns
    -------
    Me : (24, 24) consistent mass matrix
    """
    k = 1.0 / np.sqrt(3.0)
    gauss_pts = np.array(
        [
            [-k, -k, -k],
            [k, -k, -k],
            [k, k, -k],
            [-k, k, -k],
            [-k, -k, k],
            [k, -k, k],
            [k, k, k],
            [-k, k, k],
        ]
    )

    Me = np.zeros((24, 24))

    for gpt in gauss_pts:
        xi, eta, mu = gpt
        N, dNdxi, dNdeta, dNdmu = _hexa8_shape_functions(xi, eta, mu)

        # Jacobian determinant
        J = np.array([dNdxi, dNdeta, dNdmu]) @ coords
        detJ = np.linalg.det(J)

        # N-matrix (3x24): interpolation of displacements
        Ni = np.zeros((3, 24))
        ix = np.arange(0, 24, 3)
        Ni[0, ix] = N
        Ni[1, ix + 1] = N
        Ni[2, ix + 2] = N

        Me += rho * (Ni.T @ Ni) * detJ

    return Me


def hexa8_volume(coords: np.ndarray) -> float:
    """Compute the volume of a CHEXA8 element via 2x2x2 Gauss quadrature.

    Parameters
    ----------
    coords : (8, 3)
        Nodal coordinates.

    Returns
    -------
    volume : float
    """
    k = 1.0 / np.sqrt(3.0)
    gauss_pts = np.array(
        [
            [-k, -k, -k],
            [k, -k, -k],
            [k, k, -k],
            [-k, k, -k],
            [-k, -k, k],
            [k, -k, k],
            [k, k, k],
            [-k, k, k],
        ]
    )

    volume = 0.0
    for gpt in gauss_pts:
        xi, eta, mu = gpt
        _, dNdxi, dNdeta, dNdmu = _hexa8_shape_functions(xi, eta, mu)
        J = np.array([dNdxi, dNdeta, dNdmu]) @ coords
        volume += np.linalg.det(J)
    return volume


def hexa8_geometric_stiffness(
    coords: np.ndarray,
    stress: np.ndarray,
) -> np.ndarray:
    """Compute the 24x24 CHEXA8 geometric stiffness matrix.

    Parameters
    ----------
    coords : (8, 3)
        Nodal coordinates.
    stress : (6,)
        Element stress [sxx, syy, szz, sxy, syz, sxz] from preload.

    Returns
    -------
    KDe : (24, 24) geometric stiffness matrix
    """
    k = 1.0 / np.sqrt(3.0)
    gauss_pts = np.array(
        [
            [-k, -k, -k],
            [k, -k, -k],
            [k, k, -k],
            [-k, k, -k],
            [-k, -k, k],
            [k, -k, k],
            [k, k, k],
            [-k, k, k],
        ]
    )
    ndof = 24
    nnodes = 8
    KDe = np.zeros((ndof, ndof))

    sxx, syy, szz, sxy, syz, sxz = stress
    sigma = np.array(
        [
            [sxx, sxy, sxz],
            [sxy, syy, syz],
            [sxz, syz, szz],
        ]
    )

    for gpt in gauss_pts:
        xi, eta, mu = gpt
        N, dNdxi, dNdeta, dNdmu = _hexa8_shape_functions(xi, eta, mu)
        J = np.array([dNdxi, dNdeta, dNdmu]) @ coords
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        dN_dxyz = Jinv @ np.array([dNdxi, dNdeta, dNdmu])

        for i in range(nnodes):
            for j in range(nnodes):
                val = (dN_dxyz[:, i] @ sigma @ dN_dxyz[:, j]) * detJ
                KDe[3 * i, 3 * j] += val
                KDe[3 * i + 1, 3 * j + 1] += val
                KDe[3 * i + 2, 3 * j + 2] += val

    return KDe


# ============================================================================
# HEXA20 (20-node serendipity quadratic hexahedron)
# ============================================================================


def _hexa20_shape_functions(
    xi: float,
    eta: float,
    mu: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute shape functions and derivatives for 20-node hex.

    Parameters
    ----------
    xi, eta, mu : float
        Natural coordinates in [-1, 1].

    Returns
    -------
    N : (20,) shape functions
    dNdnat : (3, 20) derivatives w.r.t. (xi, eta, mu)
    """
    # Corner node natural coordinates
    xi_c = np.array([-1, 1, 1, -1, -1, 1, 1, -1], dtype=float)
    eta_c = np.array([-1, -1, 1, 1, -1, -1, 1, 1], dtype=float)
    mu_c = np.array([-1, -1, -1, -1, 1, 1, 1, 1], dtype=float)

    N = np.zeros(20)
    dNdxi_arr = np.zeros(20)
    dNdeta_arr = np.zeros(20)
    dNdmu_arr = np.zeros(20)

    # Corner nodes (1-8): serendipity quadratic formula
    for i in range(8):
        xii = xi_c[i]
        etai = eta_c[i]
        mui = mu_c[i]
        xp = 1.0 + xii * xi
        ep = 1.0 + etai * eta
        mp = 1.0 + mui * mu
        s = xii * xi + etai * eta + mui * mu - 2.0

        N[i] = 0.125 * xp * ep * mp * s
        dNdxi_arr[i] = 0.125 * (xii * ep * mp * s + xp * ep * mp * xii)
        dNdeta_arr[i] = 0.125 * (xp * etai * mp * s + xp * ep * mp * etai)
        dNdmu_arr[i] = 0.125 * (xp * ep * mui * s + xp * ep * mp * mui)

    # Midside nodes on edges parallel to xi (nodes 9, 11, 17, 19 -> indices 8,10,16,18)
    # eta_m, mu_m for these nodes
    mid_xi = [(8, -1, -1), (10, 1, -1), (16, -1, 1), (18, 1, 1)]
    for idx, etai, mui in mid_xi:
        ep = 1.0 + etai * eta
        mp = 1.0 + mui * mu
        N[idx] = 0.25 * (1.0 - xi**2) * ep * mp
        dNdxi_arr[idx] = 0.25 * (-2.0 * xi) * ep * mp
        dNdeta_arr[idx] = 0.25 * (1.0 - xi**2) * etai * mp
        dNdmu_arr[idx] = 0.25 * (1.0 - xi**2) * ep * mui

    # Midside nodes on edges parallel to eta (nodes 10, 12, 18, 20 -> indices 9,11,17,19)
    mid_eta = [(9, -1, -1), (11, 1, -1), (17, -1, 1), (19, 1, 1)]
    for idx, xii, mui in mid_eta:
        xp = 1.0 + xii * xi
        mp = 1.0 + mui * mu
        N[idx] = 0.25 * xp * (1.0 - eta**2) * mp
        dNdxi_arr[idx] = 0.25 * xii * (1.0 - eta**2) * mp
        dNdeta_arr[idx] = 0.25 * xp * (-2.0 * eta) * mp
        dNdmu_arr[idx] = 0.25 * xp * (1.0 - eta**2) * mui

    # Midside nodes on edges parallel to mu (nodes 13-16 -> indices 12-15)
    mid_mu = [(12, -1, -1), (13, 1, -1), (14, 1, 1), (15, -1, 1)]
    for idx, xii, etai in mid_mu:
        xp = 1.0 + xii * xi
        ep = 1.0 + etai * eta
        N[idx] = 0.25 * xp * ep * (1.0 - mu**2)
        dNdxi_arr[idx] = 0.25 * xii * ep * (1.0 - mu**2)
        dNdeta_arr[idx] = 0.25 * xp * etai * (1.0 - mu**2)
        dNdmu_arr[idx] = 0.25 * xp * ep * (-2.0 * mu)

    dNdnat = np.array([dNdxi_arr, dNdeta_arr, dNdmu_arr])
    return N, dNdnat


def hexa20_stiffness(coords: np.ndarray, D: np.ndarray) -> np.ndarray:
    """Compute the 60x60 CHEXA20 element stiffness matrix.

    Parameters
    ----------
    coords : (20, 3)
        Nodal coordinates.
    D : (6, 6)
        Material constitutive matrix.

    Returns
    -------
    Ke : (60, 60) element stiffness matrix
    """
    from .tetra import _solid_B_matrix

    k = np.sqrt(3.0 / 5.0)
    pts_1d = np.array([-k, 0.0, k])
    wts_1d = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])

    ndof = 60
    Ke = np.zeros((ndof, ndof))

    for i in range(3):
        for j in range(3):
            for m in range(3):
                xi, eta, mu = pts_1d[i], pts_1d[j], pts_1d[m]
                w = wts_1d[i] * wts_1d[j] * wts_1d[m]

                N, dNdnat = _hexa20_shape_functions(xi, eta, mu)
                J = dNdnat @ coords
                detJ = np.linalg.det(J)
                Jinv = np.linalg.inv(J)
                dN_dxyz = Jinv @ dNdnat
                B = _solid_B_matrix(dN_dxyz, 20)
                Ke += (B.T @ D @ B) * detJ * w

    return Ke


def hexa20_mass(coords: np.ndarray, rho: float) -> np.ndarray:
    """Compute the 60x60 CHEXA20 consistent mass matrix.

    Parameters
    ----------
    coords : (20, 3)
        Nodal coordinates.
    rho : float
        Material density.

    Returns
    -------
    Me : (60, 60) consistent mass matrix
    """
    k = np.sqrt(3.0 / 5.0)
    pts_1d = np.array([-k, 0.0, k])
    wts_1d = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])

    ndof = 60
    Me = np.zeros((ndof, ndof))

    for i in range(3):
        for j in range(3):
            for m in range(3):
                xi, eta, mu = pts_1d[i], pts_1d[j], pts_1d[m]
                w = wts_1d[i] * wts_1d[j] * wts_1d[m]

                N, dNdnat = _hexa20_shape_functions(xi, eta, mu)
                J = dNdnat @ coords
                detJ = np.linalg.det(J)

                Ni = np.zeros((3, ndof))
                ix = np.arange(0, ndof, 3)
                Ni[0, ix] = N
                Ni[1, ix + 1] = N
                Ni[2, ix + 2] = N
                Me += rho * (Ni.T @ Ni) * detJ * w

    return Me


def hexa20_geometric_stiffness(
    coords: np.ndarray,
    stress: np.ndarray,
) -> np.ndarray:
    """Compute the 60x60 CHEXA20 geometric stiffness matrix.

    Parameters
    ----------
    coords : (20, 3)
        Nodal coordinates.
    stress : (6,)
        Element stress [sxx, syy, szz, sxy, syz, sxz] from preload.

    Returns
    -------
    KDe : (60, 60) geometric stiffness matrix
    """
    k = np.sqrt(3.0 / 5.0)
    pts_1d = np.array([-k, 0.0, k])
    wts_1d = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])

    ndof = 60
    nnodes = 20
    KDe = np.zeros((ndof, ndof))

    sxx, syy, szz, sxy, syz, sxz = stress
    sigma = np.array(
        [
            [sxx, sxy, sxz],
            [sxy, syy, syz],
            [sxz, syz, szz],
        ]
    )

    for i in range(3):
        for j in range(3):
            for m in range(3):
                xi, eta, mu = pts_1d[i], pts_1d[j], pts_1d[m]
                w = wts_1d[i] * wts_1d[j] * wts_1d[m]

                N, dNdnat = _hexa20_shape_functions(xi, eta, mu)
                J = dNdnat @ coords
                detJ = np.linalg.det(J)
                Jinv = np.linalg.inv(J)
                dN_dxyz = Jinv @ dNdnat

                for ii in range(nnodes):
                    for jj in range(nnodes):
                        val = (dN_dxyz[:, ii] @ sigma @ dN_dxyz[:, jj]) * detJ * w
                        KDe[3 * ii, 3 * jj] += val
                        KDe[3 * ii + 1, 3 * jj + 1] += val
                        KDe[3 * ii + 2, 3 * jj + 2] += val

    return KDe
