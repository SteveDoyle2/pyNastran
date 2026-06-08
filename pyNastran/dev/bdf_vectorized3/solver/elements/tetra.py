"""CTETRA element stiffness, mass, and geometric stiffness matrices.

4-node linear and 10-node quadratic isoparametric tetrahedral solid elements.
Uses volume-coordinate-based Gauss quadrature.

DOF ordering: 3 translational DOF per node (u, v, w).

References
----------
Bathe, K.J. (2014). "Finite Element Procedures", 2nd ed., Chapter 5.
Zienkiewicz, O.C. & Taylor, R.L. (2000). "The Finite Element Method", Vol. 1.
MYSTRAN Source: EMG5/TETRA.f90, EMG7/SHP3DT.f90, EMG7/ORDER_TETRA.f90
"""

from __future__ import annotations

import numpy as np


def _tetra4_shape_functions(
    xi: float,
    eta: float,
    zeta: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute shape functions and derivatives for 4-node tet.

    Parameters
    ----------
    xi, eta, zeta : float
        Volume coordinates in [0, 1] with xi + eta + zeta <= 1.

    Returns
    -------
    N : (4,) shape functions
    dNdnat : (3, 4) derivatives w.r.t. (xi, eta, zeta)
    """
    N = np.array([1.0 - xi - eta - zeta, xi, eta, zeta])
    dNdnat = np.array(
        [
            [-1.0, 1.0, 0.0, 0.0],
            [-1.0, 0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0, 1.0],
        ]
    )
    return N, dNdnat


def _tetra10_shape_functions(
    xi: float,
    eta: float,
    zeta: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute shape functions and derivatives for 10-node tet.

    Parameters
    ----------
    xi, eta, zeta : float
        Volume coordinates in [0, 1] with xi + eta + zeta <= 1.

    Returns
    -------
    N : (10,) shape functions
    dNdnat : (3, 10) derivatives w.r.t. (xi, eta, zeta)
    """
    phi = 1.0 - xi - eta - zeta

    N = np.array(
        [
            phi * (2.0 * phi - 1.0),
            xi * (2.0 * xi - 1.0),
            eta * (2.0 * eta - 1.0),
            zeta * (2.0 * zeta - 1.0),
            4.0 * xi * phi,
            4.0 * xi * eta,
            4.0 * eta * phi,
            4.0 * zeta * phi,
            4.0 * xi * zeta,
            4.0 * eta * zeta,
        ]
    )

    dNdxi = np.array(
        [
            -(4.0 * phi - 1.0),
            4.0 * xi - 1.0,
            0.0,
            0.0,
            4.0 * (phi - xi),
            4.0 * eta,
            -4.0 * eta,
            -4.0 * zeta,
            4.0 * zeta,
            0.0,
        ]
    )
    dNdeta = np.array(
        [
            -(4.0 * phi - 1.0),
            0.0,
            4.0 * eta - 1.0,
            0.0,
            -4.0 * xi,
            4.0 * xi,
            4.0 * (phi - eta),
            -4.0 * zeta,
            0.0,
            4.0 * zeta,
        ]
    )
    dNdzeta = np.array(
        [
            -(4.0 * phi - 1.0),
            0.0,
            0.0,
            4.0 * zeta - 1.0,
            -4.0 * xi,
            0.0,
            -4.0 * eta,
            4.0 * (phi - zeta),
            4.0 * xi,
            4.0 * eta,
        ]
    )
    dNdnat = np.array([dNdxi, dNdeta, dNdzeta])
    return N, dNdnat


# Gauss quadrature points for tetrahedra
_TETRA4_GAUSS = (
    np.array([[0.25, 0.25, 0.25]]),
    np.array([1.0 / 6.0]),
)

_alpha = 0.58541020
_beta = 0.13819660
_TETRA10_GAUSS = (
    np.array(
        [
            [_alpha, _beta, _beta],
            [_beta, _alpha, _beta],
            [_beta, _beta, _alpha],
            [_beta, _beta, _beta],
        ]
    ),
    np.array([1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0]),
)


def _solid_B_matrix(dN_dxyz: np.ndarray, nnodes: int) -> np.ndarray:
    """Build the 6 x (3*nnodes) strain-displacement matrix.

    Parameters
    ----------
    dN_dxyz : (3, nnodes)
        Shape function derivatives in physical coords.
    nnodes : int
        Number of element nodes.

    Returns
    -------
    B : (6, 3*nnodes)
    """
    ndof = 3 * nnodes
    B = np.zeros((6, ndof))
    ix = np.arange(0, ndof, 3)
    Nx = dN_dxyz[0]
    Ny = dN_dxyz[1]
    Nz = dN_dxyz[2]
    B[0, ix] = Nx
    B[1, ix + 1] = Ny
    B[2, ix + 2] = Nz
    B[3, ix] = Ny
    B[3, ix + 1] = Nx
    B[4, ix + 1] = Nz
    B[4, ix + 2] = Ny
    B[5, ix] = Nz
    B[5, ix + 2] = Nx
    return B


def tetra4_stiffness(coords: np.ndarray, D: np.ndarray) -> np.ndarray:
    """Compute the 12x12 CTETRA4 element stiffness matrix.

    Parameters
    ----------
    coords : (4, 3)
        Nodal coordinates [x, y, z] for each of the 4 nodes.
    D : (6, 6)
        Material constitutive matrix (Voigt notation).

    Returns
    -------
    Ke : (12, 12) element stiffness matrix
    """
    gauss_pts, weights = _TETRA4_GAUSS
    ndof = 12
    Ke = np.zeros((ndof, ndof))

    for gpt, w in zip(gauss_pts, weights):
        xi, eta, zeta = gpt
        N, dNdnat = _tetra4_shape_functions(xi, eta, zeta)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        dN_dxyz = Jinv @ dNdnat
        B = _solid_B_matrix(dN_dxyz, 4)
        Ke += (B.T @ D @ B) * detJ * w

    return Ke


def tetra10_stiffness(coords: np.ndarray, D: np.ndarray) -> np.ndarray:
    """Compute the 30x30 CTETRA10 element stiffness matrix.

    Parameters
    ----------
    coords : (10, 3)
        Nodal coordinates [x, y, z] for each of the 10 nodes.
    D : (6, 6)
        Material constitutive matrix (Voigt notation).

    Returns
    -------
    Ke : (30, 30) element stiffness matrix
    """
    gauss_pts, weights = _TETRA10_GAUSS
    ndof = 30
    Ke = np.zeros((ndof, ndof))

    for gpt, w in zip(gauss_pts, weights):
        xi, eta, zeta = gpt
        N, dNdnat = _tetra10_shape_functions(xi, eta, zeta)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        dN_dxyz = Jinv @ dNdnat
        B = _solid_B_matrix(dN_dxyz, 10)
        Ke += (B.T @ D @ B) * detJ * w

    return Ke


def tetra4_mass(coords: np.ndarray, rho: float) -> np.ndarray:
    """Compute the 12x12 CTETRA4 consistent mass matrix.

    Parameters
    ----------
    coords : (4, 3)
        Nodal coordinates.
    rho : float
        Material density.

    Returns
    -------
    Me : (12, 12) consistent mass matrix
    """
    gauss_pts, weights = _TETRA4_GAUSS
    ndof = 12
    Me = np.zeros((ndof, ndof))

    for gpt, w in zip(gauss_pts, weights):
        xi, eta, zeta = gpt
        N, dNdnat = _tetra4_shape_functions(xi, eta, zeta)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)

        Ni = np.zeros((3, ndof))
        ix = np.arange(0, ndof, 3)
        Ni[0, ix] = N
        Ni[1, ix + 1] = N
        Ni[2, ix + 2] = N
        Me += rho * (Ni.T @ Ni) * detJ * w

    return Me


def tetra10_mass(coords: np.ndarray, rho: float) -> np.ndarray:
    """Compute the 30x30 CTETRA10 consistent mass matrix.

    Parameters
    ----------
    coords : (10, 3)
        Nodal coordinates.
    rho : float
        Material density.

    Returns
    -------
    Me : (30, 30) consistent mass matrix
    """
    gauss_pts, weights = _TETRA10_GAUSS
    ndof = 30
    Me = np.zeros((ndof, ndof))

    for gpt, w in zip(gauss_pts, weights):
        xi, eta, zeta = gpt
        N, dNdnat = _tetra10_shape_functions(xi, eta, zeta)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)

        Ni = np.zeros((3, ndof))
        ix = np.arange(0, ndof, 3)
        Ni[0, ix] = N
        Ni[1, ix + 1] = N
        Ni[2, ix + 2] = N
        Me += rho * (Ni.T @ Ni) * detJ * w

    return Me


def tetra4_geometric_stiffness(
    coords: np.ndarray,
    stress: np.ndarray,
) -> np.ndarray:
    """Compute the 12x12 CTETRA4 geometric stiffness matrix.

    Parameters
    ----------
    coords : (4, 3)
        Nodal coordinates.
    stress : (6,)
        Element stress [sxx, syy, szz, sxy, syz, sxz] from preload.

    Returns
    -------
    KDe : (12, 12) geometric stiffness matrix
    """
    gauss_pts, weights = _TETRA4_GAUSS
    ndof = 12
    KDe = np.zeros((ndof, ndof))

    sxx, syy, szz, sxy, syz, sxz = stress
    sigma = np.array(
        [
            [sxx, sxy, sxz],
            [sxy, syy, syz],
            [sxz, syz, szz],
        ]
    )

    for gpt, w in zip(gauss_pts, weights):
        xi, eta, zeta = gpt
        N, dNdnat = _tetra4_shape_functions(xi, eta, zeta)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        dN_dxyz = Jinv @ dNdnat

        # G matrix: (3*nnodes) columns, arranged as [dN1/dx, dN1/dy, dN1/dz, dN2/dx, ...]
        # but for geometric stiffness we use: KD = integral(G^T * S_hat * G * detJ * w)
        # where S_hat is block-diagonal 3x3 sigma repeated nnodes times
        # Simpler: G(3, 3*4) with G[:, 3*i:3*i+3] = dNi/dxyz * I_3
        nnodes = 4
        G = np.zeros((3, ndof))
        for i in range(nnodes):
            G[0, 3 * i] = dN_dxyz[0, i]
            G[1, 3 * i + 1] = dN_dxyz[1, i]
            G[2, 3 * i + 2] = dN_dxyz[2, i]

        # ANSYS/MYSTRAN-style: sum over each pair of nodes
        # KD[3i:3i+3, 3j:3j+3] = dNi^T @ sigma @ dNj * detJ * w
        for i in range(nnodes):
            for j in range(nnodes):
                val = (dN_dxyz[:, i] @ sigma @ dN_dxyz[:, j]) * detJ * w
                KDe[3 * i, 3 * j] += val
                KDe[3 * i + 1, 3 * j + 1] += val
                KDe[3 * i + 2, 3 * j + 2] += val

    return KDe


def tetra10_geometric_stiffness(
    coords: np.ndarray,
    stress: np.ndarray,
) -> np.ndarray:
    """Compute the 30x30 CTETRA10 geometric stiffness matrix.

    Parameters
    ----------
    coords : (10, 3)
        Nodal coordinates.
    stress : (6,)
        Element stress [sxx, syy, szz, sxy, syz, sxz] from preload.

    Returns
    -------
    KDe : (30, 30) geometric stiffness matrix
    """
    gauss_pts, weights = _TETRA10_GAUSS
    ndof = 30
    nnodes = 10
    KDe = np.zeros((ndof, ndof))

    sxx, syy, szz, sxy, syz, sxz = stress
    sigma = np.array(
        [
            [sxx, sxy, sxz],
            [sxy, syy, syz],
            [sxz, syz, szz],
        ]
    )

    for gpt, w in zip(gauss_pts, weights):
        xi, eta, zeta = gpt
        N, dNdnat = _tetra10_shape_functions(xi, eta, zeta)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        dN_dxyz = Jinv @ dNdnat

        for i in range(nnodes):
            for j in range(nnodes):
                val = (dN_dxyz[:, i] @ sigma @ dN_dxyz[:, j]) * detJ * w
                KDe[3 * i, 3 * j] += val
                KDe[3 * i + 1, 3 * j + 1] += val
                KDe[3 * i + 2, 3 * j + 2] += val

    return KDe
