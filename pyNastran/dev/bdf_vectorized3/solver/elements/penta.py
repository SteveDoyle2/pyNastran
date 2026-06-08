"""CPENTA element stiffness, mass, and geometric stiffness matrices.

6-node linear and 15-node quadratic isoparametric pentahedral (wedge) elements.
Uses combined triangular + Gauss quadrature.

DOF ordering: 3 translational DOF per node (u, v, w).

References
----------
Bathe, K.J. (2014). "Finite Element Procedures", 2nd ed., Chapter 5.
Zienkiewicz, O.C. & Taylor, R.L. (2000). "The Finite Element Method", Vol. 1.
MYSTRAN Source: EMG5/PENTA.f90, EMG7/SHP3DP.f90, EMG7/ORDER_TRIA.f90
"""

from __future__ import annotations

import numpy as np

from .tetra import _solid_B_matrix


def _penta6_shape_functions(
    xi: float,
    eta: float,
    zeta: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute shape functions and derivatives for 6-node wedge.

    Parameters
    ----------
    xi, eta : float
        Triangular coordinates in [0, 1] with xi + eta <= 1.
    zeta : float
        Axial coordinate in [-1, 1].

    Returns
    -------
    N : (6,) shape functions
    dNdnat : (3, 6) derivatives w.r.t. (xi, eta, zeta)
    """
    phi = 1.0 - xi - eta

    N = np.array(
        [
            0.5 * phi * (1.0 - zeta),
            0.5 * xi * (1.0 - zeta),
            0.5 * eta * (1.0 - zeta),
            0.5 * phi * (1.0 + zeta),
            0.5 * xi * (1.0 + zeta),
            0.5 * eta * (1.0 + zeta),
        ]
    )

    dNdxi = np.array(
        [
            -0.5 * (1.0 - zeta),
            0.5 * (1.0 - zeta),
            0.0,
            -0.5 * (1.0 + zeta),
            0.5 * (1.0 + zeta),
            0.0,
        ]
    )
    dNdeta = np.array(
        [
            -0.5 * (1.0 - zeta),
            0.0,
            0.5 * (1.0 - zeta),
            -0.5 * (1.0 + zeta),
            0.0,
            0.5 * (1.0 + zeta),
        ]
    )
    dNdzeta = np.array(
        [
            -0.5 * phi,
            -0.5 * xi,
            -0.5 * eta,
            0.5 * phi,
            0.5 * xi,
            0.5 * eta,
        ]
    )
    dNdnat = np.array([dNdxi, dNdeta, dNdzeta])
    return N, dNdnat


def _penta15_shape_functions(
    xi: float,
    eta: float,
    zeta: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute shape functions and derivatives for 15-node wedge.

    Parameters
    ----------
    xi, eta : float
        Triangular coordinates in [0, 1] with xi + eta <= 1.
    zeta : float
        Axial coordinate in [-1, 1].

    Returns
    -------
    N : (15,) shape functions
    dNdnat : (3, 15) derivatives w.r.t. (xi, eta, zeta)

    Notes
    -----
    Node ordering: 1-6 corners (bottom tri 1-3, top tri 4-6),
    7-9 bottom edge midsides, 10-12 vertical edge midsides,
    13-15 top edge midsides.
    """
    phi = 1.0 - xi - eta
    zm = 1.0 - zeta
    zp = 1.0 + zeta

    N = np.zeros(15)
    # Corner nodes (quadratic in zeta, linear in triangle)
    N[0] = 0.5 * phi * zm * (2.0 * phi - 2.0 - zeta)
    N[1] = 0.5 * xi * zm * (2.0 * xi - 2.0 - zeta)
    N[2] = 0.5 * eta * zm * (2.0 * eta - 2.0 - zeta)
    N[3] = 0.5 * phi * zp * (2.0 * phi - 2.0 + zeta)
    N[4] = 0.5 * xi * zp * (2.0 * xi - 2.0 + zeta)
    N[5] = 0.5 * eta * zp * (2.0 * eta - 2.0 + zeta)
    # Bottom edge midsides
    N[6] = 2.0 * xi * phi * zm
    N[7] = 2.0 * xi * eta * zm
    N[8] = 2.0 * eta * phi * zm
    # Vertical edge midsides
    N[9] = phi * (1.0 - zeta**2)
    N[10] = xi * (1.0 - zeta**2)
    N[11] = eta * (1.0 - zeta**2)
    # Top edge midsides
    N[12] = 2.0 * xi * phi * zp
    N[13] = 2.0 * xi * eta * zp
    N[14] = 2.0 * eta * phi * zp

    dNdxi = np.zeros(15)
    dNdxi[0] = 0.5 * zm * (-(2.0 * phi - 2.0 - zeta) + phi * (-2.0))
    dNdxi[1] = 0.5 * zm * ((2.0 * xi - 2.0 - zeta) + xi * 2.0)
    dNdxi[2] = 0.0
    dNdxi[3] = 0.5 * zp * (-(2.0 * phi - 2.0 + zeta) + phi * (-2.0))
    dNdxi[4] = 0.5 * zp * ((2.0 * xi - 2.0 + zeta) + xi * 2.0)
    dNdxi[5] = 0.0
    dNdxi[6] = 2.0 * (phi - xi) * zm
    dNdxi[7] = 2.0 * eta * zm
    dNdxi[8] = -2.0 * eta * zm
    dNdxi[9] = -(1.0 - zeta**2)
    dNdxi[10] = 1.0 - zeta**2
    dNdxi[11] = 0.0
    dNdxi[12] = 2.0 * (phi - xi) * zp
    dNdxi[13] = 2.0 * eta * zp
    dNdxi[14] = -2.0 * eta * zp

    dNdeta = np.zeros(15)
    dNdeta[0] = 0.5 * zm * (-(2.0 * phi - 2.0 - zeta) + phi * (-2.0))
    dNdeta[1] = 0.0
    dNdeta[2] = 0.5 * zm * ((2.0 * eta - 2.0 - zeta) + eta * 2.0)
    dNdeta[3] = 0.5 * zp * (-(2.0 * phi - 2.0 + zeta) + phi * (-2.0))
    dNdeta[4] = 0.0
    dNdeta[5] = 0.5 * zp * ((2.0 * eta - 2.0 + zeta) + eta * 2.0)
    dNdeta[6] = -2.0 * xi * zm
    dNdeta[7] = 2.0 * xi * zm
    dNdeta[8] = 2.0 * (phi - eta) * zm
    dNdeta[9] = -(1.0 - zeta**2)
    dNdeta[10] = 0.0
    dNdeta[11] = 1.0 - zeta**2
    dNdeta[12] = -2.0 * xi * zp
    dNdeta[13] = 2.0 * xi * zp
    dNdeta[14] = 2.0 * (phi - eta) * zp

    dNdzeta = np.zeros(15)
    dNdzeta[0] = 0.5 * phi * (-(2.0 * phi - 2.0 - zeta) + zm * (-1.0))
    dNdzeta[1] = 0.5 * xi * (-(2.0 * xi - 2.0 - zeta) + zm * (-1.0))
    dNdzeta[2] = 0.5 * eta * (-(2.0 * eta - 2.0 - zeta) + zm * (-1.0))
    dNdzeta[3] = 0.5 * phi * ((2.0 * phi - 2.0 + zeta) + zp * 1.0)
    dNdzeta[4] = 0.5 * xi * ((2.0 * xi - 2.0 + zeta) + zp * 1.0)
    dNdzeta[5] = 0.5 * eta * ((2.0 * eta - 2.0 + zeta) + zp * 1.0)
    dNdzeta[6] = -2.0 * xi * phi
    dNdzeta[7] = -2.0 * xi * eta
    dNdzeta[8] = -2.0 * eta * phi
    dNdzeta[9] = -2.0 * phi * zeta
    dNdzeta[10] = -2.0 * xi * zeta
    dNdzeta[11] = -2.0 * eta * zeta
    dNdzeta[12] = 2.0 * xi * phi
    dNdzeta[13] = 2.0 * xi * eta
    dNdzeta[14] = 2.0 * eta * phi

    dNdnat = np.array([dNdxi, dNdeta, dNdzeta])
    return N, dNdnat


# Gauss quadrature: triangular (1-point) x axial (2-point) for PENTA6
_k1 = 1.0 / np.sqrt(3.0)
_PENTA6_GAUSS = (
    np.array(
        [
            [1.0 / 3.0, 1.0 / 3.0, -_k1],
            [1.0 / 3.0, 1.0 / 3.0, _k1],
        ]
    ),
    np.array([0.5, 0.5]),  # tri_weight(1/2) * axial_weight(1.0) = 0.5 each
)

# Gauss quadrature: triangular (3-point) x axial (2-point) for PENTA6 full
# Triangular weights (1/6 each, sum to 1/2 = reference triangle area)
# Axial Gauss weights (1.0 each for 2-point rule over [-1,1])
# Combined weight = tri_weight * axial_weight
_PENTA6_FULL_GAUSS = (
    np.array(
        [
            [1.0 / 6.0, 1.0 / 6.0, -_k1],
            [2.0 / 3.0, 1.0 / 6.0, -_k1],
            [1.0 / 6.0, 2.0 / 3.0, -_k1],
            [1.0 / 6.0, 1.0 / 6.0, _k1],
            [2.0 / 3.0, 1.0 / 6.0, _k1],
            [1.0 / 6.0, 2.0 / 3.0, _k1],
        ]
    ),
    np.array([1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]),
)

# PENTA15: triangular (3-point) x axial (3-point Simpson/Gauss)
_k3 = np.sqrt(3.0 / 5.0)
_w3 = np.array([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])
_tri3_pts = np.array([[1.0 / 6.0, 1.0 / 6.0], [2.0 / 3.0, 1.0 / 6.0], [1.0 / 6.0, 2.0 / 3.0]])
_tri3_w = np.array([1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0])
_zeta3 = np.array([-_k3, 0.0, _k3])

_penta15_pts = []
_penta15_wts = []
for _ti in range(3):
    for _zi in range(3):
        _penta15_pts.append([_tri3_pts[_ti, 0], _tri3_pts[_ti, 1], _zeta3[_zi]])
        _penta15_wts.append(_tri3_w[_ti] * _w3[_zi])

_PENTA15_GAUSS = (np.array(_penta15_pts), np.array(_penta15_wts))


def penta6_stiffness(coords: np.ndarray, D: np.ndarray) -> np.ndarray:
    """Compute the 18x18 CPENTA6 element stiffness matrix.

    Parameters
    ----------
    coords : (6, 3)
        Nodal coordinates.
    D : (6, 6)
        Material constitutive matrix.

    Returns
    -------
    Ke : (18, 18) element stiffness matrix
    """
    gauss_pts, weights = _PENTA6_FULL_GAUSS
    ndof = 18
    Ke = np.zeros((ndof, ndof))

    for gpt, w in zip(gauss_pts, weights):
        xi, eta, zeta = gpt
        N, dNdnat = _penta6_shape_functions(xi, eta, zeta)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        dN_dxyz = Jinv @ dNdnat
        B = _solid_B_matrix(dN_dxyz, 6)
        Ke += (B.T @ D @ B) * detJ * w

    return Ke


def penta15_stiffness(coords: np.ndarray, D: np.ndarray) -> np.ndarray:
    """Compute the 45x45 CPENTA15 element stiffness matrix.

    Parameters
    ----------
    coords : (15, 3)
        Nodal coordinates.
    D : (6, 6)
        Material constitutive matrix.

    Returns
    -------
    Ke : (45, 45) element stiffness matrix
    """
    gauss_pts, weights = _PENTA15_GAUSS
    ndof = 45
    Ke = np.zeros((ndof, ndof))

    for gpt, w in zip(gauss_pts, weights):
        xi, eta, zeta = gpt
        N, dNdnat = _penta15_shape_functions(xi, eta, zeta)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        dN_dxyz = Jinv @ dNdnat
        B = _solid_B_matrix(dN_dxyz, 15)
        Ke += (B.T @ D @ B) * detJ * w

    return Ke


def penta6_mass(coords: np.ndarray, rho: float) -> np.ndarray:
    """Compute the 18x18 CPENTA6 consistent mass matrix.

    Parameters
    ----------
    coords : (6, 3)
        Nodal coordinates.
    rho : float
        Material density.

    Returns
    -------
    Me : (18, 18) consistent mass matrix
    """
    gauss_pts, weights = _PENTA6_FULL_GAUSS
    ndof = 18
    Me = np.zeros((ndof, ndof))

    for gpt, w in zip(gauss_pts, weights):
        xi, eta, zeta = gpt
        N, dNdnat = _penta6_shape_functions(xi, eta, zeta)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)

        Ni = np.zeros((3, ndof))
        ix = np.arange(0, ndof, 3)
        Ni[0, ix] = N
        Ni[1, ix + 1] = N
        Ni[2, ix + 2] = N
        Me += rho * (Ni.T @ Ni) * detJ * w

    return Me


def penta15_mass(coords: np.ndarray, rho: float) -> np.ndarray:
    """Compute the 45x45 CPENTA15 consistent mass matrix.

    Parameters
    ----------
    coords : (15, 3)
        Nodal coordinates.
    rho : float
        Material density.

    Returns
    -------
    Me : (45, 45) consistent mass matrix
    """
    gauss_pts, weights = _PENTA15_GAUSS
    ndof = 45
    Me = np.zeros((ndof, ndof))

    for gpt, w in zip(gauss_pts, weights):
        xi, eta, zeta = gpt
        N, dNdnat = _penta15_shape_functions(xi, eta, zeta)
        J = dNdnat @ coords
        detJ = np.linalg.det(J)

        Ni = np.zeros((3, ndof))
        ix = np.arange(0, ndof, 3)
        Ni[0, ix] = N
        Ni[1, ix + 1] = N
        Ni[2, ix + 2] = N
        Me += rho * (Ni.T @ Ni) * detJ * w

    return Me


def penta6_geometric_stiffness(
    coords: np.ndarray,
    stress: np.ndarray,
) -> np.ndarray:
    """Compute the 18x18 CPENTA6 geometric stiffness matrix.

    Parameters
    ----------
    coords : (6, 3)
        Nodal coordinates.
    stress : (6,)
        Element stress [sxx, syy, szz, sxy, syz, sxz] from preload.

    Returns
    -------
    KDe : (18, 18) geometric stiffness matrix
    """
    gauss_pts, weights = _PENTA6_FULL_GAUSS
    ndof = 18
    nnodes = 6
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
        N, dNdnat = _penta6_shape_functions(xi, eta, zeta)
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


def penta15_geometric_stiffness(
    coords: np.ndarray,
    stress: np.ndarray,
) -> np.ndarray:
    """Compute the 45x45 CPENTA15 geometric stiffness matrix.

    Parameters
    ----------
    coords : (15, 3)
        Nodal coordinates.
    stress : (6,)
        Element stress [sxx, syy, szz, sxy, syz, sxz] from preload.

    Returns
    -------
    KDe : (45, 45) geometric stiffness matrix
    """
    gauss_pts, weights = _PENTA15_GAUSS
    ndof = 45
    nnodes = 15
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
        N, dNdnat = _penta15_shape_functions(xi, eta, zeta)
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
