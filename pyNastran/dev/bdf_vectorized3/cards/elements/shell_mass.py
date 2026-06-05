"""Consistent mass matrices for shell elements via Gauss quadrature.

The consistent mass matrix for a shell element is:
    M[3i+a, 3j+b] = delta_ab * integral(Ni * Nj * mass_per_area * |J| dA)

where Ni are 2D shape functions, J is the Jacobian mapping parametric to
physical surface coordinates, and mass_per_area = rho*t + nsm.

For shells with 6 DOF/node (3 translational + 3 rotational), only the
translational DOFs are populated here. Rotational inertia requires
additional thickness integration and is not included.

TODO: rotational DOF mass (drilling and bending inertia)
"""
from __future__ import annotations
import numpy as np


# --- 2D Shape functions and derivatives ---

def _shape_tri3(xi: float, eta: float) -> np.ndarray:
    """3-node linear triangle. L1=xi, L2=eta, L3=1-xi-eta. Returns (3,)."""
    return np.array([xi, eta, 1.0 - xi - eta])


def _dshape_tri3(xi: float, eta: float) -> np.ndarray:
    """Derivatives for 3-node triangle. Returns (3, 2): dN/dxi, dN/deta."""
    dN = np.zeros((3, 2))
    dN[0, 0] = 1.0;  dN[0, 1] = 0.0
    dN[1, 0] = 0.0;  dN[1, 1] = 1.0
    dN[2, 0] = -1.0; dN[2, 1] = -1.0
    return dN


def _shape_tri6(xi: float, eta: float) -> np.ndarray:
    """6-node quadratic triangle. Returns (6,)."""
    L1 = xi
    L2 = eta
    L3 = 1.0 - xi - eta
    N = np.zeros(6)
    N[0] = L1 * (2*L1 - 1)
    N[1] = L2 * (2*L2 - 1)
    N[2] = L3 * (2*L3 - 1)
    N[3] = 4*L1*L2
    N[4] = 4*L2*L3
    N[5] = 4*L1*L3
    return N


def _dshape_tri6(xi: float, eta: float) -> np.ndarray:
    """Derivatives for 6-node triangle. Returns (6, 2)."""
    L1 = xi
    L2 = eta
    L3 = 1.0 - xi - eta
    dN = np.zeros((6, 2))
    # dN/dxi (dL1/dxi=1, dL3/dxi=-1)
    dN[0, 0] = 4*L1 - 1
    dN[1, 0] = 0.0
    dN[2, 0] = -(4*L3 - 1)
    dN[3, 0] = 4*L2
    dN[4, 0] = -4*L2
    dN[5, 0] = 4*(L3 - L1)
    # dN/deta (dL2/deta=1, dL3/deta=-1)
    dN[0, 1] = 0.0
    dN[1, 1] = 4*L2 - 1
    dN[2, 1] = -(4*L3 - 1)
    dN[3, 1] = 4*L1
    dN[4, 1] = 4*(L3 - L2)
    dN[5, 1] = -4*L1
    return dN


def _shape_quad4(xi: float, eta: float) -> np.ndarray:
    """4-node bilinear quad. xi,eta in [-1,1]. Returns (4,)."""
    N = np.zeros(4)
    N[0] = (1 - xi) * (1 - eta) / 4.
    N[1] = (1 + xi) * (1 - eta) / 4.
    N[2] = (1 + xi) * (1 + eta) / 4.
    N[3] = (1 - xi) * (1 + eta) / 4.
    return N


def _dshape_quad4(xi: float, eta: float) -> np.ndarray:
    """Derivatives for 4-node quad. Returns (4, 2)."""
    dN = np.zeros((4, 2))
    dN[0, 0] = -(1 - eta) / 4.;  dN[0, 1] = -(1 - xi) / 4.
    dN[1, 0] = (1 - eta) / 4.;   dN[1, 1] = -(1 + xi) / 4.
    dN[2, 0] = (1 + eta) / 4.;   dN[2, 1] = (1 + xi) / 4.
    dN[3, 0] = -(1 + eta) / 4.;  dN[3, 1] = (1 - xi) / 4.
    return dN


def _shape_quad8(xi: float, eta: float) -> np.ndarray:
    """8-node serendipity quad. xi,eta in [-1,1]. Returns (8,)."""
    corners = [(-1,-1),(1,-1),(1,1),(-1,1)]
    midsides = [(0,-1),(1,0),(0,1),(-1,0)]

    N = np.zeros(8)
    for i, (xi_i, eta_i) in enumerate(corners):
        xp = 1 + xi_i * xi
        ep = 1 + eta_i * eta
        N[i] = xp * ep * (xi_i*xi + eta_i*eta - 1) / 4.

    for i, (xi_i, eta_i) in enumerate(midsides):
        if abs(xi_i) < 0.5:  # xi=0 edge
            N[4+i] = (1 - xi**2) * (1 + eta_i*eta) / 2.
        else:  # eta=0 edge
            N[4+i] = (1 + xi_i*xi) * (1 - eta**2) / 2.
    return N


def _dshape_quad8(xi: float, eta: float) -> np.ndarray:
    """Derivatives for 8-node serendipity quad. Returns (8, 2)."""
    corners = [(-1,-1),(1,-1),(1,1),(-1,1)]
    midsides = [(0,-1),(1,0),(0,1),(-1,0)]

    dN = np.zeros((8, 2))
    for i, (xi_i, eta_i) in enumerate(corners):
        xp = 1 + xi_i * xi
        ep = 1 + eta_i * eta
        s = xi_i*xi + eta_i*eta - 1
        dN[i, 0] = xi_i * ep * s / 4. + xp * ep * xi_i / 4.
        dN[i, 1] = eta_i * xp * s / 4. + xp * ep * eta_i / 4.

    for i, (xi_i, eta_i) in enumerate(midsides):
        if abs(xi_i) < 0.5:  # xi=0 edge
            dN[4+i, 0] = -2*xi * (1 + eta_i*eta) / 2.
            dN[4+i, 1] = (1 - xi**2) * eta_i / 2.
        else:  # eta=0 edge
            dN[4+i, 0] = xi_i * (1 - eta**2) / 2.
            dN[4+i, 1] = (1 + xi_i*xi) * (-2*eta) / 2.
    return dN


# --- 2D Gauss quadrature rules ---

def _tri_gauss_3pt() -> tuple[np.ndarray, np.ndarray]:
    """3-point triangle rule (degree 2). Sum of weights = 1/2 (reference area)."""
    pts = np.array([
        [1./6., 1./6.],
        [2./3., 1./6.],
        [1./6., 2./3.],
    ])
    wts = np.array([1./6., 1./6., 1./6.])
    return pts, wts


def _quad_gauss_2x2() -> tuple[np.ndarray, np.ndarray]:
    """2x2 Gauss-Legendre for quad [-1,1]^2."""
    g = 1.0 / np.sqrt(3.0)
    pts = np.array([[-g,-g],[g,-g],[g,g],[-g,g]])
    wts = np.array([1.0, 1.0, 1.0, 1.0])
    return pts, wts


def _quad_gauss_3x3() -> tuple[np.ndarray, np.ndarray]:
    """3x3 Gauss-Legendre for quad [-1,1]^2."""
    g = np.sqrt(3./5.)
    pts_1d = np.array([-g, 0.0, g])
    wts_1d = np.array([5./9., 8./9., 5./9.])
    pts = np.zeros((9, 2))
    wts = np.zeros(9)
    idx = 0
    for i in range(3):
        for j in range(3):
            pts[idx] = [pts_1d[i], pts_1d[j]]
            wts[idx] = wts_1d[i] * wts_1d[j]
            idx += 1
    return pts, wts


# --- Common 2D surface mass matrix integrator ---

def _consistent_mass_matrix_2d(
    all_nodes: np.ndarray,
    mass_per_area: np.ndarray,
    gauss_pts: np.ndarray,
    gauss_wts: np.ndarray,
    shape_func,
    dshape_func,
) -> np.ndarray:
    """Compute consistent mass matrix for shell elements via 2D Gauss quadrature.

    Parameters
    ----------
    all_nodes : (nelements, nnodes, 3) node coordinates in physical space
    mass_per_area : (nelements,) mass per unit area (rho*t + nsm)
    gauss_pts : (npts, 2) parametric coordinates
    gauss_wts : (npts,) weights
    shape_func : callable(xi, eta) -> (nnodes,)
    dshape_func : callable(xi, eta) -> (nnodes, 2)

    Returns
    -------
    M : (nelements, ndof, ndof) where ndof = 3*nnodes (translational DOFs only)
    """
    nelements = all_nodes.shape[0]
    nnodes = all_nodes.shape[1]
    ndof = 3 * nnodes

    M = np.zeros((nelements, ndof, ndof), dtype=all_nodes.dtype)

    for i_gp in range(len(gauss_wts)):
        xi, eta = gauss_pts[i_gp]
        w = gauss_wts[i_gp]

        N = shape_func(xi, eta)       # (nnodes,)
        dN = dshape_func(xi, eta)     # (nnodes, 2)

        # 2D Jacobian: maps (xi, eta) -> (x, y, z)
        # J = (nelements, 3, 2): J[:, i, j] = sum_k dN[k,j] * nodes[:, k, i]
        J = np.einsum('kj,eki->eij', dN, all_nodes)  # (nelements, 3, 2)

        # |J| = norm of cross product of the two tangent vectors
        # tangent1 = J[:, :, 0], tangent2 = J[:, :, 1]
        t1 = J[:, :, 0]  # (nelements, 3)
        t2 = J[:, :, 1]  # (nelements, 3)
        cross = np.cross(t1, t2, axis=1)  # (nelements, 3)
        det_J = np.linalg.norm(cross, axis=1)  # (nelements,)

        factor = mass_per_area * w * det_J  # (nelements,)

        # NNT outer product (same for all elements at this gauss point)
        NNT = np.outer(N, N)  # (nnodes, nnodes)

        for a in range(3):
            rows = np.arange(a, ndof, 3)
            cols = np.arange(a, ndof, 3)
            M[:, rows[:, None], cols[None, :]] += factor[:, None, None] * NNT[None, :, :]

    return M


# --- Public API ---

def consistent_mass_ctria3(all_nodes: np.ndarray, mass_per_area: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 3-node linear triangles.

    Parameters
    ----------
    all_nodes : (nelements, 3, 3)
    mass_per_area : (nelements,)

    Returns
    -------
    M : (nelements, 9, 9)
    """
    pts, wts = _tri_gauss_3pt()
    return _consistent_mass_matrix_2d(all_nodes, mass_per_area, pts, wts,
                                      _shape_tri3, _dshape_tri3)


def consistent_mass_ctria6(all_nodes: np.ndarray, mass_per_area: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 6-node quadratic triangles.

    Parameters
    ----------
    all_nodes : (nelements, 6, 3)
    mass_per_area : (nelements,)

    Returns
    -------
    M : (nelements, 18, 18)
    """
    pts, wts = _tri_gauss_3pt()
    return _consistent_mass_matrix_2d(all_nodes, mass_per_area, pts, wts,
                                      _shape_tri6, _dshape_tri6)


def consistent_mass_cquad4(all_nodes: np.ndarray, mass_per_area: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 4-node bilinear quads.

    Parameters
    ----------
    all_nodes : (nelements, 4, 3)
    mass_per_area : (nelements,)

    Returns
    -------
    M : (nelements, 12, 12)
    """
    pts, wts = _quad_gauss_2x2()
    return _consistent_mass_matrix_2d(all_nodes, mass_per_area, pts, wts,
                                      _shape_quad4, _dshape_quad4)


def consistent_mass_cquad8(all_nodes: np.ndarray, mass_per_area: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 8-node serendipity quads.

    Uses 3x3 Gauss quadrature for accuracy with quadratic shape functions.

    Parameters
    ----------
    all_nodes : (nelements, 8, 3)
    mass_per_area : (nelements,)

    Returns
    -------
    M : (nelements, 24, 24)
    """
    pts, wts = _quad_gauss_3x3()
    return _consistent_mass_matrix_2d(all_nodes, mass_per_area, pts, wts,
                                      _shape_quad8, _dshape_quad8)
