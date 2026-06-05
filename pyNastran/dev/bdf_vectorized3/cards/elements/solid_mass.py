"""Consistent mass matrices for solid elements via Gauss quadrature.

The consistent mass matrix for a solid element is:
    M[3i+a, 3j+b] = rho * delta_ab * integral(Ni * Nj * det(J) dV)

where Ni are shape functions, J is the Jacobian matrix mapping parametric
to physical coordinates, and rho is the material density.

This gives a full (non-diagonal) mass matrix that exactly represents the
kinetic energy for the assumed displacement field.
"""
from __future__ import annotations
import numpy as np

from .solid_volume import (
    # Shape functions
    _shape_ctetra4, _shape_ctetra10,
    _shape_cpenta6, _shape_cpenta15,
    _shape_chexa8, _shape_chexa20,
    _shape_cpyram5, _shape_cpyram13,
    # Shape function derivatives (for Jacobian)
    _dshape_ctetra10, _dshape_cpenta15,
    _dshape_chexa20, _dshape_cpyram13,
    # Quadrature rules
    _tet_gauss_4pt, _penta_gauss_6pt, _gauss_pts_1d, _pyram_gauss_8pt,
)


def _consistent_mass_matrix(
    all_nodes: np.ndarray,
    rho: np.ndarray,
    gauss_pts: np.ndarray,
    gauss_wts: np.ndarray,
    shape_func,
    dshape_func,
    negate_detJ: bool = False,
) -> np.ndarray:
    """Compute consistent mass matrix via Gauss quadrature.

    Parameters
    ----------
    all_nodes : (nelements, nnodes, 3)
    rho : (nelements,) density per element
    gauss_pts : (npts, 3)
    gauss_wts : (npts,)
    shape_func : callable(xi, eta, zeta) -> (nnodes,)
    dshape_func : callable(xi, eta, zeta) -> (nnodes, 3)
    negate_detJ : bool
        If True, negate det(J) (needed for tet parametric convention)

    Returns
    -------
    M : (nelements, ndof, ndof) where ndof = 3*nnodes
    """
    nelements = all_nodes.shape[0]
    nnodes = all_nodes.shape[1]
    ndof = 3 * nnodes

    M = np.zeros((nelements, ndof, ndof), dtype=all_nodes.dtype)

    for i_gp in range(len(gauss_wts)):
        xi, eta, zeta = gauss_pts[i_gp]
        w = gauss_wts[i_gp]

        N = shape_func(xi, eta, zeta)      # (nnodes,)
        dN = dshape_func(xi, eta, zeta)    # (nnodes, 3)

        # Jacobian: J[e, i, j] = sum_k dN[k,j] * nodes[e,k,i]
        J = np.einsum('kj,eki->eji', dN, all_nodes)
        det_J = np.linalg.det(J)  # (nelements,)
        if negate_detJ:
            det_J = -det_J

        # Scalar integrand factor: rho * w * |det_J|
        factor = rho * w * det_J  # (nelements,)

        # NNT[i,j] = N[i]*N[j] -- outer product of shape functions
        NNT = np.outer(N, N)  # (nnodes, nnodes) -- same for all elements at this gauss point

        # Expand to DOF space: M[3i+a, 3j+a] += factor * NNT[i,j] for a=0,1,2
        for a in range(3):
            # indices for DOF direction a
            rows = np.arange(a, ndof, 3)  # [a, a+3, a+6, ...]
            cols = np.arange(a, ndof, 3)
            # NNT_expanded[ii, jj] for this direction
            M[:, rows[:, None], cols[None, :]] += factor[:, None, None] * NNT[None, :, :]

    return M


# --- Linear element derivatives (for Jacobian in mass matrix) ---

def _dshape_ctetra4(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape function derivatives for 4-node linear tet. Returns (4, 3)."""
    dN = np.zeros((4, 3))
    dN[0, 0] = 1.0   # dN1/dxi
    dN[1, 1] = 1.0   # dN2/deta
    dN[2, 2] = 1.0   # dN3/dzeta
    dN[3, :] = -1.0  # dN4/d* = -1
    return dN


def _dshape_cpenta6(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape function derivatives for 6-node linear wedge. Returns (6, 3)."""
    L1 = xi
    L2 = eta
    L3 = 1.0 - xi - eta
    zm = (1 - zeta) / 2.
    zp = (1 + zeta) / 2.
    dN = np.zeros((6, 3))
    # dN/dxi
    dN[0, 0] = zm;   dN[1, 0] = 0;    dN[2, 0] = -zm
    dN[3, 0] = zp;   dN[4, 0] = 0;    dN[5, 0] = -zp
    # dN/deta
    dN[0, 1] = 0;    dN[1, 1] = zm;   dN[2, 1] = -zm
    dN[3, 1] = 0;    dN[4, 1] = zp;   dN[5, 1] = -zp
    # dN/dzeta
    dN[0, 2] = -L1/2.; dN[1, 2] = -L2/2.; dN[2, 2] = -L3/2.
    dN[3, 2] = L1/2.;  dN[4, 2] = L2/2.;  dN[5, 2] = L3/2.
    return dN


def _dshape_chexa8(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape function derivatives for 8-node linear hex. Returns (8, 3)."""
    signs = [(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
             (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]
    dN = np.zeros((8, 3))
    for i, (si, ei, zi) in enumerate(signs):
        dN[i, 0] = si * (1 + ei*eta) * (1 + zi*zeta) / 8.
        dN[i, 1] = ei * (1 + si*xi) * (1 + zi*zeta) / 8.
        dN[i, 2] = zi * (1 + si*xi) * (1 + ei*eta) / 8.
    return dN


def _dshape_cpyram5(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape function derivatives for 5-node linear pyramid. Returns (5, 3)."""
    h = 1e-7
    dN = np.zeros((5, 3))
    dN[:, 0] = (_shape_cpyram5(xi+h, eta, zeta) - _shape_cpyram5(xi-h, eta, zeta)) / (2*h)
    dN[:, 1] = (_shape_cpyram5(xi, eta+h, zeta) - _shape_cpyram5(xi, eta-h, zeta)) / (2*h)
    dN[:, 2] = (_shape_cpyram5(xi, eta, zeta+h) - _shape_cpyram5(xi, eta, zeta-h)) / (2*h)
    return dN


# --- Public API ---

def consistent_mass_ctetra4(all_nodes: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 4-node linear tetrahedra.

    Parameters
    ----------
    all_nodes : (nelements, 4, 3)
    rho : (nelements,)

    Returns
    -------
    M : (nelements, 12, 12)
    """
    pts, wts = _tet_gauss_4pt()
    return _consistent_mass_matrix(
        all_nodes, rho, pts, wts,
        _shape_ctetra4, _dshape_ctetra4, negate_detJ=True)


def consistent_mass_ctetra10(all_nodes: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 10-node quadratic tetrahedra.

    Parameters
    ----------
    all_nodes : (nelements, 10, 3)
    rho : (nelements,)

    Returns
    -------
    M : (nelements, 30, 30)
    """
    pts, wts = _tet_gauss_4pt()
    return _consistent_mass_matrix(
        all_nodes, rho, pts, wts,
        _shape_ctetra10, _dshape_ctetra10, negate_detJ=True)


def consistent_mass_cpenta6(all_nodes: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 6-node linear wedge elements.

    Parameters
    ----------
    all_nodes : (nelements, 6, 3)
    rho : (nelements,)

    Returns
    -------
    M : (nelements, 18, 18)
    """
    pts, wts = _penta_gauss_6pt()
    return _consistent_mass_matrix(
        all_nodes, rho, pts, wts,
        _shape_cpenta6, _dshape_cpenta6)


def consistent_mass_cpenta15(all_nodes: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 15-node quadratic wedge elements.

    Parameters
    ----------
    all_nodes : (nelements, 15, 3)
    rho : (nelements,)

    Returns
    -------
    M : (nelements, 45, 45)
    """
    pts, wts = _penta_gauss_6pt()
    return _consistent_mass_matrix(
        all_nodes, rho, pts, wts,
        _shape_cpenta15, _dshape_cpenta15)


def consistent_mass_chexa8(all_nodes: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 8-node linear hexahedra.

    Parameters
    ----------
    all_nodes : (nelements, 8, 3)
    rho : (nelements,)

    Returns
    -------
    M : (nelements, 24, 24)
    """
    pts_1d, wts_1d = _gauss_pts_1d(2)
    pts = np.zeros((8, 3))
    wts = np.zeros(8)
    idx = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                pts[idx] = [pts_1d[i], pts_1d[j], pts_1d[k]]
                wts[idx] = wts_1d[i] * wts_1d[j] * wts_1d[k]
                idx += 1
    return _consistent_mass_matrix(
        all_nodes, rho, pts, wts,
        _shape_chexa8, _dshape_chexa8)


def consistent_mass_chexa20(all_nodes: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 20-node serendipity hexahedra.

    Uses 3x3x3 Gauss quadrature (27 points) for accuracy with quadratic functions.

    Parameters
    ----------
    all_nodes : (nelements, 20, 3)
    rho : (nelements,)

    Returns
    -------
    M : (nelements, 60, 60)
    """
    pts_1d, wts_1d = _gauss_pts_1d(3)
    pts = np.zeros((27, 3))
    wts = np.zeros(27)
    idx = 0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                pts[idx] = [pts_1d[i], pts_1d[j], pts_1d[k]]
                wts[idx] = wts_1d[i] * wts_1d[j] * wts_1d[k]
                idx += 1
    return _consistent_mass_matrix(
        all_nodes, rho, pts, wts,
        _shape_chexa20, _dshape_chexa20)


def consistent_mass_cpyram5(all_nodes: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 5-node linear pyramid.

    Parameters
    ----------
    all_nodes : (nelements, 5, 3)
    rho : (nelements,)

    Returns
    -------
    M : (nelements, 15, 15)
    """
    pts, wts = _pyram_gauss_8pt()
    return _consistent_mass_matrix(
        all_nodes, rho, pts, wts,
        _shape_cpyram5, _dshape_cpyram5)


def consistent_mass_cpyram13(all_nodes: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Consistent mass matrix for 13-node quadratic pyramid.

    Parameters
    ----------
    all_nodes : (nelements, 13, 3)
    rho : (nelements,)

    Returns
    -------
    M : (nelements, 39, 39)
    """
    pts, wts = _pyram_gauss_8pt()
    return _consistent_mass_matrix(
        all_nodes, rho, pts, wts,
        _shape_cpyram13, _dshape_cpyram13)
