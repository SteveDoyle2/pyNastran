"""Vectorized volume calculations for solid elements.

Linear elements use closed-form formulas.
Higher-order elements use Gauss quadrature with isoparametric shape functions.
"""
from __future__ import annotations
import numpy as np


# --- Gauss quadrature points and weights ---

def _gauss_pts_1d(n: int) -> tuple[np.ndarray, np.ndarray]:
    """1D Gauss-Legendre points and weights for n=2 or n=3."""
    if n == 2:
        g = 1.0 / np.sqrt(3.0)
        return np.array([-g, g]), np.array([1.0, 1.0])
    if n == 3:
        g = np.sqrt(3.0 / 5.0)
        return np.array([-g, 0.0, g]), np.array([5./9., 8./9., 5./9.])
    raise ValueError(n)


# --- Common isoparametric volume integration ---

def _volume_isoparametric(
    all_nodes: np.ndarray,
    gauss_pts: np.ndarray,
    gauss_wts: np.ndarray,
    dshape_func,
) -> np.ndarray:
    """Compute volume via Gauss quadrature for isoparametric elements.

    Parameters
    ----------
    all_nodes : (nelements, nnodes, 3) array of node coordinates
    gauss_pts : (npts, ndim) array of parametric gauss point coordinates
    gauss_wts : (npts,) array of gauss weights
    dshape_func : callable(xi, eta, zeta) -> (nnodes, 3) shape function derivatives
        dN/dxi, dN/deta, dN/dzeta as columns

    Returns
    -------
    volume : (nelements,) array
    """
    nelements = all_nodes.shape[0]
    volume = np.zeros(nelements, dtype=all_nodes.dtype)

    for i_gp in range(len(gauss_wts)):
        xi, eta, zeta = gauss_pts[i_gp]
        w = gauss_wts[i_gp]
        # dN_dxi: (nnodes, 3) - derivatives of shape functions wrt (xi, eta, zeta)
        dN = dshape_func(xi, eta, zeta)

        # Jacobian: J[i,j] = sum_k (dN_k/dxi_j) * x_k_i
        # all_nodes is (nelements, nnodes, 3), dN is (nnodes, 3)
        # J = (nelements, 3, 3): J[:, i, j] = sum_k dN[k, j] * all_nodes[:, k, i]
        J = np.einsum('kj,eki->eji', dN, all_nodes)

        det_J = np.linalg.det(J)
        volume += w * det_J

    return volume


# --- 20-node CHEXA shape function derivatives ---

def _dshape_chexa20(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape function derivatives for 20-node serendipity hexahedron.

    Node ordering follows Nastran convention:
    Corners: 1-8 (standard hex ordering)
    Midside: 9-20 (edges of hex)

    Returns (20, 3) array: dN/dxi, dN/deta, dN/dzeta
    """
    # corner nodes in (xi, eta, zeta) space
    corners = np.array([
        [-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1],
        [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1],
    ], dtype=float)
    # midside nodes: edges 1-2, 2-3, 3-4, 4-1, 5-6, 6-7, 7-8, 8-5, 1-5, 2-6, 3-7, 4-8
    midsides = np.array([
        [0, -1, -1], [+1, 0, -1], [0, +1, -1], [-1, 0, -1],
        [0, -1, +1], [+1, 0, +1], [0, +1, +1], [-1, 0, +1],
        [-1, -1, 0], [+1, -1, 0], [+1, +1, 0], [-1, +1, 0],
    ], dtype=float)

    dN = np.zeros((20, 3))
    # Corner nodes: N_i = (1/8)(1+xi_i*xi)(1+eta_i*eta)(1+zeta_i*zeta)(xi_i*xi+eta_i*eta+zeta_i*zeta-2)
    for i in range(8):
        xi_i, eta_i, zeta_i = corners[i]
        xp = 1 + xi_i * xi
        ep = 1 + eta_i * eta
        zp = 1 + zeta_i * zeta
        s = xi_i * xi + eta_i * eta + zeta_i * zeta - 2

        dN[i, 0] = (1./8.) * xi_i * ep * zp * s + (1./8.) * xp * ep * zp * xi_i
        dN[i, 1] = (1./8.) * eta_i * xp * zp * s + (1./8.) * xp * ep * zp * eta_i
        dN[i, 2] = (1./8.) * zeta_i * xp * ep * s + (1./8.) * xp * ep * zp * zeta_i

    # Midside nodes on xi-direction edges (xi=0): nodes 9,11 (eta=±1, zeta=±1)
    # N = (1/4)(1-xi^2)(1+eta_i*eta)(1+zeta_i*zeta)
    for i in range(12):
        xi_i, eta_i, zeta_i = midsides[i]
        if abs(xi_i) < 0.5:  # xi=0 edge
            ep = 1 + eta_i * eta
            zp = 1 + zeta_i * zeta
            dN[8+i, 0] = (1./4.) * (-2*xi) * ep * zp
            dN[8+i, 1] = (1./4.) * (1 - xi**2) * eta_i * zp
            dN[8+i, 2] = (1./4.) * (1 - xi**2) * ep * zeta_i
        elif abs(eta_i) < 0.5:  # eta=0 edge
            xp = 1 + xi_i * xi
            zp = 1 + zeta_i * zeta
            dN[8+i, 0] = (1./4.) * xi_i * (1 - eta**2) * zp
            dN[8+i, 1] = (1./4.) * xp * (-2*eta) * zp
            dN[8+i, 2] = (1./4.) * xp * (1 - eta**2) * zeta_i
        else:  # zeta=0 edge
            xp = 1 + xi_i * xi
            ep = 1 + eta_i * eta
            dN[8+i, 0] = (1./4.) * xi_i * ep * (1 - zeta**2)
            dN[8+i, 1] = (1./4.) * xp * eta_i * (1 - zeta**2)
            dN[8+i, 2] = (1./4.) * xp * ep * (-2*zeta)

    return dN


def _dshape_ctetra10(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape function derivatives for 10-node quadratic tetrahedron.

    Parametric coords: L1=xi, L2=eta, L3=zeta, L4=1-xi-eta-zeta
    Corners: N1=L1(2L1-1), N2=L2(2L2-1), N3=L3(2L3-1), N4=L4(2L4-1)
    Midside: N5=4*L1*L2, N6=4*L2*L3, N7=4*L1*L3, N8=4*L1*L4, N9=4*L2*L4, N10=4*L3*L4

    Returns (10, 3): dN/dxi, dN/deta, dN/dzeta
    """
    L1 = xi
    L2 = eta
    L3 = zeta
    L4 = 1.0 - xi - eta - zeta

    dN = np.zeros((10, 3))
    # dL1/dxi=1, dL2/deta=1, dL3/dzeta=1, dL4/d*=-1

    # N1 = L1(2L1-1): dN1/dxi = 4L1-1
    dN[0, 0] = 4*L1 - 1;  dN[0, 1] = 0.0;        dN[0, 2] = 0.0
    # N2 = L2(2L2-1): dN2/deta = 4L2-1
    dN[1, 0] = 0.0;       dN[1, 1] = 4*L2 - 1;   dN[1, 2] = 0.0
    # N3 = L3(2L3-1): dN3/dzeta = 4L3-1
    dN[2, 0] = 0.0;       dN[2, 1] = 0.0;         dN[2, 2] = 4*L3 - 1
    # N4 = L4(2L4-1): dN4/d* = -1*(2L4-1) + L4*(-2) = -(4L4-1)
    dN[3, 0] = -(4*L4-1); dN[3, 1] = -(4*L4-1);  dN[3, 2] = -(4*L4-1)

    # N5 = 4*L1*L2: dN5/dxi=4L2, dN5/deta=4L1
    dN[4, 0] = 4*L2;      dN[4, 1] = 4*L1;       dN[4, 2] = 0.0
    # N6 = 4*L2*L3: dN6/deta=4L3, dN6/dzeta=4L2
    dN[5, 0] = 0.0;       dN[5, 1] = 4*L3;       dN[5, 2] = 4*L2
    # N7 = 4*L1*L3: dN7/dxi=4L3, dN7/dzeta=4L1
    dN[6, 0] = 4*L3;      dN[6, 1] = 0.0;        dN[6, 2] = 4*L1
    # N8 = 4*L1*L4: dN8/dxi=4(L4-L1), dN8/deta=-4L1, dN8/dzeta=-4L1
    dN[7, 0] = 4*(L4-L1); dN[7, 1] = -4*L1;      dN[7, 2] = -4*L1
    # N9 = 4*L2*L4: dN9/dxi=-4L2, dN9/deta=4(L4-L2), dN9/dzeta=-4L2
    dN[8, 0] = -4*L2;     dN[8, 1] = 4*(L4-L2);  dN[8, 2] = -4*L2
    # N10 = 4*L3*L4: dN10/dxi=-4L3, dN10/deta=-4L3, dN10/dzeta=4(L4-L3)
    dN[9, 0] = -4*L3;     dN[9, 1] = -4*L3;      dN[9, 2] = 4*(L4-L3)

    return dN


def _dshape_cpenta15(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape function derivatives for 15-node quadratic wedge/pentahedron.

    Parametric coords: L1=xi, L2=eta, L3=1-xi-eta (triangular), zeta in [-1,1].
    Zeta interpolation uses quadratic 1D functions:
      h1 = z(z-1)/2 (bottom), h2 = 1-z^2 (mid), h3 = z(z+1)/2 (top)

    Corners 1-3 at zeta=-1, corners 4-6 at zeta=+1.
    Midside 7-9 on triangle edges at zeta=-1, 10-12 at zeta=+1, 13-15 vertical edges.

    Returns (15, 3): dN/dxi, dN/deta, dN/dzeta
    """
    L1 = xi
    L2 = eta
    L3 = 1.0 - xi - eta
    z = zeta

    # Quadratic 1D interpolation in zeta
    h1 = z * (z - 1) / 2.    # = 1 at z=-1, 0 at z=0, 0 at z=+1
    h3 = z * (z + 1) / 2.    # = 0 at z=-1, 0 at z=0, 1 at z=+1
    h2 = 1.0 - z**2          # = 0 at z=-1, 1 at z=0, 0 at z=+1
    dh1 = z - 0.5            # dh1/dz
    dh3 = z + 0.5            # dh3/dz
    dh2 = -2.0 * z           # dh2/dz

    dN = np.zeros((15, 3))

    # Corner bottom (i=1,2,3): N_i = L_i*(2L_i - 1)*h1
    # N1: L1*(2L1-1)*h1
    dN[0, 0] = (4*L1 - 1)*h1;     dN[0, 1] = 0.0;              dN[0, 2] = L1*(2*L1-1)*dh1
    # N2: L2*(2L2-1)*h1
    dN[1, 0] = 0.0;               dN[1, 1] = (4*L2 - 1)*h1;    dN[1, 2] = L2*(2*L2-1)*dh1
    # N3: L3*(2L3-1)*h1;  dL3/dxi=-1, dL3/deta=-1
    dN[2, 0] = -(4*L3 - 1)*h1;    dN[2, 1] = -(4*L3 - 1)*h1;   dN[2, 2] = L3*(2*L3-1)*dh1

    # Corner top (i=4,5,6): N_i = L_i*(2L_i - 1)*h3
    # N4: L1*(2L1-1)*h3
    dN[3, 0] = (4*L1 - 1)*h3;     dN[3, 1] = 0.0;              dN[3, 2] = L1*(2*L1-1)*dh3
    # N5: L2*(2L2-1)*h3
    dN[4, 0] = 0.0;               dN[4, 1] = (4*L2 - 1)*h3;    dN[4, 2] = L2*(2*L2-1)*dh3
    # N6: L3*(2L3-1)*h3
    dN[5, 0] = -(4*L3 - 1)*h3;    dN[5, 1] = -(4*L3 - 1)*h3;   dN[5, 2] = L3*(2*L3-1)*dh3

    # Midside on bottom triangle edges: N = 4*Li*Lj*h1
    # N7 = 4*L1*L2*h1
    dN[6, 0] = 4*L2*h1;           dN[6, 1] = 4*L1*h1;          dN[6, 2] = 4*L1*L2*dh1
    # N8 = 4*L2*L3*h1
    dN[7, 0] = -4*L2*h1;          dN[7, 1] = 4*(L3-L2)*h1;     dN[7, 2] = 4*L2*L3*dh1
    # N9 = 4*L1*L3*h1
    dN[8, 0] = 4*(L3-L1)*h1;      dN[8, 1] = -4*L1*h1;         dN[8, 2] = 4*L1*L3*dh1

    # Midside on top triangle edges: N = 4*Li*Lj*h3
    # N10 = 4*L1*L2*h3
    dN[9, 0] = 4*L2*h3;           dN[9, 1] = 4*L1*h3;          dN[9, 2] = 4*L1*L2*dh3
    # N11 = 4*L2*L3*h3
    dN[10, 0] = -4*L2*h3;         dN[10, 1] = 4*(L3-L2)*h3;    dN[10, 2] = 4*L2*L3*dh3
    # N12 = 4*L1*L3*h3
    dN[11, 0] = 4*(L3-L1)*h3;     dN[11, 1] = -4*L1*h3;        dN[11, 2] = 4*L1*L3*dh3

    # Midside on vertical edges: N_i = L_i*h2 = L_i*(1-z^2)
    # N13 = L1*h2
    dN[12, 0] = h2;               dN[12, 1] = 0.0;             dN[12, 2] = L1*dh2
    # N14 = L2*h2
    dN[13, 0] = 0.0;              dN[13, 1] = h2;              dN[13, 2] = L2*dh2
    # N15 = L3*h2
    dN[14, 0] = -h2;              dN[14, 1] = -h2;             dN[14, 2] = L3*dh2

    return dN


def _shape_cpyram13(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape functions for 13-node quadratic pyramid.

    Parametric space: xi,eta in [-(1-zeta), (1-zeta)], zeta in [0,1].
    Uses collapsed-coordinate formulation with u=xi/(1-zeta), v=eta/(1-zeta).
    """
    N = np.zeros(13)
    t = zeta
    r = 1.0 - t
    if r < 1e-14:
        N[4] = 1.0
        return N

    u = xi / r
    v = eta / r

    # 5-node linear pyramid functions
    L1 = (1 - u) * (1 - v) * r / 4.
    L2 = (1 + u) * (1 - v) * r / 4.
    L3 = (1 + u) * (1 + v) * r / 4.
    L4 = (1 - u) * (1 + v) * r / 4.
    L5 = t

    # Base midside (serendipity on collapsed quad)
    M6 = (1 - u**2) * (1 - v) * r / 2.
    M7 = (1 + u) * (1 - v**2) * r / 2.
    M8 = (1 - u**2) * (1 + v) * r / 2.
    M9 = (1 - u) * (1 - v**2) * r / 2.

    # Lateral midside
    M10 = 4 * L1 * L5
    M11 = 4 * L2 * L5
    M12 = 4 * L3 * L5
    M13 = 4 * L4 * L5

    # Serendipity corrections to corner and apex nodes
    N[0] = L1 - (M6 + M9) / 2. - M10 / 2.
    N[1] = L2 - (M6 + M7) / 2. - M11 / 2.
    N[2] = L3 - (M7 + M8) / 2. - M12 / 2.
    N[3] = L4 - (M8 + M9) / 2. - M13 / 2.
    N[4] = L5 - (M10 + M11 + M12 + M13) / 2.

    N[5] = M6
    N[6] = M7
    N[7] = M8
    N[8] = M9
    N[9] = M10
    N[10] = M11
    N[11] = M12
    N[12] = M13

    return N


def _dshape_cpyram13(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape function derivatives for 13-node quadratic pyramid.

    Uses central finite differences of _shape_cpyram13 for robustness,
    since the collapsed-coordinate analytical derivatives are complex.

    Returns (13, 3): dN/dxi, dN/deta, dN/dzeta
    """
    h = 1e-7
    dN = np.zeros((13, 3))
    dN[:, 0] = (_shape_cpyram13(xi + h, eta, zeta) - _shape_cpyram13(xi - h, eta, zeta)) / (2*h)
    dN[:, 1] = (_shape_cpyram13(xi, eta + h, zeta) - _shape_cpyram13(xi, eta - h, zeta)) / (2*h)
    dN[:, 2] = (_shape_cpyram13(xi, eta, zeta + h) - _shape_cpyram13(xi, eta, zeta - h)) / (2*h)
    return dN


# --- Shape function evaluators (for mass matrix / field interpolation) ---

def _shape_chexa20(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape functions for 20-node serendipity hexahedron. Returns (20,) array."""
    corners = np.array([
        [-1, -1, -1], [+1, -1, -1], [+1, +1, -1], [-1, +1, -1],
        [-1, -1, +1], [+1, -1, +1], [+1, +1, +1], [-1, +1, +1],
    ], dtype=float)
    midsides = np.array([
        [0, -1, -1], [+1, 0, -1], [0, +1, -1], [-1, 0, -1],
        [0, -1, +1], [+1, 0, +1], [0, +1, +1], [-1, 0, +1],
        [-1, -1, 0], [+1, -1, 0], [+1, +1, 0], [-1, +1, 0],
    ], dtype=float)

    N = np.zeros(20)
    for i in range(8):
        xi_i, eta_i, zeta_i = corners[i]
        xp = 1 + xi_i * xi
        ep = 1 + eta_i * eta
        zp = 1 + zeta_i * zeta
        s = xi_i * xi + eta_i * eta + zeta_i * zeta - 2
        N[i] = (1./8.) * xp * ep * zp * s

    for i in range(12):
        xi_i, eta_i, zeta_i = midsides[i]
        if abs(xi_i) < 0.5:
            N[8+i] = (1./4.) * (1 - xi**2) * (1 + eta_i*eta) * (1 + zeta_i*zeta)
        elif abs(eta_i) < 0.5:
            N[8+i] = (1./4.) * (1 + xi_i*xi) * (1 - eta**2) * (1 + zeta_i*zeta)
        else:
            N[8+i] = (1./4.) * (1 + xi_i*xi) * (1 + eta_i*eta) * (1 - zeta**2)
    return N


def _shape_ctetra10(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape functions for 10-node quadratic tetrahedron. Returns (10,) array."""
    L1 = xi
    L2 = eta
    L3 = zeta
    L4 = 1.0 - xi - eta - zeta
    N = np.zeros(10)
    N[0] = L1 * (2*L1 - 1)
    N[1] = L2 * (2*L2 - 1)
    N[2] = L3 * (2*L3 - 1)
    N[3] = L4 * (2*L4 - 1)
    N[4] = 4*L1*L2
    N[5] = 4*L2*L3
    N[6] = 4*L1*L3
    N[7] = 4*L1*L4
    N[8] = 4*L2*L4
    N[9] = 4*L3*L4
    return N


def _shape_cpenta15(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape functions for 15-node quadratic wedge. Returns (15,) array."""
    L1 = xi
    L2 = eta
    L3 = 1.0 - xi - eta
    z = zeta
    h1 = z * (z - 1) / 2.
    h3 = z * (z + 1) / 2.
    h2 = 1.0 - z**2

    N = np.zeros(15)
    N[0] = L1 * (2*L1 - 1) * h1
    N[1] = L2 * (2*L2 - 1) * h1
    N[2] = L3 * (2*L3 - 1) * h1
    N[3] = L1 * (2*L1 - 1) * h3
    N[4] = L2 * (2*L2 - 1) * h3
    N[5] = L3 * (2*L3 - 1) * h3
    N[6] = 4*L1*L2 * h1
    N[7] = 4*L2*L3 * h1
    N[8] = 4*L1*L3 * h1
    N[9] = 4*L1*L2 * h3
    N[10] = 4*L2*L3 * h3
    N[11] = 4*L1*L3 * h3
    N[12] = L1 * h2
    N[13] = L2 * h2
    N[14] = L3 * h2
    return N


# _shape_cpyram13 is defined above (near _dshape_cpyram13)


# --- Shape functions for linear elements ---

def _shape_ctetra4(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape functions for 4-node linear tetrahedron. Returns (4,) array."""
    return np.array([xi, eta, zeta, 1.0 - xi - eta - zeta])


def _shape_cpenta6(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape functions for 6-node linear wedge. Returns (6,) array."""
    L1 = xi
    L2 = eta
    L3 = 1.0 - xi - eta
    zm = (1 - zeta) / 2.
    zp = (1 + zeta) / 2.
    return np.array([L1*zm, L2*zm, L3*zm, L1*zp, L2*zp, L3*zp])


def _shape_chexa8(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape functions for 8-node linear hexahedron. Returns (8,) array."""
    N = np.zeros(8)
    signs = [(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
             (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]
    for i, (si, ei, zi) in enumerate(signs):
        N[i] = (1 + si*xi) * (1 + ei*eta) * (1 + zi*zeta) / 8.
    return N


def _shape_cpyram5(xi: float, eta: float, zeta: float) -> np.ndarray:
    """Shape functions for 5-node linear pyramid. Returns (5,) array."""
    t = zeta
    r = 1.0 - t
    if r < 1e-14:
        N = np.zeros(5)
        N[4] = 1.0
        return N
    u = xi / r
    v = eta / r
    N = np.zeros(5)
    N[0] = (1 - u) * (1 - v) * r / 4.
    N[1] = (1 + u) * (1 - v) * r / 4.
    N[2] = (1 + u) * (1 + v) * r / 4.
    N[3] = (1 - u) * (1 + v) * r / 4.
    N[4] = t
    return N


# --- Gauss quadrature points for tetrahedra ---

def _tet_gauss_4pt() -> tuple[np.ndarray, np.ndarray]:
    """4-point Gauss quadrature for tetrahedron (degree 2).

    Points in volume coordinates (L1, L2, L3), L4 = 1-L1-L2-L3.
    """
    a = 0.1381966011250105
    b = 0.5854101966249685
    pts = np.array([
        [a, a, a],
        [b, a, a],
        [a, b, a],
        [a, a, b],
    ])
    wts = np.array([1./24., 1./24., 1./24., 1./24.])  # sum = 1/6 (tet reference volume)
    return pts, wts


def _penta_gauss_6pt() -> tuple[np.ndarray, np.ndarray]:
    """6-point Gauss quadrature for wedge (triangular prism).

    3 triangle points x 2 line points.
    Triangle in (xi, eta), zeta in [-1, 1].
    """
    # 3-point triangle rule (degree 2)
    tri_pts = np.array([
        [1./6., 1./6.],
        [2./3., 1./6.],
        [1./6., 2./3.],
    ])
    tri_wts = np.array([1./6., 1./6., 1./6.])  # sum = 1/2 (triangle area)

    # 2-point Gauss on [-1, 1]
    g = 1.0 / np.sqrt(3.0)
    line_pts = np.array([-g, g])
    line_wts = np.array([1.0, 1.0])

    pts = np.zeros((6, 3))
    wts = np.zeros(6)
    idx = 0
    for i in range(3):
        for j in range(2):
            pts[idx, 0] = tri_pts[i, 0]
            pts[idx, 1] = tri_pts[i, 1]
            pts[idx, 2] = line_pts[j]
            wts[idx] = tri_wts[i] * line_wts[j]
            idx += 1
    return pts, wts


def _pyram_gauss_8pt() -> tuple[np.ndarray, np.ndarray]:
    """8-point Gauss quadrature for pyramid.

    Uses 2x2x2 Gauss-Legendre mapped to pyramid parametric space [0,1] in zeta.
    """
    g = 1.0 / np.sqrt(3.0)
    pts_1d = np.array([-g, g])
    wts_1d = np.array([1.0, 1.0])

    # Map zeta from [-1,1] to [0,1]: zeta = (1+z)/2, dz = 2*dzeta
    pts = []
    wts = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                z_ref = (1 + pts_1d[k]) / 2.  # zeta in [0,1]
                pts.append([pts_1d[i], pts_1d[j], z_ref])
                # weight includes Jacobian of mapping from cube to pyramid: (1-zeta)^2
                # and the 1/2 from the [-1,1]->[0,1] mapping
                wts.append(wts_1d[i] * wts_1d[j] * wts_1d[k] * (1 - z_ref)**2 * 0.5)
    return np.array(pts), np.array(wts)


# --- Public higher-order volume functions ---

def volume_ctetra10(all_nodes: np.ndarray) -> np.ndarray:
    """Volume of 10-node quadratic tetrahedra via 4-point Gauss quadrature.

    Parameters
    ----------
    all_nodes : (nelements, 10, 3)
    """
    pts, wts = _tet_gauss_4pt()
    # Negate to match volume_ctetra sign convention (parametric L1,L2,L3 form
    # a left-handed system relative to standard Nastran node ordering)
    return -_volume_isoparametric(all_nodes, pts, wts, _dshape_ctetra10)


def volume_cpenta15(all_nodes: np.ndarray) -> np.ndarray:
    """Volume of 15-node quadratic wedge via 6-point Gauss quadrature.

    Parameters
    ----------
    all_nodes : (nelements, 15, 3)
    """
    pts, wts = _penta_gauss_6pt()
    return _volume_isoparametric(all_nodes, pts, wts, _dshape_cpenta15)


def volume_chexa20(all_nodes: np.ndarray) -> np.ndarray:
    """Volume of 20-node serendipity hexahedra via 2x2x2 Gauss quadrature.

    Parameters
    ----------
    all_nodes : (nelements, 20, 3)
    """
    pts_1d, wts_1d = _gauss_pts_1d(2)
    npts = 8
    pts = np.zeros((npts, 3))
    wts = np.zeros(npts)
    idx = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                pts[idx] = [pts_1d[i], pts_1d[j], pts_1d[k]]
                wts[idx] = wts_1d[i] * wts_1d[j] * wts_1d[k]
                idx += 1
    return _volume_isoparametric(all_nodes, pts, wts, _dshape_chexa20)


def volume_cpyram13(all_nodes: np.ndarray) -> np.ndarray:
    """Volume of 13-node quadratic pyramid via 8-point Gauss quadrature.

    Parameters
    ----------
    all_nodes : (nelements, 13, 3)
    """
    pts, wts = _pyram_gauss_8pt()
    return _volume_isoparametric(all_nodes, pts, wts, _dshape_cpyram13)


# --- Linear element volume functions ---

def volume_ctetra(n1: np.ndarray, n2: np.ndarray,
                  n3: np.ndarray, n4: np.ndarray) -> np.ndarray:
    """is this signed volume?"""
    #volume = -np.dot(n1 - n4, np.cross(n2 - n4, n3 - n4)) / 6.
    n14 = n1 - n4
    n24 = n2 - n4
    n34 = n3 - n4
    n234 = np.cross(n24, n34, axis=1)

    # dot product
    volume = np.einsum("ij, ij->i", n14, n234) / -6.
    return volume

def volume_cpyram(n1: np.ndarray, n2: np.ndarray,
                  n3: np.ndarray, n4: np.ndarray,
                  n5: np.ndarray) -> np.ndarray:
    vpyramid1 = volume_ctetra(n1, n2, n3, n5)
    vpyramid2 = volume_ctetra(n1, n3, n4, n5)
    volume = vpyramid1 + vpyramid2
    return volume

def volume_cpenta(n1: np.ndarray, n2: np.ndarray,
                  n3: np.ndarray, n4: np.ndarray,
                  n5: np.ndarray, n6: np.ndarray) -> np.ndarray:
    nelement = n1.shape[0]
    length1 = np.linalg.norm(n1 - n4, axis=1)
    length2 = np.linalg.norm(n2 - n5, axis=1)
    length3 = np.linalg.norm(n3 - n6, axis=1)
    assert len(length1) == nelement
    avg_length = (length1 + length2 + length3) / 3

    #c1 = (n1 + n2 + n3) / 3.
    #c2 = (n4 + n5 + n6) / 3.
    #avg_length = np.linalg.norm(c1-c2, axis=1)
    #assert len(avg_length) == nelements

    a1 = np.cross(n2-n1, n3-n1, axis=1)
    a2 = np.cross(n5-n4, n6-n4, axis=1)
    assert a1.shape == (nelement, 3)

    # need to divide A1 and A2 by 0.5
    # additional 0.5 factor for the average area
    a_avg = 0.25 * (np.linalg.norm(a1, axis=1) + np.linalg.norm(a2, axis=1))
    assert len(a_avg) == nelement
    #vpyramid1 = volume4_array(n1, n2, n3, n4)
    v_triangular_prism = avg_length * a_avg
    assert len(v_triangular_prism) == nelement
    return v_triangular_prism

def volume_chexa(n1: np.ndarray, n2: np.ndarray,
                 n3: np.ndarray, n4: np.ndarray,
                 n5: np.ndarray, n6: np.ndarray,
                 n7: np.ndarray, n8: np.ndarray) -> np.ndarray:
    nelement = n1.shape[0]
    #vpyramid1 = volume4_array(n1, n2, n3, n5)
    #vpyramid2 = volume4_array(n1, n3, n4, n5)
    #vpyramid = vpyramid1 + vpyramid2

    # https://www.osti.gov/servlets/purl/632793/
    #volume = (
        #det3(x7 - x0, x1 - x0, x3 - x5) +
        #det3(x7 - x0, x4 - x0, x5 - x6) +
        #det3(x7 - x0, x2 - x0, x6 - x3)
    #) / 6.
    #  swap points
    # x2 / x3
    # x6 / x7
    def det3(a, b, c):
        stack = np.dstack([a, b, c])
        d3 = np.linalg.det(stack)
        #except FloatingPointError:
            # recasting it to a float64 array gets rid of the underflow
            #if stack.dtype.name == 'float64':
                #stack2 = stack.astype('float64')
                #d3 = np.linalg.det(stack2)
            #else:
                #raise NotImplementedError(stack.dtype.name)
            #abs_stack = np.abs(stack-stack2)
            #print(f'dstack = {abs_stack.max():.3e}')
        return d3

    #volume = (
        #det3(x6 - x0, x1 - x0, x2 - x5) +
        #det3(x6 - x0, x4 - x0, x5 - x7) +
        #det3(x6 - x0, x3 - x0, x7 - x2)
    #) / 6.
    # add 1
    n71 = n7 - n1
    v1 = det3(n71, n2 - n1, n3 - n6)
    v2 = det3(n71, n5 - n1, n6 - n8)
    v3 = det3(n71, n4 - n1, n8 - n3)
    assert len(v1) == nelement
    #volume = (
        #det3(n71, n2 - n1, n3 - n6) +
        #det3(n71, n5 - n1, n6 - n8) +
        #det3(n71, n4 - n1, n8 - n3)
    #) / 6.
    volume = (v1 + v2 + v3) / 6.
    return volume
