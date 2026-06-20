"""Beam element matrices for CBAR/CBEAM: K, M, KD, PG.

Implements the Timoshenko beam (reduces to Euler-Bernoulli when k1=k2=inf).

DOF order (element local):
    [u1x, u1y, u1z, r1x, r1y, r1z, u2x, u2y, u2z, r2x, r2y, r2z]
     0    1    2    3    4    5    6    7    8    9    10   11

References
----------
- Przemieniecki, "Theory of Matrix Structural Analysis", 1968, Ch. 5.
- NX Nastran Theoretical Manual, Sections 3.1 (BAR/BEAM).
- Bathe, "Finite Element Procedures", 2014, Section 5.4.
"""
from itertools import count
import numpy as np
from numpy.polynomial.legendre import leggauss
from pyNastran.f06.errors import FatalError


# Cache 3-point Gauss-Legendre quadrature (used by every element)
_GAUSS3_XI, _GAUSS3_W = leggauss(3)
_GAUSS3_S = (1.0 + _GAUSS3_XI) / 2.0  # mapped to [0, 1]


def timoshenko_stiffness(
    A: float,
    E: float,
    G: float,
    L: float,
    Iy: float,
    Iz: float,
    J: float,
    k1: float,
    k2: float,
    pa: int = 0,
    pb: int = 0,
    s1: float = 0.0,
    s2: float = 0.0,
) -> np.ndarray:
    """12x12 Timoshenko beam stiffness in element local coordinates.

    Uses the flexibility-based method with 3-point Gauss-Legendre quadrature,
    matching NX Nastran KGG at machine precision. Supports shear relief
    coefficients S1/S2 (PBEAM).

    Parameters
    ----------
    A : float
        Cross-sectional area.
    E : float
        Young's modulus.
    G : float
        Shear modulus.
    L : float
        Element length.
    Iy : float
        Second moment of area about local y (bending in xz-plane).
    Iz : float
        Second moment of area about local z (bending in xy-plane).
    J : float
        Torsional constant.
    k1 : float
        Shear correction factor for plane 1 (xy-plane, Iz). Use >= 1e8 for E-B.
    k2 : float
        Shear correction factor for plane 2 (xz-plane, Iy). Use >= 1e8 for E-B.
    pa, pb : int
        Pin flags at ends A and B.
    s1 : float
        Shear relief coefficient for plane 1 (xy-plane).
    s2 : float
        Shear relief coefficient for plane 2 (xz-plane).

    Returns
    -------
    K : np.ndarray, shape (12, 12)
    """
    n_gauss = 3
    s_pts = _GAUSS3_S
    weights = _GAUSS3_W
    dx_dxi = L / 2.0

    # Axial flexibility
    f_ax = 0.0
    for ig in range(n_gauss):
        f_ax += weights[ig] * dx_dxi / (E * A)
    k_ax = 1.0 / f_ax

    # Torsion flexibility
    if G * J == 0.0:
        raise FatalError(f'G*J=0 (G={G}, J={J}); torsional stiffness is undefined')
    f_tx = 0.0
    for ig in range(n_gauss):
        f_tx += weights[ig] * dx_dxi / (G * J)
    k_tx = 1.0 / f_tx

    # Bending flexibility integrals
    f11_y = 0.0
    f22_y = 0.0
    f12_y = 0.0
    f11_z = 0.0
    f22_z = 0.0
    f12_z = 0.0

    for ig in range(n_gauss):
        s = s_pts[ig]
        wt = weights[ig] * dx_dxi
        t = s  # x/L

        phi_A = 1.0 - t
        phi_B = -t

        # Bending compliance
        f11_y += wt * phi_A * phi_A / (E * Iy)
        f22_y += wt * phi_B * phi_B / (E * Iy)
        f12_y += wt * phi_A * phi_B / (E * Iy)

        f11_z += wt * phi_A * phi_A / (E * Iz)
        f22_z += wt * phi_B * phi_B / (E * Iz)
        f12_z += wt * phi_A * phi_B / (E * Iz)

        # Shear compliance with relief: Q = V - (S/L)*M
        # xz-plane (k2, shear relief s2)
        if k2 > 0.0:
            inv_ksGA_y = 1.0 / (k2 * G * A)
            q_A_y = (1.0 - s2 * (1.0 - t)) / L
            q_B_y = (1.0 + s2 * t) / L
            f11_y += wt * q_A_y * q_A_y * inv_ksGA_y
            f22_y += wt * q_B_y * q_B_y * inv_ksGA_y
            f12_y += wt * q_A_y * q_B_y * inv_ksGA_y

        # xy-plane (k1, shear relief s1)
        if k1 > 0.0:
            inv_ksGA_z = 1.0 / (k1 * G * A)
            q_A_z = (1.0 - s1 * (1.0 - t)) / L
            q_B_z = (1.0 + s1 * t) / L
            f11_z += wt * q_A_z * q_A_z * inv_ksGA_z
            f22_z += wt * q_B_z * q_B_z * inv_ksGA_z
            f12_z += wt * q_A_z * q_B_z * inv_ksGA_z

    # Invert flexibility -> moment stiffness
    F_y = np.array([[f11_y, f12_y], [f12_y, f22_y]])
    K_mm_y = np.linalg.inv(F_y)
    F_z = np.array([[f11_z, f12_z], [f12_z, f22_z]])
    K_mm_z = np.linalg.inv(F_z)

    # Transform to DOF stiffness
    T_xz = np.array([[1.0 / L, 1.0 / L], [-1.0, 0.0], [-1.0 / L, -1.0 / L], [0.0, -1.0]])
    T_xy = np.array([[1.0 / L, 1.0 / L], [1.0, 0.0], [-1.0 / L, -1.0 / L], [0.0, 1.0]])

    K_bxz = T_xz @ K_mm_y @ T_xz.T
    K_bxy = T_xy @ K_mm_z @ T_xy.T

    # Assemble 12x12
    K = np.zeros((12, 12), dtype="float64")

    K[0, 0] = k_ax
    K[0, 6] = -k_ax
    K[6, 0] = -k_ax
    K[6, 6] = k_ax

    K[3, 3] = k_tx
    K[3, 9] = -k_tx
    K[9, 3] = -k_tx
    K[9, 9] = k_tx

    bxy = [1, 5, 7, 11]
    bxz = [2, 4, 8, 10]
    for i_local, i1_global, i2_global in zip(count(), bxy, bxz):
        for j_local, j1_global, j2_global in zip(count(), bxy, bxz):
            K[i1_global, j1_global] = K_bxy[i_local, j_local]
            K[i2_global, j2_global] = K_bxz[i_local, j_local]

    # Pin flags
    if pa != 0:
        for ch in str(pa):
            i = int(ch) - 1
            K[i, :] = 0.0
            K[:, i] = 0.0
    if pb != 0:
        for ch in str(pb):
            i = int(ch) + 5
            K[i, :] = 0.0
            K[:, i] = 0.0
    return K


def consistent_mass(
    A: float,
    L: float,
    rho: float,
    Iy: float,
    Iz: float,
    J: float,
    k1: float,
    k2: float,
    nsm: float = 0.0,
) -> np.ndarray:
    """12x12 consistent Timoshenko beam mass matrix in element local coords.

    Includes translational inertia, rotary inertia, and coupling.

    Parameters
    ----------
    A : float
        Cross-sectional area.
    L : float
        Element length.
    rho : float
        Material density.
    Iy, Iz : float
        Second moments of area.
    J : float
        Polar moment (used for rotary inertia about x).
    k1, k2 : float
        Shear correction factors.
    nsm : float
        Non-structural mass per unit length.

    Returns
    -------
    M : np.ndarray, shape (12, 12)
    """
    m_per_L = rho * A + nsm
    m_total = m_per_L * L

    L2 = L * L

    # Shear parameters
    if k1 > 0.0:
        phi_z = 12.0 * Iz / (k1 * A * L2)
    else:
        phi_z = 0.0
    if k2 > 0.0:
        phi_y = 12.0 * Iy / (k2 * A * L2)
    else:
        phi_y = 0.0

    M = np.zeros((12, 12), dtype="float64")

    # Axial: Nastran CBAR COUPMASS uses m/420*[[175,35],[35,175]]
    M[0, 0] = M[6, 6] = 175.0 * m_total / 420.0
    M[0, 6] = M[6, 0] = 35.0 * m_total / 420.0

    # Torsional: linear interpolation of θ -> rho*Ip*L/6*[[2,1],[1,2]]
    Ip = rho * (Iy + Iz)
    Ip_L = Ip * L
    M[3, 3] = M[9, 9] = Ip_L / 3.0
    M[3, 9] = M[9, 3] = Ip_L / 6.0

    # Bending in xy-plane (Fy/Mz DOFs 1,5,7,11) using Iz
    # Przemieniecki Table 5.3 with Timoshenko shear correction
    # DOF order: [uy1, θz1, uy2, θz2] -> indices [1, 5, 7, 11]
    dz = (1.0 + phi_z) ** 2
    rz = Iz / A  # radius of gyration squared

    a1 = (13.0 / 35.0 + 7.0 / 10.0 * phi_z + phi_z**2 / 3.0 + 6.0 / 5.0 * rz / L2) / dz
    b1 = (
        11.0 / 210.0 + 11.0 / 120.0 * phi_z + phi_z**2 / 24.0 + (1.0 / 10.0 - phi_z / 2.0) * rz / L
    ) / dz
    c1 = (9.0 / 70.0 + 3.0 / 10.0 * phi_z + phi_z**2 / 6.0 - 6.0 / 5.0 * rz / L2) / dz
    d1 = (
        13.0 / 420.0 + 3.0 / 40.0 * phi_z + phi_z**2 / 24.0 - (1.0 / 10.0 - phi_z / 2.0) * rz / L
    ) / dz
    e1 = (
        1.0 / 105.0
        + phi_z / 60.0
        + phi_z**2 / 120.0
        + (2.0 / 15.0 + phi_z / 6.0 + phi_z**2 / 3.0) * rz / L2
    ) / dz
    f1 = (
        1.0 / 140.0
        + phi_z / 60.0
        + phi_z**2 / 120.0
        + (1.0 / 30.0 + phi_z / 6.0 - phi_z**2 / 6.0) * rz / L2
    ) / dz

    # Mass submatrix [uy1, θz1, uy2, θz2]:
    #   [a1,    b1*L,   c1,    -d1*L ]
    #   [b1*L,  e1*L², d1*L,  -f1*L²]
    #   [c1,    d1*L,   a1,    -b1*L ]
    #   [-d1*L, -f1*L², -b1*L, e1*L² ]
    M[1, 1] = M[7, 7] = a1 * m_total
    M[1, 7] = M[7, 1] = c1 * m_total
    M[1, 5] = M[5, 1] = b1 * m_total * L
    M[1, 11] = M[11, 1] = -d1 * m_total * L
    M[7, 5] = M[5, 7] = d1 * m_total * L
    M[7, 11] = M[11, 7] = -b1 * m_total * L
    M[5, 5] = M[11, 11] = e1 * m_total * L2
    M[5, 11] = M[11, 5] = -f1 * m_total * L2

    # Bending in xz-plane (Fz/My DOFs 2,4,8,10) using Iy
    # DOF order: [uz1, θy1, uz2, θy2] -> indices [2, 4, 8, 10]
    # Sign convention: positive Fz -> negative θy (right-hand rule)
    # So the coupling terms flip sign relative to the xy-plane
    dy_sq = (1.0 + phi_y) ** 2
    ry = Iy / A

    a2 = (13.0 / 35.0 + 7.0 / 10.0 * phi_y + phi_y**2 / 3.0 + 6.0 / 5.0 * ry / L2) / dy_sq
    b2 = (
        11.0 / 210.0 + 11.0 / 120.0 * phi_y + phi_y**2 / 24.0 + (1.0 / 10.0 - phi_y / 2.0) * ry / L
    ) / dy_sq
    c2 = (9.0 / 70.0 + 3.0 / 10.0 * phi_y + phi_y**2 / 6.0 - 6.0 / 5.0 * ry / L2) / dy_sq
    d2 = (
        13.0 / 420.0 + 3.0 / 40.0 * phi_y + phi_y**2 / 24.0 - (1.0 / 10.0 - phi_y / 2.0) * ry / L
    ) / dy_sq
    e2 = (
        1.0 / 105.0
        + phi_y / 60.0
        + phi_y**2 / 120.0
        + (2.0 / 15.0 + phi_y / 6.0 + phi_y**2 / 3.0) * ry / L2
    ) / dy_sq
    f2 = (
        1.0 / 140.0
        + phi_y / 60.0
        + phi_y**2 / 120.0
        + (1.0 / 30.0 + phi_y / 6.0 - phi_y**2 / 6.0) * ry / L2
    ) / dy_sq

    # Mass submatrix [uz1, θy1, uz2, θy2]:
    #   [a2,    -b2*L,  c2,     d2*L ]
    #   [-b2*L, e2*L², -d2*L, -f2*L²]
    #   [c2,    -d2*L,  a2,     b2*L ]
    #   [d2*L, -f2*L²,  b2*L,  e2*L²]
    M[2, 2] = M[8, 8] = a2 * m_total
    M[2, 8] = M[8, 2] = c2 * m_total
    M[2, 4] = M[4, 2] = -b2 * m_total * L
    M[2, 10] = M[10, 2] = d2 * m_total * L
    M[8, 4] = M[4, 8] = -d2 * m_total * L
    M[8, 10] = M[10, 8] = b2 * m_total * L
    M[4, 4] = M[10, 10] = e2 * m_total * L2
    M[4, 10] = M[10, 4] = -f2 * m_total * L2

    return M


def geometric_stiffness(
    A: float,
    E: float,
    G: float,
    L: float,
    Iy: float,
    Iz: float,
    J: float,
    k1: float,
    k2: float,
    P: float,
) -> np.ndarray:
    """12x12 geometric (differential) stiffness for a beam under axial load P.

    Uses the EB formula P/(30L)*[36, 3L, ...] which matches NX Nastran KDGG
    at machine precision. NX computes forces from Timoshenko Ke but uses EB Kg.

    Parameters
    ----------
    A : float
        Cross-sectional area.
    E : float
        Young's modulus.
    G : float
        Shear modulus.
    L : float
        Element length.
    Iy, Iz : float
        Second moments of area.
    J : float
        Torsional constant (unused but kept for interface consistency).
    k1, k2 : float
        Shear correction factors (unused — EB formula used for Nastran match).
    P : float
        Axial force in the element (positive = tension).

    Returns
    -------
    KD : np.ndarray, shape (12, 12)

    Notes
    -----
    Validated against NX Nastran 2412 (2026-05-27). The EB formula gives
    machine-precision agreement with Nastran KDGG output.
    """
    KD = np.zeros((12, 12), dtype="float64")
    L2 = L * L
    Ip = Iy + Iz

    c = P / (30.0 * L)

    # Bending in xy-plane (DOFs 1,5,7,11)
    KD[1, 1] = 36.0 * c
    KD[1, 5] = 3.0 * L * c
    KD[1, 7] = -36.0 * c
    KD[1, 11] = 3.0 * L * c
    KD[5, 1] = 3.0 * L * c
    KD[5, 5] = 4.0 * L2 * c
    KD[5, 7] = -3.0 * L * c
    KD[5, 11] = -L2 * c
    KD[7, 1] = -36.0 * c
    KD[7, 5] = -3.0 * L * c
    KD[7, 7] = 36.0 * c
    KD[7, 11] = -3.0 * L * c
    KD[11, 1] = 3.0 * L * c
    KD[11, 5] = -L2 * c
    KD[11, 7] = -3.0 * L * c
    KD[11, 11] = 4.0 * L2 * c

    # Bending in xz-plane (DOFs 2,4,8,10)
    KD[2, 2] = 36.0 * c
    KD[2, 4] = -3.0 * L * c
    KD[2, 8] = -36.0 * c
    KD[2, 10] = -3.0 * L * c
    KD[4, 2] = -3.0 * L * c
    KD[4, 4] = 4.0 * L2 * c
    KD[4, 8] = 3.0 * L * c
    KD[4, 10] = -L2 * c
    KD[8, 2] = -36.0 * c
    KD[8, 4] = 3.0 * L * c
    KD[8, 8] = 36.0 * c
    KD[8, 10] = 3.0 * L * c
    KD[10, 2] = -3.0 * L * c
    KD[10, 4] = -L2 * c
    KD[10, 8] = 3.0 * L * c
    KD[10, 10] = 4.0 * L2 * c

    # Torsion geometric stiffness
    c_tor = P * Ip / (A * L)
    KD[3, 3] = c_tor
    KD[3, 9] = -c_tor
    KD[9, 3] = -c_tor
    KD[9, 9] = c_tor

    return KD


def thermal_load_beam(
    A: float,
    E: float,
    alpha: float,
    L: float,
    ihat: np.ndarray,
    jhat: np.ndarray,
    khat: np.ndarray,
    dT: float,
) -> np.ndarray:
    """12-component thermal load vector for a uniform temperature beam.

    Parameters
    ----------
    A : float
        Cross-sectional area.
    E : float
        Young's modulus.
    alpha : float
        Coefficient of thermal expansion.
    L : float
        Element length.
    ihat, jhat, khat : np.ndarray, shape (3,)
        Element local axes.
    dT : float
        Temperature change (uniform).

    Returns
    -------
    PG : np.ndarray, shape (12,)
        Thermal load vector in basic (global) coordinates.
    """
    # Thermal axial force = E*A*alpha*dT (compressive in element, applied at nodes)
    F_thermal = E * A * alpha * dT

    # In element local coords: equal and opposite axial forces
    Fe = np.zeros(12, dtype="float64")
    Fe[0] = -F_thermal  # node 1 (reaction to thermal expansion)
    Fe[6] = F_thermal  # node 2

    # Transform to basic
    T = np.vstack([ihat, jhat, khat])
    z = np.zeros((3, 3), dtype="float64")
    Teb = np.block([
        [T, z, z, z],
        [z, T, z, z],
        [z, z, T, z],
        [z, z, z, T],
    ])
    PG = Teb.T @ Fe
    return PG


def beam_transforms(ihat, jhat, khat):
    neid = len(ihat)
    Tebs = np.full((neid, 12, 12), dtype="float64")
    for i in range(4):
        r = i * 3
        Teb[:, r, r] = ihat[:, 0]
        Teb[:, r, r + 1] = ihat[:, 1]
        Teb[:, r, r + 2] = ihat[:, 2]
        Teb[:, r + 1, r] = jhat[:, 0]
        Teb[:, r + 1, r + 1] = jhat[:, 1]
        Teb[:, r + 1, r + 2] = jhat[:, 2]
        Teb[:, r + 2, r] = khat[0]
        Teb[:, r + 2, r + 1] = khat[:, 1]
        Teb[:, r + 2, r + 2] = khat[:, 2]
    return Tebs

def beam_transform(
    ihat: np.ndarray,
    jhat: np.ndarray,
    khat: np.ndarray,
) -> np.ndarray:
    """Build the 12x12 coordinate transformation from element to basic.

    Parameters
    ----------
    ihat, jhat, khat : np.ndarray, shape (3,)
        Unit vectors of the element coordinate system.

    Returns
    -------
    Teb : np.ndarray, shape (12, 12)
    """
    Teb = np.zeros((12, 12), dtype="float64")
    for i in range(4):
        r = i * 3
        Teb[r, r] = ihat[0]
        Teb[r, r + 1] = ihat[1]
        Teb[r, r + 2] = ihat[2]
        Teb[r + 1, r] = jhat[0]
        Teb[r + 1, r + 1] = jhat[1]
        Teb[r + 1, r + 2] = jhat[2]
        Teb[r + 2, r] = khat[0]
        Teb[r + 2, r + 1] = khat[1]
        Teb[r + 2, r + 2] = khat[2]
    return Teb


def recover_beam_force(
    Ke: np.ndarray,
    Teb: np.ndarray,
    q_basic: np.ndarray,
) -> np.ndarray:
    """Recover element internal forces from global displacements.

    Parameters
    ----------
    Ke : np.ndarray, shape (12, 12)
        Element stiffness in element local coordinates.
    Teb : np.ndarray, shape (12, 12)
        Transformation from element to basic coordinates.
    q_basic : np.ndarray, shape (12,)
        Nodal displacements in basic coordinates.

    Returns
    -------
    Fe : np.ndarray, shape (12,)
        Element forces in element local coordinates.
        [Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2]
    """
    q_element = Teb @ q_basic
    Fe = Ke @ q_element
    return Fe


def beam_stress_at_points(
    Fe: np.ndarray,
    A: float,
    I1: float,
    I2: float,
    J: float,
    c1: float,
    c2: float,
    d1: float,
    d2: float,
    e1: float,
    e2: float,
    f1: float,
    f2: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute beam longitudinal stress at C,D,E,F recovery points.

    Parameters
    ----------
    Fe : np.ndarray, shape (12,)
        Element forces in local coordinates.
    A : float
        Cross-sectional area.
    I1 : float
        Moment of inertia about axis 1 (for bending in plane 1).
    I2 : float
        Moment of inertia about axis 2 (for bending in plane 2).
    J : float
        Torsional constant.
    c1, c2, d1, d2, e1, e2, f1, f2 : float
        Stress recovery point coordinates (y, z offsets).

    Returns
    -------
    stress_a : np.ndarray, shape (4,)
        Longitudinal stress at end A [C, D, E, F].
    stress_b : np.ndarray, shape (4,)
        Longitudinal stress at end B [C, D, E, F].

    Notes
    -----
    sigma = P/A + My*z/Iy - Mz*y/Iz
    The sign convention follows Nastran: point (c1,c2) means y=c1, z=c2.
    sigma = Fx/A + Mz*c1/Iz + My*c2/Iy  (Nastran convention)
    """
    # Safe division
    inv_I1 = 1.0 / I1 if I1 > 0.0 else 0.0
    inv_I2 = 1.0 / I2 if I2 > 0.0 else 0.0

    # force2stress maps [Fx, Fy, Fz, Mx, My, Mz] -> stress at each point
    # sigma = Fx/A + My*z/Iy + Mz*y/Iz  (Nastran: My=Fe[4], Mz=Fe[5])
    # But Nastran convention: c1=y coord, c2=z coord of recovery point
    # sigma = Fx/A + Mz*(c1)/Iz + My*(c2)/Iy
    #       = Fx/A + 0*Fy + 0*Fz + 0*Mx + c2/I1*My + c1/I2*Mz
    # where I1=Iy (bending about y, uses z-coord) and I2=Iz (bending about z, uses y-coord)
    # NX convention: sigma = P/A - M1*c1/I1 + M2*c2/I2
    # where M1=Mz (plane 1), M2=My (plane 2)
    force2stress = np.array([
        [1.0 / A, 0.0, 0.0, 0.0, c2 * inv_I1, -c1 * inv_I2],
        [1.0 / A, 0.0, 0.0, 0.0, d2 * inv_I1, -d1 * inv_I2],
        [1.0 / A, 0.0, 0.0, 0.0, e2 * inv_I1, -e1 * inv_I2],
        [1.0 / A, 0.0, 0.0, 0.0, f2 * inv_I1, -f1 * inv_I2],
    ])

    stress_a = -force2stress @ Fe[:6]
    stress_b = force2stress @ Fe[6:]
    return stress_a, stress_b


# =============================================================================
# Lumped mass matrix
# =============================================================================


def lumped_mass(
    A: float,
    L: float,
    rho: float,
    Iy: float,
    Iz: float,
    J: float,
    nsm: float = 0.0,
) -> np.ndarray:
    """12x12 lumped (diagonal) beam mass matrix in element local coords.

    Parameters
    ----------
    A : float
        Cross-sectional area.
    L : float
        Element length.
    rho : float
        Material density.
    Iy, Iz : float
        Second moments of area.
    J : float
        Polar moment.
    nsm : float
        Non-structural mass per unit length.

    Returns
    -------
    M : np.ndarray, shape (12, 12)
    """
    m_total = (rho * A + nsm) * L
    Ip = Iy + Iz

    M = np.zeros((12, 12), dtype="float64")

    # Translational: half mass at each node in x, y, z
    for dof in [0, 1, 2]:
        M[dof, dof] = m_total / 2.0
        M[6 + dof, 6 + dof] = m_total / 2.0

    # Torsional rotary inertia
    M[3, 3] = rho * Ip * L / 2.0
    M[9, 9] = rho * Ip * L / 2.0

    # Bending rotary inertia (for Timoshenko beams)
    M[4, 4] = rho * Iy * L / 2.0
    M[10, 10] = rho * Iy * L / 2.0
    M[5, 5] = rho * Iz * L / 2.0
    M[11, 11] = rho * Iz * L / 2.0

    return M


# =============================================================================
# PLOAD1: Consistent distributed/point load vector
# =============================================================================


def _timoshenko_shape_functions(t: float, Phi: float) -> np.ndarray:
    """Phi-dependent displacement shape functions on [0, 1].

    DOFs: [w1, theta1, w2, theta2] where theta is the cross-section rotation.
    When Phi=0, reduces to standard cubic Hermite.

    Parameters
    ----------
    t : float
        Parametric coordinate in [0, 1] (= x/L).
    Phi : float
        Shear deformation parameter = 12*E*I / (ks*G*A*L^2).

    Returns
    -------
    N : np.ndarray, shape (4,)
    """
    d = 1.0 + Phi
    N1 = (1.0 - 3.0 * t**2 + 2.0 * t**3 + Phi * (1.0 - t)) / d
    N2 = (t - 2.0 * t**2 + t**3 + Phi * t * (1.0 - t) / 2.0) / d
    N3 = (3.0 * t**2 - 2.0 * t**3 + Phi * t) / d
    N4 = (-t**2 + t**3 + Phi * t * (t - 1.0) / 2.0) / d
    return np.array([N1, N2, N3, N4])


def beam_pg_distributed(
    qx: float,
    qy: float,
    qz: float,
    L: float,
    E: float,
    G: float,
    A: float,
    Iy: float,
    Iz: float,
    k1: float,
    k2: float,
    x_start: float = 0.0,
    x_end: float | None = None,
    qx_end: float | None = None,
    qy_end: float | None = None,
    qz_end: float | None = None,
    n_gauss: int = 5,
) -> np.ndarray:
    """Consistent nodal load vector for distributed load on a beam.

    Uses Timoshenko shape functions with Gauss-Legendre quadrature.
    Supports uniform and linearly varying (trapezoidal) loads.

    Parameters
    ----------
    qx, qy, qz : float
        Distributed load intensity at x_start (force/length) in local coords.
    L : float
        Element length.
    E, G : float
        Elastic and shear moduli.
    A, Iy, Iz : float
        Section properties.
    k1, k2 : float
        Shear correction factors.
    x_start, x_end : float
        Load application range. Default: full span [0, L].
    qx_end, qy_end, qz_end : float or None
        Load intensity at x_end for trapezoidal load. None = uniform.
    n_gauss : int
        Gauss points for integration.

    Returns
    -------
    fe : np.ndarray, shape (12,)
        Equivalent nodal force vector in element local coordinates.
    """
    if x_end is None:
        x_end = L
    if qx_end is None:
        qx_end = qx
    if qy_end is None:
        qy_end = qy
    if qz_end is None:
        qz_end = qz

    # Phi parameters
    Phi_y = 12.0 * E * Iy / (k2 * G * A * L**2) if k2 > 0.0 else 0.0
    Phi_z = 12.0 * E * Iz / (k1 * G * A * L**2) if k1 > 0.0 else 0.0

    # Gauss points on [x_start, x_end]
    load_L = x_end - x_start
    if load_L <= 0.0:
        return np.zeros(12, dtype="float64")

    xi_pts, w_pts = leggauss(n_gauss)
    x_pts = (xi_pts + 1.0) / 2.0 * load_L + x_start
    w_scaled = w_pts * load_L / 2.0

    fe = np.zeros(12, dtype="float64")

    for k in range(n_gauss):
        x = x_pts[k]
        t = x / L
        wk = w_scaled[k]

        frac = (x - x_start) / load_L
        qx_k = qx * (1.0 - frac) + qx_end * frac
        qy_k = qy * (1.0 - frac) + qy_end * frac
        qz_k = qz * (1.0 - frac) + qz_end * frac

        # Axial: linear shape functions
        fe[0] += wk * qx_k * (1.0 - t)
        fe[6] += wk * qx_k * t

        # xy-plane (qy): Timoshenko shape functions with Phi_z
        N_xy = _timoshenko_shape_functions(t, Phi_z)
        fe[1] += wk * qy_k * N_xy[0]
        fe[5] += wk * qy_k * L * N_xy[1]
        fe[7] += wk * qy_k * N_xy[2]
        fe[11] += wk * qy_k * L * N_xy[3]

        # xz-plane (qz): ry = -dw/dx convention
        N_xz = _timoshenko_shape_functions(t, Phi_y)
        fe[2] += wk * qz_k * N_xz[0]
        fe[4] += wk * qz_k * (-L * N_xz[1])
        fe[8] += wk * qz_k * N_xz[2]
        fe[10] += wk * qz_k * (-L * N_xz[3])

    return fe


def beam_pg_point(
    Px: float,
    Py: float,
    Pz: float,
    x_load: float,
    L: float,
    E: float,
    G: float,
    A: float,
    Iy: float,
    Iz: float,
    k1: float,
    k2: float,
) -> np.ndarray:
    """Consistent nodal forces for a point load at position x_load.

    Parameters
    ----------
    Px, Py, Pz : float
        Point load components in local coordinates.
    x_load : float
        Position along the beam (0 = node i, L = node j).
    L : float
        Element length.
    E, G, A, Iy, Iz : float
        Material and section properties.
    k1, k2 : float
        Shear correction factors.

    Returns
    -------
    fe : np.ndarray, shape (12,)
        Equivalent nodal force vector in element local coordinates.
    """
    Phi_y = 12.0 * E * Iy / (k2 * G * A * L**2) if k2 > 0.0 else 0.0
    Phi_z = 12.0 * E * Iz / (k1 * G * A * L**2) if k1 > 0.0 else 0.0

    t = x_load / L
    fe = np.zeros(12, dtype="float64")

    # Axial
    fe[0] = Px * (1.0 - t)
    fe[6] = Px * t

    # xy-plane (Py)
    N_xy = _timoshenko_shape_functions(t, Phi_z)
    fe[1] = Py * N_xy[0]
    fe[5] = Py * L * N_xy[1]
    fe[7] = Py * N_xy[2]
    fe[11] = Py * L * N_xy[3]

    # xz-plane (Pz): ry = -dw/dx convention
    N_xz = _timoshenko_shape_functions(t, Phi_y)
    fe[2] = Pz * N_xz[0]
    fe[4] = Pz * (-L * N_xz[1])
    fe[8] = Pz * N_xz[2]
    fe[10] = Pz * (-L * N_xz[3])

    return fe
