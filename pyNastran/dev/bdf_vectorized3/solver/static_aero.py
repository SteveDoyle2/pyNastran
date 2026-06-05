"""Static aeroelastic solver stubs (SOL 144 framework).

Provides the matrix algebra flow for static trim analysis following the
NX Nastran DMAP sequence. All aerodynamic matrices are currently stubbed
as zeros — the actual DLM/spline implementations must be connected.

Matrix definitions
------------------
AJJ : (nj, nj) Aerodynamic Influence Coefficient matrix
SKJ : (nk, nj) Integration matrix (aero forces from pressures)
D1JK : (nj, nk) Substantial differentiation matrix (real part)
D2JK : (nj, nk) Substantial differentiation matrix (imag part)
QKK : (nk, nk) Aerodynamic force matrix in k-set
GKA : (nk, na) Spline matrix (structural DOF -> aero boxes)
WKK : (nk, nk) W-matrix (downwash weighting)
WTFACT : (nk, nk) W-factor weighting
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def partition_matrix(Maa: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Partition a-set matrix into L (flexible) and R (rigid) subsets.

    Parameters
    ----------
    Maa : (na, na)
        Matrix in analysis set.

    Returns
    -------
    MLL, MLR, MRR, MRL : partitioned submatrices
    """
    # TODO: implement actual L/R partitioning based on SUPORT cards
    n = Maa.shape[0]
    nh = n // 2
    MLL = Maa[:nh, :nh]
    MLR = Maa[:nh, nh:]
    MRR = Maa[nh:, nh:]
    MRL = Maa[nh:, :nh]
    return MLL, MLR, MRR, MRL


def pfaero(
    model: BDF,
    nj: int = 0,
    nk: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute aero matrices independent of the structural model.

    This builds the aerodynamic force matrix QKK from the AIC matrix (AJJ),
    integration matrix (SKJ), and downwash matrices (D1JK, D2JK).

    Parameters
    ----------
    model : BDF
        The model (used to extract AERO, CAERO, SPLINE cards in future).
    nj : int
        Number of aero DOFs (boxes). If 0, uses a placeholder size of 3.
    nk : int
        Number of structural DOFs in k-set. If 0, uses nj.

    Returns
    -------
    SKJ : (nk, nj) integration matrix
    D1JK : (nj, nk) real downwash derivative matrix
    AJJi : (nj, nj) inverse of AIC matrix
    QKK : (nk, nk) aero force matrix
    """
    if nj == 0:
        nj = 3
    if nk == 0:
        nk = nj

    # TODO: replace stubs with actual DLM computation
    AJJ = np.zeros((nj, nj))
    SKJ = np.zeros((nk, nj))
    D1JK = np.zeros((nj, nk))
    D2JK = np.zeros((nj, nk))
    WKK = np.zeros((nk, nk))
    WTFACT = np.zeros((nk, nk))
    SKJF = np.zeros((nk, nj))
    K2JK = np.zeros((nj, nk))

    # AIC matrix inverse (would come from DLM panel method)
    # Guard against singular stub
    AJJ_reg = AJJ + np.eye(nj) * 1e-12
    AJJi = np.linalg.inv(AJJ_reg)

    # Weighted integration
    WSKJ = WKK @ SKJ

    # Aerodynamic force matrix: QKK = SKJ @ AJJ^{-1} @ (D1JK + i*D2JK)
    # For static (real) case, only D1JK matters
    QKK = SKJ @ AJJi @ (D1JK + 1j * D2JK)

    # Strip theory / mach box alternate path (not used here)
    SKJ1 = WTFACT @ SKJF

    return SKJ, D1JK, AJJi, QKK.real


def aestat_rs(
    model: BDF,
    Kaa: np.ndarray,
    Maa: np.ndarray,
    TR: np.ndarray,
    TRX: np.ndarray,
    q_inf: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Static aeroelastic trim solution (SOL 144 framework).

    Implements the trim equation sequence from the NX Nastran DMAP:
    1. Partition KAA into flexible (L) and rigid (R) sets via SUPORT
    2. Build aero force matrix QAA from splined QKK
    3. Form aeroelastic stiffness KALL = KLL - q*QLL
    4. Solve for trim unknowns

    Parameters
    ----------
    model : BDF
        The model.
    Kaa : (na, na)
        Stiffness in analysis set.
    Maa : (na, na)
        Mass in analysis set.
    TR : (nr, nr)
        Rigid body transformation matrix.
    TRX : (nr, nx)
        Transform from accelerations to aero extra points.
    q_inf : float
        Dynamic pressure.

    Returns
    -------
    HP : np.ndarray
        Perturbation in support point deformations due to inertia relief.
    HP0 : np.ndarray
        Perturbation in support point deformations due to external loads.
    """
    na = Kaa.shape[0]
    nj = na  # placeholder sizing

    SKJ, D1JK, AJJi, QKK = pfaero(model, nj=nj, nk=na)

    # Partition stiffness and mass
    KLL, KLR, KRR, KRL = partition_matrix(Kaa)
    MLL, MLR, MRR, MRL = partition_matrix(Maa)

    nr = KLL.shape[0]
    nx = TRX.shape[1]

    nk = na  # structural k-set size (placeholder = a-set)
    nx = TRX.shape[1]

    # Placeholder spline and aero matrices (sized to partitioned sets)
    DJX = np.zeros((nj, nx))
    GKA = np.zeros((nk, na))
    GTKL = np.zeros((nk, nr))
    QLL = np.zeros((nr, nr))
    QKKS = np.zeros((nk, nk))
    QKX = np.zeros((nk, nx))
    WKK = np.zeros((nk, nk))
    SRKT = np.zeros((nk, na))
    KALLA = np.zeros((nr, nr))
    fi_q = np.zeros(nj)
    u_l = np.zeros(nr)
    u_x = np.zeros(nx)
    PL = np.zeros(nr)
    WGJ = np.zeros(nj)
    FAJE = np.zeros(nj)

    # Static condensation
    KLLi = np.linalg.inv(KLL + np.eye(nr) * 1e-12)
    D = -KLLi @ KLR

    # Rigid body mass
    mr = D.T @ MLL @ D + MRL @ D + D.T @ MLR + MRR

    MSLR = MLL @ D + MLR

    # Aeroelastic stiffness
    QAA = GKA.T @ QKKS @ GKA
    QAX = GKA.T @ QKX

    KAAA = -q_inf * QAA
    KALL = KLL - q_inf * QLL

    # Aeroelastic stiffness partitioned
    KSALL = KALL + np.eye(nr) * 1e-12
    KALX = np.zeros((nr, TRX.shape[1]))
    KRXA = np.zeros((MRR.shape[0], na))
    KLXA = np.zeros((nr, na))

    KSALLi = np.linalg.inv(KSALL)

    # Downwash calculation
    # GTKL: (nk, nr), so GTKL.T: (nr, nk) — but we need u_k in k-set (nk,)
    u_k = GTKL @ u_l  # (nk,)
    w_j = D1JK @ u_k + DJX @ u_x + WGJ  # (nj,)

    # Stability derivatives
    RINT = q_inf * SRKT.T @ (WKK @ SKJ @ AJJi @ w_j + SKJ @ fi_q)  # (na,)

    # Trim solution
    mri = np.linalg.inv(mr + np.eye(mr.shape[0]) * 1e-12)
    KALLi = np.linalg.inv(KALL + np.eye(KALL.shape[0]) * 1e-12)

    HP = mri @ (D.T @ MLL + MRL) @ KALLA @ ((MLL @ D + MLR) @ TR.T @ TRX + KALX)
    HP0 = -mri @ (D.T @ MLL + MRL) @ KALLi @ PL

    return HP, HP0
