import numpy as np
from pyNastran.bdf import read_bdf, BDF

def partition_matrix(Maa):
    str(Maa)
    MLL = MLR = MRR = MRL = np.zeros((3, 3))
    return MLL, MLR, MRR, MRL

def pfaero(model: BDF):
    """performs aero calcs that are independnt of the structural model"""
    str(model)
    # ???
    AJJ = np.zeros((3, 3))
    SKL = np.zeros((3, 3))
    SKJ = np.zeros((3, 3))
    SKJF = np.zeros((3, 3))
    #KAA = np.zeros((3, 3))
    DJK = np.zeros((3, 3))
    WTFACT = np.zeros((3, 3))
    QKJ = np.zeros((3, 3))
    WKK = np.zeros((3, 3))
    D1JK = np.zeros((3, 3))
    D2JK = np.zeros((3, 3))
    DJKB = np.zeros((3, 3))
    K2JK = np.zeros((3, 3))
    DKJB = np.zeros((3, 3))
    #q_inf = 1.0

    # GI - spline matrix
    unused_GI = np.zeros((3, 3))

    # AIC matirx
    AJJ = np.zeros((3, 3))
    AJJT = AJJ.T

    AJJi = np.linalg.inv(AJJ)

    WSKJ = WKK @ SKL
    unused_KKSJ = WSKJ @ AJJi @ DJK
    unused_KX = WSKJ @ AJJi

    SKJ1 = WTFACT @ SKJF
    unused_KDJB = D1JK + 1j * K2JK

    # SKJ - integration matrix
    unused_QKK = SKJ @ AJJi @ (D1JK + 1j * D2JK)

    is_strip_theory_or_mach_box = False
    if not is_strip_theory_or_mach_box:
        unused_KJ = AJJi @ SKJ1
        unused_QKK = QKJ @ DJKB
        unused_KK = SKJ1 @ AJJT @ DKJB
    return SKJ, D1JK, AJJi

def aestat_rs(model: BDF, Kaa, Maa, TR, TRX):
    """static aero"""
    SKJ, D1JK, AJJi = pfaero(model)
    DJX = np.zeros((3, 3))
    KAA = np.zeros((3, 3))
    unused_KAAX = np.zeros((3, 3))
    KSALL = np.zeros((3, 3))
    KALX = np.zeros((3, 3))
    KRXA = np.zeros((3, 3))
    KLXA = np.zeros((3, 3))
    GKA = np.zeros((3, 3))
    GTKL = np.zeros((3, 3))
    QLL = np.zeros((3, 3))
    QKKS = np.zeros((3, 3))
    QKX = np.zeros((3, 3))
    WKK = np.zeros((3, 3))
    SRKT = np.zeros((3, 3))
    KALLA = np.zeros((3, 3))
    fi_q = np.zeros(3)
    u_l = np.zeros(3)
    u_x = np.zeros(3)
    PL = np.zeros(3)
    WGJ = np.zeros(3)
    FAJE = np.zeros(3)
    q_inf = 1.0

    KLL, KLR, unused_KRR, unused_KRL = partition_matrix(Kaa)
    MLL, MLR, MRR, MRL = partition_matrix(Maa)
    KLLi = np.linalg.inv(KLL)
    D = -KLLi @ KLR
    mr = D.T @ MLL @ D + MRL @ D + D.T @ MLR + MRR

    #ZZX = mr @ TR.T @ TRX

    MSLR = MLL @ D + MLR
    unused_MSRR = MRL @ D + D.T @ MSLR + MRR
    unused_M1RR = D.T @ MRL + MRR
    unused_M1RL = D.T @ MLL + MRL
    unused_RFSOP = TR.T @ TRX

    QAA = GKA.T @ QKKS @ GKA
    QAX = GKA.T @ QKX

    KAAA = -q_inf * QAA
    unused_KAAX = -q_inf * QAX
    unused_KSAA1 = KAA + KAAA

    # partition SUPORTs using KSAA1, KALX, KAAA
    #KLLA, KLRA, KRRA, KRLA = partition_matrix(KAAA)
    KALL = KLL - q_inf * QLL

    # TRX - boolean matrix that selects accelerations from the
    #       aerodynamic extra points
    # TR - transforms acceleratiosnf from aero ref point to supported
    #      DOFs form ALX
    KSALLi = np.linalg.inv(KSALL)
    unused_ALX = KSALLi @ KALX

    # WGJ - user input downwash vector
    # FAJE - user input pressure coefficient vector
    # PSA - external load vector
    # PZ - loads for trim calculation
    # IPZ - restrained elastic dimensional intercepts
    # IPZF2 - Unrestrained elastic dimensional intercepts
    # RINT - Rigid, unsplined dimensional intercepts
    # INTZ - Rigid, splined dimensional intercepts
    # HP0 - Perturbation in the support point deformations relative to mean axes due to external loads

    u_k = GTKL.T @ u_l

    # total downwash velocity
    w_j = D1JK @ u_k + DJX @ u_x + WGJ

    # pressure on aero elements
    unused_FFAJ = q_inf * AJJi @ w_j + q_inf * FAJE


    #  statibility
    QKX = WKK @ SKJ @ AJJi @ DJX
    unused_RSTAB = q_inf * SRKT.T @ QKX
    RINT = q_inf * SRKT.T @ (WKK @ SKJ @ AJJi @ w_j + SKJ @ fi_q)
    unused_KSAZX = D.T @ KLXA + KRXA
    unused_INTZ = GKA.T @ RINT

    mri = np.linalg.inv(mr)
    KALLi = np.linalg.inv(KALL)
    HP = mri @ (D.T @ MLL + MRL) @ KALLA @ ((MLL @ D + MLR) @ TR.T @ TRX + KALX)
    HP0 = -mri @ (D.T @ MLL + MRL) @ KALLi @ PL
    return HP, HP0
