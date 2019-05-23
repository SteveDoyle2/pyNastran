#pylint: disable=C0301,W0612,C0111,R0201,C0103,W0613,R0914
import numpy as np

def merge(amatrix, bmatrix):
    return amatrix + bmatrix

def run(Gka, Wkk, Skj, AJJI, DJK,
        K, KALX, KARX, KALL, KALLI, KAXL, Ka_rl, Ka_rr, Ka_lli,
        K2JE, KARL, KARZX, KSALL, KSALLI, KAA, Ka_lr, K2RZX,
        K2JK, K1JE,
        ik, A, B,
        MLL, MLR, mr, MRR, MRL,
        QPKI, Qke, Qkx,
        PHIA,
        WGJ,
        FAJE,
        D, PL, PR, ZXX, TR, TRX, INTZ,
        bref, cref, qbar, machs,
        ARLR, M4RR, MSRR, AMLR, M5RR, K4LX, TMP2, IPZ,
        SRKT, SKJ, SKL, SKJF,
        DJX, D1JE, D2JE, DJ2K, DJKB, D2JK,
        GTKL,
        WTFACT,
        Ajj, Ajx, ALX, AJJT,
        UX, UL, UINTL,

        NDIM, PHAI, wj, f_j_e,
        subcases, only_steady, only_unsteady):
    WKK = Wkk
    GKA = Gka
    AJJ = Ajj
    Ka_ll = KALL
    Djx = DJX
    Ka_rx = KARX
    Ka_lx = KALX

    q = 1.
    S = 1.
    phi_ai = PHAI

    # KAA = structural stiffness matrix
    # MAA = structural mass matrix
    # PA = vector of applied loads

    # eq 1-64
    Qaa = Gka.T @ Wkk @ Skj @ AJJI @ DJK @ GKA

    # eq 1-65
    QAX = GKA.T @ WKK @ SKJ @ AJJI @ DJX


    # eq 1-75
    ZXX = (
        mr @ TR.T @ TRX
        - (D.T @ Ka_ll + Ka_rl) @ Ka_lli @ (MLL @ D + MLR) @ TR.T @ TRX
        - (D.T @ KALX + KARX) @ (D.T @ KALL + KARL) @ KALLI @ KALX
    )
    PZ = D.T @ PL + PR - (D.T @ KALL + KARL) @ KALLI @ PL

    # eq 1-78
    KRZX = ZXX - mr @ TR.T @ TRX

    # eq 1-80
    NDIM = np.ones((6, 6), dtype='float32')
    NDIM[3, 3] /= bref
    NDIM[4, 4] /= cref
    NDIM[5, 5] /= bref

    # eq 1-88
    K2RR = -(D.T @ Ka_ll + Ka_rl) @ ARLR + (D.T @ Ka_lr + Ka_rr)
    MSRR = mr
    KAZL = D.T @ Ka_ll + Ka_rl
    KARZX = KAZL - KAXL * ALX
    IPZ = INTZ - (K.T @ Ka_ll + Ka_rl) @ UINTL

    # eq 1-89
    M5RR = -K2RR @ M4RR + MSRR
    MIRR = -KAZL @ AMLR + M5RR
    KR1ZX = -K2RR @ K4LX + KARZX
    IPZF = K2RR @ TMP2 + IPZ

    # eq 1-90
    IPZF1 = MIRR**(-1) @ IPZF # solve?
    IPZF2 = MSRR @ IPZF1
    KR2ZX = -MIRR**(-1) @ KR1ZX # solve?
    Z1ZX = MSRR @ K2RZX

    # eq 1-91
    aero_coeffs = 1./(q*S) * NDIM @ TR @ Z1ZX
    aero_coeffs0 = 1./(q*S) * NDIM @ TR @ IPZF2


    RSTAB = qbar @ SRKT.T @ Qkx
    AjjI = Ajj**(-1)
    Qkx = Wkk @ Skj @ AjjI @ Djx # solve?

    # eq 1-92
    RINT = qbar @ SRKT.T @ (Wkk @ Skj @ AjjI @ wj + Skj @ (f_j_e / qbar))

    # eq 1-93
    KSAZX = D.T @ (Ka_lx + Ka_rx)

    # Qii = the generalized aerodynamic matrix
    # phi_ai = matrix of i-set normal mode vectors in the physical a-set
    # Gka = spline matrix from Eq 1-22 reduced to the a-set
    # ue = vector of extra point displacements
    # D1JE =
    # D2JE =
    # wj = downwash
    # WTFACT = weighting matrix Wkk from Eq 1-21

    # eq 1-94
    INTZ = Gka.T @ RINT


    # eq 1-103
    #wj = D1JE + i

    # eq 1-104
    Qie = phi_ai.T @ Gka.T @ WTFACT @ Qke

    # eq 1-105
    Qke = WTFACT @ Skj @ AjjI @ (D1JE + ik * K2JE)


    # eq 1-108
    # A is symmetric nhdpts + 2 matrix
    C = A**(-1) @ B

    #--------------------------------------------------------------
    # PFAERO (section 5.3)

    # 1. Read in DMI
    # 2. process aero model geometry (APD)
    # 4. print USET data if requested (TABPRT)
    # 5. form spline matrix (GI)

    # 7.
    if not only_unsteady:
        # 8. form static aero matrices that are only a function of geometry (ADG)
        ADG = 'ADG'

        # 9. loop on the static aeroelastic subcases
        for subcase in subcases:
            # 10. loop on number of machs per subcase.  For trim, this is one..
            #     for divergence, this is nmachs on the DIVERG card
            # 11. determine mach on the current loop (AELOOP)
            for mach in machs:
                # 12. if aero data exists, skip to 18

                # 13. generate aero matrices (AMG)
                AMG = 'AMG'

                # 14. apply weighting factors (if any)
                WSKJ = WKK @ SKL

                # 16.
                AJJI = AJJ ** (-1)
                QKKS = WSKJ @ AJJI @ DJK

                # 17.
                QKX = WSKJ @ AJJI @ DJK
                 # 18. next mach
            # 19. next subcase


    # 20. skip if only_steady
    if not only_steady:
        # 21. do we need to generate PARAML?
        # 22. loop on machs/reduced frequency pairs
        # 23. determine mach and reduced frequency
        # 24. if the required aero data exists, skip to 35
        # 25. create aero matrices (AMG)
        AMG = 'AMG'

        # 26. apply weight factors (if any)
        SKJ1 = WTFACT @ SKJF

        # 27.
        DKJB = K1JK = ik * D2JK

        # 30.
        QKJ = AJJI @ SKJ.T

        # 31.
        QKK = QKJ @ DJKB

        # 32. if there are user input downwash matrices due to extra points
        QKE = QKJ @ (D1JE + ik * D2JE)

        # 33.
        goto = 35

        if goto == 34:
            # 34.
            QKK = SKJ1 @ AJJT @ DKJB

    # 36.
    #return

    #--------------------------------------------------------------
    # AESTATRS

    # 1.
    MSLR = MLL @ D + MLR
    MSRR = MRL @ D + D.T @ MSLR + MRR
    M1RR = D.T @ MRL + MRR
    M1RR = D.T @ MLL + MRL

    # 2.
    RFSOP = TR.T @ TRX

    # 7.
    QAA = GKA.T @ QKKS @ GKA
    QAX = GKA.T @ QKX
    KAAA = -qbar * QAA
    KAAX = -qbar * QAX
    KSAA1 = KAA + KAAA

    # 9.
    KsALLI = KSALL ** (-1)
    ALX = KSALLI @ KALX

    # KRZX = restrained elastic dimensional derivatives
    # Z1ZX = unrestrained elastic dimensional derivatives
    # RSTAB = rigid, splined dimensional derivatives
    # KSAZX = rigid, splined dimensional derivatives
    # ZZX = solution matrix for trim calculations
    # HP = perturbation in support point deformations relative to mean axis due to aerodyanmic extra points

    # 19.
    UK = GTKL.T @ UL

    # 20. total downwash is computed
    WJ = K1JK @ UK + DJX @ UX + WGJ

    # 21. compute pressures on aero elements
    FFAJ = qbar @ AJJI @ WJ + qbar * FAJE

    # 22. compute forces on aero elements
    PAK = qbar * WSKJ @ AJJ @ WJ + qbar @ SKJ @ FAJE

    #--------------------------------------------------------------
    # DIVERGRS

    # 1.
    GKL = GTKL.T

    # 8.
    QLL = GKL.T @ QKK @ GKL

    #--------------------------------------------------------------
    # FLUTTER


    # 19.
    GPKI = GKA @ PHAI

    # 24.
    QKI = QKK @ GPKI
    QII = GPKI @ QKI

    # 26.
    QIE = GPKI @ QKE

    # 27. form QHH

    # 28. form QKH

    #--------------------------------------------------------------
    # MRFREQRS

    # 3. form:
    GPKI = GKA @ PHIA

    # 7. form:
    QKI = QKK @ GPKI
    QKK = GPKI.T @ QKI

    # 9. form:
    QIE = QPKI.T @ QKE

    # 10. form:
    QKH = merge(QII, QIE)

    # 11. form:
    QKH = merge(QKI, QKE)

    # 12. append QHH & QKH onto QHHA and QKHA using SDR1


def thin_plate_spline(C, wS, mesh, aero_points, node_list, D=1.0):
    """spline function"""
    piD16 = np.pi * D * 16.

    nnodes = len(node_list)
    npoints = len(aero_points.keys())
    Cws = np.linalg.inv(C) * wS  # Cws matrix, P matrix

    wa = {}
    i = 0
    for iaero, aero_node in sorted(aero_points.items()):
        xK = np.zeros(nnodes+3, 'd')
        #nodeI = mesh.Node(iNode)

        xa, ya, za = aero_node

        xK[0] = 1.
        xK[1] = xa
        xK[2] = ya

        j = 3
        for jNode in node_list:
            sNode = mesh.Node(jNode)
            (xs, ys, zs) = sNode.get_position()

            Rij2 = (xa-xs)**2. + (ya-ys)**2  # Rij^2
            if Rij2 == 0.:
                xK[j] = 0.
            else:
                Kij = Rij2 * np.log(Rij2) / piD16  # natural log
                xK[j] = Kij
            j += 1

        wai = xK * Cws
        wa[iaero] = wai[0, 0]
        #print("w[%s]=%s" % (iAero, wi[0, 0]))
        i += 1


    #P = solve(C, wS)
    #C*P = wS
    #P = C^-1*wS
    return Cws, wa

#run()
