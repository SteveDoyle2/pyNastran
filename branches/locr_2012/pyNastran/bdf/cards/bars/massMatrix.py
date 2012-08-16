from sympy import Symbol, Matrix, integrate

z = Symbol('z')
n = Symbol('n')
b = Symbol('b')

p = Symbol('p')
A = Symbol('A')
L = Symbol('L')
t = Symbol('t')
v = Symbol('V')


def bar():
    N1 = (1 - z) / 2
    N2 = (1 + z) / 2
    NT = Matrix([N1, N2])
    pdV = p * A * L / 2
    M = makeM(pdV, NT)
    print "Mbar = \n", M


def truss():
    N1 = (1 - z) / 2
    N2 = (1 + z) / 2
    NT = Matrix([[N1, 0, N2, 0],
                 [0, N1, 0, N2]])
    pdV = p * A * L / 2
    M = makeM(pdV, NT)
    print "Mtruss = \n", M


def quad():
    N1 = (1 - z) * (1 - n) / 4
    N2 = (1 + z) * (1 - n) / 4
    N3 = (1 + z) * (1 + n) / 4
    N4 = (1 - z) * (1 + n) / 4

    N = Matrix([[N1, 0, N2, 0, N3, 0, N4, 0],
                [0, N1, 0, N2, 0, N3, 0, N4]])
    NT = N.transpose()
    pdV = p * A * t / 4  # 4 refers to number of nodes??
    Jacobian = Matrix([[0, 0, 0, 0],  # not done
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0]])
    factorI = Jacobian
    M = makeM(pdV, NT, factorI, levels=2)
    print "Mquad = \n", M


def tet4():
    N1 = z
    N2 = n
    N3 = b
    N4 = 1 - z - n - b

    X = Matrix([[1, x, y, z]])
    X2 = Matrix([[1, x1, y1, z1],
                 [1, x2, y2, z2],
                 [1, x3, y3, z3],
                 [1, x4, y4, z4]])
    N = X * X2.inv()
    N1 = N[0, 0]
    N2 = N[0, 1]
    N3 = N[0, 2]
    N4 = N[0, 3]

    N = Matrix([[N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0, 0],
                [0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0],
                [0, 0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4]])
    NT = N.transpose()
    #pdV = p*v
    pdV = 3
    factorI = 1
    M = makeM(pdV, NT, factorI, levels=3)
    print "Mtet = \n", M


def makeM(pdV, NT, factorI=1, levels=1):
    N = NT.transpose()

    #print "N = \n",N

    print "size(NT) = ", NT.shape
    print "size(N) = ", N.shape

    NtN = NT * N
    B = []
    print "NtN = \n", NtN
    print "size(NtN) = ", NtN.shape

    M = integrate(NtN * factorI, z)
    Mp1 = M.subs(z, 1)
    Mm1 = M.subs(z, -1)
    M2 = Mp1 - Mm1
    if levels >= 2:
        M3 = integrate(M2, n)
        M3p1 = M3.subs(n, 1)
        M3m1 = M3.subs(n, -1)
        M4 = M3p1 - M3m1
        M2 = M4
        if levels >= 3:
            print "M4 = ", M4
            M5 = integrate(M4, b)
            M5p1 = M5.subs(b, 1)
            M5m1 = M5.subs(b, -1)
            M6 = M5p1 - M5m1
            M2 = M6
            print "M6 = ", M6
        ###
    ###
    print "pdV = ", pdV
    MM = pdV * M2
    MM.simplify()
    return MM

#bar()
#truss()
tet4()
