from numpy import matrix, cos, sin, dot, cross
from numpy.linalg import solve


def T(alphaR):
    ca = cos(alphaR)
    sa = sin(alphaR)

    t = matrix([[ca, sa, 0.],
                [-sa, ca, 0.],
                [0., 0., 1.]])
    return t


def makeT(v1, v2):
    t = matrix([[ca, sa, 0.],
                [-sa, ca, 0.],
                [0., 0., 1.]])
    return t


def buildGlobalStiffness(elements):
    nElements = len(elements)

    Kg = matrix(zeros((nElements, nElements), 'd'))

    for element in elements:
        nodes = element.nodes
        K = element.getStiffness()

        # put in global stiffness matrix
        for (i, iNode) in enumerate(nodes):  # row
            for (j, jNode) in enumerate(nodes):  # column
                Kg[iNode, jNode] = K[i, j]
            ###
        ###
    ###
    return Kg

# B - strain displacement matrix
# Uo = 1/2*d.T*B.T*E*B
# Ke = integral(B.T*E*B,dV)


class tri(object):
    def __init__(self):
        pass

    def B(self):
        """
        http://www.me.mtu.edu/~mavable/MEEM4405/Plane.pdf
        """
        A = 1 / 2 * A
        m = matrix(([y23, 0, y31, 0, y12, 0],
                    [0, x32, 0, x13, 0, x21],
                    [x32, y23, x13, y31, x21, y21]))

#def applyBoundaryConditons():
#    pass


def doProblem(elements):
    Kg = buildGlobalStiffness()     # global  stiffness (Kg)
    Kr = applyBoundaryConditions()  # reduced stiffness (Kr)

    #F = K*x
    F = getForces()
    x = solve(Kr, F)  # deflections

    (stress, strain) = recoverStressStrain(elements, x)


def frame3d():
    """
    http://www.ce.ufl.edu/~mih/courses/CES4141/Notes%2049%20-%20Sstan%20Three%20Dimensional%20Example.pdf
    """
    p1 = array([0., 0., 0.]) * 12.
    p2 = array([10., 0., 0.]) * 12.
    p3 = array([10., 10., 0.]) * 12.
    p4 = array([0., 10., 0.]) * 12.

    p5 = array([0., 0., 10.]) * 12.
    p6 = array([10., 0., 10.]) * 12.
    p7 = array([10., 10., 10.]) * 12.
    p8 = array([0., 10., 10.]) * 12.

    Ac = 4.  # in^2
    Jc = 60.  # in4
    Ixc = 650.  # in4
    Iyc = 65.  # in4

    Ab = 3.2
    Jb = 43.
    Ixb = 450.
    Iyb = 32.

    c1 = rod(p5, p1, 5, 1, Ac, Jc, Ixc, Iyc)
    c2 = rod(p6, p2, 6, 2, Ac, Jc, Ixc, Iyc)
    c3 = rod(p7, p3, 7, 3, Ac, Jc, Ixc, Iyc)
    c4 = rod(p8, p4, 8, 4, Ac, Jc, Ixc, Iyc)

    b1 = rod(p1, p2, 1, 2, Ab, Jb, Ixb, Iyb)
    b2 = rod(p2, p3, 2, 3, Ab, Jb, Ixb, Iyb)
    b3 = rod(p3, p4, 3, 4, Ab, Jb, Ixb, Iyb)
    b4 = rod(p4, p1, 4, 1, Ab, Jb, Ixb, Iyb)

    b5 = rod(p5, p6, 5, 6, Ab, Jb, Ixb, Iyb)
    b6 = rod(p6, p7, 6, 7, Ab, Jb, Ixb, Iyb)
    b7 = rod(p7, p8, 7, 8, Ab, Jb, Ixb, Iyb)
    b8 = rod(p8, p5, 8, 5, Ab, Jb, Ixb, Iyb)


class rod(object):
    def __init__(pA, pB, a, b, A, J, Ix, Iy):
        self.pA = pA
        self.pB = pB
        self.a = a
        self.b = b
        self.A = A
        self.J = J
        self.Ix = Ix
        self.Iy = Iy


def Stiffness(self, r, A, E):
    L = norm(r)
    K = A * E / L * matrix([[1., -1.], [-1., 1.]])  # rod
    return K


def truss():
    """
    @code
       ^y
       |
       1
       | \
    L1 |   \
       2-----3 ->x
          L2

    L1 = 4 ft
    L2 = 3 ft
    @endcode
    """
    p1 = array([4., 0., 0.]) * 12.
    p2 = array([0., 0., 0.]) * 12.
    p3 = array([0., 3., 0.]) * 12.
    A = 3.  # in^2
    E = 29000  # ksi

    r12 = r1 - r2
    r13 = r1 - r3
    r23 = r2 - r3


def fKx(K, x):
    f = solve(K, x)
    return f


def beamStiffness(r, A, E, I):
    Ke = matrix(zeros((6, 6), 'd'))
    AE = A * E
    EI = E * I

    if 1:
        Ke[0, 0] = AE / L
        Ke[3, 0] = -AE / L
        Ke[0, 3] = -AE / L
        Ke[3, 3] = AE / L

        Ke[1, 1] = 12 * EI / L ** 3
        Ke[1, 2] = 6 * EI / L ** 2
        Ke[2, 1] = Ke[1, 2]  # 6*EI/L**2
        Ke[2, 2] = 4 * EI / L

        Ke[1, 4] = -Ke[1, 1]  # -12*EI/L**3
        Ke[1, 5] = Ke[1, 2]  # 6*EI/L**2
        Ke[2, 4] = -Ke[1, 2]  # -6*EI/L**2
        Ke[2, 5] = 2 * EI / L

        Ke[4, 1] = -Ke[1, 4]  # -12*EI/L**3
        Ke[4, 2] = Ke[2, 4]  # -6*EI/L**2
        Ke[5, 1] = Ke[1, 2]  # 6*EI/L**2
        Ke[5, 2] = Ke[2, 5]  # 2*EI/L

        Ke[4, 4] = Ke[1, 1]  # 12*EI/L**3
        Ke[4, 5] = Ke[2, 4]  # -6*EI/L**2
        Ke[5, 4] = Ke[2, 4]  # -6*EI/L**2
        Ke[5, 5] = Ke[2, 2]  # 4*EI/L
    else:
        Ke[0, 0] = AE
        Ke[3, 0] = -AE
        Ke[0, 3] = -AE
        Ke[3, 3] = AE

        Ke[1, 1] = 12 * EI / L ** 2
        Ke[1, 2] = 6 * EI / L
        Ke[2, 1] = Ke[1, 2]  # 6*EI/L**2
        Ke[2, 2] = 4 * EI

        Ke[1, 4] = -Ke[1, 1]  # -12*EI/L**3
        Ke[1, 5] = Ke[1, 2]  # 6*EI/L**2
        Ke[2, 4] = -Ke[1, 2]  # -6*EI/L**2
        Ke[2, 5] = 2 * EI

        Ke[4, 1] = -Ke[1, 4]  # -12*EI/L**3
        Ke[4, 2] = Ke[2, 4]  # -6*EI/L**2
        Ke[5, 1] = Ke[1, 2]  # 6*EI/L**2
        Ke[5, 2] = Ke[2, 5]  # 2*EI/L

        Ke[4, 4] = Ke[1, 1]  # 12*EI/L**3
        Ke[4, 5] = Ke[2, 4]  # -6*EI/L**2
        Ke[5, 4] = Ke[2, 4]  # -6*EI/L**2
        Ke[5, 5] = Ke[2, 2]  # 4*EI/L

        Ke = Ke / L
    ###
    return Ke
