import numpy as np
from numpy import matrix, cos, sin
from numpy.linalg import solve, norm  # type: ignore


def T(alphaR):
    ca = cos(alphaR)
    sa = sin(alphaR)

    t = np.array([[ca, sa, 0.],
                  [-sa, ca, 0.],
                  [0., 0., 1.]])
    return t


def makeT(v1, v2):
    t = np.array([[ca, sa, 0.],
                  [-sa, ca, 0.],
                  [0., 0., 1.]])
    return t


def build_global_stiffness(elements):
    nElements = len(elements)

    Kg = np.zeros((nElements, nElements), 'd')

    for element in elements:
        nodes = element.nodes
        K = element.get_stiffness_matrix()

        # put in global stiffness matrix
        for (i, inode) in enumerate(nodes):  # row
            for (j, jnode) in enumerate(nodes):  # column
                Kg[inode, jnode] = K[i, j]
    return Kg

# B - strain displacement matrix
# Uo = 1/2*d.T*B.T*E*B
# Ke = integral(B.T*E*B,dV)


class Tri:
    def __init__(self):
        pass

    def B(self):
        """
        http://www.me.mtu.edu/~mavable/MEEM4405/Plane.pdf
        """
        A = 1.0
        A = 1 / 2 * A
        m = matrix(([y23, 0, y31, 0, y12, 0],
                    [0, x32, 0, x13, 0, x21],
                    [x32, y23, x13, y31, x21, y21]))

#def applyBoundaryConditions():
#    pass


def do_problem(elements):
    Kg = build_global_stiffness()     # global  stiffness (Kg)
    Kr = apply_boundary_conditions()  # reduced stiffness (Kr)

    #F = K*x
    F = get_forces()
    x = solve(Kr, F)  # deflections

    (stress, strain) = recover_stress_strain(elements, x)


def frame3d():
    """
    http://www.ce.ufl.edu/~mih/courses/CES4141/Notes%2049%20-%20Sstan%20Three%20Dimensional%20Example.pdf
    """
    p1 = np.array([0., 0., 0.]) * 12.
    p2 = np.array([10., 0., 0.]) * 12.
    p3 = np.array([10., 10., 0.]) * 12.
    p4 = np.array([0., 10., 0.]) * 12.

    p5 = np.array([0., 0., 10.]) * 12.
    p6 = np.array([10., 0., 10.]) * 12.
    p7 = np.array([10., 10., 10.]) * 12.
    p8 = np.array([0., 10., 10.]) * 12.

    Ac = 4.  # in^2
    Jc = 60.  # in4
    Ixc = 650.  # in4
    Iyc = 65.  # in4

    Ab = 3.2
    Jb = 43.
    Ixb = 450.
    Iyb = 32.

    c1 = Rod(p5, p1, 5, 1, Ac, Jc, Ixc, Iyc)
    c2 = Rod(p6, p2, 6, 2, Ac, Jc, Ixc, Iyc)
    c3 = Rod(p7, p3, 7, 3, Ac, Jc, Ixc, Iyc)
    c4 = Rod(p8, p4, 8, 4, Ac, Jc, Ixc, Iyc)

    b1 = Rod(p1, p2, 1, 2, Ab, Jb, Ixb, Iyb)
    b2 = Rod(p2, p3, 2, 3, Ab, Jb, Ixb, Iyb)
    b3 = Rod(p3, p4, 3, 4, Ab, Jb, Ixb, Iyb)
    b4 = Rod(p4, p1, 4, 1, Ab, Jb, Ixb, Iyb)

    b5 = Rod(p5, p6, 5, 6, Ab, Jb, Ixb, Iyb)
    b6 = Rod(p6, p7, 6, 7, Ab, Jb, Ixb, Iyb)
    b7 = Rod(p7, p8, 7, 8, Ab, Jb, Ixb, Iyb)
    b8 = Rod(p8, p5, 8, 5, Ab, Jb, Ixb, Iyb)


class Rod:
    def __init__(self, pA, pB, a, b, A, J, Ix, Iy):
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
    K = A * E / L * np.array([[1., -1.],
                              [-1., 1.]])  # rod
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
    p1 = np.array([4., 0., 0.]) * 12.
    p2 = np.array([0., 0., 0.]) * 12.
    p3 = np.array([0., 3., 0.]) * 12.
    A = 3.  # in^2
    E = 29000  # ksi

    r12 = r1 - r2
    r13 = r1 - r3
    r23 = r2 - r3


def fKx(K, x):
    f = solve(K, x)
    return f


def beam_stiffness(r, A, E, I):
    Ke = np.zeros((6, 6), dtype='float64')
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
    return Ke
