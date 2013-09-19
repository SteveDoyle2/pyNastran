from sympy import Symbol, solve
#from sympy import Matrix

a = Symbol('a')
b = Symbol('b')
c = Symbol('c')


class FEM(object):
    def solveABC(self, N, k, av, bv=0., cv=0.):
        N2 = self.subABC(N, av, bv, cv)
        #N2 = N.subs(a,av).subs(b,bv).subs(c,cv)
        print "N  = ", N
        print "Nb = ", N2
        k = solve(N2, k)
        print "k  = ", k[0]
        print ""
        return k[0]

    def subABC(self, N, av, bv=0, cv=0):
        N2 = N.subs(a, av).subs(b, bv).subs(c, cv)
        return N2


class CTRIA3(FEM):
    def N(self):
        k1 = Symbol('k1')
        k2 = Symbol('k2')
        k3 = Symbol('k3')

        N1 = k1 * a - 1
        N2 = k2 * b - 1
        N3 = k3 * c - 1
        # N1+N2+N3 = 1
        #N3 = 1-a-b

        k1 = self.solveABC(N1, k1, 1, 0, 0)
        k2 = self.solveABC(N2, k2, 0, 1, 0)
        k3 = self.solveABC(N3, k3, 0, 0, 1)

        # corners
        N1 = k1 * a
        N2 = k2 * b

        #N3 = c
        # N1+N2+N3 = 1
        #N3 = 1-a-b
        N3 = k3 * (1 - N1 - N2)

        print "CTRIA3"
        print "N1 = %s" % (N1)
        print "N2 = %s" % (N2)
        print "N3 = %s" % (N3)


class CTRIA4(FEM):
    def N(self):
        """
            3
         x     x

        1   4    2

        3-2 -> a=0
        4      a=1/2
        1      a=1

        1-3 -> b=0
        4      b=1/2
        2      b=1

        1-4-2 -> c=0
        3        c=1

        p1 = (1,0,0)
        p2 = (0,1,0)
        p3 = (0,0,1)
        p4 = (0.5,0.5,0)
        """
        k1 = Symbol('k1')
        k2 = Symbol('k2')
        k3 = Symbol('k3')
        k4 = Symbol('k4')
        k5 = Symbol('k5')
        k6 = Symbol('k6')

        # corners
        N1 = a + k1 * a * b
        N2 = b + k2 * a * b
        N3 = c + k3 * a * b

        # midside
        #N4 = (1-a-b-c-1)
        L13 = b
        L23 = a
        N4 = k4 * L13 * L23 - 1

        k1 = self.solveABC(N1, k1, 0.5, 0.5, 0)
        k2 = self.solveABC(N2, k2, 0.5, 0.5, 0)
        k3 = self.solveABC(N3, k3, 0.5, 0.5, 1)
        k4 = self.solveABC(N4, k4, 0.5, 0.5, 0)

        # corners
        N1 = a + k1 * a * b
        N2 = b + k2 * a * b
        N3 = c + k3 * a * b

        # midside
        N4 = k4 * L13 * L23

        #p1 = (1,0,0)
        assert self.subABC(N1, 0, 0, 0) == 0, self.subABC(N1, 0, 0, 0)
        assert self.subABC(N1, 1, 0, 0) == 1, self.subABC(N1, 1, 0, 0)
        assert self.subABC(N1, 0, 1, 0) == 0, self.subABC(N1, 0, 1, 0)
        assert self.subABC(N1, 0, 0, 1) == 0, self.subABC(N1, 0, 0, 1)

        #p2 = (0,1,0)
        assert self.subABC(N2, 0, 0, 0) == 0, self.subABC(N2, 0, 0, 0)
        assert self.subABC(N2, 1, 0, 0) == 0, self.subABC(N2, 1, 0, 0)
        assert self.subABC(N2, 0, 1, 0) == 1, self.subABC(N2, 0, 1, 0)
        assert self.subABC(N2, 0, 0, 1) == 0, self.subABC(N2, 0, 0, 1)

        #p3 = (0,0,1)
        assert self.subABC(N3, 0, 0, 0) == 0, self.subABC(N3, 0, 0, 0)
        assert self.subABC(N3, 1, 0, 0) == 0, self.subABC(N3, 1, 0, 0)
        assert self.subABC(N3, 0, 1, 0) == 0, self.subABC(N3, 0, 1, 0)
        assert self.subABC(N3, 0, 0, 1) == 1, self.subABC(N3, 0, 0, 1)

        #p4 = (0.5,0.5,0)
        assert N4.subs(a, 0.5).subs(b, 0.5) == 1, self.subABC(N4, 0.5, 0.5, 0)

        print "CTRIA4"
        print "N1 = %s" % (N1)
        print "N2 = %s" % (N2)
        print "N3 = %s" % (N3)
        print "N4 = %s" % (N4)


class CTRIA5(FEM):
    def N(self):
        """
            3
         x     5

        1   4    2

        3-2-5 -> a=0
        4        a=1/2
        1        a=1

        1-3   -> b=0
        4-5      b=1/2
        2        b=1

        1-4-2 -> c=0
        5        c=1/2
        3        c=1

        p1 = (1,0,0)
        p2 = (0,1,0)
        p3 = (0,0,1)
        p4 = (0.5,0.5,  0)
        p5 = (0,  0.5,0.5)

        """
        k1 = Symbol('k1')
        k2 = Symbol('k2')
        k3 = Symbol('k3')
        k4 = Symbol('k4')
        k5 = Symbol('k5')
        k6 = Symbol('k6')
        k7 = Symbol('k7')
        k8 = Symbol('k8')

        #c = 1-a-b

        # corners
        N1 = a + k1 * a * b + k6 * b * c
        N2 = b + k2 * a * b + k7 * b * c
        N3 = c + k3 * a * b + k8 * b * c

        # midside
        #N4 = (1-a-b-c-1)
        N4 = k4 * a * b - 1  # from CTRIA6 node 4
        N5 = k5 * b * c - 1  # from CTRIA6 node 5

        k1 = self.solveABC(N1, k1, 0.5, 0.5, 0)
        k2 = self.solveABC(N2, k2, 0.5, 0.5, 0)
        k3 = self.solveABC(N3, k3, 0.5, 0.5, 0)

        k4 = self.solveABC(N4, k4, 0.5, 0.5, 0)
        k5 = self.solveABC(N5, k5, 0, 0.5, 0.5)

        # corners
        N1 = a + k1 * a * b + k6 * b * c
        N2 = b + k2 * a * b + k7 * b * c
        N3 = c + k3 * a * b + k8 * b * c

        # midside
        N4 = k4 * a * b
        N5 = k5 * b * c

        #p1 = (1,0,0)
        #p2 = (0,1,0)
        #p3 = (0,0,1)
        #p4 = (0.5,0.5,  0)
        #p5 = (0,  0.5,0.5)

        #p1 = (1,0,0)
        #assert self.subABC(N1,0,0,0) == 0,self.subABC(N1,0,0,0) # ???
        #self.verifyN(N1,[1,0,0,0,0])
        #self.verifyN(N2,[0,1,0,0,0])
        #self.verifyN(N3,[0,0,1,0,0])
        #self.verifyN(N4,[0,0,0,1,0])
        #self.verifyN(N5,[0,0,0,0,1])

        print "CTRIA4"
        print "N1 = %s" % (N1)
        print "N2 = %s" % (N2)
        print "N3 = %s" % (N3)
        print "N4 = %s" % (N4)
        print "N5 = %s" % (N5)

    def verifyN(self, N, V):
        v = [1, 0, 0]
        assert self.subABC(N, 1, 0, 0) == V[0], 'f(%s)=%s = %s' % (v, N, self.subABC(N, 1, 0, 0))
        v = [0, 1, 0]
        assert self.subABC(N, 0, 1, 0) == V[1], 'f(%s)=%s = %s' % (v, N, self.subABC(N, 0, 1, 0))

        v = [0, 0, 1]
        assert self.subABC(N, 0, 0, 1) == V[2], 'f(%s)%s = %s' % (v, N, self.subABC(N, 0, 0, 1))

        v = [.5, .5, 0]
        assert self.subABC(N, 0.5, 0.5, 0) == V[3], 'f(%s)%s = %s' % (v, N, self.subABC(N, 0.5, 0.5, 0))

        v = [0, .5, .5]
        assert self.subABC(N, 0, 0.5, 0.5) == V[4], 'f(%s)%s = %s' % (v, N, self.subABC(N, 0, 0.5, 0.5))


class CTRIA6(FEM):
    def __init__(self):
        pass

    def N(self):
        """
        @code
        N1 = a**2
        N2 = b**2
        N3 = c**3
        N4 = 4*a*b
        N5 = 4*b*c
        N6 = 4*a*c

            3
         6     5

        1   4    2

        3-5-2 -> a=0
        6-4      a=1/2
        1        a=1

        1-6-3 -> b=0
        4-5      b=1/2
        2        b=1

        1-4-2 -> c=0
        5-6      c=1/2
        3        c=1
        @endcode
        """
        k1 = Symbol('k1')
        k2 = Symbol('k2')
        k3 = Symbol('k3')
        k4 = Symbol('k4')
        k5 = Symbol('k5')
        k6 = Symbol('k6')

        print k1
        L12 = c
        L13 = b
        L23 = a
        L45 = b - 1 / 2
        L46 = a - 1 / 2
        L56 = c - 1 / 2

        # corners
        N1 = k1 * L23 * L46 - 1
        N2 = k2 * L13 * L45 - 1
        N3 = k3 * L12 * L56 - 1

        # midside
        N4 = k4 * L13 * L23 - 1
        N5 = k5 * L12 * L13 - 1
        N6 = k6 * L12 * L23 - 1

        print "N1 = %s" % (N1)
        k1 = self.solveABC(N1, k1, 1, 0, 0)
        k2 = self.solveABC(N2, k2, 0, 1, 0)
        k3 = self.solveABC(N3, k3, 0, 0, 1)
        k4 = self.solveABC(N4, k4, 0.5, 0.5, 0)
        k5 = self.solveABC(N5, k5, 0, 0.5, 0.5)
        k6 = self.solveABC(N6, k6, 0.5, 0, 0.5)
        #print "k1 = ",k1

        # corners
        N1 = k1 * L23 * L46
        N2 = k2 * L13 * L45
        N3 = k3 * L12 * L56

        # midside
        N4 = k4 * L13 * L23
        N5 = k5 * L12 * L13
        N6 = k6 * L12 * L23

        print "CTRIA6"
        print "N1 = %s" % (N1)
        print "N2 = %s" % (N2)
        print "N3 = %s" % (N3)
        print "N4 = %s" % (N4)
        print "N5 = %s" % (N5)
        print "N6 = %s" % (N6)

    def solveABC(self, N, k, av, bv=0., cv=0.):
        N2 = N.subs(a, av).subs(b, bv).subs(c, cv)
        print "N  = ", N
        print "Nb = ", N2
        k = solve(N2, k)
        print "k  = ", k  # [0]
        print ""
        return k[0]


class CQUAD4(FEM):
    def N(self):
        """
        4----3    p1=(-1,-1)
        |    |    p2=( 1,-1)
        |    |    p3=( 1, 1)
        1----2    p4=(-1, 1)

        1-4 -> a=-1
        2-3    a=1

        1-2 -> b=-1
        3-4    b=1
        """
        k1 = Symbol('k1')
        k2 = Symbol('k2')
        k3 = Symbol('k3')
        k4 = Symbol('k4')

        L12 = b + 1
        L23 = a - 1
        L34 = b - 1
        L14 = a + 1

        # corners
        N1 = k1 * L34 * L23 - 1
        N2 = k2 * L14 * L34 - 1
        N3 = k3 * L14 * L12 - 1
        N4 = k4 * L12 * L23 - 1

        k1 = self.solveABC(N1, k1, -1, -1)
        k2 = self.solveABC(N2, k2, 1, -1)
        k3 = self.solveABC(N3, k3, 1, 1)
        k4 = self.solveABC(N4, k4, -1, 1)

        # corners
        N1 = k1 * L34 * L23
        N2 = k2 * L14 * L34
        N3 = k3 * L14 * L12
        N4 = k4 * L12 * L23

        print "CQUAD4"
        print "N1 = %s" % (N1)
        print "N2 = %s" % (N2)
        print "N3 = %s" % (N3)
        print "N4 = %s" % (N4)


class CQUAD9(FEM):  # really a CQUAD8
    def N(self):
        """
        4--7--3    p1=(-1,-1)
        |     |    p2=( 1,-1)
        8  9  6    p3=( 1, 1)
        |     |    p4=(-1, 1)
        1--5--2    p5=(

        1-8-4 -> a=-1
        5-9-7    a=0
        2-6-3    a=1

        1-5-2 -> b=-1
        8-9-6    b=0
        3-7-4    b=1
        """
        k1 = Symbol('k1')
        k2 = Symbol('k2')
        k3 = Symbol('k3')
        k4 = Symbol('k4')
        k5 = Symbol('k5')
        k6 = Symbol('k6')
        k7 = Symbol('k7')
        k8 = Symbol('k8')
        k9 = Symbol('k9')

        L12 = b + 1
        L23 = a - 1
        L34 = b - 1
        L14 = a + 1
        L68 = b
        L57 = a

        # corners
        N1 = k1 * L34 * L23 * L68 * L57 - 1
        N2 = k2 * L14 * L34 * L68 * L57 - 1
        N3 = k3 * L14 * L12 * L68 * L57 - 1
        N4 = k4 * L12 * L23 * L68 * L57 - 1

        #midside
        N5 = k5 * L14 * L23 * L34 * L68 - 1
        N6 = k6 * L14 * L57 * L34 * L12 - 1
        N7 = k7 * L14 * L23 * L68 * L12 - 1
        N8 = k8 * L34 * L12 * L57 * L23 - 1

        #center
        N9 = k9 * L34 * L12 * L14 * L23 - 1

        k1 = self.solveABC(N1, k1, -1, -1)
        k2 = self.solveABC(N2, k2, 1, -1)
        k3 = self.solveABC(N3, k3, 1, 1)
        k4 = self.solveABC(N4, k4, -1, 1)

        #1-8-4 -> a=-1
        #5-9-7    a=0
        #2-6-3    a=1

        #1-5-2 -> b=-1
        #8-9-6    b=0
        #3-7-4    b=1

        k5 = self.solveABC(N5, k5, 0, -1)
        k6 = self.solveABC(N6, k6, 1, 0)
        k7 = self.solveABC(N7, k7, 0, 1)
        k8 = self.solveABC(N8, k8, -1, 0)
        k9 = self.solveABC(N9, k9, 0, 0)

        # corners
        N1 = k1 * L34 * L23 * L68 * L57
        N2 = k2 * L14 * L34 * L68 * L57
        N3 = k3 * L14 * L12 * L68 * L57
        N4 = k4 * L12 * L23 * L68 * L57

        #midside
        N5 = k5 * L14 * L23 * L34 * L68
        N6 = k6 * L14 * L57 * L34 * L12
        N7 = k7 * L14 * L23 * L68 * L12
        N8 = k8 * L34 * L12 * L57 * L23

        #center
        N9 = k9 * L34 * L12 * L14 * L23

        print "CQUAD8"
        print "N1 = %s" % (N1)
        print "N2 = %s" % (N2)
        print "N3 = %s" % (N3)
        print "N4 = %s" % (N4)

        print "N5 = %s" % (N5)
        print "N6 = %s" % (N6)
        print "N7 = %s" % (N7)
        print "N8 = %s" % (N8)
        print "N9 = %s" % (N9)


obj = CTRIA3()
#obj = CTRIA6()
#obj.N()
#obj = CQUAD4()
#obj = CQUAD9()

#obj = CTRIA4()
#obj = CTRIA5()
obj.N()
