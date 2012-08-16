from sympy import Symbol, simplify, Rational, Matrix

## The purpose of this code is to automatically generate
## the moment of inertia equations.  only the locations of
## The nodes of the elements are required, as well as the
## triangular sections that are generated.
## For example, [1,a,b] means a positive triangular section (+1)
## with the 3 nodes at a, b, and <0,0>.
## This is easily understood if the triangle creates an area.
## Some are positive, some are negative, but if you sum them all
## up you get the total area of the section.  This allows you
## to generate the moments of inertia for complex sections.


def bar():
    """
         y
         ^
         |
    1----|----2
    |    |    |
    |    |    | h
    |    o----------> x
    |         |
    |         |
    4---------3
         b

    b = w1
    h = h1
    """
    print "---Bar---"
    #A = b*h
    h = Symbol('h')
    b = Symbol('b')
    half = Rational(1, 2)
    A = b * h * half

    p1 = [-b * half, h * half]
   #p2 = [ b*half, h*half]
   #p3 = [ b*half,-h*half]
   #p4 = [-b*half,-h*half]

    p2 = addx(p1, b)
    p3 = suby(p2, h)
    p4 = subx(p3, b)
    print "p1 = ", p1
    print "p2 = ", p2
    print "p3 = ", p3
    print "p4 = ", p4
    Atri = half * p1[0] * p1[0]
    print "Atri = ", Atri

    sections = [[1, p1, p2], [1, p2, p3], [1, p3, p4], [1, p4, p1]]
    Ixx(sections)


def cross():
    """
           y
    b1     ^      b3
           |
          b2
       a---p1--b
       |   |   |               h1
    k--l   |   c-----d
    |      |         |
    |      o---------p2-> x    h2
    |                |
    j--i       f-----e
       |       |               h3
       |       |
       h-------g
    """
    #given
    #k-l = 0.5 D1
    #k-j = D4
    #h-g = D2
    #b-g = D3
    print "---Cross---"

    D1 = Symbol('D1')
    D2 = Symbol('D2')
    D3 = Symbol('D3')
    D4 = Symbol('D4')
    o2 = Rational(1, 2)
    # bases
    b2 = D2
    b1 = b3 = D1 * o2

    # heights
    h2 = D4
    h1 = h3 = (D3 - D4) * o2

    p1 = [0, D3 * o2]
    p2 = [(D2 + D1) * o2, 0]

    a = subx(p1, D2 * o2)
    b = addx(p1, D2 * o2)
    c = suby(b, h1)
    d = addx(c, b3)
    e = suby(d, h2)
    f = subx(e, b3)
    g = suby(f, h3)
    h = subx(g, b2)
    i = addy(h, h3)
    j = subx(i, b1)
    k = addy(j, h2)
    l = addx(k, b1)
    A = addy(l, h1)
    #print "a = ",a
    #print "b = ",b

    #print "d = ",d
    #print "e = ",e

    #print "g = ",g
    #print "h = ",h

    #print "k = ",k
    #print "j = ",j

    #print "A = ",A
    sections = [[1, a, b], [1, b, c], [1, c, d], [1, d, e],
                [1, e, f], [1, f, g], [1, g, h], [1, h, i],
                [1, i, j], [1, j, k], [1, k, l], [1, l, a]]
    Ixx(sections)


def box():
    """
          y
          ^
          |
      w1  | w2  w3
          |
    a-----p1------b
    |     |       |  h1
    |  e--|----f  |
    |  |  |    |  |
    |  |  |    |  |  h2
    |  |  o-------|--> x
    |  |       |  |
    |  h-------g  |
    |             |  h3
    d-------------c

    # given
    #a-d     = D2
    #(f-g)_y = D3
    #d-c     = D1
    #(g-c)_x = D4
    """
    print "---Box---"

    D1 = Symbol('D1')
    D2 = Symbol('D2')
    D3 = Symbol('D3')
    D4 = Symbol('D4')

    # height
    h1 = h3 = D3
    h2 = D2 - h1 - h3

    # bases
    w1 = w3 = D4
    w2 = D1 - w1 - w3

    o2 = Rational(1, 2)

    # outer
    p1 = [0, D2 * o2]
    a = subx(p1, D1 * o2)
    b = addx(p1, D1 * o2)
    c = suby(b, D2)
    d = subx(c, D1)

    # inner
    e = suby(p1, D3)
    f = addx(e, w2)
    g = suby(f, h2)
    h = subx(g, w2)

    sections = [[1, a, b], [1, b, c], [1, c, d], [1, d, a],
                [-1, e, f], [-1, f, g], [-1, g, h], [-1, h, e]]
    Ixx(sections)


def Box1():
    """
          y
          ^
          |
      w1  | w2  w3
          |
    a-----p1------b
    |     |       |  h1
    |  e--|----f  |
    |  |  |    |  |
    |  |  |    |  |  h2
    |  |  o-------|--> x
    |  |       |  |
    |  h-------g  |
    |             |  h3
    d-------------c

    # given
    #a-b = D1
    #a-d = D2

    #(b-f)_y = D3
    #(g-c)_y = D4
    #(g-c)_x = D5
    #(d-h)_x = D6

    @warning inertia is wrong b/c 0 is assumed to be <0,0> when it should
    be the output of the first inertia calculation
    """
    print "---Box1---"

    D1 = Symbol('D1')
    D2 = Symbol('D2')
    D3 = Symbol('D3')
    D4 = Symbol('D4')
    D5 = Symbol('D5')
    D6 = Symbol('D6')

    # height
    h1 = D3
    h2 = D2 - D3 - D4
    h3 = D4

    # bases
    w1 = D6
    w2 = D1 - D5 - D6
    w3 = D5

    o2 = Rational(1, 2)

    # outer
    p1 = [0, D2 * o2]
    a = subx(p1, D1 * o2)
    b = addx(p1, D1 * o2)
    c = suby(b, D2)
    d = subx(c, D1)

    # inner
    e = suby(p1, D3)
    f = addx(e, w2)
    g = suby(f, h2)
    h = subx(g, w2)

    sections = [[1, a, b], [1, b, c], [1, c, d], [1, d, a],
                [-1, e, f], [-1, f, g], [-1, g, h], [-1, h, e]]
    Ixx(sections)


def T():
    """
           y
           ^
           |
    a---k-p1-l---b
    |      |     |  h1
    |      |     |
    h---g  | d---c
        |  | |
        |  o----------> x
        |    |
        |    |  h2
    j   f----e   i
      w1  w2   w3

    # given
    a-b     = D1
    (b-e)_y = D2
    a-h     = D3
    f-e     = D4
    @todo Area seems off
    @todo Inertias should be relative to the CG
    @warning wrong...
    """
    print "---T---"

    D1 = Symbol('D1')
    D2 = Symbol('D2')
    D3 = Symbol('D3')
    D4 = Symbol('D4')

    o2 = Rational(1, 2)
    w2 = D4
    h1 = D3
    h2 = D2 - h1
    w1 = w3 = (D1 - D4) * o2

    p1 = [0, 0]
    a = subx(p1, D1 * o2)
    b = addx(p1, D1 * o2)
    c = suby(b, h1)
    d = subx(c, w3)
    e = suby(d, h2)
    f = subx(e, w2)
    g = addy(f, h2)
    h = subx(g, w1)
    A = addy(h, h1)

    i = suby(c, h2)
    j = subx(i, D1)
    k = subx(p1, D4 * o2)
    l = addx(p1, D4 * o2)
    #print "a = ",a
    #print "A = ",A

    sections = [[1, a, b], [1, b, c], [1, d, c], [1, d, e], [1, f, e],
                [1, f, g], [1, g, h], [1, h, a],
                [-1, l, d], [-1, d, g], [-1, g, k]]
    Ixx(sections)


def Ixx(sections):
    r"""
    from http://en.wikipedia.org/wiki/Second_moment_of_area
    \f[ J_{xx} = \frac{1}{12} \sum_{i = 1}^{n-1} ( y_i^2 + y_i y_{i+1} + y_{i+1}^2 ) a_i \f]
    \f[ J_{yy} = \frac{1}{12} \sum_{i = 1}^{n-1} ( x_i^2 + x_i x_{i+1} + x_{i+1}^2 ) a_i \f]
    \f[ J_{xy} = \frac{1}{24} \sum_{i = 1}^{n-1} ( x_i y_{i+1} + 2 x_i y_i + 2 x_{i+1} y_{i+1} + x_{i+1} y_i ) a_i \f]
    """
    h = Symbol('h')
    b = Symbol('b')
    Ixx = 0
    Iyy = 0
    Ixy = 0
    half = Rational(1, 2)
    CG = [0, 0]
    Area = 0
    o3 = Rational(1, 3)
    for i, section in enumerate(sections):
        (sign, p1, p2) = section
        #print "p%s = %s" %(i,p1)
        #print "p%s = %s" %(i+1,p2)
        #Atri = half * abs(p1[0]*p2[1])
        AA = Matrix([p1, p2])
        #print "AA = \n",AA
        detAA = AA.det()
        msgAA = '%s' % (detAA)
        #print "msgAA = ",msgAA
        if '-' in msgAA:
            detAA *= Rational(-1, 1)
        Atri = detAA
        #Atri = half*b*h
        #Atri = (p1[0]*p2[1] - p2[0]*p1[1])*half
        ixxi = p1[1] * p1[1] + p1[1] * p2[1] + p2[1] * p2[1]
        iyyi = p1[0] * p1[0] + p1[0] * p2[0] + p2[0] * p2[0]
        ixyi = (p1[0] * p1[1] + 2 * p1[0] * p1[1] + 2 * p2[0] * p2[1] + p2[0] * p1[1])
        Ixxi = ixxi * Atri
        Iyyi = iyyi * Atri
        Ixyi = ixyi * Atri
        #print "Atri = ",Atri
        #print "Ixxi = ",Ixxi
        Ixx += Ixxi
        Iyy += Iyyi
        Ixy += Ixyi
        #print "Ixx = ",simplify(Ixx),'\n'
        CG = [CG[0] + p1[0] * o3 + p2[0] * o3, CG[1] + p1[1] * o3 + p2[1] * o3]
        Area += Atri
    CG = [CG[0] / Area, CG[1] / Area]
    #print "Ixx = ",simplify(Ixx)
    o12 = Rational(1, 12)
    o24 = Rational(1, 24)
    Ixx = Ixx * o12
    Iyy = Iyy * o12
    Ixy = Ixy * o24
    J = Ixx + Iyy
    print "Area = ", simplify(Area * half)
    print "CG   = ", simplify(CG)
    print "Ixx  = ", simplify(Ixx)
    print "Iyy  = ", simplify(Iyy)
    print "Ixy  = ", simplify(Ixy)
    print "J    = ", simplify(J), '\n'


def add(p1, p2, s1=1, s2=1):
    x = p1[0] * s1 + p2[0] * s2
    y = p1[1] * s1 + p2[1] * s2
    return [simplify(x), simplify(y)]


def sub(p1, p2):
    return add(p1, p2, s2=-1)


def subx(p1, p2x):
    return sub(p1, [p2x, 0])


def suby(p1, p2y):
    return sub(p1, [0, p2y])


def addx(p1, p2x):
    return add(p1, [p2x, 0])


def addy(p1, p2y):
    return add(p1, [0, p2y])

#bar()
#cross()
#box()
T()
