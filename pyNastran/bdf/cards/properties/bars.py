# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All bar properties are defined in this file.  This includes:
 *   PBAR
 *   PBARL

All bars are LineProperty objects.
Multi-segment beams are IntegratedLineProperty objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from itertools import count
from numpy import pi, array

from pyNastran.bdf.field_writer_8 import set_blank_if_default, set_default_if_blank
from pyNastran.bdf.cards.baseCard import Property
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string, string_or_blank, blank,
    integer_or_double, double_string_or_blank, fields, integer_double_string_or_blank)
from pyNastran.utils.mathematics import integrate_line, integrate_positive_line
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

def IyyBeam(b, h):
    return 1 / 12. * b * h ** 3


def IBeam(b, h):
    f = 1 / 12. * b * h
    Iyy = f * h * h  # 1/12.*b*h**3
    Izz = f * b * b  # 1/12.*h*b**3
    Iyz = 0.
    return (Iyy, Izz, Iyz)


def IBeamOffset(b, h, y, z):
    A = b * h
    f = 1. / 12. * A

    Iyy = f * h * h  # 1/12.*b*h**3
    Izz = f * b * b  # 1/12.*h*b**3
    Iyz = 0.

    Iyy += A * y * y
    Izz += A * z * z
    Iyz += A * y * z
    return (Iyy, Izz, Iyz)


def getInertiaRectangular(sections):
    """
    Calculates the moment of inertia for a section about the CG.

    :param sections: [[b,h,y,z]_1,...] y,z is the centroid
                     (x in the direction of the beam,
                      y right, z up)
    :returns: interiaParameters list of [Area, Iyy, Izz, Iyz]

    .. seealso:: http://www.webs1.uidaho.edu/mindworks/Machine_Design/Posters/PDF/Moment%20of%20Inertia.pdf
    """
    As = []
    Ax = 0.
    Ay = 0.
    for section in sections:
        (b, h, x, y) = section
        A = b * h
        As.append(A)
        Ax += A * x
        Ay += A * y

    xCG = Ax / A
    yCG = Ay / A
    Axx = 0.
    Ayy = 0.
    Axy = 0.
    for (i, section) in enumerate(sections):
        (b, h, x, y) = section
        #A = b*h
        #As.append(A)
        Axx += As[i] * (x - xCG) ** 2
        Ayy += As[i] * (y - yCG) ** 2
        Axy += As[i] * (x - xCG) * (y - yCG)
    Ixx = Axx / A
    Iyy = Ayy / A
    Ixy = Axy / A
    return (A, Ixx, Iyy, Ixy)


class LineProperty(Property):
    def __init__(self, card, data):
        self.Type = None
        self.dim = None
        self.A = None
        self.i1 = None
        self.i2 = None
        self.i12 = None
        self.j = None
        self.nsm = None
        Property.__init__(self, card, data)

    #def D_bending(self):
        #pass

    #def D_axial(self):
        #pass

    #def D_thermal(self):
        #pass

    #def D_shear(self):
        #pass

    def Rho(self):
        return self.mid.rho

    def Area(self):
        return self.A

    def Nsm(self):
        return self.nsm

    def J(self):
        return self.j

    def I11(self):
        return self.i1

    def I22(self):
        return self.i2

    def E(self):
        return self.mid.E

    def G(self):
        return self.mid.G

    def Nu(self):
        return self.mid.nu

    def CA_Section(self, iFace, iStart, dims):
        """
        ::

          ---msg1---
          H1=0.1
          W1=0.05

          ---msg2---
          Face_1 = geompy.MakeFaceHW(H1, W1, 1)
          geompy.addToStudy( Face_1, 'Face_1' )

          ---msg---
          H1=0.1
          W1=0.05
          Face_1 = geompy.MakeFaceHW(H1, W1, 1)
          geompy.addToStudy( Face_1, 'Face_1' )
        """
        msg1 = ''
        msg2 = 'Face_%s = geompy.MakeFaceHW(' % (iFace + 1)
        for (i, dim) in enumerate(dims):
            msg1 += 'D%s = %s\n' % (iStart + i, dim)
            msg2 += 'D%s,' % (iStart + i)
        msg2 += '1)\n'
        msg2 += "geompy.addToStudy(Face_%i, 'Face_%i')\n" % (iFace, iFace)
        return msg1 + msg2

    def IAreaL(self, dim):
        if self.Type == 'ROD':
            R = dim[0]
            A = pi * R ** 2
            Iyy = A * R ** 2 / 4.
            Izz = Iyy
            Iyz = 0.
        elif self.Type == 'TUBE':
            R1 = dim[0]
            R2 = dim[1]
            A1 = pi * R1 ** 2
            Iyy1 = A1 * R1 ** 2 / 4.
            A2 = pi * R2 ** 2
            Iyy2 = A2 * R2 ** 2 / 4.
            A = A1 - A2
            Iyy = Iyy1 - Iyy2
            Izz = Iyy
            Iyz = 0.

        elif self.Type == 'I':
            # |  ------------
            # |  |    A     | d5
            # |  ------------
            # |     >| |<--d3
            # |      |B|           "I" beam
            # | d1   | |
            # |      | |
            # |   ----------
            # |   |   C    |  d5
            # |   ----------
            sections = []
            h1 = dim[5]  # d2
            w1 = dim[2]
            y1 = dim[0] / 2. - h1
            sections.append([w1, h1, 0., y1])

            h3 = dim[4]
            w3 = dim[1]
            y3 = -dim[0] / 2. + h3
            sections.append([w3, h3, 0., y1])

            h2 = dim[0] - h1 - h3
            w2 = dim[3]  # d1
            sections.append([w2, h2, 0., 0.])

            (A, Iyy, Izz, Iyz) = getInertiaRectangular(sections)
            assert Iyz == 0.

        elif self.Type == 'BAR':
            #: *-------*
            #: |       |
            #: |  BAR  |h1
            #: |       |
            #: *-------*
            #:    w1
            #: I_{xx}=\frac{bh^3}{12}
            #: I_{yy}=\frac{hb^3}{12}
            h1 = dim[1]
            w1 = dim[0]
            A = h1 * w1
            Iyy = 1 / 12. * w1 * h1 ** 3
            Izz = 1 / 12. * h1 * w1 ** 3
            Iyz = 0.  #: .. todo:: is the Ixy of a bar 0 ???

        else:
            msg = 'Type=%s is not supported for %s class...' % (self.Type,
                                                                self.type)
            raise NotImplementedError(msg)
        return (A, Iyy, Izz, Iyz)

    def I1(self):
        I = self.I1_I2_I12()
        return I[0]

    def I2(self):
        I = self.I1_I2_I12()
        return I[1]

    def I12(self):
        try:
            I = self.I1_I2_I12()
        except:
            print(str(self))
            raise
        return I[2]

    def I1_I2_I12(self):
        r"""
        ::

          BAR
              2
              ^
              |
          *---|--*
          |   |  |
          |   |  |
          |h  *-----------> 1
          |      |
          |   b  |
          *------*

        .. math:: I_1 = \frac{1}{12} b h^3

        .. math:: I_2 = \frac{1}{12} h b^3
        """
        dim = self.dim
        if self.Type == 'ROD':
            R = dim[0]
            A = pi * R ** 2
            I1 = A * R ** 2 / 4.
            I2 = I1
            I12 = 0.
        elif self.Type == 'BAR':
            I1 = 1 / 12. * dim[0] * dim[1] ** 3
            I2 = 1 / 12. * dim[1] * dim[0] ** 3
            I12 = 0.
        else:
            msg = 'I1_I2_I12; Type=%s is not supported for %s class...' % (self.Type, self.type)
            raise NotImplementedError(msg)
        return(I1, I2, I12)

def _bar_areaL(class_name, Type, dim):
    """
    Area(x) method for the PBARL and PBEAML classes (pronounced **Area-L**)

    :param self:   the object pointer
    :param dim:    a list of the dimensions associated with **Type**
    :returns Area: Area of the given cross section defined
                   by **self.Type**

    .. note:: internal method
    """
    if Type == 'ROD':
        """
        This is a circle if you couldn't tell...
          __C__
         /     \
        |       |
        F       D
        |       |
         \     /
          \_E_/

        Radius = dim1
        """
        A = pi * dim[0] ** 2
    elif Type == 'TUBE':
        A = pi * (dim[0] ** 2 - dim[1] ** 2)
    elif Type == 'I':
        h1 = dim[5]
        w1 = dim[2]

        h3 = dim[4]
        w3 = dim[1]

        h2 = dim[0] - h1 - h3
        w2 = dim[3]
        A = h1 * w1 + h2 * w2 + h3 * w3
    elif Type == 'L':
        """

         D4
        F--C      ^
        |  |      |
        |  |      |
        |  |      | D2
        |  +---+  |
        |   D3 |  |
        E------D  |

        <------> D1
        """
        (d1, d2, d3, d4) = dim
        A1 = d1 * d3

        h2 = (d2 - d3)
        A2 = h2 * d4
        A = A1 + A2
    elif Type == 'CHAN':
        h1 = dim[3]
        w1 = dim[0]

        h3 = h1
        w3 = w1
        h2 = dim[1] - h1 - h3
        w2 = dim[2]
        A = h1 * w1 + h2 * w2 + h3 * w3
    elif Type == 'T':
        h1 = dim[2]
        w1 = dim[0]

        h2 = dim[1] - h1
        w2 = dim[3]
        A = h1 * w1 + h2 * w2
    elif Type == 'BOX':
        h1 = dim[2]
        w1 = dim[0]

        h2 = dim[1] - 2 * h1
        w2 = dim[3]
        A = 2 * (h1 * w1 + h2 * w2)
    elif Type == 'BAR':
        """
        <------> D1

        F------C  ^
        |      |  |
        |      |  | D2
        |      |  |
        E------D  |
        """
        h1 = dim[1]
        w1 = dim[0]
        A = h1 * w1
    elif Type == 'CROSS':
        h1 = dim[2]
        w1 = dim[1]

        h2 = dim[3]
        w2 = dim[0]
        A = h1 * w1 + h2 * w2
    elif Type == 'H':
        h1 = dim[2]
        w1 = dim[1]

        h2 = dim[3]
        w2 = dim[0]
        A = h1 * w1 + h2 * w2
    elif Type == 'T1':
        h1 = dim[0]
        w1 = dim[2]

        h2 = dim[3]
        w2 = dim[1]
        A = h1 * w1 + h2 * w2
    elif Type == 'I1':
        h2 = dim[2]
        w2 = dim[1]

        h1 = dim[3] - h2
        w1 = dim[0] + w2
        A = h1 * w1 + h2 * w2
    elif Type == 'CHAN1':
        h2 = dim[2]
        w2 = dim[1]

        h1 = dim[3] - h2
        w1 = dim[0] + w2
        A = h1 * w1 + h2 * w2
    elif Type == 'Z':
        h2 = dim[2]
        w2 = dim[1]

        h1 = dim[3] - h2
        w1 = dim[0]
        A = h1 * w1 + h2 * w2
    elif Type == 'CHAN2':
        h2 = dim[1]
        w2 = dim[3]

        h1 = dim[2] - h2
        w1 = dim[0] * 2
        A = h1 * w1 + h2 * w2

    elif Type == 'T2':
        h1 = dim[3]
        w1 = dim[1]

        h2 = h1 - dim[2]
        w2 = dim[0]
        A = h1 * w1 + h2 * w2
    elif Type == 'BOX1':
        h1 = dim[2]  # top
        w1 = dim[0]

        h2 = dim[3]  # btm
        A1 = (h1 + h2) * w1

        h3 = dim[1] - h1 - h2  # left
        w3 = dim[5]

        w4 = dim[4]  # right
        A2 = h3 * (w3 + w4)
        A = A1 + A2
    elif Type == 'HEXA':
        hBox = dim[2]
        wBox = dim[1]

        wTri = dim[0]
        A = hBox * wBox - wTri * hBox
    elif Type == 'HAT':
        w = dim[1]      # constant width (note h is sometimes w)
        h1 = w           # horizontal lower bar
        w1 = dim[3]

        h2 = dim[0] - 2 * w  # vertical bar
        w2 = w

        h3 = w           # half of top bar
        w3 = dim[2] / 2.

        A = 2 * (h1 * w1 + h2 * w2 + h3 * w3)  # symmetrical box
    elif Type == 'HAT1':
        w = dim[3]

        h0 = dim[4]         # btm bar
        w0 = dim[0] / 2.

        h2 = dim[1] - h0 - 2 * w  # vertical bar
        w2 = w

        h3 = w              # top bar
        w3 = dim[2] / 2.

        h1 = w              # upper, horizontal lower bar (see HAT)
        w1 = w0 - w3

        A = 2 * (h0 * w0 + h1 * w1 + h2 * w2 + h3 * w3)

    elif Type == 'DBOX':
        #
        #  |--2------5----
        #  |     |       |
        #  1     3       6
        #  |     |       |
        #  |--4--|---7---|
        #

        #0,1,2,6,11
        #1,2,3,7,12

        hTotal = dim[11]
        wTotal = dim[0]

        h2 = dim[6]
        w2 = dim[3]

        h4 = dim[7]
        w4 = w2

        h1 = hTotal - h2 - h4
        w1 = dim[3]

        h5 = dim[8]
        w5 = wTotal - w2

        h7 = dim[9]
        w7 = w5

        h6 = hTotal - h5 - h7
        w6 = dim[5]

        h3 = (h1 + h6) / 2.
        w3 = dim[4]

        A = (h1 * w1 + h2 * w2 + h3 * w3 + h4 * w4 +
             h5 * w5 + h6 * w6 + h7 * w7)
    else:
        msg = 'areaL; Type=%s is not supported for %s class...' % (Type,
                                                            class_name)
        raise NotImplementedError(msg)
    return A

class IntegratedLineProperty(LineProperty):
    def __init__(self, card, data):
        self.xxb = None
        self.A = None
        self.j = None
        self.i1 = None
        self.i2 = None
        self.i12 = None
        LineProperty.__init__(self, card, data)

    def Area(self):
        A = integrate_positive_line(self.xxb, self.A)
        return A

    def J(self):
        J = integrate_positive_line(self.xxb, self.j)
        return J

    def I11(self):
        i1 = integrate_positive_line(self.xxb, self.i1)
        return i1

    def I22(self):
        i2 = integrate_positive_line(self.xxb, self.i2)
        return i2

    def I12(self):
        i12 = integrate_line(self.xxb, self.i12)
        return i12

    def Nsm(self):
        #print("xxb = ",self.xxb)
        #print("nsm = ",self.nsm)
        nsm = integrate_positive_line(self.xxb, self.nsm)
        return nsm


class PBAR(LineProperty):
    type = 'PBAR'

    def __init__(self, card=None, data=None, comment=''):
        """
        .. todo::
            support solution 600 default
            do a check for mid -> MAT1      for structural
            do a check for mid -> MAT4/MAT5 for thermal

        +------+-----+-----+-----+----+----+----+-----+-----+
        | PBAR | PID | MID | A   | I1 | I2 | J  | NSM |     |
        +------+-----+-----+-----+----+----+----+-----+-----+
        |      | C1  | C2  | D1  | D2 | E1 | E2 | F1  | F2  |
        +------+-----+-----+-----+----+----+----+-----+-----+
        |      | K1  | K2  | I12 |    |    |    |     |     |
        +------+-----+-----+-----+----+----+----+-----+-----+
        """
        LineProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: property ID -> use Pid()
            self.pid = integer(card, 1, 'pid')
            #: material ID -> use Mid()
            self.mid = integer(card, 2, 'mid')
            #: Area -> use Area()
            self.A = double_or_blank(card, 3, 'A', 0.0)
            #: I1 -> use I1()
            self.i1 = double_or_blank(card, 4, 'I1', 0.0)
            #: I2 -> use I2()
            self.i2 = double_or_blank(card, 5, 'I2', 0.0)

            #: Polar Moment of Inertia J -> use J()
            #: default=1/2(I1+I2) for SOL=600, otherwise 0.0
            #: .. todo:: support SOL 600 default
            self.j = double_or_blank(card, 6, 'J', 0.0)
            #: nonstructral mass -> use Nsm()
            self.nsm = double_or_blank(card, 7, 'nsm', 0.0)

            self.C1 = double_or_blank(card, 9, 'C1', 0.0)
            self.C2 = double_or_blank(card, 10, 'C2', 0.0)
            self.D1 = double_or_blank(card, 11, 'D1', 0.0)
            self.D2 = double_or_blank(card, 12, 'D2', 0.0)
            self.E1 = double_or_blank(card, 13, 'E1', 0.0)
            self.E2 = double_or_blank(card, 14, 'E2', 0.0)
            self.F1 = double_or_blank(card, 15, 'F1', 0.0)
            self.F2 = double_or_blank(card, 16, 'F2', 0.0)

            if self.i1 < 0.:
                raise RuntimeError('I1=%r must be greater than or equal to 0.0' % self.i1)
            if self.i2 < 0.:
                raise RuntimeError('I2=%r must be greater than or equal to 0.0' % self.i2)
            if self.j < 0.:
                raise RuntimeError('J=%r must be greater than or equal to 0.0' % self.j)

            #: I12 -> use I12()
            self.i12 = double_or_blank(card, 19, 'I12', 0.0)

            if self.A == 0.0:
                # K1/K2 must be blank
                #: default=infinite; assume 1e8
                self.K1 = blank(card, 17, 'K1')
                #: default=infinite; assume 1e8
                self.K2 = blank(card, 18, 'K2')
            elif self.i12 != 0.0:
                # K1 / K2 are ignored
                self.K1 = None
                self.K2 = None
            else:
                #: default=infinite; assume 1e8
                self.K1 = double_or_blank(card, 17, 'K1', 1e8)
                #: default=infinite; assume 1e8
                self.K2 = double_or_blank(card, 18, 'K2', 1e8)

            assert len(card) <= 20, 'len(PBAR card) = %i' % len(card)
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.A = data[2]
            self.i1 = data[3]
            self.i2 = data[4]
            self.j = data[5]

            self.nsm = data[6]
            #self.fe  = data[7] #: .. todo:: not documented....
            self.C1 = data[8]
            self.C2 = data[9]
            self.D1 = data[10]
            self.D2 = data[11]
            self.E1 = data[12]
            self.E2 = data[13]
            self.F1 = data[14]
            self.F2 = data[15]
            self.K1 = data[16]
            self.K2 = data[17]
            self.i12 = data[18]

    def _verify(self, xref=False):
        pid = self.Pid()
        mid = self.Mid()
        A = self.Area()
        J = self.J()
        #c = self.c
        nsm = self.Nsm()
        mpa = self.MassPerLength()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(A, float), 'pid=%r' % A
        assert isinstance(J, float), 'cid=%r' % J
        #assert isinstance(c, float), 'c=%r' % c
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        assert isinstance(mpa, float), 'mass_per_length=%r' % mpa

    def MassPerLength(self):
        r"""
        Gets the mass per length :math:`\frac{m}{L}` of the CBAR.

        .. math:: \frac{m}{L} = \rho A + nsm
        """
        A = self.Area()
        rho = self.Rho()
        nsm = self.Nsm()
        return rho * A + nsm

    def cross_reference(self, model):
        msg = ' which is required by PBAR mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)

    def Area(self):
        """
        Gets the area :math:`A` of the CBAR.
        """
        return self.A

    #def Nsm(self):
    #    return self.nsm

    #def J(self):
    #    return self.j

    def I11(self):
        return self.i1

    def I22(self):
        return self.i2

    def I12(self):
        return self.i12

    def write_code_aster(self):  # PBAR
        a = self.Area()
        iy = self.I11()
        iz = self.I22()
        j = self.J()
        msg = ''
        msg += "    POUTRE=_F(GROUP_MA='P%s', # PBAR\n" % (self.pid)
        msg += "              SECTION='GENERALE',\n"
        msg += "              CARA=('A','IY','IZ','JX')\n"
        msg += "              VALE=(%g,  %g,  %g,  %g,)\n" % (a, iy, iz, j)
        msg += "              ORIENTATION=(\n"
        msg += "                    CARA=('VECT_Y'),\n"
        msg += "                    VALE=(1.0,0.0,0.0,),),\n"
        return msg

    def raw_fields(self):
        list_fields = ['PBAR', self.pid, self.Mid(), self.A, self.i1, self.i2,
                       self.j, self.nsm, None, self.C1, self.C2, self.D1, self.D2,
                       self.E1, self.E2, self.F1, self.F2, self.K1, self.K2,
                       self.i12]
        return list_fields

    def repr_fields(self):
        #A  = set_blank_if_default(self.A,0.0)
        i1 = set_blank_if_default(self.i1, 0.0)
        i2 = set_blank_if_default(self.i2, 0.0)
        i12 = set_blank_if_default(self.i12, 0.0)
        j = set_blank_if_default(self.j, 0.0)
        nsm = set_blank_if_default(self.nsm, 0.0)

        C1 = set_blank_if_default(self.C1, 0.0)
        C2 = set_blank_if_default(self.C2, 0.0)

        D1 = set_blank_if_default(self.D1, 0.0)
        D2 = set_blank_if_default(self.D2, 0.0)

        E1 = set_blank_if_default(self.E1, 0.0)
        E2 = set_blank_if_default(self.E2, 0.0)

        F1 = set_blank_if_default(self.F1, 0.0)
        F2 = set_blank_if_default(self.F2, 0.0)

        K1 = set_blank_if_default(self.K1, 1e8)
        K2 = set_blank_if_default(self.K2, 1e8)

        list_fields = ['PBAR', self.pid, self.Mid(), self.A, i1, i2, j, nsm,
                       None, C1, C2, D1, D2, E1, E2, F1, F2, K1, K2, i12]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class PBARL(LineProperty):
    """
    .. todo:: doesnt support user-defined types
    """
    type = 'PBARL'
    validTypes = {
        "ROD": 1,
        "TUBE": 2,
        "I": 6,
        "CHAN": 4,
        "T": 4,
        "BOX": 4,
        "BAR": 2,
        "CROSS": 4,
        "H": 4,
        "T1": 4,
        "I1": 4,
        "CHAN1": 4,
        "Z": 4,
        "CHAN2": 4,
        "T2": 4,
        "BOX1": 6,
        "HEXA": 3,
        "HAT": 4,
        "HAT1": 5,
        "DBOX": 10,  # was 12
    }  # for GROUP="MSCBML0"

    def __init__(self, card=None, data=None, comment=''):
        LineProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Material ID
            self.mid = integer(card, 2, 'mid')
            self.group = string_or_blank(card, 3, 'group', 'MSCBMLO')
            #: Section Type (e.g. 'ROD', 'TUBE', 'I', 'H')
            self.Type = string(card, 4, 'Type')

            ndim = self.validTypes[self.Type]
            j = 9 + ndim + 1

            dims = []
            #dim_old = None  ## TODO: is there a default?
            for i in range(ndim):
                dim = double_or_blank(card, 9 + i, 'dim%i' % (i + 1))
                dims.append(dim)

            #: dimension list
            self.dim = dims
            assert len(dims) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(dims))
            assert len(dims) == len(self.dim), 'PBARL ndim=%s len(dims)=%s' % (ndim, len(self.dim))

            #: non-structural mass
            self.nsm = double_or_blank(card, 9 + ndim + 1, 'nsm', 0.0)

        else:
            self.pid = data[0]
            self.mid = data[1]
            self.group = data[2].strip()
            self.Type = data[3].strip()
            self.dim = list(data[4:-1])
            self.nsm = data[-1]
            #print("group = %r" % self.group)
            #print("Type  = %r" % self.Type)
            #print("dim = ",self.dim)
            #print(str(self))
            #print("*PBARL = ",data)
            #raise NotImplementedError('not finished...')
        if self.Type not in self.validTypes:
            msg = ('Invalid PBARL Type, Type=%s '
                   'validTypes=%s' % (self.Type, self.validTypes.keys()))
            raise RuntimeError(msg)

        if len(self.dim) != self.validTypes[self.Type]:
            msg = 'dim=%s len(dim)=%s Type=%s len(dimType)=%s' % (
                self.dim, len(self.dim), self.Type,
                self.validTypes[self.Type])
            raise RuntimeError(msg)

        assert None not in self.dim

    def cross_reference(self, model):
        msg = ' which is required by PBARL mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)

    def _verify(self, xref=False):
        pid = self.Pid()
        mid = self.Mid()
        A = self.Area()
        try:
            J = self.J()
        except NotImplementedError:
            msg = "J is not implemented for pid.type=%s pid.Type=%s" % (self.type, self.Type)
            print(msg)
            J = 0.0
        nsm = self.Nsm()
        mpl = self.MassPerLength()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(A, float), 'pid=%r' % A
        assert isinstance(J, float), 'cid=%r' % J
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        assert isinstance(mpl, float), 'mass_per_length=%r' % mpl

    def Area(self):
        """
        Gets the area :math:`A` of the CBAR.
        """
        return _bar_areaL('PBARL', self.Type, self.dim)

    def Nsm(self):
        """
        Gets the non-structural mass :math:`nsm` of the CBAR.
        """
        return self.nsm

    def MassPerLength(self):
        r"""
        Gets the mass per length :math:`\frac{m}{L}` of the CBAR.

        .. math:: \frac{m}{L} = A \rho + nsm
        """
        rho = self.Rho()
        area = self.Area()
        nsm = self.Nsm()
        return area * rho + nsm

    def I11(self):
        return self.I1()
        if self.Type in ['ROD']:
            assert len(self.dim) == 1, 'dim=%r' % self.dim
            r = self.dim[0]
            #Ix = pi*r**4/4.
            #J = pi*r**4/2.
            (Ix, Iy, Ixy) = self.I1_I2_I12()

        elif self.Type in ['BAR']:
            assert len(self.dim) == 2, 'dim=%r' % self.dim
            b, h = self.dim
            #Ix = self.I11()
            #Iy = self.I22()
            (Ix, Iy, Ixy) = self.I1_I2_I12()
            #print
            #J = Ix + Iy
        else:
            msg = 'I11 for Type=%r dim=%r on PBARL is not supported' % (self.Type, self.dim)
            raise NotImplementedError(msg)
        return Ix
        raise RuntimeError('who is calling this...')

    #def I12(self):
        #return self.I12()

    def _points(self, Type, dim):
        if Type == 'BAR':  # origin ar center
            (d1, d2) = dim
            Area = d1 * d2
            y1 = d2 / 2.
            x1 = d1 / 2.
            points = [  # start at upper right, go clockwise
                [x1, y1],    # p1
                [x1, y1],    # p2
                [-x1, -y1],  # p3
                [-x1, -y1],  # p4
            ]
        elif Type == 'CROSS':
            (d1, d2, d3, d4) = dim  # origin at center
            x1 = d2 / 2.
            x2 = d2 / 2. + d1
            y1 = -d3 / 2.
            y2 = -d4 / 2.
            y3 = d4 / 2.
            y4 = d3 / 2.
            points = [  # start at top right, go clockwise, down first
                [x1, y4],  # p1
                [x1, y3],  # p2
                [x2, y3],  # p3
                [x2, y2],  # p4
                [x1, y2],  # p5
                [x1, y1],  # p6

                [-x1, y1],  # p7
                [-x1, y2],  # p8
                [-x2, y2],  # p9
                [-x2, y3],  # p10
                [-x1, y3],  # p11
                [-x1, y4],  # p12
            ]
            Area = d2*d3 + 2*d1*d4
        elif Type in ['HEXA']:
            (d1, d2, d3) = dim
            x1 = d2 / 2. - d1
            x2 = d2 / 2.
            y1 = 0.
            y2 = d3 / 2.
            y3 = d3
            points = [  # start at upper center, go clockwise, diagonal down right
                [x1, y3],   # p1
                [x2, y2],   # p2
                [x1, y1],   # p3
                [x1, y1],   # p4
                [-x2, y2],  # p5
                [-x1, y3],  # p6
            ]
            Area = d1 * (d2 + 2 * d3)

        elif Type in ['I']:
            (d1, d2, d3) = dim
            raise NotImplementedError('PBARL Type=%r' % Type)

        elif Type in ['H']:
            (d1, d2, d3, d4) = dim
            x1 = d1 / 2.
            x2 = (d1 + d2) / 2.
            y3 = d4 / 2.
            y4 = d3 / 2.
            y1 = -y4
            y2 = -y3
            points = [  # start at top of H in dip, go clockwise, up first
                [x1, y3],  # p1 # right side
                [x1, y4],  # p2
                [x2, y4],  # p3
                [x2, y1],  # p4
                [x1, y1],  # p5
                [x1, y2],  # p6

                [-x1, y2],  # p7 # left side
                [-x1, y1],  # p8
                [-x2, y1],  # p9
                [-x2, y4],  # p10
                [-x1, y4],  # p11
                [-x1, y3],  # p12
            ]
            Area = d2 * d3 + d1 * d4

        elif Type in ['T2']:
            d1, d2, d3, d4 = dim  # check origin, y3 at bottom, x1 innner
            x1 = d4 / 2.
            x2 = d1 / 2.
            y1 = -d3 / 2.
            y2 = d3 / 2.
            y3 = -d3 / 2.
            points = [  # start at upper right, go clockwise
                [x1, y3],   # p1
                [x1, y2],   # p2
                [x2, y2],   # p3
                [x2, y1],   # p4
                [-x2, y1],  # p5
                [-x2, y2],  # p6
                [-x1, y2],  # p7
                [-x1, y3]   # p8
            ]
            Area = d1*d3 + (d2-d3)*d4
        else:
            msg = '_points for Type=%r dim=%r on PBARL is not supported' % (self.Type, self.dim)
            raise NotImplementedError(msg)
        return array(points), Area

    def J(self):
        if self.Type in ['ROD']:
            r = self.dim[0]
            Ixx = pi * r**4 / 4.
            Iyy = Ixx
            Ixy = 0.
        elif self.Type in ['TUBE']:
            rout, rin = self.dim
            #rin = rout - 2*t
            Ixx = pi * (rout**4 - rin**4) / 4.
            Iyy = Ixx
            Ixy = 0.
        elif self.Type in ['TUBE2']:
            rout, t = self.dim
            rin = rout - 2*t
            Ixx = pi * (rout**4 - rin**4) / 4.
            Iyy = Ixx
            Ixy = 0.
        elif self.Type in ['BOX']:
            (d1, d2, d3, d4) = self.dim
            hout = d2
            wout = d1
            hin = d2 - 2. * d3
            win = d1 - 2. * d4
            points, Area = self._points('BAR', [hout, wout])
            yi = points[0, :-1]
            yip1 = points[0, 1:]
            xi = points[1, :-1]
            xip1 = points[1, 1:]

            #: .. seealso:: http://en.wikipedia.org/wiki/Area_moment_of_inertia
            ai = xi*yip1 - xip1*yi
            Ixx1 = 1/12 * sum((yi**2 + yi * yip1 + yip1**2)*ai)
            Iyy1 = 1/12 * sum((xi**2 + xi * xip1 + xip1**2)*ai)
            #Ixy1 = 1/24*sum((xi*yip1 + 2*xi*yi + 2*xip1*yip1 + xip1*yi)*ai)

            points, Area = self._points('BAR', [hin, win])
            yi = points[0, :-1]
            yip1 = points[0, 1:]
            xi = points[1, :-1]
            xip1 = points[1, 1:]

            #: .. seealso:: http://en.wikipedia.org/wiki/Area_moment_of_inertia
            ai = xi*yip1 - xip1*yi
            Ixx2 = 1/12*sum((yi**2 + yi*yip1+yip1**2)*ai)
            Iyy2 = 1/12*sum((xi**2 + xi*xip1+xip1**2)*ai)
            #Ixy2 = 1/24*sum((xi*yip1 + 2*xi*yi + 2*xip1*yip1 + xip1*yi)*ai)

            Ixx = Ixx1 - Ixx2
            Iyy = Iyy1 - Iyy2
            #Ixy = Ixy1 - Ixy2

        #elif self.Type in ['BAR']:
            #assert len(self.dim) == 2, 'dim=%r' % self.dim
            #b, h = self.dim
            #(Ix, Iy, Ixy) = self.I1_I2_I12()
            #J = Ix + Iy
        elif self.Type in ['BAR', 'CROSS', 'HEXA', 'T2', 'H']:
            points, Area = self._points(self.Type, self.dim)
            yi = points[0, :-1]
            yip1 = points[0, 1:]

            xi = points[1, :-1]
            xip1 = points[1, 1:]

            #: .. seealso:: http://en.wikipedia.org/wiki/Area_moment_of_inertia
            ai = xi*yip1 - xip1*yi
            Ixx = 1/12*sum((yi**2 + yi*yip1+yip1**2)*ai)
            Iyy = 1/12*sum((xi**2 + xi*xip1+xip1**2)*ai)
            #Ixy = 1/24*sum((xi*yip1 + 2*xi*yi + 2*xip1*yip1 + xip1*yi)*ai)
        else:
            msg = 'J for Type=%r dim=%r on PBARL is not supported' % (self.Type, self.dim)
            raise NotImplementedError(msg)

        #: .. seealso:: http://en.wikipedia.org/wiki/Perpendicular_axis_theorem
        J = Ixx + Iyy
        return J

    def I22(self):
        return self.I2()

    def write_code_aster(self, iCut=0, iFace=0, iStart=0):  # PBARL
        msg = '# BAR Type=%s pid=%s\n' % (self.type, self.pid)
        msg2 = ''
        msg += self.CA_Section(iFace, iStart, self.dim)
        iFace += 1
        iStart += len(self.dim)

        msg += self.CA_Section(iFace, iStart, self.dim)
        iFace += 1
        msg2 += 'Cut_%s = geompy.MakeCut(Face_%i, Face_%i)\n' % (
            iCut + 1, iFace + 1, iFace + 2)
        msg2 += "geompy.addToStudy(Cut_%i,  'Cut_%i')\n" % (
            iCut + 1, iCut + 1)
        iStart += len(self.dim)
        return msg + msg2, iCut, iFace, iStart

    def raw_fields(self):
        list_fields = ['PBARL', self.pid, self.Mid(), self.group, self.Type,
                       None, None, None, None] + self.dim + [self.nsm]
        return list_fields

    def repr_fields(self):
        group = set_blank_if_default(self.group, 'MSCBMLO')
        list_fields = ['PBARL', self.pid, self.Mid(), group, self.Type, None,
                       None, None, None] + self.dim + [self.nsm]
        return list_fields

    def write_card(self, size, is_double):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class PBEAM3(LineProperty):  # not done, cleanup
    type = 'PBEAM3'

    def __init__(self, card=None, data=None, comment=''):
        LineProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.pid = integer(card, 1, 'pid')
            self.mid = integer(card, 2, 'mid')

            self.A = double(card, 3, 'A')
            self.iz = double(card, 4, 'Iz')
            self.iy = double(card, 5, 'Iy')
            self.iyz = double_or_blank(card, 6, 'Iyz', 0.0)
            self.j = double_or_blank(card, 7, 'J', self.iy + self.iz)
            self.nsm = double_or_blank(card, 8, 'nsm', 0.0)

            self.cy = double(card, 9, 'cy')
            self.cz = double(card, 10, 'cz')
            self.dy = double(card, 11, 'dy')
            self.dz = double(card, 12, 'dz')
            self.ey = double(card, 13, 'ey')
            self.dz = double(card, 14, 'ez')
            self.fy = double(card, 15, 'fy')
            self.fz = double(card, 16, 'fz')
            # more...
        else:
            raise NotImplementedError(data)

    def Nsm(self):
        """
        Gets the non-structural mass :math:`nsm`.
        .. warning:: nsm field not supported fully on PBEAM3 card
        """
        return self.nsm

    def cross_reference(self, model):
        msg = ' which is required by PBEAM3 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)

    def repr_fields(self):
        """.. todo:: not done"""
        list_fields = ['PBEAM3', self.pid, self.Mid(), ]  # other
        return list_fields


class PBEND(LineProperty):
    type = 'PBEND'

    def __init__(self, card=None, data=None, comment=''):
        LineProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.pid = integer(card, 1, 'pid')
            self.mid = integer(card, 2, 'mid')

            value3 = integer_or_double(card, 3, 'A_FSI')
            if isinstance(value3, float):
                self.beamType = 1
                #: Area of the beam cross section
                self.A = double(card, 3, 'A')

                #: Area moments of inertia in planes 1 and 2.
                self.i1 = double(card, 4, 'I1')
                self.i2 = double(card, 5, 'I2')

                #: Torsional stiffness :math:`J`
                self.j = double(card, 6, 'J')

                # line2
                #: The r,z locations from the geometric centroid for stress
                #: data recovery.
                self.c1 = double(card, 9, 'c1')
                self.c2 = double(card, 10, 'c2')
                self.d1 = double(card, 11, 'd1')
                self.d2 = double(card, 12, 'd2')
                self.e1 = double(card, 13, 'e1')
                self.e2 = double(card, 14, 'e2')
                self.f1 = double(card, 15, 'f1')
                self.f2 = double(card, 16, 'f2')

                # line 3
                #: Shear stiffness factor K in K*A*G for plane 1.
                self.k1 = double(card, 17, 'k1')
                #: Shear stiffness factor K in K*A*G for plane 2.
                self.k2 = double(card, 18, 'k2')

                #: Nonstructural mass per unit length.
                self.nsm = double(card, 19, 'nsm')

                #: Radial offset of the geometric centroid from points GA and GB.
                self.rc = double(card, 20, 'rc')

                #: Offset of the geometric centroid in a direction perpendicular
                #: to the plane of points GA and GB and vector v
                self.zc = double(card, 21, 'zc')

                #: Radial offset of the neutral axis from the geometric
                #: centroid, positive is toward the center of curvature
                self.deltaN = double(card, 22, 'deltaN')

            elif isinstance(value3, integer_types):  # alternate form
                self.beamType = 2
                #: Flag selecting the flexibility and stress intensification
                #: factors. See Remark 3. (Integer = 1, 2, or 3)
                self.fsi = integer(card, 3, 'fsi')
                assert self.fsi in [1, 2, 3]

                #: Mean cross-sectional radius of the curved pipe
                self.rm = double(card, 4, 'rm')

                #: Wall thickness of the curved pipe
                self.t = double(card, 5, 't')

                #: Internal pressure
                self.p = double(card, 6, 'p')

                # line3
                # Non-structural mass :math:`nsm`
                self.nsm = double(card, 11, 'nsm')
                self.rc = double(card, 12, 'rc')
                self.zc = double(card, 13, 'zc')
            else:
                raise RuntimeError('Area/FSI on CBEND must be defined...')

            #: Bend radius of the line of centroids
            self.rb = double_or_blank(card, 7, 'rb')

            #: Arc angle :math:`\theta_B` of element  (optional)
            self.thetab = double_or_blank(card, 8, 'thetab')
            assert len(card) <= 23, 'len(PBEND card) = %i' % len(card)

        else:
            raise NotImplementedError(data)

    #def Nsm(self):
        #""".. warning:: nsm field not supported fully on PBEND card"""
        #raise RuntimeError(self.nsm[0])
        #return self.nsm

    def cross_reference(self, model):
        msg = ' which is required by PBEND mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)

    def repr_fields(self):
        list_fields = ['PBEND', self.pid, self.Mid(), ]  # other
        if self.beamType == 1:
            list_fields += [self.A, self.i1, self.i2, self.j, self.rb,
                            self.thetab, self.c1, self.c2, self.d1, self.d2,
                            self.e1, self.e2, self.f1, self.f2, self.k1, self.k2,
                            self.nsm, self.rc, self.zc, self.deltaN]
        elif self.beamType == 2:
            list_fields += [self.fsi, self.rm, self.t, self.p, self.rb,
                            self.thetab, None, None, self.nsm, self.rc, self.zc]
        else:
            raise ValueError('only beamType=1 and 2 supported')
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)
