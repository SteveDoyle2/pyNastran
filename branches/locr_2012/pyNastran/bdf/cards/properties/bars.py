# pylint: disable=C0103,R0902,R0904,R0914
"""

All beam properties are defined in this file.  This includes:
 *   PBEAM
 *   PBEAML
 *   PBAR
 *   PBARL
 *   PROD
 *   PTUBE

All beams are LineProperty objects.
Multi-segment beams are IntegratedLineProperty objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from itertools import izip, count
from numpy import zeros, pi

from pyNastran.bdf.fieldWriter import (set_blank_if_default,
                                       set_default_if_blank)
from pyNastran.general.mathematics import integrate_line, integrate_positive_line
from pyNastran.bdf.cards.baseCard import Property


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
    calculates the moment of inertia for a section about the CG
    @param sections
      [[b,h,y,z]_1,...] y,z is the centroid (x in the direction of the beam,
      y right, z up)
    @retval interiaParameters
        list of [Area,Iyy,Izz,Iyz]
    @see http://www.webs1.uidaho.edu/mindworks/Machine_Design/Posters/PDF/Moment%20of%20Inertia.pdf
    
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
    type = 'LineProperty'

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

    def D_bending(self):
        pass

    def D_axial(self):
        pass

    def D_thermal(self):
        pass

    def D_shear(self):
        pass

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

    def I12(self):
        return self.i12

    def E(self):
        return self.mid.E

    def G(self):
        return self.mid.G

    def Nu(self):
        return self.mid.nu

    def CA_Section(self, iFace, iStart, dims):
        """
        @code
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
        @endcode
        """
        msg1 = ''
        msg2 = 'Face_%s = geompy.MakeFaceHW(' % (iFace + 1)
        for (i, dim) in enumerate(dims):
            msg1 += 'D%s = %s\n' % (iStart + i, dim)
            msg2 += 'D%s,' % (iStart + i)
        msg2 += '1)\n'
        msg2 += "geompy.addToStudy(Face_%i, 'Face_%i')\n" % (iFace, iFace)
        return (msg1 + msg2)

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
            sections = []
            h1 = dim[5]  # d2
            w1 = dim[2]                         # |  ------------
            y1 = dim[0] / 2. - h1                   # |  |    A     | d5
            sections.append([w1, h1, 0., y1])   # |  ------------
                                                # |     >| |<--d3
            h3 = dim[4]                         # |      |B|
            w3 = dim[1]                         # | d1   | |
            y3 = -dim[0] / 2. + h3                  # |      | |
            sections.append([w3, h3, 0., y1])   # |   ----------
                                                # |   |   C    |  d5
            h2 = dim[0] - h1 - h3                   # |   ----------
            w2 = dim[3]  # d1
            sections.append([w2, h2, 0., 0.])

            (A, Iyy, Izz, Iyz) = getInertiaRectangular(sections)
            assert Iyz == 0.

        elif self.Type == 'BAR':  # *-------*
            h1 = dim[1]  # |       |
            w1 = dim[0]  # |       |h1
            A = h1 * w1  # |       |
            ## *-------*          I_{xx}=\frac{bh^3}{12}
            Iyy = 1 / 12. * w1 * h1 ** 3
            ## w1             I_{yy}=\frac{hb^3}{12}
            Izz = 1 / 12. * h1 * w1 ** 3  
            Iyz = 0.  ## @todo is the Ixy of a bar 0 ???

        else:
            msg = 'Type=%s is not supported for %s class...' % (self.Type,
                                                                self.type)
            raise NotImplementedError(msg)
        return (A, Iyy, Izz, Iyz)

    def I1(self):
        return self.I1_I2()[0]

    def I2(self):
        return self.I1_I2()[1]

    def I1_I2(self):
        """
        @code
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
        I1 = 1/12*b*h^3
        I2 = 1/12*h*b^3
        @endcode
        """
        dim = self.dim
        if self.Type == 'BAR':
            I1 = 1 / 12. * dim[0] * dim[1] ** 3
            I2 = 1 / 12. * dim[1] * dim[0] ** 3
        return(I1, I2)

    def areaL(self, dim):
        """
        Area method for the PBARL and PBEAML classes (pronounced Area-L)
        @param self the object pointer
        @param dim a list of the dimensions associated with self.Type
        @retval Area of the given cross section
        """
        try:
            if   self.Type == 'ROD':
                A = pi * dim[0] ** 2
            elif self.Type == 'TUBE':
                A = pi * (dim[0] ** 2 - dim[1] ** 2)
            elif self.Type == 'I':
                h1 = dim[5]
                w1 = dim[2]

                h3 = dim[4]
                w3 = dim[1]

                h2 = dim[0] - h1 - h3
                w2 = dim[3]
                A = h1 * w1 + h2 * w2 + h3 * w3
            elif self.Type == 'CHAN':
                h1 = dim[3]
                w1 = dim[0]

                h3 = h1
                w3 = w1
                h2 = dim[1] - h1 - h3
                w2 = dim[2]
                A = h1 * w1 + h2 * w2 + h3 * w3
            elif self.Type == 'T':
                h1 = dim[2]
                w1 = dim[0]

                h2 = dim[1] - h1
                w2 = dim[3]
                A = h1 * w1 + h2 * w2
            elif self.Type == 'BOX':
                h1 = dim[2]
                w1 = dim[0]

                h2 = dim[1] - 2 * h1
                w2 = dim[3]
                A = 2 * (h1 * w1 + h2 * w2)
            elif self.Type == 'BAR':
                h1 = dim[1]
                w1 = dim[0]
                A = h1 * w1
            elif self.Type == 'CROSS':
                h1 = dim[2]
                w1 = dim[1]

                h2 = dim[3]
                w2 = dim[0]
                A = h1 * w1 + h2 * w2
            elif self.Type == 'H':
                h1 = dim[2]
                w1 = dim[1]

                h2 = dim[3]
                w2 = dim[0]
                A = h1 * w1 + h2 * w2
            elif self.Type == 'T1':
                h1 = dim[0]
                w1 = dim[2]

                h2 = dim[3]
                w2 = dim[1]
                A = h1 * w1 + h2 * w2
            elif self.Type == 'I1':
                h2 = dim[2]
                w2 = dim[1]

                h1 = dim[3] - h2
                w1 = dim[0] + w2
                A = h1 * w1 + h2 * w2
            elif self.Type == 'CHAN1':
                h2 = dim[2]
                w2 = dim[1]

                h1 = dim[3] - h2
                w1 = dim[0] + w2
                A = h1 * w1 + h2 * w2
            elif self.Type == 'Z':
                h2 = dim[2]
                w2 = dim[1]

                h1 = dim[3] - h2
                w1 = dim[0]
                A = h1 * w1 + h2 * w2
            elif self.Type == 'CHAN2':
                h2 = dim[1]
                w2 = dim[3]

                h1 = dim[2] - h2
                w1 = dim[0] * 2
                A = h1 * w1 + h2 * w2

            elif self.Type == 'T2':
                h1 = dim[3]
                w1 = dim[1]

                h2 = h1 - dim[2]
                w2 = dim[0]
                A = h1 * w1 + h2 * w2
            elif self.Type == 'BOX1':
                h1 = dim[2]  # top
                w1 = dim[0]

                h2 = dim[3]  # btm
                A1 = (h1 + h2) * w1

                h3 = dim[1] - h1 - h2  # left
                w3 = dim[5]

                w4 = dim[4]  # right
                A2 = h3 * (w3 + w4)
                A = A1 + A2
            elif self.Type == 'HEXA':
                hBox = dim[2]
                wBox = dim[1]

                wTri = dim[0]
                A = hBox * wBox - wTri * hBox
            elif self.Type == 'HAT':
                w = dim[1]      # constant width (note h is sometimes w)
                h1 = w           # horizontal lower bar
                w1 = dim[3]

                h2 = dim[0] - 2 * w  # vertical bar
                w2 = w

                h3 = w           # half of top bar
                w3 = dim[2] / 2.

                A = 2 * (h1 * w1 + h2 * w2 + h3 * w3)  # symmetrical box
            elif self.Type == 'HAT1':
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

            elif self.Type == 'DBOX':
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
                msg = 'Type=%s is not supported for %s class...' % (self.Type,
                                                                    self.type)
                raise NotImplementedError(msg)
        except IndexError as e:
            msg = 'There was an error extracting fields'
            msg += ' from a %s dim=%s for a %s' % (self.Type, dim, self.type)
            msg += '-' * 80 + '\n'
            msg += 'Traceback:\n%s' % (e.message)
            raise IndexError(msg)
        return A


class IntegratedLineProperty(LineProperty):
    type = 'IntegratedLineProperty'

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
        #print "xxb = ",self.xxb
        #print "nsm = ",self.nsm
        nsm = integrate_positive_line(self.xxb, self.nsm)
        return nsm


class PROD(LineProperty):
    type = 'PROD'

    def __init__(self, card=None, data=None):
        LineProperty.__init__(self, card, data)

        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.A = card.field(3)
            self.j = card.field(4, 0.0)
            self.c = card.field(5, 0.0)
            self.nsm = card.field(6, 0.0)
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.A = data[2]
            self.j = data[3]
            self.c = data[4]
            self.nsm = data[5]
        ###

    #def Radius(self):
        #"""assumes circular cross section - probably will remove this"""
        #return (self.A/pi)**0.5

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def writeCodeAster(self, iCut, iFace, iStart):  # PROD
        msg = ''
        msg += "    POUTRE=_F(GROUP_MA='P%s', # PROD\n" % (self.pid)
        msg += "              SECTION='CERCLE',  # circular section\n"
        msg += "              CARA=('R')   # radius\n"
        msg += "              VALE=(%g),),\n" % (self.Radius())

        #msg += "              SECTION='GENERALE',\n"
        #msg += "              CARA=('A','IY','IZ','JX')\n"
        #msg += "              VALE=(%g,  %g,  %g,  %g),\n"  %(self.Area(),self.I11(),self.I22(),self.J())
        #msg += "                    CARA='VECT_Y'),\n"
        #msg += "                    VALE=(1.0,0.0,0.0,),),\n"
        return (msg, iCut, iFace, iStart)

    def rawFields(self):
        fields = ['PROD', self.pid, self.Mid(), self.A, self.j, self.c,
                  self.nsm]
        return fields

    def reprFields(self):
        j = set_blank_if_default(self.j, 0.0)
        c = set_blank_if_default(self.c, 0.0)
        nsm = set_blank_if_default(self.nsm, 0.0)
        fields = ['PROD', self.pid, self.Mid(), self.A, j, c, nsm]
        return fields


class PTUBE(LineProperty):
    type = 'PTUBE'

    def __init__(self, card=None, data=None):
        LineProperty.__init__(self, card, data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.OD1 = card.field(3)
            self.t = card.field(4, self.OD1 / 2.)
            self.nsm = card.field(5, 0.0)
            self.OD2 = card.field(6, self.OD1)
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.OD1 = data[2]
            self.t = data[3]
            self.nsm = data[4]
            self.OD2 = self.OD1
            #self.OD2 = data[5]  ## @note quirk to this one...

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def Area(self):
        A = (self.area1() + self.area2()) / 2.

        #A1 = pi*D1^2/4 - pi*((D1-2t)^2)/4
        #A2 = pi*D2^2/4 - pi*((D2-2t)^2)/4
        #A = A1+A2

        #A = pi*D1*t/2 + pi*D2*t/2 - pi*t
        #A = pi*t*(D1/2 + D2/2 - t)
        #A = pi*t*( (D1+D2)/2.-t )

        #A2 = pi*t*( (D1+D2)/2.-t )

        #if A != A2:
            #msg = 'AREA method has problem in PTUBE Aold=%s != Anew=%s' %(A,A2)
            #raise RuntimeError(msg)
        return A

    def area1(self):
        """@todo remove after verifying formula..."""
        Dout = self.OD1
        Din = Dout - 2 * self.t
        A1 = pi / 4. * (Dout * Dout - Din * Din)
        return A1

    def area2(self):
        """@todo remove after verifying formula..."""
        Dout = self.OD2
        Din = Dout - 2 * self.t
        A2 = pi / 4. * (Dout * Dout - Din * Din)
        return A2

    def massMatrix(self):
        """@todo not done"""
        m = zeros(6, 6)
        m[0, 0] = 1.
        return m

    def rawFields(self):
        fields = ['PTUBE', self.pid, self.Mid(), self.OD1, self.t, self.nsm,
                  self.OD2]
        return fields

    def reprFields(self):
        t = set_blank_if_default(self.t, self.OD1 / 2.)
        nsm = set_blank_if_default(self.nsm, 0.0)
        OD2 = set_blank_if_default(self.OD2, self.OD1)
        fields = ['PTUBE', self.pid, self.Mid(), self.OD1, t, nsm, OD2]
        return fields


class PBAR(LineProperty):
    type = 'PBAR'

    def __init__(self, card=None, data=None):
        """
        @todo
            support solution 600 default
            do a check for mid -> MAT1      for structural
            do a check for mid -> MAT4/MAT5 for thermal
            
        """
        LineProperty.__init__(self, card, data)
        if card:
            ## property ID -> use Pid()
            self.pid = card.field(1)
            ## material ID -> use Mid()
            self.mid = card.field(2)
            ## Area -> use Area()
            self.A = card.field(3, 0.0)
            ## Izz -> use Izz()
            self.i1 = card.field(4, 0.0)
            ## Iyy -> use Iyy()
            self.i2 = card.field(5, 0.0)
            ## Polar Moment of Inertia J -> use J()
            # default=1/2(I1+I2) for SOL=600,otherwise 0.0
            self.j = card.field(6, 0.0)
            ## nonstructral mass -> use Nsm()
            self.nsm = card.field(7, 0.0)

            self.C1 = card.field(9, 0.0)
            self.C2 = card.field(10, 0.0)
            self.D1 = card.field(11, 0.0)
            self.D2 = card.field(12, 0.0)
            self.E1 = card.field(13, 0.0)
            self.E2 = card.field(14, 0.0)
            self.F1 = card.field(15, 0.0)
            self.F2 = card.field(16, 0.0)

            ## default=infinite
            self.K1 = card.field(17, 1e8)
            ## default=infinite
            self.K2 = card.field(18, 1e8)
            self.i12 = card.field(19, 0.0)
            if self.A == 0.0:
                assert self.K1 is None
                assert self.K2 is None
            ###
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.A = data[2]
            self.i1 = data[3]
            self.i2 = data[4]
            self.j = data[5]

            self.nsm = data[6]
            #self.fe  = data[7] ## @todo not documented....
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
        ###

    def MassPerLength(self):
        r"""
        \f[ \frac{m}{L} = \rho A+nsm \f]
        """
        A = self.Area()
        rho = self.Rho()
        nsm = self.Nsm()
        return rho * A + nsm

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def Area(self):
        return self.A

    #def Nsm(self):
    #    return self.nsm

    #def J(self):
    #    return self.j

    #def Izz(self):
    #    return self.i1

    #def Iyy(self):
    #    return self.i2

    #def Iyz(self):
    #    return self.i12

    def writeCodeAster(self):  # PBAR
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
        return (msg)

    def rawFields(self):
        fields = ['PBAR', self.pid, self.Mid(), self.A, self.i1, self.i2,
                  self.j, self.nsm, None, self.C1, self.C2, self.D1, self.D2,
                  self.E1, self.E2, self.F1, self.F2, self.K1, self.K2,
                  self.i12]
        return fields

    def reprFields(self):
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

        #print "line3 = ",line3

        fields = ['PBAR', self.pid, self.Mid(), self.A, i1, i2, j, nsm, None,
                  C1, C2, D1, D2, E1, E2, F1, F2,
                  K1, K2, i12]

        return fields


class PBARL(LineProperty):
    """
    @warning doesnt support user-defined types
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

    def __init__(self, card=None, data=None):
        LineProperty.__init__(self, card, data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.group = card.field(3, 'MSCBMLO')
            self.Type = card.field(4)

            #nDim = len(self.dim)-1
            nDim = self.validTypes[self.Type]
            self.dim = card.fields(9, 9 + nDim + 1)
            self.nsm = card.field(9 + nDim + 1, 0.0)

            #self.dim = fields(9)
            if nDim > 0:
                self.nsm = set_default_if_blank(self.dim.pop(), 0.0)
            else:
                self.nsm = 0.0
            ###
            assert isinstance(self.nsm, float)
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.group = data[2].strip()
            self.Type = data[3].strip()
            self.dim = list(data[4:-1])
            self.nsm = data[-1]
            #print "group = |%s|" %(self.group)
            #print "Type  = |%s|" %(self.Type)
            #print "dim = ",self.dim
            #print str(self)
            #print "*PBARL = ",data
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
        self.mid = model.Material(self.mid)

    def Area(self):
        return self.areaL(self.dim)

    def Nsm(self):
        return self.nsm

    def MassPerLength(self):
        """
        @code
        mass = L*(Area*rho+nsm)
        mass/L = Area*rho+nsm
        @endcode
        """
        rho = self.Rho()
        area = self.Area()
        nsm = self.Nsm()
        return area * rho + nsm

    def I11(self):
        return None

    def I12(self):
        return None

    def J(self):
        return None

    def I22(self):
        return None

    def writeCodeAster(self, iCut=0, iFace=0, iStart=0):  # PBARL
        msg = '# BAR Type=%s pid=%s\n' % (self.type, self.pid)
        msg2 = ''
        (msg) += self.CA_Section(iFace, iStart, self.dim)
        iFace += 1
        iStart += len(self.dim)

        (msg) += self.CA_Section(iFace, iStart, self.dim)
        iFace += 1
        msg2 += 'Cut_%s = geompy.MakeCut(Face_%i, Face_%i)\n' % (
            iCut + 1, iFace + 1, iFace + 2)
        msg2 += "geompy.addToStudy(Cut_%i,  'Cut_%i')\n" % (
            iCut + 1, iCut + 1)
        iStart += len(self.dim)
        return (msg + msg2, iCut, iFace, iStart)

    def rawFields(self):
        fields = ['PBARL', self.pid, self.Mid(), self.group, self.Type, None,
                  None, None, None, ] + self.dim + [self.nsm]
        return fields

    def reprFields(self):
        group = set_blank_if_default(self.group, 'MSCBMLO')
        fields = ['PBARL', self.pid, self.Mid(), group, self.Type, None, None,
                  None, None, ] + self.dim + [self.nsm]
        return fields


class PBCOMP(LineProperty):
    type = 'PBCOMP'

    def __init__(self, card=None, data=None):
        LineProperty.__init__(self, card, data)
        if card:
            ## Property ID
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.A = card.field(3, 0.0)
            self.i1 = card.field(4, 0.0)
            self.i2 = card.field(5, 0.0)
            self.i12 = card.field(6, 0.0)
            self.j = card.field(7, 0.0)
            self.nsm = card.field(8, 0.0)
            self.k1 = card.field(9, 1.0)
            self.k2 = card.field(10, 1.0)
            self.m1 = card.field(11, 0.0)
            self.m2 = card.field(12, 0.0)
            self.n1 = card.field(13, 0.0)
            self.n2 = card.field(14, 0.0)
            self.symopt = card.field(15, 0)
            self.y = []
            self.z = []
            self.c = []
            self.mids = []

            fields = card.fields(17)
            nfields = len(fields)
            nrows = nfields // 8
            if nfields % 8 > 0:
                nrows += 1
            
            for row in xrange(nrows):
                i = 8*row + 17
                yi = card.field(i)
                zi = card.field(i+1)
                ci = card.field(i+2, 0.0)
                mid = card.field(i+3, self.mid)
                self.y.append(yi)
                self.z.append(zi)
                self.c.append(ci)
                self.mids.append(mid)
            
    def MassPerLength(self):
        return self.nsm+self.mid.Rho()*self.A

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def rawFields(self):
        fields = ['PBCOMP', self.pid, self.Mid(), self.A, self.i1, self.i2,
             self.i12, self.j, self.nsm, self.k1, self.k2, self.m1, self.m2,
             self.n1, self.n2, self.symopt, None]
        for (yi, zi, ci, mid) in zip(self.y,self.z,self.c,self.mids):
            fields += [yi, zi, ci, mid, None, None, None, None]
        return fields
                
    def reprFields(self):
        area = set_blank_if_default(self.A, 0.0)
        j = set_blank_if_default(self.j, 0.0)
        i1 = set_blank_if_default(self.i1, 0.0)
        i2 = set_blank_if_default(self.i2, 0.0)
        i12 = set_blank_if_default(self.i12, 0.0)
        nsm = set_blank_if_default(self.nsm, 0.0)

        k1 = set_blank_if_default(self.k1, 1.0)
        k2 = set_blank_if_default(self.k2, 1.0)

        m1 = set_blank_if_default(self.m1, 0.0)
        m2 = set_blank_if_default(self.m2, 0.0)

        n1 = set_blank_if_default(self.n1, 0.0)
        n2 = set_blank_if_default(self.n2, 0.0)

        symopt = set_blank_if_default(self.symopt, 0)
        
        fields = ['PBCOMP', self.pid, self.Mid(), area, i1, i2, i12, j, nsm,
                  k1, k2, m1, m2, n1, n2, symopt, None]

        for (yi, zi, ci, mid) in zip(self.y,self.z,self.c,self.mids):
            ci = set_blank_if_default(ci, 0.0)
            fields += [yi, zi, ci, mid, None, None, None, None]
        return fields

class PBEAM(IntegratedLineProperty):
    type = 'PBEAM'

    def __init__(self, card=None, data=None):
        """
        @todo fix 0th entry of self.so, self.xxb
        """
        IntegratedLineProperty.__init__(self, card, data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            #print "pid = ",self.pid

            nFields = card.nFields() - 1  # -1 for PBEAM field

            self.so = ['YES']  # @todo what are these values (so[0],xxb[0])???
            self.xxb = [0.]
            self.A = [card.field(3)]
            self.i1 = [card.field(4, 0.0)]
            self.i2 = [card.field(5, 0.0)]
            self.i12 = [card.field(6, 0.0)]
            self.j = [card.field(7, 0.0)]
            self.nsm = [card.field(8, 0.0)]
            self.c1 = [card.field(9, 0.0)]
            self.c2 = [card.field(10, 0.0)]
            self.d1 = [card.field(11, 0.0)]
            self.d2 = [card.field(12, 0.0)]
            self.e1 = [card.field(13, 0.0)]
            self.e2 = [card.field(14, 0.0)]
            self.f1 = [card.field(15, 0.0)]
            self.f2 = [card.field(16, 0.0)]

            #fields = card.fields(0)
            #print "fieldsPBEAM = ",fields
            #fieldsMid = fields[16:]
            #print "fieldsMid = ",fieldsMid

            #fields = card.fields(9)
            #print ""
            #print "  nFields = ",nFields
            #nFields = card.nFields()-16 # 17+16 (leading + trailing fields)
            # counting continuation cards
            nMajor = nFields // 16
            nLeftover = nFields % 16
            #print "  nMajor=%s nLeftover=%s" %(nMajor,nLeftover)
            if nLeftover == 0:
                nMajor -= 1

            if nMajor == 0:
                nMajor = 1

            x = (nMajor) * 16 + 1
            if card.field(x) in ['YES', 'YESA', 'NO']:  # there is no footer
                nMajor += 1
                x += 16

            #print "  nMajor = ",nMajor
            for nRepeated in xrange(1, nMajor):
                #print "  adding a major"
                ## field 17 is the first possible so
                nStart = nRepeated * 16 + 1
                propFields = card.fields(nStart, nStart + 16)
                #print "propFields = ",propFields

                #print "  so = ",propFields[0]
                if propFields[0] not in [None, 'YES', 'YESA', 'NO']:
                    msg = "SO=%s for PBEAM pid=%s" % (propFields[0], self.pid)
                    raise RuntimeError(msg)
                self.so.append(propFields[0])
                self.xxb.append(propFields[1])
                self.A.append(propFields[2])
                self.i1.append(set_default_if_blank(propFields[3], 0.0))
                self.i2.append(set_default_if_blank(propFields[4], 0.0))
                self.i12.append(set_default_if_blank(propFields[5], 0.0))
                self.j.append(set_default_if_blank(propFields[6], 0.0))
                self.nsm.append(set_default_if_blank(propFields[7], 0.0))
                self.c1.append(set_default_if_blank(propFields[8], 0.0))
                self.c2.append(set_default_if_blank(propFields[9], 0.0))
                self.d1.append(set_default_if_blank(propFields[10], 0.0))
                self.d2.append(set_default_if_blank(propFields[11], 0.0))
                self.e1.append(set_default_if_blank(propFields[12], 0.0))
                self.e2.append(set_default_if_blank(propFields[13], 0.0))
                self.f1.append(set_default_if_blank(propFields[14], 0.0))
                self.f2.append(set_default_if_blank(propFields[15], 0.0))
            #print("nRepeated = %s" %(nRepeated))

            # footer fields
            self.k1 = card.field(x, 1.0)

            if self.k1 in ['YES', 'YESA', 'NO']:
                msg = 'error reading PBEAM card pid=%s' % (self.pid)
                raise RuntimeError(msg)
            #print "  k1 = ",self.k1
            self.k2 = card.field(x + 1, 1.0)
            self.s1 = card.field(x + 2, 0.0)
            self.s2 = card.field(x + 3, 0.0)
            self.nsia = card.field(x + 4, 0.0)
            self.nsib = card.field(x + 5, self.nsia)
            self.cwa = card.field(x + 6, 0.0)
            self.cwb = card.field(x + 7, self.cwa)

            self.m1a = card.field(x + 8, 0.0)
            self.m2a = card.field(x + 9, self.m1a)
            self.m1b = card.field(x + 10, 0.0)
            self.m2b = card.field(x + 11, self.m1b)
            self.n1a = card.field(x + 12, 0.0)
            self.n2a = card.field(x + 13, self.n1a)
            self.n1b = card.field(x + 14, 0.0)
            self.n2b = card.field(x + 15, self.n1b)
        else:
            raise NotImplementedError('not supported')
        ###

    #def Area(self):
    #    """@warning area field not supported fully on PBEAM card"""
    #    #raise RuntimeError(self.A[0])
    #    return self.A[0]

    #def Nsm(self):
    #    """@warning nsm field not supported fully on PBEAM card"""
    #    #raise RuntimeError(self.nsm[0])
    #    return self.nsm[0]

    def MassPerLength(self):
        """
        mass = L*(Area*rho+nsm)
        mass/L = Area*rho+nsm
        """
        rho = self.Rho()
        massPerLs = []
        for (area, nsm) in izip(self.A, self.nsm):
            massPerLs.append(area * rho + nsm)
        massPerL = integrate_positive_line(self.xxb, massPerLs)
        return massPerL

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def writeCodeAster(self):  # PBEAM
        a = self.Area()
        iy = self.I11()
        iz = self.I22()
        j = self.J()
        msg = ''
        msg += "    POUTRE=_F(GROUP_MA='P%s', # PBEAM\n" % (self.pid)
        msg += "              SECTION='GENERALE',\n"
        msg += "              CARA=('A','IY','IZ','JX'), # area, moments of inertia\n"
        msg += "              VALE=(%g,  %g,  %g,  %g),\n" % (a, iy, iz, j)

        msg += "              ORIENTATION=_F( \n"
        msg += "                  CARA=('VECT_Y'), # direction of beam ???\n"  # @todo is this correct
        msg += "                  VALE=(1.0,0.0,0.0,)"

        if [self.n1a, self.n1b] != [0., 0.]:
            msg += "              \n),\n"
            msg += "              CARA=('AX','AY'), # shear centers\n"
            msg += "              VALE=(%g, %g),\n" % (self.n1a, self.n1b)
            msg += "             ),\n"
        else:
            msg += " )),\n"
        return (msg)

    def rawFields(self):
        fields = ['PBEAM', self.pid, self.Mid()]

        for (so, xxb, A, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2,
             f1, f2) in izip(self.so, self.xxb, self.A, self.i1, self.i2,
                             self.i12, self.j, self.nsm, self.c1, self.c2, self.d1, self.d2,
                             self.e1, self.e2, self.f1, self.f2):
            fields += [so, xxb, A, i1, i2, i12, j, nsm, c1, c2, d1, d2,
                       e1, e2, f1, f2]
        ###
        footer = [self.k1, self.k2, self.s1, self.s2, self.nsia, self.nsib,
                  self.cwa, self.cwb, self.m1a, self.m2a, self.m1b, self.m2b,
                  self.n1a, self.n2a, self.n1b, self.n2b]
        fields += footer
        return fields

    def reprFields(self):
        fields = ['PBEAM', self.pid, self.Mid()]

        i = 0
        for (so, xxb, A, i1, i2, i12, j, nsm, c1, c2, d1, d2, e1, e2, f1, f2) in izip(
            self.so, self.xxb, self.A, self.i1, self.i2, self.i12, self.j, self.nsm,
            self.c1, self.c2, self.d1, self.d2, self.e1, self.e2, self.f1, self.f2):

            i1 = set_blank_if_default(i1, 0.0)
            i2 = set_blank_if_default(i2, 0.0)
            i12 = set_blank_if_default(i12, 0.0)
            j = set_blank_if_default(j, 0.0)

            nsm = set_blank_if_default(nsm, 0.0)
            c1 = set_blank_if_default(c1, 0.0)
            d1 = set_blank_if_default(d1, 0.0)
            e1 = set_blank_if_default(e1, 0.0)
            f1 = set_blank_if_default(f1, 0.0)

            c2 = set_blank_if_default(c2, 0.0)
            d2 = set_blank_if_default(d2, 0.0)
            e2 = set_blank_if_default(e2, 0.0)
            f2 = set_blank_if_default(f2, 0.0)

            #print "  i = ",i
            if i == 0:  # the 1st 2 fields aren't written
                fields += [A, i1, i2, i12, j, nsm, c1, c2, d1, d2,
                                    e1, e2, f1, f2]
            else:
                fields += [so, xxb, A, i1, i2, i12, j, nsm, c1, c2, d1, d2,
                                    e1, e2, f1, f2]
            i += 1
        k1 = set_blank_if_default(self.k1, 1.0)
        k2 = set_blank_if_default(self.k2, 1.0)
        s1 = set_blank_if_default(self.s1, 0.0)
        s2 = set_blank_if_default(self.s2, 0.0)

        nsia = set_blank_if_default(self.nsia, 0.0)
        nsib = set_blank_if_default(self.nsib, self.nsia)

        cwa = set_blank_if_default(self.cwa, 0.0)
        cwb = set_blank_if_default(self.cwb, self.cwa)

        m1a = set_blank_if_default(self.m1a, 0.0)
        m2a = set_blank_if_default(self.m2a, self.m1a)
        m1b = set_blank_if_default(self.m1b, 0.0)
        m2b = set_blank_if_default(self.m2b, self.m1b)

        #print "m1a=%s m2a=%s" %(m1a,m2a)
        #print "m1b=%s m2b=%s" %(m1b,m2b)

        n1a = set_blank_if_default(self.n1a, 0.0)
        n2a = set_blank_if_default(self.n2a, self.n1a)
        n1b = set_blank_if_default(self.n1b, 0.0)
        n2b = set_blank_if_default(self.n2b, self.n1b)
        #print "n1a=%s n2a=%s" %(n1a,n2a)
        #print "n1b=%s n2b=%s" %(n1b,n2b)

        footer = [k1, k2, s1, s2, nsia, nsib, cwa, cwb,
                  m1a, m2a, m1b, m2b, n1a, n2a, n1b, n2b]

        #if footer == [self.k1,None,None,None,None,None,None,None,   None,None,None,None,None,None,None,None]:
        #    footer = []
        fields += footer
        #print fields
        return fields


class PBEAML(IntegratedLineProperty):
    type = 'PBEAML'
    validTypes = {
        "ROD": 1,
        "TUBE": 2,
        "L": 4,
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

    def __init__(self, card=None, data=None):
        IntegratedLineProperty.__init__(self, card, data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)
            self.group = card.field(3, 'MSCBMLO')
            self.Type = card.field(4)

            nDim = self.validTypes[self.Type]
            #nAll = nDim+1

            #nDim = len(self.dim)-1
            dimAll = card.fields(9)

            self.dim = []
            Dim = []
            self.xxb = [0.]
            self.so = ['YES']
            self.nsm = []

            j = 0  # the dimension counter
            n = 0  # there is no SO or XXB for the first section (n=0)
            i = 9  # the index in the card
            #print "dimAll = ",dimAll,nDim,i
            for (ii, dim) in enumerate(dimAll):  ## ii is the counter starting from 9
                if j < nDim:  # the first block, n=0 ???
                    #print "*1"
                    Dim.append(dim)
                    j += 1
                    i += 1
                else:  # other blocks, n>0 ???
                    #print "*2",
                    #print "dim = ",Dim
                    if isinstance(dim, unicode):
                        raise RuntimeError('nsm is a string...nsm=|%s|; fields=%s' % (dim, card.fields()))
                    self.nsm.append(dim)
                    if n > 0:
                        so = card.field(i + 1, 'YES')
                        xxb = card.field(i + 2, 1.0)
                        self.so.append(so)  # dimAll[i+1]
                        self.xxb.append(xxb)  # dimAll[i+2]
                    #j =
                    n += 1
                    i += 2
                    #print "Dim = ",Dim
                    self.dim.append(Dim)
                    Dim = []
                    j = 0
                ###
                #print("i=%s ii=%s j=%s Dim=%s" %(i, ii, j, Dim))
            ###
            if j <= nDim:  # if the last field is blank
                #print "DimB = ",Dim
                self.dim.append(Dim)
                self.nsm.append(0.0)
                if isinstance(self.nsm[0], unicode):
                    msg = 'nsm is a string...nsm=|%s|' % (self.nsm)
                    raise RuntimeError(msg)

            ###
            #print("nsm = %s" %(self.nsm))
            #print self

    def MassPerLength(self):
        """
        mass = L*(Area*rho+nsm)
        mass/L = Area*rho+nsm
        """
        rho = self.Rho()
        massPerLs = []
        for (dim, nsm) in izip(self.dim, self.nsm):
            a = self.areaL(dim)
            try:
                massPerLs.append(a * rho + nsm)
            except:
                msg = "PBEAML a*rho+nsm a=%s rho=%s nsm=%s" % (a, rho, nsm)
                raise RuntimeError(msg)
        massPerL = integrate_positive_line(self.xxb, massPerLs)
        return massPerL

    def Area(self):
        Areas = []
        for dim in self.dim:
            Areas.append(self.areaL(dim))
        A = integrate_line(self.xxb, Areas)
        return A

    #def Mid(self):
    #    return self.mid

    def cross_reference(self, model):
        """
        @warning For structural problems, PBEAML entries must reference a MAT1 material entry
        @warning For heat-transfer problems, the MID must reference a MAT4 or MAT5 material entry.
        @note what happens when there are 2 subcases?
        """
        self.mid = model.Material(self.mid)

    def verify(self, model, iSubcase):
        if model.is_thermal_solution(iSubcase):
            assert self.mid.type in ['MAT4', 'MAT5']
        else:
            assert self.mid.type in ['MAT1']
        ###

    def _J(self):
        j = []
        for dims in self.dim:
            pass
            #print "dims = ",dims
            #IAreaL()
        return j

    def J(self):
        #Js = self._J()
        #j = integrate_positive_line(self.xxb, Js)
        j = None
        return j

    def I11(self):
        #i1 = integrate_positive_line(self.xxb,self.i1)
        i1 = None
        return i1

    def I22(self):
        #i2 = integrate_positive_line(self.xxb,self.i2)
        i2 = None
        return i2

    def I12(self):
        #i12 = integrate_line(self.xxb,self.i12)
        i12 = None
        return i12

    def writeCodeAster(self, iCut=0, iFace=0, iStart=0):  # PBEAML
        msg = ''

        msg2 = 'Cut_%s = geompy.MakeCut(' % (iCut + 1)
        for (xxb, so, dim, nsm) in izip(self.xxb, self.so, self.dim, self.nsm):
            msg += self.CA_Section(iFace, iStart, self.dim)
            msg2 += 'Face_%i, ' % (iFace + 1)
            iFace += 1
            iStart += len(self.dim)

        msg2 = msg2[-2:]
        msg2 += ')\n'

        msg2 += "geompy.addToStudy(Cut_%i,  'Cut_%i')\n" % (
            iCut + 1, iCut + 1)
        iCut += 1
        return (msg + msg2, iCut, iFace, iStart)

    def rawFields(self):
        fields = ['PBEAML', self.pid, self.Mid(), self.group,
            self.Type, None, None, None, None]
        #print "self.nsm = ",self.nsm
        #print "xxb=%s so=%s dim=%s nsm=%s" %(self.xxb,self.so,self.dim,self.nsm)
        for (i, xxb, so, dim, nsm) in izip(count(), self.xxb, self.so, self.dim, self.nsm):
            if i == 0:
                fields += dim + [nsm]
            else:
                fields += [xxb, so] + dim + [nsm]
            ###
        ###
        #print self.printCard(fields)
        #raise NotImplementedError('verify PBEAML...')
        return fields

    def reprFields(self):
        group = set_blank_if_default(self.group, 'MSCBMLO')
        fields = self.rawFields()
        fields[3] = group
        return fields


class PBEAM3(LineProperty):  # not done, cleanup
    type = 'PBEAM3'

    def __init__(self, card=None, data=None):
        LineProperty.__init__(self, card, data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)

            self.A = card.field(3)
            self.iz = card.field(4)
            self.iy = card.field(5)
            self.iyz = card.field(6, 0.0)
            self.j = card.field(7, self.iy + self.iz)
            self.nsm = card.field(8, 0.0)

            self.cy = card.field(9)
            self.cz = card.field(10)
            self.dy = card.field(11)
            self.dz = card.field(12)
            self.ey = card.field(13)
            self.dz = card.field(14)
            self.fy = card.field(15)
            self.fz = card.field(16)

            # more...
        ###
        else:
            raise NotImplementedError('not implemented...')
        ###

    def Nsm(self):
        """@warning nsm field not supported fully on PBEAM3 card"""
        return self.nsm

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def reprFields(self):
        raise NotImplementedError('not done...')
        fields = ['PBEAM3', self.pid, self.Mid(), ]  # other
        return fields


class PBEND(LineProperty):
    type = 'PBEND'

    def __init__(self, card=None, data=None):
        LineProperty.__init__(self, card, data)
        if card:
            self.pid = card.field(1)
            self.mid = card.field(2)

            value3 = card.field(3)
            if isinstance(value3, float):
                self.beamType = 1
                ## Area of the beam cross section
                self.A = card.field(3)
                ## Area moments of inertia in planes 1 and 2.
                self.i1 = card.field(4)
                self.i2 = card.field(5)
                ## Torsional stiffness
                self.j = card.field(6)

                # line2
                ## The r,z locations from the geometric centroid for stress
                ## data recovery.
                self.c1 = card.field(9)
                self.c2 = card.field(10)
                self.d1 = card.field(11)
                self.d2 = card.field(12)
                self.e1 = card.field(13)
                self.e2 = card.field(14)
                self.f1 = card.field(15)
                self.f2 = card.field(16)

                # line 3
                ## Shear stiffness factor K in K*A*G for plane 1 and plane 2.
                self.k1 = card.field(17)
                self.k2 = card.field(18)
                ## Radial offset of the neutral axis from the geometric
                ## centroid, positive is toward the center of curvature
                self.deltaN = card.field(22)

            elif isinstance(value3, int):
                self.beamType = 2
                ## Flag selecting the flexibility and stress intensification
                ## factors. See Remark 3. (Integer = 1, 2, or 3)
                self.fsi = card.field(3)
                ## Mean cross-sectional radius of the curved pipe
                self.rm = card.field(4)
                ## Wall thickness of the curved pipe
                self.t = card.field(5)
                ## Internal pressure
                self.p = card.field(6)

                # line3
                self.nsm = card.field(19)
                self.rc = card.field(20)
                self.zc = card.field(21)
            else:
                raise RuntimeError('Area/FSI on CBEND must be defined...')
            ## Bend radius of the line of centroids
            self.rb = card.field(7)
            ## Arc angle of element  (optional)
            self.thetab = card.field(8)

            ## Nonstructural mass per unit length.
            self.nsm = card.field(19)
            ## Radial offset of the geometric centroid from points GA and GB
            self.rc = card.field(20)
            ## Offset of the geometric centroid in a direction perpendicular to the
            ## plane of points GA and GB and vector v
            self.zc = card.field(21)
        ###
        else:
            raise NotImplementedError('PBEND')
        ###

    #def Nsm(self):
        #"""@warning nsm field not supported fully on PBEND card"""
        #raise RuntimeError(self.nsm[0])
        #return self.nsm

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def reprFields(self):
        fields = ['PBEND', self.pid, self.Mid(), ]  # other
        if self.beamType == 1:
            fields += [self.A, self.i1, self.i2, self.j, self.rb, self.thetab,
                       self.c1, self.c2, self.d1, self.d2, self.e1, self.e2, self.f1, self.f2,
                       self.k1, self.k2, self.nsm, self.rc, self.zc, self.deltaN]
        elif self.beamType == 2:
            fields += [self.fsi, self.rm, self.t, self.p, self.rb, self.thetab,
                       None, None, self.nsm, self.rc, self.zc]
        else:
            raise ValueError('only beamType=1 and 2 supported')
        return fields
