"""
All bar properties are defined in this file.  This includes:
 *   PBAR
 *   PBARL
 *   PBEAM3
 *   PBEND
 *   PBRSECT

All bars are LineProperty objects.
Multi-segment beams are IntegratedLineProperty objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from six import integer_types, string_types
from numpy import pi, array
import numpy as np

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, double_or_blank, string, string_or_blank,
    blank, integer_or_double, #integer_or_blank,
)
from pyNastran.utils.mathematics import integrate_unit_line, integrate_positive_unit_line
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.utils import to_fields
from pyNastran.utils import float_types


def Iyy_beam(b, h):
    return 1 / 12. * b * h ** 3


def I_beam(b, h):
    f = 1 / 12. * b * h
    Iyy = f * h * h  # 1/12.*b*h**3
    Izz = f * b * b  # 1/12.*h*b**3
    Iyz = 0.
    return (Iyy, Izz, Iyz)


def I_beam_offset(b, h, y, z):
    A = b * h
    f = 1. / 12. * A

    Iyy = f * h * h  # 1/12.*b*h**3
    Izz = f * b * b  # 1/12.*h*b**3
    Iyz = 0.

    Iyy += A * y * y
    Izz += A * z * z
    Iyz += A * y * z
    return (Iyy, Izz, Iyz)


def get_inertia_rectangular(sections):
    """
    Calculates the moment of inertia for a section about the CG.

    Parameters
    ----------
    sections : [[b,h,y,z]_1,...]
        [[b,h,y,z]_1,...] y,z is the centroid
        (x in the direction of the beam, y right, z up)

    Returns
    -------
    interia_parameters : List[Area, Iyy, Izz, Iyz]
        the inertia parameters

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

    xcg = Ax / A
    ycg = Ay / A
    Axx = 0.
    Ayy = 0.
    Axy = 0.
    for (i, section) in enumerate(sections):
        (b, h, x, y) = section
        #A = b*h
        #As.append(A)
        Axx += As[i] * (x - xcg) ** 2
        Ayy += As[i] * (y - ycg) ** 2
        Axy += As[i] * (x - xcg) * (y - ycg)
    Ixx = Axx / A
    Iyy = Ayy / A
    Ixy = Axy / A
    return (A, Ixx, Iyy, Ixy)


class LineProperty(Property):
    def __init__(self):
        self.beam_type = None
        self.dim = None
        self.A = None
        self.i1 = None
        self.i2 = None
        self.i12 = None
        self.j = None
        self.nsm = None
        Property.__init__(self)
        self.mid_ref = None

    #def D_bending(self):
        #pass

    #def D_axial(self):
        #pass

    #def D_thermal(self):
        #pass

    #def D_shear(self):
        #pass

    def Rho(self):
        return self.mid_ref.rho

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
        return self.mid_ref.E

    def G(self):
        return self.mid_ref.G

    def Nu(self):
        return self.mid_ref.nu

    def IAreaL(self, dim):
        if self.beam_type == 'ROD':
            R = dim[0]
            A = pi * R ** 2
            Iyy = A * R ** 2 / 4.
            Izz = Iyy
            Iyz = 0.
        elif self.beam_type == 'TUBE':
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
        elif self.beam_type == 'TUBE2':
            R1 = dim[0]
            t = dim[1]
            R2 = R1 - t
            A1 = pi * R1 ** 2
            Iyy1 = A1 * R1 ** 2 / 4.
            A2 = pi * R2 ** 2
            Iyy2 = A2 * R2 ** 2 / 4.
            A = A1 - A2
            Iyy = Iyy1 - Iyy2
            Izz = Iyy
            Iyz = 0.

        elif self.beam_type == 'I':
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
            #y3 = -dim[0] / 2. + h3
            sections.append([w3, h3, 0., y1])

            h2 = dim[0] - h1 - h3
            w2 = dim[3]  # d1
            sections.append([w2, h2, 0., 0.])

            (A, Iyy, Izz, Iyz) = get_inertia_rectangular(sections)
            assert Iyz == 0.

        elif self.beam_type == 'BAR':
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
            msg = 'beam_type=%s is not supported for %s class...' % (
                self.beam_type, self.type)
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
        if self.beam_type == 'ROD':
            R = dim[0]
            A = pi * R ** 2
            I1 = A * R ** 2 / 4.
            I2 = I1
            I12 = 0.
        elif self.beam_type == 'BAR':
            I1 = 1 / 12. * dim[0] * dim[1] ** 3
            I2 = 1 / 12. * dim[1] * dim[0] ** 3
            I12 = 0.
        else:
            msg = 'I1_I2_I12; beam_type=%s is not supported for %s class...' % (
                self.beam_type, self.type)
            raise NotImplementedError(msg)
        return(I1, I2, I12)

def _bar_areaL(class_name, beam_type, dim, prop):
    """
    Area(x) method for the PBARL and PBEAML classes (pronounced **Area-L**)

    Parameters
    ----------
    dim : List[float]
        a list of the dimensions associated with **beam_type**

    Returns
    -------
    Area : float
        Area of the given cross section defined by **self.beam_type**

    Notes
    -----
    internal method
    """
    if beam_type == 'ROD':
        # This is a circle if you couldn't tell...
        #   __C__
        #  /     \
        # |       |
        # F       D
        # |       |
        #  \     /
        #   \_E_/
        #
        # Radius = dim1
        A = pi * dim[0] ** 2
    elif beam_type == 'TUBE':
        r_outer, r_inner = dim
        A = pi * (r_outer ** 2 - r_inner ** 2)
        assert r_outer > r_inner, 'TUBE; r_outer=%s r_inner=%s' % (r_outer, r_inner)
    elif beam_type == 'TUBE2':
        r_outer, thick = dim
        r_inner = r_outer - thick
        assert r_outer > r_inner, 'TUBE2; r_outer=%s r_inner=%s' % (r_outer, r_inner)
        A = pi * (r_outer ** 2 - r_inner ** 2)

    #I = (DIM3*DIM6)+(DIM2*DIM5) + ((DIM1-(DIM5+DIM6))*DIM4)
    elif beam_type == 'I':
        # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
        #  <----d3---->
        #
        #  1----------2      ^
        #  |    A     |  d6  |
        # 12-11---4---3      |
        #     |   |          |
        #     |   |          |  d1
        #     | B | <--- d4  |
        #     |   |          |
        #  9-10---5---6      |
        #  |    C     |  d5  |
        #  8----------7      v
        #
        #  <----d2---->
        #
        #
        # h1 = hA = d6
        # h2 = hB = d1 - d6 - d5
        # h3 = hC = d5
        # w1 = wA = d3
        # w2 = wB = d4
        # w3 = wC = d5
        # A = A1 + A2 + A3
        dim1, dim2, dim3, dim4, dim5, dim6 = dim

        h = dim1
        a = dim2
        b = dim3
        tw = dim4
        ta = dim5
        tb = dim6
        hw = h - (ta + tb)
        hf = h - 0.5 * (ta + tb)
        assert hw > 0, 'hw=%s' % hw
        assert hf > 0, 'hf=%s' % hf

        A = ta * a + hw * tw + b * tb
        #yc = (0.5 * hw * (hw + ta)*tw + hf*tb*b) / A
        #ys = tb * hf * b**3/(tb * b**3 + ta * a**a)
        #yna = yc - ys
        #i1 = (
        #    1/12 * (h*tb**3 + a*ta**3 + tw*hw**3)
        #    + (hf-yc)**2*b*tb+yc**2*a*ta
        #    + (yc - 0.5*(hw+ta))**2*hw*tw
        #)
        #i2 = (b**3 * tb + ta * a**3 + hw*tw**3)/12.
        #i12 = 0.
        #j = 1/3 * (tb**3*b + ta**3*a + tw**3*hf)

        assert dim1 > dim5+dim6, 'I; required: dim1 > dim5+dim6; dim1=%s dim5=%s dim6=%s\n%s' % (dim1, dim5, dim6, prop)
        #THE SUM OF THE FLANGE THICKNESSES,DIM5 + DIM6,
        #CAN NOT BE GREATER THAT THE WEB HEIGHT, DIM1.

    #I1 = (DIM2*DIM3)+((DIM4-DIM3)*(DIM1+DIM2))
    elif beam_type == 'I1':
        #
        #  d1/2   d2  d1/2
        #  <---><---><--->
        #
        #  1--------------2     ^
        #  |      A       |     |
        # 12---11---4-----3     |
        #       |   |        ^  |
        #       |   |        |  |  d4
        #       | B |     d3 |  |
        #       |   |        v  |
        #  9----10---5----6     |
        #  |      C       |     |
        #  8--------------7     v
        #
        #  <----d2---->
        #
        #
        # h1 = hA = (d4 - d3) / 2
        # h2 = hB = d3
        #
        # w1 = wA = d1 + d2
        # w2 = d2
        #
        # A3 = A1
        # A = A1 + A2 + A3
        w1 = dim[0] + dim[1]
        h1 = (dim[3] - dim[2]) / 2.

        w2 = dim[1]
        h2 = dim[2]

        A = 2. * (h1 * w1) + h2 * w2

    elif beam_type == 'L':
        # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
        #
        #  D4
        # F---C      ^
        # |   |      |
        # | 2 |      |
        # |   |      | D2
        # +---+---+  |
        # |1   D3 |  |
        # E-------D  v
        #
        # <------> D1
        #
        (dim1, dim2, dim3, dim4) = dim
        t1 = dim3
        t2 = dim4
        b = dim1 - 0.5 * t2
        h = dim2 - 0.5 * t1
        h2 = dim2 - t1
        #b1 = dim1 - t2
        A = (b + 0.5 * t2) * t1 + h2 * t2
        #yc = t2*h2 * (h2 + t1) / (2 * A)
        #zc = t1*b1 * (b1 + t2) / (2 * A)
        #i1 = (
        #    t1 ** 3 * (b + 0.5 * t2) / 12. +
        #    t1 * (b + 0.5 * t2) * yc ** 2 +
        #    t2 * h ** 3 / 12 +
        #    h2 * t2 * (0.5 * (h2 + t1) - yc) ** 2
        #)
        #i2 = (
        #    t2 ** 3 * h2 / 12 +
        #    t1*(b + 0.5 * t2) ** 3 / 12 +
        #    t1*(b + 0.5 * t2) * (0.5 * b1 - zc) ** 2  # zc is z2 in the docs...
        #)
        #i12 = (
        #    zc * yc * t1 * t2
        #    - b1 * t1 * yc * (0.5 * (b1 + t2) - zc)
        #    - h2 * t2 * zc * (0.5 * (h2 + t1) - yc)
        #)
        #j = 1/3 * (t1**3*b + t2**3 * h)
        #k1 = h2 * t2 / A
        #k2 = b1 * t1 / A
        #yna = yc
        #zna = zc

    #CHAN = 2*(DIM1*DIM4) + (DIM2-(2*DIM4))*DIM3
    elif beam_type == 'CHAN':
        #
        # +-+----+   ^  ^
        # | |    |   |  | d4
        # | +----+   |  v
        # | |        |
        # | |<-- d3  |
        # | |        | d2
        # | |        |
        # | +----+   |
        # | |    |   |
        # +-+----+   v
        #
        # <--d1-->
        #

        dim1, dim2, dim3, dim4 = dim
        tweb = dim3
        tf = dim4

        bf = dim1
        d = dim2 - 2 * tf

        # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
        tw = dim3
        tf = dim4
        #b = dim1 - 0.5 * tw
        h = dim2 - tf
        bf = dim1 - tw
        #hw = dim2 - 2 * tf

        #A = 2 * tf * bf + (h + tf) * tweb  # I think tt is tf...
        #zc = bf * tf * (bf + tw) / A
        #zs = b**2 * tf / (2 * b * tw + h * tf / 3)
        #i1 = (
            #h ** 2 * tf * bf / 2
            #+ bf * tf ** 3 / 6
            #+ (h + tf) ** 3 * tw / 12
        #)
        #i2 = (
            #(h + tf) * tw**3/12
            #+ bf**3 * tf / 6
            #+ 0.5 * (bf + tw) ** 2 * bf * tf
            #- zc ** 2 * A
        #)
        #j = (2 * b * tf **3 + h * tw ** 3) / 3

        # warping coefficient for the cross section relative to the shear center
        #iw = tf * b**3 * h**2 / 12 * (2 * tw * h + 3 * tf * b) / (tw * h + 6 * tf * b)
        #k1 = tw * hw / A
        #k2 = 2 * tf * bf / A
        #zna = zc + zs
        hweb = dim2 - 2 * tf
        A1 = 2 * tf * bf + hweb * tweb

        #A2 = 2. * bf * tf + (dim2 - 2. * dim4) * tweb
        A2 = 2. * tf * bf + (dim2 - 2. * dim4) * tweb
        A0 = A1
        assert np.allclose(A0, A2), 'A0=%s A1=%s A2=%s' % (A0, A2, A2)
        assert np.allclose(A1, A2), 'A0=%s A1=%s A2=%s' % (A0, A1, A2)
        A = A1

    #CHAN1 = DIM2*DIM3 + (DIM4-DIM3)*(DIM1+DIM2)
    elif beam_type == 'CHAN1':
        h2 = dim[2]
        w2 = dim[1]

        h1 = dim[3] - h2
        w1 = dim[0] + w2
        A = h1 * w1 + h2 * w2
        dim1, dim2, dim3, dim4 = dim
        A2 = dim2 * dim3 + (dim4 - dim3) * (dim1 + dim2)
        assert np.allclose(A, A2), 'A=%s A2=%s' % (A, A2)

    #CHAN2 = 2*(DIM1*DIM3)+((DIM4-(2*DIM1))*DIM2)
    elif beam_type == 'CHAN2':
        #  d1        d1
        # <-->      <-->
        # +--+      +--+        ^
        # |  |      |  |        |
        # |  |      |  |        | d3
        # +--+------+--+  ^     |
        # |            |  | d2  |
        # +------------+  v     v
        #
        # <-----d4----->
        #
        dim1, dim2, dim3, dim4 = dim
        hupper = dim3 - dim2
        hlower = dim2
        wlower = dim4
        wupper = dim1
        A = 2 * hupper * wupper + hlower * wlower

        h2 = dim2
        w2 = dim4

        h1 = dim3 - h2
        w1 = dim1 * 2
        A2 = h1 * w1 + h2 * w2
        A3 = 2 * dim1 * dim3 + (dim4 - 2 * dim1) * dim2
        assert np.allclose(A, A2), 'A=%s A2=%s A3=%s' % (A, A2, A3)
        assert np.allclose(A, A3), 'A=%s A2=%s A3=%s' % (A, A2, A3)

    #T = (DIM1*DIM3)+((DIM2-DIM3)*DIM4)
    elif beam_type == 'T':
        # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
        #           ^ y
        #           |
        #           |
        # ^    +-----------+ ^
        # | d3 |           | |  ----> z
        # v    +---+  +----+ |
        #          |  |      |
        #          |  |      | d2
        #          |  |      |
        #          |  |      |
        #          +--+      v
        #
        #          <--> d4
        #
        dim1, dim2, dim3, dim4 = dim
        d = dim1
        tf = dim3
        tw = dim4
        h = dim2 - 0.5 * tf
        hw = dim2 - tf
        h1 = dim[2]
        w1 = dim[0]

        A = d*tf + hw*tw
        #yna = hw*tw*(hw+tf)/(2*A)
        #i1 = (
        #    (d*tf**3 + tw*hw**3) / 12. +
        #    hw*tw*(yna + 0.5 * (hw +tf))**2 + d*tf*yna**2
        #)
        #i2 = (tf*d**3 + hw*tw**3)/12.
        #i12 = 0.
        #j = 1/3 * (tf**3*d + tw**3 * h)

    #T1 = (DIM1*DIM3)+(DIM2*DIM4)
    elif beam_type == 'T1':
        h1 = dim[0]
        w1 = dim[2]

        h2 = dim[3]
        w2 = dim[1]
        A = h1 * w1 + h2 * w2

    #T2 = (DIM1*DIM3)+((DIM2-DIM3)*DIM4)
    elif beam_type == 'T2':
        h1 = dim[3]
        w1 = dim[1]

        h2 = h1 - dim[2]
        w2 = dim[0]
        A = h1 * w1 + h2 * w2

    #BOX = 2*(DIM1*DIM3)+2*((DIM2-(2*DIM3))*DIM4)
    elif beam_type == 'BOX':
        #
        # +----------+ ^     ^
        # |          | | d3  |
        # |  +----+  | v     |
        # |  |    |  |       |  d2
        # |  +----+  |       |
        # |          |       |
        # +----------+       v
        #          <-> d4
        # <----d1---->
        #
        dim1, dim2, dim3, dim4 = dim
        #h1 = dim3
        #w1 = dim1

        #h2 = dim2 - 2 * h1
        #w2 = dim4
        assert dim1 > 2.*dim4, 'BOX; required: dim1 > 2*dim4; dim1=%s dim4=%s\n%s' % (dim1, dim4, prop)
        assert dim2 > 2.*dim3, 'BOX; required: dim2 > 2*dim3; dim2=%s dim3=%s\n%s' % (dim2, dim3, prop)
        #A = 2 * (h1 * w1 + h2 * w2)

        # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
        b = dim1
        h = dim2
        t1 = dim3
        t2 = dim4
        bi = b - 2 * t2
        hi = h - 2 * t1
        A = b * h - bi * hi
        #i1 = (b*h**3 * bi*hi**3) / 12.
        #i2 = (h*b**3 * hi*bi**3) / 12.
        #i12 = 0.
        #j = (2*t2*t1*(b-t2)**2*(h-t1)**2)/(b*t2+h*t1-t2**2-t1**2)

    #BOX1 = (DIM2*DIM6)+(DIM2*DIM5)+((DIM1-DIM5-DIM6)*DIM3)+((DIM1-DIM5-DIM6)*DIM4)
    elif beam_type == 'BOX1':
        h1 = dim[2]  # top
        w1 = dim[0]

        h2 = dim[3]  # btm
        A1 = (h1 + h2) * w1

        h3 = dim[1] - h1 - h2  # left
        w3 = dim[5]

        w4 = dim[4]  # right
        A2 = h3 * (w3 + w4)
        A = A1 + A2

    #BAR   = DIM1*DIM2
    elif beam_type == 'BAR':
        #per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
        # <------> D1
        #
        # F------C  ^
        # |      |  |
        # |      |  | D2
        # |      |  |
        # E------D  v
        dim1, dim2 = dim
        b = dim1
        h = dim2
        A = b * h
        #i1 = b*h**3/12.
        #i2 = b**3*h/12.
        #i12 = 0.
        #J = b*h**3*(1/3. - 0.21*h/b*(1-h**4/(12*b**4)))

    #CROSS = (DIM2*DIM3)+2*((0.5*DIM1)*DIM4)
    elif beam_type == 'CROSS':
        h1 = dim[2]
        w1 = dim[1]

        h2 = dim[3]
        w2 = dim[0]
        A = h1 * w1 + h2 * w2

    #H = 2*((0.5*DIM2)*DIM3)+(DIM1*DIM4)
    elif beam_type == 'H':
        h1 = dim[2]
        w1 = dim[1]

        h2 = dim[3]
        w2 = dim[0]
        A = h1 * w1 + h2 * w2

    #Z = (DIM2*DIM3)+((DIM4-DIM3)*(DIM1+DIM2))
    elif beam_type == 'Z':
        #dim1, dim2, dim3, dim4 = dim
        h2 = dim[2]
        w2 = dim[1]

        h1 = dim[3] - h2
        w1 = dim[0]
        A = h1 * w1 + h2 * w2

    #HEXA = ((DIM2-(2*DIM1))*DIM3)+(DIM3*DIM1)
    elif beam_type == 'HEXA':
        #     _______
        #   /        \     ^
        #  /          \    |
        # *            *   |  d3
        #  \          /    |
        #   \        /     |
        #    \______/      v
        #           |d1|
        # <-----d2---->
        #
        dim1, dim2, dim3 = dim
        hbox = dim3
        wbox = dim2
        wtri = dim1
        A = hbox * wbox - 2. * wtri * hbox
        #print('hbox=%s wbox=%s hbox*wbox=%s 2*wtri*hbox=%s A=%s' % (
            #hbox, wbox, hbox*wbox, 2*wtri*hbox, A))

    #HAT = (DIM2*DIM3)+2*((DIM1-DIM2)*DIM2)+2*(DIM2*DIM4)
    elif beam_type == 'HAT':
        #
        #        <--------d3------->
        #
        #        +-----------------+              ^
        #   d4   |        A        |   d4         |
        # <----> +-d2-+-------+-d2-+ <---->       |
        #        | B  |       |  B |              | d1
        # +------+----+       +----+------+       |
        # |     C     |       |     C     | t=d2  |
        # +-----------+       +-----------+       v
        dim1, dim2, dim3, dim4 = dim
        assert dim3 > 2.*dim2, 'HAT; required: dim3 > 2*dim2; dim2=%s dim3=%s; delta=%s\n%s' % (dim2, dim3, dim3-2*dim2, prop)
        #DIM3, CAN NOT BE LESS THAN THE SUM OF FLANGE
            #THICKNESSES, 2*DIM2
        t = dim[1]
        wa = dim[2]
        hb = dim[0] - 2. * t
        wc = dim[3] + t
        A = wa * t + (2. * wc * t) + (2. * hb * t)

    #HAT1 = (DIM1*DIM5)+(DIM3*DIM4)+((DIM1-DIM3)*DIM4)+2*((DIM2-DIM5-DIM4)*DIM4)
    elif beam_type == 'HAT1':
        # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
        w = dim[3]

        h0 = dim[4]         # btm bar
        w0 = dim[0] / 2.

        h2 = dim[1] - h0 - 2 * w  # vertical bar
        w2 = w

        h3 = w              # top bar
        w3 = dim[2] / 2.

        dim1, dim2, dim3, dim4, dim5 = dim
        assert dim2 > dim4+dim5, 'HAT1; required: dim2 > dim4+dim5; dim2=%s dim4=%s; dim5=%s\n%s' % (dim2, dim4, dim5, prop)
        #*DIM4+DIM5, CAN NOT BE LARGER THAN THE HEIGHT OF
            #THE HAT, DIM2.


        h1 = w              # upper, horizontal lower bar (see HAT)
        w1 = w0 - w3
        A = 2. * (h0 * w0 + h1 * w1 + h2 * w2 + h3 * w3)
    #DBOX = ((DIM2*DIM3)-((DIM2-DIM7-DIM8)*(DIM3-((0.5*DIM5)+DIM4)))) +
    #       (((DIM1-DIM3)*DIM2)-((DIM2-(DIM9+DIM10))*(DIM1-DIM3-(0.5*DIM5)-DIM6)))
    elif beam_type == 'DBOX':
        #
        #  |--2------5----
        #  |     |       |
        #  1     3       6
        #  |     |       |
        #  |--4--|---7---|
        #

        #0,1,2,6,11
        #1,2,3,7,12

        htotal = dim[1]
        wtotal = dim[0]

        h2 = dim[6]
        w2 = dim[3]

        h4 = dim[7]
        w4 = w2

        h1 = htotal - h2 - h4
        w1 = dim[3]

        h5 = dim[8]
        w5 = wtotal - w2

        h7 = dim[9]
        w7 = w5

        h6 = htotal - h5 - h7
        w6 = dim[5]

        h3 = (h1 + h6) / 2.
        w3 = dim[4]

        A = (h1 * w1 + h2 * w2 + h3 * w3 + h4 * w4 +
             h5 * w5 + h6 * w6 + h7 * w7)
    else:
        msg = 'areaL; beam_type=%s is not supported for %s class...' % (
            beam_type, class_name)
        raise NotImplementedError(msg)
    assert A > 0, 'beam_type=%r dim=%r A=%s\n%s' % (beam_type, dim, A, prop)
    #A = 1.
    return A

class IntegratedLineProperty(LineProperty):
    def __init__(self):
        self.xxb = None
        self.A = None
        self.j = None
        self.i1 = None
        self.i2 = None
        self.i12 = None
        LineProperty.__init__(self)

    def Area(self):
        A = integrate_positive_unit_line(self.xxb, self.A)
        return A

    def J(self):
        J = integrate_positive_unit_line(self.xxb, self.j)
        return J

    def I11(self):
        i1 = integrate_positive_unit_line(self.xxb, self.i1)
        return i1

    def I22(self):
        i2 = integrate_positive_unit_line(self.xxb, self.i2)
        return i2

    def I12(self):
        i12 = integrate_unit_line(self.xxb, self.i12)
        return i12

    def Nsm(self):
        #print("xxb = ",self.xxb)
        #print("nsm = ",self.nsm)
        nsm = integrate_positive_unit_line(self.xxb, self.nsm)
        return nsm


class PBAR(LineProperty):
    """
    Defines the properties of a simple beam element (CBAR entry).

    +------+-----+-----+-----+----+----+----+-----+-----+
    |   1  |  2  |  3  |  4  |  5 |  6 |  7 |  8  |  9  |
    +======+=====+=====+=====+====+====+====+=====+=====+
    | PBAR | PID | MID |  A  | I1 | I2 | J  | NSM |     |
    +------+-----+-----+-----+----+----+----+-----+-----+
    |      | C1  | C2  | D1  | D2 | E1 | E2 | F1  | F2  |
    +------+-----+-----+-----+----+----+----+-----+-----+
    |      | K1  | K2  | I12 |    |    |    |     |     |
    +------+-----+-----+-----+----+----+----+-----+-----+

    .. todo::
        support solution 600 default
        do a check for mid -> MAT1      for structural
        do a check for mid -> MAT4/MAT5 for thermal
    """
    type = 'PBAR'
    pname_fid_map = {
        # 1-based
        4 : 'A', 'A' : 'A',
        5 : 'i1', 'I1' : 'i1',
        6 : 'i2', 'I2' : 'i2',
        7 : 'j', 'J' : 'j',
        10 : 'c1',
        11 : 'c2',
        12 : 'd1',
        13 : 'd2',
        14 : 'e1',
        15 : 'e2',
        16 : 'f1',
        17 : 'f2',
        18 : 'k1',
        19 : 'k1',
        20 : 'i12', 'I12' : 'i12',
    }

    def __init__(self, pid, mid, A=0., i1=0., i2=0., i12=0., j=0., nsm=0.,
                 c1=0., c2=0., d1=0., d2=0., e1=0., e2=0., f1=0., f2=0.,
                 k1=1.e8, k2=1.e8, comment=''):
        """
        Creates a PBAR card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        area : float
            area
        i1, i2, i12, j : float
            moments of inertia
        nsm : float; default=0.
            nonstructural mass per unit length
        c1/c2, d1/d2, e1/e2, f1/f2 : float
           the y/z locations of the stress recovery points
           c1 - point C.y
           c2 - point C.z

        k1 / k2 : float; default=1.e8
            Shear stiffness factor K in K*A*G for plane 1/2.
        comment : str; default=''
            a comment for the card
        """
        LineProperty.__init__(self)
        if comment:
            self.comment = comment
        #: property ID -> use Pid()
        self.pid = pid
        #: material ID -> use Mid()
        self.mid = mid
        #: Area -> use Area()
        self.A = A
        #: I1 -> use I1()
        self.i1 = i1
        #: I2 -> use I2()
        self.i2 = i2

        #: I12 -> use I12()
        self.i12 = i12

        #: Polar Moment of Inertia J -> use J()
        #: default=1/2(I1+I2) for SOL=600, otherwise 0.0
        #: .. todo:: support SOL 600 default
        self.j = j

        #: nonstructral mass -> use Nsm()
        self.nsm = nsm

        self.c1 = c1
        self.c2 = c2
        self.d1 = d1
        self.d2 = d2
        self.e1 = e1
        self.e2 = e2
        self.f1 = f1
        self.f2 = f2

        # K1/K2 must be blank
        #: default=infinite; assume 1e8
        self.k1 = k1
        #: default=infinite; assume 1e8
        self.k2 = k2
        self.mid_ref = None

    def validate(self):
        if self.i1 < 0.:
            raise ValueError('I1=%r must be greater than or equal to 0.0' % self.i1)
        if self.i2 < 0.:
            raise ValueError('I2=%r must be greater than or equal to 0.0' % self.i2)
        if self.j < 0.:
            raise ValueError('J=%r must be greater than or equal to 0.0' % self.j)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBAR card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        A = double_or_blank(card, 3, 'A', 0.0)
        i1 = double_or_blank(card, 4, 'I1', 0.0)
        i2 = double_or_blank(card, 5, 'I2', 0.0)

        j = double_or_blank(card, 6, 'J', 0.0)
        nsm = double_or_blank(card, 7, 'nsm', 0.0)

        c1 = double_or_blank(card, 9, 'C1', 0.0)
        c2 = double_or_blank(card, 10, 'C2', 0.0)
        d1 = double_or_blank(card, 11, 'D1', 0.0)
        d2 = double_or_blank(card, 12, 'D2', 0.0)
        e1 = double_or_blank(card, 13, 'E1', 0.0)
        e2 = double_or_blank(card, 14, 'E2', 0.0)
        f1 = double_or_blank(card, 15, 'F1', 0.0)
        f2 = double_or_blank(card, 16, 'F2', 0.0)

        i12 = double_or_blank(card, 19, 'I12', 0.0)

        if A == 0.0:
            k1 = blank(card, 17, 'K1')
            k2 = blank(card, 18, 'K2')
        elif i12 != 0.0:
            # K1 / K2 are ignored
            k1 = None
            k2 = None
        else:
            #: default=infinite; assume 1e8
            k1 = double_or_blank(card, 17, 'K1', 1e8)
            #: default=infinite; assume 1e8
            k2 = double_or_blank(card, 18, 'K2', 1e8)

        assert len(card) <= 20, 'len(PBAR card) = %i\ncard=%s' % (len(card), card)
        return PBAR(pid, mid, A, i1, i2, i12, j, nsm,
                    c1, c2, d1, d2, e1, e2,
                    f1, f2, k1, k2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        pid = data[0]
        mid = data[1]
        A = data[2]
        i1 = data[3]
        i2 = data[4]
        j = data[5]

        nsm = data[6]
        #self.fe  = data[7] #: .. todo:: not documented....
        c1 = data[8]
        c2 = data[9]
        d1 = data[10]
        d2 = data[11]
        e1 = data[12]
        e2 = data[13]
        f1 = data[14]
        f2 = data[15]
        k1 = data[16]
        k2 = data[17]
        i12 = data[18]
        if k1 == 0.:
            k1 = None
        if k2 == 0.:
            k2 = None
        return PBAR(pid, mid, A, i1, i2, i12, j, nsm,
                    c1, c2, d1, d2, e1, e2,
                    f1, f2, k1, k2, comment=comment)

    def _verify(self, xref):
        pid = self.pid
        mid = self.Mid()
        A = self.Area()
        J = self.J()
        #c = self.c
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(A, float), 'pid=%r' % A
        assert isinstance(J, float), 'cid=%r' % J
        #assert isinstance(c, float), 'c=%r' % c
        if xref:
            nsm = self.Nsm()
            mpa = self.MassPerLength()
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
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PBAR mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def uncross_reference(self):
        self.mid = self.Mid()
        self.mid_ref = None

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

    def raw_fields(self):
        list_fields = ['PBAR', self.pid, self.Mid(), self.A, self.i1, self.i2,
                       self.j, self.nsm, None, self.c1, self.c2, self.d1, self.d2,
                       self.e1, self.e2, self.f1, self.f2, self.k1, self.k2,
                       self.i12]
        return list_fields

    def repr_fields(self):
        #A  = set_blank_if_default(self.A,0.0)
        i1 = set_blank_if_default(self.i1, 0.0)
        i2 = set_blank_if_default(self.i2, 0.0)
        i12 = set_blank_if_default(self.i12, 0.0)
        j = set_blank_if_default(self.j, 0.0)
        nsm = set_blank_if_default(self.nsm, 0.0)

        c1 = set_blank_if_default(self.c1, 0.0)
        c2 = set_blank_if_default(self.c2, 0.0)

        d1 = set_blank_if_default(self.d1, 0.0)
        d2 = set_blank_if_default(self.d2, 0.0)

        e1 = set_blank_if_default(self.e1, 0.0)
        e2 = set_blank_if_default(self.e2, 0.0)

        f1 = set_blank_if_default(self.f1, 0.0)
        f2 = set_blank_if_default(self.f2, 0.0)

        k1 = set_blank_if_default(self.k1, 1e8)
        k2 = set_blank_if_default(self.k2, 1e8)

        list_fields = ['PBAR', self.pid, self.Mid(), self.A, i1, i2, j, nsm,
                       None, c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, i12]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBARL(LineProperty):
    """
    .. todo:: doesnt support user-defined types

    +-------+------+------+-------+------+------+------+------+------+
    |   1   |   2  |   3  |   4   |   5  |   6  |   7  |   8  |   9  |
    +=======+======+======+=======+======+======+======+======+======+
    | PBARL | PID  | MID  | GROUP | TYPE |      |      |      |      |
    +-------+------+------+-------+------+------+------+------+------+
    |       | DIM1 | DIM2 | DIM3  | DIM4 | DIM5 | DIM6 | DIM7 | DIM8 |
    +-------+------+------+-------+------+------+------+------+------+
    |       | DIM9 | etc. |  NSM  |      |      |      |      |      |
    +-------+------+------+-------+------+------+------+------+------+
    """
    type = 'PBARL'
    valid_types = {
        "ROD": 1,
        "TUBE": 2,
        "TUBE2": 2,
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
    #pname_fid_map = {
        #12 : 'DIM1',
        #13 : 'DIM2',
        #14 : 'DIM3',
        #15 : 'DIM3',
    #}
    def update_by_pname_fid(self, pname_fid, value):
        if isinstance(pname_fid, string_types) and pname_fid.startswith('DIM'):
            num = int(pname_fid[3:])
            self.dim[num - 1] = value
        else:
            raise NotImplementedError('PBARL Type=%r name=%r has not been implemented' % (
                self.Type, self.pname_fid))

    def __init__(self, pid, mid, Type, dim, group='MSCBML0', nsm=0., comment=''):
        """
        Creates a PBARL card, which defines A, I1, I2, I12, and J using
        dimensions rather than explicit values.

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        Type : str
            type of the bar
            {ROD, TUBE, TUBE2, I, CHAN, T, BOX, BAR, CROSS, H, T1, I1,
            CHAN1, Z, CHAN2, T2, BOX1, HEXA, HAT, HAT1, DBOX}
        dim : List[float]
            dimensions for cross-section corresponding to Type;
            the length varies
        group : str; default='MSCBML0'
            this parameter can lead to a very broken deck with a very
            bad error message; don't touch it!
        nsm : float; default=0.
           non-structural mass
        comment : str; default=''
            a comment for the card

        The shear center and neutral axis do not coincide when:
           - Type = I and dim2 != dim3
           - Type = CHAN, CHAN1, CHAN2
           - Type = T
           - Type = T1, T2
           - Type = BOX1
           - Type = HAT, HAT1
           - Type = DBOX
        """
        LineProperty.__init__(self)
        if comment:
            self.comment = comment

        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        self.group = group
        #: Section Type (e.g. 'ROD', 'TUBE', 'I', 'H')
        self.Type = Type
        self.dim = dim
        #: non-structural mass
        self.nsm = nsm
        #ndim = self.valid_types[Type]
        #assert len(dim) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(dim))
        #self.validate()
        #area = self.Area()
        #assert area > 0, 'Type=%s dim=%s A=%s\n%s' % (self.Type, self.dim, area, str(self))
        self.mid_ref = None

    def validate(self):
        if self.Type not in self.valid_types:
            msg = ('Invalid PBARL Type, Type=%s '
                   'valid_types=%s' % (self.Type, self.valid_types.keys()))
            raise ValueError(msg)

        ndim = self.valid_types[self.Type]
        if not isinstance(self.dim, list):
            msg = 'PBARL pid=%s; dim must be a list; type=%r' % (self.pid, type(self.dim))
            raise TypeError(msg)
        if len(self.dim) != ndim:
            msg = 'dim=%s len(dim)=%s Type=%s len(dimType)=%s' % (
                self.dim, len(self.dim), self.Type,
                self.valid_types[self.Type])
            raise RuntimeError(msg)

        assert len(self.dim) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(self.dim))
        if not isinstance(self.group, string_types):
            raise TypeError('Invalid group; pid=%s group=%r' % (self.pid, self.group))
        #if self.group != 'MSCBML0':
            #msg = 'Invalid group; pid=%s group=%r expected=[MSCBML0]' % (self.pid, self.group)
            #raise ValueError(msg)

        assert None not in self.dim

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBARL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        group = string_or_blank(card, 3, 'group', 'MSCBML0')
        Type = string(card, 4, 'Type')

        ndim = cls.valid_types[Type]

        dim = []
        for i in range(ndim):
            dimi = double(card, 9 + i, 'ndim=%s; dim%i' % (ndim, i + 1))
            dim.append(dimi)

        #: dimension list
        assert len(dim) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(dim))
        #assert len(dims) == len(self.dim), 'PBARL ndim=%s len(dims)=%s' % (ndim, len(self.dim))

        nsm = double_or_blank(card, 9 + ndim, 'nsm', 0.0)
        return PBARL(pid, mid, Type, dim, group=group, nsm=nsm, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        pid = data[0]
        mid = data[1]
        group = data[2].strip()
        Type = data[3].strip()
        dim = list(data[4:-1])
        nsm = data[-1]
        return PBARL(pid, mid, Type, dim, group=group, nsm=nsm, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PBARL mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def uncross_reference(self):
        self.mid = self.Mid()
        self.mid_ref = None

    def _verify(self, xref):
        pid = self.pid
        mid = self.Mid()
        A = self.Area()
        try:
            J = self.J()
        except NotImplementedError:
            msg = "J is not implemented for pid.type=%s pid.Type=%s" % (self.type, self.Type)
            print(msg)
            J = 0.0
        nsm = self.Nsm()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(A, float), 'pid=%r' % A
        assert isinstance(J, float), 'cid=%r' % J
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        if xref:
            mpl = self.MassPerLength()
            assert isinstance(mpl, float), 'mass_per_length=%r' % mpl

    @property
    def Type(self):
        """gets Type"""
        return self.beam_type
    @Type.setter
    def Type(self, beam_type):
        """sets Type"""
        self.beam_type = beam_type

    def Area(self):
        """
        Gets the area :math:`A` of the CBAR.
        """
        return _bar_areaL('PBARL', self.beam_type, self.dim, self)

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
        #if self.beam_type in ['ROD']:
            #assert len(self.dim) == 1, 'dim=%r' % self.dim
            #r = self.dim[0]
            ##Ix = pi*r**4/4.
            ##J = pi*r**4/2.
            #(Ix, Iy, Ixy) = self.I1_I2_I12()

        #elif self.beam_type in ['BAR']:
            #assert len(self.dim) == 2, 'dim=%r' % self.dim
            #b, h = self.dim
            ##Ix = self.I11()
            ##Iy = self.I22()
            #(Ix, Iy, Ixy) = self.I1_I2_I12()
            ##print
            ##J = Ix + Iy
        #else:
            #msg = 'I11 for beam_type=%r dim=%r on PBARL is not supported' % (
                #self.beam_type, self.dim)
            #raise NotImplementedError(msg)
        #return Ix

    #def I12(self):
        #return self.I12()

    def _points(self, beam_type, dim):
        if beam_type == 'BAR':  # origin ar center
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
        elif beam_type == 'CROSS':
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
        elif beam_type in ['HEXA']:
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

        elif beam_type in ['I']:
            (d1, d2, d3) = dim
            raise NotImplementedError('PBARL beam_type=%r' % beam_type)

        elif beam_type in ['H']:
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

        elif beam_type in ['T2']:
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
            msg = '_points for beam_type=%r dim=%r on PBARL is not supported' % (self.beam_type, self.dim)
            raise NotImplementedError(msg)
        return array(points), Area

    def J(self):
        if self.beam_type in ['ROD']:
            r = self.dim[0]
            Ixx = pi * r**4 / 4.
            Iyy = Ixx
            unused_Ixy = 0.
        elif self.beam_type in ['TUBE']:
            rout, rin = self.dim
            #rin = rout - 2*t
            Ixx = pi * (rout**4 - rin**4) / 4.
            Iyy = Ixx
            unused_Ixy = 0.
        elif self.beam_type in ['TUBE2']:
            rout, t = self.dim
            rin = rout - 2*t
            Ixx = pi * (rout**4 - rin**4) / 4.
            Iyy = Ixx
            unused_Ixy = 0.
        elif self.beam_type in ['BOX']:
            (d1, d2, d3, d4) = self.dim
            hout = d2
            wout = d1
            hin = d2 - 2. * d3
            win = d1 - 2. * d4
            points, unused_Area = self._points('BAR', [hout, wout])
            yi = points[0, :-1]
            yip1 = points[0, 1:]
            xi = points[1, :-1]
            xip1 = points[1, 1:]

            #: .. seealso:: http://en.wikipedia.org/wiki/Area_moment_of_inertia
            ai = xi*yip1 - xip1*yi
            Ixx1 = 1/12 * sum((yi**2 + yi * yip1 + yip1**2)*ai)
            Iyy1 = 1/12 * sum((xi**2 + xi * xip1 + xip1**2)*ai)
            #Ixy1 = 1/24*sum((xi*yip1 + 2*xi*yi + 2*xip1*yip1 + xip1*yi)*ai)

            points, unused_Area = self._points('BAR', [hin, win])
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

        #elif self.beam_type in ['BAR']:
            #assert len(self.dim) == 2, 'dim=%r' % self.dim
            #b, h = self.dim
            #(Ix, Iy, Ixy) = self.I1_I2_I12()
            #J = Ix + Iy
        elif self.beam_type in ['BAR', 'CROSS', 'HEXA', 'T2', 'H']:
            points, unused_Area = self._points(self.beam_type, self.dim)
            yi = points[0, :-1]
            yip1 = points[0, 1:]

            xi = points[1, :-1]
            xip1 = points[1, 1:]

            #: .. seealso:: http://en.wikipedia.org/wiki/Area_moment_of_inertia
            ai = xi*yip1 - xip1*yi
            Ixx = 1/12 * sum((yi**2 + yi*yip1+yip1**2)*ai)
            Iyy = 1/12 * sum((xi**2 + xi*xip1+xip1**2)*ai)
            #Ixy = 1/24*sum((xi*yip1 + 2*xi*yi + 2*xip1*yip1 + xip1*yi)*ai)
        elif self.beam_type == 'I':
            # http://www.efunda.com/designstandards/beams/SquareIBeam.cfm
            # d - outside height
            # h - inside height
            # b - base
            # t - l thickness
            # s - web thickness
            #(b, d, t, s) = self.dim
            #h = d - 2 * s
            #cx = b / 2.
            #cy = d / 2.
            (d, b1, b2, t, s1, s2) = self.dim
            if b1 != b2:
                msg = 'J for beam_type=%r dim=%r on PBARL b1 != b2 is not supported' % (
                    self.beam_type, self.dim)
                raise NotImplementedError(msg)
            if s1 != s2:
                msg = 'J for beam_type=%r dim=%r on PBARL s1 != s2 is not supported' % (
                    self.beam_type, self.dim)
                raise NotImplementedError(msg)
            h = d - b1 - b2
            s = s1
            b = b1
            Ixx = (b*d**3-h**3*(b-t)) / 12.
            Iyy = (2.*s*b**3 + h*t**3) / 12.
        #elif self.beam_type == 'T': # test
            # http://www.amesweb.info/SectionalPropertiesTabs/SectionalPropertiesTbeam.aspx
            # http://www.amesweb.info/SectionalPropertiesTabs/SectionalPropertiesTbeam.aspx
            # d - outside height
            # h - inside height
            # b - base
            # t - l thickness
            # s - web thickness
            #(b, d, t, s) = self.dim
            #h = d - 2 * s
            #(b, d, s, t) = self.dim
            #if b1 != b2:
                #msg = 'J for beam_type=%r dim=%r on PBARL b1 != b2 is not supported' % (
                    #self.beam_type, self.dim)
                #raise NotImplementedError(msg)
            #if s1 != s2:
                #msg = 'J for beam_type=%r dim=%r on PBARL s1 != s2 is not supported' % (
                    #self.beam_type, self.dim)
                #raise NotImplementedError(msg)
            #h = d - b1 - b2
            #s = s1
            #b = b1

            # http://www.engineersedge.com/material_science/moment-inertia-gyration-6.htm
            #y = d**2*t+s**2*(b-t)/(2*(b*s+h*t))
            #Ixx = (t*y**3 + b*(d-y)**3 - (b-t)*(d-y-s)**3)/3.
            #Iyy = t**3*(h-s)/12. + b**3*s/12.
            #A = b*s + h*t

        elif self.beam_type == 'C':
            # http://www.efunda.com/math/areas/squarechannel.cfm
            # d - outside height
            # h - inside height
            # b - base
            # t - l thickness
            # s - web thickness
            (b, d, t, s) = self.dim
            h = d - 2 * s
            #cx = (2.*b**2*s + h*t**2)/(2*b*d - 2*h*(b-t))
            #cy = d / 2.
            Ixx = (b * d**3 - h **3 * (b-t)) / 12.
            #Iyx = (2.*s*b**3 + h*t**3)/3 - A*cx**2
        else:
            msg = 'J for beam_type=%r dim=%r on PBARL is not supported' % (self.beam_type, self.dim)
            raise NotImplementedError(msg)

        #: .. seealso:: http://en.wikipedia.org/wiki/Perpendicular_axis_theorem
        J = Ixx + Iyy
        return J

    def I22(self):
        return self.I2()

    def raw_fields(self):
        list_fields = ['PBARL', self.pid, self.Mid(), self.group, self.beam_type,
                       None, None, None, None] + self.dim + [self.nsm]
        return list_fields

    def repr_fields(self):
        group = set_blank_if_default(self.group, 'MSCBML0')
        ndim = self.valid_types[self.beam_type]
        assert len(self.dim) == ndim, 'PBARL ndim=%s len(dims)=%s' % (ndim, len(self.dim))
        list_fields = ['PBARL', self.pid, self.Mid(), group, self.beam_type, None,
                       None, None, None] + self.dim + [self.nsm]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBRSECT(LineProperty):
    """
    not done
    """
    type = 'PBRSECT'

    def __init__(self, pid, mid, form, options, comment=''):
        LineProperty.__init__(self)
        if comment:
            self.comment = comment

        #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        self.form = form

        self.nsm = 0.
        self.t = None
        self.outp = None
        self.brp1 = None
        for key, value in options.items():
            key = key.upper()
            #if key == 'NSM':
                #self.nsm = float(value)
            #elif key == 'OUTP':
                #self.outp = int(value)
            #elif key == 'BRP(1)':
                #self.brp1 = int(value)
            #elif key == 'T':
                #self.t = float(value)
            #else:
            raise NotImplementedError('PBRSECT.pid=%s key=%r value=%r' % (pid, key, value))

        self._validate_input()
        self.mid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBRSECT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : List[str]
            this card is special and is not a ``BDFCard`` like other cards
        comment : str; default=''
            a comment for the card
        """
        line0 = card[0]
        if '\t' in line0:
            line0 = line0.expandtabs()

        bdf_card = BDFCard(to_fields([line0], 'PBMSECT'))
        unused_line0_eq = line0[16:]
        lines_joined = ''.join(card[1:]).replace(' ', '')

        if lines_joined:
            fields = lines_joined.split(',')
            slines = [field.split('=') for field in fields]
            #C:\MSC.Software\MSC.Nastran\msc20051\nast\tpl\zbr3.dat
            #slines = [
                #[u'OUTP', u'201'],
                #[u'T', u'1.0'],
                #[u'BRP', u'202'],
                #[u'T(11)', u'[1.2'],
                #[u'PT', u'(202'], [u'224)]'],
                #[u'T(12)', u'[1.2'],
                #[u'PT', u'(224'],
                #[u'205)]'],
            #]
            try:
                options = {key : value for (key, value) in slines}
            except:
                print('PBRSECT slines=%s' % slines)
                raise
        else:
            options = {}

        pid = integer(bdf_card, 1, 'pid')
        mid = integer(bdf_card, 2, 'mid')
        form = string_or_blank(bdf_card, 3, 'form')
        assert form in ['GS', 'OP', 'CP'], 'pid=%s form=%r' % (pid, form)

        return PBRSECT(pid, mid, form, options, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        #pid = data[0]
        #mid = data[1]
        #group = data[2].strip()
        #beam_type = data[3].strip()
        #dim = list(data[4:-1])
        #nsm = data[-1]
        #print("group = %r" % self.group)
        #print("beam_type  = %r" % self.beam_type)
        #print("dim = ",self.dim)
        #print(str(self))
        #print("*PBARL = ",data)
        raise NotImplementedError('not finished...')
        #return PBRSECT(pid, mid, group, beam_type, dim, nsm, comment=comment)

    def _validate_input(self):
        pass

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PBRSECT mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def uncross_reference(self):
        self.mid = self.Mid()
        self.mid_ref = None

    def _verify(self, xref):
        pid = self.pid
        mid = self.Mid()
        #A = self.Area()
        #J = self.J()
        #nsm = self.Nsm()
        #mpl = self.MassPerLength()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid
        #assert isinstance(A, float), 'pid=%r' % A
        #assert isinstance(J, float), 'cid=%r' % J
        #assert isinstance(nsm, float), 'nsm=%r' % nsm
        #assert isinstance(mpl, float), 'mass_per_length=%r' % mpl

    def Area(self):
        """
        Gets the area :math:`A` of the CBAR.
        """
        return 0.
        #raise NotImplementedError('Area is not implemented for PBRSECT')

    def Nsm(self):
        """
        Gets the non-structural mass :math:`nsm` of the CBAR.
        """
        return 0.
        #raise NotImplementedError('Nsm is not implemented for PBRSECT')

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
        raise NotImplementedError('I11 is not implemented for PBRSECT')

    #def I12(self):
        #return self.I12()

    def J(self):
        raise NotImplementedError('J is not implemented for PBRSECT')

    def I22(self):
        raise NotImplementedError('I22 is not implemented for PBRSECT')

    def raw_fields(self):
        """not done..."""
        list_fields = ['PBRSECT', self.pid, self.Mid(), self.form]
                       #None, None, None, None] + self.dim + [self.nsm]
        return list_fields

    def repr_fields(self):
        """not done..."""
        list_fields = ['PBRSECT', self.pid, self.Mid(), self.form]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBEAM3(LineProperty):  # not done, cleanup
    type = 'PBEAM3'

    def __init__(self, pid, mid, A, iz, iy, iyz, j, nsm=0.,
                 cy=0., cz=0., dy=0., dz=0., ey=0., ez=0., fy=0., fz=0., comment=''):
        LineProperty.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.mid = mid

        self.A = A
        self.iz = iz
        self.iy = iy
        self.iyz = iyz
        self.j = j
        self.nsm = nsm

        self.cy = cy
        self.cz = cz
        self.dy = dy
        self.dz = dz
        self.ey = ey
        self.ez = ez
        self.fy = fy
        self.fz = fz
        self.mid_ref = None

    def add_card(self, card, comment=''):
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')

        A = double(card, 3, 'A')
        iz = double(card, 4, 'Iz')
        iy = double(card, 5, 'Iy')
        iyz = double_or_blank(card, 6, 'Iyz', 0.0)
        j = double_or_blank(card, 7, 'J', self.iy + self.iz)
        nsm = double_or_blank(card, 8, 'nsm', 0.0)

        cy = double(card, 9, 'cy')
        cz = double(card, 10, 'cz')

        dy = double(card, 11, 'dy')
        dz = double(card, 12, 'dz')

        ey = double(card, 13, 'ey')
        ez = double(card, 14, 'ez')

        fy = double(card, 15, 'fy')
        fz = double(card, 16, 'fz')
        # more...
        return PBEAM3(pid, mid, A, iz, iy, iyz, j, nsm=nsm,
                      cy=cy, cz=cz, dy=dy, dz=dz, ey=ey, ez=ez, fy=fy, fz=fz, comment=comment)

    def add_op2_data(self, data, comment=''):
        if comment:
            self.comment = comment
        raise NotImplementedError(data)

    def Nsm(self):
        """
        Gets the non-structural mass :math:`nsm`.
        .. warning:: nsm field not supported fully on PBEAM3 card
        """
        return self.nsm

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PBEAM3 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def uncross_reference(self):
        self.mid = self.Mid()
        self.mid_ref = None

    def repr_fields(self):
        """.. todo:: not done"""
        list_fields = [
            'PBEAM3', self.pid, self.Mid(), self.A, self.iz, self.iy, self.iyz, self.j, self.nsm,
            self.cy, self.cz, self.dy, self.dz, self.ey, self.ez, self.fy, self.fz]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBEND(LineProperty):
    """
    MSC/NX Option A

    +-------+------+-------+-----+----+----+--------+----+--------+
    |   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
    +=======+======+=======+=====+====+====+========+====+========+
    | PBEND | PID  |  MID  | A   | I1 | I2 |   J    | RB | THETAB |
    +-------+------+-------+-----+----+----+--------+----+--------+
    |       |  C1  |  C2   | D1  | D2 | E1 |   E2   | F1 |   F2   |
    +-------+------+-------+-----+----+----+--------+----+--------+
    |       |  K1  |  K2   | NSM | RC | ZC | DELTAN |    |        |
    +-------+------+-------+-----+----+----+--------+----+--------+

    MSC Option B

    +-------+------+-------+-----+----+----+--------+----+--------+
    |   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
    +=======+======+=======+=====+====+====+========+====+========+
    | PBEND | PID  |  MID  | FSI | RM | T  |   P    | RB | THETAB |
    +-------+------+-------+-----+----+----+--------+----+--------+
    |       |      |       | NSM | RC | ZC |        |    |        |
    +-------+------+-------+-----+----+----+--------+----+--------+

    NX Option B

    +-------+------+-------+-----+----+----+--------+----+--------+
    |   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
    +=======+======+=======+=====+====+====+========+====+========+
    | PBEND | PID  |  MID  | FSI | RM | T  |   P    | RB | THETAB |
    +-------+------+-------+-----+----+----+--------+----+--------+
    |       | SACL | ALPHA | NSM | RC | ZC | FLANGE |    |        |
    +-------+------+-------+-----+----+----+--------+----+--------+
    |       |  KX  |  KY   | KZ  |    | SY |   SZ   |    |        |
    +-------+------+-------+-----+----+----+--------+----+--------+
    """
    type = 'PBEND'

    def __init__(self, pid, mid, beam_type, A, i1, i2, j,
                 c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                 nsm, rc, zc, delta_n, fsi, rm, t, p,
                 rb, theta_b, comment=''):
        LineProperty.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.mid = mid
        self.beam_type = beam_type
        self.A = A
        self.i1 = i1
        self.i2 = i2
        self.j = j
        self.c1 = c1
        self.c2 = c2
        self.d1 = d1
        self.d2 = d2
        self.e1 = e1
        self.e2 = e2
        self.f1 = f1
        self.f2 = f2
        self.k1 = k1
        self.k2 = k2
        self.nsm = nsm
        self.rc = rc
        self.zc = zc
        self.delta_n = delta_n
        self.fsi = fsi
        self.rm = rm
        self.t = t
        self.p = p
        self.rb = rb
        self.theta_b = theta_b
        self.mid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PBEND card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')

        value3 = integer_or_double(card, 3, 'Area/FSI')
        #print("PBEND: area/fsi=%s" % value3)

        # MSC/NX option A
        A = None
        i1 = None
        i2 = None
        j = None
        c1 = None
        c2 = None
        d1 = None
        d2 = None
        e1 = None
        e2 = None
        f1 = None
        f2 = None
        k1 = None
        k2 = None
        delta_n = None

        # MSC option B
        rm = None
        t = None
        p = None

        # NX option B
        #sacl = None
        #alpha = None
        #flange = None
        #kx = None
        #ky = None
        #kz = None
        #sy = None
        #sz = None
        if isinstance(value3, float):
            fsi = 0
            beam_type = 1
            #: Area of the beam cross section
            A = double(card, 3, 'A')

            #: Area moments of inertia in planes 1 and 2.
            i1 = double(card, 4, 'I1')
            i2 = double(card, 5, 'I2')

            #: Torsional stiffness :math:`J`
            j = double(card, 6, 'J')

            # line2
            #: The r,z locations from the geometric centroid for stress
            #: data recovery.
            c1 = double_or_blank(card, 9, 'c1', 0.)
            c2 = double_or_blank(card, 10, 'c2', 0.)
            d1 = double_or_blank(card, 11, 'd1', 0.)
            d2 = double_or_blank(card, 12, 'd2', 0.)
            e1 = double_or_blank(card, 13, 'e1', 0.)
            e2 = double_or_blank(card, 14, 'e2', 0.)
            f1 = double_or_blank(card, 15, 'f1', 0.)
            f2 = double_or_blank(card, 16, 'f2', 0.)

            # line 3
            #: Shear stiffness factor K in K*A*G for plane 1.
            k1 = double_or_blank(card, 17, 'k1')
            #: Shear stiffness factor K in K*A*G for plane 2.
            k2 = double_or_blank(card, 18, 'k2')

            #: Nonstructural mass per unit length.
            nsm = double_or_blank(card, 19, 'nsm', 0.)

            #: Radial offset of the geometric centroid from points GA and GB.
            rc = double_or_blank(card, 20, 'rc', 0.)

            #: Offset of the geometric centroid in a direction perpendicular
            #: to the plane of points GA and GB and vector v
            zc = double_or_blank(card, 21, 'zc', 0.)

            #: Radial offset of the neutral axis from the geometric
            #: centroid, positive is toward the center of curvature
            delta_n = double_or_blank(card, 22, 'delta_n', 0.)

        elif isinstance(value3, integer_types):  # alternate form
            beam_type = 2
            #: Flag selecting the flexibility and stress intensification
            #: factors. See Remark 3. (Integer = 1, 2, or 3)
            fsi = integer(card, 3, 'fsi')
            if fsi in [1, 2, 3]:
                # assuming MSC
                #: Mean cross-sectional radius of the curved pipe
                rm = double(card, 4, 'rm')

                #: Wall thickness of the curved pipe
                t = double(card, 5, 't')

                #: Internal pressure
                p = double_or_blank(card, 6, 'p')

                # line3
                # Non-structural mass :math:`nsm`
                nsm = double_or_blank(card, 11, 'nsm', 0.)
                rc = double_or_blank(card, 12, 'rc', 0.)
                zc = double_or_blank(card, 13, 'zc', 0.)
            elif fsi in [4, 5, 6]:
                # Non-structural mass :math:`nsm`
                nsm = double_or_blank(card, 11, 'nsm', 0.)
                rc = double_or_blank(card, 12, 'rc', 0.)
                zc = double_or_blank(card, 13, 'zc', 0.)

                #sacl = double_or_blank(card, 9, 'sacl')
                #alpha = double_or_blank(card, 10, 'alpha', 0.)
                #flange = integer_or_blank(card, 15, 'flange', 0)
                #kx = double_or_blank(card, 18, 'kx', 1.0)
                #ky = double_or_blank(card, 19, 'ky', 1.0)
                #kz = double_or_blank(card, 20, 'kz', 1.0)
                #sy = double_or_blank(card, 22, 'sy', 1.0)
                #sz = double_or_blank(card, 23, 'sz', 1.0)
            else:
                assert fsi in [1, 2, 3, 4, 5, 6], 'pid=%s fsi=%s\ncard:%s' % (pid, fsi, card)
        else:
            raise RuntimeError('Area/FSI on CBEND must be defined...')
        assert fsi in [0, 1, 2, 3, 4, 5, 6], 'pid=%s fsi=%s\ncard:%s' % (pid, fsi, card)

        #: Bend radius of the line of centroids
        rb = double_or_blank(card, 7, 'rb')

        #: Arc angle :math:`\theta_B` of element  (optional)
        theta_b = double_or_blank(card, 8, 'thetab')
        assert len(card) <= 23, 'len(PBEND card) = %i\ncard=%s' % (len(card), card)
        return PBEND(pid, mid, beam_type, A, i1, i2, j, c1, c2, d1, d2,
                     e1, e2, f1, f2, k1, k2, nsm,
                     rc, zc, delta_n, fsi, rm, t,
                     p, rb, theta_b, comment=comment)

    def validate(self):
        """card checking method"""
        if self.delta_n is not None and not isinstance(self.delta_n, float_types):
            raise RuntimeError('delta_n=%r must be None or a float; type=%s; fsi=%s\n%s' % (
                self.delta_n, type(self.delta_n), self.fsi, str(self)))

    #def Nsm(self):
        #""".. warning:: nsm field not supported fully on PBEND card"""
        #raise RuntimeError(self.nsm[0])
        #return self.nsm

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PBEND mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def uncross_reference(self):
        self.mid = self.Mid()
        self.mid_ref = None

    def MassPerLength(self):
        """m/L = rho*A + nsm"""
        rho = self.mid_ref.Rho()
        return self.A * rho + self.nsm

    def raw_fields(self):
        return self.repr_fields()

    def repr_fields(self):
        list_fields = ['PBEND', self.pid, self.Mid(), ]  # other
        if self.beam_type == 1:
            list_fields += [self.A, self.i1, self.i2, self.j, self.rb,
                            self.theta_b, self.c1, self.c2, self.d1, self.d2,
                            self.e1, self.e2, self.f1, self.f2, self.k1, self.k2,
                            self.nsm, self.rc, self.zc, self.delta_n]
            #print("beam_type=0 I1=%s I2=%s; J=%s RM=%s T=%s P=%s" % (
                #self.i1, self.i2, self.j, self.rm, self.t, self.p), list_fields)
        elif self.beam_type == 2:
            list_fields += [self.fsi, self.rm, self.t, self.p, self.rb,
                            self.theta_b, None, None, self.nsm, self.rc, self.zc]
        elif self.beam_type == 0:
            # dunno
            list_fields += [self.A, self.i1, self.i2, self.j, self.rb,
                            self.theta_b, self.c1, self.c2, self.d1, self.d2,
                            self.e1, self.e2, self.f1, self.f2, self.k1, self.k2,
                            self.nsm, self.rc, self.zc, self.delta_n]
            #print("beam_type=0 I1=%s I2=%s; J=%s RM=%s T=%s P=%s" % (
                #self.i1, self.i2, self.j, self.rm, self.t, self.p), list_fields)
        else:
            raise ValueError('only beam_type=1 and 2 supported; beam_type/fsi=%s' % self.beam_type)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
