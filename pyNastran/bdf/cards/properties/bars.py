# pylint: disable=C0103,R0914,R0902,R0913
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
from __future__ import annotations
from itertools import count
from typing import List, Tuple, Union, Any, TYPE_CHECKING
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
from pyNastran.bdf.bdf_interface.utils import to_fields
from pyNastran.utils.numpy_utils import float_types
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


def Iyy_beam(b, h):
    """gets the Iyy for a solid beam"""
    return 1 / 12. * b * h ** 3


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


def _IAreaL(prop, dim):
    # type: (Any, List[float]) -> Tuple[float, float, float, float]
    beam_type = prop.beam_type
    if beam_type == 'ROD':
        A, I1, I2, I12 = rod_section(prop.type, prop.beam_type, dim, prop)
    elif beam_type == 'TUBE':
        A, I1, I2, I12 = tube_section(prop.type, prop.beam_type, dim, prop)
    elif beam_type == 'TUBE2':
        A, I1, I2, I12 = tube2_section(prop.type, prop.beam_type, dim, prop)

    elif beam_type == 'I':
        # |  ------------
        # |  |    A     | d6
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

    elif beam_type == 'BAR':
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
        msg = 'beam_type=%r is not supported for %s class...' % (
            beam_type, prop.type)
        raise NotImplementedError(msg)
    return (A, Iyy, Izz, Iyz)

class LineProperty(Property):
    def __init__(self):
        self.beam_type = None
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

    def Area(self) -> float:
        """gets area"""
        return self.A

    def Nsm(self) -> float:
        """gets nonstructural mass per unit length"""
        return self.nsm

    def J(self) -> float:
        """gets J"""
        return self.j

    def I11(self) -> float:
        """gets I11"""
        return self.i1

    def I22(self) -> float:
        """gets I22"""
        return self.i2

    def Rho(self) -> float:
        """gets the material density"""
        return self.mid_ref.rho

    def E(self) -> float:
        """gets the material Young's ratio"""
        return self.mid_ref.E

    def G(self) -> float:
        """gets the material Shear ratio"""
        return self.mid_ref.G

    def Nu(self) -> float:
        """gets the material Poisson's ratio"""
        return self.mid_ref.nu


def A_I1_I2_I12(prop, beam_type: str,
                dim: List[float]) -> Tuple[float, float, float, float]:
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
    if beam_type == 'ROD':
        A, I1, I2, I12 = rod_section(prop.type, beam_type, dim, prop)
    elif beam_type == 'TUBE':
        A, I1, I2, I12 = tube_section(prop.type, beam_type, dim, prop)

    elif beam_type == 'BAR':
        A, I1, I2, I12 = bar_section(prop.type, beam_type, dim, prop)
    elif beam_type == 'BOX':
        A, I1, I2, I12 = box_section(prop.type, beam_type, dim, prop)
    elif beam_type == 'L':
        A, I1, I2, I12 = l_section(prop.type, beam_type, dim, prop)
    else:  # pragma: no cover
        msg = 'A_I1_I2_I12; beam_type=%s is not supported for %s class...' % (
            beam_type, prop.type)
        raise NotImplementedError(msg)
    assert A == prop.Area(), prop
    return A, I1, I2, I12

def _bar_areaL(class_name: str, beam_type: str, dim: List[float],
               prop: Any) -> float:
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
        A, I1, I2, I12 = rod_section(class_name, beam_type, dim, prop)
    elif beam_type == 'TUBE':
        A, I1, I2, I12 = tube_section(class_name, beam_type, dim, prop)
    elif beam_type == 'TUBE2':
        A, I1, I2, I12 = tube2_section(class_name, beam_type, dim, prop)

    elif beam_type == 'I':
        A, I1, I2, I12 = i_section(class_name, beam_type, dim, prop)
    elif beam_type == 'I1':
        A, I1, I2, I12 = i1_section(class_name, beam_type, dim, prop)

    elif beam_type == 'L':
        A, I1, I2, I12 = l_section(class_name, beam_type, dim, prop)

    elif beam_type == 'CHAN':
        A, I1, I2, I12 = chan_section(class_name, beam_type, dim, prop)
    elif beam_type == 'CHAN1':
        A, I1, I2, I12 = chan1_section(class_name, beam_type, dim, prop)
    elif beam_type == 'CHAN2':
        A, I1, I2, I12 = chan2_section(class_name, beam_type, dim, prop)

    elif beam_type == 'T':
        A, I1, I2, I12 = t_section(class_name, beam_type, dim, prop)
    elif beam_type == 'T1':
        A, I1, I2, I12 = t1_section(class_name, beam_type, dim, prop)
    elif beam_type == 'T2':
        A, I1, I2, I12 = t2_section(class_name, beam_type, dim, prop)

    elif beam_type == 'BOX':
        A, I1, I2, I12 = box_section(class_name, beam_type, dim, prop)
    elif beam_type == 'BOX1':
        A, I1, I2, I12 = box1_section(class_name, beam_type, dim, prop)

    elif beam_type == 'BAR':
        A, I1, I2, I12 = bar_section(class_name, beam_type, dim, prop)

    elif beam_type == 'CROSS':
        A, I1, I2, I12 = cross_section(class_name, beam_type, dim, prop)

    elif beam_type == 'H':
        A, I1, I2, I12 = h_section(class_name, beam_type, dim, prop)
    elif beam_type == 'Z':
        A, I1, I2, I12 = zee_section(class_name, beam_type, dim, prop)

    elif beam_type == 'HEXA':
        A, I1, I2, I12 = hexa_section(class_name, beam_type, dim, prop)
    elif beam_type == 'HAT':
        A, I1, I2, I12 = hat_section(class_name, beam_type, dim, prop)
    elif beam_type == 'HAT1':
        A, I1, I2, I12 = hat1_section(class_name, beam_type, dim, prop)
    elif beam_type == 'DBOX':
        A, I1, I2, I12 = dbox_section(class_name, beam_type, dim, prop)
    else:  # pragma: no cover
        msg = 'areaL; beam_type=%s is not supported for %s class...' % (
            beam_type, class_name)
        raise NotImplementedError(msg)
    assert A > 0, 'beam_type=%r dim=%r A=%s\n%s' % (beam_type, dim, A, prop)
    return A

def rod_section(class_name, beam_type, dim, prop):
    # This is a circle if you couldn't tell...
    #   __C__
    #  /     \
    # |       |
    # F       D
    # |       |
    #  \     /
    #   \_E_/
    #
    R = dim[0]
    A = pi * R ** 2
    I1 = A * R ** 2 / 4.
    I2 = I1
    I12 = 0.
    return A, I1, I2, I12

def tube_section(class_type, beam_type, dim, prop):
    r_outer, r_inner = dim
    assert r_outer > r_inner, 'TUBE; r_outer=%s r_inner=%s' % (r_outer, r_inner)
    A = pi * (r_outer ** 2 - r_inner ** 2)
    I1 = pi * (r_outer ** 4 - r_inner ** 4) / 4.
    I2 = I1
    I12 = 0.
    return A, I1, I2, I12

def tube2_section(class_type, beam_type, dim, prop):
    r_outer, thick = dim
    r_inner = r_outer - thick
    assert r_outer > r_inner, 'TUBE2; r_outer=%s r_inner=%s' % (r_outer, r_inner)

    # same as tube
    A = pi * (r_outer ** 2 - r_inner ** 2)
    I1 = pi * (r_outer ** 4 - r_inner ** 4) / 4.
    I2 = I1
    I12 = 0.
    return A, I1, I2, I12

def bar_section(class_type, beam_type, dim, prop):
    #per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
    # <------> D1
    #
    # F------C  ^
    # |      |  |
    # |      |  | D2
    # |      |  |
    # E------D  v
    b, h = dim
    A = b * h
    I1 = b * h ** 3 / 12.
    I2 = h * b ** 3 / 12.
    I12 = 0.
    #J = b*h**3*(1/3. - 0.21*h/b*(1-h**4/(12*b**4)))
    return A, I1, I2, I12

def box_section(class_name, beam_type, dim, prop):
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
    assert dim1 > 2.*dim4, 'BOX; required: dim1 > 2*dim4; dim1=%s dim4=%s\n%s' % (dim1, dim4, prop)
    assert dim2 > 2.*dim3, 'BOX; required: dim2 > 2*dim3; dim2=%s dim3=%s\n%s' % (dim2, dim3, prop)

    # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
    b = dim1
    h = dim2
    t1 = dim3
    t2 = dim4
    bi = b - 2 * t2
    hi = h - 2 * t1
    A = b * h - bi * hi

    I1 = (b * h ** 3 - bi * hi ** 3) / 12
    I2 = (h * b ** 3 - bi * hi ** 3) / 12
    I12 = 0.
    #j = (2*t2*t1*(b-t2)**2*(h-t1)**2)/(b*t2+h*t1-t2**2-t1**2)
    return A, I1, I2, I12

def box1_section(class_name, beam_type, dim, prop):
    dim1, dim2, dim3, dim4, dim5, dim6 = dim
    h1 = dim3  # top
    w1 = dim1

    h2 = dim4  # btm
    A1 = (h1 + h2) * w1

    h3 = dim2 - h1 - h2  # left
    w3 = dim6

    w4 = dim5  # right
    A2 = h3 * (w3 + w4)
    A = A1 + A2
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def cross_section(class_type, beam_type, dim, prop):
    dim1, dim2, dim3, dim4 = dim
    h1 = dim3
    w1 = dim2

    h2 = dim4
    w2 = dim1
    A = h1 * w1 + h2 * w2
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def h_section(class_name, beam_type, dim, prop):
    dim1, dim2, dim3, dim4 = dim
    h1 = dim3
    w1 = dim2

    h2 = dim4
    w2 = dim1
    A = h1 * w1 + h2 * w2
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def t_section(class_name, beam_type, dim, prop):
    # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
    #       ^ y
    #       |
    #       |
    #  <----d1=d--->
    #  +-----------+ ^     ^
    #  | 1         | |     | d3=tf ----> z
    #  +---+  +----+ |     v
    #      |  |      |     ^
    #      |  |      | d2  |
    #      |  |      |     | hw
    #      |  |      |     |
    #      +--+      v     v
    #
    #      <--> d4=tw
    #
    dim1, dim2, dim3, dim4 = dim
    d = dim1 # bf
    tf = dim3
    tw = dim4
    h = dim2 - 0.5 * tf
    hw = dim2 - tf
    #print('tf=%s tw=%s bf=%s' % (tf, tw, bf))
    #print('h=%s hw=%s' % (h, hw))

    A = d*tf + hw*tw
    yna = hw*tw*(hw+tf)/(2*A)
    i1 = (
       (d*tf**3 + tw*hw**3) / 12. +
       hw*tw*(yna + 0.5 * (hw +tf))**2 + d*tf*yna**2
    )
    i2 = (tf*d**3 + hw*tw**3)/12.
    i12 = 0.
    #j = 1/3 * (tf**3*d + tw**3 * h)
    return A, i1, i2, i12

def t1_section(class_name, beam_type, dim, prop):
    dim1, dim2, dim3, dim4 = dim
    h1 = dim1
    w1 = dim3

    h2 = dim4
    w2 = dim2
    A = h1 * w1 + h2 * w2
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def t2_section(class_name, beam_type, dim, prop):
    #       ^ y
    #       |
    #       |
    #      +--+      ^     v
    #      |  |      |     ^
    #      |1 |      | d2  |
    #      |  |      |     | hw=bweb
    #      |  |      |     |
    #  +---+  +----+ |     v
    #  | 2 |  |    | |     | d3=tf ----> z
    #  +---+--+----+ v     v
    #  <----d1=d--->
    #
    #      <--> d4=tweb
    #
    dim1, dim2, dim3, dim4 = dim
    ball = dim1
    hall = dim2
    tflange = dim3
    tweb = dim4
    # bweb = hall - tflange
    bflange = ball - tweb
    #print(f'hall={hall} ball={ball} tweb={tweb} bflange={bflange}')
    assert hall - tflange > 0, f'bweb={hall-tflange} hall(dim2)={hall} tflange(dim3)={tflange}'
    assert ball - tweb > 0, f'bflange={bflange} ball(dim1)={ball} tweb(dim4)={tweb}'
    A = tweb * hall + tflange * bflange
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def zee_section(class_name, beam_type, dim, prop):
    #  bflange bweb
    # <--d1-><--d2-->
    # +-------------+            ^
    # |             |            |
    # +------+      |  ^         |
    #        |      |  |         |
    #        |      |  | d3=hweb |  d4 = hall
    #        |      |  v         |
    #        |      +---------+  |
    #        |                |  |
    #        +----------------+  v
    #
    dim1, dim2, dim3, dim4 = dim
    bflange = dim1
    bweb = dim2
    hweb = dim3
    hall = dim4
    tflange = hall - hweb # actually 2*tflange
    A = bflange * tflange + hall * bweb
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def hexa_section(class_name, beam_type, dim, prop):
    #     _______
    #   /        \     ^
    #  /          \    |
    # *            *   |  d3 = h
    #  \          /    |
    #   \        /     |
    #    \______/      v
    #           |d1| = wtri
    # <-----d2----> = wall
    #
    dim1, dim2, dim3 = dim
    hbox = dim3
    wall = dim2
    wtri = dim1
    wbox = wall - 2 * wtri
    assert wbox > 0, 'wbox=%s' % (wbox)
    A = hbox * wbox + hbox * wtri
    #print('hbox=%s wbox=%s hbox*wbox=%s 2*wtri*hbox=%s A=%s' % (
        #hbox, wbox, hbox*wbox, 2*wtri*hbox, A))
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def hat_section(class_name, beam_type, dim, prop):
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
    t = dim2
    wa = dim3
    assert wa > 2.*t, f'HAT; required: dim3 > 2*dim2; dim2={t} dim3={wa}; delta={wa-2*t}\n{prop}'
    #DIM3, CAN NOT BE LESS THAN THE SUM OF FLANGE
        #THICKNESSES, 2*DIM2
    hb = dim1 - 2. * t
    wc = dim4 + t
    assert hb > 0., f'HAT; required hb=dim1-2*dim2={hb} dim1={dim1} dim2={dim2}'
    A = wa * t + (2. * wc * t) + (2. * hb * t)
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def hat1_section(class_name, beam_type, dim, prop):
    # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
    #
    #        <-----d3=wmid----->
    #
    #        +-----------------+              ^
    #   d4   |                 |              |
    #        +----+-------+-d4-+  d4=t        |
    #        |    |       |    |              | d2=hall
    # +------+----+   x   +----+------+       |
    # |           |       |           | t=d2  |
    # +-----------+-------+-----------+       |  ^
    # |                               |       |  |
    # |                               |       |  | d5=hbox
    # +-------------------------------+       v  v
    #
    # <------------d1=wbox------------>
    dim1, dim2, dim3, dim4, dim5 = dim
    wbox = dim1
    hall = dim2
    unused_wmid = dim3
    t = dim4
    hbox = dim5
    hweb = hall - hbox - t
    assert hweb > 0, hweb

    Abox = hbox * wbox
    # we do a first pass assuming the hat is flat and then
    # add the leftover web sections
    Ahat = wbox * t + 2 * hweb * t
    A = Abox + Ahat

    b1 = dim1
    h = dim2
    b2 = dim3
    t2 = dim4
    t1 = dim5
    A0 = b1 * t1 + t2 * (b1 - b2) + 2 * t2 * (h - t1) + t2 * (b2 - 2 * t2)
    assert np.allclose(A, A0), 'A=%s A0=%s' % (A, A0)
    #assert dim2 > dim4+dim5, 'HAT1; required: dim2 > dim4+dim5; dim2=%s dim4=%s; dim5=%s\n%s' % (dim2, dim4, dim5, prop)
    #*DIM4+DIM5, CAN NOT BE LARGER THAN THE HEIGHT OF
        #THE HAT, DIM2.


    #h1 = w  # upper, horizontal lower bar (see HAT)
    #w1 = w0 - w3
    #A = 2. * (h0 * w0 + h1 * w1 + h2 * w2 + h3 * w3)
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def i_section(class_name, beam_type, dim, prop):
    # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf

    #  <----d3=b---->
    #
    #  1------------2         ^
    #  |      B     |  d6=tb  |
    # 12---11---4---3         |
    #       |   |             |
    #       |   |             |  d1=h
    #       | C | <--- d4=tw  |
    #       |   |             |
    #  9---10---5---6         |
    #  |      A     |  d5=ta  |
    #  8------------7         v
    #
    #  <----d2=a---->
    #
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
    yc = (0.5 * hw * (hw + ta)*tw + hf*tb*b) / A
    #ys = tb * hf * b**3/(tb * b**3 + ta * a**a)
    #yna = yc - ys
    i1 = (
       1 / 12 * (h*tb**3 + a*ta**3 + tw*hw**3)
       + (hf-yc)**2*b*tb+yc**2*a*ta
       + (yc - 0.5*(hw+ta))**2*hw*tw
    )
    i2 = (b**3 * tb + ta * a**3 + hw*tw**3) / 12.
    i12 = 0.
    #j = 1/3 * (tb**3*b + ta**3*a + tw**3*hf)
    #assert np.allclose(i1+i2, j), 'I1+I2=%s+%s=%s J=%s' % (i1, i2, i1 + i2, j)

    assert dim1 > dim5+dim6, 'I; required: dim1 > dim5+dim6; dim1=%s dim5=%s dim6=%s\n%s' % (dim1, dim5, dim6, prop)
    #THE SUM OF THE FLANGE THICKNESSES,DIM5 + DIM6,
    #CAN NOT BE GREATER THAT THE WEB HEIGHT, DIM1.
    return A, i1, i2, i12

def i1_section(class_name, beam_type, dim, prop):
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
    dim1, dim2, dim3, dim4 = dim
    w1 = dim1 + dim2
    h1 = (dim4 - dim3) / 2.

    w2 = dim2
    h2 = dim3

    A = 2. * (h1 * w1) + h2 * w2
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def l_section(class_type, beam_type, dim, prop):
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
    dim1, dim2, dim3, dim4 = dim
    t1 = dim3
    t2 = dim4
    b = dim1 - 0.5 * t2
    h = dim2 - 0.5 * t1
    h2 = dim2 - t1
    b1 = dim1 - t2
    A = (b + 0.5 * t2) * t1 + h2 * t2

    yc = t2*h2 * (h2 + t1) / (2 * A)
    zc = t1*b1 * (b1 + t2) / (2 * A)
    i1 = (
        t1 ** 3 * (b + 0.5 * t2) / 12. +
        t1 * (b + 0.5 * t2) * yc ** 2 +
        t2 * h ** 3 / 12 +
        h2 * t2 * (0.5 * (h2 + t1) - yc) ** 2
    )
    i2 = (
        t2 ** 3 * h2 / 12. +
        t1*(b + 0.5 * t2) ** 3 / 12. +
        t1*(b + 0.5 * t2) * (0.5 * b1 - zc) ** 2  # zc is z2 in the docs...
    )
    i12 = (
        zc * yc * t1 * t2
        - b1 * t1 * yc * (0.5 * (b1 + t2) - zc)
        - h2 * t2 * zc * (0.5 * (h2 + t1) - yc)
    )
    #j = 1/3 * (t1**3*b + t2**3 * h)
    #k1 = h2 * t2 / A
    #k2 = b1 * t1 / A
    #yna = yc
    #zna = zc
    return A, i1, i2, i12

def chan_section(class_name, beam_type, dim, prop):
    #
    # +-+----+   ^  ^
    # | |    |   |  | d4 = tf
    # | +----+   |  v
    # | |        |         ^
    # | |<-- d3  |         |
    # | |    tw  | d2 = h  | hw
    # | |        |         |
    # | +----+   |         v
    # | |    |   |
    # +-+----+   v
    #   <-bf->
    # <--d1-->
    #

    dim1, dim2, dim3, dim4 = dim
    tw = tweb = dim3
    tf = tflange = dim4
    bflange_tweb = dim1

    bflange = bflange_tweb - tweb
    hweb_all = dim2
    hweb = dim2 - 2. * tf
    #print(f'tflange={tflange} bflage={bflange} hweb_all={hweb_all} hweb={hweb} tweb={tweb}')
    A0 = 2 * (tflange * bflange) + tweb * hweb_all

    d = dim2 - 2 * tf

    # per https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf
    #b = dim1 - 0.5 * tw
    #h = dim2 - tf
    #bf = dim1 - tw
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
    A1 = 2 * tf * dim1 + hweb * tweb
    A = A1
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def chan1_section(class_name, beam_type, dim, prop):
    #
    # +--------+  ^
    # |        |  |
    # +---+----+  |
    # |   |       |         ^
    # |   |       |         |
    # |   |       | d4 = h  | d3 = hw = h2
    # |   |       |         |
    # +---+----+  |         v
    # |        |  |
    # +--------+  v
    # <---w1--->
    # <--><---->
    # dim2 dim1
    # tw    bf
    dim1, dim2, dim3, dim4 = dim
    tweb = dim2
    bflange = dim1
    h_all = dim4
    h_inner = dim3

    h1 = h_all - h_inner  # 2*tflange
    assert h_all > h_inner, f'h_all(dim4)={h_all} h_inner(dim3)={h_inner}'
    tflange = h1 / 2.
    w_all = bflange + tweb  # bf + tw
    A = h1 * w_all + h_inner * tweb
    #print(f'h_all={h_all} h_inner={h_inner} tweb={tweb} tflange={tflange} bflange={bflange}')
    A2 = dim2 * dim3 + (h_all - dim3) * (dim1 + dim2)
    assert np.allclose(A, A2), 'A=%s A2=%s' % (A, A2)
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def chan2_section(class_name, beam_type, dim, prop):
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
    I1 = I2 = I12 = None
    return A, I1, I2, I12

def dbox_section(class_name, beam_type, dim, prop):
    #DBOX = ((DIM2*DIM3)-((DIM2-DIM7-DIM8)*(DIM3-((0.5*DIM5)+DIM4)))) +
    #       (((DIM1-DIM3)*DIM2)-((DIM2-(DIM9+DIM10))*(DIM1-DIM3-(0.5*DIM5)-DIM6)))
    #
    #  |--2--|---5---|
    #  |     |       |
    #  1     3       6
    #  |     |       |
    #  |--4--|---7---|
    #
    dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, dim9, dim10 = dim

    #0,1,2,6,11
    #1,2,3,7,12

    htotal = dim2
    wtotal = dim1

    wleft = dim3 - dim4 - dim5 / 2.
    wright = dim1 - dim3 - dim6 - dim5 / 2.
    assert wleft > 0, wleft
    #assert wright > 0, wright

    #A1 = dim4 * htotal
    #A2 = dim7 * wleft
    #A3 = dim5 * htotal
    #A4 = dim8 * wleft
    #A5 = dim9 * wright
    #A6 = dim6 * htotal
    #A7 = dim10 * wright
    #A = A1 + A2 + A3 + A4 + A5 + A6 + A7

    A1 = dim4 * htotal
    A2 = wleft * dim7
    A3 = dim5 * htotal
    A4 = wleft * dim8
    A5 = wright * dim9
    A6 = dim6 * htotal
    A7 = wright * dim10
    A = A1 + A2 + A3 + A4 + A5 + A6 + A7

    z1 = dim4 / 2
    z2 = dim4 + wleft / 2
    z3 = dim3
    z4 = z2
    z5 = dim3 + dim5 / 2 + wright / 2
    z6 = dim1 - dim6 / 2
    z7 = z5

    y1 = htotal / 2
    y2 = htotal - dim7 / 2
    y3 = y1
    y4 = dim8 / 2
    y5 = htotal - dim9 / 2.
    y6 = y1
    y7 = dim10 / 2

    cg_y = (A1 * y1 + A2 * y2 + A3 * y3 + A4 * y4 + A5 * y5 + A6 * y6 + A7 * y7) / A
    cg_z = (A1 * z1 + A2 * z2 + A3 * z3 + A4 * z4 + A5 * z5 + A6 * z6 + A7 * z7) / A

    dy1 = y1 - cg_y
    dy2 = y2 - cg_y
    dy3 = y3 - cg_y
    dy4 = y4 - cg_y
    dy5 = y5 - cg_y
    dy6 = y6 - cg_y
    dy7 = y7 - cg_y

    dz1 = z1 - cg_z
    dz2 = z2 - cg_z
    dz3 = z3 - cg_z
    dz4 = z4 - cg_z
    dz5 = z5 - cg_z
    dz6 = z6 - cg_z
    dz7 = z7 - cg_z
    assert dz1 < 0., dz1
    assert dy2 > 0., dy2

    # inertia taken about the centroid
    # 1/12 * b * h^3 = 1/12 * A * h^2
    Iy1 = A1 * htotal ** 2
    Iy2 = A2 * dim7 ** 2
    Iy3 = A3 * htotal ** 2
    Iy4 = A4 * dim8 ** 2
    Iy5 = A5 * dim9 ** 2
    Iy6 = A6 * htotal ** 2
    Iy7 = A7 * dim10 ** 2
    dI1 = A1 * dy1 + A2 * dy2 + A3 * dy3 + A4 * dy4 + A5 * dy5 + A6 * dy6 + A7 * dy7
    I1 = (Iy1 + Iy2 + Iy3 + Iy4 + Iy5 + Iy6 + Iy7) / 12. + dI1

    # 1/12 * b^3 * h = 1/12 * A * b^2
    Iz1 = A1 * dim4 ** 2
    Iz2 = A2 * wleft ** 2
    Iz3 = A3 * dim5 ** 2
    Iz4 = A4 * wleft ** 2
    Iz5 = A5 * wright ** 2
    Iz6 = A6 * dim6 ** 2
    Iz7 = A7 * wright ** 2
    dI2 = A1 * dz1 + A2 * dz2 + A3 * dz3 + A4 * dz4 + A5 * dz5 + A6 * dz6 + A7 * dz7
    I2 = (Iz1 + Iz2 + Iz3 + Iz4 + Iz5 + Iz6 + Iz7) / 12. + dI2

    I12 = None
    return A, I1, I2, I12


class IntegratedLineProperty(LineProperty):
    def __init__(self):
        self.xxb = None
        self.A = None
        self.j = None
        self.i1 = None
        self.i2 = None
        self.i12 = None
        LineProperty.__init__(self)

    def Area(self) -> float:
        A = integrate_positive_unit_line(self.xxb, self.A)
        return A

    def J(self) -> float:
        J = integrate_positive_unit_line(self.xxb, self.j)
        return J

    def I11(self) -> float:
        i1 = integrate_positive_unit_line(self.xxb, self.i1)
        return i1

    def I22(self) -> float:
        i2 = integrate_positive_unit_line(self.xxb, self.i2)
        return i2

    def I12(self) -> float:
        i12 = integrate_unit_line(self.xxb, self.i12)
        return i12

    def Nsm(self) -> float:
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

    @classmethod
    def export_to_hdf5(cls, h5_file, model, pids):
        """exports the properties in a vectorized way"""
        #comments = []
        mids = []
        A = []
        J = []
        I = []

        c = []
        d = []
        e = []
        f = []
        k = []
        nsm = []
        for pid in pids:
            prop = model.properties[pid]
            #comments.append(prop.comment)
            mids.append(prop.mid)
            A.append(prop.A)
            I.append([prop.i1, prop.i2, prop.i12])
            J.append(prop.j)

            c.append([prop.c1, prop.c2])
            d.append([prop.d1, prop.d2])
            e.append([prop.e1, prop.e2])
            f.append([prop.f1, prop.f2])

            ki = []
            if prop.k1 is None:
                ki.append(np.nan)
            else:
                ki.append(prop.k1)
            if prop.k2 is None:
                ki.append(np.nan)
            else:
                ki.append(prop.k2)

            k.append(ki)
            nsm.append(prop.nsm)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('mid', data=mids)
        h5_file.create_dataset('A', data=A)
        h5_file.create_dataset('J', data=J)
        h5_file.create_dataset('I', data=I)
        h5_file.create_dataset('c', data=c)
        h5_file.create_dataset('d', data=d)
        h5_file.create_dataset('e', data=e)
        h5_file.create_dataset('f', data=f)
        h5_file.create_dataset('k', data=k)
        h5_file.create_dataset('nsm', data=nsm)
        #h5_file.create_dataset('_comment', data=comments)

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

        return PBAR(pid, mid, A=A, i1=i1, i2=i2, i12=i12, j=j, nsm=nsm,
                    c1=c1, c2=c2, d1=d1, d2=d2, e1=e1, e2=e2,
                    f1=f1, f2=f2, k1=k1, k2=k2, comment=comment)

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

    def MassPerLength(self) -> float:
        r"""
        Gets the mass per length :math:`\frac{m}{L}` of the CBAR.

        .. math:: \frac{m}{L} = \rho A + nsm

        """
        A = self.Area()
        rho = self.Rho()
        nsm = self.Nsm()
        return rho * A + nsm

    def get_cdef(self):
        cdef = np.array([
                [self.c1, self.c2],
                [self.d1, self.d2],
                [self.e1, self.e2],
                [self.f1, self.f2],
            ], dtype='float64')
        return cdef

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PBAR mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def Area(self) -> float:
        """Gets the area :math:`A` of the CBAR."""
        return self.A

    #def Nsm(self) -> float:
    #    return self.nsm

    #def J(self) -> float:
       #return self.j

    def I11(self) -> float:
        """gets the section I11 moment of inertia"""
        return self.i1

    def I22(self) -> float:
        """gets the section I22 moment of inertia"""
        return self.i2

    def I12(self) -> float:
        """gets the section I12 moment of inertia"""
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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
    _properties = ['Type', 'valid_types']
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

        # approximate
        #'I', 'CHAN', 'T', 'CHAN1', 'T1', 'CHAN2', 'T2', 'L' and 'BOX1'.
        'L' : 4,
    }  # for GROUP="MSCBML0"
    #pname_fid_map = {
        #12 : 'DIM1',
        #13 : 'DIM2',
        #14 : 'DIM3',
        #15 : 'DIM3',
    #}
    def update_by_pname_fid(self, pname_fid, value):
        if isinstance(pname_fid, str) and pname_fid.startswith('DIM'):
            num = int(pname_fid[3:])
            self.dim[num - 1] = value
        else:
            raise NotImplementedError('PBARL Type=%r name=%r has not been implemented' % (
                self.Type, pname_fid))

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        Type = 'ROD'
        dim = [1.]
        return PBARL(pid, mid, Type, dim, group='MSCBML0', nsm=0., comment='')

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
            keys = list(self.valid_types.keys())
            msg = ('Invalid PBARL Type, Type=%s '
                   'valid_types=%s' % (self.Type, ', '.join(sorted(keys))))
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
        if not isinstance(self.group, str):
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

        try:
            ndim = cls.valid_types[Type]
        except KeyError:
            keys = list(cls.valid_types.keys())
            raise KeyError('%r is not a valid PBARL type\nallowed_types={%s}' % (
                Type, ', '.join(sorted(keys))))

        # PBARL
        # 9. For DBOX section, the default value for DIM5 to DIM10 are
        #    based on the following rules:
        #     a. DIM5, DIM6, DIM7 and DIM8 have a default value of
        #        DIM4 if not provided.
        #     b. DIM9 and DIM10 have a default value of DIM6 if not
        #        provided.
        dim = []
        if Type == 'DBOX':
            for i in range(ndim):
                if i in [4, 5, 6, 7]:
                    dim4 = dim[3]
                    dimi = double_or_blank(card, 9 + i, 'ndim=%s; dim%i' % (ndim, i + 1),
                                           default=dim4)
                elif i in [8, 9]:
                    dim6 = dim[5]
                    dimi = double_or_blank(card, 9 + i, 'ndim=%s; dim%i' % (ndim, i + 1),
                                           default=dim6)
                else:
                    dimi = double(card, 9 + i, 'ndim=%s; dim%i' % (ndim, i + 1))
                dim.append(dimi)
        else:
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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PBARL mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
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


    def get_cdef(self):
        r"""
        these axes are backwards...

            ^ y
        +---|---+
        |   |   |
        |   +------> z
        |       |
        +-------+
        """
        if self.Type in ['ROD', 'TUBE']:
            R = self.dims[0]
            cdef = np.array([
                [R, 0.],
                [0., R],
                [-R, 0.],
                [0., -R],
            ], dtype='float64')
        elif self.Type in ['BAR', 'BOX']:
            height, width = self.dims
            cdef = np.array([
                [width, height],
                [width, -height],
                [-width, -height],
                [-width, height],
            ], dtype='float64') / 2.
        else:
            raise NotImplementedError(self.Type)
        return cdef

    @property
    def Type(self):
        """gets Type"""
        return self.beam_type
    @Type.setter
    def Type(self, beam_type):
        """sets Type"""
        self.beam_type = beam_type

    def Area(self) -> float:
        """Gets the area :math:`A` of the CBAR."""
        return _bar_areaL('PBARL', self.beam_type, self.dim, self)

    def Nsm(self) -> float:
        """Gets the non-structural mass :math:`nsm` of the CBAR."""
        return self.nsm

    def MassPerLength(self) -> float:
        r"""
        Gets the mass per length :math:`\frac{m}{L}` of the CBAR.

        .. math:: \frac{m}{L} = A \rho + nsm
        """
        rho = self.Rho()
        area = self.Area()
        nsm = self.Nsm()
        return area * rho + nsm

    def I1(self) -> float:
        """gets the section I1 moment of inertia"""
        I1 = A_I1_I2_I12(self, self.beam_type, self.dim)
        return I1

    def I2(self) -> float:
        """gets the section I2 moment of inertia"""
        I2 = A_I1_I2_I12(self, self.beam_type, self.dim)[2]
        return I2

    def I12(self) -> float:
        """gets the section I12 moment of inertia"""
        try:
            I12 = A_I1_I2_I12(self, self.beam_type, self.dim)
        except:
            print(str(self))
            raise
        return I12

    #def I1_I2_I12(self):
        #"""gets the section I1, I2, I12 moment of inertia"""
        #return I1_I2_I12(prop, prop.dim)

    def I11(self) -> float:
        return self.I1()

    def I22(self) -> float:
        return self.I2()

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
        elif beam_type == 'HEXA':
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

        elif beam_type == 'I':
            (d1, d2, d3) = dim
            raise NotImplementedError('PBARL beam_type=%r' % beam_type)

        elif beam_type == 'H':
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

        elif beam_type == 'T2':
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
            msg = '_points for beam_type=%r dim=%r on PBARL is not supported' % (beam_type, self.dim)
            raise NotImplementedError(msg)
        return array(points), Area

    def J(self) -> float:
        beam_type = self.beam_type
        if beam_type == 'ROD':
            A, I1, I2, I12 = rod_section(self.type, beam_type, self.dim, self)
        elif beam_type == 'TUBE':
            A, I1, I2, I12 = tube_section(self.type, beam_type, self.dim, self)
        elif beam_type == 'TUBE2':
            A, I1, I2, I12 = tube2_section(self.type, beam_type, self.dim, self)
        elif beam_type == 'BOX':
            A, I1, I2, I12 = box_section(self.type, beam_type, self.dim, self)

        #elif beam_type in ['BAR']:
            #assert len(self.dim) == 2, 'dim=%r' % self.dim
            #b, h = self.dim
            #(A, Ix, Iy, Ixy) = self.A_I1_I2_I12()
            #J = Ix + Iy
        elif beam_type in ['BAR', 'CROSS', 'HEXA', 'T2', 'H']:
            points, unused_Area = self._points(beam_type, self.dim)
            yi = points[0, :-1]
            yip1 = points[0, 1:]

            xi = points[1, :-1]
            xip1 = points[1, 1:]

            #: .. seealso:: http://en.wikipedia.org/wiki/Area_moment_of_inertia
            ai = xi*yip1 - xip1*yi
            I1 = 1/12 * sum((yi**2 + yi*yip1+yip1**2)*ai)
            I2 = 1/12 * sum((xi**2 + xi*xip1+xip1**2)*ai)
            #I12 = 1/24*sum((xi*yip1 + 2*xi*yi + 2*xip1*yip1 + xip1*yi)*ai)
        elif beam_type == 'I':
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
                    beam_type, self.dim)
                raise NotImplementedError(msg)
            if s1 != s2:
                msg = 'J for beam_type=%r dim=%r on PBARL s1 != s2 is not supported' % (
                    beam_type, self.dim)
                raise NotImplementedError(msg)
            h = d - b1 - b2
            s = s1
            b = b1
            I1 = (b*d**3-h**3*(b-t)) / 12.
            I2 = (2.*s*b**3 + h*t**3) / 12.
        #elif beam_type == 'T': # test
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
                    #beam_type, self.dim)
                #raise NotImplementedError(msg)
            #if s1 != s2:
                #msg = 'J for beam_type=%r dim=%r on PBARL s1 != s2 is not supported' % (
                    #beam_type, self.dim)
                #raise NotImplementedError(msg)
            #h = d - b1 - b2
            #s = s1
            #b = b1

            # http://www.engineersedge.com/material_science/moment-inertia-gyration-6.htm
            #y = d**2*t+s**2*(b-t)/(2*(b*s+h*t))
            #I1 = (t*y**3 + b*(d-y)**3 - (b-t)*(d-y-s)**3)/3.
            #I2 = t**3*(h-s)/12. + b**3*s/12.
            #A = b*s + h*t

        elif beam_type == 'C':
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
            I1 = (b * d**3 - h **3 * (b-t)) / 12.
            #I12 = (2.*s*b**3 + h*t**3)/3 - A*cx**2
        else:  # pragma: no cover
            msg = 'J for beam_type=%r dim=%r on PBARL is not supported' % (beam_type, self.dim)
            raise NotImplementedError(msg)

        #: .. seealso:: http://en.wikipedia.org/wiki/Perpendicular_axis_theorem
        J = I1 + I2
        return J

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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBRSECT(LineProperty):
    """
    not done
    """
    type = 'PBRSECT'

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 2
        form = 'FORM'
        options = [('OUTP', 10)]
        return PBRSECT(pid, mid, form, options, comment='')

    def _finalize_hdf5(self, encoding):
        self.brps = {key : value for key, value in zip(*self.brps)}
        self.ts = {key : value for key, value in zip(*self.ts)}
        self.inps = {key : value for key, value in zip(*self.inps)}

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

        # int : int
        self.brps = {}
        self.inps = {}

        # int : floats
        self.ts = {}
        assert isinstance(options, list), options
        for key_value in options:
            try:
                key, value = key_value
            except ValueError:
                print(key_value)
                raise
            key = key.upper()

            if key == 'NSM':
                self.nsm = float(value)
            elif 'INP' in key:
                if key.startswith('INP('):
                    assert key.endswith(')'), 'key=%r' % key
                    key_id = int(key[4:-1])
                    self.inps[key_id] = int(value)
                else:
                    self.inps[0] = int(value)
            elif key == 'OUTP':
                self.outp = int(value)

            elif key.startswith('BRP'):
                if key.startswith('BRP('):
                    assert key.endswith(')'), 'key=%r' % key
                    key_id = int(key[4:-1])
                    self.brps[key_id] = int(value)
                else:
                    self.brps[0] = int(value)

            elif key.startswith('T('):
                index, out = split_arbitrary_thickness_section(key, value)
                self.ts[index] = out
            elif key == 'T':
                self.ts[1] = float(value)

            #if key == 'NSM':
                #self.nsm = float(value)
            #elif key == 'OUTP':
                #self.outp = int(value)
            #elif key == 'BRP(1)':
                #self.brp1 = int(value)
            #elif key == 'T':
                #self.t = float(value)
            else:
                raise NotImplementedError('PBRSECT.pid=%s key=%r value=%r' % (pid, key, value))

        self.mid_ref = None
        self.brps_ref = {}
        self.outp_ref = None

    def validate(self):
        assert self.form in ['GS', 'OP', 'CP'], 'pid=%s form=%r' % (self.pid, self.form)

        #assert self.outp is not None, 'form=%s outp=%s' % (self.form, self.outp)
        if self.form == 'GS':
            assert len(self.inps) > 0, 'form=%s inps=%s' % (self.form, self.inps)
            assert len(self.brps) == 0, 'form=%s brps=%s' % (self.form, self.brps)
            assert len(self.ts) == 0, 'form=%s ts=%s' % (self.form, self.ts)
        elif self.form in ['OP', 'CP']:
            assert len(self.inps) == 0, 'form=%s inps=%s' % (self.form, self.inps)
            assert len(self.brps) >= 0, 'form=%s brps=%s' % (self.form, self.brps)
            assert len(self.ts) >= 0, 'form=%s ts=%s' % (self.form, self.ts)

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
        lines_joined = ','.join(card[1:]).replace(' ', '').replace(',,', ',')

        if lines_joined:
            fields = get_beam_sections(lines_joined)
            options = [field.split('=', 1) for field in fields]
            #C:\MSC.Software\MSC.Nastran\msc20051\nast\tpl\zbr3.dat
            #options = [
                #[u'OUTP', u'201'],
                #[u'T', u'1.0'],
                #[u'BRP', u'202'],
                #[u'T(11)', u'[1.2'],
                #[u'PT', u'(202'], [u'224)]'],
                #[u'T(12)', u'[1.2'],
                #[u'PT', u'(224'],
                #[u'205)]'],
            #]
        else:
            options = []

        pid = integer(bdf_card, 1, 'pid')
        mid = integer(bdf_card, 2, 'mid')
        form = string_or_blank(bdf_card, 3, 'form')

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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PBMSECT mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        if self.outp is not None:
            self.outp_ref = model.Set(self.outp)
            self.outp_ref.cross_reference_set(model, 'Point', msg=msg)

        if len(self.brps):
            for key, brpi in self.brps.items():
                brpi_ref = model.Set(brpi, msg=msg)
                brpi_ref.cross_reference_set(model, 'Point', msg=msg)
                self.brps_ref[key] = brpi_ref

    def plot(self, model, figure_id=1, show=False):
        """
        Plots the beam section

        Parameters
        ----------
        model : BDF()
            the BDF object
        figure_id : int; default=1
            the figure id
        show : bool; default=False
            show the figure when done

        """
        class_name = self.__class__.__name__
        form_map = {
            'GS' : 'General Section',
            'OP' : 'Open Profile',
            'CP' : 'Closed Profile',
        }
        formi = ' form=%s' % form_map[self.form]
        plot_arbitrary_section(
            model, self,
            self.inps, self.ts, self.brps_ref, self.nsm, self.outp_ref,
            figure_id=figure_id,
            title=class_name + ' pid=%s' % self.pid + formi,
            show=show)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None
        self.outp_ref = None
        self.brps_ref = {}

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

    def Area(self) -> float:
        """Gets the area :math:`A` of the CBAR."""
        return 0.
        #raise NotImplementedError('Area is not implemented for PBRSECT')

    def Nsm(self) -> float:
        """Gets the non-structural mass :math:`nsm` of the CBAR."""
        return 0.
        #raise NotImplementedError('Nsm is not implemented for PBRSECT')

    def MassPerLength(self) -> float:
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        end = write_arbitrary_beam_section(self.inps, self.ts, self.brps, self.nsm, self.outp)
        out = self.comment + print_card_8(card) + end
        return out

    def __repr__(self):
        return self.write_card()


class PBEAM3(LineProperty):  # not done, cleanup; MSC specific card
    """
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |    1   |     2    |    3    |    4    |     5    |    6    |    7    |     8    |     9    |
    +========+==========+=========+=========+==========+=========+=========+==========+==========+
    | PBEAM3 |    PID   |   MID   |   A(A)  |   IZ(A)  |  IY(A)  |  IYZ(A) |   J(A)   |  NSM(A)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |   CY(A)  |  CZ(A)  |  DY(A)  |   DZ(A)  |  EY(A)  |  EZ(A)  |   FY(A)  |   FZ(A)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |   SO(B)  |         |   A(B)  |   IZ(B)  |  IY(B)  |  IYZ(B) |   J(B)   |  NSM(B)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |   CY(B)  |  CZ(B)  |  DY(B)  |   DZ(B)  |  EY(B)  |  EZ(B)  |   FY(B)  |   FZ(B)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |   SO(C)  |         |   A(C)  |   IZ(C)  |  IY(C)  |  IYZ(C) |   J(C)   |  NSM(C)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |   CY(C)  |  CZ(C)  |  DY(C)  |   DZ(C)  |  EY(C)  |  EZ(C)  |   FY(C)  |   FZ(C)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |    KY    |   KZ    |  NY(A)  |   NZ(A)  |  NY(B)  |  NZ(B)  |   NY(C)  |   NZ(C)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |  MY(A)   |  MZ(A)  |  MY(B)  |   MZ(B)  |  MY(C)  |  MZ(C)  |  NSIY(A) |  NSIZ(A) |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        | NSIYZ(A) | NSIY(B) | NSIZ(B) | NSIYZ(B) | NSIY(C) | NSIZ(C) | NSIYZ(C) |   CW(A)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |   CW(B)  |  CW(C)  |  STRESS |          |         |         |          |          |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |   WC(A)  |  WYC(A) |  WZC(A) |   WD(A)  |  WYD(A) |  WZD(A) |   WE(A)  |  WYE(A)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |  WZE(A)  |  WF(A)  |  WYF(A) |  WZF(A)  |  WC(B)  |  WYC(B) |  WZC(B)  |   WD(B)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |  WYD(B)  |  WZD(B) |  WE(B)  |  WYE(B)  |  WZE(B) |  WF(B)  |  WYF(B)  |  WZF(B)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |   WC(C)  |  WYC(C) |  WZC(C) |   WD(C)  |  WYD(C) |  WZD(C) |  WE(C)   |  WYE(C)  |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+
    |        |  WZE(C)  |  WF(C)  |  WYF(C) |  WZF(C)  |         |         |          |          |
    +--------+----------+---------+---------+----------+---------+---------+----------+----------+

    """
    type = 'PBEAM3'

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 2
        A = 3.
        iz = 4.
        iy = 5.
        return PBEAM3(pid, mid, A, iz, iy,
                      iyz=None, j=None, nsm=0.,
                      so=None, cy=None, cz=None, dy=None, dz=None,
                      ey=None, ez=None, fy=None, fz=None,
                      ky=1., kz=1.,
                      ny=None, nz=None, my=None, mz=None,
                      nsiy=None, nsiz=None, nsiyz=None,
                      cw=None, stress='GRID',
                      w=None, wy=None, wz=None, comment='')

    def __init__(self, pid, mid, A, iz, iy, iyz=None, j=None, nsm=0.,
                 so=None,
                 cy=None, cz=None, dy=None, dz=None, ey=None, ez=None, fy=None, fz=None,
                 ky=1., kz=1.,
                 ny=None, nz=None, my=None, mz=None,
                 nsiy=None, nsiz=None, nsiyz=None,
                 cw=None, stress='GRID',
                 w=None, wy=None, wz=None,
                 comment=''):
        """
        Creates a PBEAM3 card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        A : List[float]
            areas for ABC
        iz / iy / iyz : List[float]
            area moment of inertias for ABC
        iyz : List[float]; default=None -> [0., 0., 0.]
            area moment of inertias for ABC
        j : List[float]; default=None
            polar moment of inertias for ABC
            None -> iy + iz from section A for ABC
        so : List[str]; default=None
            None -> ['YES', 'YESA', 'YESA']
        cy / cz / dy / dz / ey / ez / fy / fz : List[float]; default=[0., 0., 0.]
            stress recovery loctions for ABC
        ny / nz : List[float]
            Local (y, z) coordinates of neutral axis for ABC
        my / mz : List[float]
            Local (y, z) coordinates of nonstructural mass center of gravity for ABC
        nsiy / nsiz / nsiyz : List[float]
            Nonstructural mass moments of inertia per unit length about
            local y and z-axes, respectively, with regard to the nonstructural mass
            center of gravity for ABC
        cw : List[float]
            warping coefficients for ABC
        stress : str; default='GRID'
            Location selection for stress, strain and force output.
        w : (4, 3) float numpy array; default=None
            Values of warping function at stress recovery points
            None : array of 0.0
        wy / wz : (4, 3) float numpy array; default=None
            Gradients of warping function in the local (y, z) coordinate
            system at stress recovery points
            None : array of 0.0

        """
        LineProperty.__init__(self)
        if comment:
            self.comment = comment

        def _pbeam3_station_a(values):
            """used by A, Iz, Iy"""
            if isinstance(values, float_types):
                values = [values] * 3
            station_a_value = values[0]
            for i, value in zip(count(), values[1:]):
                if value is None:
                    values[i + 1] = station_a_value
            return values

        def _pbeam3_station_a_default(values, default_value):
            """used by Iyz, nsm"""
            if values is None:
                values = [default_value] * 3
            elif isinstance(values, float_types):
                values = [values] * 3

            station_a_value = default_value if values[0] is None else values[0]
            for i, value in zip(count(), values[1:]):
                if value is None:
                    values[i + 1] = station_a_value
            return values

        self.pid = pid
        self.mid = mid

        self.A = _pbeam3_station_a(A)
        self.iz = _pbeam3_station_a(iz)
        self.iy = _pbeam3_station_a(iy)
        self.iyz = _pbeam3_station_a_default(iyz, 0.0)
        self.j = _pbeam3_station_a_default(j, self.iy[0] + self.iz[0])
        self.nsm = _pbeam3_station_a_default(nsm, 0.0)

        if so is None:
            so = ['YES', 'YESA', 'YESA']
        elif isinstance(so, str):
            so = [so] * 3
        self.so = so

        def _pbeam3_default_list(values, default):
            """used by Cy/Cz, Dy/Dz, Ey/Ez, Fy/Fz"""
            if values is None:
                values = [default] * 3
            elif isinstance(values, float_types):
                values = [values] * 3
            for i, value in zip(count(), values):
                if value is None:
                    values[i] = default
            return values

        self.cy = _pbeam3_default_list(cy, 0.)
        self.cz = _pbeam3_default_list(cz, 0.)
        self.dy = _pbeam3_default_list(dy, 0.)
        self.dz = _pbeam3_default_list(dz, 0.)
        self.ey = _pbeam3_default_list(ey, 0.)
        self.ez = _pbeam3_default_list(ez, 0.)
        self.fy = _pbeam3_default_list(fy, 0.)
        self.fz = _pbeam3_default_list(fz, 0.)

        self.ky = ky
        self.kz = kz

        def _pbeam3_default_list_station_a(values, default):
            """used by Ny/Nz, My/Mz, NSIy/NSIz/NSIyz, Cw"""
            if values is None:
                values = [default] * 3
            elif isinstance(values, float_types):
                values = [values] * 3
            if values[0] is None:
                values[0] = default

            station_a_value = values[0]
            for i, value in zip(count(), values[1:]):
                if value is None:
                    values[i + 1] = station_a_value
            return values

        self.ny = _pbeam3_default_list_station_a(ny, 0.)
        self.nz = _pbeam3_default_list_station_a(nz, 0.)
        self.my = _pbeam3_default_list_station_a(my, 0.)
        self.mz = _pbeam3_default_list_station_a(mz, 0.)

        self.nsiy = _pbeam3_default_list_station_a(nsiy, 0.)
        self.nsiz = _pbeam3_default_list_station_a(nsiz, 0.)
        self.nsiyz = _pbeam3_default_list_station_a(nsiyz, 0.)

        self.cw = _pbeam3_default_list_station_a(cw, 0.)
        self.stress = stress

        if w is None:
            w = np.zeros((3, 4), dtype='float64')
        if wy is None:
            wy = np.zeros((3, 4), dtype='float64')
        if wz is None:
            wz = np.zeros((3, 4), dtype='float64')

        self.w = w
        self.wy = wy
        self.wz = wz
        self.mid_ref = None

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
        #PID MID A(A) IZ(A) IY(A) IYZ(A) J(A) NSM(A)
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')

        area = [double(card, 3, 'A')]
        iz = [double(card, 4, 'Iz')]
        iy = [double(card, 5, 'Iy')]
        iyz = [double_or_blank(card, 6, 'Iyz', 0.0)]
        j = [double_or_blank(card, 7, 'J', iy[0] + iz[0])]
        nsm = [double_or_blank(card, 8, 'nsm', 0.0)]

        #CY(A) CZ(A) DY(A) DZ(A) EY(A) EZ(A) FY(A) FZ(A)
        cy = [double_or_blank(card, 9, 'cy', default=0.)]
        cz = [double_or_blank(card, 10, 'cz', default=0.)]

        dy = [double_or_blank(card, 11, 'dy', default=0.)]
        dz = [double_or_blank(card, 12, 'dz', default=0.)]

        ey = [double_or_blank(card, 13, 'ey', default=0.)]
        ez = [double_or_blank(card, 14, 'ez', default=0.)]

        fy = [double_or_blank(card, 15, 'fy', default=0.)]
        fz = [double_or_blank(card, 16, 'fz', default=0.)]

        #SO(B)        A(B) IZ(B) IY(B) IYZ(B)  J(B) NSM(B)
        #CY(B) CZ(B) DY(B) DZ(B) EY(B)  EZ(B) FY(B)  FZ(B)

        #SO(C)        A(C) IZ(C) IY(C) IYZ(C) J(C)  NSM(C)
        #CY(C) CZ(C) DY(C) DZ(C) EY(C)  EZ(C) FY(C)  FZ(C)

        so = ['YES']
        locations = ['B', 'C']
        for i, location in enumerate(locations):
            offset = 17 + i * 16
            so.append(string_or_blank(card, offset, 'SO_%s' % location, default='YESA'))

            area.append(double_or_blank(card, offset + 2, 'area_%s' % location, default=area[0]))
            iz.append(double_or_blank(card, offset + 3, 'Iz', default=iz[0]))
            iy.append(double_or_blank(card, offset + 4, 'Iy', default=iy[0]))
            iyz.append(double_or_blank(card, offset + 5, 'Iyz', default=iyz[0]))
            j.append(double_or_blank(card, offset + 6, 'J', default=j[0]))
            nsm.append(double_or_blank(card, offset + 7, 'nsm', default=nsm[0]))

            cy.append(double_or_blank(card, offset + 8, 'cy', default=0.))
            cz.append(double_or_blank(card, offset + 9, 'cz', default=0.))

            dy.append(double_or_blank(card, offset + 10, 'dy', default=0.))
            dz.append(double_or_blank(card, offset + 11, 'dz', default=0.))

            ey.append(double_or_blank(card, offset + 12, 'ey', default=0.))
            ez.append(double_or_blank(card, offset + 13, 'ez', default=0.))

            fy.append(double_or_blank(card, offset + 14, 'fy', default=0.))
            fz.append(double_or_blank(card, offset + 15, 'fz', default=0.))

        #KY       KZ      NY(A)   NZ(A)    NY(B)   NZ(B)   NY(C)    NZ(C)
        #MY(A)    MZ(A)   MY(B)   MZ(B)    MY(C)   MZ(C)   NSIY(A)  NSIZ(A)
        #NSIYZ(A) NSIY(B) NSIZ(B) NSIYZ(B) NSIY(C) NSIZ(C) NSIYZ(C) CW(A)
        #CW(B)    CW(C)   STRESS

        ifield = 49
        ky = double_or_blank(card, ifield, 'Ky', default=1.0)
        kz = double_or_blank(card, ifield + 1, 'Kz', default=1.0)
        ifield += 2

        locations = ['A', 'B', 'C']
        ny = []
        nz = []
        for i, location in enumerate(locations):
            if i == 0:
                nyi = double_or_blank(card, ifield, 'NY(%s)' % location, default=0.0)
                nzi = double_or_blank(card, ifield + 1, 'NZ(%s)' % location, default=0.0)
            else:
                nyi = double_or_blank(card, ifield, 'NY(%s)' % location, default=ny[0])
                nzi = double_or_blank(card, ifield + 1, 'NZ(%s)' % location, default=nz[0])
            ny.append(nyi)
            nz.append(nzi)
            ifield += 2

        my = []
        mz = []
        for i, location in enumerate(locations):
            if i == 0:
                myi = double_or_blank(card, ifield, 'MY(%s)' % location, default=0.0)
                mzi = double_or_blank(card, ifield + 1, 'MZ(%s)' % location, default=0.0)
            else:
                myi = double_or_blank(card, ifield, 'MY(%s)' % location, default=my[0])
                mzi = double_or_blank(card, ifield + 1, 'MZ(%s)' % location, default=mz[0])
            my.append(myi)
            mz.append(mzi)
            ifield += 2

        nsiy = []
        nsiz = []
        nsiyz = []
        for i, location in enumerate(locations):
            if i == 0:
                nsiyi = double_or_blank(card, ifield, 'NSIY(%s)' % location, default=0.0)
                nsizi = double_or_blank(card, ifield + 1, 'NSIZ(%s)' % location, default=0.0)
                nsiyzi = double_or_blank(card, ifield + 2, 'NSIYZ(%s)' % location, default=0.0)
            else:
                nsiyi = double_or_blank(card, ifield, 'NSIY(%s)' % location, default=nsiy[0])
                nsizi = double_or_blank(card, ifield + 1, 'NSIZ(%s)' % location, default=nsiz[0])
                nsiyzi = double_or_blank(card, ifield + 2, 'NSIYZ(%s)' % location, default=nsiyz[0])
            nsiy.append(nsiyi)
            nsiz.append(nsizi)
            nsiyz.append(nsiyzi)
            ifield += 3

        cw = []
        for location in locations:
            cwi = double_or_blank(card, ifield, 'CW(%s)' % location, default=0.0)
            cw.append(cwi)
            ifield += 1
        stress = string_or_blank(card, ifield, 'STRESS', default='GRID')
        ifield += 6


        # WC(A)  WYC(A) WZC(A) WD(A)  WYD(A) WZD(A) WE(A)  WYE(A)
        # WZE(A) WF(A)) WYF(A) WZF(A) WC(B)  WYC(B) WZC(B) WD(B)
        # WYD(B) WZD(B) WE(B)  WYE(B) WZE(B) WF(B)  WYF(B) WZF(B)
        # WC(C)  WYC(C) WZC(C) WD(C)  WYD(C) WZD(C) WE(C)  WYE(C)
        # WZE(C) WF(C)  WYF(C) WZF(C)

        spots = ('C', 'D', 'E', 'F')
        w = np.zeros((3, 4), dtype='float64')
        wy = np.zeros((3, 4), dtype='float64')
        wz = np.zeros((3, 4), dtype='float64')
        for iloc, location in enumerate(locations):
            for ispot, spot in enumerate(spots):
                wi = double_or_blank(card, ifield, 'W%s(%s)' % (spot, location), default=0.0)
                wyi = double_or_blank(card, ifield + 1, 'WY%s(%s)' % (spot, location), default=0.0)
                wzi = double_or_blank(card, ifield + 2, 'WZ%s(%s)' % (spot, location), default=0.0)
                w[iloc, ispot] = wi
                wy[iloc, ispot] = wyi
                wz[iloc, ispot] = wzi
                ifield += 3

        return PBEAM3(pid, mid, area, iz, iy, iyz, j, nsm=nsm,
                      so=so,
                      cy=cy, cz=cz, dy=dy, dz=dz, ey=ey, ez=ez, fy=fy, fz=fz,
                      ky=ky, kz=kz,
                      ny=ny, nz=nz, my=my, mz=mz,
                      nsiy=nsiy, nsiz=nsiz, nsiyz=nsiyz,
                      cw=cw, stress=stress,
                      w=w, wy=wy, wz=wz,
                      comment=comment)

    def add_op2_data(self, data, comment=''):
        if comment:
            self.comment = comment
        raise NotImplementedError(data)

    def Nsm(self) -> List[float]:
        """
        Gets the non-structural mass :math:`nsm`.
        .. warning:: nsm field not supported fully on PBEAM3 card
        """
        return self.nsm

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PBEAM3 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def raw_fields(self):
        list_fields = ['PBEAM3', self.pid, self.Mid()]
        for (i, soi, ai, iz, iy, iyz, j, nsm, cy, cz, dy, dz, ey, ez, fy, fz) in zip(
                count(), self.so, self.A, self.iz, self.iy, self.iyz, self.j, self.nsm,
                self.cy, self.cz, self.dy, self.dz, self.ey, self.ez, self.fy, self.fz):
            if i == 0:
                list_fields += [
                    ai, iz, iy, iyz, j, nsm,
                    cy, cz, dy, dz, ey, ez, fy, fz]
            else:
                list_fields += [
                    soi, None,
                    ai, iz, iy, iyz, j, nsm,
                    cy, cz, dy, dz, ey, ez, fy, fz]

        #KY KZ NY(A) NZ(A) NY(B) NZ(B) NY(C) NZ(C)
        #MY(A) MZ(A) MY(B) MZ(B) MY(C) MZ(C) NSIY(A) NSIZ(A)
        #NSIYZ(A) NSIY(B) NSIZ(B) NSIYZ(B) NSIY(C) NSIZ(C) NSIYZ(C) CW(A)
        #CW(B) CW(C) STRESS
        list_fields += [self.ky, self.kz]
        for ny, nz in zip(self.ny, self.nz):
            list_fields += [ny, nz]
        for my, mz in zip(self.my, self.mz):
            list_fields += [my, mz]
        for nsiy, nsiz, nsiyz in zip(self.nsiy, self.nsiz, self.nsiyz):
            list_fields += [nsiy, nsiz, nsiyz]

        list_fields += self.cw
        list_fields += [self.stress, None, None, None, None, None]

        #WC(A) WYC(A) WZC(A) WD(A) WYD(A) WZD(A) WE(A) WYE(A)
        #WZE(A) WF(A)) WYF(A) WZF(A)
        for w, wy, wz in zip(self.w, self.wy, self.wz):
            for wi, wyi, wzi in zip(w, wy, wz):
                list_fields += [wi, wyi, wzi]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 1
        fsi = 1
        rm = 0.1
        t = 0.01
        return cls.add_beam_type_2(pid, mid,
                                   fsi, rm, t, p=None, rb=None, theta_b=None,
                                   nsm=0., rc=0., zc=0., comment='')

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
    def add_beam_type_1(cls, pid, mid,
                        A, i1, i2, j,
                        rb=None, theta_b=None,
                        c1=0., c2=0., d1=0., d2=0., e1=0., e2=0., f1=0., f2=0.,
                        k1=None, k2=None,
                        nsm=0., rc=0., zc=0., delta_n=0., comment=''):
        """
        +-------+------+-------+-----+----+----+--------+----+--------+
        |   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
        +=======+======+=======+=====+====+====+========+====+========+
        | PBEND | PID  |  MID  | A   | I1 | I2 |   J    | RB | THETAB |
        +-------+------+-------+-----+----+----+--------+----+--------+
        |       |  C1  |  C2   | D1  | D2 | E1 |   E2   | F1 |   F2   |
        +-------+------+-------+-----+----+----+--------+----+--------+
        |       |  K1  |  K2   | NSM | RC | ZC | DELTAN |    |        |
        +-------+------+-------+-----+----+----+--------+----+--------+

        Parameters
        ----------
        A : float
            cross-sectional area
        i1, i2 : float
            area moments of inertia for plane 1/2
        j : float
            torsional stiffness
        rb : float; default=None
            bend radius of the line of centroids
        theta_b : float; default=None
            arc angle of element (degrees)
        c1, c2, d1, d2, e1, e2, f1, f2 : float; default=0.0
            the r/z locations from the geometric centroid for stress recovery
        k1, k2 : float; default=None
            Shear stiffness factor K in K*A*G for plane 1 and plane 2
        nsm : float; default=0.
            nonstructural mass per unit length???
        zc : float; default=None
            Offset of the geometric centroid in a direction perpendicular to
            the plane of points GA and GB and vector v.
        delta_n : float; default=None
            Radial offset of the neutral axis from the geometric centroid,
            positive is toward the center of curvature
        """
        beam_type = 1
        fsi = None
        rm = None
        t = None
        p = None
        return PBEND(pid, mid, beam_type, A, i1, i2, j,
                     c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                     nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b, comment=comment)

    @classmethod
    def add_beam_type_2(cls, pid, mid,
                        fsi, rm, t, p=None, rb=None, theta_b=None,
                        nsm=0., rc=0., zc=0., comment=''):
        """
        +-------+------+-------+-----+----+----+--------+----+--------+
        |   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
        +=======+======+=======+=====+====+====+========+====+========+
        | PBEND | PID  |  MID  | FSI | RM | T  |   P    | RB | THETAB |
        +-------+------+-------+-----+----+----+--------+----+--------+
        |       |      |       | NSM | RC | ZC |        |    |        |
        +-------+------+-------+-----+----+----+--------+----+--------+

        Parameters
        ----------
        fsi : int
            Flag selecting the flexibility and stress intensification
            factors. See Remark 3. (Integer = 1, 2, or 3)
        rm : float
            Mean cross-sectional radius of the curved pipe
        t : float
            Wall thickness of the curved pipe
        p : float; default=None
            Internal pressure
        rb : float; default=None
            bend radius of the line of centroids
        theta_b : float; default=None
            arc angle of element (degrees)
        nsm : float; default=0.
            nonstructural mass per unit length???
        rc : float; default=None
            Radial offset of the geometric centroid from points GA and GB.
        zc : float; default=None
            Offset of the geometric centroid in a direction perpendicular
            to the plane of points GA and GB and vector v
        """
        beam_type = 2
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
        return PBEND(pid, mid, beam_type, A, i1, i2, j,
                     c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                     nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b, comment=comment)

    #@classmethod
    #def add_beam_type_3(cls, pid, mid,
                        #fsi, rm, t, p, rb, theta_b,
                        ##sacl, alpha, nsm, rc, zc, flange, kx, ky, kz, sy, sz,
                        #comment=''):
        #"""
        #+-------+------+-------+-----+----+----+--------+----+--------+
        #|   1   |   2  |   3   |  4  |  5 |  6 |   7    |  7 |    8   |
        #+=======+======+=======+=====+====+====+========+====+========+
        #| PBEND | PID  |  MID  | FSI | RM | T  |   P    | RB | THETAB |
        #+-------+------+-------+-----+----+----+--------+----+--------+
        #|       | SACL | ALPHA | NSM | RC | ZC | FLANGE |    |        |
        #+-------+------+-------+-----+----+----+--------+----+--------+
        #|       |  KX  |  KY   | KZ  |    | SY |   SZ   |    |        |
        #+-------+------+-------+-----+----+----+--------+----+--------+
        #"""
        #beam_type = 3
        #A = None
        #i1 = None
        #i2 = None
        #j = None
        #c1 = None
        #c2 = None
        #d1 = None
        #d2 = None
        #e1 = None
        #e2 = None
        #f1 = None
        #f2 = None
        #k1 = None
        #k2 = None
        #rc = None
        #zc = None
        #nsm = None
        #delta_n = None
        #return PBEND(pid, mid, beam_type, A, i1, i2, j,
                     #c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                     #nsm, rc, zc, delta_n, fsi, rm, t, p, rb, theta_b, comment=comment)

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

        elif isinstance(value3, int):  # alternate form
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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by PBEND mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def MassPerLength(self) -> float:
        """m/L = rho*A + nsm"""
        rho = self.mid_ref.Rho()
        assert isinstance(self.A, float), self.get_stats()
        return self.A * rho + self.nsm

    def raw_fields(self) -> List[Union[str, float, int, None]]:
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


def split_arbitrary_thickness_section(key, value):
    # type: (str, Union[str, float, List[int]]) -> Tuple[int, Union[float, List[int]]]
    """
    Helper method for PBRSECT/PBMSECT

    >>> key = 'T(11)'
    >>> value = '[1.2,PT=(123,204)]'
    >>> index, out = split_arbitrary_thickness_section(key, value)
    >>> index
    11
    >>> out
    [1.2, [123, 204]]
    """
    assert key.endswith(')'), 'key=%r' % key
    # T(3), CORE(3)
    key_id = key[:-1].split('(', 1)[1]
    key_id = int(key_id)

    if isinstance(value, (int, float)):
        return key_id, value

    value = value.replace(' ', '')
    if 'PT' in value:
        bracketed_values = value.strip('[]')
        sline = bracketed_values.split(',', 1)
        thicknessi = float(sline[0])
        pt_value = sline[1].split('=')
        assert pt_value[0] == 'PT', pt_value
        points = pt_value[1].strip('()').split(',')
        assert len(points) == 2, pt_value
        int_points = [int(pointi) for pointi in points]
        out = [thicknessi, int_points]
    else:
        out = float(value)
    return key_id, out


def get_beam_sections(line: str) -> List[str]:
    """
    Splits a PBRSECT/PBMSECT line

    >>> line = 'OUTP=10,BRP=20,T=1.0,T(11)=[1.2,PT=(123,204)], NSM=0.01'
    >>> sections = get_beam_sections(line)
    >>> sections
    ['OUTP=10', 'BRP=20', 'T=1.0', 'T(11)=[1.2,PT=(123,204)', 'NSM=0.01'], sections
    """
    line = line.replace(' ', '')
    words = []
    i0 = None
    nopen_parantheses = 0
    nopen_brackets = 0
    i = 0
    while i < len(line):
        char = line[i]
        if char == '(':
            nopen_parantheses += 1
        elif char == ')':
            nopen_parantheses -= 1
        elif char == '[':
            nopen_brackets += 1
        elif char == ']':
            nopen_brackets -= 1

        elif nopen_parantheses == 0 and nopen_brackets == 0 and char == ',':
            word = line[i0:i].strip(',')
            words.append(word)
            i0 = i
        i += 1
    word = line[i0:].strip(',')
    if word:
        words.append(word)
    return words

def write_arbitrary_beam_section(inps, ts, branch_paths, nsm, outp_id, core=None):
    """writes the PBRSECT/PBMSECT card"""
    end = ''
    for key, dicts in [('INP', inps), ('T', ts), ('BRP', branch_paths), ('CORE', core)]:
        if dicts is None:
            continue
        # dicts = {int index : int/float value}
        for index, value1 in sorted(dicts.items()):
            if index == 0:
                if isinstance(value1, list):
                    for value1i in value1:
                        end += '        %s=%s,\n' % (key, value1i)
                else:
                    end += '        %s=%s,\n' % (key, value1)
            else:
                if isinstance(value1, list):
                    if len(value1) == 2:
                        assert len(value1) == 2, value1
                        thicknessi = value1[0]
                        points = value1[1]
                        end += '        %s(%s)=[%s, PT=(%s,%s)],\n' % (
                            key, index, thicknessi, points[0], points[1])
                    else:
                        assert len(value1) == 1, value1
                        end += '        %s=%s,\n' % (key, value1[0])
                else:
                    end += '        %s(%s)=%s,\n' % (key, index, value1)

    for key, value in [('NSM', nsm), ('outp', outp_id),]:
        if value:
            end += '        %s=%s,\n' % (key, value)
    if end:
        end = end[:-2] + '\n'
    return end

def plot_arbitrary_section(model, self,
                           inps, ts, branch_paths, nsm, outp_ref,
                           figure_id=1, title='', show=False):
    """helper for PBRSECT/PBMSECT"""
    import matplotlib.pyplot as plt
    if ts:
        try:
            ts2 = {1 : ts[1]}
        except KeyError:
            print('ts =', ts)
            print(self)
            raise
        for key, value in ts.items():
            if key == 1:
                ts2[key] = value
            else:
                thickness_value, section = value
                p1, p2 = section
                p1, p2 = min(p1, p2), max(p1, p2)
                ts2[(p1, p2)] = thickness_value
        ts = ts2

    def _plot_rectangles(ax, sections, xy_dict, ts):
        """helper method for ``plot_arbitrary_section``"""
        for section in sections:
            #p1, p2 = section
            out = xy_dict[section]
            (x1, x2, y1, y2) = out
            dy = y2 - y1
            dx = x2 - x1
            length = np.sqrt(dy**2 + dx**2)
            angle = np.arctan2(dy, dx)
            angled = np.degrees(angle)

            thickness = ts.get(section, ts[1])
            assert isinstance(thickness, float), thickness

            width = thickness

            dx_width = +width / 2. * np.sin(angle)
            dy_width = -width / 2. * np.cos(angle)
            xy = (x1+dx_width, y1+dy_width)

            rect_height = length
            rect_width = width

            #print('dxy_width = (%.2f,%.2f)' % (dx_width, dy_width))
            #print('p1,2=(%s, %s) xy=(%.2f,%.2f) t=%s width=%s height=%s angled=%s\n' % (
                #p1, p2, xy[0], xy[1], thickness, rect_width, rect_height, angled))

            rect = plt.Rectangle(xy, rect_height, rect_width, angle=angled,
                                 fill=True) #, alpha=1.2+0.15*i)
            ax.add_patch(rect)
            #break

    def add_to_sections(sections, xy, points, x, y):
        """helper method for ``plot_arbitrary_section``"""
        for i in range(len(points)-1):
            p1 = points[i]
            p2 = points[i+1]
            x1 = x[i]
            x2 = x[i+1]
            y1 = y[i]
            y2 = y[i+1]
            p1, p2 = min(p1, p2), max(p1, p2)
            sections.add((p1, p2))
            xy[(p1, p2)] = (x1, x2, y1, y2)

    fig = plt.figure(figure_id)
    ax = fig.add_subplot(111, aspect='equal')
    #print('outp:\n%s' % outp_ref)
    out_points = outp_ref.ids
    if self.form == 'CP' and out_points[0] != out_points[-1]:
        out_points = out_points + [out_points[0]]

    #out_points_ref = outp_ref.ids_ref
    #print('out_points =', out_points)
    #out_xyz = np.array([point.get_position() for point in out_points_ref])
    out_xyz = np.array([model.points[point_id].get_position()
                        for point_id in out_points])
    #print('out_xyz:\n%s' % out_xyz)
    sections = set()
    x = out_xyz[:, 0]
    y = out_xyz[:, 1]
    xy = {}
    add_to_sections(sections, xy, out_points, x, y)
    #print('x=%s y=%s' % (x, y))
    ax.plot(x, y, '-o', label='OUTP')
    all_points = {point_id : (xi, yi)
                  for point_id, xi, yi in zip(out_points, x, y)}
    #print('out_points =', out_points)
    #print('all_points =', all_points)
    #plt.show()

    for key, unused_thickness in ts.items():
        if key == 1:
            continue
        #print(thickness)
        section = key
        #print(model.points)
        #thickness_value, section = thickness
        out_xyz = np.array([model.points[point_id].get_position()
                            for point_id in section])
        x = out_xyz[:, 0]
        y = out_xyz[:, 1]
        #print('adding t=%s section %s' % (thickness, str(section)))
        #print(sections)

        add_to_sections(sections, xy, section, x, y)
        #print(sections)

    if branch_paths:
        for key, brp_set_ref in branch_paths.items():
            brp_points = brp_set_ref.ids
            brp_points_ref = brp_set_ref.ids_ref
            brp_xyz = np.array([point.get_position() for point in brp_points_ref])
            #print('branch = %s' % brp_points)
            x = brp_xyz[:, 0]
            y = brp_xyz[:, 1]
            ax.plot(x, y, '-o', label='BRP(%i)' % key)
            add_to_sections(sections, xy, brp_points, x, y)

            for point_id, xi, yi in zip(brp_points, x, y):
                all_points[point_id] = (xi, yi)

    #print('xy =', xy)

    if inps:
        #print('inps! = %s' % inps)
        for key, inp in inps.items():
            if isinstance(inp, int):
                inp = [inp]
            for inpi in inp:
                inp_ref = model.Set(inpi)
                #inp_ref.cross_reference_set(model, 'Point', msg='')
                inp_points = inp_ref.ids
                #if inp_points[0] != inp_points[-1]:
                    #inp_points = inp_points + [inp_points[0]]
                #print('inp_points = %s' % inp_points)
                #inp_points_ref = inp_ref.ids_ref
                inp_xyz = np.array([model.points[point_id].get_position()
                                    for point_id in inp_points])
                #inp_xyz = np.array([point.get_position() for point in inp_points_ref])
                #print('inp_xyz:\n%s' % (inp_xyz[:, :2]))
                x = inp_xyz[:, 0]
                y = inp_xyz[:, 1]
                ax.plot(x, y, '--x', label='INP')


    if ts:
        _plot_rectangles(ax, sections, xy, ts)

    #print('all_points =', all_points)
    for point_id, xy in sorted(all_points.items()):
        ax.annotate(str(point_id), xy=xy)

    ax.grid(True)
    ax.set_title(title)
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    ax.legend()
    if show:
        plt.show()
