# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All beam properties are defined in this file.  This includes:
 *   PBEAM
 *   PBEAML
 *   PBAR
 *   PBARL

All beams are Property objects.
Multi-segment beams are IntegratedLineProperty objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from six.moves import zip
from itertools import count
from numpy import pi, array

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Property
from pyNastran.bdf.bdfInterface.assign_type import (integer,
    double, double_or_blank)
from pyNastran.utils.mathematics import integrate_line, integrate_positive_line
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

class PROD(Property):
    type = 'PROD'

    def __init__(self):
        Property.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.pid = integer(card, 1, 'pid')
        self.mid = integer(card, 2, 'mid')
        self.A = double(card, 3, 'A')
        self.j = double_or_blank(card, 4, 'J', 0.0)
        self.c = double_or_blank(card, 5, 'c', 0.0)
        self.nsm = double_or_blank(card, 6, 'nsm', 0.0)
        assert len(card) <= 7, 'len(PROD card) = %i' % len(card)

    def add_op2_data(self, data):
        self.pid = data[0]
        self.mid = data[1]
        self.A = data[2]
        self.j = data[3]
        self.c = data[4]
        self.nsm = data[5]

    def _verify(self, xref=False):
        pid = self.Pid()
        mid = self.Mid()
        A = self.Area()
        J = self.J()
        c = self.c
        nsm = self.Nsm()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(A, float), 'pid=%r' % A
        assert isinstance(J, float), 'cid=%r' % J
        assert isinstance(c, float), 'c=%r' % c
        assert isinstance(nsm, float), 'nsm=%r' % nsm

    def Area(self):
        return self.A

    def J(self):
        return self.j

    def C(self):
        return self.c

    def Nsm(self):
        return self.nsm

    def E(self):
        return self.mid_ref.E()

    def G(self):
        return self.mid_ref.G()

    def Rho(self):
        return self.mid_ref.Rho()

    def cross_reference(self, model):
        msg = ' which is required by PROD mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

    def write_code_aster(self, icut, iface, istart):  # PROD
        msg = ''
        msg += "    POUTRE=_F(GROUP_MA='P%s', # PROD\n" % (self.pid)
        msg += "              SECTION='CERCLE',  # circular section\n"
        msg += "              CARA=('R')   # radius\n"
        msg += "              VALE=(%g),),\n" % (self.Radius())

        msg += "              SECTION='GENERALE',\n"
        msg += "              CARA=('A', 'JX')\n"
        msg += "              VALE=(%g, %g),\n"  %(self.Area(), self.J())
        msg += "                    CARA='VECT_Y'),\n"
        msg += "                    VALE=(1.0,0.0,0.0,),),\n"
        return (msg, icut, iface, istart)

    def raw_fields(self):
        list_fields = ['PROD', self.pid, self.Mid(), self.A, self.j, self.c,
                       self.nsm]
        return list_fields

    def repr_fields(self):
        j = set_blank_if_default(self.j, 0.0)
        c = set_blank_if_default(self.c, 0.0)
        nsm = set_blank_if_default(self.nsm, 0.0)
        list_fields = ['PROD', self.pid, self.Mid(), self.A, j, c, nsm]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PTUBE(Property):
    type = 'PTUBE'

    def __init__(self):
        Property.__init__(self)

    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        self.pid = integer(card, 1, 'pid')
        self.mid = integer(card, 2, 'mid')
        self.OD1 = double(card, 3, 'OD1')
        self.t = double_or_blank(card, 4, 't')
        if self.t is None:
            self.t = self.OD1 / 2.
        self.nsm = double_or_blank(card, 5, 'nsm', 0.0)
        self.OD2 = double_or_blank(card, 6, 'OD2', self.OD1)
        assert len(card) <= 7, 'len(PTUBE card) = %i' % len(card)

    def add_op2_data(self, data):
        self.pid = data[0]
        self.mid = data[1]
        self.OD1 = data[2]
        self.t = data[3]
        self.nsm = data[4]
        self.OD2 = self.OD1
        #self.OD2 = data[5]  #: .. note:: quirk to this one...

    def _verify(self, xref=False):
        pid = self.Pid()
        mid = self.Mid()
        A = self.Area()
        nsm = self.Nsm()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(A, float), 'A=%r' % A
        assert isinstance(nsm, float), 'nsm=%r' % nsm

    def Nsm(self):
        """
        Gets the non-structural mass :math:`nsm` of the CTUBE.
        """
        return self.nsm

    def cross_reference(self, model):
        msg = ' which is required by PTUBE mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

    def Rho(self):
        """
        Gets the density :math:`\rho` of the CTUBE.
        """
        return self.mid_ref.Rho()

    def MassPerLength(self):
        r"""
        Gets the mass per length :math:`\frac{m}{L}` of the CTUBE.

        .. math:: \frac{m}{L} = (A \rho) nsm
        """
        return self.Area() * self.Rho() + self.nsm

    def E(self):
        return self.mid_ref.E()

    def G(self):
        return self.mid_ref.G()

    def J(self):
        Dout = self.OD1
        if self.t == 0.0:
            return pi / 8. * Dout**4
        Din = Dout - 2 * self.t
        return pi / 8. * (Dout**4 - Din**2)

    def Area(self):
        r"""
        Gets the area :math:`A` of the CTUBE.

        .. math:: A_1 = \pi \frac{d_1^2}{4} - \pi {(D_1-2t)^2}{4}

        .. math:: A_2 = \pi \frac{d_2^2}{4} - \pi {(D_2-2t)^2}{4}

        .. math:: A = A_1 + A_2
        """
        A = (self._area1() + self._area2()) / 2.

        #A1 = pi*D1^2/4 - pi*((D1-2t)^2)/4
        #A2 = pi*D2^2/4 - pi*((D2-2t)^2)/4
        #A = A1 + A2

        #A = pi*D1*t/2 + pi*D2*t/2 - pi*t
        #A = pi*t*(D1/2 + D2/2 - t)
        #A = pi*t*( (D1+D2)/2.-t )

        #A2 = pi*t*( (D1+D2)/2.-t )

        #if A != A2:
            #msg = ('AREA method has problem in PTUBE '
            #       'Aold=%s != Anew=%s' %(A, A2))
            #raise RuntimeError(msg)
        return A

    def _area1(self):
        """Gets the Area of Section 1 of the CTUBE."""
        Dout = self.OD1
        if self.t == 0:
            return pi / 4. * Dout**2
        Din = Dout - 2 * self.t
        A1 = pi / 4. * (Dout * Dout - Din * Din)
        return A1

    def _area2(self):
        """Gets the Area of Section 2 of the CTUBE."""
        Dout = self.OD2
        if self.t == 0:
            return pi / 4. * Dout**2
        Din = Dout - 2 * self.t
        A2 = pi / 4. * (Dout * Dout - Din * Din)
        return A2

    #def massMatrix(self):
        #"""
        #.. todo:: not done
        #"""
        #m = zeros(6, 6)
        #m[0, 0] = 1.
        #return m

    def raw_fields(self):
        list_fields = ['PTUBE', self.pid, self.Mid(), self.OD1, self.t,
                       self.nsm, self.OD2]
        return list_fields

    def repr_fields(self):
        t = set_blank_if_default(self.t, self.OD1 / 2.)
        nsm = set_blank_if_default(self.nsm, 0.0)
        OD2 = set_blank_if_default(self.OD2, self.OD1)
        list_fields = ['PTUBE', self.pid, self.Mid(), self.OD1, t, nsm, OD2]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
