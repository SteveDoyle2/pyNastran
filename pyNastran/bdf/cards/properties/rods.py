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
from __future__ import annotations
from typing import TYPE_CHECKING
from numpy import pi

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

class PROD(Property):
    """
    +------+-----+-----+-----+-----+-----+-----+
    |  1   |  2  |  3  |  4  |  5  |  6  |  7  |
    +======+=====+=====+=====+=====+=====+=====+
    | PROD | PID | MID |  A  |  J  |  C  | NSM |
    +------+-----+-----+-----+-----+-----+-----+
    | PROD |  1  |  2  | 2.0 | 3.0 | 0.5 | 1.0 |
    +------+-----+-----+-----+-----+-----+-----+
    """
    type = 'PROD'
    pname_fid_map = {
        4 : 'A', 'A' : 'A',
        5 : 'J', 'J' : 'j',
        6 : 'C', 'C' : 'c',
    }

    @classmethod
    def export_to_hdf5(cls, h5_file, model, pids):
        """exports the properties in a vectorized way"""
        #comments = []
        mids = []
        A = []
        J = []
        c = []
        nsm = []
        for pid in pids:
            prop = model.properties[pid]
            #comments.append(prop.comment)
            mids.append(prop.mid)
            A.append(prop.A)
            J.append(prop.j)
            c.append(prop.c)
            nsm.append(prop.nsm)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('mid', data=mids)
        h5_file.create_dataset('A', data=A)
        h5_file.create_dataset('J', data=J)
        h5_file.create_dataset('c', data=c)
        h5_file.create_dataset('nsm', data=nsm)
        #h5_file.create_dataset('_comment', data=comments)

    def __init__(self, pid, mid, A, j=0., c=0., nsm=0., comment=''):
        """
        Creates a PROD card

        Parameters
        ----------
        pid : int
           property id
        mid : int
           material id
        A : float
           area
        J : float; default=0.
           polar moment of inertia
        c : float; default=0.
           stress factor
        nsm : float; default=0.
           nonstructural mass per unit length
        comment : str; default=''
            a comment for the card
        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.mid = mid
        self.A = A
        self.j = j
        self.c = c
        self.nsm = nsm
        self.mid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PROD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        A = double(card, 3, 'A')
        j = double_or_blank(card, 4, 'J', 0.0)
        c = double_or_blank(card, 5, 'c', 0.0)
        nsm = double_or_blank(card, 6, 'nsm', 0.0)
        assert len(card) <= 7, 'len(PROD card) = %i\ncard=%s' % (len(card), card)
        return PROD(pid, mid, A, j, c, nsm, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PROD card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        pid = data[0]
        mid = data[1]
        A = data[2]
        j = data[3]
        c = data[4]
        nsm = data[5]
        return PROD(pid, mid, A, j, c, nsm, comment=comment)

    def _verify(self, xref):
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
        if self.mid_ref is None:
            raise RuntimeError('Material mid=%i has not been cross referenced.\n%s' % (self.mid, str(self)))
        return self.mid_ref.E()

    def G(self):
        if self.mid_ref is None:
            raise RuntimeError('Material mid=%i has not been cross referenced.\n%s' % (self.mid, str(self)))
        return self.mid_ref.G()

    def Rho(self):
        if self.mid_ref is None:
            raise RuntimeError('Material mid=%i has not been cross referenced.\n%s' % (self.mid, str(self)))
        return self.mid_ref.Rho()

    def MassPerLength(self):
        r"""
            Gets the mass per length :math:`\frac{m}{L}` of the CBEAM.

            .. math:: \frac{m}{L} = A \rho + nsm
            """
        rho = self.Rho()
        area = self.A
        nsm = self.nsm
        return area * rho + nsm

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PROD mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PROD mid=%s' % self.mid
        self.mid_ref = model.safe_material(self.mid, self.pid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)


class PTUBE(Property):
    """
    +-------+------+-----+------+------+------+-----+
    |   1   |  2   |  3  |   4  |  5   |  6   |  7  |
    +=======+======+=====+======+======+======+=====+
    | PTUBE | PID  | MID |  OD  |   T  |  NSM | OD2 |
    +-------+------+-----+------+------+------+-----+
    | PTUBE |  2   |  6  | 6.29 | 0.25 |      |     |
    +-------+------+-----+------+------+------+-----+
    """
    type = 'PTUBE'
    pname_fid_map = {
        4 : 'OD1', 'OD' : 'OD1',
        5 : 't', 'T' : 't',
        7 : 'OD2', 'OD2' : 'OD2',
    }

    @classmethod
    def export_to_hdf5(cls, h5_file, model, pids):
        """exports the properties in a vectorized way"""
        #comments = []
        mids = []
        OD = []
        t = []
        nsm = []
        for pid in pids:
            prop = model.properties[pid]
            #comments.append(prop.comment)
            mids.append(prop.mid)
            OD.append([prop.OD1, prop.OD2])
            t.append(prop.t)
            nsm.append(prop.nsm)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('pid', data=pids)
        h5_file.create_dataset('mid', data=mids)
        h5_file.create_dataset('OD', data=OD)
        h5_file.create_dataset('t', data=t)
        h5_file.create_dataset('nsm', data=nsm)
        #h5_file.create_dataset('_comment', data=comments)

    def __init__(self, pid, mid, OD1, t=None, nsm=0., OD2=None, comment=''):
        """
        Adds a PTUBE card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        OD1 : float
            outer diameter at End A
        t : float; default=None -> OD1/2.
            thickness
        nsm : float; default=0.
            non-structural mass per unit length
        OD2 : float; default=None -> OD1
            outer diameter at End B
        comment : str; default=''
            a comment for the card
        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        if t is None:
            t = OD1 / 2.
        if OD2 is None:
            OD2 = OD1
        self.pid = pid
        self.mid = mid
        self.OD1 = OD1
        self.OD2 = OD2
        self.t = t
        self.nsm = nsm
        self.mid_ref = None

    def validate(self):
        assert self.OD1 >= 0., 'PTUBE: pid=%s OD1=%s' % (self.pid, self.OD1)
        assert self.OD2 >= 0., 'PTUBE: pid=%s OD2=%s' % (self.pid, self.OD2)
        assert self.t >= 0., 'PTUBE: pid=%s t=%s' % (self.pid, self.t)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PTUBE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        OD1 = double(card, 3, 'OD1')
        t = double_or_blank(card, 4, 't', OD1/2.)
        nsm = double_or_blank(card, 5, 'nsm', 0.0)
        OD2 = double_or_blank(card, 6, 'OD2', OD1)
        assert len(card) <= 7, 'len(PTUBE card) = %i\ncard=%s' % (len(card), card)
        return PTUBE(pid, mid, OD1, t, nsm, OD2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PTUBE card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        pid = data[0]
        mid = data[1]
        OD1 = data[2]
        t = data[3]
        nsm = data[4]
        OD2 = OD1
        #OD2 = data[5]  #: .. note:: quirk to this one...
        return PTUBE(pid, mid, OD1, t, nsm, OD2, comment=comment)

    def _verify(self, xref):
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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PTUBE mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by PTUBE mid=%s' % self.mid
        self.mid_ref = model.safe_material(self.mid, self.pid, xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def Rho(self):
        """
        Gets the density :math:`\rho` of the CTUBE.
        """
        if self.mid_ref is None:
            raise RuntimeError('Material mid=%i has not been cross referenced.\n%s' % (self.mid, str(self)))
        return self.mid_ref.Rho()

    def MassPerLength(self):
        r"""
        Gets the mass per length :math:`\frac{m}{L}` of the CTUBE.

        .. math:: \frac{m}{L} = (A \rho) nsm
        """
        return self.Area() * self.Rho() + self.nsm

    def E(self):
        if self.mid_ref is None:
            raise RuntimeError('Material mid=%i has not been cross referenced.\n%s' % (self.mid, str(self)))
        return self.mid_ref.E()

    def G(self):
        if self.mid_ref is None:
            raise RuntimeError('Material mid=%i has not been cross referenced.\n%s' % (self.mid, str(self)))
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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
