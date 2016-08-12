# pylint: disable=R0904,R0902,E1101,E1103,C0111,C0302,C0103,W0101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types

import numpy as np
from numpy.linalg import norm

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Element
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_double_or_blank, double_or_blank,
    integer_string_or_blank, string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class LineElement(Element):  # CBAR, CBEAM, CBEAM3, CBEND
    def __init__(self):
        Element.__init__(self)

    def C(self):
        """torsional constant"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.C()

    def Area(self):
        """returns the area of the element face"""
        raise NotImplementedError('implement self.Area() for %s' % self.type)

    def E(self):
        """returns the Young's Modulus, :math:`E`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.mid_ref.E()

    def G(self):
        """returns the Shear Modulus, :math:`G`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.mid_ref.G()

    def J(self):
        """returns the Polar Moment of Inertia, :math:`J`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.J()

    def I11(self):
        """returns the Moment of Inertia, :math:`I_{11}`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.I11()

    def I22(self):
        """returns the Moment of Inertia, :math:`I_{22}`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.I22()

    def I12(self):
        """returns the Moment of Inertia, :math:`I_{12}`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.I12()

    def Nu(self):
        """Get Poisson's Ratio, :math:`\nu`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.mid_ref.nu

    def Rho(self):
        """Get the material density, :math:`\rho`"""
        #print(str(self.pid), type(self.pid))
        #raise NotImplementedError('implement self.Rho() for %s' % self.type)
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.mid_ref.rho

    def Nsm(self):
        """Placeholder method for the non-structural mass, :math:`nsm`"""
        raise NotImplementedError('implement self.Area() for %s' % self.type)

    def MassPerLength(self):
        """
        Get the mass per unit length, :math:`\frac{m}{L}`
        """
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.MassPerLength()

    def Mass(self):
        r"""
        Get the mass of the element.

        .. math:: m = \left( \rho A + nsm \right) L
        """
        L = self.Length()
        mass = L * self.MassPerLength()
        #try:
            #mass = (self.Rho() * self.Area() + self.Nsm()) * L
        #except TypeError:
            #msg = 'TypeError on eid=%s pid=%s:\n' % (self.eid, self.Pid())
            #msg += 'rho = %s\narea = %s\nnsm = %s\nL = %s' % (self.Rho(),
            #                                                  self.Area(),
            #                                                  self.Nsm(), L)
            #raise TypeError(msg)

        return mass

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        #self.g0 = model.nodes[self.g0]
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()
        del self.nodes_ref, self.pid_Ref

    def Length(self):
        r"""
        Gets the length, :math:`L`, of the element.

        .. math:: L = \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  }
        """
        L = norm(self.nodes_ref[1].get_position() - self.nodes_ref[0].get_position())
        return L

    def get_edge_ids(self):
        """
        Return the edge IDs
        """
        node_ids = self.node_ids
        return [(node_ids[0], node_ids[1])]


class CBAROR(object):
    type = 'CBAROR'
    def __init__(self):
        self.n = 0

    def add(self, card=None, data=None, comment=''):
        if self.n == 1:
            raise RuntimeError('only one CBAROR is allowed')
        self.n = 1
        if comment:
            self._comment = comment

        self.property_id = integer_or_blank(card, 2, 'pid')

        # x / g0
        field5 = integer_double_or_blank(card, 5, 'g0_x1', 0.0)
        if isinstance(field5, integer_types):
            self.is_g0 = True
            self.g0 = field5
            self.x = [0., 0., 0.]
        elif isinstance(field5, float):
            self.is_g0 = False
            self.g0 = None
            self.x = np.array([field5,
                               double_or_blank(card, 6, 'x2', 0.0),
                               double_or_blank(card, 7, 'x3', 0.0)],
                              dtype='float64')
        self.offt = string_or_blank(card, 8, 'offt', 'GGG')
        assert len(card) <= 9, 'len(CBAROR card) = %i\ncard=%s' % (len(card), card)


class CBAR(LineElement):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |   9  |
    +=======+=====+=====+=====+=====+=====+=====+=====+======+
    | CBAR  | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B  |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |   9  |
    +=======+=====+=====+=====+=====+=====+=====+=====+======+
    | CBAR  | EID | PID | GA  | GB  | G0  |     |     | OFFT |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B  |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+

    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    |   1   |   2   |  3  |   4   |   5   |    6   |   7   |   8   |   9   |
    +=======+=======+=====+=======+=======+========+=======+=======+=======+
    |  CBAR | 2     |  39 | 7     | 6     |  105   |       |       |  GGG  |
    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    |       |       | 513 | 0.0+0 | 0.0+0 |    -9. | 0.0+0 | 0.0+0 |   -9. |
    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    """
    type = 'CBAR'
    aster_type = 'CBAR'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb',
        8:'offt', 9:'pa', 10:'pb',
    }

    def _update_field_helper(self, n, value):
        if n == 11:
            self.wa[0] = value
        elif n == 12:
            self.wa[1] = value
        elif n == 13:
            self.wa[2] = value
        elif n == 14:
            self.wb[0] = value
        elif n == 15:
            self.wb[1] = value
        elif n == 16:
            self.wb[2] = value
        else:
            if self.g0 is not None:
                if n == 5:
                    self.g0 = value
                else:
                    raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))
            else:
                if n == 5:
                    self.x[0] = value
                elif n == 6:
                    self.x[1] = value
                elif n == 7:
                    self.x[2] = value
                else:
                    raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, pid, ga, gb,
                 x, g0, offt='GGG',
                 pa=0, pb=0, wa=None, wb=None, comment=''):
        LineElement.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.x = x
        self.g0 = g0
        self.ga = ga
        self.gb = gb
        self.offt = offt
        self.pa = pa
        self.pb = pb
        self.wa = wa
        self.wb = wb
        self._validate_input()

    def _validate_input(self):
        if isinstance(self.offt, integer_types):
            assert self.offt in [1, 2], 'invalid offt; offt=%i' % self.offt
            raise NotImplementedError('invalid offt; offt=%i' % self.offt)
        elif not isinstance(self.offt, string_types):
            raise SyntaxError('invalid offt expected a string of length 3 '
                              'offt=%r; Type=%s' % (self.offt, type(self.offt)))

        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)

        msg = 'invalid offt parameter of %s...offt=%s' % (self.type, self.offt)
        # B,G,O
        assert self.offt[0] in ['G', 'B'], msg
        assert self.offt[1] in ['G', 'O', 'E'], msg
        assert self.offt[2] in ['G', 'O', 'E'], msg

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        x, g0 = cls._init_x_g0(card, eid)

        # doesn't exist in NX nastran
        offt = integer_string_or_blank(card, 8, 'offt', 'GGG')
        #print('cls.offt = %r' % (cls.offt))

        pa = integer_or_blank(card, 9, 'pa', 0)
        pb = integer_or_blank(card, 10, 'pb', 0)

        wa = np.array([double_or_blank(card, 11, 'w1a', 0.0),
                       double_or_blank(card, 12, 'w2a', 0.0),
                       double_or_blank(card, 13, 'w3a', 0.0)], dtype='float64')

        wb = np.array([double_or_blank(card, 14, 'w1b', 0.0),
                       double_or_blank(card, 15, 'w2b', 0.0),
                       double_or_blank(card, 16, 'w3b', 0.0)], dtype='float64')
        assert len(card) <= 17, 'len(CBAR card) = %i\ncard=%s' % (len(card), card)
        return CBAR(eid, pid, ga, gb, x, g0,
                    offt, pa, pb, wa, wb, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        #: .. todo:: verify
        #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,g0]]
        #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,x1,x2,x3]]

        main = data[0]
        flag = data[1][0]
        if flag in [0, 1]:
            g0 = None
            x = np.array([data[1][1],
                          data[1][2],
                          data[1][3]], dtype='float64')
        else:
            g0 = data[1][1]
            x = None

        eid = main[0]
        pid = main[1]
        ga = main[2]
        gb = main[3]
        #self.offt = str(data[4]) # GGG
        offt = 'GGG'  #: .. todo:: offt can be an integer; translate to char
        pa = main[4]
        pb = main[5]

        wa = np.array([main[6], main[7], main[8]], dtype='float64')
        wb = np.array([main[9], main[10], main[11]], dtype='float64')
        return CBAR(eid, pid, ga, gb, x, g0,
                    offt, pa, pb, wa, wb, comment=comment)

    def _verify(self, xref=False):
        eid = self.eid
        pid = self.Pid()
        edges = self.get_edge_ids()
        mid = self.Mid()
        nsm = self.Nsm()
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        if xref:  # True
            A = self.Area()
            mpl = self.MassPerLength()
            L = self.Length()
            mass = self.Mass()
            assert isinstance(A, float), 'eid=%s A=%r' % (eid, A)
            assert isinstance(L, float), 'eid=%s L=%r' % (eid, L)
            assert isinstance(mpl, float), 'eid=%s mass_per_length=%r' % (eid, mpl)
            assert isinstance(mass, float), 'eid=%s mass=%r' % (eid, mass)
            assert L > 0.0, 'eid=%s L=%s' % (eid, L)

    def Mid(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.Mid()

    def Area(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        A = self.pid_ref.Area()
        assert isinstance(A, float)
        return A

    def J(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        j = self.pid_ref.J()
        assert isinstance(j, float), 'J=%r for CBAR eid=%s pid=%s pidType=%s' % (j, self.eid, self.pid_ref.pid, self.pid_ref.type)
        return j

    def Length(self):
        # TODO: consider w1a and w1b in the length formulation
        L = norm(self.gb_ref.get_position() - self.ga_ref.get_position())
        assert isinstance(L, float)
        return L

    def Nsm(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        nsm = self.pid_ref.Nsm()
        assert isinstance(nsm, float)
        return nsm

    def I1(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.I1()

    def I2(self):
        return self.pid_ref.I2()

    def Centroid(self):
        return (self.ga_ref.get_position() + self.gb_ref.get_position()) / 2.

    @classmethod
    def _init_x_g0(cls, card, eid):
        field5 = integer_double_or_blank(card, 5, 'g0_x1', 0.0)
        if isinstance(field5, integer_types):
            g0 = field5
            x = None
        elif isinstance(field5, float):
            g0 = None
            x = np.array([field5,
                          double_or_blank(card, 6, 'x2', 0.0),
                          double_or_blank(card, 7, 'x3', 0.0)], dtype='float64')
            if norm(x) == 0.0:
                msg = 'G0 vector defining plane 1 is not defined.\n'
                msg += 'G0 = %s\n' % g0
                msg += 'X  = %s\n' % x
                raise RuntimeError(msg)
        else:
            msg = ('field5 on %s (G0/X1) is the wrong type...id=%s field5=%s '
                   'type=%s' % (cls.type, eid, field5, type(field5)))
            raise RuntimeError(msg)
        return x, g0

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        #if self.g0:
        #    self.x = nodes[self.g0].get_position() - nodes[self.ga].get_position()
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.ga_ref = self.ga
        self.gb = model.Node(self.gb, msg=msg)
        self.gb_ref = self.gb
        self.nodes = model.Nodes([self.ga.nid, self.gb.nid], msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid
        if model.is_nx:
            assert self.offt == 'GGG', 'NX only support offt=GGG; offt=%r' % self.offt

    def uncross_reference(self):
        self.pid = self.Pid()
        self.ga = self.Ga()
        self.gb = self.Gb()
        del self.pid_ref, self.ga_ref, self.gb_ref

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        else:
            return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, integer_types):
            return self.gb
        else:
            return self.gb_ref.nid

    def getX_G0_defaults(self):
        if self.g0 is not None:
            return (self.g0, None, None)
        else:
            #print('x =', self.x)
            #print('g0 =', self.g0)
            #x1 = set_blank_if_default(self.x[0], 0.0)
            #x2 = set_blank_if_default(self.x[1], 0.0)
            #x3 = set_blank_if_default(self.x[2], 0.0)
            return list(self.x)

    def get_orientation_vector(self, xyz):
        if self.g0:
            v = xyz[self.g0] - xyz[self.Ga()]
        else:
            v = self.x
        assert self.offt == 'GGG', self.offt
        return v

    def nodeIDs(self):
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    @property
    def nodes(self):
        return [self.ga, self.gb]

    @nodes.setter
    def nodes(self, values):
        self.ga = values[0]
        self.gb = values[1]

    def raw_fields(self):
        """
        .. todo:: not perfectly accurate b/c ???
        """
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = set_blank_if_default(self.offt, 'GGG')
        list_fields = ['CBAR', self.eid, self.Pid(), self.Ga(), self.Gb(), x1, x2,
                       x3, offt, self.pa, self.pb] + list(self.wa) + list(self.wb)
        return list_fields

    def repr_fields(self):
        pa = set_blank_if_default(self.pa, 0)
        pb = set_blank_if_default(self.pb, 0)

        w1a = set_blank_if_default(self.wa[0], 0.0)
        w2a = set_blank_if_default(self.wa[1], 0.0)
        w3a = set_blank_if_default(self.wa[2], 0.0)

        w1b = set_blank_if_default(self.wb[0], 0.0)
        w2b = set_blank_if_default(self.wb[1], 0.0)
        w3b = set_blank_if_default(self.wb[2], 0.0)
        (x1, x2, x3) = self.getX_G0_defaults()

        # offt doesn't exist in NX nastran
        offt = set_blank_if_default(self.offt, 'GGG')

        list_fields = ['CBAR', self.eid, self.Pid(), self.Ga(), self.Gb(), x1, x2,
                       x3, offt, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def write_card_16(self, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_16(card)


class CBEAM3(LineElement):  # was CBAR
    """
    Defines a three-node beam element
    """
    type = 'CBEAM3'

    def __init__(self, eid, pid, ga, gb, gc, x, g0,
                 wa, wb, wc, tw, s, comment=''):
        LineElement.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.ga = ga
        self.gb = gb
        self.gc = gc
        self.x = x
        self.g0 = g0
        self.wa = wa
        self.wb = wb
        self.wc = wc
        self.tw = tw
        self.s = s

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        gc = integer(card, 5, 'gc')

        x, g0 = cls._init_x_g0(card, eid)

        wa = np.array([double_or_blank(card, 9, 'w1a', 0.0),
                       double_or_blank(card, 10, 'w2a', 0.0),
                       double_or_blank(card, 11, 'w3a', 0.0)], dtype='float64')

        wb = np.array([double_or_blank(card, 12, 'w1b', 0.0),
                       double_or_blank(card, 13, 'w2b', 0.0),
                       double_or_blank(card, 14, 'w3b', 0.0)], dtype='float64')

        wc = np.array([double_or_blank(card, 15, 'w1c', 0.0),
                       double_or_blank(card, 16, 'w2c', 0.0),
                       double_or_blank(card, 17, 'w3c', 0.0)], dtype='float64')

        tw = np.array([double_or_blank(card, 18, 0., 'twa'),
                       double_or_blank(card, 19, 0., 'twb'),
                       double_or_blank(card, 20, 0., 'twc')], dtype='float64')

        s = np.array([integer_or_blank(card, 21, 'sa'),
                      integer_or_blank(card, 22, 'sb'),
                      integer_or_blank(card, 23, 'sc')], dtype='float64')
        assert len(card) <= 24, 'len(CBEAM3 card) = %i\ncard=%s' % (len(card), card)
        return CBEAM3(eid, pid, ga, gb, gc, x, g0,
                      wa, wb, wc, tw, s, comment='')

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.ga_ref = self.ga
        self.gb = model.Node(self.gb, msg=msg)
        self.gb_ref = self.gb
        self.gc = model.Node(self.gc, msg=msg)
        self.gc_ref = self.gc
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid

    def uncross_reference(self):
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.gc = self.Gc()
        self.pid = self.Pid()
        del self.ga_ref, self.gb_ref, self.gc_ref, self.pid_ref

    def Length(self):
        """
        # TODO: consider w1a and w1b in the length formulation
        # TODO: add gc to length formula
        """
        L = norm(self.gb_ref.get_position() - self.ga_ref.get_position())
        assert isinstance(L, float)
        return L

    def Area(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        A = self.pid_ref.Area()
        assert isinstance(A, float)
        return A

    def Nsm(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        nsm = self.pid_ref.Nsm()
        assert isinstance(nsm, float)
        return nsm

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        else:
            return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, integer_types):
            return self.gb
        else:
            return self.gb_ref.nid

    def Gc(self):
        if isinstance(self.gc, integer_types):
            return self.gc
        else:
            return self.gc_ref.nid

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb(), self.Gc()]

    def raw_fields(self):
        (x1, x2, x3) = self.getX_G0_defaults()
        (ga, gb, gc) = self.node_ids
        list_fields = ['CBEAM3', self.eid, self.Pid(), ga, gb, gc, x1, x2, x3] + \
                  list(self.wa) + list(self.wb) + list(self.wc) + list(self.tw) + list(self.s)
        return list_fields

    def repr_fields(self):
        w1a = set_blank_if_default(self.wa[0], 0.0)
        w2a = set_blank_if_default(self.wa[1], 0.0)
        w3a = set_blank_if_default(self.wa[2], 0.0)
        w1b = set_blank_if_default(self.wb[0], 0.0)
        w2b = set_blank_if_default(self.wb[1], 0.0)
        w3b = set_blank_if_default(self.wb[2], 0.0)
        w1c = set_blank_if_default(self.wc[0], 0.0)
        w2c = set_blank_if_default(self.wc[1], 0.0)
        w3c = set_blank_if_default(self.wc[2], 0.0)

        twa = set_blank_if_default(self.tw[0], 0.0)
        twb = set_blank_if_default(self.tw[1], 0.0)
        twc = set_blank_if_default(self.tw[2], 0.0)

        (x1, x2, x3) = self.getX_G0_defaults()
        (ga, gb, gc) = self.node_ids
        list_fields = ['CBEAM3', self.eid, self.Pid(), ga, gb, gc, x1, x2, x3,
                       w1a, w2a, w3a, w1b, w2b, w3b, w1c, w2c, w3c,
                       twa, twb, twc, self.s[0], self.s[1], self.s[2]]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def _verify(self, xref=False):
        edges = self.get_edge_ids()

class CBEND(LineElement):
    type = 'CBEND'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', 8:'geom',
    }

    def _update_field_helper(self, n, value):
        if self.g0 is not None:
            if n == 5:
                self.g0 = value
            else:
                raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))
        else:
            if n == 5:
                self.x[0] = value
            elif n == 6:
                self.x[1] = value
            elif n == 7:
                self.x[2] = value
            else:
                raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, eid, pid, ga, gb, g0, x, geom, comment=''):
        LineElement.__init__(self)
        if comment:
            self._comment = comment
        self.eid = eid
        self.pid = pid
        self.ga = ga
        self.gb = gb
        self.g0 = g0
        self.x = x
        self.geom = geom
        assert self.geom in [1, 2, 3, 4], 'geom is invalid geom=%r' % self.geom
        self._validate_input()

    @classmethod
    def add_card(cls, card, comment=''):
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')
        x1_g0 = integer_double_or_blank(card, 5, 'x1_g0', 0.0)
        if isinstance(x1_g0, integer_types):
            g0 = x1_g0
            x = None
        elif isinstance(x1_g0, float):
            g0 = None
            x = np.array([double_or_blank(card, 5, 'x1', 0.0),
                          double_or_blank(card, 6, 'x2', 0.0),
                          double_or_blank(card, 7, 'x3', 0.0)], dtype='float64')
            if norm(x) == 0.0:
                msg = 'G0 vector defining plane 1 is not defined.\n'
                msg += 'G0 = %s\n' % g0
                msg += 'X  = %s\n' % x
                raise RuntimeError(msg)
        else:
            raise ValueError('invalid x1Go=%r on CBEND' % x1_g0)
        geom = integer(card, 8, 'geom')

        assert len(card) == 9, 'len(CBEND card) = %i\ncard=%s' % (len(card), card)
        return CBEND(eid, pid, ga, gb, g0, x, geom, comment=comment)

    #def add_op2_data(self, data, comment=''):
        #if comment:
            #self._comment = comment
        #raise NotImplementedError(data)

    def Length(self):
        # TODO: consider w1a and w1b in the length formulation
        L = norm(self.gb_ref.get_position() - self.ga_ref.get_position())
        assert isinstance(L, float)
        return L

        #prop = self.pid_ref
        #bend_radius = prop.rb
        #theta_bend = prop.thetab
        #length_oa = None
        if self.geom == 1:
            #The center of curvature lies on the line AO
            #(or its extension) or vector .
            pass
        elif self.geom == 2:
            # The tangent of centroid arc at end A is
            # parallel to line AO or vector . Point O (or
            # vector) and the arc must be on the
            # same side of the chord .
            pass
        elif self.geom == 3:
            # The bend radius (RB) is specified on the
            # PBEND entry: Points A, B, and O (or
            # vector ) define a plane parallel or
            # coincident with the plane of the element
            # arc. Point O (or vector ) lies on the
            # opposite side of line AB from the center of
            # the curvature.
            pass
        elif self.geom == 4:
            # THETAB is specified on the PBEND entry.
            # Points A, B, and O (or vector ) define a
            # plane parallel or coincident with the plane
            # of the element arc. Point O (or vector )
            # lies on the opposite side of line AB from the
            # center of curvature.
            pass
        else:
            raise RuntimeError('geom=%r is not supported on the CBEND' % geom)
        return L

    def _validate_input(self):
        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    @property
    def nodes(self):
        return [self.ga, self.gb]

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        else:
            return self.ga_ref.nid

    def Gb(self):
        if isinstance(self.gb, integer_types):
            return self.gb
        else:
            return self.gb_ref.nid

    def Area(self):
        if isinstance(self.pid, integer_types):
            msg = 'Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self))
            raise RuntimeError(msg)
        return self.pid_ref.Area()

    def _verify(self, xref):
        edges = self.get_edge_ids()

    def raw_fields(self):
        (x1, x2, x3) = self.getX_G0_defaults()
        list_fields = ['CBEND', self.eid, self.Pid(), self.Ga(), self.Gb(),
                       x1, x2, x3, self.geom]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        #self.g0 = model.nodes[self.g0]
        self.ga_ref = self.ga
        self.gb_ref = self.gb
        self.pid_ref = self.pid

    def uncross_reference(self):
        node_ids = self.node_ids
        self.ga = node_ids[0]
        self.gb = node_ids[1]
        self.pid = self.Pid()
        del self.ga_ref, self.gb_ref, self.pid_ref

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        else:
            return self.comment + print_card_16(card)
