# pylint: disable=R0904,R0902,E1101,E1103,C0111,C0302,C0103,W0101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types, integer_types

from numpy import array
from numpy.linalg import norm

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    integer_double_or_blank, double_or_blank, string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class LineElement(Element):  # CBAR, CBEAM, CBEAM3, CBEND
    def __init__(self, card, data):
        Element.__init__(self, card, data)

    def C(self):
        """torsional constant"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.C()

    def Area(self):
        """returns the area of the element face"""
        raise NotImplementedError('implement self.Area() for %s' % self.type)

    def E(self):
        """returns the Young's Modulus, :math:`E`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.E()

    def G(self):
        """returns the Shear Modulus, :math:`G`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.G()

    def J(self):
        """returns the Polar Moment of Inertia, :math:`J`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.J()

    def I11(self):
        """returns the Moment of Inertia, :math:`I_{11}`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.I11()

    def I22(self):
        """returns the Moment of Inertia, :math:`I_{22}`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.I22()

    def I12(self):
        """returns the Moment of Inertia, :math:`I_{12}`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.I12()

    def Nu(self):
        """Get Poisson's Ratio, :math:`\nu`"""
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.nu

    def Rho(self):
        """Get the material density, :math:`\rho`"""
        #print(str(self.pid), type(self.pid))
        #raise NotImplementedError('implement self.Rho() for %s' % self.type)
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.rho

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
        return self.pid.MassPerLength()

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
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        #self.g0 = model.nodes[self.g0]

    def uncross_reference(self):
        self.nodes = self.node_ids
        self.pid = self.Pid()

    def Length(self):
        r"""
        Gets the length, :math:`L`, of the element.

        .. math:: L = \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  }

        :param self: the object pointer
        """
        L = norm(self.nodes[1].get_position() - self.nodes[0].get_position())
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
            self.x = array([field5,
                            double_or_blank(card, 6, 'x2', 0.0),
                            double_or_blank(card, 7, 'x3', 0.0)], dtype='float64')
        self.offt = string_or_blank(card, 8, 'offt', 'GGG')
        assert len(card) <= 9, 'len(CBAROR card) = %i' % len(card)


class CBAR(LineElement):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    | CBAR  | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B  |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    | CBAR  | EID | PID | GA  | GB  | G0  |     |     | OFFT |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B  |
    +-------+-----+-----+-----+-----+-----+-----+-----+------+

    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    |  CBAR | 2     |  39 | 7     | 6     |  105   |       |       |  GGG  |
    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    |       |       | 513 | 0.0+0 | 0.0+0 |    -9. | 0.0+0 | 0.0+0 |   -9. |
    +-------+-------+-----+-------+-------+--------+-------+-------+-------+
    """
    type = 'CBAR'
    asterType = 'CBAR'
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

    def __init__(self, card=None, data=None, comment=''):
        LineElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.ga = integer(card, 3, 'ga')
            self.gb = integer(card, 4, 'gb')
            self.initX_G0(card)

            self.offt = string_or_blank(card, 8, 'offt', 'GGG')
            #print('self.offt = |%s|' % (self.offt))

            self.pa = integer_or_blank(card, 9, 'pa', 0)
            self.pb = integer_or_blank(card, 10, 'pb', 0)

            self.wa = array([double_or_blank(card, 11, 'w1a', 0.0),
                             double_or_blank(card, 12, 'w2a', 0.0),
                             double_or_blank(card, 13, 'w3a', 0.0)], dtype='float64')

            self.wb = array([double_or_blank(card, 14, 'w1b', 0.0),
                             double_or_blank(card, 15, 'w2b', 0.0),
                             double_or_blank(card, 16, 'w3b', 0.0)], dtype='float64')
            assert len(card) <= 17, 'len(CBAR card) = %i' % len(card)
        else:  #: .. todo:: verify
            #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,g0]]
            #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,x1,x2,x3]]

            main = data[0]
            flag = data[1][0]
            if flag in [0, 1]:
                self.g0 = None
                self.x = array([data[1][1],
                                data[1][2],
                                data[1][3]], dtype='float64')
            else:
                self.g0 = data[1][1]
                self.x = None

            self.eid = main[0]
            self.pid = main[1]
            self.ga = main[2]
            self.gb = main[3]
            #self.offt = str(data[4]) # GGG
            self.offt = 'GGG'  #: .. todo:: offt can be an integer; translate to char
            self.pa = main[4]
            self.pb = main[5]

            self.wa = array([main[6], main[7], main[8]], dtype='float64')
            self.wb = array([main[9], main[10], main[11]], dtype='float64')

        if not isinstance(self.offt, string_types):
            raise SyntaxError('invalid offt expected a string of length 3 '
                              'offt=|%r|; Type=%s' % (self.offt, type(self.offt)))

        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)

        msg = 'invalid offt parameter of %s...offt=%s' % (self.type, self.offt)
        # B,G,O
        assert self.offt[0] in ['G', 'B'], msg
        assert self.offt[1] in ['G', 'O', 'E'], msg
        assert self.offt[2] in ['G', 'O', 'E'], msg

    def _verify(self, xref=False):
        pid = self.Pid()
        edges = self.get_edge_ids()
        if xref:  # True
            mid = self.Mid()
            A = self.Area()
            nsm = self.Nsm()
            mpl = self.MassPerLength()
            L = self.Length()
            mass = self.Mass()
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(A, float), 'A=%r' % A
        assert isinstance(L, float), 'L=%r' % L
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        assert isinstance(mpl, float), 'mass_per_length=%r' % mpl
        assert isinstance(mass, float), 'nass=%r' % mass

    def Mid(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.Mid()

    def Area(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        A = self.pid.Area()
        assert isinstance(A, float)
        return A

    def J(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        j = self.pid.J()
        assert isinstance(j, float), 'J=%r for CBAR eid=%s pid=%s pidType=%s' % (j, self.eid, self.pid.pid, self.pid.type)
        return j

    def Length(self):
        # TODO: consider w1a and w1b in the length formulation
        L = norm(self.gb.get_position() - self.ga.get_position())
        assert isinstance(L, float)
        return L

    def Nsm(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        nsm = self.pid.Nsm()
        assert isinstance(nsm, float)
        return nsm

    def I1(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.I1()

    def I2(self):
        return self.pid.I2()

    def Centroid(self):
        return (self.ga.get_position() + self.gb.get_position()) / 2.

    def initX_G0(self, card):
        field5 = integer_double_or_blank(card, 5, 'g0_x1', 0.0)
        if isinstance(field5, integer_types):
            self.g0 = field5
            self.x = None
        elif isinstance(field5, float):
            self.g0 = None
            self.x = array([field5,
                            double_or_blank(card, 6, 'x2', 0.0),
                            double_or_blank(card, 7, 'x3', 0.0)], dtype='float64')
            if norm(self.x) == 0.0:
                msg = 'G0 vector defining plane 1 is not defined.\n'
                msg += 'G0 = %s\n' % self.g0
                msg += 'X  = %s\n' % self.x
                raise RuntimeError(msg)
        else:
            msg = ('field5 on %s (G0/X1) is the wrong type...id=%s field5=%s '
                   'type=%s' % (self.type, self.eid, field5, type(field5)))
            raise RuntimeError(msg)

    def cross_reference(self, model):
        #if self.g0:
        #    self.x = nodes[self.g0].get_position() - nodes[self.ga].get_position()
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.Ga(), msg=msg)
        self.gb = model.Node(self.Gb(), msg=msg)
        self.pid = model.Property(self.Pid(), msg=msg)

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        else:
            return self.ga.nid

    def Gb(self):
        if isinstance(self.gb, integer_types):
            return self.gb
        else:
            return self.gb.nid

    def getX_G0_defaults(self):
        if self.g0:
            return (self.g0, None, None)
        else:
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
        offt = set_blank_if_default(self.offt, 'GGG')
        list_fields = ['CBAR', self.eid, self.Pid(), self.Ga(), self.Gb(), x1, x2,
                       x3, offt, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class CBEAM3(CBAR):
    """
    Defines a three-node beam element
    """
    type = 'CBEAM3'

    def __init__(self, card=None, data=None, comment=''):
        LineElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.ga = integer(card, 3, 'ga')
            self.gb = integer(card, 4, 'gb')
            self.gc = integer(card, 5, 'gc')

            self.initX_G0(card)

            self.wa = array([double_or_blank(card, 9, 'w1a', 0.0),
                             double_or_blank(card, 10, 'w2a', 0.0),
                             double_or_blank(card, 11, 'w3a', 0.0)], dtype='float64')

            self.wb = array([double_or_blank(card, 12, 'w1b', 0.0),
                             double_or_blank(card, 13, 'w2b', 0.0),
                             double_or_blank(card, 14, 'w3b', 0.0)], dtype='float64')

            self.wc = array([double_or_blank(card, 15, 'w1c', 0.0),
                             double_or_blank(card, 16, 'w2c', 0.0),
                             double_or_blank(card, 17, 'w3c', 0.0)], dtype='float64')

            self.tw = array([double_or_blank(card, 18, 0., 'twa'),
                             double_or_blank(card, 19, 0., 'twb'),
                             double_or_blank(card, 20, 0., 'twc')], dtype='float64')

            self.s = array([integer_or_blank(card, 21, 'sa'),
                            integer_or_blank(card, 22, 'sb'),
                            integer_or_blank(card, 23, 'sc')], dtype='float64')
            assert len(card) <= 24, 'len(CBEAM3 card) = %i' % len(card)
        else:
            raise NotImplementedError(data)

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)
        self.gc = model.Node(self.gc, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def Length(self):
        """
        .. math:: L = g_b - g_a
        """
        L = norm(self.gb.get_position() - self.ga.get_position())
        return L

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

    def __init__(self, card=None, data=None, comment=''):
        LineElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.ga = integer(card, 3, 'ga')
            self.gb = integer(card, 4, 'gb')
            x1Go = integer_double_or_blank(card, 5, 'x1_g0', 0.0)
            if isinstance(x1Go, integer_types):
                self.g0 = x1Go
                self.x = None
            elif isinstance(x1Go, float):
                self.g0 = None
                self.x = array([double_or_blank(card, 5, 'x1', 0.0),
                                double_or_blank(card, 6, 'x2', 0.0),
                                double_or_blank(card, 7, 'x3', 0.0)], dtype='float64')
                if norm(self.x) == 0.0:
                    msg = 'G0 vector defining plane 1 is not defined.\n'
                    msg += 'G0 = %s\n' % self.g0
                    msg += 'X  = %s\n' % self.x
                    raise RuntimeError(msg)
            else:
                raise ValueError('invalid x1Go=|%s| on CBEND' % x1Go)
            self.geom = integer(card, 8, 'geom')

            assert len(card) == 9, 'len(CBEND card) = %i' % len(card)
            assert self.geom in [1, 2, 3, 4], 'geom is invalid geom=|%s|' % self.geom
        else:
            raise NotImplementedError(data)

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
            return self.ga.nid

    def Gb(self):
        if isinstance(self.gb, integer_types):
            return self.gb
        else:
            return self.gb.nid

    def Area(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.Area()

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
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        #self.pid = model.Property(self.pid, msg=msg)
        #self.g0 = model.nodes[self.g0]

    def write_card(self, size, is_double):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        else:
            return self.comment + print_card_16(card)
