# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All mass elements are defined in this file.  This includes:

 * CMASS1
 * CMASS2
 * CMASS3
 * CMASS4
 * CONM1
 * CONM2

All mass elements are PointMassElement and Element objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import zeros, array

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element #, BaseCard
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
                                       double_or_blank)
from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16


class PointElement(Element):
    def __init__(self, card, data):
        Element.__init__(self, card, data)


class PointMassElement(PointElement):
    def __init__(self, card, data):
        self.mass = None
        PointElement.__init__(self, card, data)

    def Mass(self):
        return self.mass


# class PointMass(BaseCard):
#     def __init__(self, card, data):
#         self.mass = None
#
#     def Mass(self):
#         return self.mass


class CMASS1(PointMassElement):
    """
    Defines a scalar mass element.
    CMASS2 EID M   G1 C1 G2 C2
    CMASS1 EID PID G1 C1 G2 C2
    """
    type = 'CMASS1'
    _field_map = {
        1: 'eid', 2:'pid', 3:'g1', 4:'c1', 5:'g2', 6:'c2',
    }

    def __init__(self, card=None, data=None, comment=''):
        PointMassElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.g1 = integer_or_blank(card, 3, 'g1')
            self.c1 = integer_or_blank(card, 4, 'c1')
            self.g2 = integer_or_blank(card, 5, 'g2')
            self.c2 = integer_or_blank(card, 6, 'c2')
            assert len(card) <= 7, 'len(CMASS1 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.g1 = data[2]
            self.g2 = data[3]
            self.c1 = data[4]
            self.c2 = data[5]

    def Eid(self):
        return self.eid

    def Mass(self):
        return self.pid.mass

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        mass = self.Mass()
        c1 = self.c1
        c2 = self.c2
        #self.nodes

        assert isinstance(eid, int), 'eid=%r' % eid
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mass, float), 'mass=%r' % mass
        assert c1 is None or isinstance(c1, int), 'c1=%r' % c1
        assert c2 is None or isinstance(c2, int), 'c2=%r' % c2

    def cross_reference(self, model):
        msg = ' which is required by CMASS1 eid=%s' % self.eid
        if isinstance(self.g1, int):
            self.g1 = model.Node(self.g1, msg=msg)
        if isinstance(self.g2, int):
            self.g2 = model.Node(self.g2, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def G1(self):
        if isinstance(self.g1, int):
            return self.g1
        elif self.g1 is None:
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, int):
            return self.g2
        elif self.g2 is None:
            return self.g2
        return self.g2.nid

    def Centroid(self):
        """
        Centroid is assumed to be c=(g1+g2)/2.
        If g2 is blank, then the centroid is the location of g1.
        """
        f = 0.
        p1 = array([0., 0., 0.])
        p2 = array([0., 0., 0.])
        if self.g1 is not None:
            p1 = self.g1.Position()
            f += 1.
        if self.g2 is not None:
            p2 = self.g2.Position()
            f += 1.
        c = (p1 + p2) / f
        return c

    def nodeIDs(self):
        g1 = self.G1()
        g2 = self.G2()
        nodes = []
        if g1:
            nodes.append(g1)
        if g2:
            nodes.append(g2)
        return nodes

    def Pid(self):
        if isinstance(self.pid, int):
            return self.pid
        return self.pid.pid

    def rawFields(self):
        fields = ['CMASS1', self.eid, self.Pid(), self.g1, self.c1,
                  self.g2, self.c2]
        return fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


class CMASS2(PointMassElement):
    """
    Defines a scalar mass element without reference to a property entry.
    CMASS2 EID M G1 C1 G2 C2
    """
    type = 'CMASS2'
    _field_map = {
        1: 'eid', 2:'pid', 3:'g1', 4:'c1', 5:'g2', 6:'c2',
    }

    def __init__(self, card=None, data=None, comment=''):
        PointMassElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.mass = double_or_blank(card, 2, 'mass', 0.)
            self.g1 = integer_or_blank(card, 3, 'g1')
            self.c1 = integer_or_blank(card, 4, 'c1')
            self.g2 = integer_or_blank(card, 5, 'g2')
            self.c2 = integer_or_blank(card, 6, 'c2')
            assert len(card) <= 7, 'len(CMASS2 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.mass = data[1]
            self.g1 = data[2]
            self.g2 = data[3]
            self.c1 = data[4]
            self.c2 = data[5]

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        mass = self.Mass()
        c1 = self.c1
        c2 = self.c2
        #self.nodes

        assert isinstance(eid, int), 'eid=%r' % eid
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mass, float), 'mass=%r' % mass
        assert c1 is None or isinstance(c1, int), 'c1=%r' % c1
        assert c2 is None or isinstance(c2, int), 'c2=%r' % c2

    def Eid(self):
        return self.eid

    def nodeIDs(self):
        g1 = self.G1()
        g2 = self.G2()
        nodes = []
        if g1:
            nodes.append(g1)
        if g2:
            nodes.append(g2)
        return nodes

    def Mass(self):
        return self.mass

    def Centroid(self):
        """
        Centroid is assumed to be c=(g1+g2)/2.
        If g2 is blank, then the centroid is the location of g1.
        """
        f = 0.
        p1 = array([0., 0., 0.])
        p2 = array([0., 0., 0.])
        if self.g1 is not None:
            p1 = self.g1.Position()
            f += 1.
        if self.g2 is not None:
            p2 = self.g2.Position()
            f += 1.
        c = (p1 + p2) / f
        return c

    def cross_reference(self, model):
        msg = ' which is required by CMASS2 eid=%s' % self.eid
        if isinstance(self.g1, int):
            self.g1 = model.Node(self.g1, msg=msg)
        if isinstance(self.g2, int):
            self.g2 = model.Node(self.g2, msg=msg)

    def G1(self):
        if isinstance(self.g1, int):
            return self.g1
        elif self.g1 is None:
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, int):
            return self.g2
        elif self.g2 is None:
            return self.g2
        return self.g2.nid

    def rawFields(self):
        fields = ['CMASS2', self.eid, self.mass, self.G1(),
                  self.c1, self.G2(), self.c2]
        return fields

    def reprFields(self):
        mass = set_blank_if_default(self.mass, 0.)
        fields = ['CMASS2', self.eid, mass, self.G1(), self.c1,
                  self.G2(), self.c2]
        return fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


class CMASS3(PointMassElement):
    """
    Defines a scalar mass element that is connected only to scalar points.
    CMASS3 EID PID S1 S2
    """
    type = 'CMASS3'
    _field_map = {
        1: 'eid', 2:'pid', 3:'s1', 4:'s2',
    }

    def __init__(self, card=None, data=None, comment=''):
        PointMassElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.s1 = integer_or_blank(card, 3, 's1')
            self.s2 = integer_or_blank(card, 4, 's2')
            assert len(card) <= 5, 'len(CMASS3 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.s1 = data[2]
            self.s2 = data[3]
        assert self.s1 != self.s2

    def Eid(self):
        return self.eid

    def Mass(self):
        return self.pid.mass

    def isSameCard(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.Pid(), self.s1, self.s2]
        fields2 = [elem.eid, elem.Pid(), elem.s1, elem.s2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def cross_reference(self, model):
        """
        Links up the propertiy ID
        """
        msg = ' which is required by CMASS3 eid=%s' % self.eid
        #self.s1 = model.Node(self.s1, msg=msg)
        #self.s2 = model.Node(self.s2, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def rawFields(self):
        fields = ['CMASS3', self.eid, self.Pid(), self.s1, self.s2]
        return fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


class CMASS4(PointMassElement):
    """
    Defines a scalar mass element that is connected only to scalar points,
    without reference to a property entry
    CMASS4 EID M S1 S2
    """
    type = 'CMASS4'
    _field_map = {
        1: 'eid', 2:'mass', 3:'s1', 4:'s2',
    }

    def __init__(self, card=None, data=None, comment=''):
        PointMassElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.mass = double_or_blank(card, 2, 'mass', 0.)
            self.s1 = integer(card, 3, 's1')
            self.s2 = integer_or_blank(card, 4, 's2', 0)
            assert len(card) <= 5, 'len(CMASS4 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.mass = data[1]
            self.s1 = data[2]
            self.s2 = data[3]
        assert self.s1 != self.s2

    def Eid(self):
        return self.eid

    def Mass(self):
        return self.mass

    def isSameCard(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.Pid(), self.s1, self.s2]
        fields2 = [elem.eid, elem.Pid(), elem.s1, elem.s2]
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def cross_reference(self, model):
        """
        """
        #self.s1 = model.Node(self.s1)
        #self.s2 = model.Node(self.s2)
        pass

    def rawFields(self):
        fields = ['CMASS4', self.eid, self.mass, self.s1, self.s2]
        return fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


class CONM1(PointMassElement):
    type = 'CONM1'
    _field_map = {
        1: 'eid', 2:'nid', 3:'cid',
    }
    def _update_field_helper(self, n, value):
        m = self.massMatrix
        if n == 4:
            m[0, 0] = value
        elif n == 5:
            m[1, 0] = value
        elif n == 6:
            m[1, 1] = value
        elif n == 7:
            m[2, 0] = value
        elif n == 8:
            m[2, 1] = value
        elif n == 9:
            m[2, 2] = value
        elif n == 10:
            m[3, 0] = value
        elif n == 11:
            m[3, 1] = value
        elif n == 12:
            m[3, 2] = value
        elif n == 13:
            m[3, 3] = value
        elif n == 14:
            m[4, 0] = value
        elif n == 15:
            m[4, 1] = value
        elif n == 16:
            m[4, 2] = value
        elif n == 17:
            m[4, 3] = value
        elif n == 18:
            m[4, 4] = value
        elif n == 19:
            m[5, 0] = value
        elif n == 20:
            m[5, 1] = value
        elif n == 21:
            m[5, 2] = value
        elif n == 22:
            m[5, 3] = value
        elif n == 23:
            m[5, 4] = value
        elif n == 24:
            m[5, 5] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        """
        Concentrated Mass Element Connection, General Form
        Defines a 6 x 6 symmetric mass matrix at a geometric grid point::

          CONM1 EID G CID M11 M21 M22 M31 M32
                M33 M41 M42 M43 M44 M51 M52 M53
                M54 M55 M61 M62 M63 M64 M65 M66

        ::

          [M] = [M11 M21 M31 M41 M51 M61]
                [    M22 M32 M42 M52 M62]
                [        M33 M43 M53 M63]
                [            M44 M54 M64]
                [    Sym         M55 M65]
                [                    M66]
        """
        PointMassElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        m = zeros((6, 6))
        if card:
            #self.nids  = [ card[1] ]
            #del self.nids
            #self.pid = None
            self.eid = integer(card, 1, 'eid')
            self.nid = integer(card, 2, 'nid')
            self.cid = integer_or_blank(card, 3, 'cid', 0)

            m[0, 0] = double_or_blank(card, 4, 'M11', 0.)
            m[1, 0] = double_or_blank(card, 5, 'M21', 0.)
            m[1, 1] = double_or_blank(card, 6, 'M22', 0.)
            m[2, 0] = double_or_blank(card, 7, 'M31', 0.)
            m[2, 1] = double_or_blank(card, 8, 'M32', 0.)
            m[2, 2] = double_or_blank(card, 9, 'M33', 0.)
            m[3, 0] = double_or_blank(card, 10, 'M41', 0.)
            m[3, 1] = double_or_blank(card, 11, 'M42', 0.)
            m[3, 2] = double_or_blank(card, 12, 'M43', 0.)
            m[3, 3] = double_or_blank(card, 13, 'M44', 0.)
            m[4, 0] = double_or_blank(card, 14, 'M51', 0.)
            m[4, 1] = double_or_blank(card, 15, 'M52', 0.)
            m[4, 2] = double_or_blank(card, 16, 'M53', 0.)
            m[4, 3] = double_or_blank(card, 17, 'M54', 0.)
            m[4, 4] = double_or_blank(card, 18, 'M55', 0.)
            m[5, 0] = double_or_blank(card, 19, 'M61', 0.)
            m[5, 1] = double_or_blank(card, 20, 'M62', 0.)
            m[5, 2] = double_or_blank(card, 21, 'M63', 0.)
            m[5, 3] = double_or_blank(card, 22, 'M64', 0.)
            m[5, 4] = double_or_blank(card, 23, 'M65', 0.)
            m[5, 5] = double_or_blank(card, 24, 'M66', 0.)
            assert len(card) <= 25, 'len(CONM1 card) = %i' % len(card)
        else:
            (eid, nid, cid, m1, m2a, m2b, m3a, m3b, m3c, m4a, m4b, m4c, m4d,
             m5a, m5b, m5c, m5d, m5e, m6a, m6b, m6c, m6d, m6e, m6f) = data
            self.eid = eid
            self.nid = nid
            self.cid = cid
            m[0, 0] = m1   # M11
            m[1, 0] = m2a  # M21
            m[1, 1] = m2b  # M22
            m[2, 0] = m3a  # M31
            m[2, 1] = m3b  # M32
            m[2, 2] = m3c  # M33
            m[3, 0] = m4a  # M41
            m[3, 1] = m4b  # M42
            m[3, 2] = m4c  # M43
            m[3, 3] = m4d  # M44
            m[4, 0] = m5a  # M51
            m[4, 1] = m5b  # M52
            m[4, 2] = m5c  # M53
            m[4, 3] = m5d  # M54
            m[4, 4] = m5e  # M55
            m[5, 0] = m6a  # M61
            m[5, 1] = m6b  # M62
            m[5, 2] = m6c  # M63
            m[5, 3] = m6d  # M64
            m[5, 4] = m6e  # M65
            m[5, 5] = m6f  # M66
        self.massMatrix = m

    def _verify(self, xref=False):
        eid = self.Eid()
        assert isinstance(eid, int), 'eid=%r' % eid

    def Eid(self):
        return self.eid

    def Mass(self):
        return 0.0

    def nodeIDs(self):
        return [self.Nid()]

    def Nid(self):
        if isinstance(self.nid, int):
            return self.nid
        return self.nid.nid

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def cross_reference(self, model):
        msg = ' which is required by CONM1 eid=%s' % self.eid
        self.nid = model.Node(self.nid, msg=msg)
        self.cid = model.Coord(self.cid, msg=msg)

    def MassMatrix(self):
        return self.massMatrix

    def rawFields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        nid = self.Nid()
        m = self.massMatrix
        list_fields = ['CONM1', self.eid, nid, cid, m[0, 0], m[1, 0], m[1, 1],
                  m[2, 0], m[2, 1], m[2, 2], m[3, 0], m[3, 1], m[3, 2],
                  m[3, 3], m[4, 0], m[4, 1], m[4, 2], m[4, 3], m[4, 4],
                  m[5, 0], m[5, 1], m[5, 2], m[5, 3], m[5, 4], m[5, 5]]
        return list_fields

    def reprFields(self):
        list_fields = self.rawFields()
        list_fields2 = list_fields[0:4]
        for field in list_fields[4:]:
            val = set_blank_if_default(field, 0.)
            list_fields2.append(val)
        return list_fields2

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


class CONM2(PointMassElement):
    type = 'CONM2'
    _field_map = {
        1: 'eid', 2:'nid', 3:'cid', 4:'mass',
    }
    def _update_field_helper(self, n, value):
        if n == 5:
            self.X[0] = value
        elif n == 6:
            self.X[1] = value
        elif n == 7:
            self.X[2] = value

        elif n == 9:
            self.I[0] = value
        elif n == 10:
            self.I[1] = value
        elif n == 11:
            self.I[2] = value
        elif n == 12:
            self.I[3] = value
        elif n == 13:
            self.I[4] = value
        elif n == 14:
            self.I[5] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        """
        :param self: the CONM2 object
        :param eid:  element ID
        :param nid:  node ID
        :param cid:  coordinate frame of the offset (-1=absolute coordinates)
        :param X:    offset vector relative to nid
        :param I:    mass moment of inertia matrix about the CG

        +-------+--------+-------+-------+---------+------+------+------+------+
        | CONM2 | EID    | NID   |  CID  |  MASS   |  X1  |  X2  |  X3  |      |
        +-------+--------+-------+-------+---------+------+------+------+------+
        |       |   I11  |   I21 |  I22  |   I31   |  I32 |  I33 |
        +-------+--------+-------+-------+---------+------+------+

        +-------+--------+-------+-------+---------+
        | CONM2 | 501274 | 11064 |       | 132.274 |
        +-------+--------+-------+-------+---------+
        """
        PointMassElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element identification number. (0 < Integer < 100,000,000)
            self.eid = integer(card, 1, 'eid')

            #: Grid point identification number. (Integer > 0)
            self.nid = integer(card, 2, 'nid')
            #: Coordinate system identification number.
            #: For CID of -1; see X1, X2, X3 below.
            #: (Integer > -1; Default = 0)
            self.cid = integer_or_blank(card, 3, 'cid', 0)

            #: Mass value. (Real)
            self.mass = double_or_blank(card, 4, 'mass', 0.)

            #: Offset distances from the grid point to the center of gravity of
            # the mass in the coordinate system defined in field 4, unless
            # CID = -1, in which case X1, X2, X3 are the coordinates, not
            # offsets, of the center of gravity of the mass in the basic
            #: coordinate system. (Real)
            self.X = array([double_or_blank(card, 5, 'x1', 0.0),
                            double_or_blank(card, 6, 'x2', 0.0),
                            double_or_blank(card, 7, 'x3', 0.0)])

            #: Mass moments of inertia measured at the mass center of gravity in
            # the coordinate system defined by field 4. If CID = -1, the basic
            # coordinate system is implied. (Real)
            self.I = array([double_or_blank(card, 9, 'I11', 0.0),
                            double_or_blank(card, 10, 'I21', 0.0),
                            double_or_blank(card, 11, 'I22', 0.0),
                            double_or_blank(card, 12, 'I31', 0.0),
                            double_or_blank(card, 13, 'I32', 0.0),
                            double_or_blank(card, 14, 'I33', 0.0)])
            assert len(card) <= 15, 'len(CONM2 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.nid = data[1]
            self.cid = data[2]
            self.mass = data[3]
            self.X = array(data[4:7])
            self.I = array(data[7:])

    def _verify(self, xref=False):
        eid = self.Eid()
        nid = self.Nid()
        cid = self.Cid()
        mass = self.Mass()
        c = self.Centroid()

        assert isinstance(eid, int), 'eid=%r' % eid
        assert isinstance(nid, int), 'nid=%r' % nid
        assert isinstance(cid, int), 'cid=%r' % cid
        assert isinstance(mass, float), 'mass=%r' % mass
        for i in range(3):
            assert isinstance(c[i], float), 'centroid[%i]=%r' % (i, c[i])

    def Eid(self):
        return self.eid

    def Mass(self):
        return self.mass

    def Inertia(self):
        """
        Returns the 3x3 inertia matrix
        .. warning:: doesnt handle offsets or coordinate systems
        """
        I = self.I
        A = [[ I[0], -I[1], -I[3]],
             [-I[1],  I[2], -I[4]],
             [-I[3], -I[4],  I[5]]]
        if self.Cid() in [0, -1]:
            return A
        else:
            # transform to global
            (dx, matrix) = self.cid.transformToGlobal(self.X)
            raise NotImplementedError('CONM2 intertia method is not implemented.')
            A2 = A * matrix
            return A2  # correct for offset using dx???

    def Centroid(self):
        """
        This method seems way more complicated than it needs to be thanks
        to all these little caveats that don't seem to be supported.
        """
        if self.Cid() == 0:
            # no transform needed
            X2 = self.nid.Position() + self.X
        elif self.Cid() == -1:
            # case X1, X2, X3 are the coordinates, not offsets, of the center of gravity of
            # the mass in the basic coordinate system.

            # 4. If CID = -1, offsets are internally computed as the difference between the grid
            # point location and X1, X2, X3.
            # this statement is not supported...
            return self.X
        else:
            # Offset distances from the grid point to the center of gravity of the mass
            # in the coordinate system

            # If CID > 0, then X1, X2, and X3 are defined by a local Cartesian system, even
            # if CID references a spherical or cylindrical coordinate system. This is similar
            # to the manner in which displacement coordinate systems are defined.
            # this statement is not supported...

            # convert self.X into the global frame
            x, matrix = self.cid.transformToGlobal(self.X)

            # self.X is an offset
            dx = x - self.cid.origin

            # the actual position of the CONM2
            X2 = self.nid.Position() + dx
        return X2

    def cross_reference(self, model):
        msg = ' which is required by CONM2 eid=%s' % self.eid
        if self.Cid() != -1:
            self.cid = model.Coord(self.cid, msg=msg)
        self.nid = model.Node(self.nid, msg=msg)

    def nodeIDs(self):
        return [self.Nid()]

    def Nid(self):
        if isinstance(self.nid, int):
            return self.nid
        return self.nid.nid

    def Cid(self):
        if isinstance(self.cid, int):
            return self.cid
        return self.cid.cid

    def writeCodeAster(self):
        msg = ''
        msg += "    DISCRET=_F(\n"
        msg += "             'CARA='M_T_D_N'\n"
        msg += "              NOEUD=N%s\n" % self.Nid()
        msg += "              VALE=%g),\n" % self.mass
        return msg

    def rawFields(self):
        list_fields = (['CONM2', self.eid, self.Nid(), self.Cid(), self.mass] +
                  list(self.X) + [None] + list(self.I))
        return list_fields

    def reprFields(self):
        I = []
        for i in self.I:
            if i == 0.:
                I.append(None)
            else:
                I.append(i)
        X = []
        for x in self.X:
            if x == 0.:
                X.append(None)
            else:
                X.append(x)

        cid = set_blank_if_default(self.Cid(), 0)
        list_fields = (['CONM2', self.eid, self.Nid(), cid, self.mass] + X +
                  [None] + I)
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)
