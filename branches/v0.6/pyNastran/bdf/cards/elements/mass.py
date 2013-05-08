## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
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
        assert c2 is None or isinstance(c1, int), 'c1=%r' % c1
        assert c2 is None or isinstance(c2, int), 'c2=%r' % c2

    def cross_reference(self, model):
        msg = ' which is required by CMASS1 eid=%s' % self.eid
        #self.g1 = model.Node(self.g1, msg=msg)
        #self.g2 = model.Node(self.g2, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def Pid(self):
        if isinstance(self.pid, int):
            return self.pid
        return self.pid.pid

    def rawFields(self):
        fields = ['CMASS1', self.eid, self.Pid(), self.g1, self.c1,
                  self.g2, self.c2]
        return fields


class CMASS2(PointMassElement):
    """
    Defines a scalar mass element without reference to a property entry.
    CMASS2 EID M G1 C1 G2 C2
    """
    type = 'CMASS2'

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


class CMASS3(PointMassElement):
    """
    Defines a scalar mass element that is connected only to scalar points.
    CMASS3 EID PID S1 S2
    """
    type = 'CMASS3'

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


class CMASS4(PointMassElement):
    """
    Defines a scalar mass element that is connected only to scalar points,
    without reference to a property entry
    CMASS4 EID M S1 S2
    """
    type = 'CMASS4'

    def __init__(self, card=None, data=None, comment=''):
        print("***")
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


class CONM1(PointMassElement):
    type = 'CONM1'

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


class CONM2(PointMassElement):
    """
    :param self: the CONM2 object
    :param eid:  element ID
    :param nid:  node ID
    :param cid:  coordinate frame of the offset (-1=absolute coordinates)
    :param X:    offset vector
    :param I:    mass moment of inertia matrix about the CG
    
    ::

      CONM2    501274  11064          132.274
    """
    type = 'CONM2'

    def __init__(self, card=None, data=None, comment=''):
        PointMassElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.nid = integer(card, 2, 'nid')
            self.cid = integer_or_blank(card, 3, 'cid', 0)
            self.mass = double_or_blank(card, 4, 'mass', 0.)
            self.X = array([double_or_blank(card, 5, 'x1', 0.0),
                            double_or_blank(card, 6, 'x2', 0.0),
                            double_or_blank(card, 7, 'x3', 0.0)])

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
            self.X = data[4:7]
            self.I = data[7:]

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
            assert isinstance(c[i], float), 'centroid[%i]=%r' %(i, c[i])

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
        A = [[I[0], I[1], I[3]],
             [I[1], I[2], I[4]],
             [I[3], I[4], I[5]]]
        if self.Cid() == -1:
            return A
        else:
            # transform to global
            (dx, matrix) = self.cid.transformToGlobal(self.X)
            raise NotImplementedError('CONM2 intertia method is not implemented.')
            A2 = A * matrix
            return A2  # correct for offset using dx???

    def Centroid(self):
        if self.Cid() == 0:
            X2 = self.nid.Position() + self.X
        elif self.Cid() == -1:
            return self.X
        else:
            dx, matrix = self.cid.transformToGlobal(self.X)
            #print "dx = ",dx
            #print "X  = ",self.nid.Position()
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