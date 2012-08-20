# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import zeros, array

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element, BaseCard


class PointElement(Element):
    def __init__(self, card=None, data=None):
        Element.__init__(self, card, data)


class PointMassElement(PointElement):
    def __init__(self, card=None, data=None):
        self.mass = None
        PointElement.__init__(self, card, data)

    def Mass(self):
        return self.mass


class PointMass(BaseCard):
    def __init__(self, card=None, data=None):
        self.mass = None

    def Mass(self):
        return self.mass


class CMASS1(PointMass):
    """
    Defines a scalar mass element.
    CMASS2 EID M   G1 C1 G2 C2
    CMASS1 EID PID G1 C1 G2 C2
    """
    type = 'CMASS1'

    def __init__(self, card=None, data=None):
        PointMass.__init__(self, card, data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2, self.eid)
            self.g1 = card.field(3)
            self.c1 = card.field(4)
            self.g2 = card.field(5)
            self.c2 = card.field(6)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.g1 = data[2]
            self.g2 = data[3]
            self.c1 = data[4]
            self.c2 = data[5]
        ###

    def Mass(self):
        return self.pid.mass

    def cross_reference(self, model):
        #self.g1 = model.Node(self.g1)
        #self.g2 = model.Node(self.g2)
        self.pid = model.Property(self.pid)

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

    def __init__(self, card=None, data=None):
        PointMassElement.__init__(self, card, data)
        if card:
            self.eid = card.field(1)
            self.mass = card.field(2, 0.)
            self.g1 = card.field(3)
            self.c1 = card.field(4)
            self.g2 = card.field(5)
            self.c2 = card.field(6)
        else:
            self.eid = data[0]
            self.mass = data[1]
            self.g1 = data[2]
            self.g2 = data[3]
            self.c1 = data[4]
            self.c2 = data[5]
        ###

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

    def cross_reference(self, mesh):
        if isinstance(self.g1, int):
            self.g1 = mesh.Node(self.g1)
        if isinstance(self.g2, int):
            self.g2 = mesh.Node(self.g2)
        ###

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
        #print "cmass2 fields = ",fields
        return fields


class CMASS3(PointMassElement):
    """
    Defines a scalar mass element that is connected only to scalar points.
    CMASS3 EID PID S1 S2
    """
    type = 'CMASS3'

    def __init__(self, card=None, data=None):
        PointMass.__init__(self, card, data)

        if card:
            self.eid = card.field(1)
            self.pid = card.field(2, self.eid)
            self.s1 = card.field(3)
            self.s2 = card.field(4)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.s1 = data[2]
            self.s2 = data[3]
        ###

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

    def cross_reference(self, mesh):
        """
        links up the propertiy ID
        """
        #self.s1 = mesh.Node(self.s1)
        #self.s2 = mesh.Node(self.s2)
        self.pid = mesh.Property(self.pid)

    def rawFields(self):
        fields = ['CMASS3', self.eid, self.Pid(), self.s1, self.s2]
        return fields


class CMASS4(PointMassElement):
    """
    Defines a scalar mass element that is connected only to scalar points, without reference to a property entry
    CMASS4 EID M S1 S2
    """
    type = 'CMASS4'

    def __init__(self, card=None, data=None):
        PointMass.__init__(self, card, data)

        if card:
            self.eid = card.field(1)
            self.mass = card.field(2, 0.)
            self.s1 = card.field(3)
            self.s2 = card.field(4)
        else:
            self.eid = data[0]
            self.mass = data[1]
            self.s1 = data[2]
            self.s2 = data[3]
        ###

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

    def cross_reference(self, mesh):
        """
        """
        #self.s1 = mesh.Node(self.s1)
        #self.s2 = mesh.Node(self.s2)
        pass

    def rawFields(self):
        fields = ['CMASS4', self.eid, self.mass, self.s1, self.s2]
        return fields


class CONM1(PointMass):
    type = 'CONM1'

    def __init__(self, card=None, data=None):
        """
        Concentrated Mass Element Connection, General Form
        Defines a 6 x 6 symmetric mass matrix at a geometric grid point

        CONM1 EID G CID M11 M21 M22 M31 M32
              M33 M41 M42 M43 M44 M51 M52 M53
              M54 M55 M61 M62 M63 M64 M65 M66

        [M] = [M11 M21 M31 M41 M51 M61]
              [    M22 M32 M42 M52 M62]
              [        M33 M43 M53 M63]
              [            M44 M54 M64]
              [    Sym         M55 M65]
              [                    M66]
        """
        PointMass.__init__(self, card, data)
        m = zeros((6, 6))
        if card:
            #self.nids  = [ card[1] ]
            #del self.nids
            #self.pid = None
            self.eid = card.field(1)
            self.nid = card.field(2)
            self.cid = card.field(3, 0)

            m[0, 0] = card.field(4, 0.)   # M11
            m[1, 0] = card.field(5, 0.)   # M21
            m[1, 1] = card.field(6, 0.)   # M22
            m[2, 0] = card.field(7, 0.)   # M31
            m[2, 1] = card.field(8, 0.)   # M32
            m[2, 2] = card.field(9, 0.)   # M33
            m[3, 0] = card.field(10, 0.)  # M41
            m[3, 1] = card.field(11, 0.)  # M42
            m[3, 2] = card.field(12, 0.)  # M43
            m[3, 3] = card.field(13, 0.)  # M44
            m[4, 0] = card.field(14, 0.)  # M51
            m[4, 1] = card.field(15, 0.)  # M52
            m[4, 2] = card.field(16, 0.)  # M53
            m[4, 3] = card.field(17, 0.)  # M54
            m[4, 4] = card.field(18, 0.)  # M55
            m[5, 0] = card.field(19, 0.)  # M61
            m[5, 1] = card.field(20, 0.)  # M62
            m[5, 2] = card.field(21, 0.)  # M63
            m[5, 3] = card.field(22, 0.)  # M64
            m[5, 4] = card.field(23, 0.)  # M65
            m[5, 5] = card.field(24, 0.)  # M66
        ###
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
        ###
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
        self.nid = model.Node(self.nid)
        self.cid = model.Coord(self.cid)

    def MassMatrix(self):
        return self.massMatrix

    def rawFields(self):
        cid = set_blank_if_default(self.Cid(), 0)
        nid = self.Nid()
        m = self.massMatrix
        fields = ['CONM1', self.eid, nid, cid, m[0, 0], m[1, 0], m[1, 1], m[2, 0], m[2, 1],
                  m[2, 2], m[3, 0], m[3, 1], m[3, 2], m[3, 3], m[4, 0], m[4, 1], m[4, 2],
                  m[4, 3], m[4, 4], m[5, 0], m[5, 1], m[5, 2], m[5, 3], m[5, 4], m[5, 5]]
        return fields

    def reprFields(self):
        fields = self.rawFields()
        fields2 = fields[0:4]
        for field in fields[4:]:
            val = set_blank_if_default(field, 0.)
            #print "field=%s value=%s" %(field,val)
            fields2.append(val)
        #print "fields  = ",fields
        #print "fields2 = ",fields2
        return fields2


class CONM2(PointMassElement):  ## @todo not done
    """
    @param self the CONM2 object
    @param eid element ID
    @param nid node ID
    @param cid coordinate frame of the offset (-1=absolute coordinates)
    @param X offset vector
    @param I mass moment of inertia matrix about the CG
    """
    type = 'CONM2'
    # 'CONM2    501274  11064          132.274'

    def __init__(self, card=None, data=None):
        PointMassElement.__init__(self, card, data)
        if card:
            #self.nids  = [ card[1] ]
            #del self.nids
            #self.pid = None
            self.eid = card.field(1)
            self.nid = card.field(2)
            self.cid = card.field(3, 0)
            self.mass = card.field(4, 0.)
            self.X = array(card.fields(5, 8, [0., 0., 0.]))
            self.I = card.fields(9, 15, [0.] * 6)
        else:
            self.eid = data[0]
            self.nid = data[1]
            self.cid = data[2]
            self.mass = data[3]
            self.X = data[4:7]
            self.I = data[7:]
        ###

    def Mass(self):
        return self.mass

    def Inertia(self):
        """
        Returns the 3x3 inertia matrix
        @warning doesnt handle offsets or coordinate systems
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
            raise NotImplementedError()
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
        if self.Cid() != -1:
            self.cid = model.Coord(self.cid)
        self.nid = model.Node(self.nid)

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
        msg += "              NOEUD=N%s\n" % (self.Nid())
        msg += "              VALE=%g),\n" % (self.mass)
        return msg

    def rawFields(self):
        #print "X=%r" %(self.X)
        #print "I=%r" %(self.I)
        fields = ['CONM2', self.eid, self.Nid(), self.Cid(),
                  self.mass] + list(self.X) + [None] + list(self.I)
        return fields

    def reprFields(self):
        I = []
        for i in self.I:
            if i == 0.:
                I.append(None)
            else:
                I.append(i)
            ###
        ###
        X = []
        for x in self.X:
            if x == 0.:
                X.append(None)
            else:
                X.append(x)
            ###
        ###

        cid = set_blank_if_default(self.Cid(), 0)
        mass = set_blank_if_default(self.mass, 0.)
        fields = ['CONM2', self.eid, self.Nid(), cid, mass] + X + [None] + I
        return fields
