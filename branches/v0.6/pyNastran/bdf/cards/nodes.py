# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import array

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru, collapse_thru
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string)
from pyNastran.bdf.fieldWriter import print_card_8

class Ring(BaseCard):
    """Generic Ring base class"""
    def __init__(self, card, data):
        assert card is None or data is None


class Node(BaseCard):
    """Generic Node base class"""
    def __init__(self, card, data):
        assert card is None or data is None

    def cross_reference(self, model):
        msg = '%s hasnt implemented a cross_reference method' % self.type
        raise NotImplementedError(msg)


class RINGAX(Ring):
    """
    Defines a ring for conical shell problems.

    +-------+-----+-----+-----+----+-----+-----+------+-----+
    |   1   |  2  |  3  |  4  |  5 |  6  |  7  |  8   |  9  |
    +=======+=====+=====+=====+====+=====+=====+======+=====+
    |RINGAX | MID |     |  R  |  Z |     |     | PS   |     |
    +-------+-----+-----+-----+----+-----+-----+------+-----+
    """
    type = 'RINGAX'
    _field_map = {1: 'mid', 3:'R', 4:'z', 7:'ps'}

    def __init__(self, card=None, data=None, comment=''):  # this card has missing fields
        Ring.__init__(self, card, data)
        if comment:
            self._comment = comment
        self.nid = integer(card, 1, 'nid')
        blank(card, 2, 'blank')
        self.R = double(card, 3, 'R')
        self.z = double(card, 4, 'z')
        blank(card, 5, 'blank')
        blank(card, 6, 'blank')
        self.ps = integer_or_blank(card, 7, 'ps')
        assert len(card) <= 8, 'len(RINGAX card) = %i' % len(card)

    def Position(self):
        return array([0., 0., 0.])

    def rawFields(self):
        list_fields = ['RINGAX', self.nid, None, self.R, self.z, None,
                  None, self.ps]
        return list_fields

    def write_bdf(self, f, method):
        card = self.reprFields()
        f.write(method(card))


class SPOINT(Node):
    type = 'SPOINT'

    def __init__(self, nid, comment=''):
        Node.__init__(self, card=None, data=None)
        if comment:
            self._comment = comment
        self.nid = nid
        assert isinstance(nid, int), nid

    def cross_reference(self, model):
        pass

    def Position(self):
        return array([0., 0., 0.])

    def rawFields(self):
        if isinstance(self.nid, int):
            list_fields = ['SPOINT', self.nid]
        else:
            list_fields = ['SPOINT'] + collapse_thru(self.nid)
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + print_card_8(card)


class SPOINTs(Node):
    """
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |   1    |  2  |  3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +========+=====+======+=====+=====+=====+=====+=====+=====+
    | SPOINT | ID1 | THRU | ID2 |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    | SPOINT | ID1 | ID1  | ID3 | ID4 | ID5 | ID6 | ID7 | ID8 |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    |        | ID8 | etc. |     |     |     |     |     |     |
    +--------+-----+------+-----+-----+-----+-----+-----+-----+
    """
    type = 'SPOINT'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            fields = []
            for i in range(1, len(card)):
                field = integer_or_string(card, i, 'ID%i' % i)
                fields.append(field)
        else:
            fields = data
            assert isinstance(data, list), data
            assert isinstance(data[0], int), data
        self.spoints = set(expand_thru(fields))

    def nDOF(self):
        return len(self.spoints)

    def addSPoints(self, sList):
        #print('old=%s new=%s' % (self.spoints, sList))
        self.spoints = self.spoints.union(set(sList))

    def cross_reference(self, model):
        pass

    def createSPOINTi(self):
        spoints = []
        for nid in self.spoints:
            spoints.append(SPOINT(nid))
        return spoints

    def rawFields(self):
        #print("SPOINTi")
        spoints = list(self.spoints)
        spoints.sort()
        #print("self.spoints = %s" % self.spoints)
        spoints = collapse_thru(spoints)
        list_fields = ['SPOINT'] + spoints
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + print_card_8(card)


class GRDSET(Node):
    """
    Defines default options for fields 3, 7, 8, and 9 of all GRID entries.

    +--------+-----+----+----+----+----+----+----+------+
    |    1   |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
    +========+=====+====+====+====+====+====+====+======+
    | GRDSET |     | CP |    |    |    | CD | PS | SEID |
    +--------+-----+----+----+----+----+----+----+------+
    """
    type = 'GRDSET'
    _field_map = {1: 'nid', 2:'cp', 6:'cd', 7:'ps', 8:'seid'}

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        #: Grid point coordinate system
        blank(card, 1, 'blank')
        self.cp = integer_or_blank(card, 2, 'cp', 0)
        blank(card, 3, 'blank')
        blank(card, 4, 'blank')
        blank(card, 5, 'blank')

        #: Analysis coordinate system
        self.cd = integer_or_blank(card, 6, 'cd', 0)

        #: Default SPC constraint on undefined nodes
        self.ps = str(integer_or_blank(card, 7, 'ps', ''))

        #: Superelement ID
        self.seid = integer_or_blank(card, 8, 'seid', 0)
        assert len(card) <= 9, 'len(GRDSET card) = %i' % len(card)

    def cross_reference(self, model):
        msg = ' which is required by the GRDSET'
        self.cp = model.Coord(self.cp, msg=msg)
        self.cd = model.Coord(self.cd, msg=msg)
        #self.seid = model.SuperElement(self.seid, msg)

    def Cd(self):
        if isinstance(self.cd, int):
            return self.cd
        else:
            return self.cd.cid

    def Cp(self):
        if isinstance(self.cp, int):
            return self.cp
        else:
            return self.cp.cid

    def SEid(self):
        if isinstance(self.seid, int):
            return self.seid
        else:
            return self.seid.seid

    def _verify(self, xref=False):
        cp = self.Cp()
        seid = self.SEid()
        cd = self.Cd()
        ps = self.Ps()
        assert isinstance(cp, int), 'cp=%r' % cp
        assert isinstance(cd, int), 'cd=%r' % cd
        assert isinstance(seid, int), 'seid=%r' % seid

    def rawFields(self):
        list_fields = ['GRDSET', None, self.Cp(), None, None, None,
                  self.Cd(), self.ps, self.SEid()]
        return list_fields

    def reprFields(self):
        cp = set_blank_if_default(self.Cp(), 0)
        cd = set_blank_if_default(self.Cd(), 0)
        ps = set_blank_if_default(self.ps, 0)
        seid = set_blank_if_default(self.SEid(), 0)
        list_fields = ['GRDSET', None, cp, None, None, None, cd, ps, seid]
        return list_fields

    def write_bdf(self, f, method):
        card = self.reprFields()
        f.write(print_card_8(card))


class GRIDB(Node):
    type = 'GRIDB'
    _field_map = {1: 'nid', 4:'phi', 6:'cd', 7:'ps', 8:'idf'}

    def __init__(self, card=None, data=None, comment=''):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            self.nid = integer(card, 1, 'nid')
            self.phi = double(card, 4, 'phi')
            self.cd = integer(card, 6, 'cd')
            self.ps = integer(card, 7, 'ps')
            self.idf = integer(card, 8, 'idf')
        else:
            self.nid = data[0]
            self.phi = data[1]
            self.cd = data[2]
            self.ps = data[3]
            self.idf = data[4]

        assert self.nid > 0, 'nid=%s' % self.nid
        assert self.phi >= 0, 'phi=%s' % self.phi
        assert self.cd >= 0, 'cd=%s' % self.cd
        assert self.ps >= 0, 'ps=%s' % self.ps
        assert self.idf >= 0, 'idf=%s' % self.idf

    def _verify(self, xref=False):
        pass

    def Cd(self):
        if isinstance(self.cd, int):
            return self.cd
        else:
            return self.cd.cid

    def rawFields(self):
        list_fields = ['GRIDB', self.nid, None, None, self.phi, None,
                  self.Cd(), self.ps, self.idf]
        return list_fields

    def reprFields(self):
        #phi = set_blank_if_default(self.phi,0.0)
        cd = set_blank_if_default(self.Cd(), 0)
        ps = set_blank_if_default(self.ps, 0)
        idf = set_blank_if_default(self.idf, 0)
        list_fields = ['GRIDB', self.nid, None, None, self.phi, None, cd, ps,
                       idf]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


class GRID(Node):
    """
    +------+-----+----+----+----+----+----+----+------+
    |   1  |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9   |
    +======+=====+====+====+====+====+====+====+======+
    | GRID | NID | CP | X1 | X2 | X3 | CD | PS | SEID |
    +------+-----+----+----+----+----+----+----+------+
    """
    type = 'GRID'
    _field_map = {1: 'nid', 2:'cp', 6:'cd', 7:'ps', 8:'seid'}

    def _update_field_helper(self, n, value):
        if n == 3:
            self.xyz[0] = value
        elif n == 4:
            self.xyz[1] = value
        elif n == 5:
            self.xyz[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            #: Node ID
            self.nid = integer(card, 1, 'nid')

            #: Grid point coordinate system
            self.cp = integer_or_blank(card, 2, 'cp', 0)

            x = double_or_blank(card, 3, 'x1', 0.)
            y = double_or_blank(card, 4, 'x2', 0.)
            z = double_or_blank(card, 5, 'x3', 0.)
            #: node location in local frame
            self.xyz = array([x, y, z])

            #: Analysis coordinate system
            self.cd = integer_or_blank(card, 6, 'cd', 0)

            #: SPC constraint
            self.ps = str(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            self.seid = integer_or_blank(card, 8, 'seid', 0)
            assert len(card) <= 9, 'len(GRID card) = %i' % len(card)
        else:
            self.nid = data[0]
            self.cp = data[1]
            self.xyz = array(data[2:5])
            self.cd = data[5]
            self.ps = data[6]
            self.seid = data[7]
            if self.ps == 0:
                self.ps = ''
            assert len(self.xyz) == 3

        assert self.nid > 0, 'nid=%s' % (self.nid)
        assert self.cp >= 0, 'cp=%s' % (self.cp)
        assert self.cd >= -1, 'cd=%s' % (self.cd)
        assert self.seid >= 0, 'seid=%s' % (self.seid)

    def Nid(self):
        return self.nid

    def Ps(self):
        return self.ps

    def Cd(self):
        if isinstance(self.cd, int):
            return self.cd
        else:
            return self.cd.cid

    def Cp(self):
        if isinstance(self.cp, int):
            return self.cp
        else:
            return self.cp.cid

    def SEid(self):
        if isinstance(self.seid, int):
            return self.seid
        else:
            return self.seid.seid

    #def SEid(self):
        #return self.seid

    def _verify(self, xref):
        nid = self.Nid()
        cp = self.Cp()
        cd = self.Cd()
        xyz = self.xyz
        ps = self.Ps()
        seid = self.SEid()
        assert isinstance(nid, int), 'nid=%r' % nid
        assert isinstance(cp, int), 'cp=%r' % cp
        assert isinstance(cd, int), 'cd=%r' % cd
        if xref:
            pos_xyz = self.Position()

    def nDOF(self):
        return 6

    def UpdatePosition(self, model, xyz, cid):
        self.xyz = xyz
        msg = ' which is required by GRID nid=%s' % self.nid
        self.cp = model.Coord(cid, msg=msg)
        #assert cid == 0

    def Position(self, debug=False):
        """
        Gets the point in the global XYZ coordinate system.

        :param self:  the object pointer
        :param debug: developer debug
        """
        p, matrix = self.cp.transformToGlobal(self.xyz, debug=debug)
        return p

    def PositionWRT(self, model, cid, debug=False):
        """
        Gets the point which started in some arbitrary local coordinate
        system and returns it in the desired coordinate system

        :param self:  the object pointer
        :param model: the BDF model object
        :param cid:   the desired coordinate ID (int)
        :param debug: developer debug
        """
        if cid == self.Cp():
            return self.xyz
        #coordA = model.Coord(cid, msg=msg)
        # converting the xyz point arbitrary->global
        p, matrixDum = self.cp.transformToGlobal(self.xyz, debug=debug)
        #print "wrt = ",p
        msg = ' which is required by %s nid=%s' % (self.type, self.nid)
        coordB = model.Coord(cid, msg=msg)

        # a matrix global->local matrix is found
        pdum, matrix = coordB.transformToGlobal(
            array([1., 0., 0]), debug=debug)
        p2 = coordB.transformToLocal(p, matrix, debug=debug)
        return p2

    def cross_reference(self, model, grdset=None):
        """
        the gridset object will only update the fields that have not been set
        """
        #print str(self)
        if grdset:  # update using a gridset object
            if not self.cp:
                self.cp = grdset.cp
            if not self.cd:
                self.cd = grdset.cd
            if not self.ps:
                self.ps = grdset.ps
            if not self.seid:
                self.seid = grdset.seid
        msg = ' which is required by %s nid=%s' % (self.type, self.nid)
        self.cp = model.Coord(self.cp, msg=msg)
        if self.cd != -1:
            self.cd = model.Coord(self.cd, msg=msg)
        #self.xyzGlobal = coord.transformToGlobal(self.xyz)

    def rawFields(self):
        list_fields = (['GRID', self.nid, self.Cp()] + list(self.xyz) +
                       [self.Cd(), self.ps, self.SEid()])
        return list_fields

    def reprFields(self):
        cp = set_blank_if_default(self.Cp(), 0)
        cd = set_blank_if_default(self.Cd(), 0)
        seid = set_blank_if_default(self.SEid(), 0)
        list_fields = ['GRID', self.nid, cp] + list(self.xyz) + [cd, self.ps,
                                                                 seid]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


class POINT(Node):
    """
    +-------+-----+----+----+----+----+----+----+-----+
    |   1   |  2  | 3  | 4  | 5  | 6  |  7 | 8  |  9  |
    +=======+=====+====+====+====+====+====+====+=====+
    | POINT | NID | CP | X1 | X2 | X3 |    |    |     |
    +-------+-----+----+----+----+----+----+----+-----+
    """
    type = 'POINT'
    _field_map = {1: 'nid', 2:'cp'}

    def __init__(self, card=None, data=None, comment=''):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            #: Node ID
            self.nid = integer(card, 1, 'nid')

            #: Grid point coordinate system
            self.cp = integer_or_blank(card, 2, 'cp', 0)

            x = double_or_blank(card, 3, 'x1', 0.)
            y = double_or_blank(card, 4, 'x2', 0.)
            z = double_or_blank(card, 5, 'x3', 0.)

            #: node location in local frame
            self.xyz = array([x, y, z])

            #: Analysis coordinate system
            self.cd = integer_or_blank(card, 6, 'cd', 0)

            #: SPC constraint
            self.ps = str(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            self.seid = integer_or_blank(card, 8, 'seid', 0)
            assert len(card) <= 9, 'len(POINT card) = %i' % len(card)
        else:
            #print data
            self.nid = data[0]
            self.cp = data[1]
            self.xyz = array(data[2:5])
            assert len(self.xyz) == 3

        assert self.nid > 0, 'nid=%s' % (self.nid)
        assert self.cp >= 0, 'cp=%s' % (self.cp)

    def nDOF(self):
        return 6

    def UpdatePosition(self, model, xyz, cid):
        self.xyz = xyz
        self.cp = model.Coord(cid)
        #assert cid == 0

    def Position(self, debug=False):
        """
        Gets the point in the global XYZ coordinate system.

        :param self:  the object pointer
        :param debug: developer debug
        """
        p, matrix = self.cp.transformToGlobal(self.xyz, debug=debug)
        return p

    def PositionWRT(self, model, cid, debug=False):
        """
        Gets the location of the GRID which started in some arbitrary
        local coordinate system and returns it in the desired coordinate
        system.

        :param self:  the object pointer
        :param model: the BDF model object
        :param cid:   the desired coordinate ID (int)
        :param debug: developer debug (default=True)

        :returns position: the position of the GRID in an arbitrary
                           coordinate system
        """
        if cid == self.Cp():
            return self.xyz
        #coordA = model.Coord(cid)
        # converting the xyz point arbitrary->global
        p, matrixDum = self.cp.transformToGlobal(self.xyz, debug=debug)
        coordB = model.Coord(cid)

        # a matrix global->local matrix is found
        pdum, matrix = coordB.transformToGlobal(array([1., 0., 0]),
                                                debug=debug)
        p2 = coordB.transformToLocal(p, matrix, debug=debug)
        return p2

    def Cp(self):
        if isinstance(self.cp, int):
            return self.cp
        else:
            return self.cp.cid

    def cross_reference(self, model):
        self.cp = model.Coord(self.cp)

    def rawFields(self):
        list_fields = ['POINT', self.nid, self.Cp()] + list(self.xyz)
        return list_fields

    def reprFields(self):
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = ['POINT', self.nid, cp] + list(self.xyz)
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)
