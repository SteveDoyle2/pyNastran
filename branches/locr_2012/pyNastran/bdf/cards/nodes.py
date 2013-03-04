# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from numpy import array

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru, collapse_thru
from pyNastran.bdf.format import (integer, integer_or_blank,
                                  double, double_or_blank)


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

    def Cp(self):
        if isinstance(self.cp, int):
            return self.cp
        else:
            return self.cp.cid

    def Seid(self):
        if isinstance(self.seid, int):
            return self.seid
        else:
            return self.seid.seid

    def Cd(self):
        if isinstance(self.cd, int):
            return self.cd
        else:
            return self.cd.cid


class RINGAX(Ring):
    """
    Defines a ring for conical shell problems

    @code
    RINGAX ID R    Z    PS
    RINGAX 3  2.0 -10.0 162
    @endcode
    """
    type = 'RINGAX'

    def __init__(self, card=None, data=None, comment=''):  # this card has missing fields
        Ring.__init__(self, card, data)
        if comment:
            self._comment = comment
        self.nid = integer(card, 1, 'nid')
        self.R = double(card, 3, 'R')
        self.z = double(card, 4, 'z')
        self.ps = integer_or_blank(card, 7, 'ps')

    def Position(self):
        return array([0., 0., 0.])

    def rawFields(self):
        list_fields = ['RINGAX', self.nid, None, self.R, self.z, None,
                  None, self.ps]
        return list_fields


class SPOINT(Node):
    type = 'SPOINT'

    def __init__(self, nid, comment=''):
        Node.__init__(self, card=None, data=None)
        if comment:
            self._comment = comment
        self.nid = nid

    def cross_reference(self, model):
        pass

    def Position(self):
        return array([0., 0., 0.])

    def rawFields(self):
        if isinstance(self.nid, int):
            list_fields = ['SPOINT'] + self.nid
        else:
            list_fields = ['SPOINT'] + collapse_thru(self.nid)
        return list_fields


class SPOINTs(Node):
    """
    @code
    SPOINT ID1 ID2 ID3 ID4 ID5 ID6 ID7 ID8
    @endcode

    or

    @code
    SPOINT ID1 THRU ID2
    SPOINT 5   THRU 649
    @endcode
    """
    type = 'SPOINT'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)
        #nFields = card.nFields()

        if card:
            fields = integer(card, 1, 'ID')
        else:
            fields = data
        self.spoints = set(expand_thru(fields))
        #i = 0
        #while i<nFields: # =1 ???
        #    if fields[i]=='THRU':
        #        self.spoints += [fields[i-1],fields[i]+1]
        #        i+=1
        #    else:
        #        self.spoints.append(fields[i])
        #    i+=1

    def nDOF(self):
        return len(self.spoints)

    def addSPoints(self, sList):
        #print('old=%s new=%s' %(self.spoints,sList))
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
        #print("self.spoints = %s" %(self.spoints))
        spoints = collapse_thru(spoints)
        list_fields = ['SPOINT'] + spoints
        return list_fields

    def reprFields(self):
        return self.rawFields()


class GRDSET(Node):
    """
    Defines default options for fields 3, 7, 8, and 9 of all GRID entries.
    """
    type = 'GRDSET'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        ## Grid point coordinate system
        self.cp = integer_or_blank(card, 2, 'cp', 0)

        ## Analysis coordinate system
        self.cd = integer_or_blank(card, 6, 'cd', 0)

        ## Default SPC constraint on undefined nodes
        self.ps = str(integer_or_blank(card, 7, 'ps', ''))

        ## Superelement ID
        self.seid = integer_or_blank(card, 8, 'seid', 0)

    def cross_reference(self, model):
        self.cp = model.Coord(self.cp)
        self.cd = model.Coord(self.cd)
        #self.seid = model.SuperElement(self.seid)

    def rawFields(self):
        list_fields = ['GRDSET', None, self.Cp(), None, None, None,
                  self.Cd(), self.ps, self.Seid()]
        return list_fields

    def reprFields(self):
        cp = set_blank_if_default(self.Cp(), 0)
        cd = set_blank_if_default(self.Cd(), 0)
        ps = set_blank_if_default(self.ps, 0)
        seid = set_blank_if_default(self.Seid(), 0)
        list_fields = ['GRDSET', None, cp, None, None, None, cd, ps, seid]
        return list_fields


class GRIDB(Node):
    type = 'GRIDB'

    def __init__(self, card=None, data=None, comment=''):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            raise NotImplementedError('GRIDB data')
        else:
            print(data)
            self.nid = data[0]
            self.phi = data[1]
            self.cd = data[2]
            self.ps = data[3]
            self.idf = data[4]

        assert self.nid > 0, 'nid=%s' % (self.nid)
        assert self.phi >= 0, 'phi=%s' % (self.phi)
        assert self.cd >= 0, 'cd=%s' % (self.cd)
        assert self.ps >= 0, 'ps=%s' % (self.ps)
        assert self.idf >= 0, 'idf=%s' % (self.idf)

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


class GRID(Node):
    type = 'GRID'

    def __init__(self, card=None, data=None, comment=''):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            ## Node ID
            self.nid = integer(card, 1, 'nid')

            ## Grid point coordinate system
            self.cp = integer_or_blank(card, 2, 'cp', 0)
            
            x = double_or_blank(card, 3, 'x1', 0.)
            y = double_or_blank(card, 4, 'x2', 0.)
            z = double_or_blank(card, 5, 'x3', 0.)
            #xyz = card.fields(3, 6, [0., 0., 0.])
            ## node location in local frame
            self.xyz = array([x, y, z])

            ## Analysis coordinate system
            self.cd = integer_or_blank(card, 6, 'cd', 0)

            ## SPC constraint
            self.ps = str(integer_or_blank(card, 7, 'ps', ''))

            ## Superelement ID
            self.seid = integer_or_blank(card, 8, 'seid', 0)
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

    def nDOF(self):
        return 6

    def UpdatePosition(self, model, xyz, cid):
        self.xyz = xyz
        self.cp = model.Coord(cid)
        #assert cid == 0

    def Position(self, debug=False):
        """
        returns the point in the global XYZ coordinate system
        @param self the object pointer
        @param debug developer debug
        """
        p, matrix = self.cp.transformToGlobal(self.xyz, debug=debug)
        return p

    def PositionWRT(self, model, cid, debug=False):
        """
        returns the point which started in some arbitrary local coordinate
        system and returns it in the desired coordinate system
        @param self the object pointer
        @param model the BDF model object
        @param cid the desired coordinate ID (int)
        @param debug developer debug
        """
        if cid == self.Cp():
            return self.xyz
        #coordA = model.Coord(cid)
        # converting the xyz point arbitrary->global
        p, matrixDum = self.cp.transformToGlobal(self.xyz, debug=debug)
        #print "wrt = ",p
        coordB = model.Coord(cid)

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
        self.cp = model.Coord(self.cp)
        if self.cd != -1:
            self.cd = model.Coord(self.cd)
        #self.xyzGlobal = coord.transformToGlobal(self.xyz)

    def rawFields(self):
        list_fields = (['GRID', self.nid, self.Cp()] + list(self.xyz) +
                       [self.Cd(), self.ps, self.Seid()])
        return list_fields

    def reprFields(self):
        cp = set_blank_if_default(self.Cp(), 0)
        cd = set_blank_if_default(self.Cd(), 0)
        seid = set_blank_if_default(self.Seid(), 0)
        list_fields = ['GRID', self.nid, cp] + list(self.xyz) + [cd, self.ps,
                                                                 seid]
        return list_fields


class POINT(Node):
    type = 'POINT'

    def __init__(self, card=None, data=None, comment=''):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        if comment:
            self._comment = comment
        Node.__init__(self, card, data)

        if card:
            ## Node ID
            self.nid = integer(card, 1, 'nid')

            ## Grid point coordinate system
            self.cp = integer_or_blank(card, 2, 'cp', 0)

            x = double_or_blank(card, 3, 'x1', 0.)
            y = double_or_blank(card, 4, 'x2', 0.)
            z = double_or_blank(card, 5, 'x3', 0.)
            #xyz = card.fields(3, 6, [0., 0., 0.])
            ## node location in local frame
            self.xyz = array([x, y, z])

            ## Analysis coordinate system
            self.cd = integer_or_blank(card, 6, 'cd', 0)

            ## SPC constraint
            self.ps = str(integer_or_blank(card, 7, 'ps', ''))

            ## Superelement ID
            self.seid = integer_or_blank(card, 8, 'seid', 0)
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
        returns the point in the global XYZ coordinate system
        @param self the object pointer
        @param debug developer debug
        """
        p, matrix = self.cp.transformToGlobal(self.xyz, debug=debug)
        return p

    def PositionWRT(self, model, cid, debug=False):
        """
        returns the point which started in some arbitrary local coordinate
        system and returns it in the desired coordinate system
        @param self the object pointer
        @param model the BDF model object
        @param cid the desired coordinate ID (int)
        @param debug developer debug
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

    def cross_reference(self, model):
        self.cp = model.Coord(self.cp)

    def rawFields(self):
        list_fields = ['POINT', self.nid, self.Cp()] + list(self.xyz)
        return list_fields

    def reprFields(self):
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = ['POINT', self.nid, cp] + list(self.xyz)
        return list_fields
