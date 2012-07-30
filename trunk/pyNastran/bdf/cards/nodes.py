# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from numpy import array

from pyNastran.bdf.fieldWriter import setBlankIfDefault
from pyNastran.bdf.cards.baseCard import BaseCard, expandThru, collapseThru

class Ring(BaseCard):
    """Generic Ring base class"""
    def __init__(self, card, data):
        assert card is None or data is None

class Node(BaseCard):
    """Generic Node base class"""
    def __init__(self, card, data):
        assert card is None or data is None

    def crossReference(self, model):
        msg = '%s hasnt implemented a crossReference method' %(self.type)
        raise NotImplementedError(msg)

    def Cp(self):
        if isinstance(self.cp, int):
            return self.cp
        else:
            return self.cp.cid
        ###

    def Seid(self):
        if isinstance(self.seid, int):
            return self.seid
        else:
            return self.seid.seid
        ###
        
    def Cd(self):
        if isinstance(self.cd, int):
            return self.cd
        else:
            return self.cd.cid
        ###

class RINGAX(Ring):
    """
    Defines a ring for conical shell problems
    RINGAX ID R    Z    PS
    RINGAX 3  2.0 -10.0 162
    """
    type = 'RINGAX'
    def __init__(self, card=None, data=None): # this card has missing fields
        Ring.__init__(self, card, data)
        self.nid = card.field(1)
        self.R   = card.field(3)
        self.z   = card.field(4)
        self.ps  = card.field(7)

    def Position(self):
        return array([0., 0., 0.])

    def rawFields(self):
        fields = ['RINGAX', self.nid, None, self.R, self.z, None, None, self.ps]
        return fields
        
    
class SPOINT(Node):
    type = 'SPOINT'
    def __init__(self, nid):
        Node.__init__(self, card=None, data=None)
        self.nid = nid

    def crossReference(self, model):
        pass

    def Position(self):
        return array([0., 0., 0.])

    def rawFields(self):
        """
        @todo support THRU in output
        """
        #print("SPOINT")
        if isinstance(self.nid, int):
            fields = ['SPOINT']+[self.nid]
        else:
            #print "self.nid = ",self.nid
            fields = ['SPOINT']+self.nid
        return fields

class SPOINTs(Node):
    """
    SPOINT ID1 ID2 ID3 ID4 ID5 ID6 ID7 ID8
    or
    SPOINT ID1 THRU ID2
    SPOINT 5   THRU 649
    """
    type = 'SPOINT'
    def __init__(self, card=None, data=None):
        Node.__init__(self, card, data)
        #nFields = card.nFields()
        
        if card:
            fields  = card.fields(1)
        else:
            fields = data
        self.spoints = set(expandThru(fields))
        #i = 0
        #while i<nFields: # =1 ???
        #    if fields[i]=='THRU':
        #        self.spoints += [fields[i-1],fields[i]+1]
        #        i+=1
        #    else:
        #        self.spoints.append(fields[i])
        #    i+=1
        ###

    def nDOF(self):
        return len(self.spoints)

    def addSPoints(self, sList):
        #print('old=%s new=%s' %(self.spoints,sList))
        self.spoints = self.spoints.union(set(sList))
        
    def crossReference(self, model):
        pass

    def createSPOINTi(self):
        spoints = []
        for nid in self.spoints:
            spoints.append(SPOINT(nid))
        ###
        return spoints

    def rawFields(self):
        #print("SPOINTi")
        spoints = list(self.spoints)
        spoints.sort()
        #print("self.spoints = %s" %(self.spoints))
        spoints = collapseThru(spoints)
        fields = ['SPOINT']+spoints
        return fields

    def reprFields(self):
        return self.rawFields()

class GRDSET(Node):
    """
    Defines default options for fields 3, 7, 8, and 9 of all GRID entries.
    """
    type = 'GRDSET'
    def __init__(self, card=None, data=None):
        ## Grid point coordinate system
        self.cp  = card.field(2, 0)
        
        ## Analysis coordinate system
        self.cd   = card.field(6, 0)
        
        ## Default SPC constraint on undefined nodes
        self.ps   = str(card.field(7, ''))
        
        ## Superelement ID
        self.seid = card.field(8, 0)

    def crossReference(self, model):
        self.cp = model.Coord(self.cp)
        self.cd = model.Coord(self.cd)
        #self.seid = model.SuperElement(self.seid)

    def rawFields(self):
        fields = ['GRDSET', None, self.Cp(), None, None, None, self.Cd(), self.ps, self.Seid()]
        return fields

    def reprFields(self):
        cp   = setBlankIfDefault(self.Cp(), 0)
        cd   = setBlankIfDefault(self.Cd(), 0)
        ps   = setBlankIfDefault(self.ps, 0)
        seid = setBlankIfDefault(self.Seid(), 0)
        fields = ['GRDSET', None, cp, None, None, None, cd, ps, seid]
        return fields

class GRIDB(Node):
    type = 'GRIDB'
    def __init__(self, card=None, data=None):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        Node.__init__(self, card, data)

        if card:
            raise Exception('not implemented...')
        else:
            print(data)
            self.nid = data[0]
            self.phi = data[1]
            self.cd  = data[2]
            self.ps  = data[3]
            self.idf = data[4]
        ###
        assert self.nid > 0,  'nid=%s'  %(self.nid)
        assert self.phi >= 0, 'phi=%s'  %(self.phi)
        assert self.cd  >= 0, 'cd=%s'   %(self.cd)
        assert self.ps  >= 0, 'ps=%s'   %(self.ps)
        assert self.idf >= 0, 'idf=%s'  %(self.idf)

    def rawFields(self):
        fields = ['GRIDB', self.nid, None, None, self.phi, None, self.Cd(), self.ps, self.idf]
        return fields

    def reprFields(self):
        #phi = setBlankIfDefault(self.phi,0.0)
        cd  = setBlankIfDefault(self.Cd(), 0)
        ps  = setBlankIfDefault(self.ps, 0)
        idf = setBlankIfDefault(self.idf, 0)
        fields = ['GRIDB', self.nid, None, None, self.phi, None, cd, ps, idf]
        return fields

class GRID(Node):
    type = 'GRID'
    def __init__(self, card=None, data=None):
        """
        if coming from a BDF object, card is used
        if coming from the OP2, data is used
        """
        Node.__init__(self, card, data)

        if card:
            ## Node ID
            self.nid = int(card.field(1))

            ## Grid point coordinate system
            self.cp = card.field(2, 0)

            xyz = card.fields(3, 6, [0., 0., 0.])
            #displayCard(card)
            #print "xyz = ",xyz
            self.xyz = array(xyz)

            ## Analysis coordinate system
            self.cd = card.field(6, 0)

            ## SPC constraint
            self.ps = str(card.field(7, ''))

            ## Superelement ID
            self.seid = card.field(8, 0)

            #print "xyz = ",self.xyz
            #print "cd = ",self.cd
            #print "ps = ",self.ps
        else:
            #print data
            self.nid  = data[0]
            self.cp   = data[1]
            self.xyz  = array(data[2:5])
            self.cd   = data[5]
            self.ps   = data[6]
            self.seid = data[7]
            if self.ps==0:
                self.ps=''
            assert len(self.xyz)==3
        ###
        assert self.nid  >   0,  'nid=%s' %(self.nid)
        assert self.cp   >=  0,  'cp=%s'  %(self.cp)
        assert self.cd   >= -1,  'cd=%s'  %(self.cd)
        assert self.seid >=  0,'seid=%s'  %(self.seid)

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
        p,matrix = self.cp.transformToGlobal(self.xyz, debug=debug)
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
        if cid==self.Cp():
            return self.xyz
        #coordA = model.Coord(cid)
        # converting the xyz point arbitrary->global
        p,matrixDum = self.cp.transformToGlobal(self.xyz, debug=debug)
        #print "wrt = ",p
        coordB = model.Coord(cid)
        
        # a matrix global->local matrix is found
        pdum,matrix = coordB.transformToGlobal(array([1.,0.,0]), debug=debug)
        p2          = coordB.transformToLocal(p, matrix, debug=debug)
        return p2

    def crossReference(self, model, grdset=None):
        """
        the gridset object will only update the fields that have not been set
        """
        #print str(self)
        if grdset: # update using a gridset object
            if not self.cp:   self.cp   = grdset.cp
            if not self.cd:   self.cd   = grdset.cd
            if not self.ps:   self.ps   = grdset.ps
            if not self.seid: self.seid = grdset.seid
        self.cp = model.Coord(self.cp)
        if self.cd != -1:
            self.cd = model.Coord(self.cd)
        #self.xyzGlobal = coord.transformToGlobal(self.xyz)

    def rawFields(self):
        fields = ['GRID', self.nid, self.Cp()]+list(self.xyz)+[self.Cd(), self.ps, self.Seid()]
        return fields

    def reprFields(self):
        cp   = setBlankIfDefault(self.Cp(), 0)
        cd   = setBlankIfDefault(self.Cd(), 0)
        seid = setBlankIfDefault(self.Seid(), 0)
        fields = ['GRID', self.nid, cp]+list(self.xyz)+[cd, self.ps, seid]
        return fields

