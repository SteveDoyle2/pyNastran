# pylint: disable=R0904,R0902
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from itertools import izip

from pyNastran.bdf.fieldWriter import (printCard, set_default_if_blank,
                                       is_same)
from pyNastran.bdf.bdfInterface.BDF_Card import BDFCard


class BaseCard(BDFCard):
    #def __init__(self,card):
    #    pass

    def writeCodeAster(self):
        return ('# skipping %s  because writeCodeAster is not implemented\n'
                % self.type)

    def writeCodeAsterLoad(self, model, gridWord='node'):
        return '# skipping %s (lid=%s) because writeCodeAsterLoad is not implemented\n' % (self.type, self.lid)

    #def verify(self, model, iSubcase):
        #"""
        #this method checks performs checks on the cards such as
        #that the PBEAML has a proper material type
        #"""
        #pass

    def isSameFields(self, fields1, fields2):
        for (field1, field2) in izip(fields1, fields2):
            if not is_same(field1, field2):
                return False
        return True

    #def Is(self, typeCheck):
        #"""retruns True if the card type is the same as the object"""
        #if self.type == typeCheck:
        #    return True
        #return False

    #def removeTrailingNones(self, fields):
        #"""removes blank fields at the end of a card object"""
        #self._wipeEmptyFields(fields)

    def printCard(self, fields, tol=0.):
        """prints a card object"""
        return printCard(fields, tol)

    def cross_reference(self, model):
        #self.mid = model.Material(self.mid)
        msg = "%s needs to implement the 'cross_reference' method" % (self.type)
        raise NotImplementedError(msg)

    def buildTableLines(self, fields, nStart=1, nEnd=0):
        """
        builds a table of the form:
        @code
        'DESVAR' DVID1 DVID2 DVID3 DVID4 DVID5 DVID6 DVID7
                 DVID8 -etc.-
        'UM'     VAL1  VAL2  -etc.-
        @endcode
        
        and then pads the rest of the fields with None's
        @param fields the fields to enter, including DESVAR
        @param nStart the number of blank fields at the start of the line (default=1)
        @param nEnd the number of blank fields at the end of the line (default=0)

        @note will be used for DVPREL2, RBE1, RBE3
        @warning only works for small field format???
        """
        fieldsOut = []
        n = 8 - nStart - nEnd

        # pack all the fields into a list.  Only consider an entry as isolated
        for (i, field) in enumerate(fields):
            fieldsOut.append(field)
            if i > 0 and i % n == 0:  # beginning of line
                #print "i=%s" %(i)
                #pad = [None]*(i+j)
                #fieldsOut += pad
                fieldsOut += [None] * (nStart + nEnd)

        # make sure they're aren't any trailing None's (from a new line)
        fieldsOut = self._wipeEmptyFields(fieldsOut)
        #print "fieldsOut = ",fieldsOut,len(fieldsOut)

        # push the next key (aka next fields[0]) onto the next line
        nSpaces = 8 - (len(fieldsOut)) % 8  # puts UM onto next line
        #print "nSpaces[%s] = %s max=%s" %(fields[0],nSpaces,nSpaceMax)
        if nSpaces < 8:
            fieldsOut += [None] * (nSpaces)
        #print ""
        return fieldsOut

    def isSameCard(self, card, debug=False):
        fields1 = self.rawFields()
        fields2 = card.rawFields()
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def printRawFields(self):
        """A card's raw fields include all defaults for all fields"""
        fields = self.rawFields()
        return self.printCard(fields)

    def reprFields(self):
        return self.rawFields()

    def __repr__(self):  # ,tol=1e-8
        """
        Prints a card in the simplest way possible
        (default values are left blank).
        """
        #print "tol = ",tol
        #self.rawFields()
        fields = self.reprFields()
        try:
            return self.printCard(fields)
        except:
            print('problem printing %s card' % (self.type))
            print("fields = ", fields)
            raise


def Mid(self):
    #print str(self)
    if isinstance(self.mid, int):
        return self.mid
    else:
        return self.mid.mid
    ###


class Property(BaseCard):
    def __init__(self, card, data):
        assert card is None or data is None
        pass

    def Mid(self):
        return Mid(self)

    def isSameCard(self, prop, debug=False):
        if self.type != prop.type:
            return False
        fields1 = self.rawFields()
        fields2 = prop.rawFields()
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)


class Material(BaseCard):
    """Base Material Class"""
    def __init__(self, card, data):
        pass
        #self.type = card[0]

    def isSameCard(self, mat, debug=False):
        if self.type != mat.type:
            return False
        fields1 = self.rawFields()
        fields2 = mat.rawFields()
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def cross_reference(self, model):
        pass

    def Mid(self):
        return self.mid


class Element(BaseCard):
    pid = 0  # CONM2, rigid

    def __init__(self, card, data):
        assert card is None or data is None
        ## the list of node IDs for an element (default=None)
        self.nodes = None
        #self.nids = []
        pass

    def isSameCard(self, element, debug=False):
        if self.type != element.type:
            return False
        fields1 = self.rawFields()
        fields2 = element.rawFields()
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def Pid(self):
        """returns the property ID of an element"""
        if isinstance(self.pid, int):
            return self.pid
        else:
            return self.pid.pid

    def nodePositions(self, nodes=None):
        """returns the positions of multiple node objects"""
        if not nodes:
            nodes = self.nodes
        return [node.Position() for node in nodes]

    def nodeIDs(self, nodes=None, allowEmptyNodes=False, msg=''):
        """returns nodeIDs for repr functions"""
        try:
            if not nodes:
                nodes = self.nodes

            if allowEmptyNodes:
                nodes2 = []
                for i, node in enumerate(nodes):
                    if node == 0 or node is None:
                        nodes2.append(None)
                    elif isinstance(node, int):
                        nodes2.append(node)
                    else:
                        nodes2.append(node.nid)
                return nodes2
            else:
                if isinstance(nodes[0], int):
                    nodeIDs = [node for node in nodes]
                else:
                    nodeIDs = [node.nid for node in nodes]
                ###
                assert 0 not in nodeIDs, 'nodeIDs = %s' % (nodeIDs)
                return nodeIDs
        except:
            print("nodes=%s allowEmptyNodes=%s\nmsg=%s" % (
                nodes, allowEmptyNodes, msg))
            raise

    def prepareNodeIDs(self, nids, allowEmptyNodes=False):
        """verifies all node IDs exist and that they're integers"""
        self.nodes = []
        for nid in nids:
            if isinstance(nid, int):
                self.nodes.append(int(nid))
            elif nid is None and allowEmptyNodes:
                self.nodes.append(nid)
            else:  # string???
                self.nodes.append(int(nid))
                #raise RuntimeError('this element may not have missing nodes...nids=%s allowEmptyNodes=False' %(nids))
            ###

    #def Normal(self,a,b):
    #    """finds the unit normal vector of 2 vectors"""
    #    return Normal(a,b)

    def CentroidTriangle(self, n1, n2, n3, debug=False):
        if debug:
            print("n1=%s \nn2=%s \nn3=%s" % (n1, n2, n3))
        centroid = (n1 + n2 + n3) / 3.
        return centroid

    def Centroid(self):
        msg = 'Centroid not implemented in the %s class' % (
            self.__class__.__name__)
        raise NotImplementedError(msg)

    def Length(self):
        msg = 'Length not implemented in the %s class' % (
            self.__class__.__name__)
        raise NotImplementedError(msg)

    def Area(self):
        msg = 'Area not implemented in the %s class' % (
            self.__class__.__name__)
        raise NotImplementedError(msg)

    def Volume(self):
        msg = 'Volume not implemented in the %s class' % (
            self.__class__.__name__)
        raise NotImplementedError(msg)

    def Mass(self):
        msg = 'Mass not implemented in the %s class' % (
            self.__class__.__name__)
        raise NotImplementedError(msg)

    def B(self):
        msg = 'B matrix not implemented in the %s class' % (
            self.__class__.__name__)
        raise NotImplementedError(msg)

    def D(self):
        msg = 'D matrix not implemented in the %s class' % (
            self.__class__.__name__)
        raise NotImplementedError(msg)

    def Jacobian(self):
        msg = 'Jacobian not implemented for %s' % (
            self.__class__.__name__)
        raise NotImplementedError(msg)

    def stiffnessMatrix(self):
        msg = 'stiffnessMatrix not implemented in the %s class' % (
            self.__class__.__name__)
        raise NotImplementedError(msg)

    def massMatrix(self):
        msg = 'massMatrix not implemented in the %s class' % (
            self.__class__.__name__)
        raise NotImplementedError(msg)


def expand_thru(fields):
    """
    expands a list of values of the form [1,5,THRU,9,13]
    to be [1,5,6,7,8,9,13]
    """
    if len(fields) == 1:
        return fields
    out = []
    nFields = len(fields)
    i = 0
    while(i < nFields):
        if fields[i] == 'THRU':
            for j in xrange(fields[i - 1], fields[i + 1] + 1):
                out.append(j)
            i += 2
        else:
            out.append(fields[i])
            i += 1
    return list(set(out))


def expand_thru_by(fields):
    """
    expands a list of values of the form [1,5,THRU,9,BY,2,13]
    to be [1,5,7,9,13]
    @todo not tested
    @note used for QBDY3, ???
    """
    if len(fields) == 1:
        return fields
    out = []
    nFields = len(fields)
    i = 0
    by = 1
    while(i < nFields):
        if fields[i] == 'THRU':
            by = 1
            byCase = False
            if i + 2 < nFields and fields[i + 2] == 'BY':
                by = fields[i + 3]
            else:
                by = 1
                byCase = True
            minValue = fields[i - 1]
            maxValue = fields[i + 1]
            maxR = int((maxValue - minValue) // by + 1) # max range value

            for j in xrange(0, maxR):  # +1 is to include final point
                value = minValue + by * j
                out.append(value)

            if byCase:  # null/standard case
                i += 2
            else:     # BY case
                i += 3
        else:
            out.append(fields[i])
            i += 1
    return list(set(out))


def expand_thru_exclude(self, fields):
    """
    expands a list of values of the form [1,5,THRU,11,EXCEPT,7,8,13]
    to be [1,5,6,9,10,11,13]
    @warning hasnt been tested
    """
    fieldsOut = []
    nFields = len(fields)
    for i in xrange(nFields):
        if fields[i] == 'THRU':
            storedList = []
            for j in xrange(fields[i - 1], fields[i + 1]):
                storedList.append(fields[j])

        elif fields[i] == 'EXCLUDE':
            storedSet = set(storedList)
            while fields[i] < max(storedList):
                storedSet.remove(fields[i])
            storedList = list(storedSet)
        else:
            if storedList:
                fieldsOut += storedList
            fieldsOut.append(fields[i])
    return fieldsOut


def collapse_thru_by(fields):
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    fields2 = build_thru(packs)
    #assert fields == expand_thru_by(fields2)  # why doesn't this work?
    return fields2

def collapse_thru_by_float(fields):
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    fields2 = build_thru_float(packs)
    #assert fields == expand_thru_by(fields2)  # why doesn't this work?
    return fields2


def collapse_thru(fields):
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    fields2 = build_thru(packs, maxDV=1)
    #assert fields == expand_thru_by(fields2), fields2  # why doesn't this work?
    return fields2


def condense(valueList):
    """
    Builds a list of packs (list of 3 values representing the first, last,
    and delta values for condensing a SET card.
    @see build_thru
    """
    if len(valueList) == 1:
        return [[valueList[0], valueList[0], 1]]
    valueList.sort()
    packs = []
    
    dvOld = None
    firstVal = valueList[0]
    lastVal = firstVal

    for val in valueList[1:]:
        dv = val - lastVal
        
        # sets up the first item of the pack
        if dvOld is None:
            dvOld = dv
        
        # fill up the pack
        if dvOld == dv:
            lastVal = val
        else:
            packs.append([firstVal,lastVal,dvOld])
            lastVal = val
            dvOld = None
            firstVal = val

    # fills the last pack
    if dvOld == dv:
        packs.append([firstVal,val,dv])
    else:
        packs.append([firstVal, val, dvOld])
    return packs


def build_thru(packs, maxDV=None):
    """
    Takes a pack [1,7,2] and converts it into fields used by a SET card.
    The values correspond to the first value, last value, and delta in the
    list.  This means that [1,1001,2] represents 500 values.
    [1,1001,1] represents 1001 values and will be written as [1,THRU,1001]..

    @param packs
      list of packs (list of 3 values: [first, last, delta] )
    @param maxDV
      integer defining the max allowable delta between two values
      (default=None; no limit)
    """
    fields = []
    for (firstVal, lastVal, dv) in packs:
        if firstVal == lastVal:
            fields.append(firstVal)
        elif dv == 1:
            if lastVal - firstVal > 2:
                fields.append(firstVal)
                fields.append('THRU')
                fields.append(lastVal)
            elif lastVal - firstVal == 2:
                fields.append(firstVal)
                fields.append(firstVal + 1)
                fields.append(lastVal)
            else:
                fields.append(firstVal)
                fields.append(lastVal)
        else:
            if maxDV is None:
                if lastVal - firstVal > 4*dv:
                    fields.append(firstVal)
                    fields.append('THRU')
                    fields.append(lastVal)
                    fields.append('BY')
                    fields.append(dv)
                else:
                    for v in xrange(firstVal, lastVal + dv, dv):
                        fields.append(v)
            else:
                for v in xrange(firstVal, lastVal + dv, dv):
                    fields.append(v)
    return fields

def build_thru_float(packs, maxDV=None):
    """
    Takes a pack [1,7,2] and converts it into fields used by a SET card.
    The values correspond to the first value, last value, and delta in the
    list.  This means that [1,1001,2] represents 500 values.
    [1,1001,1] represents 1001 values and will be written as [1,THRU,1001]..

    @param packs
      list of packs (list of 3 values: [first, last, delta] )
    @param maxDV
      integer defining the max allowable delta between two values
      (default=None; no limit)
    """
    fields = []
    for (firstVal, lastVal, dv) in packs:
        if lastVal - firstVal > 4*dv:
            fields.append(firstVal)
            fields.append('THRU')
            fields.append(lastVal)
            fields.append('BY')
            fields.append(dv)
        else:
            nv = int(round((lastVal-firstVal)/dv))+1
            for i in xrange(nv):
                v = firstVal + i*dv
                fields.append(v)
    return fields
