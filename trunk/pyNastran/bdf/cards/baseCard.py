# pylint: disable=R0904,R0902,C0111,C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import izip

from pyNastran.bdf.fieldWriter import print_card, is_same
                               #print_card_8, set_default_if_blank, print_card
#from pyNastran.bdf.fieldWriter16 import print_card_16
from pyNastran.bdf.bdfInterface.assign_type import interpret_value
from pyNastran.bdf.bdfInterface.BDF_Card import wipe_empty_fields


class BaseCard(object):
    def __init__(self):
        pass

    def comment(self):
        if hasattr(self, '_comment'):
            return '%s' % self._comment
        return ''

    def _test_update_fields(self):
        n = 1
        while 1:
            try:
                self.update_field(n, 1.0)  # dummy updating the field
            except IndexError:
                return
            except KeyError:
                return

    def update_field(self, n, value):
        try:
            key_name = self._field_map[n]
            setattr(self, key_name, value)
        except KeyError:
            self._update_field_helper(n, value)

    def get_field(self, n):
        try:
            key_name = self._field_map[n]
            value = getattr(self, key_name)
        except KeyError:
            value = self._get_field_helper(n)
        return value

    def writeCodeAster(self):
        return ('# skipping %s  because writeCodeAster is not implemented\n'
                % self.type)

    def writeCodeAsterLoad(self, model, gridWord='node'):
        return ('# skipping %s (lid=%s) because writeCodeAsterLoad is '
                'not implemented\n' % (self.type, self.lid))

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        print('# skipping _verify (type=%s) because _verify is '
              'not implemented\n' % self.type)

    def isSameFields(self, fields1, fields2):
        for (field1, field2) in izip(fields1, fields2):
            if not is_same(field1, field2):
                return False
        return True

    def buildTableLines(self, fields, nStart=1, nEnd=0):
        """
        Builds a table of the form:

        +--------+-------+--------+-------+--------+-------+-------+-------+
        | DESVAR | DVID1 | DVID2  | DVID3 |DVID4   | DVID5 | DVID6 | DVID7 |
        +--------+-------+--------+-------+--------+-------+-------+-------+
        |        | DVID8 | -etc.- |       |        |
        +--------+-------+--------+-------+--------+
        |        |  UM   | VAL1   | VAL2  | -etc.- |
        +--------+-------+--------+-------+--------+

        and then pads the rest of the fields with None's

        :param fields: the fields to enter, including DESVAR
        :type fields:  list of values
        :param nStart: the number of blank fields at the start of the
                       line (default=1)
        :param nStart: int
        :param nEnd:   the number of blank fields at the end of the
                       line (default=0)
        :param nEnd:   int

        .. note:: will be used for DVPREL2, RBE1, RBE3
        .. warning:: only works for small field format???
        """
        fieldsOut = []
        n = 8 - nStart - nEnd

        # pack all the fields into a list.  Only consider an entry as isolated
        for (i, field) in enumerate(fields):
            fieldsOut.append(field)
            if i > 0 and i % n == 0:  # beginning of line
                #print("i=%s" %(i))
                #pad = [None]*(i+j)
                #fieldsOut += pad
                fieldsOut += [None] * (nStart + nEnd)

        # make sure they're aren't any trailing None's (from a new line)
        fieldsOut = wipe_empty_fields(fieldsOut)

        # push the next key (aka next fields[0]) onto the next line
        nSpaces = 8 - (len(fieldsOut)) % 8  # puts UM onto next line
        if nSpaces < 8:
            fieldsOut += [None] * nSpaces
        return fieldsOut

    def isSameCard(self, card):
        if self.type != card.type:
            return False
        fields1 = self.rawFields()
        fields2 = card.rawFields()
        return self.isSameFields(fields1, fields2)

    def printRawFields(self, size=8):
        """A card's raw fields include all defaults for all fields"""
        list_fields = self.rawFields()
        return print_card(list_fields, size=size)

    def reprFields(self):
        return self.rawFields()

    def print_card(self, size=8):
        list_fields = self.reprFields()
        return self.comment() + print_card(list_fields, size=size)

    def repr_card(self, size=8):
        list_fields = self.reprFields()
        return print_card(list_fields, size=size)

    def __repr__(self):
        """
        Prints a card in the simplest way possible
        (default values are left blank).
        """
        try:
            return self.print_card()
        except:
            print('problem printing %s card' % self.type)
            fields = self.reprFields()
            print("fields = ", fields)
            raise


class Property(BaseCard):
    def __init__(self, card, data):
        assert card is None or data is None

    def Pid(self):
        """
        returns the property ID of an property

        :param self:  the Property pointer
        :returns pid: the Property ID
        :type pid:    int
        """
        return self.pid

    def Mid(self):
        """
        returns the material ID of an element

        :param self:  the Property pointer
        :returns mid: the Material ID
        :type mid:    int
        """
        if isinstance(self.mid, int):
            return self.mid
        else:
            return self.mid.mid

    def cross_reference(self, model):
        msg = ' which is required by %s pid=%s' % (self.type, self.pid)
        self.mid = model.Material(self.mid, msg)


class Material(BaseCard):
    """Base Material Class"""
    def __init__(self, card, data):
        BaseCard.__init__(self)

    def cross_reference(self, model):
        pass

    def Mid(self):
        """
        returns the material ID of an element

        :param self:  the Property pointer
        :returns mid: the Material ID
        :type mid:    int
        """
        return self.mid


class Element(BaseCard):
    pid = 0  # CONM2, rigid

    def __init__(self, card, data):
        BaseCard.__init__(self)
        assert card is None or data is None
        #: the list of node IDs for an element (default=None)
        self.nodes = None
        #self.nids = []

    def Pid(self):
        """
        returns the property ID of an element
        :param self:  the Element pointer
        :returns pid: the Property ID
        :type pid:    int
        """
        if isinstance(self.pid, int):
            return self.pid
        else:
            return self.pid.pid

    def nodePositions(self, nodes=None):
        """returns the positions of multiple node objects"""
        if not nodes:
            nodes = self.nodes

        positions = []
        for node in nodes:
            if node is not None:
                positions.append(node.Position())
            else:
                positions.append(None)
        return positions

    def _nodeIDs(self, nodes=None, allowEmptyNodes=False, msg=''):
        """returns nodeIDs for repr functions"""
        try:
            if not nodes:
                nodes = self.nodes

            if allowEmptyNodes:
                nodes2 = []
                for i, node in enumerate(nodes):
                    #print("node=%r type=%r" % (node, type(node)))
                    if node == 0 or node is None:
                        nodes2.append(None)
                    elif isinstance(node, int):
                        nodes2.append(node)
                    else:
                        nodes2.append(node.nid)
                return nodes2
            else:
                try:
                    nodeIDs = []
                    for i, node in enumerate(nodes):
                        #print("node=%r type=%r" % (node, type(node)))
                        if isinstance(node, int):
                            nodeIDs.append(node)
                        else:
                            nodeIDs.append(node.nid)

                    #if isinstance(nodes[0], int):
                        #nodeIDs = [node for node in nodes]
                    #else:
                        #nodeIDs = [node.nid for node in nodes]
                except:
                    print("type=%s nodes=%s allowEmptyNodes=%s\nmsg=%s" % (
                          self.type, nodes, allowEmptyNodes, msg))
                    raise
                assert 0 not in nodeIDs, 'nodeIDs = %s' % nodeIDs
                return nodeIDs
        except:
            print("type=%s nodes=%s allowEmptyNodes=%s\nmsg=%s" % (
                  self.type, nodes, allowEmptyNodes, msg))
            raise

    def prepareNodeIDs(self, nids, allowEmptyNodes=False):
        """Verifies all node IDs exist and that they're integers"""
        self.nodes = []
        for nid in nids:
            if isinstance(nid, int):
                self.nodes.append(nid)
            elif nid is None and allowEmptyNodes:
                self.nodes.append(None)
            else:  # string???
                #self.nodes.append(int(nid))
                raise RuntimeError('this element may not have missing '
                                   'nodes...nids=%s allowEmptyNodes=False'
                                   % nids)

    @property  # I think this means you can just call it as an attribute...
    def faces(self):
        """
        Gets the faces of the element

        :returns: dictionary with face number as the keys and
                  a list of nodes (integer pointers) as the values.
        ..note::  The order of the nodes are consistent with ANSYS numbering.

        Example
        =======
        >>> print element.faces  # TODO: is this right???
        """
        faces = {}
        nodes = self.nodeIDs()
        if self.type.find('HEX') != -1:  # CHEXA
            faces[1] = [nodes[0], nodes[1], nodes[2], nodes[3]]
            faces[2] = [nodes[0], nodes[1], nodes[5], nodes[4]]
            faces[3] = [nodes[1], nodes[2], nodes[6], nodes[5]]
            faces[4] = [nodes[2], nodes[3], nodes[7], nodes[6]]
            faces[5] = [nodes[3], nodes[0], nodes[4], nodes[7]]
            faces[6] = [nodes[4], nodes[5], nodes[6], nodes[7]]
        elif self.type.find('TET') != -1:  # CTETRA
            faces[1] = [nodes[0], nodes[1], nodes[2]]
            faces[2] = [nodes[0], nodes[1], nodes[3]]
            faces[3] = [nodes[1], nodes[2], nodes[3]]
            faces[4] = [nodes[2], nodes[0], nodes[3]]
        elif self.type.find('PYR') != -1:  # PYRAMID
            faces[1] = [nodes[0], nodes[1], nodes[2], nodes[3]]
            faces[2] = [nodes[0], nodes[1], nodes[4]]
            faces[3] = [nodes[1], nodes[2], nodes[4]]
            faces[4] = [nodes[2], nodes[3], nodes[4]]
            faces[5] = [nodes[3], nodes[0], nodes[4]]
        elif self.type.find('PEN') != -1:  # CPENTA
            faces[1] = [nodes[0], nodes[1], nodes[2]]
            faces[2] = [nodes[3], nodes[4], nodes[5]]
            faces[3] = [nodes[0], nodes[1], nodes[4], nodes[3]]
            faces[4] = [nodes[1], nodes[2], nodes[5], nodes[4]]
            faces[5] = [nodes[2], nodes[0], nodes[3], nodes[5]]
        elif self.type.find('QUAD') != -1:  # CQUADx
            # both sides
            faces[1] = [nodes[0], nodes[1], nodes[2], nodes[3]]
            faces[2] = [nodes[1], nodes[0], nodes[3], nodes[2]]
        elif self.type.find('TRI') != -1:  # CTRIAx
            # both sides
            faces[1] = [nodes[0], nodes[1], nodes[2]]
            faces[2] = [nodes[1], nodes[0], nodes[2]]
        else:
            faces = None
        return faces

    def nodes2face(self, nodes):
        """
        returns the face number that matches the list of nodes input

        :param nodes:   list of nodes
        :returns faces: the face number as an integer

        ..warning:: It's assumed you have the nodes in the proper order.
        """
        assert isinstance(nodes, list), 'nodes=%s' % str(nodes)
        face = None
        for i in self.faces.keys():
            chck = self.faces[i]
            if nodes == chck:
                face = i
        return face


def expand_thru(fields, set_fields=True, sort_fields=False):
    """
    Expands a list of values of the form [1,5,THRU,9,13]
    to be [1,5,6,7,8,9,13]

    :param fields:      the fields to expand
    :param set_fields:  should the fields be converted to a set and then back
                        to a list? This is useful for [2, 'THRU' 5, 1] (default=True)
    :param sort_fields: should the fields be sorted at the end? (default=False)
    """
    # ..todo:  should this be removed...is the field capitalized when read in?
    fields = [field.upper() if isinstance(field, basestring) else field for field in fields]

    if isinstance(fields, int):
        return [fields]
    if len(fields) == 1:
        return [int(fields[0])]
    out = []
    nFields = len(fields)
    i = 0
    while(i < nFields):
        if fields[i] == 'THRU':
            istart = int(fields[i - 1])
            iend = int(fields[i + 1])
            for j in xrange(istart, iend + 1): # adding 1 to iend for the xrange offset
                out.append(j)
            i += 2
        else:
            out.append(int(fields[i]))
            i += 1

    if set_fields:
        out = list(set(out))
    if sort_fields:
        out.sort()
    return out


def expand_thru_by(fields, set_fields=True, sort_fields=False):
    """
    Expands a list of values of the form [1,5,THRU,9,BY,2,13]
    to be [1,5,7,9,13]

    :param fields:      the fields to expand
    :param set_fields:  should the fields be converted to a set and then back
                        to a list? This is useful for [2, 'THRU' 5, 1] (default=True)
    :param sort_fields: should the fields be sorted at the end? (default=False)

    .. todo:: not tested
    .. note:: used for QBDY3, ???
    """
    # ..todo:  should this be removed...is the field capitalized when read in?
    fields = [field.upper() if isinstance(field, basestring) else field for field in fields]

    if len(fields) == 1:
        return [interpret_value(fields[0])]
    out = []
    nFields = len(fields)
    i = 0
    by = 1
    while(i < nFields):
        if fields[i] == 'THRU':
            by = 1
            byCase = False
            if i + 2 < nFields and fields[i + 2] == 'BY':
                by = interpret_value(fields[i + 3])
            else:
                by = 1
                byCase = True
            minValue = interpret_value(fields[i - 1])
            maxValue = interpret_value(fields[i + 1])
            maxR = int((maxValue - minValue) // by + 1)  # max range value

            for j in xrange(0, maxR):  # +1 is to include final point
                value = minValue + by * j
                out.append(value)

            if byCase:  # null/standard case
                i += 2
            else:     # BY case
                i += 3
        else:
            out.append(interpret_value(fields[i]))
            i += 1

    if set_fields:
        out = list(set(out))
    if sort_fields:
        out.sort()
    return out


def expand_thru_exclude(self, fields):
    """
    Expands a list of values of the form [1,5,THRU,11,EXCEPT,7,8,13]
    to be [1,5,6,9,10,11,13]

    .. warning:: hasnt been tested
    """
    # ..todo:  should this be removed...is the field capitalized when read in?
    fields = [interpret_value(field.upper()) if isinstance(field, basestring) else field for field in fields]

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


def collapse_thru_by(fields, get_packs=False):
    """
    :param fields:    the list of fields to collapse
    :param get_packs: get the list of packs so "special" formatting can be done

    fields              packs
    [1, 2, 3...150]  -> [1, 150, 1]
    [1, 3, 5...150]  -> [1, 150, 2]
    """
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    if get_packs:
        return packs
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


def collapse_thru_packs(fields):
    assert 'THRU' not in fields, fields
    fields.sort()
    packs = condense(fields)
    singles, doubles = build_thru_packs(packs, maxDV=1)

    #assert fields == expand_thru_by(fields2), fields2  # why doesn't this work?
    return singles, doubles


def condense(valueList):
    """
    Builds a list of packs (list of 3 values representing the first, last,
    and delta values for condensing a SET card.

    .. seealso:: build_thru
    """
    if len(valueList) == 0:
        return []
    if len(valueList) == 1:
        return [[valueList[0], valueList[0], 1]]
    valueList.sort()
    packs = []

    dvOld = None
    firstVal = valueList[0]
    lastVal = firstVal

    for val in valueList[1:]:
        try:
            dv = val - lastVal
        except TypeError:
            print("lastVal=%r val=%r" % (lastVal, val))
            print("valueList=%r" % valueList)
            raise

        # sets up the first item of the pack
        if dvOld is None:
            dvOld = dv

        # fill up the pack
        if dvOld == dv:
            lastVal = val
        else:
            packs.append([firstVal, lastVal, dvOld])
            lastVal = val
            dvOld = None
            firstVal = val

    # fills the last pack
    if dvOld == dv:
        packs.append([firstVal, val, dv])
    else:
        packs.append([firstVal, val, dvOld])
    return packs


def build_thru_packs(packs, maxDV=1):
    """
    # invalid
    SET1,4000, 1, 3, THRU, 10, 20, THRU, 30

    # valid
    SET1,4000, 1
    SET1,4000, 3,  THRU, 10
    SET1,4000, 20, THRU, 30

    returns
      singles = [1]
      doubles = [[3, 'THRU', 10], [20, 'THRU', 30]]
    """
    singles = []
    doubles = []
    for (firstVal, lastVal, by) in packs:
        if firstVal == lastVal:
            singles.append(firstVal)
        else:
            if lastVal - firstVal > 2:
                double = [firstVal, 'THRU', lastVal]
                doubles.append(double)
            else:
                for val in xrange(lastVal, firstVal+1, by):
                    singles.append(val)
    return singles, doubles


def build_thru(packs, maxDV=None):
    """
    Takes a pack [1,7,2] and converts it into fields used by a SET card.
    The values correspond to the first value, last value, and delta in the
    list.  This means that [1,1001,2] represents 500 values.
    [1,1001,1] represents 1001 values and will be written as [1,THRU,1001]..

    :param packs: list of packs (list of 3 values: [first, last, delta] )
    :param maxDV: integer defining the max allowable delta between two values
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
                if lastVal - firstVal > 4 * dv:
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

    :param packs: list of packs (list of 3 values: [first, last, delta] )
    :param maxDV: integer defining the max allowable delta between two values
                  (default=None; no limit)
    """
    fields = []
    for (firstVal, lastVal, dv) in packs:
        if lastVal - firstVal > 4 * dv:
            fields.append(firstVal)
            fields.append('THRU')
            fields.append(lastVal)
            fields.append('BY')
            fields.append(dv)
        else:
            nv = int(round((lastVal - firstVal) / dv)) + 1
            for i in xrange(nv):
                v = firstVal + i * dv
                fields.append(v)
    return fields