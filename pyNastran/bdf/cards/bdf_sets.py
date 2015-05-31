from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import zip, range

from pyNastran.bdf.cards.baseCard import (BaseCard, _node_ids, expand_thru,
    collapse_thru, collapse_thru_packs)
from pyNastran.bdf.field_writer_8 import print_card_8, print_int_card_blocks
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    components, components_or_blank, fields, integer_or_string, string,
    string_or_blank, integer_string_or_blank)

"""
All constraint cards are defined in this file.  This includes:

* sets
  * SET1, SET3, RADSET # ??? RADSET
* asets - aset, aset1
* bsets - bset, bset1
* csets - cset, cset1
* qsets - qset, qset1
* usets - uset, uset1  # USET 1 is not supported

The superelement sets start with SE:
* se_bsets - sebset, sebset1
* se_csets - secset, secset1
* se_qsets - seqset, seqset1
* se_usets - seuset, seuset1
*se_sets
  * SESET
  * SEQSEP

 #* Set
 #* SetSuper

+------------+-----------------+
| Entry Type | Equivalent Type |
+------------+-----------------+
|  SEQSETi   | QSETi           |
+------------+-----------------+
|  SESUP     | SUPORT          |
+------------+-----------------+
|  SECSETi   | CSETi           |
+------------+-----------------+
|  SEBSETi   | BSETi           |
+------------+-----------------+
"""

class Set(BaseCard):
    """Generic Class all SETx cards inherit from"""

    def __init__(self, card, data):
        #:  Unique identification number. (Integer > 0)
        self.sid = None
        #:  list of IDs in the SETx
        self.IDs = None

    def cleanIDs(self):
        """eliminates duplicate IDs from self.IDs and sorts self.IDs"""
        self.IDs = list(set(self.IDs))
        self.IDs.sort()

    def SetIDs(self):
        """gets the IDs of the SETx"""
        return collapse_thru(self.IDs)

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def __repr__(self):
        return self.comment() + print_card_8(self.repr_fields())

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment() + print_card_8(card)


class SetSuper(Set):
    """Generic Class all Superelement SETx cards inherit from."""
    def __init__(self, card, data):
        #:  Superelement identification number. Must be a primary superelement.
        #:  (Integer >= 0)
        self.seid = None
        #:  list of IDs in the SESETx
        self.IDs = None


class ABCQSet(Set):
    """
    Generic Class ASET, BSET, CSET, QSET cards inherit from.

    Defines degrees-of-freedom in the analysis set (A-set)

    +------+-----+----+-----+------+-----+----+-----+----+
    | ASET | ID1 | C1 | ID2 |  C2  | ID3 | C3 | ID4 | C4 |
    +------+-----+----+-----+------+-----+----+-----+----+
    | ASET | 16  |  2 |  23 | 3516 |  1  | 4  |     |    |
    +------+-----+----+-----+------+-----+----+-----+----+
    """
    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment

        #:  Identifiers of grids points. (Integer > 0)
        self.IDs = []
        self.components = []

        nterms = len(card) // 2
        for n in range(nterms):
            i = n * 2 + 1
            ID = integer(card, i, 'ID' + str(n))
            component = components(card, i + 1, 'component' + str(n))
            self.IDs.append(ID)
            self.components.append(component)

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = [self.type]  # ASET, BSET
        for (ID, comp) in zip(self.IDs, self.components):
            list_fields += [ID, comp]
        return list_fields

    def __repr__(self):
        list_fields = self.raw_fields()
        return self.comment() + print_card_8(list_fields)


class ASET(ABCQSet):
    """
    Defines degrees-of-freedom in the analysis set (A-set).

    +------+-----+----+-----+------+-----+----+-----+----+
    | ASET | ID1 | C1 | ID2 |  C2  | ID3 | C3 | ID4 | C4 |
    +------+-----+----+-----+------+-----+----+-----+----+
    | ASET | 16  |  2 |  23 | 3516 |  1  | 4  |     |    |
    +------+-----+----+-----+------+-----+----+-----+----+
    """
    type = 'ASET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)

    def cross_reference(self, model):
        pass

class BSET(ABCQSet):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.

    +------+-----+----+-----+------+-----+----+-----+----+
    | BSET | ID1 | C1 | ID2 |  C2  | ID3 | C3 | ID4 | C4 |
    +------+-----+----+-----+------+-----+----+-----+----+
    | BSET | 16  |  2 |  23 | 3516 |  1  | 4  |     |    |
    +------+-----+----+-----+------+-----+----+-----+----+
    """
    type = 'BSET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)


class CSET(ABCQSet):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.

    +------+-----+----+-----+------+-----+----+-----+----+
    | CSET | ID1 | C1 | ID2 |  C2  | ID3 | C3 | ID4 | C4 |
    +------+-----+----+-----+------+-----+----+-----+----+
    | CSET | 16  |  2 |  23 | 3516 |  1  | 4  |     |    |
    +------+-----+----+-----+------+-----+----+-----+----+
    """
    type = 'CSET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)


class QSET(ABCQSet):
    """
    Defines generalized degrees-of-freedom (q-set) to be used for dynamic
    reduction or component mode synthesis.

    +------+-----+----+-----+------+-----+----+-----+----+
    | QSET | ID1 | C1 | ID2 |  C2  | ID3 | C3 | ID4 | C4 |
    +------+-----+----+-----+------+-----+----+-----+----+
    | QSET | 16  | 2  | 23  | 3516 |  1  |  4 |     |    |
    +------+-----+----+-----+------+-----+----+-----+----+
    """
    type = 'QSET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)


class ABQSet1(Set):
    """
    Generic Class ASET1, BSET1, QSET1 cards inherit from.

    Defines degrees-of-freedom in the analysis set (a-set).

    +--=----+-----+-----+------+------+-----+-----+-----+-----+
    | ASET1 |  C  | ID1 | ID2  | ID3  | ID4 | ID5 | ID6 | ID7 |
    +-------+-----+-----+------+------+-----+-----+-----+-----+
    |       | ID8 | ID9 |      |      |     |     |     |     |
    +-------+-----+-----+------+------+-----+-----+-----+-----+
    | ASET1 |  C  | ID1 | THRU | ID2  |     |     |     |     |
    +-------+-----+-----+------+------+-----+-----+-----+-----+
    """

    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment
        #:  Component number. (Integer zero or blank for scalar points or any
        #:  unique combination of the Integers 1 through 6 for grid points with
        #:  no embedded blanks.)
        self.components = components_or_blank(card, 1, 'components', 0)

        nfields = len(card)
        IDs = []
        i = 1
        for ifield in range(2, nfields):
            ID = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if ID:
                i += 1
                IDs.append(ID)
        #IDs = fields(integer_or_string, card, 'ID', i=2, j=nfields)
        #:  Identifiers of grids points. (Integer > 0)
        self.IDs = expand_thru(IDs)

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = [self.type, self.components] + collapse_thru(self.IDs)
        return list_fields

    def __repr__(self):
        list_fields = self.raw_fields()
        return self.comment() + print_card_8(list_fields)


class ASET1(ABQSet1):
    """
    Defines degrees-of-freedom in the analysis set (a-set)

    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    | ASET1 |  C  | ID1 |  ID2 | ID3 | ID4 | ID5 | ID6 | ID7 |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    |       | ID8 | ID9 |      |     |     |     |     |     |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    | ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    """
    type = 'ASET1'

    def __init__(self, card=None, data=None, comment=''):
        ABQSet1.__init__(self, card, data, comment)

    def cross_reference(self, model):
        pass


class BSET1(ABQSet1):
    type = 'BSET1'

    def __init__(self, card=None, data=None, comment=''):
        ABQSet1.__init__(self, card, data, comment)


class CSET1(Set):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.

    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    | CSET1 |  C  | ID1 |  ID2 | ID3 | ID4 | ID5 | ID6 | ID7 |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    |       | ID8 | ID9 |      |     |     |     |     |     |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    | CSET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    | CSET1 |  ,, | ALL |      |     |     |     |     |     |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    """
    type = 'CSET1'

    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment

        #:  Identifiers of grids points. (Integer > 0)
        self.IDs = []
        if string_or_blank(card, 2, 'C') == 'ALL':
            self.components = '123456'
        else:
            self.components = components(card, 1, 'components')

            IDs2 = []
            ii = 1
            for ifield in range(2, len(card)):
                integer_or_string(card, ifield, 'ID' % ii)
                ii += 1
            IDs = fields(integer_or_string, 'ID', i=2, j=len(card))
            assert IDs2 == IDs
            self.IDs = expand_thru(IDs)

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['CSET1', self.components] + collapse_thru(self.IDs)
        return list_fields

    def __repr__(self):
        list_fields = self.raw_fields()
        return self.comment() + print_card_8(list_fields)


class QSET1(ABQSet1):
    """
    Defines generalized degrees-of-freedom (q-set) to be used for dynamic
    reduction or component mode synthesis.
    """
    type = 'QSET1'

    def __init__(self, card=None, data=None, comment=''):
        ABQSet1.__init__(self, card, data, comment)

    def cross_reference(self, model):
        pass


class SET1(Set):
    """
    Defines a list of structural grid points or element identification
    numbers.

    +------+--------+--------+-----+------+-----+-----+------+-----+
    | SET1 | SID    |   ID1  | ID2 | ID3  | ID4 | ID5 | ID6  | ID7 |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    |      |  ID8   | -etc.- |     |      |     |     |      |     |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    | SET1 |  3     |   31   | 62  |  93  | 124 | 16  |  17  | 18  |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    |      |   19   |        |     |      |     |     |      |     |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    | SET1 |  6     |  29    | 32  | THRU | 50  | 61  | THRU | 70  |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    |      |   17   |  57    |     |      |     |     |      |     |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    """
    type = 'SET1'

    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment
        #:  Unique identification number. (Integer > 0)
        self.sid = integer(card, 1, 'sid')

        self.IDs = []
        IDs = []
        i = 1
        for ifield in range(2, len(card)):
            ID = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if ID:
                i += 1
                IDs.append(ID)
        #IDs = fields(integer_or_string, card, 'ID', i=2, j=len(card))

        self.isSkin = False
        i = 0
        if isinstance(IDs[0], str) and IDs[0] == 'SKIN':
            self.isSkin = True
            i += 1

        #:  List of structural grid point or element identification numbers.
        #:  (Integer > 0 or 'THRU'; for the 'THRU' option, ID1 < ID2 or 'SKIN';
        #:  in field 3)
        self.IDs = expand_thru(IDs[i:])
        self.cleanIDs()

    def __eq__(self, set1):
        self.cleanIDs()
        set1.cleanIDs()
        if self.IDs == set1.IDs:
            return True
        return False

    def symmetric_difference(self, set1):
        s1 = set(self.IDs)
        s2 = set(set1.IDs)
        return s1.symmetric_difference(s2)

    def add_set(self, set1):
        self.IDs += set1.IDs
        self.cleanIDs()

    def IsSkin(self):
        return self.isSkin

    def raw_fields(self):
        skin = []
        if self.isSkin:
            skin = ['SKIN']

        return ['SET1', self.sid] + skin + self.IDs

    def write_card(self, size=8, is_double=False):
        skin = []
        if self.isSkin:
            skin = ['SKIN']

        #return self.comment() + print_card_8(['SET1', self.sid] + skin + self.IDs)

        field_packs = []
        singles, doubles = collapse_thru_packs(self.IDs)
        if singles:
            field_packs.append(['SET1', self.sid] + skin + singles)
        if doubles:
            for pack in doubles:
                field_packs.append(['SET1', self.sid] + skin + pack)

        msg = []
        for field_pack in field_packs:
            msg.append(print_card_8(field_pack))
        return ''.join(msg)


class SET3(Set):
    """
    Defines a list of grids, elements or points.

    +------+-----+-------+-----+-----+-----+-----+-----+-----+
    | SET3 | SID |  DES  | ID1 | ID2 | ID3 | ID4 | ID5 | ID6 |
    +------+-----+-------+-----+-----+-----+-----+-----+-----+
    |      | ID7 |  ID8  | etc |     |     |     |     |     |
    +------+-----+-------+-----+-----+-----+-----+-----+-----+

    Example
    +------+-----+-------+-----+----+
    | SET3 |  1  | POINT | 11  | 12 |
    +------+-----+-------+-----+----+
    """
    type = 'SET3'

    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment
        #:  Unique identification number. (Integer > 0)
        self.sid = integer(card, 1, 'sid')

        #:  Set description (Character). Valid options are 'GRID', 'ELEM',
        #:  'POINT' and 'PROP'.
        self.desc = string(card, 2, 'desc')
        assert self.desc in ['GRID', 'POINT', 'ELEM', 'PROP'], 'desc = %r' % self.desc

        #:  Identifiers of grids points, elements, points or properties.
        #:  (Integer > 0)
        self.IDs = []

        IDs = fields(integer_or_string, card, 'ID', i=3, j=len(card))
        self.IDs = expand_thru(IDs)
        self.cleanIDs()

    def union(self, s2):
        assert self.desc == s2.desc, 'self.desc=%r set2.desc=%r' % (self.desc, s2.desc)
        ids1 = set(self.IDs)
        ids2 = set(s2.IDs)
        self.IDs = list(ids1.union(ids2))

    def symmetric_difference(self, s2):
        ids1 = set(self.IDs)
        ids2 = set(s2.IDs)
        return ids1.symmetric_difference(ids2)

    def IsGrid(self):
        if self.desc == 'GRID':
            return True
        return False

    def IsPoint(self):
        if self.desc == 'POINT':
            return True
        return False

    def IsProperty(self):
        if self.desc == 'PROP':
            return True
        return False

    def IsElement(self):
        if self.desc == 'ELEM':
            return True
        return False

    def raw_fields(self):
        """Gets the "raw" card without any processing as a list for printing"""
        list_fields = ['SET3', self.sid, self.desc] + self.SetIDs()
        return list_fields

    def __repr__(self):
        fields_blocks = [
            'SET3',
            [[self.sid, self.desc], False], # these are not all integers
            [self.SetIDs(), True], # these are all integers
        ]
        return self.comment() + print_int_card_blocks(fields_blocks)


class SESET(SetSuper):
    """
    Defines interior grid points for a superelement.
    """
    type = 'SESET'

    def __init__(self, card=None, data=None, comment=''):
        SetSuper.__init__(self, card, data)
        if comment:
            self._comment = comment
        self.seid = integer_or_blank(card, 1, 'seid', 0)
        #:  Grid or scalar point identification number.
        #:  (0 < Integer < 1000000; G1 < G2)
        self.IDs = []

        IDs = fields(integer_or_string, card, 'ID', i=2, j=len(card))
        self.IDs = expand_thru(IDs)
        self.cleanIDs()

    def add_SESET_Object(self, seset):
        self.IDs += seset.IDs
        self.cleanIDs()

    def raw_fields(self):
        return ['SESET', self.seid] + collapse_thru(self.IDs)

    def __repr__(self):
        thruFields = collapse_thru(self.IDs)
        #list_fields = ['SESET', self.seid]

        cards = []
        while 'THRU' in thruFields:
            iThru = thruFields.index('THRU')
            card = print_card_8(['SESET', self.seid] +
                                thruFields[iThru - 1:iThru + 2])
            cards.append(card)
            thruFields = thruFields[0:iThru - 1]+thruFields[iThru + 2:]
        if thruFields:
            card = print_card_8(['SESET', self.seid] + thruFields)
            cards.append(card)
        return ''.join(cards)

    def cross_reference(self, model):
        pass


class SEBSET(ABCQSet):
    """
    Defines boundary degrees-of-freedom to be fixed (b-set) during generalized
    dynamic reduction or component mode calculations.

    +--------+------+-----+------+-----+----+-----+----+
    | SEBSET | SEID | ID1 |  C1  | ID2 | C2 | ID3 | C3 |
    +--------+------+-----+------+-----+----+-----+----+
    | SEBSET |  C   | ID1 | THRU | ID2 |    |     |    |
    +--------+------+-----+------+-----+----+-----+----+
    """
    type = 'SEBSET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)

class SEBSET1(ABQSet1):
    type = 'SEBSET1'

    def __init__(self, card=None, data=None, comment=''):
        ABQSet1.__init__(self, card, data, comment)

    def cross_reference(self, model):
        pass

class SECSET(ABCQSet):
    type = 'SECSET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)

class SECSET1(ABQSet1):
    type = 'SECSET1'

    def __init__(self, card=None, data=None, comment=''):
        ABQSet1.__init__(self, card, data, comment)

    def cross_reference(self, model):
        pass


class SEQSET(ABCQSet):
    type = 'SEQSET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)

class SEQSET1(ABQSet1):
    type = 'SEQSET1'

    def __init__(self, card=None, data=None, comment=''):
        ABQSet1.__init__(self, card, data, comment)

    def cross_reference(self, model):
        pass


class SEQSEP(SetSuper):  # not integrated...is this an SESET ???
    """
    Used with the CSUPER entry to define the correspondence of the
    exterior grid points between an identical or mirror-image superelement
    and its primary superelement.
    """
    type = 'SEQSEP'

    def __init__(self, card=None, data=None, comment=''):
        SetSuper.__init__(self, card, data)
        if comment:
            self._comment = comment
        #: Identification number for secondary superelement. (Integer >= 0).
        self.ssid = integer(card, 1, 'ssid')
        #: Identification number for the primary superelement. (Integer >= 0).
        self.psid = integer(card, 2, 'psid')
        #: Exterior grid point identification numbers for the primary
        #: superelement. (Integer > 0)
        self.IDs = []

        IDs = fields(integer_or_string, card, i=3, j=len(card))
        self.IDs = expand_thru(IDs)
        self.cleanIDs()

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['SEQSEP', self.ssid, self.psid] + self.SetIDs()
        return list_fields


class RADSET(Set):  # not integrated
    """
    Specifies which radiation cavities are to be included for
    radiation enclosure analysis.

    +--------+----------+----------+----------+----------+----------+----------+----------+
    | RADSET | ICAVITY1 | ICAVITY2 | ICAVITY3 | ICAVITY4 | ICAVITY5 | ICAVITY6 | ICAVITY7 |
    +--------+----------+----------+----------+----------+----------+----------+----------+
    |        | ICAVITY8 | ICAVITY9 | -etc.-   |          |          |          |          |
    +--------+----------+----------+----------+----------+----------+----------+----------+
    """
    type = 'RADSET'

    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment
        self.seid = integer(card, 1, 'seid')
        #: Grid or scalar point identification number.
        #: (0 < Integer < 1000000; G1 < G2)
        self.IDs = []

        IDs = fields(integer_or_string, card, 'ID', i=2, j=len(card))
        self.IDs = expand_thru(IDs)
        self.cleanIDs()

    def addRadsetObject(self, radset):
        self.IDs += radset.IDs
        self.cleanIDs()

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['RADSET', self.seid] + self.SetIDs()
        return list_fields


class USET(Set):
    """
    Defines a degrees-of-freedom set.

    +------+-------+-----+------+-----+----+-----+----+
    | USET | SNAME | ID1 |  C1  | ID2 | C2 | ID3 | C3 |
    +------+-------+-----+------+-----+----+-----+----+
    | USET |  JUNK | ID1 | THRU | ID2 |    |     |    |
    +------+-------+-----+------+-----+----+-----+----+
    """
    type = 'USET'
    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment

        self.name = string(card, 1, 'name')

        #:  Identifiers of grids points. (Integer > 0)
        self.components = []
        self.IDs = []

        nsets = (len(card) - 1) // 2
        for n in range(nsets):
            i = n * 2 + 2
            ID = integer(card, i, 'node_id' + str(n))
            component = components(card, i + 1, 'component' + str(n))
            self.components.append(component)
            self.IDs.append(ID)

    def raw_fields(self):
        """
        gets the "raw" card without any processing as a list for printing
        """
        list_fields = ['USET', self.name]
        for (component, ID) in zip(self.components, self.node_ids):
            list_fields += [ID, component]
        return list_fields

    def cross_reference(self, model):
        msg = ' which is required by %s name=%s' % (self.type, self.name)
        self.IDs = model.Nodes(self.IDs, allowEmptyNodes=True, msg=msg)

    @property
    def node_ids(self):
        msg = ' which is required by %s name=%s' % (self.type, self.name)
        return _node_ids(self, self.IDs, allowEmptyNodes=True, msg=msg)

