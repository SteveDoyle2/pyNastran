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
+============+=================+
|  SEQSETi   | QSETi           |
+------------+-----------------+
|  SESUP     | SUPORT          |
+------------+-----------------+
|  SECSETi   | CSETi           |
+------------+-----------------+
|  SEBSETi   | BSETi           |
+------------+-----------------+
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types
from six.moves import zip, range

from pyNastran.utils import integer_types
from pyNastran.bdf.cards.base_card import (
    BaseCard, _node_ids, expand_thru
)
from pyNastran.bdf.cards.collpase_card import collapse_thru, condense, build_thru_packs
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_or_string,
    parse_components, components_or_blank as fcomponents_or_blank,
    fields, string, integer_string_or_blank)


class Set(BaseCard):
    """Generic Class all SETx cards inherit from"""

    def __init__(self):
        #:  Unique identification number. (Integer > 0)
        self.sid = None
        #:  list of IDs in the SETx
        self.ids = []

    def clean_ids(self):
        """eliminates duplicate IDs from self.IDs and sorts self.IDs"""
        self.ids = list(set(self.ids))
        self.ids.sort()

    #def cleanIDs(self):
        #self.clean_ids()

    #def SetIDs(self):
        #"""gets the IDs of the SETx"""
        #return collapse_thru(self.ids)

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def __repr__(self):
        return self.comment + print_card_8(self.repr_fields())

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SetSuper(Set):
    """Generic Class all Superelement SETx cards inherit from."""
    def __init__(self):
        Set.__init__(self)
        #:  Superelement identification number. Must be a primary superelement.
        #:  (Integer >= 0)
        self.seid = None
        #:  list of IDs in the SESETx
        self.ids = None


class ABCQSet(Set):
    """
    Generic Class ASET, BSET, CSET, QSET cards inherit from.

    Defines degrees-of-freedom in the analysis set (A-set)

    +------+-----+----+-----+------+-----+----+-----+----+
    |   1  |   2 |  3 |  4  |   5  |  6  |  7 |   8 |  9 |
    +======+=====+====+=====+======+=====+====+=====+====+
    | ASET | ID1 | C1 | ID2 |  C2  | ID3 | C3 | ID4 | C4 |
    +------+-----+----+-----+------+-----+----+-----+----+
    | ASET | 16  |  2 |  23 | 3516 |  1  | 4  |     |    |
    +------+-----+----+-----+------+-----+----+-----+----+
    """
    type = 'ABCQSet'
    def __init__(self, ids, components, comment=''):
        Set.__init__(self)
        if comment:
            self.comment = comment
        #:  Identifiers of grids points. (Integer > 0)
        self.ids = ids
        self.components = components

    def validate(self):
        assert isinstance(self.ids, list), type(self.ids)
        assert isinstance(self.components, list), type(self.components)
        assert len(self.ids) == len(self.components), 'len(ids)=%s len(components)=%s' % (len(self.ids), len(self.components))

    @classmethod
    def add_card(cls, card, comment=''):
        ids = []
        components = []

        nterms = len(card) // 2
        for n in range(nterms):
            i = n * 2 + 1
            idi = integer(card, i, 'ID' + str(n))
            component = parse_components(card, i + 1, 'component' + str(n))
            ids.append(idi)
            components.append(component)
        return cls(ids, components, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        ids = [data[0]]
        components = [data[1]]
        return cls(ids, components, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s' % self.type
        self.ids = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.ids_ref = self.ids

    def uncross_reference(self):
        self.ids = self.node_ids
        del self.ids_ref

    @property
    def node_ids(self):
        msg = ' which is required by %s' % self.type
        return _node_ids(self, self.ids, allow_empty_nodes=True, msg=msg)

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = [self.type]  # ASET, BSET
        for (idi, comp) in zip(self.node_ids, self.components):
            list_fields += [idi, comp]
        return list_fields

    def __repr__(self):
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)


class SuperABCQSet(Set):
    """
    Generic Class ASET, BSET, CSET, QSET cards inherit from.

    Defines degrees-of-freedom in the analysis set (A-set)

    +--------+------+-----+----+-----+------+-----+-----+-----+
    |   1    |  2   |  3  | 4  |  5  |  6   |  7  |  8  |  9  |
    +========+======+=====+====+=====+======+=====+=====+=====+
    | SEBSET | SEID | ID1 | C1 | ID2 |  C2  | ID3 | C3  |     |
    +--------+------+-----+----+-----+------+-----+-----+-----+
    | SEBSET | 100  | 16  |  2 |  23 | 3516 |  1  | 4   |     |
    +--------+------+-----+----+-----+------+-----+-----+-----+
    """
    type = 'SuperABCQSet'
    def __init__(self, seid, ids, components, comment=''):
        Set.__init__(self)
        if comment:
            self.comment = comment

        self.seid = seid
        #:  Identifiers of grids points. (Integer > 0)
        self.ids = ids
        self.components = components

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        ids = []
        components = []

        nterms = len(card) // 2
        for n in range(nterms):
            i = n * 2 + 2
            idi = integer(card, i, 'ID' + str(n))
            component = parse_components(card, i + 1, 'component' + str(n))
            ids.append(idi)
            components.append(component)
        return cls(seid, ids, components, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s seid=%s' % (self.type, self.seid)
        self.ids = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.ids_ref = self.ids

    def uncross_reference(self):
        self.ids = self.node_ids
        del self.ids_ref

    @property
    def node_ids(self):
        msg = ' which is required by %s seid=%s' % (self.type, self.seid)
        return _node_ids(self, self.ids, allow_empty_nodes=True, msg=msg)

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = [self.type, self.seid]  # SEASET, SEBSET
        for (idi, comp) in zip(self.node_ids, self.components):
            list_fields += [idi, comp]
        return list_fields

    def __repr__(self):
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)


class ASET(ABCQSet):
    """
    Defines degrees-of-freedom in the analysis set (A-set).

    +------+-----+----+-----+------+-----+----+-----+----+
    |  1   |  2  | 3  |  4  |  5   |  6  |  7 |  8  | 9  |
    +======+=====+====+=====+======+=====+====+=====+====+
    | ASET | ID1 | C1 | ID2 |  C2  | ID3 | C3 | ID4 | C4 |
    +------+-----+----+-----+------+-----+----+-----+----+
    | ASET | 16  |  2 |  23 | 3516 |  1  | 4  |     |    |
    +------+-----+----+-----+------+-----+----+-----+----+
    """
    type = 'ASET'

    def __init__(self, ids, components, comment=''):
        ABCQSet.__init__(self, ids, components, comment)

class BSET(ABCQSet):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.

    +------+-----+----+-----+------+-----+----+-----+----+
    |  1   |  2  | 3  |  4  |  5   |  6  |  7 |  8  | 9  |
    +======+=====+====+=====+======+=====+====+=====+====+
    | BSET | ID1 | C1 | ID2 |  C2  | ID3 | C3 | ID4 | C4 |
    +------+-----+----+-----+------+-----+----+-----+----+
    | BSET | 16  |  2 |  23 | 3516 |  1  | 4  |     |    |
    +------+-----+----+-----+------+-----+----+-----+----+
    """
    type = 'BSET'

    def __init__(self, ids, components, comment=''):
        ABCQSet.__init__(self, ids, components, comment)


class CSET(ABCQSet):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.

    +------+-----+----+-----+------+-----+----+-----+----+
    |  1   |  2  | 3  |  4  |  5   |  6  |  7 |  8  | 9  |
    +======+=====+====+=====+======+=====+====+=====+====+
    | CSET | ID1 | C1 | ID2 |  C2  | ID3 | C3 | ID4 | C4 |
    +------+-----+----+-----+------+-----+----+-----+----+
    | CSET | 16  |  2 |  23 | 3516 |  1  | 4  |     |    |
    +------+-----+----+-----+------+-----+----+-----+----+
    """
    type = 'CSET'

    def __init__(self, ids, components, comment=''):
        ABCQSet.__init__(self, ids, components, comment)


class QSET(ABCQSet):
    """
    Defines generalized degrees-of-freedom (q-set) to be used for dynamic
    reduction or component mode synthesis.

    +------+-----+----+-----+------+-----+----+-----+----+
    |  1   |  2  | 3  |  4  |  5   |  6  |  7 |  8  | 9  |
    +======+=====+====+=====+======+=====+====+=====+====+
    | QSET | ID1 | C1 | ID2 |  C2  | ID3 | C3 | ID4 | C4 |
    +------+-----+----+-----+------+-----+----+-----+----+
    | QSET | 16  | 2  | 23  | 3516 |  1  |  4 |     |    |
    +------+-----+----+-----+------+-----+----+-----+----+
    """
    type = 'QSET'

    def __init__(self, ids, components, comment=''):
        ABCQSet.__init__(self, ids, components, comment)


class ABQSet1(Set):
    """
    Generic Class ASET1, BSET1, QSET1 cards inherit from.

    Defines degrees-of-freedom in the analysis set (a-set).

    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    |   1   |  2  |  3  |   4  |  5  |  6  |  7  |  8  |  9  |
    +=======+=====+=====+======+=====+=====+=====+=====+=====+
    | xSET1 |  C  | ID1 |  ID2 | ID3 | ID4 | ID5 | ID6 | ID7 |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    |       | ID8 | ID9 |      |     |     |     |     |     |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    | xSET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    """
    type = 'ABQSet1'
    def __init__(self, components, ids, comment=''):
        Set.__init__(self)
        if comment:
            self.comment = comment

        #:  Component number. (Integer zero or blank for scalar points or any
        #:  unique combination of the Integers 1 through 6 for grid points with
        #:  no embedded blanks.)
        self.components = components

        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)


    @classmethod
    def add_card(cls, card, comment=''):
        components = fcomponents_or_blank(card, 1, 'components', 0)

        nfields = len(card)
        ids = []
        i = 1
        for ifield in range(2, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        return cls(components, ids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        components = str(data[0])
        thru_flag = data[1]
        if thru_flag == 0:
            ids = data[2:]
        elif thru_flag == 1:
            assert len(data) == 4, data
            #ids = [data[2], 'THRU', data[3]]
            ids = list(range(data[2], data[3]+1))
        else:
            raise NotImplementedError('thru_flag=%s data=%s' % (thru_flag, data))
        return cls(components, ids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s' % self.type
        self.ids = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.ids_ref = self.ids

    def uncross_reference(self):
        self.ids = self.node_ids
        del self.ids_ref

    #@property
    #def node_ids(self):
        #return self.get_ids()

    @property
    def node_ids(self):
        msg = ' which is required by %s' % self.type
        return _node_ids(self, self.ids, allow_empty_nodes=True, msg=msg)

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = [self.type, self.components] + collapse_thru(self.node_ids)
        return list_fields

    def __repr__(self):
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)


class SuperABQSet1(Set):
    """
    Generic Class SEBSET1, SEQSET1 cards inherit from.

    Defines degrees-of-freedom in the analysis set (a-set).

    +----------+------+-----+------+------+-----+-----+-----+-----+
    |    1     |  2   |  3  |   4  |   5  |  6  |  7  |  8  |  9  |
    +==========+======+=====+======+======+=====+=====+=====+=====+
    | SEBSET1  | SEID |  C  | ID1  | ID2  | ID3 | ID4 | ID5 | ID6 |
    +----------+------+-----+------+------+-----+-----+-----+-----+
    |          | ID7  | ID9 |      |      |     |     |     |     |
    +----------+------+-----+------+------+-----+-----+-----+-----+
    | SEBSET1  | SEID |  C  | ID1  | THRU | ID2 |     |     |     |
    +----------+------+-----+------+------+-----+-----+-----+-----+
    """
    type = 'SuperABQSet1'
    def __init__(self, seid, components, ids, comment=''):
        Set.__init__(self)
        if comment:
            self.comment = comment
        self.seid = seid

        #:  Component number. (Integer zero or blank for scalar points or any
        #:  unique combination of the Integers 1 through 6 for grid points with
        #:  no embedded blanks.)
        self.components = components

        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)
        #print('ids =', self.ids)
        assert None not in self.ids

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        components = fcomponents_or_blank(card, 2, 'components', 0)

        nfields = len(card)
        ids = []
        i = 1
        for ifield in range(3, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        ids = expand_thru(ids)
        return cls(seid, components, ids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        seid, components, nids = data
        #assert None not in components, 'Type=%s components=%s' % (self.type, components)
        assert None not in nids, 'Type=%s nids=%s' % (self.type, nids)
        assert -1 not in nids, 'nids=%s' % (nids.tolist())
        assert 0 not in nids, 'nids=%s' % (nids.tolist())
        return cls(seid, components, nids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s seid=%s' % (self.type, self.seid)
        self.ids = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.ids_ref = self.ids

    def uncross_reference(self):
        self.ids = self.node_ids
        del self.ids_ref

    @property
    def node_ids(self):
        msg = ' which is required by %s seid=%s' % (self.type, self.seid)
        return _node_ids(self, self.ids, allow_empty_nodes=True, msg=msg)

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = [self.type, self.seid, self.components] + collapse_thru(self.node_ids)
        return list_fields

    def __repr__(self):
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)


class ASET1(ABQSet1):
    """
    Defines degrees-of-freedom in the analysis set (a-set)

    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    |   1   |  2  |  3  |   4  |  5  |  6  |  7  |  8  |  9  |
    +=======+=====+=====+==================+=====+=====+=====+
    | ASET1 |  C  | ID1 |  ID2 | ID3 | ID4 | ID5 | ID6 | ID7 |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    |       | ID8 | ID9 |      |     |     |     |     |     |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    | ASET1 |  C  | ID1 | THRU | ID2 |     |     |     |     |
    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    """
    type = 'ASET1'

    def __init__(self, components, ids, comment=''):
        ABQSet1.__init__(self, components, ids, comment)


class BSET1(ABQSet1):
    type = 'BSET1'

    def __init__(self, components, ids, comment=''):
        ABQSet1.__init__(self, components, ids, comment)


class CSET1(Set):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.

    +-------+-----+-----+------+-----+-----+-----+-----+-----+
    |   1   |  2  |  3  |   4  |  5  |  6  |  7  |  8  |  9  |
    +=======+=====+=====+==================+=====+=====+=====+
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

    def __init__(self, ids, components, comment=''):
        Set.__init__(self)
        if comment:
            self.comment = comment
        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)
        self.components = components

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CSET1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        if integer_string_or_blank(card, 2, 'C') == 'ALL':
            components = '123456'
        else:
            components = parse_components(card, 1, 'components')

        ids = []
        id_count = 1
        for ifield in range(2, len(card)):
            idi = integer_or_string(card, ifield, 'ID%i' % id_count)
            ids.append(idi)
            id_count += 1
        return CSET1(ids, components, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CSET1'
        self.ids = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.ids_ref = self.ids

    def uncross_reference(self):
        self.ids = self.node_ids
        del self.ids_ref

    @property
    def node_ids(self):
        msg = ' which is required by CSET1'
        return _node_ids(self, self.ids, allow_empty_nodes=True, msg=msg)

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['CSET1', self.components] + collapse_thru(self.node_ids)
        return list_fields

    def __repr__(self):
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)


class QSET1(ABQSet1):
    """
    Defines generalized degrees-of-freedom (q-set) to be used for dynamic
    reduction or component mode synthesis.
    """
    type = 'QSET1'

    def __init__(self, components, ids, comment=''):
        ABQSet1.__init__(self, components, ids, comment)


class SET1(Set):
    """
    Defines a list of structural grid points or element identification
    numbers.

    +------+--------+--------+-----+------+-----+-----+------+-----+
    |  1   |    2   |    3   |  4  |   5  |  6  |  7  |   8  |  9  |
    +======+========+========+=====+======+=====+=====+======+=====+
    | SET1 |  SID   |   ID1  | ID2 | ID3  | ID4 | ID5 | ID6  | ID7 |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    |      |  ID8   | -etc.- |     |      |     |     |      |     |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    | SET1 |   3    |   31   | 62  |  93  | 124 | 16  |  17  | 18  |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    |      |   19   |        |     |      |     |     |      |     |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    | SET1 |   6    |   29   | 32  | THRU | 50  | 61  | THRU | 70  |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    |      |   17   |  57    |     |      |     |     |      |     |
    +------+--------+--------+-----+------+-----+-----+------+-----+
    """
    type = 'SET1'

    def __init__(self, sid, ids, is_skin=False, comment=''):
        Set.__init__(self)
        if comment:
            self.comment = comment
        #:  Unique identification number. (Integer > 0)
        self.sid = sid

        #:  List of structural grid point or element identification numbers.
        #:  (Integer > 0 or 'THRU'; for the 'THRU' option, ID1 < ID2 or 'SKIN';
        #:  in field 3)
        self.ids = expand_thru(ids)
        self.clean_ids()

        self.is_skin = is_skin
        self.xref_type = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SET1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')

        ids = fields(integer_or_string, card, 'ID', i=2, j=len(card))

        is_skin = False
        i = 0
        if len(ids) > 0:
            if isinstance(ids[0], string_types) and ids[0] == 'SKIN':
                is_skin = True
                i += 1
        else:
            assert len(card) > 2, card
        return SET1(sid, ids[i:], is_skin=is_skin, comment=comment)

    #def __eq__(self, set1):
        #assert self.type == set1.type, 'type=%s set1.type=%s' % (self.type, set1.type)
        #self.clean_ids()
        #set1.clean_ids()
        #if self.get_IDs() == set1.get_IDs():
            #return True
        #return False

    def symmetric_difference(self, set1):
        ids1 = set(self.get_ids())
        ids2 = set(set1.get_ids())
        return ids1.symmetric_difference(ids2)

    def add_set(self, set1):
        self.ids += set1.get_ids()
        self.clean_ids()

    def raw_fields(self):
        skin = []
        if self.is_skin:
            skin = ['SKIN']
        ids = self.get_ids()
        return ['SET1', self.sid] + skin + self.get_ids()

    def cross_reference(self, model, xref_type, msg='', allow_empty_nodes=False):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        xref_type : str
            {'Node'}
        allow_empty_nodes : bool; default=False
            do all nodes need to exist?

        SPLINEx, ACMODL, PANEL, AECOMP, XYOUTPUT

        - nodes
          - SPLINEx (all nodes must exist)
          - PANEL (all nodes must exist)
          - XYOUTPUT (missing nodes ignored)
          - AECOMP
          - ACMODL (optional)
        - elements
          - ACMODL (optional)
        """
        msg = ' which is required by SET1 sid=%s%s' % (self.sid, msg)
        if xref_type == 'Node':
            self.ids = model.Nodes(self.get_ids(), msg=msg)
        elif xref_type == 'Point':
            self.ids = model.Points(self.get_ids(), msg=msg)
        else:
            raise NotImplementedError("xref_type=%r and must be ['Node']" % xref_type)
        self.xref_type = xref_type

    def uncross_reference(self):
        if self.xref_type in ['Node', 'Point']:
            self.ids = self.get_ids()
            del self.ids_ref
            self.xref_type = None
        else:
            raise NotImplementedError("xref_type=%r and must be ['Node']" % self.xref_type)

    def get_ids(self):
        if self.xref_type is None:
            ids = self.ids
        elif self.xref_type in ['Node', 'Point']:
            ids = [node if isinstance(node, integer_types) else node.nid
                   for node in self.ids]
        else:
            raise NotImplementedError("xref_type=%r and must be ['Node']" % self.xref_type)
        return ids

    def write_card(self, size=8, is_double=False):
        skin = []
        if self.is_skin:
            skin = ['SKIN']

        # checked in NX 2014 / MSC 2005.1
        return self.comment + print_card_8(['SET1', self.sid] + skin + self.get_ids())

        # I thought this worked in the new MSC Nastran...
        # Doesn't work in NX 2014 / MSC 2005.1 (multiple duplicate sids).
        # It may work with one sid, with singles and doubles on one card.
        #field_packs = []
        #singles, doubles = collapse_thru_packs(self.get_ids())
        #if singles:
            #field_packs.append(['SET1', self.sid] + skin + singles)
        #if doubles:
            #for pack in doubles:
                #field_packs.append(['SET1', self.sid] + skin + pack)

        #msg = []
        #for field_pack in field_packs:
            #msg.append(print_card_8(field_pack))
        #return ''.join(msg)


class SET3(Set):
    """
    Defines a list of grids, elements or points.

    SET3 entries are referenced by:
    - NX
      - ACMODL
      - PANEL
    - MSC
      - PBMSECT
      - PBRSECT
      - RFORCE
        - ELEM only (SOL 600)
      - DEACTEL
        - ELEM only (SOL 400)
      - RBAR, RBAR1, RBE1, RBE2, RBE2GS, RBE3, RROD,
        RSPLINE, RSSCON, RTRPLT and RTRPLT1
         - RBEin / RBEex only
      - ELSIDi / XELSIDi
         - ELEM only
      - NDSIDi
         - GRID only

    +------+-----+-------+-----+-----+-----+-----+-----+-----+
    |   1  |  2  |   3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +======+=====+=======+=====+=====+=====+=====+=====+=====+
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
    valid_descs = ['GRID', 'POINT', 'ELEMENT', 'PROP', 'RBEIN', 'RBEEX']

    def __init__(self, sid, desc, ids, comment=''):
        Set.__init__(self)
        if comment:
            self.comment = comment
        #:  Unique identification number. (Integer > 0)
        self.sid = sid

        #:  Set description (Character). Valid options are 'GRID', 'ELEM',
        #:  'POINT' and 'PROP'.
        if desc == 'ELEM':
            desc = 'ELEMENT'
        self.desc = desc

        #:  Identifiers of grids points, elements, points or properties.
        #:  (Integer > 0)
        self.ids = expand_thru(ids)
        self.xref_type = None

    def validate(self):
        if self.desc not in self.valid_descs:
            msg = 'desc=%r; valid_descs=[%s]' % (self.desc, ', '.join(self.valid_descs))
            raise ValueError(msg)

    def get_ids(self):
        if self.xref_type is None:
            ids = self.ids
        elif self.xref_type == 'Point':
            # TODO: improve this...
            ids = [point if isinstance(point, integer_types) else point.nid
                   for point in self.ids]
        else:
            raise NotImplementedError("xref_type=%r and must be ['Node']" % self.xref_type)
        return ids

    def cross_reference(self, model, xref_type, msg=''):
        msg = ' which is required by SET3 sid=%s%s' % (self.sid, msg)
        #if xref_type == 'Node':
            #self.ids = model.Nodes(self.get_ids(), msg=msg)
        if xref_type == 'Point':
            self.ids = model.Points(self.get_ids(), msg=msg)
        else:
            raise NotImplementedError("xref_type=%r and must be ['Point']" % xref_type)
        self.xref_type = xref_type

    def add_set(self, set3):
        self.ids += set3.get_ids()
        assert self.sid == set3.sid, 'SET3.sid=%r; existing sid=%r new=%r' % (self.sid, self.sid, set3.sid)
        assert self.desc == set3.desc, 'SET3.sid=%r; existing desc=%r new=%r' % (self.sid, self.desc, set3.desc)
        self.clean_ids()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SET3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        desc = string(card, 2, 'desc')
        ids = fields(integer_or_string, card, 'ID', i=3, j=len(card))
        return SET3(sid, desc, ids, comment=comment)

    def union(self, set3):
        assert self.type == set3.type, 'type=%r set3.type=%r' % (self.type, set3.type)
        assert self.desc == set3.desc, 'self.desc=%r set3.desc=%r' % (self.desc, set3.desc)
        ids1 = set(self.ids)
        ids2 = set(set3.ids)
        self.ids = list(ids1.union(ids2))

    def symmetric_difference(self, set3):
        assert self.type == set3.type, 'type=%r set3.type=%r' % (self.type, set3.type)
        ids1 = set(self.ids)
        ids2 = set(set3.ids)
        return ids1.symmetric_difference(ids2)

    def is_grid(self):
        if self.desc == 'GRID':
            return True
        return False

    def is_point(self):
        if self.desc == 'POINT':
            return True
        return False

    def is_property(self):
        if self.desc == 'PROP':
            return True
        return False

    def is_element(self):
        if self.desc == 'ELEMENT':
            return True
        return False

    def SetIDs(self, collapse=True):
        """gets the IDs of the SETx"""
        if collapse:
            return collapse_thru(self.ids, nthru=1)
        else:
            return self.ids

    def raw_fields(self):
        """Gets the "raw" card without any processing as a list for printing"""
        list_fields = ['SET3', self.sid, self.desc] + self.SetIDs()
        return list_fields

    def __repr__(self):
        #fields_blocks = [
            #'SET3',
            #[[self.sid, self.desc], False], # these are not all integers
            #[self.SetIDs(), True], # these are all integers
        #]
        #print(fields_blocks)
        #return self.comment + print_int_card_blocks(fields_blocks)
        msg = self.comment

        self.ids.sort()
        ids = self.get_ids()
        packs = condense(ids)
        if len(packs) == 1:
            singles, doubles = build_thru_packs(packs, max_dv=1)

            packs = collapse_thru(ids)
            for pack in doubles:
                msg += print_card_8(['SET3', self.sid, self.desc] + pack)
            if singles:
                msg += print_card_8(['SET3', self.sid, self.desc] + singles)
        else:
            msg += print_card_8(['SET3', self.sid, self.desc] + ids)
        return msg

    def write_card(self, size=8, is_double=False):
        return str(self)

class SESET(SetSuper):
    """
    Defines interior grid points for a superelement.
    """
    type = 'SESET'

    def __init__(self, seid, ids, comment=''):
        SetSuper.__init__(self)
        if comment:
            self.comment = comment
        self.seid = seid
        #:  Grid or scalar point identification number.
        #:  (0 < Integer < 1000000; G1 < G2)
        self.ids = expand_thru(ids)
        self.clean_ids()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SESET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        seid = integer_or_blank(card, 1, 'seid', 0)
        ids = fields(integer_or_string, card, 'ID', i=2, j=len(card))
        return SESET(seid, ids, comment=comment)

    def add_seset(self, seset):
        self.ids += seset.ids
        self.clean_ids()

    def raw_fields(self):
        list_fields = ['SESET', self.seid] + collapse_thru(self.ids)
        return list_fields

    def __repr__(self):
        thru_fields = collapse_thru(self.ids)
        #list_fields = ['SESET', self.seid]

        cards = []
        while 'THRU' in thru_fields:
            ithru = thru_fields.index('THRU')
            card = print_card_8(['SESET', self.seid] +
                                thru_fields[ithru - 1:ithru + 2])
            cards.append(card)
            thru_fields = thru_fields[0:ithru - 1]+thru_fields[ithru + 2:]
        if thru_fields:
            card = print_card_8(['SESET', self.seid] + thru_fields)
            cards.append(card)
        return ''.join(cards)

    def cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass


class SEBSET(SuperABCQSet):
    """
    Defines boundary degrees-of-freedom to be fixed (b-set) during generalized
    dynamic reduction or component mode calculations.

    +--------+------+-----+------+-----+----+-----+----+
    |    1   |   2  |  3  |   4  |  5  |  6 |  7  |  8 |
    +========+======+=====+======+=====+====+=====+====+
    | SEBSET | SEID | ID1 |  C1  | ID2 | C2 | ID3 | C3 |
    +--------+------+-----+------+-----+----+-----+----+
    | SEBSET |  C   | ID1 | THRU | ID2 |    |     |    |
    +--------+------+-----+------+-----+----+-----+----+
    """
    type = 'SEBSET'

    def __init__(self, seid, ids, components, comment=''):
        SuperABCQSet.__init__(self, seid, ids, components, comment)

class SEBSET1(SuperABQSet1):
    type = 'SEBSET1'

    def __init__(self, seid, components, ids, comment=''):
        SuperABQSet1.__init__(self, seid, components, ids, comment)


class SECSET(SuperABCQSet):
    type = 'SECSET'

    def __init__(self, seid, components, ids, comment=''):
        SuperABCQSet.__init__(self, seid, components, ids, comment)

class SECSET1(SuperABQSet1):
    type = 'SECSET1'

    def __init__(self, seid, components, ids, comment=''):
        SuperABQSet1.__init__(self, seid, components, ids, comment)


class SEQSET(SuperABCQSet):
    type = 'SEQSET'

    def __init__(self, seid, ids, components, comment=''):
        SuperABCQSet.__init__(self, seid, ids, components, comment)

class SEQSET1(SuperABQSet1):
    type = 'SEQSET1'

    def __init__(self, seid, components, ids, comment=''):
        SuperABQSet1.__init__(self, seid, components, ids, comment)


class SEQSEP(SetSuper):  # not integrated...is this an SESET ???
    """
    Used with the CSUPER entry to define the correspondence of the
    exterior grid points between an identical or mirror-image superelement
    and its primary superelement.
    """
    type = 'SEQSEP'

    def __init__(self, ssid, psid, ids, comment=''):
        SetSuper.__init__(self)
        if comment:
            self.comment = comment
        #: Identification number for secondary superelement. (Integer >= 0).
        self.ssid = ssid

        #: Identification number for the primary superelement. (Integer >= 0).
        self.psid = psid

        #: Exterior grid point identification numbers for the primary
        #: superelement. (Integer > 0)
        self.ids = expand_thru(ids)
        self.clean_ids()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SEQSEP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        ssid = integer(card, 1, 'ssid')
        psid = integer(card, 2, 'psid')
        ids = fields(integer_or_string, card, 'ID', i=3, j=len(card))
        return SEQSEP(ssid, psid, ids, comment=comment)

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['SEQSEP', self.ssid, self.psid] + self.get_ids()
        return list_fields


class RADSET(Set):  # not integrated
    """
    Specifies which radiation cavities are to be included for
    radiation enclosure analysis.

    +--------+----------+----------+----------+----------+----------+----------+----------+
    |   1    |     2    |     3    |     4    |     5    |     6    |     7    |     8    |
    +========+==========+==========+==========+==========+==========+==========+==========+
    | RADSET | ICAVITY1 | ICAVITY2 | ICAVITY3 | ICAVITY4 | ICAVITY5 | ICAVITY6 | ICAVITY7 |
    +--------+----------+----------+----------+----------+----------+----------+----------+
    |        | ICAVITY8 | ICAVITY9 | -etc.-   |          |          |          |          |
    +--------+----------+----------+----------+----------+----------+----------+----------+
    """
    type = 'RADSET'

    def __init__(self, seid, ids, comment=''):
        Set.__init__(self)
        if comment:
            self.comment = comment
        self.seid = seid
        #: Grid or scalar point identification number.
        #: (0 < Integer < 1000000; G1 < G2)
        self.ids = expand_thru(ids)
        self.clean_ids()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RADSET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        seid = integer(card, 1, 'seid')
        ids = fields(integer_or_string, card, 'ID', i=2, j=len(card))
        return RADSET(seid, ids, comment=comment)

    def add_radset(self, radset):
        self.ids += radset.ids
        self.clean_ids()

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['RADSET', self.seid] + self.get_ids()
        return list_fields


class USET(Set):
    """
    Defines a degrees-of-freedom set.

    +------+-------+-----+------+-----+----+-----+----+
    |  1   |   2   |  3  |   4  |  5  |  6  |  7 | 8  |
    +======+=======+=====+======+=====+====+=====+====+
    | USET | SNAME | ID1 |  C1  | ID2 | C2 | ID3 | C3 |
    +------+-------+-----+------+-----+----+-----+----+
    | USET |  JUNK | ID1 | THRU | ID2 |    |     |    |
    +------+-------+-----+------+-----+----+-----+----+
    """
    type = 'USET'
    def __init__(self, name, components, ids, comment=''):
        Set.__init__(self)
        if comment:
            self.comment = comment
        self.name = name
        #:  Identifiers of grids points. (Integer > 0)
        self.components = components
        self.ids = ids

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a USET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        name = string(card, 1, 'name')
        components = []
        ids = []

        nsets = (len(card) - 1) // 2
        for iset in range(nsets):
            i = iset * 2 + 2
            idi = integer(card, i, 'node_id' + str(iset))
            component = parse_components(card, i + 1, 'component' + str(iset))
            components.append(component)
            ids.append(idi)
        return USET(name, components, ids, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        tested by gspc1.op2

        for some reason, the setname is an integer and has bizarre rules
        that I don't understand like:
          - the setname is 1-4 characters, except if it's 'ZERO%i' % sid
            ummm...odd
        """
        sid = data[0]
        nid = data[1]
        if sid < 0:
            name = 'ZERO'
        else:
            comment = 'sid=%s (???)' % sid
            name = 'U%i' % nid
        assert nid > 0, nid
        component = str(data[2])
        for componenti in component:
            assert componenti in '0123456', component
        return USET(name, [component], [nid], comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by USET name=%s' % (self.name)
        self.ids = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.ids_ref = self.ids

    def uncross_reference(self):
        self.ids = self.node_ids
        del self.ids_ref

    @property
    def node_ids(self):
        msg = ' which is required by USET name=%s' % (self.name)
        return _node_ids(self, self.ids, allow_empty_nodes=True, msg=msg)

    def raw_fields(self):
        """
        gets the "raw" card without any processing as a list for printing
        """
        list_fields = ['USET', self.name]
        for (component, idi) in zip(self.components, self.node_ids):
            list_fields += [idi, component]
        return list_fields

class USET1(ABQSet1):
    """
    Defines degrees-of-freedom in the analysis set (a-set)

    +-------+-------+-----+------+------+-----+-----+-----+-----+
    |   1   |   2   |  3  |  4   |  5   |  6  |  7  |  8  |  9  |
    +=======+=======+=====+======+======+=====+=====+=====+=====+
    | USET1 | SNAME |  C  |  ID2 | ID3  | ID4 | ID5 | ID6 | ID7 |
    +-------+-------+-----+------+------+-----+-----+-----+-----+
    |       | ID9   |     |      |      |     |     |     |     |
    +-------+-------+-----+------+------+-----+-----+-----+-----+
    | USET1 | SNAME |  C  | ID1  | THRU | ID2 |     |     |     |
    +-------+-------+-----+------+------+-----+-----+-----+-----+
    """
    type = 'USET1'

    def __init__(self, name, components, ids, comment=''):
        Set.__init__(self)
        if comment:
            self.comment = comment
        self.name = name

        #:  Component number. (Integer zero or blank for scalar points or any
        #:  unique combination of the Integers 1 through 6 for grid points with
        #:  no embedded blanks.)
        self.components = components

        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a USET1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        name = string(card, 1, 'name')
        components = fcomponents_or_blank(card, 2, 'components', 0)

        nfields = len(card)
        ids = []
        i = 1
        for ifield in range(3, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        return USET1(name, components, ids, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by USET1 name=%s' % (self.name)
        self.ids = model.Nodes(self.node_ids, allow_empty_nodes=True, msg=msg)
        self.ids_ref = self.ids

    def uncross_reference(self):
        self.ids = self.node_ids
        del self.ids_ref

    @property
    def node_ids(self):
        msg = ' which is required by USET1 name=%s' % (self.name)
        return _node_ids(self, self.ids, allow_empty_nodes=True, msg=msg)

    def raw_fields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['USET1', self.name, self.components] + collapse_thru(self.node_ids)
        return list_fields

    def __repr__(self):
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)
