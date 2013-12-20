# pylint: disable=C0103,R0902,R0904,R0914,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import izip

from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru, collapse_thru
from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    components, components_or_blank, fields,
    integer_or_string, string, string_or_blank,
    integer_string_or_blank)

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

    def reprFields(self):
        list_fields = self.rawFields()
        return list_fields

    def __repr__(self):
        return print_card_8(self.reprFields())


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

    Defines degrees-of-freedom in the analysis set (a-set)::

      ASET ID1 C1 ID2 C2   ID3 C3 ID4 C4
      ASET 16  2  23  3516 1   4
    """
    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment

        #:  Identifiers of grids points. (Integer > 0)
        self.IDs = []
        self.components = []

        nterms = len(card) // 2
        for n in xrange(nterms):
            i = n * 2 + 1
            ID = integer(card, i, 'ID' + str(n))
            component = components(card, i + 1, 'component' + str(n))
            self.IDs.append(ID)
            self.components.append(component)

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = [self.type]  # ASET, BSET
        for (ID, comp) in izip(self.IDs, self.components):
            list_fields += [ID, comp]
        return list_fields

    def __repr__(self):
        list_fields = self.rawFields()
        return print_card_8(list_fields)


class ASET(ABCQSet):
    """
    Defines degrees-of-freedom in the analysis set (a-set).::

      ASET ID1 C1 ID2 C2   ID3 C3 ID4 C4
      ASET 16  2  23  3516 1   4
    """
    type = 'ASET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)


class BSET(ABCQSet):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.::

      BSET ID1 C1 ID2 C2   ID3 C3 ID4 C4
      BSET 16  2  23  3516 1   4
    """
    type = 'BSET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)


class CSET(ABCQSet):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.::

      CSET ID1 C1 ID2 C2   ID3 C3 ID4 C4
      CSET 16  2  23  3516 1   4
    """
    type = 'CSET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)


class QSET(ABCQSet):
    """
    Defines generalized degrees-of-freedom (q-set) to be used for dynamic
    reduction or component mode synthesis.::

      QSET ID1 C1 ID2 C2   ID3 C3 ID4 C4
      QSET 16  2  23  3516 1   4
    """
    type = 'QSET'

    def __init__(self, card=None, data=None, comment=''):
        ABCQSet.__init__(self, card, data, comment)


class ABQSet1(Set):
    """
    Generic Class ASET1, BSET1, QSET1 cards inherit from.

    Defines degrees-of-freedom in the analysis set (a-set).::

      ASET1 C ID1 ID2 ID3 ID4 ID5 ID6 ID7
      ID8 ID9
      ASET1 C ID1 'THRU' ID2
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

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = [self.type, self.components] + collapse_thru(self.IDs)
        return list_fields

    def __repr__(self):
        list_fields = self.rawFields()
        return print_card_8(list_fields)


class ASET1(ABQSet1):
    """
    Defines degrees-of-freedom in the analysis set (a-set)::

      ASET1 C ID1 ID2 ID3 ID4 ID5 ID6 ID7
      ID8 ID9
      ASET1 C ID1 'THRU' ID2
    """
    type = 'ASET1'

    def __init__(self, card=None, data=None, comment=''):
        ABQSet1.__init__(self, card, data, comment)


class BSET1(ABQSet1):
    type = 'BSET1'

    def __init__(self, card=None, data=None, comment=''):
        ABQSet1.__init__(self, card, data, comment)


class CSET1(Set):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.::

      CSET1 C ID1 ID2 ID3 ID4 ID5 ID6 ID7
      ID8 ID9
      CSET1 C ID1 'THRU' ID2
      CSET1,,'ALL'
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
            IDs = fields(integer_or_string, 'ID', i=2, j=len(card))
            self.IDs = expand_thru(IDs)

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['CSET1', self.components] + collapse_thru(self.IDs)
        return list_fields

    def __repr__(self):
        list_fields = self.rawFields()
        return print_card_8(list_fields)


class QSET1(ABQSet1):
    """
    Defines generalized degrees-of-freedom (q-set) to be used for dynamic
    reduction or component mode synthesis.
    """
    type = 'QSET1'

    def __init__(self, card=None, data=None, comment=''):
        ABQSet1.__init__(self, card, data, comment)


class SET1(Set):
    """
    Defines a list of structural grid points or element identification
    numbers.::

      SET1 SID ID1 ID2 ID3 ID4 ID5 ID6 ID7
      ID8 -etc.-
      SET1 3 31 62 93 124 16 17 18
      19
      SET1 6 29 32 THRU 50 61 THRU 70
      17 57
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

    def IsSkin(self):
        return self.isSkin

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['SET1', self.sid]
        if self.isSkin:
            list_fields.append('SKIN')
        list_fields += self.SetIDs()
        return list_fields


class SET3(Set):
    """
    Defines a list of grids, elements or points.::

      SET3 SID DES ID1 ID2 ID3 ID4 ID5 ID6
      ID7 ID8 -etc-
      SET3 1 POINT 11 12
    """
    type = 'SET1'

    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment
        #:  Unique identification number. (Integer > 0)
        self.sid = integer(card, 1, 'sid')

        #:  Set description (Character). Valid options are 'GRID', 'ELEM',
        #:  'POINT' and 'PROP'.
        self.desc = string(card, 2, 'desc')
        assert self.desc in ['GRID', 'POINT', 'ELEM', 'PROP']

        #:  Identifiers of grids points, elements, points or properties.
        #:  (Integer > 0)
        self.IDs = []

        IDs = fields(integer_or_string, card, 'ID', i=3, j=len(card))
        self.IDs = expand_thru(IDs)
        self.cleanIDs()

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

    def rawFields(self):
        """Gets the "raw" card without any processing as a list for printing"""
        list_fields = ['SET3', self.sid, self.desc] + self.SetIDs()
        return list_fields

    def __repr__(self):
        list_fields = ['SET3', self.sid, self.desc] + self.SetIDs()
        return print_card_8(list_fields)


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

    def rawFields(self):
        return ['SESET', self.seid] + collapse_thru(self.IDs)

    def __repr__(self):
        thruFields = collapse_thru(self.IDs)

        #list_fields = ['SESET', self.seid]

        #i = 0
        cards = []
        #print("thruFields", thruFields)
        while 'THRU' in thruFields:
            #print("thruFields", thruFields)
            iThru = thruFields.index('THRU')
            card = print_card_8(['SESET', self.seid] +
                                thruFields[iThru - 1:iThru + 2])
            cards.append(card)
            thruFields = thruFields[0:iThru - 1]+thruFields[iThru + 2:]
        #print("thruFields", thruFields)
        if thruFields:
            card = print_card_8(['SESET', self.seid] + thruFields)
            cards.append(card)
        #print("cards",cards)
        return ''.join(cards)

class SEBSET(Set):
    """
    Defines boundary degrees-of-freedom to be fixed (b-set) during generalized
    dynamic reduction or component mode calculations.::

      SEBSET SEID ID1 C1 ID2 C2 ID3 C3
      SEBSET C ID1 'THRU' ID2
    """
    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment

        #:  Identifiers of grids points. (Integer > 0)
        self.components = []
        self.IDs = []

        #fields = str(card.fields(1))
        nsets = (len(card) - 1 ) // 2
        for n in range(nsets):
            i = n * 2 + 1
            component = components(card, i, 'component' + str(n))
            ID = components(card, i + 1, 'ID' + str(n))
            self.components.append(component)
            self.IDs.append(ID)

    def rawFields(self):
        """
        gets the "raw" card without any processing as a list for printing
        """
        list_fields = ['SEBSET1']
        for (component, ID) in zip(self.components, self.IDs):
            list_fields += [component, ID]
        return list_fields

class SEBSET1(Set):
    """
    Defines boundary degrees-of-freedom to be fixed (b-set) during generalized
    dynamic reduction or component mode calculations.::

      SEBSET1 C ID1 ID2 ID3 ID4 ID5 ID6 ID7
      ID8 ID9
      SEBSET1 C ID1 'THRU' ID2
    """
    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment

        #:  Identifiers of grids points. (Integer > 0)
        self.IDs = []
        self.components = components(card, 1, 'components')
        IDs = fields(integer_or_string, card, 2, 'ID')
        self.IDs = expand_thru(IDs)

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['SEBSET1', self.components] + collapse_thru(self.IDs)
        return list_fields

    def __repr__(self):
        list_fields = self.rawFields()
        return print_card_8(list_fields)

class SEQSET1(Set):
    """
    Defines the generalized degrees-of-freedom of the superelement to be used in
    generalized dynamic reduction or component mode synthesis.::

      SEQSET1 C ID1 ID2 ID3 ID4 ID5 ID6 ID7
      ID8 ID9
      SEQSET1 C ID1 'THRU' ID2
    """
    def __init__(self, card=None, data=None, comment=''):
        Set.__init__(self, card, data)
        if comment:
            self._comment = comment

        #:  Identifiers of grids points. (Integer > 0)
        self.IDs = []
        self.components = components(card, 1, 'components')
        IDs = fields(integer_or_string, card, i=2, j=len(card))
        self.IDs = expand_thru(IDs)

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['SEQSET1', self.components] + collapse_thru(self.IDs)
        return list_fields

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

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['SEQSEP', self.ssid, self.psid] + self.SetIDs()
        return list_fields


class RADSET(Set):  # not integrated
    """
    Specifies which radiation cavities are to be included for
    radiation enclosure analysis.::

      RADSET ICAVITY1 ICAVITY2 ICAVITY3 ICAVITY4 ICAVITY5 ICAVITY6 ICAVITY7
      ICAVITY8 ICAVITY9 -etc.-
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

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        list_fields = ['RADSET', self.seid] + self.SetIDs()
        return list_fields
