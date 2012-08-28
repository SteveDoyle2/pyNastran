# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import izip

from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru, collapse_thru
from pyNastran.bdf.fieldWriter import print_int_card
from pyNastran.bdf.fieldWriter import printCard as print_card


class Set(BaseCard):
    """Generic Class all SETx cards inherit from"""
    type = 'Set'

    def __init__(self, card, data):
        ## Unique identification number. (Integer > 0)
        self.sid = None
        ## list of IDs in the SETx
        self.IDs = None

    def cleanIDs(self):
        """eliminates duplicate IDs from self.IDs and sorts self.IDs"""
        self.IDs = list(set(self.IDs))
        self.IDs.sort()

    def SetIDs(self):
        """gets the IDs of the SETx"""
        return collapse_thru(self.IDs)

    def reprFields(self):
        fields = self.rawFields()
        return fields

    def __repr__(self):
        return print_card(self.reprFields())


class SetSuper(Set):
    """Generic Class all Superelement SETx cards inherit from."""
    def __init__(self, card, data):
        ## Superelement identification number. Must be a primary superelement.
        ## (Integer >= 0)
        self.seid = None
        ## list of IDs in the SESETx
        self.IDs = None


class ABCQSet(Set):
    """
    Generic Class ASET, BSET, CSET, QSET cards inherit from.

    Defines degrees-of-freedom in the analysis set (a-set)
    ASET ID1 C1 ID2 C2   ID3 C3 ID4 C4
    ASET 16  2  23  3516 1   4
    """
    def __init__(self, card=None, data=None): ## @todo doesnt support data
        Set.__init__(self, card, data)

        ## Identifiers of grids points. (Integer > 0)
        self.IDs = []
        self.components = []

        fields = card.fields(1)
        for i in xrange(0, fields, 2):
            self.IDs.append(fields[i])
            self.components.append(str(fields[i + 1]))

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        fields = [self.type]  # ASET, BSET
        for (ID, comp) in izip(self.IDs, self.components):
            fields += [ID, comp]
        return fields

    def __repr__(self):
        fields = self.rawFields()
        return self.printCard(fields)


class ASET(ABCQSet):
    """
    Defines degrees-of-freedom in the analysis set (a-set).
    ASET ID1 C1 ID2 C2   ID3 C3 ID4 C4
    ASET 16  2  23  3516 1   4
    """
    type = 'ASET'

    def __init__(self, card=None, data=None):  ## @todo doesnt support data
        ABCQSet.__init__(self, card, data)


class BSET(ABCQSet):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.
    ASET ID1 C1 ID2 C2   ID3 C3 ID4 C4
    ASET 16  2  23  3516 1   4
    """
    type = 'BSET'

    def __init__(self, card=None, data=None):  ## @todo doesnt support data
        ABCQSet.__init__(self, card, data)


class CSET(ABCQSet):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.
    ASET ID1 C1 ID2 C2   ID3 C3 ID4 C4
    ASET 16  2  23  3516 1   4
    """
    type = 'CSET'

    def __init__(self, card=None, data=None):  ## @todo doesnt support data
        ABCQSet.__init__(self, card, data)


class QSET(ABCQSet):
    """
    Defines generalized degrees-of-freedom (q-set) to be used for dynamic
    reduction or component mode synthesis.
    QSET ID1 C1 ID2 C2   ID3 C3 ID4 C4
    QSET 16  2  23  3516 1   4
    """
    type = 'QSET'

    def __init__(self, card=None, data=None):  ## @todo doesnt support data
        ABCQSet.__init__(self, card, data)


class ABQSet1(Set):
    """
    Generic Class ASET1, BSET1, QSET1 cards inherit from.

    Defines degrees-of-freedom in the analysis set (a-set).
    ASET1 C ID1 ID2 ID3 ID4 ID5 ID6 ID7
    ID8 ID9
    ASET1 C ID1 'THRU' ID2
    """
    def __init__(self, card=None, data=None):  ## @todo doesnt support data
        Set.__init__(self, card, data)
        ## Component number. (Integer zero or blank for scalar points or any
        ## unique combination of the Integers 1 through 6 for grid points with
        ## no embedded blanks.)
        self.components = str(card.field(1, 0))

        ## Identifiers of grids points. (Integer > 0)
        self.IDs = expand_thru(card.fields(2))

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        fields = [self.type, self.components] + collapse_thru(self.IDs)
        return fields

    def __repr__(self):
        fields = self.rawFields()
        return self.printCard(fields)


class ASET1(ABQSet1):
    """
    Defines degrees-of-freedom in the analysis set (a-set)
    ASET1 C ID1 ID2 ID3 ID4 ID5 ID6 ID7
    ID8 ID9
    ASET1 C ID1 'THRU' ID2
    """
    type = 'ASET1'

    def __init__(self, card=None, data=None):
        ABQSet1.__init__(self, card, data)


class BSET1(ABQSet1):
    type = 'BSET1'

    def __init__(self, card=None, data=None):
        ABQSet1.__init__(self, card, data)


class CSET1(Set):
    """
    Defines analysis set (a-set) degrees-of-freedom to be fixed (b-set) during
    generalized dynamic reduction or component mode synthesis calculations.
    CSET1 C ID1 ID2 ID3 ID4 ID5 ID6 ID7
    ID8 ID9
    CSET1 C ID1 'THRU' ID2
    CSET1,,'ALL'
    """
    def __init__(self, card=None, data=None):  ## @todo doesnt support data
        Set.__init__(self, card, data)

        ## Identifiers of grids points. (Integer > 0)
        self.IDs = []
        if card.field(2) == 'ALL':
            self.components = '123456'
        else:
            self.components = str(card.field(1))
            self.IDs = expand_thru(card.fields(2))
        ###

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        fields = [self.type, self.components] + collapse_thru(self.IDs)
        return fields

    def __repr__(self):
        fields = self.rawFields()
        return self.printCard(fields)


class QSET1(ABQSet1):
    """
    Defines generalized degrees-of-freedom (q-set) to be used for dynamic
    reduction or component mode synthesis.
    """
    type = 'QSET1'

    def __init__(self, card=None, data=None):
        ABQSet1.__init__(self, card, data)


class SET1(Set):
    """
    Defines a list of structural grid points or element identification numbers.
    SET1 SID ID1 ID2 ID3 ID4 ID5 ID6 ID7
    ID8 -etc.-
    SET1 3 31 62 93 124 16 17 18
    19
    SET1 6 29 32 THRU 50 61 THRU 70
    17 57
    """
    type = 'SET1'

    def __init__(self, card=None, data=None):  ## @todo doesnt support data
        Set.__init__(self, card, data)
        ## Unique identification number. (Integer > 0)
        self.sid = card.field(1)

        ## List of structural grid point or element identification numbers.
        ## (Integer > 0 or 'THRU'; for the 'THRU' option, ID1 < ID2 or 'SKIN';
        ## in field 3)
        self.IDs = []

        fields = card.fields(2)
        self.isSkin = False
        i = 0
        if isinstance(fields[0], str) and fields[0].upper() == 'SKIN':
            self.isSkin = True
            i += 1
        self.IDs = expand_thru(fields[i:])
        self.cleanIDs()

    def IsSkin(self):
        return self.isSkin

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        fields = ['SET1', self.sid]
        if self.isSkin:
            fields.append('SKIN')
        fields += self.SetIDs()
        return fields


class SET3(Set):
    """
    Defines a list of grids, elements or points.
    SET3 SID DES ID1 ID2 ID3 ID4 ID5 ID6
    ID7 ID8 -etc-
    SET3 1 POINT 11 12
    """
    type = 'SET1'

    def __init__(self, card=None, data=None):  ## @todo doesnt support data
        Set.__init__(self, card, data)
        ## Unique identification number. (Integer > 0)
        self.sid = card.field(1)

        ## Set description (Character). Valid options are 'GRID', 'ELEM',
        ## 'POINT' and 'PROP'.
        self.desc = card.field(2).upper()
        assert self.desc in ['GRID', 'POINT', 'ELEM', 'PROP']

        ## Identifiers of grids points, elements, points or properties.
        ## (Integer > 0)
        self.IDs = []

        fields = card.fields(2)
        self.IDs = expand_thru(fields)
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
        """gets the "raw" card without any processing as a list for printing"""
        fields = ['SET3', self.sid, self.desc] + self.SetIDs()
        return fields

    def __repr__(self):
        fields = ['SET3', self.sid, self.desc] + self.SetIDs()
        return self.printCard(fields)


class SESET(SetSuper):
    """
    Defines interior grid points for a superelement.
    """
    type = 'SESET'

    def __init__(self, card=None, data=None):
        SetSuper.__init__(self, card, data)
        self.seid = card.field(1, 0)
        ## Grid or scalar point identification number.
        ## (0 < Integer < 1000000; G1 < G2)
        self.IDs = []

        self.IDs = expand_thru(card.fields(2))
        self.cleanIDs()

    def rawFields(self):
        fields = ['SESET', self.seid] + self.SetIDs()
        return fields


class SEQSEP(SetSuper):  # not integrated...is this an SESET ???
    """
    Used with the CSUPER entry to define the correspondence of the
    exterior grid points between an identical or mirror-image superelement
    and its primary superelement.
    """
    type = 'SEQSEP'

    def __init__(self, card=None, data=None):
        SetSuper.__init__(self, card, data)
        ## Identification number for secondary superelement. (Integer >= 0).
        self.ssid = card.field(1)
        ## Identification number for the primary superelement. (Integer >= 0).
        self.psid = card.field(2)
        ## Exterior grid point identification numbers for the primary
        ## superelement. (Integer > 0)
        self.IDs = []

        fields = card.fields(3)
        self.IDs = expand_thru(fields)
        self.cleanIDs()

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        fields = ['SEQSEP', self.ssid, self.psid] + self.SetIDs()
        return fields


class RADSET(Set):  # not integrated
    """
    Specifies which radiation cavities are to be included for
    radiation enclosure analysis.
    RADSET ICAVITY1 ICAVITY2 ICAVITY3 ICAVITY4 ICAVITY5 ICAVITY6 ICAVITY7
    ICAVITY8 ICAVITY9 -etc.-
    """
    type = 'RADSET'

    def __init__(self, card=None, data=None):
        Set.__init__(self, card, data)
        self.seid = card.field(1)
        ## Grid or scalar point identification number.
        ## (0 < Integer < 1000000; G1 < G2)
        self.IDs = []

        self.IDs = expand_thru(card.fields(2))
        self.cleanIDs()

    def addRadsetObject(self, radset):
        self.IDs += radset.IDs
        self.cleanIDs()

    def rawFields(self):
        """gets the "raw" card without any processing as a list for printing"""
        fields = ['RADSET', self.seid] + self.SetIDs()
        return fields
