"""
All superelements are defined in this file.  This includes:
 * CSUPER
 * CSUPEXT
 * SEBNDRY
 * SEBULK
 * SECONCT
 * SEELT
 * SEEXCLD
 * SELABEL
 * SELOC
 * SELOAD
 * SEMPLN
 * SENQSET
 * SETREE
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.base_card import (
    BaseCard, expand_thru
)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_or_string,
    string, string_or_blank, double_or_blank, integer_string_or_blank,
    exact_string_or_blank,
)


class SEBNDRY(BaseCard):
    """
    Superelement Boundary-Point Definition

    Defines a list of grid points in a partitioned superelement for the automatic boundary
    search between a specified superelement or between all other superelements in the
    model.

    SEBNDRY SEIDA SEIDB GIDA1 GIDA2 GIDA3 GIDA4 GIDA5 GIDA6
            GIDA7 GIDA8 -etc.-
    SEBNDRY 400  4  10 20 30 40
    SEBNDRY 400 ALL 10 20 30 THRU 40
    """
    type = 'SEBNDRY'
    def __init__(self, seida, seidb, ids, comment=''):
        #BaseCard.__init__()
        if comment:
            self.comment = comment
        self.seida = seida
        self.seidb = seidb

        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)
        self.ids_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid_a = integer(card, 1, 'seid_a')
        seid_b = integer_string_or_blank(card, 2, 'seid_b')
        ids = []
        i = 1
        nfields = len(card)
        for ifield in range(3, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        assert len(card) >= 3, 'len(SEBNDRY card) = %i\ncard=%s' % (len(card), card)
        return SEBNDRY(seid_a, seid_b, ids, comment=comment)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SEELT(BaseCard):
    """
    SEELT SEID EID1 EID2 EID3 EID4 EID5 EID6 EID7
    """
    type = 'SEELT'
    def __init__(self, seid, ids, comment=''):
        #BaseCard.__init__()
        if comment:
            self.comment = comment

        self.seid = seid
        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)
        self.ids_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        ids = []
        i = 1
        nfields = len(card)
        for ifield in range(2, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        assert len(card) <= 9, 'len(SEELT card) = %i\ncard=%s' % (len(card), card)
        return SEELT(seid, ids, comment=comment)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SELOAD(BaseCard):
    """
    External Superelement Load Mapping to Residual

    Maps loads from an external superelement to a specified load set for
    the residual structure.

    SELOAD LIDS0 SEID LIDSE
    SELOAD 10010 100   10
    """
    type = 'SELOC'
    def __init__(self, lid_s0, seid, lid_se, comment=''):
        #BaseCard.__init__()
        if comment:
            self.comment = comment
        self.lid_s0 = lid_s0
        self.seid = seid
        self.lid_se = lid_se
        print(self)

    @classmethod
    def add_card(cls, card, comment=''):
        lid_s0 = integer(card, 1, 'lid_s0')
        seid = integer(card, 2, 'seid')
        lid_se = integer(card, 3, 'lid_se')

        assert len(card) <= 4, 'len(SELOAD card) = %i\ncard=%s' % (len(card), card)
        return SELOAD(lid_s0, seid, lid_se, comment=comment)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SEEXCLD(BaseCard):
    """
    Partitioned Superelement Exclusion

    Defines grids that will be excluded during the attachment of a
    partitioned superelement.

    SEEXCLD SEIDA SEIDB GIDA1 GIDA2 GIDA3 GIDA4 GIDA5 GIDA6
            GIDA7 GIDA8 -etc.-
    """
    def __init__(self, seida, seidb, ids, comment=''):
        #BaseCard.__init__()
        if comment:
            self.comment = comment
        self.seida = seida
        self.seidb = seidb

        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)
        self.ids_ref = None
        print(self)

    @classmethod
    def add_card(cls, card, comment=''):
        seid_a = integer(card, 1, 'seid_a')
        seid_b = integer_string_or_blank(card, 2, 'seid_b')
        ids = []
        i = 1
        nfields = len(card)
        for ifield in range(3, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        assert len(card) >= 3, 'len(SEEXCLD card) = %i\ncard=%s' % (len(card), card)
        return SEEXCLD(seid_a, seid_b, ids, comment=comment)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SEMPLN(BaseCard):
    """
    Superelement Mirror Plane
    Defines a mirror plane for mirroring a partitioned superelement.

    SEMPLN SEID PLANE P1 P2  P3
    SEMPLN 110  PLANE 12 45 1125
    """
    type = 'SEMPLN'
    def __init__(self, seid, p1, p2, p3, comment=''):
        #BaseCard.__init__()
        if comment:
            self.comment = comment
        self.seid = seid
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'lid_s0')
        plane = string(card, 2, 'seid')
        p1 = integer(card, 3, 'p1')
        p2 = integer(card, 4, 'p2')
        p3 = integer(card, 5, 'p3')
        assert plane == 'PLANE', plane

        assert len(card) <= 6, 'len(SEMPLN card) = %i\ncard=%s' % (len(card), card)
        return SEMPLN(seid, p1, p2, p3, comment=comment)

    def raw_fields(self):
        list_fields = ['SEMPLN', self.seid, 'PLANE', self.p1, self.p2, self.p3]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SELABEL(BaseCard):
    """
    Superelement Output Label

    Defines a label or name to be printed in the superelement output headings.

    SELABEL SEID LABEL
    SELABEL  10  LEFT REAR FENDER, MODEL XYZ2000
    """
    type = 'SELABEL'
    def __init__(self, seid, label, comment=''):
        #BaseCard.__init__()
        if comment:
            self.comment = comment

        self.seid = seid
        #:  Identifiers of grids points. (Integer > 0)
        self.label = label
        print(self)

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        label = [exact_string_or_blank(card, ifield, 'label', '        ')
                 for ifield in range(2, len(card))]
        print('label =', label)
        return SELABEL(seid, label, comment=comment)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SELOC(BaseCard):
    """
    Partitioned Superelement Location

    Defines a partitioned superelement relocation by listing three noncolinear points in
    the superelement and three corresponding points not belonging to the superelement.

    SELOC SEID PA1 PA2 PA3  PB1 PB2 PB
    SELOC 110  10  100 111 1010 112 30
    """
    type = 'SELOC'
    def __init__(self, seid, ids, comment=''):
        #BaseCard.__init__()
        if comment:
            self.comment = comment

        self.seid = seid
        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)
        self.ids_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        ids = []
        i = 1
        nfields = len(card)
        for ifield in range(2, nfields):
            idi = integer(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        assert len(card) <= 8, 'len(SELOC card) = %i\ncard=%s' % (len(card), card)
        return SELOC(seid, ids, comment=comment)

    def raw_fields(self):
        list_fields = ['SELOC', self.seid] + list(self.ids)

        return list_fields
    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class SETREE(BaseCard):
    """
    Superelement Tree Definition (Alternate Form of DTI,SETREE)

    Specifies superelement reduction order.

    SETREE SEID SEUP1 SEUP2 SEUP3 SEUP4 SEUP5 SEUP6 SEUP7
           SEUP8 SEUP9 -etc.-
    SETREE 400 10 20 30 40
    """
    type = 'SETREE'
    def __init__(self, seid, ids, comment=''):
        #BaseCard.__init__()
        if comment:
            self.comment = comment
        self.seid = seid
        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)
        self.ids_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        ids = []
        i = 1
        nfields = len(card)
        for ifield in range(2, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        assert len(card) >= 3, 'len(SETREE card) = %i\ncard=%s' % (len(card), card)
        return SETREE(seid, ids, comment=comment)

    def raw_fields(self):
        list_fields = ['SETREE', self.seid] + list(self.ids)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class CSUPER(BaseCard):
    """
    Secondary Superelement Connection

    Defines the grid or scalar point connections for identical or mirror image
    superelements or superelements from an external source. These are all known as
    secondary superelements.

    CSUPER SSlD PSID GP1 GP2 GP3 GP4 GP5 GP6
           GP7 GP8 -etc.-
    """
    type = 'CSUPER'
    def __init__(self, seid, psid, ids, comment=''):
        #BaseCard.__init__()
        if comment:
            self.comment = comment

        self.seid = seid
        self.psid = psid
        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)
        self.ids_ref = None
        print(self)

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        psid = integer(card, 2, 'psid')
        ids = []
        i = 1
        nfields = len(card)
        for ifield in range(3, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        assert len(card) >= 3, 'len(CSUPER card) = %i\ncard=%s' % (len(card), card)
        return CSUPER(seid, psid, ids, comment=comment)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CSUPEXT(BaseCard):
    """
    Superelement Exterior Point Definition

    Assigns exterior points to a superelement.

    CSUPEXT SEID GP1 GP2 GP3 GP4 GP5 GP6 GP7
    """
    type = 'CSUPEXT'
    def __init__(self, seid, ids, comment=''):
        #BaseCard.__init__()
        if comment:
            self.comment = comment

        self.seid = seid
        #:  Identifiers of grids points. (Integer > 0)
        self.ids = expand_thru(ids)
        self.ids_ref = None
        print(self)

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        ids = []
        i = 1
        nfields = len(card)
        for ifield in range(2, nfields):
            idi = integer_string_or_blank(card, ifield, 'ID%i' % i)
            if idi:
                i += 1
                ids.append(idi)
        assert len(card) <= 9, 'len(CSUPEXT card) = %i\ncard=%s' % (len(card), card)
        return CSUPEXT(seid, ids, comment=comment)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SEBULK(BaseCard):
    """
    Partitional Superelement Connection

    Defines superelement boundary search options and a repeated,
    mirrored, or collector superelement.

    | SEBULK | SEID |  TYPE  | RSEID | METHOD |  TOL   | LOC | UNITNO
    | SEBULK | 14   | REPEAT |   4   |  AUTO  | 1.0E-3 |     |
    """
    type = 'SECONCT'
    def __init__(self, seid, Type, rseid,
                 method='AUTO', tol=1e-5, loc='YES', unitno=None,
                 comment=''):
        """
        Parameters
        ----------
        seid : int
            Partitioned superelement identification number.
        Type : str
            Superelement type.
            {PRIMARY, REPEAT, MIRROR, COLLCTR, EXTERNAL, EXTOP2}
        rseid : int; default=0
            Identification number of the reference superelement,
            used if TYPE = 'REPEAT' and 'MIRROR'.
        method : str; default='AUTO'
            Method to be used when searching for boundary grid points.
            {AUTO, MANUAL}
        tol : float; default=1e-5
            Location tolerance to be used when searching for boundary grid points.
        loc : str; default='YES'
            Coincident location check option for manual conection option.
            {YES, NO}
        unitno : int / None
            FORTRAN unit number for the OUTPUT2 file (applicable and
            meaningful only when TYPE='EXTOP2').

        """
        #BaseCard.__init__()
        if comment:
            self.comment = comment
        self.seid = seid
        self.Type = Type
        self.rseid = rseid
        self.method = method
        self.tol = tol
        self.loc = loc
        self.unitno = unitno

    @classmethod
    def add_card(cls, card, comment=''):
        seid = integer(card, 1, 'seid')
        Type = string(card, 2, 'Type')
        rseid = integer_or_blank(card, 3, 'rseid', 0)
        method = string_or_blank(card, 4, 'method', 'AUTO')
        tol = double_or_blank(card, 5, 'tol', 1e-5)
        loc = string_or_blank(card, 6, 'loc', 'YES')
        unitno = integer_or_blank(card, 7, 'seid')
        assert len(card) <= 8, 'len(SENQSET card) = %i\ncard=%s' % (len(card), card)
        return SEBULK(seid, Type, rseid, method=method, tol=tol,
                      loc=loc, unitno=unitno, comment=comment)

    def raw_fields(self):
        list_fields = ['SEBULK', self.seid, self.Type, self.rseid, self.method, self.tol,
                       self.loc, self.unitno]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SECONCT(BaseCard):
    """
    Partitioned Superelement Boundary-Point Connection

    Explicitly defines grid and scalar point connection procedures for a
    partitioned superelement.

    SECONCT SEIDA SEIDB TOL LOC
            GIDA1 GIDB1 GIDA2 GIDB2 GIDA3 GIDB3 -etc.-
    SECONCT 10 20 1.0E-4 YES
            1001  4001              2222  4444
    SECONCT SEIDA SEIDB TOL LOC
            GIDA1 THRU GIDA2 GIDB1 THRU GIDB2
    SECONCT 10  20
            101 THRU 110 201 THRU 210
    """
    type = 'SECONCT'
    def __init__(self, set_id, n, comment=''):
        """
        Parameters
        ----------
        SEIDA : int
            Partitioned superelement identification number.
        SEIDB : int
            Identification number of superelement for connection to SEIDA.
        TOL : float; default=1e-5
            Location tolerance to be used when searching for or checking boundary
            grid points.
        LOC : str; default='YES'
            Coincident location check option for manual connection.
            {YES, NO}
        GIDAi : int
            Identification number of a grid or scalar point in superelement SEIDA,
            which will be connected to GIDBi.
        GIDBi : int
            Identification number of a grid or scalar point in superelement SEIDB,
            which will be connected to GIDAi.

        """
        #BaseCard.__init__()
        if comment:
            self.comment = comment
        self.set_id = set_id
        self.n = n
        print(self)

    @classmethod
    def add_card(cls, card, comment=''):
        raise NotImplementedError('SECONCT')
        set_id = integer_or_string(card, 1, 'set_id')
        n = integer_or_blank(card, 2, 'n', 0)
        assert len(card) <= 3, 'len(SENQSET card) = %i\ncard=%s' % (len(card), card)
        return SENQSET(set_id, n, comment=comment)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SENQSET(BaseCard):
    """
    Superelement Internal Generalized Degree-of-Freedom

    Defines number of internally generated scalar points for superelement dynamic
    reduction.

    SENQSET SEID N
    SENQSET 110 45
    """
    type = 'CSUPER'
    def __init__(self, set_id, n, comment=''):
        """
        Parameters
        ----------
        set_id : int / str
            Partitioned superelement identification number.
            (Integer > 0 or Character='ALL')
        n : int; default=0
            Number of internally generated scalar points for dynamic
            reduction generalized coordinates (Integer > 0).
        """
        #BaseCard.__init__()
        if comment:
            self.comment = comment
        self.set_id = set_id
        self.n = n
        print(self)

    @classmethod
    def add_card(cls, card, comment=''):
        set_id = integer_or_string(card, 1, 'set_id')
        n = integer_or_blank(card, 2, 'n', 0)
        assert len(card) <= 3, 'len(SENQSET card) = %i\ncard=%s' % (len(card), card)
        return SENQSET(set_id, n, comment=comment)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)
