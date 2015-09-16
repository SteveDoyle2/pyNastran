# pylint: disable=R0902,R0904,R0914,C0302,C0111
"""
All aero cards are defined in this file.  This includes:

 * AEFACT
 * AELINK
 * AELIST
 * AEPARM
 * AESTAT
 * AESURF / AESURFS
 * AERO / AEROS
 * CSSCHD
 * CAERO1 / CAERO2 / CAERO3 / CAERO4 / CAERO5
 * FLFACT
 * FLUTTER
 * GUST
 * MKAERO1 / MKAERO2
 * PAERO1 / PAERO2 / PAERO3
 * SPLINE1 / SPLINE2 / SPLINE4 / SPLINE5
 * MNTPNT1

All cards are BaseCard objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import integer_types
from six.moves import zip, range
from itertools import count
from numpy import array, pi, linspace, zeros, arange, repeat, cos, arcsin
from numpy.linalg import norm

from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8, print_float_8
from pyNastran.bdf.cards.baseCard import BaseCard, expand_thru
from pyNastran.bdf.bdfInterface.assign_type import (fields,
    integer, integer_or_blank, double, double_or_blank, string, string_or_blank,
    integer_or_string, double_string_or_blank, blank, interpret_value, components)
from pyNastran.bdf.cards.utils import wipe_empty_fields


class AECOMP(BaseCard):
    """
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    |   1    |   2   |    3     |   4   |   5   |   6   |   7   |   8   |   9   |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    | AECOMP | NAME  | LISTTYPE | LIST1 | LIST2 | LIST3 | LIST4 | LIST5 | LIST6 |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    |        | LIST7 |  -etc.-  |       |       |       |       |       |       |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+

    Attributes
    ----------
    name : str
        The name.
    list_type : str
        {'SET1', 'AELIST', 'CAEROx'}
    lists : list[int]
        lis of values of AECOMP lists
    """
    type = 'AECOMP'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.name = string(card, 1, 'name')
            self.list_type = string(card, 2, 'list_type')
            assert self.list_type in ['SET1', 'AELIST'] or self.list_type.startswith('CAERO'), self.list_type
            j = 1
            self.lists = []
            for i in range(3, len(card)):
                list_i = integer(card, i, '%s_%i' % (self.list_type, j))
                self.lists.append(list_i)
                j += 1
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def cross_reference(self, model):
        msg = ', which is required by AECOMP name=%r' % self.name
        #return
        if self.list_type == 'SET1':
            self.lists = [model.SET1(key, msg) for key in self.lists]
        elif self.list_type == 'AELIST':
            self.lists = [model.AELIST(key, msg) for key in self.lists]
        elif self.list_type == 'CAERO':
            self.lists = [model.CAERO(key, msg) for key in self.lists]
        else:
            raise NotImplementedError(self.list_type)

    def get_lists(self):
        if self.list_type == 'SET1':
            lists = [set1 if isinstance(set1, integer_types)
                          else set1.sid for set1 in self.lists]
        elif self.list_type == 'AELIST':
            lists = [aelist if isinstance(aelist, integer_types)
                          else aelist.sid for aelist in self.lists]
        elif self.list_type == 'CAERO':
            lists = [caero if isinstance(caero, integer_types)
                          else caero.eid for caero in self.lists]
        else:
            raise NotImplementedError(self.list_type)
        return lists

    def raw_fields(self):
        list_fields = ['AECOMP', self.name, self.list_type] + self.get_lists()
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class AEFACT(BaseCard):
    """
    Defines real numbers for aeroelastic analysis.

    +--------+-----+----+--------+-----+----+----+----+----+
    | AEFACT | SID | D1 | D2     | D3  | D4 | D5 | D6 | D7 |
    +--------+-----+----+--------+-----+----+----+----+----+
    |        | D8  | D9 | -etc.- |     |    |    |    |    |
    +--------+-----+----+--------+-----+----+----+----+----+

    +--------+-----+----+--------+-----+
    | AEFACT | 97  |.3  | 0.7    | 1.0 |
    +--------+-----+----+--------+-----+
    """
    type = 'AEFACT'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Set identification number. (Unique Integer > 0)
            self.sid = integer(card, 1, 'sid')

            Di = []
            for i in range(2, len(card)):
                di = double(card, i, 'Di_%i' % (i - 1))
                Di.append(di)

            #: Number (float)
            self.Di = array(Di, dtype='float64')
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the AEFACT object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        fields = ['AEFACT', self.sid] + list(self.Di)
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AELINK(BaseCard):
    r"""
    Defines relationships between or among AESTAT and AESURF entries, such
    that:

    .. math:: u^D + \Sigma_{i=1}^n C_i u_i^I = 0.0

    +--------+-------+-------+--------+----+-------+----+-------+----+
    | AELINK | ID    | LABLD | LABL1  | C1 | LABL2 | C2 | LABL3 | C3 |
    +--------+-------+-------+--------+----+-------+----+-------+----+
    |        | LABL4 | C4    | etc.   |    |       |    |       |    |
    +--------+-------+-------+--------+----+-------+----+-------+----+

    +--------+-------+-------+-------+------+
    | AELINK | 10    | INBDA | OTBDA | -2.0 |
    +--------+-------+-------+-------+------+
    """
    type = 'AELINK'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: an ID=0 is applicable to the global subcase, ID=1 only subcase 1
            self.id = integer(card, 1, 'ID')
            #: defines the dependent variable name (string)
            self.label = string(card, 2, 'label')
            #: defines the independent variable name (string)
            self.independentLabels = []
            #: linking coefficient (real)
            self.Cis = []

            fields = [interpret_value(field) for field in card[3:]]
            assert len(fields) % 2 == 0, 'fields=%s' % fields
            for i in range(0, len(fields), 2):
                independentLabel = fields[i]
                Ci = fields[i + 1]
                self.independentLabels.append(independentLabel)
                self.Cis.append(Ci)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the AELINK object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        fields = ['AELINK', self.id, self.label]
        for (ivar, ival) in zip(self.independentLabels, self.Cis):
            fields += [ivar, ival]
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class AELIST(BaseCard):
    """
    Defines a list of aerodynamic elements to undergo the motion prescribed
    with the AESURF Bulk Data entry for static aeroelasticity.

    +---------+------+------+------+------+------+------+------+------+
    |  AELIST |  SID | E1   | E2   | E3   | E4   | E5   | E6   | E7   |
    +---------+------+------+------+------+------+------+------+------+
    |         |  E8  | etc. |      |      |      |      |      |      |
    +---------+------+------+------+------+------+------+------+------+

    +---------+------+------+------+------+------+------+------+------+
    |  AELIST |  75  | 1001 | THRU | 1075 | 1101 | THRU | 1109 | 1201 |
    +---------+------+------+------+------+------+------+------+------+
    |         | 1202 |      |      |      |      |      |      |      |
    +---------+------+------+------+------+------+------+------+------+

    Remarks
    -------

    1. These entries are referenced by the AESURF entry.
    2. When the THRU option is used, all intermediate grid points must exist.
       The word THRU may not appear in field 3 or 9 (2 or 9 for continuations).
    3. Intervening blank fields are not allowed.
    """
    type = 'AELIST'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Set identification number. (Integer > 0)
            self.sid = integer(card, 1, 'sid')

            #: List of aerodynamic boxes generated by CAERO1 entries to define a
            #: surface. (Integer > 0 or 'THRU')
            eids = fields(integer_or_string, card, 'eid', i=2, j=len(card))
            self.elements = expand_thru(eids)
            self.clean_ids()
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def cross_reference(self, model):
        pass

    def clean_ids(self):
        self.elements = list(set(self.elements))
        self.elements.sort()

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the AELIST object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['AELIST', self.sid] + self.elements
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AEPARM(BaseCard):
    """
    Defines a general aerodynamic trim variable degree-of-freedom (aerodynamic
    extra point). The forces associated with this controller will be derived
    from AEDW, AEFORCE and AEPRESS input data.

    +--------+----+--------+-------+
    | AEPARM | ID | LABEL  | UNITS |
    +--------+----+--------+-------+
    | AEPARM | 5  | THRUST | LBS   |
    +--------+----+--------+-------+
    """
    type = 'AEPARM'
    _field_map = {
        1: 'id', 2:'label', 3:'units'
    }

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.id = integer(card, 1, 'id')
            self.label = string(card, 2, 'lable')
            self.units = string(card, 3, 'units')
            assert len(card) <= 4, 'len(AEPARM card) = %i' % len(card)
        else:
            self.id = data[0]
            self.label = data[1]
            self.units = data[2]
            assert len(data) == 3, 'data = %s' % data

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the AEPARM object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['AEPARM', self.id, self.label, self.units]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class AESTAT(BaseCard):
    """
    Specifies rigid body motions to be used as trim variables in static
    aeroelasticity.

    +--------+------+--------+
    | AESTAT | ID   | LABEL  |
    +--------+------+--------+
    | AESTAT | 5001 | ANGLEA |
    +--------+------+--------+
    """
    type = 'AESTAT'

    _field_map = {
        1: 'id', 2:'label',
    }
    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.id = integer(card, 1, 'ID')
            self.label = string(card, 2, 'label')
            assert len(card) <= 3, 'len(AESTAT card) = %i' % len(card)
        else:
            self.id = data[0]
            self.label = data[1]
            assert len(data) == 2, 'data = %s' % data

    def cross_reference(self, model):
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the AESTAT object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['AESTAT', self.id, self.label]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class AESURF(BaseCard):
    """
    Specifies an aerodynamic control surface as a member of the set of
    aerodynamic extra points. The forces associated with this controller will
    be derived from rigid rotation of the aerodynamic model about the hinge
    line(s) and from AEDW, AEFORCE and AEPRESS input data. The mass properties
    of the control surface can be specified using an AESURFS entry.

    +--------+--------+-------+-------+-------+--------+--------+--------+--------+
    | AESURF |  ID    | LABEL | CID1  | ALID1 | CID2   | ALID2  | EFF    | LDW    |
    +--------+--------+-------+-------+-------+--------+--------+--------+--------+
    |        |  CREFC | CREFS | PLLIM | PULIM | HMLLIM | HMULIM | TQLLIM | TQULIM |
    +--------+--------+-------+-------+-------+--------+--------+--------+--------+
    """
    type = 'AESURF'
    _field_map = {
        1: 'aesid', 2:'label', 3:'cid1', 4:'alid1', 5:'cid2', 6:'alid2',
        7:'eff', 8:'ldw', 9:'crefc', 10:'crefs', 11:'pllim', 12:'pulim',
        13:'hmllim', 14:'hmulim', 15:'tqllim', '16':'tqulim',
    }

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Controller identification number
            self.aesid = integer(card, 1, 'aesid')
            #: Controller name.
            self.label = string(card, 2, 'label')

            #: Identification number of a rectangular coordinate system with a
            #: y-axis that defines the hinge line of the control surface
            #: component.
            self.cid1 = integer(card, 3, 'cid1')
            #: Identification of an AELIST Bulk Data entry that identifies all
            #: aerodynamic elements that make up the control surface
            #: component. (Integer > 0)
            self.alid1 = integer(card, 4, 'alid1')

            self.cid2 = integer_or_blank(card, 5, 'cid2')
            self.alid2 = integer_or_blank(card, 6, 'alid2')

            #: Control surface effectiveness. See Remark 4. (Real != 0.0;
            #: Default=1.0)
            self.eff = double_or_blank(card, 7, 'eff', 1.0)
            #: Linear downwash flag. See Remark 2.
            #: (Character, one of LDW or NOLDW; Default=LDW).
            self.ldw = string_or_blank(card, 8, 'ldw', 'LDW')
            #: Reference chord length for the control surface. (Real>0.0;
            #: Default=1.0)
            self.crefc = double_or_blank(card, 9, 'crefc', 1.0)
            #: Reference surface area for the control surface. (Real>0.0;
            #: Default=1.0)
            self.crefs = double_or_blank(card, 10, 'crefs', 1.0)
            #: Lower and upper deflection limits for the control surface in
            #: radians. (Real, Default = +/- pi/2)
            self.pllim = double_or_blank(card, 11, 'pllim', -pi / 2.)
            self.pulim = double_or_blank(card, 12, 'pulim', pi / 2.)
            #: Lower and upper hinge moment limits for the control surface in
            #: force-length units. (Real, Default = no limit) -> 1e8
            self.hmllim = double_or_blank(card, 13, 'hmllim')
            self.hmulim = double_or_blank(card, 14, 'hmulim')
            #: Set identification numbers of TABLEDi entries that provide the
            #: lower and upper deflection limits for the control surface as a
            #: function of the dynamic pressure. (Integer>0, Default = no limit)
            self.tqllim = integer_or_blank(card, 15, 'tqllim')
            self.tqulim = integer_or_blank(card, 16, 'tqulim')
            assert len(card) <= 17, 'len(AESURF card) = %i' % len(card)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def Cid1(self):
        if isinstance(self.cid1, integer_types):
            return self.cid1
        return self.cid1.cid

    def Cid2(self):
        if isinstance(self.cid2, integer_types) or self.cid2 is None:
            return self.cid2
        return self.cid2.cid

    def AELIST_id1(self):
        if isinstance(self.alid1, integer_types):
            return self.alid1
        return self.alid1.sid

    def AELIST_id2(self):
        if isinstance(self.alid2, integer_types) or self.alid2 is None:
            return self.alid2
        return self.alid2.sid

    def cross_reference(self, model):
        self.cid1 = model.Coord(self.Cid1())
        if self.cid2 is not None:
            self.cid2 = model.Coord(self.Cid2())
        self.alid1 = model.AELIST(self.AELIST_id1())
        if self.alid2:
            self.alid2 = model.AELIST(self.AELIST_id2())

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the AESURF object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['AESURF', self.aesid, self.label, self.Cid1(), self.AELIST_id1(),
                       self.Cid2(), self.AELIST_id2(), self.eff, self.ldw,
                       self.crefc, self.crefs, self.pllim, self.pulim, self.hmllim,
                       self.hmulim, self.tqllim, self.tqulim]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the AESURF object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        eff = set_blank_if_default(self.eff, 1.0)
        ldw = set_blank_if_default(self.ldw, 'LDW')
        crefc = set_blank_if_default(self.crefc, 1.0)
        crefs = set_blank_if_default(self.crefs, 1.0)

        pllim = set_blank_if_default(self.pllim, -pi / 2.)
        pulim = set_blank_if_default(self.pulim, pi / 2.)

        list_fields = ['AESURF', self.aesid, self.label, self.Cid1(), self.AELIST_id1(),
                       self.Cid2(), self.AELIST_id2(), eff, ldw, crefc, crefs,
                       pllim, pulim, self.hmllim, self.hmulim, self.tqllim,
                       self.tqulim]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AESURFS(BaseCard):  # not integrated
    """
    Optional specification of the structural nodes associated with an
    aerodynamic control surface that has been defined on an AESURF entry. The
    mass associated with these structural nodes define the control surface
    moment(s) of inertia about the hinge line(s).
    Specifies rigid body motions to be used as trim variables in static
    aeroelasticity.

    +---------+------+-------+---+-------+---+-------+
    |    1    |  2   |   3   | 4 |   5   | 6 |   7   |
    +---------+------+-------+---+-------+---+-------+
    | AESURFS | ID   | LABEL |   | LIST1 |   | LIST2 |
    +---------+------+-------+---+-------+---+-------+
    | AESURFS | 6001 | ELEV  |   | 6002  |   | 6003  |
    +---------+------+-------+---+-------+---+-------+
    """
    type = 'AESURFS'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.id = integer(card, 1, 'ID')
            self.label = string(card, 2, 'label')
            self.list1 = integer(card, 4, 'list1')
            self.list2 = integer(card, 6, 'list2')
            assert len(card) <= 7, 'len(AESURFS card) = %i' % len(card)
        else:
            self.id = data[0]
            self.label = data[1]
            self.list1 = data[2]
            self.list2 = data[3]
            assert len(data) == 4, 'data = %s' % data

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the AESURFS object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['AESURFS', self.id, self.label, None, self.list1, None,
                       self.list2]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class Aero(BaseCard):
    """Base class for AERO and AEROS cards."""
    def __init__(self, card, data):
        self.symXY = None
        self.symXZ = None

    def IsSymmetricalXY(self):
        self.deprecated('IsSymmetricalXY()', 'is_symmetric_xy()', '0.7')
        return self.is_symmetric_xy()

    def IsSymmetricalXZ(self):
        self.deprecated('IsSymmetricalXZ()', 'is_symmetric_xz()', '0.7')
        return self.is_symmetric_xz()

    def IsAntiSymmetricalXY(self):
        self.deprecated('IsAntiSymmetricalXY()', 'is_anti_symmetric_xy()', '0.7')
        return self.is_anti_symmetric_xy()

    def IsAntiSymmetricalXZ(self):
        self.deprecated('IsAntiSymmetricalXZ()', 'is_anti_symmetric_xz()', '0.7')
        return self.is_anti_symmetric_xz()

    def EnableGroundEffect(self):
        self.deprecated('EnableGroundEffect()', 'set_ground_effect(True)', '0.7')
        self.set_ground_effect(True)

    def DisableGroundEffect(self):
        self.deprecated('DisableGroundEffect()', 'set_ground_effect(False)', '0.7')
        self.set_ground_effect(False)

    def is_symmetric_xy(self):
        if self.symXY == 1:
            return True
        return False

    def is_symmetric_xz(self):
        if self.symXZ == 1:
            return True
        return False

    def is_anti_symmetric_xy(self):
        if self.symXY == -1:
            return True
        return False

    def is_anti_symmetric_xz(self):
        if self.symXY == -1:
            return True
        return False

    def set_ground_effect(self, enable):
        if enable:
            self.symXY = -1
        else:
            self.symXY = 1


class AERO(Aero):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.

    +------+-------+----------+------+--------+-------+-------+
    | 1    | 2     | 3        | 4    | 5      | 6     | 7     |
    +------+-------+----------+------+--------+-------+-------+
    | AERO | ACSID | VELOCITY | REFC | RHOREF | SYMXZ | SYMXY |
    +------+-------+----------+------+--------+-------+-------+
    | AERO | 3     | 1.3+4    | 100. |  1.-5  | 1     | -1    |
    +------+-------+----------+------+--------+-------+-------+
    """
    type = 'AERO'
    _field_map = {
        1: 'acsid', 2:'velocity', 3:'cRef', 4:'rhoRef', 5:'symXZ',
        6:'symXY',
    }

    def __init__(self, card=None, data=None, comment=''):
        Aero.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.acsid = integer_or_blank(card, 1, 'acsid', 0)
            self.velocity = double_or_blank(card, 2, 'velocity')
            self.cRef = double(card, 3, 'cRef')
            self.rhoRef = double(card, 4, 'rhoRef')
            self.symXZ = integer_or_blank(card, 5, 'symXZ', 0)
            self.symXY = integer_or_blank(card, 6, 'symXY', 0)
            assert len(card) <= 7, 'len(AERO card) = %i' % len(card)
        else:
            self.acsid = data[0]
            self.velocity = data[1]
            self.cRef = data[2]
            self.rhoRef = data[3]
            self.symXZ = data[4]
            self.symXY = data[5]
            assert len(data) == 6, 'data = %s' % data

        # T is the tabular function
        #angle = self.wg*self.t*(t-(x-self.x0)/self.V)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the AERO object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['AERO', self.acsid, self.velocity, self.cRef,
                       self.rhoRef, self.symXZ, self.symXY]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the AERO object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        symXZ = set_blank_if_default(self.symXZ, 0)
        symXY = set_blank_if_default(self.symXY, 0)
        list_fields = ['AERO', self.acsid, self.velocity, self.cRef,
                       self.rhoRef, symXZ, symXY]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AEROS(Aero):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.

    +-------+-------+-------+------+------+-------+------+-------+
    | 1     | 2     | 3     | 4    | 5    | 6     | 7    |   8   |
    +-------+-------+-------+------+------+-------+------+-------+
    | AEROS | ACSID | RCSID | REFC | REFB | REFS  |SYMXZ | SYMXY |
    +-------+-------+-------+------+------+-------+------+-------+

    +-------+-------+-------+------+------+-------+------+-------+
    | AEROS | 10    | 20    | 10.  | 100. | 1000. | 1    |       |
    +-------+-------+-------+------+------+-------+------+-------+
    """
    type = 'AEROS'
    _field_map = {
        1: 'acsid', 2:'rcsid', 3:'cRef', 4:'bRef', 5:'Sref',
        6:'symXZ', 7:'symXY',
    }

    def __init__(self, card=None, data=None, comment=''):
        Aero.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.acsid = integer_or_blank(card, 1, 'acsid', 0)
            self.rcsid = integer_or_blank(card, 2, 'rcsid', 0)
            self.cRef = double(card, 3, 'cRef')
            self.bRef = double(card, 4, 'bRef')
            self.Sref = double(card, 5, 'Sref')
            self.symXZ = integer_or_blank(card, 6, 'symXZ', 0)
            self.symXY = integer_or_blank(card, 7, 'symXY', 0)
            assert len(card) <= 8, 'len(AEROS card) = %i' % len(card)
        else:
            self.acsid = data[0]
            self.rcsid = data[1]
            self.cRef = data[2]
            self.bRef = data[3]
            self.Sref = data[4]
            self.symXZ = data[5]
            self.symXY = data[6]
            assert len(data) == 7, 'data = %s' % data

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the AEROS object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['AEROS', self.acsid, self.rcsid, self.cRef,
                       self.bRef, self.Sref, self.symXZ, self.symXY]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        :param self:
          the AEROS object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        symXZ = set_blank_if_default(self.symXZ, 0)
        symXY = set_blank_if_default(self.symXY, 0)
        list_fields = ['AEROS', self.acsid, self.rcsid, self.cRef,
                       self.bRef, self.Sref, symXZ, symXY]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CSSCHD(Aero):
    """
    Defines a scheduled control surface deflection as a function of Mach number
    and angle of attack.

    +--------+-----+-------+--------+-------+-------+
    |    1   |  2  |   3   |   4    |   5   |   6   |
    +--------+-----+-------+--------+-------+-------+
    | CSSCHD | SlD | AESID | LALPHA | LMACH | LSCHD |
    +--------+-----+-------+--------+-------+-------+
    """
    type = 'CSSCHD'
    _field_map = {
        1: 'sid', 2:'aesid', 3:'lAlpha', 4:'lMach', 5:'lSchd',
    }

    def __init__(self, card=None, data=None, comment=''):
        Aero.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.aesid = integer(card, 2, 'aesid')             # AESURF
            self.lAlpha = integer_or_blank(card, 3, 'lAlpha')  # AEFACT
            self.lMach = integer_or_blank(card, 4, 'lMach')    # AEFACT
            self.lSchd = integer(card, 5, 'lSchd')             # AEFACT
            assert len(card) <= 6, 'len(CSSCHD card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.aesid = data[1]  # AESURF
            self.lAlpha = data[2]  # AEFACT
            self.lMach = data[3]  # AEFACT
            self.lSchd = data[4]  # AEFACT

    def cross_reference(self, model):
        """
        Cross links the card

        :param self:   the CSSCHD object pointer
        :param model:  the BDF object
        :type model:   BDF()
        """
        msg = ' which is required by ASSCHD sid=%s' % self.sid
        self.aesid = model.AESurf(self.aesid, msg=msg)
        self.lAlpha = model.AEFact(self.lAlpha, msg=msg)
        self.lMach = model.AEFact(self.lMach, msg=msg)
        self.lSchd = model.AEFact(self.lSchd, msg=msg)

    def AESid(self):
        if isinstance(self.aesid, integer_types):
            return self.aesid
        return self.aesid.aesid

    def LAlpha(self):
        if isinstance(self.lAlpha, integer_types):
            return self.lAlpha
        return self.lAlpha.sid

    def LMach(self):
        if isinstance(self.lMach, integer_types):
            return self.lMach
        return self.lMach.sid

    def LSchd(self):
        if isinstance(self.lSchd, integer_types):
            return self.lSchd
        return self.lSchd.sid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the CSSCHD object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['CSSCHD', self.sid, self.AESid(), self.LAlpha(),
                       self.LMach(), self.LSchd()]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO1(BaseCard):
    """
    Defines an aerodynamic macro element (panel) in terms of two leading edge
    locations and side chords. This is used for Doublet-Lattice theory for
    subsonic aerodynamics and the ZONA51 theory for supersonic aerodynamics.

    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    | CAERO1 | EID | PID | CP | NSPAN | NCHORD |  LSPAN | LCHORD | IGID |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |        |  X1 | Y1  | Z1 | X12   | X4     | Y4     | Z4     | X43  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+

    Attributes
    ----------
    sid : int
       element id
    pid : int, PAERO1
       int : PAERO1 ID
       PAERO1 : PAERO1 object (xref)

    cp : int, CORDx
       int : coordinate system
       CORDx : Coordinate object (xref)
    nspan : int
       int > 0 : N spanwise boxes distributed evenly
       int = 0 : use lchord
    nchord : int
       int > 0 : N chordwise boxes distributed evenly
       int = 0 : use lchord
    lspan : int, AEFACT
       int > 0 : AEFACT reference for non-uniform nspan
       int = 0 : use nspan
       AEFACT : AEFACT object  (xref)
    lchord : int, AEFACT
       int > 0 : AEFACT reference for non-uniform nchord
       int = 0 : use nchord
       AEFACT : AEFACT object  (xref)
    igid : int
       Group number
    p1 : (1, 3) ndarray float
       xyz location of point 1 (leading edge; inboard)
    p4 : (1, 3) ndarray float
       xyz location of point 4 (leading edge; outboard)
    """
    type = 'CAERO1'
    _field_map = {
        1: 'sid', 2:'pid', 3:'cp', 4:'nspan', 5:'nchord',
        6:'lspan', 7:'lchord', 8:'igid', 12:'x12', 16:'x43',
    }
    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the CAERO1 card

        :param self:  the CAERO1 object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        if n == 9:
            return self.p1[0]
        elif n == 10:
            return self.p1[1]
        elif n == 11:
            return self.p1[2]

        elif n == 13:
            return self.p4[0]
        elif n == 14:
            return self.p4[1]
        elif n == 15:
            return self.p4[2]
        else:
            raise KeyError('Field %r is an invalid %s entry.' % (n, self.type))

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the CAERO1 card

        :param self:  the CAERO1 object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        if n == 9:
            self.p1[0] = value
        elif n == 10:
            self.p1[1] = value
        elif n == 11:
            self.p1[2] = value

        elif n == 13:
            self.p4[0] = value
        elif n == 14:
            self.p4[1] = value
        elif n == 15:
            self.p4[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        r"""
        ::

          1
          | \
          |   \
          |     \
          |      4
          |      |
          |      |
          2------3
        """
        if comment:
            self._comment = comment
        if card:
            #: Element identification number
            self.eid = integer(card, 1, 'eid')

            #: Property identification number of a PAERO2 entry.
            self.pid = integer(card, 2, 'pid')

            #: Coordinate system for locating point 1.
            self.cp = integer_or_blank(card, 3, 'cp', 0)

            self.nspan = integer_or_blank(card, 4, 'nspan', 0)
            self.nchord = integer_or_blank(card, 5, 'nchord', 0)

            self.lspan = integer_or_blank(card, 6, 'lspan', 0)
            self.lchord = integer_or_blank(card, 7, 'lchord', 0)

            self.igid = integer(card, 8, 'igid')

            self.p1 = array([double_or_blank(card, 9, 'x1', 0.0),
                             double_or_blank(card, 10, 'y1', 0.0),
                             double_or_blank(card, 11, 'z1', 0.0)])
            self.x12 = double_or_blank(card, 12, 'x12', 0.)

            self.p4 = array([double_or_blank(card, 13, 'x4', 0.0),
                             double_or_blank(card, 14, 'y4', 0.0),
                             double_or_blank(card, 15, 'z4', 0.0)])
            self.x43 = double_or_blank(card, 16, 'x43', 0.)
            if self.nspan == 0 and self.lspan == 0:
                msg = 'NSPAN or LSPAN must be greater than 0'
                raise ValueError(msg)
            if self.nspan != 0 and self.lspan != 0:
                msg = 'Either NSPAN or LSPAN must 0'
                raise ValueError(msg)

            if self.nchord == 0 and self.lchord == 0:
                msg = 'NCHORD or LCHORD must be greater than 0'
                raise ValueError(msg)
            if self.nchord != 0 and self.lchord != 0:
                msg = 'Either NCHORD or LCHORD must 0'
                raise ValueError(msg)

            assert len(card) <= 17, 'len(CAERO1 card) = %i' % len(card)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def _init_ids(self):
        """
        Fill `self.box_ids` with the sub-box ids. Shape is (nchord, nspan)

        Parameters
        ----------
        self : CAERO1 obj
            The CAERO1 object.
        """
        nchord, nspan = self.shape
        self.box_ids = zeros((nchord, nspan), dtype='int32')
        for ichord in range(nchord):
            for ispan in range(nspan):
                self.box_ids[ichord, ispan] = self.eid + ichord + ispan * nchord

    def Points(self):
        self.deprecated('Points()', 'get_points()', '0.7')
        return self.get_points()

    def SetPoints(self, points):
        self.deprecated('SetPoints(points)', 'set_points(points)', '0.7')
        return self.set_points(points)

    def Cp(self):
        if isinstance(self.cp, integer_types):
            return self.cp
        return self.cp.cid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid.pid

    def cross_reference(self, model):
        """
        Cross links the card

        Parameters
        ----------
        self : CAERO1 obj
            The CAERO1 object.
        model : BDF obj
            The BDF object.
        """
        msg = ' which is required by CAERO1 eid=%s' % self.eid
        self.pid = model.PAero(self.pid, msg=msg)
        self.cp = model.Coord(self.cp, msg=msg)
        if self.nchord == 0:
            assert isinstance(self.lchord, integer_types), self.lchord
            self.lchord = model.AEFact(self.lchord, msg)
        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan = model.AEFact(self.lspan, msg)
        self._init_ids()

    @property
    def min_max_eid(self):
        """
        Gets the min and max element ids of the CAERO card

        Parameters
        ----------
        self : CAERO1 obj
            The CAERO1 object.

        Returns
        -------
        min_max_eid : (2, ) list
            The [min_eid, max_eid]
        """
        nchord, nspan = self.shape
        return [self.eid, self.eid + nchord * nspan]

    def get_points(self):
        """
        Get the 4 corner points for the CAERO card

        Parameters
        ----------
        self : CAERO1 obj
            The CAERO1 object.

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame
        """
        p1 = self.cp.transform_node_to_global(self.p1)
        p4 = self.cp.transform_node_to_global(self.p4)
        p2 = p1 + array([self.x12, 0., 0.])
        p3 = p4 + array([self.x43, 0., 0.])
        return [p1, p2, p3, p4]

    @property
    def shape(self):
        if self.nchord == 0:
            x = self.lchord.Di
            nchord = len(x) - 1
        else:
            nchord = self.nchord

        if self.nspan == 0:
            y = self.lspan.Di
            nspan = len(y) - 1
        else:
            nspan = self.nspan
        if nchord < 1 or nspan < 1:
            msg = '%s eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                self.type, self.eid, self.nchord, self.nspan,
                self.lchord, self.lspan)
            raise RuntimeError(msg)
        return nchord, nspan

    def get_npanel_points_elements(self, box_ids=None):
        """
        Gets the number of sub-points and sub-elements for the CAERO card

        Parameters
        ----------
        self : CAERO1 obj
            The CAERO1 object.

        Returns
        -------
        npoints : int
            The number of nodes for the CAERO
        nelmements : int
            The number of elements for the CAERO
        """
        nchord, nspan = self.shape
        nelements = nchord * nspan
        npoints = (nchord + 1) * (nspan + 1)
        return npoints, nelements

    @property
    def xy(self):
        """
        Parameters
        ----------
        self : CAERO1 obj
            The CAERO1 object.

        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel
        """
        if self.nchord == 0:
            x = self.lchord.Di
            nchord = len(x) - 1
        else:
            nchord = self.nchord
            x = linspace(0., 1., nchord + 1)

        if self.nspan == 0:
            y = self.lspan.Di
            nspan = len(y) - 1
        else:
            nspan = self.nspan
            y = linspace(0., 1., nspan + 1)

        if nchord < 1 or nspan < 1:
            msg = '%s eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                self.type, self.eid, self.nchord, self.nspan,
                self.lchord, self.lspan)
            raise RuntimeError(msg)
        return x, y

    def panel_points_elements(self):
        """
        Gets the sub-points and sub-elements for the CAERO card

        Parameters
        ----------
        self : CAERO1 obj
            The CAERO1 object.

        Returns
        -------
        points : (nnodes,3) ndarray of floats
            the array of points
        elements : (nelements,4) ndarray of integers
            the array of point ids
        """
        p1, p2, p3, p4 = self.get_points()
        x, y = self.xy
        return points_elements_from_quad_points(p1, p2, p3, p4, x, y)

    def set_points(self, points):
        self.p1 = points[0]
        self.p2 = points[1]
        self.p3 = points[2]
        self.p4 = points[3]
        self.x12 = self.p2[0] - self.p1[0]
        self.x43 = self.p4[0] - self.p3[0]

    def uncross_reference(self):
        self.pid = self.Pid()
        self.cp = self.Cp()
        self.lchord = self.get_LChord()
        self.lspan = self.get_LSpan()


    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Parameters
        ----------
        self : CAERO 1 obj
            the CAERO1 object pointer.

        Returns
        -------
        fields : list
          the fields that define the card
        """
        lchord = self.get_LChord()
        lspan = self.get_LSpan()
        list_fields = (['CAERO1', self.eid, self.Pid(), self.Cp(), self.nspan,
                        self.nchord, lspan, lchord, self.igid, ] +
                       list(self.p1) + [self.x12] + list(self.p4) + [self.x43])
        return list_fields

    def get_LChord(self):
        if isinstance(self.lchord, integer_types):
            return self.lchord
        return self.lchord.sid

    def get_LSpan(self):
        if isinstance(self.lspan, integer_types):
            return self.lspan
        return self.lspan.sid

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Parameters
        ----------
        self : CAERO1 obj
            The CAERO1 object.

        Returns
        -------
        fields : LIST
            The fields that define the card
        """
        cp = set_blank_if_default(self.Cp(), 0)
        nspan = set_blank_if_default(self.nspan, 0)
        nchord = set_blank_if_default(self.nchord, 0)
        lchord = set_blank_if_default(self.get_LChord(), 0)
        lspan = set_blank_if_default(self.get_LSpan(), 0)
        list_fields = (['CAERO1', self.eid, self.Pid(), cp, nspan, nchord,
                        lspan, lchord, self.igid] + list(self.p1) +
                       [self.x12] + list(self.p4) + [self.x43])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO2(BaseCard):
    """
    Aerodynamic Body Connection
    Defines aerodynamic slender body and interference elements for
    Doublet-Lattice aerodynamics.
    """
    type = 'CAERO2'
    _field_map = {
        1: 'sid', 2:'pid', 3:'cp', 4:'nsb', 5:'lsb',
        6:'nint', 7:'lint', 8:'igid', 12:'x12',
    }
    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the CAERO2 card

        Parameters
        ----------
        self : CAERO2 obj
            The CAERO2 object.
        n : int
            The field number to update

        Returns
        -------
        value : int, float, None
            The value for the appropriate field
        """
        if n == 9:
            return self.p1[0]
        elif n == 10:
            return self.p1[1]
        elif n == 11:
            return self.p1[2]
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the CAERO2 card

        Parameters
        ----------
        self : CAERO2 obj
            The CAERO2 object.
        n : int
            The field number to update
        value : int, float, None
            The value for the appropriate field
        """
        if n == 9:
            self.p1[0] = value
        elif n == 10:
            self.p1[1] = value
        elif n == 11:
            self.p1[2] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        """
        ::

          1 \
          |   \
          |     \
          |      3
          |      |
          |      |
          2------4
        """
        if comment:
            self._comment = comment
        if card:
            #: Element identification number
            self.eid = integer(card, 1, 'eid')

            #: Property identification number of a PAERO2 entry.
            self.pid = integer(card, 2, 'pid')

            #: Coordinate system for locating point 1.
            self.cp = integer_or_blank(card, 3, 'cp', 0)

            #: Number of slender body elements. If NSB > 0, then NSB equal
            #: divisions are assumed; if zero or blank, specify a list of
            #: divisions in LSB. (Integer >= 0)
            self.nsb = integer_or_blank(card, 4, 'nsb', 0)

            #: Number of interference elements. If NINT > 0, then NINT equal
            #: divisions are assumed; if zero or blank, specify a list of
            #: divisions in LINT. (Integer >= 0)
            self.nint = integer_or_blank(card, 5, 'nint', 0)

            if self.nsb == 0:
                #: ID of an AEFACT Bulk Data entry for slender body division
                #: points; used only if NSB is zero or blank. (Integer >= 0)
                self.lsb = integer(card, 6, 'nsb=%s lsb' % self.nsb)
            else:
                self.lsb = blank(card, 6, 'nsb=%s lsb' % self.nsb)

            if self.nint == 0:
                #: ID of an AEFACT data entry containing a list of division
                #: points for interference elements; used only if NINT is zero
                #: or blank. (Integer > 0)
                self.lint = integer(card, 7, 'nint=%s lint' % self.nint)
            else:
                self.lint = blank(card, 7, 'nint=%s lint' % self.nint)

            #: Interference group identification. Aerodynamic elements with
            #: different IGIDs are uncoupled. (Integer >= 0)
            self.igid = integer(card, 8, 'igid')

            #: Location of point 1 in coordinate system CP
            self.p1 = array([double_or_blank(card, 9, 'x1', 0.0),
                             double_or_blank(card, 10, 'y1', 0.0),
                             double_or_blank(card, 11, 'z1', 0.0)])

            #: Length of body in the x-direction of the aerodynamic coordinate
            #: system.  (Real > 0)
            self.x12 = double_or_blank(card, 12, 'x12', 0.)
            assert len(card) <= 13, 'len(CAERO2 card) = %i' % len(card)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def Points(self):
        self.deprecated('Points()', 'get_points()', '0.7')
        return self.get_points()

    def SetPoints(self, points):
        self.deprecated('SetPoints(points)', 'set_points(points)', '0.7')
        return self.set_points(points)

    def Cp(self):
        if isinstance(self.cp, integer_types):
            return self.cp
        return self.cp.cid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid.pid

    def Lsb(self):  # AEFACT
        if isinstance(self.lsb, integer_types):
            return self.lsb
        return self.lsb.sid

    def cross_reference(self, model):
        """
        Cross links the card

        Parameters
        ----------
        self : CAERO2 obj
            The CAERO2 object.
        model : BDF obj
            The BDF object.
        """
        msg = ' which is required by CAERO2 eid=%s' % self.eid
        self.pid = model.PAero(self.pid, msg=msg)  # links to PAERO2
        self.cp = model.Coord(self.cp, msg=msg)
        #self.lsb = model.AeFact(self.lsb, msg=msg) # not added

    def get_points(self):
        (p1, matrix) = self.cp.transformToGlobal(self.p1)
        p2 = p1 + array([self.x12, 0., 0.])
        #print("x12 = %s" % self.x12)
        #print("pcaero[%s] = %s" % (self.eid, [p1,p2]))
        return [p1, p2]

    def set_points(self, points):
        self.p1 = points[0]
        self.p2 = points[1]
        x12 = self.p2 - self.p1
        self.x12 = x12[0]

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Parameters
        ----------
        self : CAERO2 obj
            The CAERO2 object.

        Returns
        -------
        fields : list
            The fields that define the card
        """
        list_fields = (['CAERO2', self.eid, self.Pid(), self.Cp(), self.nsb,
                        self.nint, self.lsb, self.lint, self.igid, ] + list(self.p1)
                       + [self.x12])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Parameters
        ----------
        self : CAERO2 obj
            The CAERO2 object.

        Returns
        -------
        fields : list
            The fields that define the card
        """
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = (['CAERO2', self.eid, self.Pid(), cp, self.nsb, self.nint,
                        self.lsb, self.lint, self.igid, ] + list(self.p1) +
                       [self.x12])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO3(BaseCard):
    type = 'CAERO3'
    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Element identification number
            self.eid = integer(card, 1, 'eid')
            #: Property identification number of a PAERO3 entry.
            self.pid = integer(card, 2, 'pid')
            #: Coordinate system for locating point 1.
            self.cp = integer_or_blank(card, 3, 'cp', 0)
            self.list_w = integer(card, 4, 'list_w')
            self.list_c1 = integer_or_blank(card, 5, 'list_c1')
            self.list_c2 = integer_or_blank(card, 6, 'list_c2')
            self.p1 = array([double_or_blank(card, 9, 'x1', 0.0),
                             double_or_blank(card, 10, 'y1', 0.0),
                             double_or_blank(card, 11, 'z1', 0.0)])
            self.x12 = double(card, 12, 'x12')
            assert self.x12 > 0., 'x12=%s' % self.x12
            self.p4 = array([double_or_blank(card, 13, 'x4', 0.0),
                             double_or_blank(card, 14, 'y4', 0.0),
                             double_or_blank(card, 15, 'z4', 0.0)])
            self.x43 = double_or_blank(card, 16, 'x43', 0.0)
            assert len(card) <= 17, 'len(CAERO3 card) = %i' % len(card)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def cross_reference(self, model):
        """
        Cross links the card

        Parameters
        ----------
        self : CAERO3 obj
            The CAERO3 object.
        model : BDF obj
            The BDF object.
        """
        msg = ' which is required by CAERO3 eid=%s' % self.eid
        self.pid = model.PAero(self.pid, msg=msg)  # links to PAERO3
        self.cp = model.Coord(self.cp, msg=msg)
        #self.list_w = model.AEFact(self.list_w, msg=msg)   # not added
        #self.list_c1 = model.AEFact(self.list_c1, msg=msg) # not added
        #self.list_c2 = model.AEFact(self.list_c2, msg=msg) # not added

    def Cp(self):
        if isinstance(self.cp, integer_types):
            return self.cp
        return self.cp.cid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid.pid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Parameters
        ----------
        self : CAERO3 obj
            The CAERO3 object.

        Returns
        -------
        fields : list
            The fields that define the card
        """
        list_fields = (['CAERO3', self.eid, self.Pid(), self.Cp(), self.list_w,
                        self.list_c1, self.list_c2, None, None] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Parameters
        ----------
        self : CAERO2 obj
            The CAERO2 object.

        Returns
        -------
        fields : list
            The fields that define the card
        """
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = (['CAERO3', self.eid, self.Pid(), cp, self.list_w,
                        self.list_c1, self.list_c2, None, None] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO4(BaseCard):
    type = 'CAERO4'
    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Element identification number
            self.eid = integer(card, 1, 'eid')
            #: Property identification number of a PAERO4 entry.
            self.pid = integer(card, 2, 'pid')
            #: Coordinate system for locating point 1.
            self.cp = integer_or_blank(card, 3, 'cp', 0)

            self.nspan = integer_or_blank(card, 4, 'nspan', 0)
            self.lspan = integer_or_blank(card, 5, 'lspan', 0)

            self.p1 = array([double_or_blank(card, 9, 'x1', 0.0),
                             double_or_blank(card, 10, 'y1', 0.0),
                             double_or_blank(card, 11, 'z1', 0.0)])
            self.x12 = double_or_blank(card, 12, 'x12', 0.)

            self.p4 = array([double_or_blank(card, 13, 'x4', 0.0),
                             double_or_blank(card, 14, 'y4', 0.0),
                             double_or_blank(card, 15, 'z4', 0.0)])
            self.x43 = double_or_blank(card, 16, 'x43', 0.)
            assert len(card) <= 17, 'len(CAERO4 card) = %i' % len(card)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def get_points(self):
        p1, matrix = self.cp.transformToGlobal(self.p1)
        p4, matrix = self.cp.transformToGlobal(self.p4)
        p2 = p1 + array([self.x12, 0., 0.])
        p3 = p4 + array([self.x43, 0., 0.])
        return [p1, p2, p3, p4]

    def cross_reference(self, model):
        """
        Cross links the card

        Parameters
        ----------
        self : CAERO4 obj
            The CAERO4 object.
        model : BDF obj
            The BDF object.
        """
        msg = ' which is required by CAERO4 eid=%s' % self.eid
        #self.pid = model.PAero(self.pid, msg=msg)  # links to PAERO4 (not added)
        self.cp = model.Coord(self.cp, msg=msg)

    def Cp(self):
        if isinstance(self.cp, integer_types):
            return self.cp
        return self.cp.cid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid.pid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Parameters
        ----------
        self : CAERO4 obj
            The CAERO4 object.

        Returns
        -------
        fields : list
            The fields that define the card
        """
        list_fields = (['CAERO4', self.eid, self.Pid(), self.Cp(), self.nspan,
                        self.lspan, None, None, None,] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Parameters
        ----------
        self : CAERO4 obj
            The CAERO4 object.

        Returns
        -------
        fields : list
            The fields that define the card
        """
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = (['CAERO4', self.eid, self.Pid(), cp, self.nspan,
                        self.lspan, None, None, None,] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO5(BaseCard):
    type = 'CAERO5'
    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Element identification number
            self.eid = integer(card, 1, 'eid')
            #: Property identification number of a PAERO5 entry.
            self.pid = integer(card, 2, 'pid')
            #: Coordinate system for locating point 1.
            self.cp = integer_or_blank(card, 3, 'cp', 0)
            self.nspan = integer_or_blank(card, 4, 'nspan', 0)
            self.lspan = integer_or_blank(card, 5, 'lspan', 0)
            self.ntheory = integer_or_blank(card, 6, 'ntheory')
            self.nthick = integer_or_blank(card, 7, 'nthick')
            # 8 - blank
            self.p1 = array([double_or_blank(card, 9, 'x1', 0.0),
                             double_or_blank(card, 10, 'y1', 0.0),
                             double_or_blank(card, 11, 'z1', 0.0)])
            self.x12 = double(card, 12, 'x12')
            assert self.x12 > 0., 'x12=%s' % self.x12
            self.p4 = array([double_or_blank(card, 13, 'x4', 0.0),
                             double_or_blank(card, 14, 'y4', 0.0),
                             double_or_blank(card, 15, 'z4', 0.0)])
            self.x43 = double_or_blank(card, 16, 'x43', 0.0)
            assert len(card) <= 17, 'len(CAERO3 card) = %i' % len(card)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

        if not (self.nspan > 0 or self.lspan > 0):
            msg = 'nspan=%r or lspan=%r must be > 0' % (self.nspan, self.lspan)
            raise ValueError(msg)
        if not (self.x12 > 0.0 or self.x43 > 0.0):
            msg = 'x12=%r or x43=%r must be > 0.0' % (self.x12, self.x43)
            raise ValueError(msg)

    def cross_reference(self, model):
        """
        Cross links the card

        Parameters
        ----------
        self : CAERO5 obj
            The CAERO5 object.
        model : BDF obj
            The BDF object.
        """
        msg = ' which is required by CAERO5 eid=%s' % self.eid
        #self.pid = model.PAero(self.pid, msg=msg)  # links to PAERO5 (not added)
        self.cp = model.Coord(self.cp, msg=msg)
        self.lspan = model.AEFact(self.lspan, msg=msg)

    def get_points(self):
        p1, matrix = self.cp.transformToGlobal(self.p1)
        p4, matrix = self.cp.transformToGlobal(self.p4)
        p2 = p1 + array([self.x12, 0., 0.])
        p3 = p4 + array([self.x43, 0., 0.])
        return [p1, p2, p3, p4]

    def get_npanel_points_elements(self):
        msg = '%s eid=%s nspan=%s lspan=%s' % (self.type,
                                               self.eid,
                                               self.nspan,
                                               self.lspan)
        if self.nspan == 0:
            y = self.lspan.Di
            nspan = len(y) - 1
        else:
            nspan = self.nspan
        assert nspan >= 1, msg

        nchord = 1
        nelements = nchord * nspan
        npoints = (nchord + 1) * (nspan + 1)
        return npoints, nelements

    def panel_points_elements(self):
        p1, p2, p3, p4 = self.get_points()

        msg = '%s eid=%s nspan=%s lspan=%s' % (self.type,
                                               self.eid,
                                               self.nspan,
                                               self.lspan)
        if self.nspan == 0:
            y = self.lspan.Di
            nspan = len(y) - 1
        else:
            nspan = self.nspan
            y = linspace(0., 1., nspan + 1)
        assert nspan >= 1, msg

        x = array([0., 1.], dtype='float64')
        assert nspan >= 1, msg

        return points_elements_from_quad_points(p1, p2, p3, p4, x, y)

    def c1_c2(self, mach):
        p1, p2, p3, p4 = self.get_points()
        #i = p2 - p1
        #ihat = i / norm(i)
        #k = cross(ihat, p4-p1)
        #khat = k / norm(k)
        #jhat = cross(khat, ihat)
        #b = self.p4 - self.p1
        L = norm(p4 - p1)

        # b = L * cos(Lambda)
        # ci = L * sin(Lambda)

        if self.ntheory == 0:
            # piston theory
            pass
        elif self.ntheory == 1:
            gamma = 1.4
            #Lambda = 0.
            c1 = 1.

            secL = 1 / cos(Lambda)
            secL2 = secL ** 2
            ma2_secL2 = mach ** 2 - secL2
            c1 = mach / ma2_secL2 ** 0.5
            c2 = (mach ** 4 * (gamma + 1) - 4 * secL2 * ma2_secL2) / (4 * ma2_secL2 ** 2)
        else:
            gamma = 1.4

            # the advance in x
            ci = p4[0] - p1[0]

            # sweep angle
            Lambda = arcsin(ci / L)

            secL = 1 / cos(Lambda)
            secL2 = secL ** 2
            ma2_secL2 = mach ** 2 - secL2
            c1 = mach / ma2_secL2 ** 0.5
            c2 = (mach ** 4 * (gamma + 1) - 4 * secL2 * ma2_secL2) / (4 * ma2_secL2 ** 2)

    def Cp(self):
        if isinstance(self.cp, integer_types):
            return self.cp
        return self.cp.cid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid.pid

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Parameters
        ----------
        self : CAERO5 obj
            The CAERO5 object.

        Returns
        -------
        fields : list
            The fields that define the card
        """
        nspan = self.nspan
        lspan = self.LSpan()
        ntheory = self.ntheory
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = (['CAERO5', self.eid, self.Pid(), cp, nspan, lspan,
                        ntheory, self.nthick, None,] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


def points_elements_from_quad_points(p1, p2, p3, p4, x, y):
    nx = x.shape[0]
    ny = y.shape[0]
    assert nx > 1
    assert ny > 1

    nelements = (nx - 1) * (ny - 1)
    npoints = nx * ny

    # shape the vectors so we can multiply them
    x = x.reshape((1, nx))
    y = y.reshape((1, ny))
    p1 = p1.reshape(1, 3)
    p2 = p2.reshape(1, 3)
    p3 = p3.reshape(1, 3)
    p4 = p4.reshape(1, 3)

    # x repeats ny times and varies slowly
    # y repeats nx times and varies quickly
    xv = repeat(x, ny, axis=1).reshape(npoints, 1)
    yv = repeat(y, nx, axis=0).reshape(npoints, 1)

    # calculate the points a and b xv% along the chord
    a = xv * p2 + (1 - xv) * p1
    b = xv * p3 + (1 - xv) * p4

    # calculate the point yv% along the span
    points = yv * b + (1 - yv) * a
    assert points.shape == (npoints, 3), "npoints=%s shape=%s" % (npoints, str(points.shape))

    # create a matrix with the point counter
    ipoints = arange(npoints, dtype='int32').reshape((nx, ny))

    # move around the CAERO quad and apply ipoints
    elements = zeros((nelements, 4), dtype='int32')
    elements[:, 0] = ipoints[:-1, :-1].ravel()  # (i,  j  )
    elements[:, 1] = ipoints[1:, :-1].ravel()   # (i+1,j  )
    elements[:, 2] = ipoints[1:, 1:].ravel()    # (i+1,j+1)
    elements[:, 3] = ipoints[:-1, 1:].ravel()   # (i,j+1  )
    return points, elements


class PAERO5(BaseCard):  ## TODO: not integrated
    def __init__(self, card=None, data=None):
        """
        +--------+-------+--------+--------+---------+-------+-------+-------+
        | PAERO5 | PID   | NALPHA | LALPHA | NXIS    | LXIS  | NTAUS | LTAUS |
        +--------+-------+--------+--------+---------+-------+-------+-------+
        |        | CAOC1 | CAOC2  | CAOC3  | CAOC4   | CAOC5 |       |       |
        +--------+-------+--------+--------+---------+-------+-------+-------+

        +--------+-------+--------+--------+---------+-------+-------+-------+
        | PAERO5 | 7001  |   1    |  702   |    1    | 701   |   1   |  700  |
        +--------+-------+--------+--------+---------+-------+-------+-------+
        |        |  0.0  |  0.0   |  5.25  | 3.99375 |  0.0  |       |       |
        +--------+-------+--------+--------+---------+-------+-------+-------+
        """
        if card:
            self.pid = integer(card, 1, 'property_id')
            self.nalpha = integer_or_blank(card, 2, 'nalpha', default=0)
            self.lalpha = integer_or_blank(card, 3, 'lalpha', default=0)

            # number of dimensionless chord coordinates in zeta ()
            self.nxis = integer_or_blank(card, 4, 'nxis', default=0)

            # ID of AEFACT that lists zeta
            self.lxis = integer_or_blank(card, 5, 'lxis', default=0)

            # number of dimensionless thickess coordinates in tau
            self.ntaus = integer_or_blank(card, 6, 'ntaus', default=0)

            # ID of AEFACT that lists thickness ratios (t/c)
            self.ltaus = integer_or_blank(card, 7, 'ltaus', default=0)

            # ca/c - control surface chord / strip chord
            caoci = []
            j = 0
            for n, i in enumerate(range(9, len(card))):
                ca = double(card, i, 'ca/ci_%i' % (n+1))
                caoci.append(ca)
            self.caoci = array(caoci, dtype='float64')
        else:
            raise NotImplementedError('data')

    @property
    def lxis_id(self):
        return self.lxis if isinstance(self.lxis, integer_types) else self.lxis.aid

    @property
    def ltaus_id(self):
        return self.ltaus if isinstance(self.ltaus, integer_types) else self.ltaus.aid

    def cross_reference(self, model):
        self.ltaus = model.AEFACT(self.lxis_id)
        self.ltaus = model.AEFACT(self.ltaus_id)

    def raw_fields(self):
        list_fields = ['PAERO5', self.pid, self.nalpha, self.lalpha, self.nxis,
                       self.lxis_id, self.ntaus, self.ltaus_id] + self.caoci
        return list_fields

    #def integrals(self):
        ## chord location
        #x = self.lxis.Di

        ## thickness
        #y = self.ltaus.Di

        ## slope of airfoil semi-thickness
        #yp = derivative1(y/2, x)

        ## x hinge
        #for xh in self.caoci:
            #I1 = integrate(yp, x, 0., 1.)
            #I2 = integrate(x * yp, x, 0., 1.)
            #I3 = integrate(x**2*yp, x, 0., 1.)
            #I4 = integrate(yp**2, x, 0., 1.)
            #I5 = integrate(x**2 * yp**2, x, 0., 1.)

            #J1 = integrate(yp, x, xh, 1.)
            #J2 = integrate(x * yp, x, xh, 1.)
            #J3 = integrate(x**2*yp, x, xh, 1.)
            #J4 = integrate(yp**2, x, xh, 1.)
            #J5 = integrate(x**2 * yp**2, x, xh, 1.)

        #return(I1, I2, I3, I4, I5,
               #J1, J2, J3, J4, J5)


class FLFACT(BaseCard):
    """
    +--------+-----+----+------+----+----+----+----+----+
    |   1    |  2  |  3 |   4  | 5  | 6  | 7  | 8  | 9  |
    +--------+-----+----+------+----+----+----+----+----+
    | FLFACT | SID | F1 | F2   | F3 | F4 | F5 | F6 | F7 |
    +--------+-----+----+------+----+----+----+----+----+
    |        | F8  | F9 | etc. |    |    |    |    |    |
    +--------+-----+----+------+----+----+----+----+----+

    +--------+-----+----+------+-----+
    |   1    |  2  |  3 |   4  | 5   |
    +--------+-----+----+------+-----+
    | FLFACT | 97  | .3 |.7    | 3.5 |
    +--------+-----+----+------+-----+

    # delta quantity approach

    +--------+-----+-------+------+-------+----+--------+
    |   1    |  2  |  3    |   4  |   5   | 6  |     7  |
    +--------+-----+-------+------+-------+----+--------+
    | FLFACT | SID | F1    | THRU | FNF   | NF | FMID   |
    +--------+-----+-------+------+-------+----+--------+
    | FLFACT | 201 | 0.200 | THRU | 0.100 | 11 | 0.1333 |
    +--------+-----+-------+------+-------+----+--------+
    """
    type = 'FLFACT'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            assert len(card) > 2, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
            field3 = double_string_or_blank(card, 3, 'THRU')
            if field3 == 'THRU':
                f1 = double(card, 2, 'f1')
                fnf = double(card, 4, 'fnf')
                nf = double(card, 5, 'nf')
                fmid = double(card, 6, 'fmid')
                assert len(card) == 7, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
                i = linspace(0, nf, nf, endpoint=False)
                self.factors = ((f1*(fnf - fmid) * (nf - 1) + fnf * (fmid - f1) * i) /
                                   ((fnf - fmid) * (nf - 1) +       (fmid - f1) * i))
            else:
                self.factors = fields(double, card, 'factors', i=2, j=len(card))

            if len(self.factors) > 1 and self.factors[1] == 'THRU':
                msg = 'embedded THRUs not supported yet on FLFACT card\n'
                raise NotImplementedError(msg)
                #(a,thru,b,n,dn) = factors
                #for i in range(
        else:
            self.sid = data[0]
            self.factors = data[1:]

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the FLFACT object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['FLFACT', self.sid] + self.factors
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class FLUTTER(BaseCard):
    """
    Defines data needed to perform flutter analysis.

    +---------+-----+--------+------+------+-------+-------+-------------+------+
    |    1    |  2  |   3    |  4   |  5   |   6   |   7   |      8      |  9   |
    +---------+-----+--------+------+------+-------+-------+-------------+------+
    | FLUTTER | SID | METHOD | DENS | MACH | RFREQ | IMETH | NVALUE/OMAX | EPS  |
    +---------+-----+--------+------+------+-------+-------+-------------+------+
    | FLUTTER | 19  | K      | 119  | 219  | 319   |     S | 5           | 1.-4 |
    +---------+-----+--------+------+------+-------+-------+-------------+------+
    """
    type = 'FLUTTER'
    _field_map = {
        1: 'sid', 2:'method', 3:'density', 4:'mach', 5:'rfreq_vel', 6:'imethod',
        8:'epsilon',
    }
    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the FLUTTER card

        :param self:  the FLUTTER object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        if n == 7:
            if self.method in ['K', 'KE']:
                value = self.nValue
            elif self.method in ['PKS', 'PKNLS']:
                value = self.omax
            else:
                value = self.nValue
            return value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the FLUTTER card

        :param self:  the FLUTTER object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        if n == 7:
            if self.method in ['K', 'KE']:
                self.nValue = value
            elif self.method in ['PKS', 'PKNLS']:
                self.omax = value
            else:
                self.nValue = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.method = string(card, 2, 'method')
            self.density = integer(card, 3, 'density')
            self.mach = integer(card, 4, 'mach')
            self.rfreq_vel = integer(card, 5, 'rfreq_vel')

            if self.method in ['K', 'KE']:
                self.imethod = string_or_blank(card, 6, 'imethod', 'L')
                self.nValue = integer_or_blank(card, 7, 'nValue')
                self.omax = None
                assert self.imethod in ['L', 'S'], 'imethod = %s' % self.imethod
            elif self.method in ['PKS', 'PKNLS']:
                self.imethod = None
                self.nValue = None
                self.omax = double_or_blank(card, 7, 'omax')
            else:
                self.nValue = integer_or_blank(card, 7, 'nValue')
                self.omax = None
                self.imethod = None

            self.epsilon = double_or_blank(card, 8, 'epsilon')  # not defined in QRG
            assert len(card) <= 9, 'len(FLUTTER card) = %i' % len(card)

        else:
            assert len(data) == 8, 'FLUTTER = %s' % data
            self.sid = data[0]
            self.method = data[1]
            self.density = data[2]
            self.mach = data[3]
            self.rfreq_vel = data[4]
            self.method = data[5]
            self.imethod = data[6]
            self.nValue = data[7]
            self.omax = data[8]
            self.epsilon = None
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

        assert self.method in ['K', 'PK', 'PKNL', 'PKS', 'PKNLS', 'KE'], 'method = %s' % self.method

    def cross_reference(self, model):
        """
        Cross links the card

        :param self:   the FLUTTER object pointer
        :param model:  the BDF object
        :type model:   BDF()
        """
        msg = ' which is required by FLUTTER sid=%s' % self.sid
        self.density = model.FLFACT(self.density, msg=msg)
        self.mach = model.FLFACT(self.density, msg=msg)
        self.rfreq_val = model.FLFACT(self.density, msg=msg)

    def get_density(self):
        if isinstance(self.density, integer_types):
            return self.density
        return self.density.sid

    def get_mach(self):
        if isinstance(self.mach, integer_types):
            return self.mach
        return self.mach.sid

    def get_rfreq_vel(self):
        if isinstance(self.rfreq_vel, integer_types):
            return self.rfreq_vel
        return self.rfreq_vel.sid

    def _get_raw_nvalue_omax(self):
        if self.method in ['K', 'KE']:
            #assert self.imethod in ['L', 'S'], 'imethod = %s' % self.imethod
            return(self.imethod, self.nValue)
        elif self.method in ['PKS', 'PKNLS']:
            return(self.imethod, self.omax)
        else:
            return(self.imethod, self.nValue)

    def _repr_nvalue_omax(self):
        if self.method in ['K', 'KE']:
            imethod = set_blank_if_default(self.imethod, 'L')
            #assert self.imethod in ['L', 'S'], 'imethod = %s' % self.imethods
            return (imethod, self.nValue)
        elif self.method in ['PKS', 'PKNLS']:
            return(self.imethod, self.omax)
        else:
            return(self.imethod, self.nValue)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the FLUTTER object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        (imethod, nValue) = self._get_raw_nvalue_omax()
        list_fields = ['FLUTTER', self.sid, self.method, self.get_density(),
                       self.get_mach(), self.get_rfreq_vel(), imethod, nValue, self.epsilon]
        return list_fields

    #def repr_fields(self):
        #(imethod, nValue) = self._reprNValueOMax()
        #list_fields = ['FLUTTER', self.sid, self.method, self.get_density(), self.get_mach(),
        #          self.get_rfreq_vel(), imethod, nValue, self.epsilon]
        #return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class GUST(BaseCard):
    """
    Defines a stationary vertical gust for use in aeroelastic response
    analysis.

    +------+-----+-------+-----+-----+------+
    |   1  |  2  |   3   |  4  |  5  |  6   |
    +------+-----+-------+-----+-----+------+
    | GUST | SID | DLOAD | WG  | X0  | V    |
    +------+-----+-------+-----+-----+------+
    | GUST | 133 | 61    | 1.0 | 0.  | 1.+4 |
    +------+-----+-------+-----+-----+------+
    """
    type = 'GUST'
    _field_map = {
        1: 'sid', 2:'dload', 3:'wg', 4:'x0', 5:'V',
    }

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.dload = integer(card, 2, 'dload')
            self.wg = double(card, 3, 'wg')
            self.x0 = double(card, 4, 'x0')
            self.V = double_or_blank(card, 4, 'V')
            assert len(card) <= 6, 'len(GUST card) = %i' % len(card)
        else:
            self.sid = data[0]
            self.dload = data[1]
            self.wg = data[2]
            self.x0 = data[3]
            self.V = data[4]
            assert len(data) == 5, 'data = %s' % data

    #def Angle(self):
        #angle = self.wg*self.t*(t-(x-self.x0)/self.V) # T is the tabular
        #return angle

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the GUST object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['GUST', self.sid, self.dload, self.wg, self.x0, self.V]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class MKAERO1(BaseCard):
    """
    Provides a table of Mach numbers (m) and reduced frequencies (k) for
    aerodynamic matrix calculation.

    +---------+----+----+----+----+----+----+----+----+
    |    1    |  2 | 3  |  4 | 5  | 6  | 7  | 8  | 9  |
    +---------+----+----+----+----+----+----+----+----+
    | MKAERO1 | m1 | m2 | m3 | m4 | m5 | m6 | m7 | m8 |
    +---------+----+----+----+----+----+----+----+----+
    |         | k1 | k2 | k3 | k4 | k5 | k6 | k7 | k8 |
    +---------+----+----+----+----+----+----+----+----+
    """
    type = 'MKAERO1'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            fields = [interpret_value(field) for field in card[1:]]
            nfields = len(fields) - 8
            self.machs = []
            self.rFreqs = []
            for i in range(1, 1 + nfields):
                self.machs.append(double_or_blank(card, i, 'mach'))
                self.rFreqs.append(double_or_blank(card, i + 8, 'rFreq'))
            self.machs = wipe_empty_fields(self.machs)
            self.v = wipe_empty_fields(self.rFreqs)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

        #print("machs  = ",self.machs)
        #print("rFreqs = ",self.rFreqs)

    def addFreqs(self, mkaero):
        self.getMach_rFreqs()
        for m in mkaero.machs:
            self.machs.append(m)
        for f in mkaero.rFreqs:
            self.rFreqs.append(f)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the MKAERO1 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        #list_fields = ['MKAERO1']
        #for (i, mach, rfreq) in zip(count(), self.machs, self.rFreqs):
        #    list_fields += [mach, rfreq]
        machs = [None] * 8
        freqs = [None] * 8
        for i, mach in enumerate(self.machs):
            machs[i] = mach
        for i, freq in enumerate(self.rFreqs):
            freqs[i] = freq
        list_fields = ['MKAERO1'] + machs + freqs
        return list_fields

    def getMach_rFreqs(self):
        return self.machs, self.rFreqs

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class MKAERO2(BaseCard):
    """
    Provides a table of Mach numbers (m) and reduced frequencies (k) for
    aerodynamic matrix calculation.

    +---------+----+----+----+----+----+----+----+----+
    |    1    | 2  | 3  | 4  | 5  | 6  | 7  | 8  | 9  |
    +---------+----+----+----+----+----+----+----+----+
    | MKAERO2 | m1 | k1 | m2 | k2 | m3 | k3 | m4 | k4 |
    +---------+----+----+----+----+----+----+----+----+
    """
    type = 'MKAERO2'

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            fields = card.fields(1)
            nfields = len(fields)
            self.machs = []
            self.rFreqs = []
            for i in range(1, 1 + nfields, 2):
                self.machs.append(double(card, i, 'mach'))
                self.rFreqs.append(double(card, i + 1, 'rFreq'))
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def addFreqs(self, mkaero):
        self.getMach_rFreqs()
        for m in mkaero.machs:
            self.machs.append(m)
        for f in mkaero.rFreqs:
            self.rFreqs.append(f)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the MKAERO2 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['MKAERO2']
        for (i, mach, rfreq) in zip(count(), self.machs, self.rFreqs):
            list_fields += [mach, rfreq]
        return list_fields

    def getMach_rFreqs(self):
        return self.machs, self.rFreqs

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class MONPNT1(BaseCard):
    type = 'MONPNT1'
    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.name = string(card, 1, 'name')

            label_fields = [labeli for labeli in card[2:8] if labeli is not None]
            self.label = ''.join(label_fields).strip()
            assert len(self.label) <= 56, self.label

            self.axes = components(card, 9, 'axes')
            self.comp = string(card, 10, 'comp')
            self.cid = integer_or_blank(card, 11, 'cid', 0)
            self.xyz = [
                double_or_blank(card, 12, 'x', default=0.0),
                double_or_blank(card, 13, 'y', default=0.0),
                double_or_blank(card, 14, 'z', default=0.0),
            ]
        else:
            raise NotImplementedError(card)

    def raw_fields(self):
        list_fields = [
            'MONPNT1', self.name, self.label.strip(), self.axes, self.comp, self.cid,
        ] + self.xyz
        return list_fields

    def write_card(self, size=8, is_double=False):
        cid = self.cid
        x, y, z = self.xyz
        msg = 'MONPNT1 %-8s%s\n' % (self.name, self.label)
        msg += '        %-8s%-8s%-8s%-8s%-8s%-8s\n' % (
            self.axes, self.comp, cid,
            print_float_8(x), print_float_8(y), print_float_8(z))
        card = self.repr_fields()
        return self.comment + msg


class PAERO1(BaseCard):
    """
    Defines associated bodies for the panels in the Doublet-Lattice method.

    +--------+-----+----+----+----+----+----+----+
    | PAERO1 | PID | B1 | B2 | B3 | B4 | B5 | B6 |
    +--------+-----+----+----+----+----+----+----+
    """
    type = 'PAERO1'
    _field_map = {1: 'pid'}

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the PAERO1 card

        :param self:  the PAERO1 object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        return self.Bi[n - 1]

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the PAERO1 card

        :param self:  the PAERO1 object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        self.Bi[n - 1] = value

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            self.pid = integer(card, 1, 'pid')
            Bi = [interpret_value(field) for field in card[2:]]
            self.Bi = []

            for bi in Bi:
                if isinstance(bi, int) and bi >= 0:
                    self.Bi.append(bi)
                elif bi is not None:
                    raise RuntimeError('invalid Bi value on PAERO1 bi=|%r|' % (bi))
                #else:
                #    pass
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def cross_reference(self, model):
        pass

    def Bodies(self):
        return self.Bi

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the PAERO1 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['PAERO1', self.pid] + self.Bi
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class PAERO2(BaseCard):
    """
    Defines the cross-sectional properties of aerodynamic bodies.

    +--------+------+--------+-------+------+------+------+------+------+
    | PAERO2 | PID  | ORIENT | WIDTH |  AR  | LRSB | LRIB | LTH1 | LTH2 |
    +--------+------+--------+-------+------+------+------+------+------+
    | THI1   | THN1 |  THI2  |  THN2 | THI3 | THN3 |      |      |      |
    +--------+------+--------+-------+------+------+------+------+------+
    """
    type = 'PAERO2'
    _field_map = {
        1: 'pid', 2:'orient', 3:'width', 4:'AR', 5:'lrsb', 6:'lrib',
        7: 'lth1', 8:'lth2',
    }

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the PAERO2 card

        :param self:  the PAERO2 object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        nnew = n - 8
        spot = nnew // 2
        i = nnew % 2
        if i == 0:
            value = self.thi[spot]
        else:
            value = self.thn[spot]
        return value

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the PAERO2 card

        :param self:  the PAERO2 object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        nnew = n - 8
        spot = nnew // 2
        i = nnew % 2
        if i == 0:
            self.thi[spot] = value
        else:
            self.thn[spot] = value

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Property identification number. (Integer > 0)
            self.pid = integer(card, 1, 'pid')
            #: Orientation flag. Type of motion allowed for bodies. Refers to
            #: the aerodynamic coordinate system of ACSID. See AERO entry.
            #: (Character = 'Z', 'Y', or 'ZY')
            self.orient = string(card, 2, 'orient')
            #: Reference half-width of body and the width of the constant width
            #: interference tube. (Real > 0.0)
            self.width = double(card, 3, 'width')
            #: Aspect ratio of the interference tube (height/width). float>0.
            self.AR = double(card, 4, 'AR')
            #: Identification number of an AEFACT entry containing a list of
            #: slender body half-widths at the end points of the slender body
            #: elements. If blank, the value of WIDTH will be used.
            #: (Integer > 0 or blank)
            self.lrsb = integer_or_blank(card, 5, 'lrsb')
            #: Identification number of an AEFACT entry containing a list of
            #: slender body half-widths at the end points of the interference
            #: elements. If blank, the value of WIDTH will be used.
            #: (Integer > 0 or blank)
            self.lrib = integer_or_blank(card, 6, 'lrib')
            #: dentification number of AEFACT entries for defining ? arrays for
            #: interference calculations. (Integer >= 0)
            self.lth1 = integer_or_blank(card, 7, 'lth1')
            self.lth2 = integer_or_blank(card, 8, 'lth2')
            self.thi = []
            self.thn = []
            fields = [interpret_value(field) for field in card[9:]]
            nfields = len(fields)
            for i in range(9, 9 + nfields, 2):
                self.thi.append(integer(card, i, 'lth'))
                self.thn.append(integer(card, i + 1, 'thn'))
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def cross_reference(self, model):
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the PAERO2 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['PAERO2', self.pid, self.orient, self.width,
                       self.AR, self.lrsb, self.lrib, self.lth1, self.lth2]
        for (thi, thn) in zip(self.thi, self.thn):
            list_fields += [thi, thn]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PAERO3(BaseCard):
    """
    Defines the number of Mach boxes in the flow direction and the location of cranks and
    control surfaces of a Mach box lifting surface.
    """
    type = 'PAERO3'
    _field_map = {
        1: 'pid', 2:'orient', 3:'width', 4:'AR',
    }

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the PAERO3 card

        :param self:  the PAERO3 object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        nnew = n - 6
        if nnew < 0:
            raise RuntimeError('field n=%i on PAERO3 is invalid' % n)
        spot = nnew // 2
        i = nnew % 2
        if i == 0:
            value = self.x[spot]
        else:
            value = self.y[spot]
        return value

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the PAERO3 card

        :param self:  the PAERO3 object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        nnew = n - 6
        if nnew < 0:
            raise RuntimeError('field n=%i on PAERO3 is invalid' % n)
        spot = nnew // 2
        i = nnew % 2
        if i == 0:
            self.x[spot] = value
        else:
            self.y[spot] = value

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Property identification number. (Integer > 0)
            self.pid = integer(card, 1, 'pid')
            self.nbox = integer(card, 2, 'nbox')
            self.ncontrol_surfaces = integer(card, 3, 'ncontrol_surfaces')
            self.x = []
            self.y = []
            nfields = card.nfields

            j = 0
            for i in range(6, nfields, 2):
                x = double(card, i, 'x%i' % j)
                y = double(card, i + 1, 'y%i' % j)
                self.x.append(x)
                self.y.append(y)
                j += 1
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def cross_reference(self, model):
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the PAERO3 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['PAERO3', self.pid, self.nbox, self.ncontrol_surfaces, None]
        for (x, y) in zip(self.x, self.y):
            list_fields += [x, y]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class Spline(BaseCard):
    def __init__(self, card, data):
        pass


class SPLINE1(Spline):
    """
    Surface Spline Methods
    Defines a surface spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by regular arrays of
    aerodynamic points.

      +---------+-------+-------+------+------+------+----+------+-------+
      | SPLINE1 | EID   | CAERO | BOX1 | BOX2 | SETG | DZ | METH | USAGE |
      +---------+-------+-------+------+------+------+----+------+-------+
      | NELEM   | MELEM |       |      |      |      |    |      |       |
      +---------+-------+-------+------+------+------+----+------+-------+
      | SPLINE1 |   3   |  111  | 115  | 122  |  14  | 0. |      |       |
      +---------+-------+-------+------+------+------+----+------+-------+
    """
    type = 'SPLINE1'
    _field_map = {
        1: 'eid', 2:'caero', 3:'box1', 4:'box2', 5:'setg', 6:'dz',
        7: 'method', 8:'usage', 9:'nelements', 10:'melements',
    }

    def __init__(self, card=None, data=None, comment=''):
        Spline.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.caero = integer(card, 2, 'caero')
            self.box1 = integer(card, 3, 'box1')
            self.box2 = integer(card, 4, 'box2')
            self.setg = integer(card, 5, 'setg')
            self.dz = double_or_blank(card, 6, 'dz', 0.0)
            self.method = string_or_blank(card, 7, 'method', 'IPS')
            self.usage = string_or_blank(card, 8, 'usage', 'BOTH')
            self.nelements = integer_or_blank(card, 9, 'nelements', 10)
            self.melements = integer_or_blank(card, 10, 'melements', 10)
            assert self.nelements > 0, 'nelements = %s' % self.nelements
            assert self.melements > 0, 'melements = %s' % self.melements
            assert len(card) <= 11, 'len(SPLINE1 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.caero = data[1]
            self.box1 = data[2]
            self.box2 = data[3]
            self.setg = data[4]
            self.dz = data[5]
            self.method = data[6]
            self.usage = data[7]
            self.nelements = data[8]
            self.melements = data[9]
            assert len(data) == 10, 'data = %s' % data

        assert self.box2 >= self.box1, 'box1=%s box2=%s' % (self.box1, self.box2)
        assert self.method in ['IPS', 'TPS', 'FPS'], 'method = %s' % self.method
        assert self.usage in ['FORCE', 'DISP', 'BOTH'], 'usage = %s' % self.usage

    def CAero(self):
        if isinstance(self.caero, integer_types):
            return self.caero
        return self.caero.eid

    def Set(self):
        if isinstance(self.setg, integer_types):
            return self.setg
        return self.setg.sid

    def cross_reference(self, model):
        """
        Cross links the card

        :param self:   the SPLINE1 object pointer
        :param model:  the BDF object
        :type model:   BDF()
        """
        msg = ' which is required by SPLINE1 eid=%s' % self.eid
        self.caero = model.CAero(self.CAero(), msg=msg)
        self.setg = model.Set(self.Set(), msg=msg)
        self.setg.cross_reference(model, 'Node')

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the SPLINE1 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['SPLINE1', self.eid, self.CAero(), self.box1, self.box2,
                       self.Set(), self.dz, self.method, self.usage, self.nelements,
                       self.melements]
        return list_fields

    def repr_fields(self):
        dz = set_blank_if_default(self.dz, 0.)
        method = set_blank_if_default(self.method, 'IPS')
        usage = set_blank_if_default(self.usage, 'BOTH')
        nelements = set_blank_if_default(self.nelements, 10)
        melements = set_blank_if_default(self.melements, 10)

        list_fields = ['SPLINE1', self.eid, self.CAero(), self.box1, self.box2,
                       self.Set(), dz, method, usage, nelements, melements]
        list_fields = wipe_empty_fields(list_fields)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SPLINE2(Spline):
    """
    Linear Spline
    Defines a surface spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by regular arrays of
    aerodynamic points.

      +---------+------+-------+-------+-------+------+----+------+-----+
      | SPLINE2 | EID  | CAERO |  ID1  |  ID2  | SETG | DZ | DTOR | CID |
      +---------+------+-------+-------+-------+------+----+------+-----+
      |         | DTHX | DTHY  | None  | USAGE |      |    |      |     |
      +---------+------+-------+-------+-------+------+----+------+-----+

      +---------+------+-------+-------+-------+------+----+------+-----+
      | SPLINE2 |   5  |   8   |  12   | 24    | 60   | 0. | 1.0  |  3  |
      +---------+------+-------+-------+-------+------+----+------+-----+
      |         |  1.  |       |       |       |      |    |      |     |
      +---------+------+-------+-------+-------+------+----+------+-----+
    """
    type = 'SPLINE2'
    _field_map = {
        1: 'eid', 2:'caero', 3:'id1', 4:'id2', 5:'setg', 6:'dz',
        7: 'dtor', 8:'cid', 9:'dthx', 10:'dthy',
    }

    def __init__(self, card=None, data=None, comment=''):
        Spline.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.caero = integer(card, 2, 'caero')
            self.id1 = integer(card, 3, 'id1')
            self.id2 = integer(card, 4, 'id2')
            self.setg = integer(card, 5, 'setg')
            self.dz = double_or_blank(card, 6, 'dz', 0.0)
            self.dtor = double_or_blank(card, 7, 'dtor', 1.0)
            self.cid = integer_or_blank(card, 8, 'cid', 0)
            self.dthx = double_or_blank(card, 9, 'dthx')
            self.dthy = double_or_blank(card, 10, 'dthy')
            assert self.id2 >= self.id1, 'id2=%s id1=%s' % (self.id2, self.id1)

            self.usage = string_or_blank(card, 12, 'usage', 'BOTH')
            assert len(card) <= 13, 'len(SPLINE2 card) = %i' % len(card)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def cross_reference(self, model):
        """
        Cross links the card

        :param self:   the SPLINE2 object pointer
        :param model:  the BDF object
        :type model:   BDF()
        """
        msg = ' which is required by SPLINE2 eid=%s' % self.eid
        self.cid = model.Coord(self.Cid(), msg=msg)
        self.caero = model.CAero(self.CAero(), msg=msg)
        self.setg = model.Set(self.Set(), msg=msg)
        self.setg.cross_reference(model, 'Node')

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid.cid

    def CAero(self):
        if isinstance(self.caero, integer_types):
            return self.caero
        return self.caero.eid

    def Set(self):
        if isinstance(self.setg, integer_types):
            return self.setg
        return self.setg.sid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the SPLINE2 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['SPLINE2', self.eid, self.CAero(), self.id1, self.id2,
                       self.Set(), self.dz, self.dtor, self.Cid(), self.dthx,
                       self.dthy, None, self.usage]
        return list_fields

    def repr_fields(self):
        dz = set_blank_if_default(self.dz, 0.)
        usage = set_blank_if_default(self.usage, 'BOTH')
        list_fields = ['SPLINE2', self.eid, self.CAero(), self.id1, self.id2,
                       self.Set(), dz, self.dtor, self.Cid(), self.dthx, self.dthy,
                       None, usage]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SPLINE3(Spline):
    type = 'SPLINE3'
    _field_map = {
        1: 'eid', 2:'caero', 3:'box_id', 5:'g1', 6:'c1',
        7: 'a1', 8:'usage',
    }

    def __init__(self, card=None, data=None, comment=''):
        Spline.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.caero = integer(card, 2, 'caero')
            self.box_idt = integer(card, 3, 'box_id')
            self.components = integer(card, 4, 'comp')
            assert self.components in [1, 2, 3, 4, 5, 6], self.components
            self.g1 = integer(card, 5, 'G1')
            self.c1 = integer(card, 5, 'C1')
            self.a1 = double(card, 5, 'A1')
            assert self.c1 in [0, 1, 2, 3, 4, 5, 6], self.components
            self.usage = string_or_blank(card, 8, 'usage', 'BOTH')
            assert self.usage in ['FORCE', 'DISP', 'BOTH'], self.usage

            nfields = len(card) - 1
            nrows = nfields // 8
            if nfields % 8:
                nrows += 1

            i = 1
            self.Gi = []
            self.ci = []
            self.ai = []
            for irow in range(1, nrows):
                j = 1 + nrows * 8
                gi = integer(card, j, 'Gi_' % i)
                ci = integer(card, j + 1, 'Ci_' % i)
                ai = double(card, j + 2, 'Ai_' % i)
                self.Gi.append(gi)
                self.ci.append(ci)
                self.ai.append(ai)
                if card[j + 3] or card[j + 4] or card[j + 5]:
                    i += 1
                    gi = integer(card, j, 'Gi_' % i)
                    ci = integer(card, j + 1, 'Ci_' % i)
                    ai = double(card, j + 2, 'Ai_' % i)
                    self.Gi.append(gi)
                    self.ci.append(ci)
                    self.ai.append(ai)
                i += 1
        else:
            raise NotImplementedError()

    def cross_reference(self, model):
        pass

    def raw_fields(self):
        list_fields = ['SPLINE3', self.eid, self.caero, self.box_id, self.g1, self.c1, self.a1, self.usage]
        for g, c, a in zip(self.Gi, self.ci, self.ai):
            list_fields += [g, c, a, None]
        return list_fields


class SPLINE4(Spline):
    """
    Surface Spline Methods
    Defines a curved surface spline for interpolating motion and/or forces for
    aeroelastic problems on general aerodynamic geometries using either the
    Infinite Plate, Thin Plate or Finite Plate splining method.

     +---------+-------+-------+--------+-----+------+----+------+-------+
     | SPLINE4 | EID   | CAERO | AELIST | --- | SETG | DZ | METH | USAGE |
     +---------+-------+-------+--------+-----+------+----+------+-------+
     | NELEM   | MELEM |       |        |     |      |    |      |       |
     +---------+-------+-------+--------+-----+------+----+------+-------+
     | SPLINE4 |   3   | 111   |   115  | --- |  14  | 0. | IPS  |       |
     +---------+-------+-------+--------+-----+------+----+------+-------+
    """
    type = 'SPLINE4'
    _field_map = {
        1: 'eid', 2:'caero', 3:'aelist', 5:'setg', 6:'dz',
        7: 'method', 8:'usage', 9:'nelements', 10:'melements',
    }

    def __init__(self, card=None, data=None, comment=''):
        Spline.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.caero = integer(card, 2, 'caero')
            self.aelist = integer(card, 3, 'aelist')
            # None
            self.setg = integer(card, 5, 'setg')
            self.dz = double_or_blank(card, 6, 'dz', 0.0)
            self.method = string_or_blank(card, 7, 'method', 'IPS')
            self.usage = string_or_blank(card, 8, 'usage', 'BOTH')
            self.nelements = integer_or_blank(card, 9, 'nelements', 10)
            self.melements = integer_or_blank(card, 10, 'melements', 10)
            assert len(card) <= 11, 'len(SPLINE4 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.caero = data[1]
            self.aelist = data[2]
            self.setg = data[3]
            self.dz = data[4]
            self.method = data[5]
            self.usage = data[6]
            self.nelements = data[7]
            self.melements = data[8]
            assert len(data) == 9, 'data = %s' % (data)

        assert self.method in ['IPS', 'TPS', 'FPS'], 'method = %s' % self.method
        assert self.usage in ['FORCE', 'DISP', 'BOTH'], 'uasge = %s' % self.usage

    def CAero(self):
        if isinstance(self.caero, integer_types):
            return self.caero
        return self.caero.eid

    def AEList(self):
        if isinstance(self.aelist, integer_types):
            return self.aelist
        return self.aelist.sid

    def Set(self):
        if isinstance(self.setg, integer_types):
            return self.setg
        return self.setg.sid

    def cross_reference(self, model):
        """
        Cross links the card

        :param self:   the SPLINE4 object pointer
        :param model:  the BDF object
        :type model:   BDF()
        """
        msg = ' which is required by SPLINE4 eid=%s' % self.eid
        self.caero = model.CAero(self.CAero(), msg=msg)
        self.setg = model.Set(self.Set(), msg=msg)
        self.aelist = model.AEList(self.aelist, msg=msg)
        self.setg.cross_reference(model, 'Node')

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the SPLINE4 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['SPLINE4', self.eid, self.CAero(), self.AEList(), None,
                       self.Set(), self.dz, self.method, self.usage, self.nelements,
                       self.melements]
        return list_fields

    def repr_fields(self):
        dz = set_blank_if_default(self.dz, 0.)
        method = set_blank_if_default(self.method, 'IPS')
        usage = set_blank_if_default(self.usage, 'BOTH')
        nelements = set_blank_if_default(self.nelements, 10)
        melements = set_blank_if_default(self.melements, 10)

        list_fields = ['SPLINE4', self.eid, self.CAero(), self.AEList(), None,
                       self.Set(), dz, method, usage, nelements, melements]
        list_fields = wipe_empty_fields(list_fields)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class SPLINE5(Spline):
    """
    Linear Spline
    Defines a 1D beam spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by irregular arrays
    of aerodynamic points. The interpolating beam supports axial rotation and
    bending in the yz-plane.

    +---------+------+-------+--------+-------+------+----+------+-----+
    | SPLINE5 | EID  | CAERO | AELIST | ---   | SETG | DZ | DTOR | CID |
    +---------+------+-------+--------+-------+------+----+------+-----+
    |         | DTHX | DTHY  | ---    | USAGE |      |    |      |     |
    +---------+------+-------+--------+-------+------+----+------+-----+
    """
    type = 'SPLINE5'
    _field_map = {
        1: 'eid', 2:'caero', 3:'aelist', 5:'setg', 6:'dz',
        7: 'dtor', 8:'cid', 9:'thx', 10:'thy', 12:'usage',
    }

    def __init__(self, card=None, data=None, comment=''):
        Spline.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.caero = integer(card, 2, 'caero')
            self.aelist = integer(card, 3, 'aelist')
            # None
            self.setg = integer(card, 5, 'setq')
            self.dz = double_or_blank(card, 6, 'dz', 0.0)
            self.dtor = double_or_blank(card, 7, 'dtor', 1.0)
            self.cid = integer_or_blank(card, 8, 'cid', 0)
            self.thx = double(card, 9, 'thx')
            self.thy = double(card, 10, 'thy')

            self.usage = string_or_blank(card, 12, 'usage', 'BOTH')
            assert len(card) <= 13, 'len(SPLINE5 card) = %i' % len(card)
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid.cid

    def CAero(self):
        if isinstance(self.caero, integer_types):
            return self.caero
        return self.caero.eid

    def AEList(self):
        if isinstance(self.aelist, integer_types):
            return self.aelist
        return self.aelist.aelist

    def Set(self):
        if isinstance(self.setg, integer_types):
            return self.setg
        return self.setg.sid

    def cross_reference(self, model):
        """
        Cross links the card

        :param self:   the SPLINE5 object pointer
        :param model:  the BDF object
        :type model:   BDF()
        """
        msg = ' which is required by SPLINE5 eid=%s' % self.eid
        self.cid = model.Coord(self.Cid(), msg=msg)
        self.caero = model.CAero(self.CAero(), msg=msg)
        self.setg = model.Set(self.Set(), msg=msg)
        self.aelist = model.AEList(self.AEList(), msg=msg)
        self.setg.cross_reference(model, 'Node')

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the SPLINE5 object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['SPLINE5', self.eid, self.CAero(), self.AEList(), None,
                       self.Set(), self.dz, self.dtor, self.Cid(), self.thx,
                       self.thy, None, self.usage]
        return list_fields

    def repr_fields(self):
        dz = set_blank_if_default(self.dz, 0.)
        usage = set_blank_if_default(self.usage, 'BOTH')
        list_fields = ['SPLINE5', self.eid, self.CAero(), self.AEList(), None,
                       self.Set(), dz, self.dtor, self.Cid(), self.thx, self.thy,
                       None, usage]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class TRIM(BaseCard):
    type = 'TRIM'
    _field_map = {
        1: 'sid', 2:'mach', 3:'q', 8:'aeqr',
    }

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the TRIM card

        :param self:  the TRIM object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        ni = 4
        for (i, label, ux) in zip(count(), self.labels, self.uxs):
            if n == ni:
                value = self.labels[i]
                return value
            elif n + 1 == ni:
                value = self.uxs[i]
                return value

            #list_fields += [label, ux]
            if i == 1:
                #list_fields += [self.aeqr]
                ni += 1
        raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the TRIM card

        :param self:  the TRIM object pointer
        :param n:     the field number to update
        :type n:      int
        :param value: the value for the appropriate field
        :type field:  varies
        """
        ni = 4
        for (i, label, ux) in zip(count(), self.labels, self.uxs):
            if n == ni:
                self.labels[i] = value
                return
            elif n + 1 == ni:
                self.uxs[i] = value
                return

            #list_fields += [label, ux]
            if i == 1:
                #list_fields += [self.aeqr]
                ni += 1
        raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        if card:
            #: Trim set identification number. (Integer > 0)
            self.sid = integer(card, 1, 'sid')
            #: Mach number. (Real > 0.0 and != 1.0)
            self.mach = double(card, 2, 'mach')
            assert self.mach >= 0.0 and self.mach != 1.0, 'mach = %s' % self.mach
            #: Dynamic pressure. (Real > 0.0)
            self.q = double(card, 3, 'q')
            assert self.q > 0.0, 'q=%s' % self.q
            #: The label identifying aerodynamic trim variables defined on an
            #: AESTAT or AESURF entry.
            self.labels = []
            #: The magnitude of the aerodynamic extra point degree-of-freedom.
            #: (Real)
            self.uxs = []
            #: Flag to request a rigid trim analysis (Real > 0.0 and < 1.0;
            #: Default = 1.0. A value of 0.0 provides a rigid trim analysis,
            #: not supported

            label = string_or_blank(card, 4, 'label1')
            if label:
                ux = double(card, 5, 'ux1')
                self.uxs.append(ux)
                self.labels.append(label)

            label = string_or_blank(card, 6, 'label2')
            if label:
                ux = double(card, 7, 'ux1')
                self.uxs.append(ux)
                self.labels.append(label)
            self.aeqr = double_or_blank(card, 8, 'aeqr')

            i = 9
            n = 3
            while i < len(card):
                label = string(card, i, 'label%i' % n)
                ux = double(card, i + 1, 'ux%i' % n)
                self.labels.append(label)
                self.uxs.append(ux)
                i += 2
        else:
            msg = '%s has not implemented data parsing' % self.type
            raise NotImplementedError(msg)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        :param self:
          the TRIM object pointer
        :returns fields:
          the fields that define the card
        :type fields:
          LIST
        """
        list_fields = ['TRIM', self.sid, self.mach, self.q]
        for (i, label, ux) in zip(count(), self.labels, self.uxs):
            list_fields += [label, ux]
            if i == 1:
                list_fields += [self.aeqr]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)
