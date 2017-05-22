# coding: utf-8
# pylint: disable=R0902,R0904,R0914,C0302,C0111
"""
All aero cards are defined in this file.  This includes:

 * AECOMP
 * AEFACT
 * AELINK
 * AELIST
 * AEPARM
 * AESTAT
 * AESURF / AESURFS
 * AERO / AEROS
 * CSSCHD
 * CAERO1 / CAERO2 / CAERO3 / CAERO4 / CAERO5
 * DIVERG
 * FLFACT
 * FLUTTER
 * GUST
 * MKAERO1 / MKAERO2
 * PAERO1 / PAERO2 / PAERO3 / PAERO4 / PAERO5
 * SPLINE1 / SPLINE2 / SPLINE3 / SPLINE4 / SPLINE5
 * MONPNT1 / MONPNT2 / MONPNT3

All cards are BaseCard objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import count
import math
from six.moves import zip, range
from six import string_types

import numpy as np

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8, print_float_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (
    fields, integer, integer_or_blank, double, double_or_blank, string,
    string_or_blank, integer_or_string, double_string_or_blank,
    interpret_value, parse_components)
from pyNastran.bdf.cards.utils import wipe_empty_fields


class AECOMP(BaseCard):
    """
    Defines a component for use in monitor point definition or external splines.

    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    |   1    |   2   |    3     |   4   |   5   |   6   |   7   |   8   |   9   |
    +========+=======+==========+=======+=======+=======+=======+=======+=======+
    | AECOMP | NAME  | LISTTYPE | LIST1 | LIST2 | LIST3 | LIST4 | LIST5 | LIST6 |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    |        | LIST7 |  -etc.-  |       |       |       |       |       |       |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+

    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    | AECOMP | WING  |  AELIST  | 1001  | 1002  |       |       |       |       |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+

    Attributes
    ----------
    name : str
        The name.
    list_type : str
        {'SET1', 'AELIST', 'CAEROx'}
    lists : list[int]
        list of values of AECOMP lists
    """
    type = 'AECOMP'
    allowed_list_types = ['SET1', 'AELIST', 'CAERO']

    def __init__(self, name, list_type, lists, comment=''):
        """
        Creates an AECOMP card

        Parameters
        ----------
        name : str
            the name of the component
        list_type : str
            One of CAERO, AELIST or CMPID for aerodynamic components and
            SET1 for structural components. Aerodynamic components are
            defined on the aerodynamic ks-set mesh while the structural
            components are defined on the g-set mesh.
        lists : List[int, int, ...]
            The identification number of either SET1, AELIST or CAEROi
            entries that define the set of grid points that comprise
            the component
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.name = name
        self.list_type = list_type
        self.lists = lists

    def validate(self):
        if not self.list_type in ['SET1', 'AELIST', 'CAERO', 'CMPID']:
            msg = 'list_type=%r not in [SET1, AELIST, CAERO, CMPID]' % self.list_type
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AECOMP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        name = string(card, 1, 'name')
        list_type = string(card, 2, 'list_type')
        j = 1
        lists = []
        for i in range(3, len(card)):
            list_i = integer(card, i, '%s_%i' % (list_type, j))
            lists.append(list_i)
            j += 1
        return AECOMP(name, list_type, lists, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by AECOMP name=%r' % self.name
        #return
        if self.list_type == 'SET1':
            self.lists = [model.SET1(key, msg) for key in self.lists]
        elif self.list_type == 'AELIST':
            self.lists = [model.AELIST(key, msg) for key in self.lists]
        elif self.list_type == 'CAERO':
            self.lists = [model.CAero(key, msg) for key in self.lists]
        #elif self.list_type == 'CMPID':
            # AEQUAD4,/AETRIA3
        else:
            raise NotImplementedError(self.list_type)
        self.lists_ref = self.lists

    def uncross_reference(self):
        self.lists = self.get_lists()
        del self.lists_ref

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
        #elif self.list_type == 'CMPID':
            # AEQUAD4,/AETRIA3
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
    |   1    |   2 |  3 |    4   |  5  |  6 |  7 | 8  |  9 |
    +========+=====+====+========+=====+====+====+====+====+
    | AEFACT | SID | D1 |   D2   | D3  | D4 | D5 | D6 | D7 |
    +--------+-----+----+--------+-----+----+----+----+----+
    |        | D8  | D9 | -etc.- |     |    |    |    |    |
    +--------+-----+----+--------+-----+----+----+----+----+

    +--------+-----+----+--------+-----+
    | AEFACT | 97  |.3  |  0.7   | 1.0 |
    +--------+-----+----+--------+-----+

    TODO: Are these defined in percentages and thus,
          should they be normalized if they are not?
    """
    type = 'AEFACT'

    def __init__(self, sid, Di, comment=''):
        """
        Creates an AEFACT card, which defines the mach, dynamic_pressure,
        velocity, and reduced frequency for an FLUTTER card

        Used in flutter (145) and gust (146) analysis.

        Parameters
        ----------
        sid : int
            unique id
        Di : List[float, ..., float]
            list of:
             - machs
             - dynamic_pressures
             - velocities
             - reduced frequency
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        #: Set identification number. (Unique Integer > 0)
        self.sid = sid
        #: Number (float)
        self.Di = np.asarray(Di, dtype='float64')

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AEFACT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')

        Di = []
        for i in range(2, len(card)):
            di = double(card, i, 'Di_%i' % (i - 1))
            Di.append(di)
        return AEFACT(sid, Di, comment=comment)

    #def cross_reference(self, model):
        #pass

    #def uncross_reference(self):
        #pass

    @property
    def data(self):
        return self.Di

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card
        """
        list_fields = ['AEFACT', self.sid] + list(self.Di)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AELINK(BaseCard):
    r"""
    Defines relationships between or among AESTAT and AESURF entries, such
    that:

    .. math:: u^D + \Sigma_{i=1}^n C_i u_i^I = 0.0

    +--------+-------+-------+--------+----+-------+----+-------+----+
    |   1    |   2   |   3   |   4    |  5 |   6   |  7 |   8   |  9 |
    +========+=======+=======+========+====+=======+====+=======+====+
    | AELINK |  ID   | LABLD | LABL1  | C1 | LABL2 | C2 | LABL3 | C3 |
    +--------+-------+-------+--------+----+-------+----+-------+----+
    |        | LABL4 |  C4   |  etc.  |    |       |    |       |    |
    +--------+-------+-------+--------+----+-------+----+-------+----+

    +--------+-------+-------+-------+------+
    | AELINK |  10   | INBDA | OTBDA | -2.0 |
    +--------+-------+-------+-------+------+
    """
    type = 'AELINK'

    def __init__(self, id, label, independent_labels, Cis, comment=''):
        """
        Creates an AELINK card, which defines an equation linking
        AESTAT and AESURF cards

        Parameters
        ----------
        id : int
            unique id
        label : str
            name of the AESURF(???) card
        independent_labels : List[str, ..., str]
            name for the AESTAT(???) cards
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        #: defines the dependent variable name (string)
        self.label = label

        #: defines the independent variable name (string)
        self.independent_labels = independent_labels

        #: linking coefficient (real)
        self.Cis = Cis

        if isinstance(id, string_types):
            if id != 'ALWAYS':
                raise RuntimeError("The only valid ID that is a string is 'ALWAYS'")
            id = 0
        #: an ID=0 is applicable to the global subcase, ID=1 only subcase 1
        self.id = id

    def validate(self):
        if len(self.independent_labels) != len(self.Cis):
            msg = 'nlabels=%s nci=%s\nindependent_labels=%s Cis=%s\n%s' % (
                len(self.independent_labels), len(self.Cis),
                self.independent_labels, self.Cis, str(self))
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AELINK card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        id = integer_or_string(card, 1, 'ID')
        label = string(card, 2, 'label')
        independent_labels = []
        Cis = []

        list_fields = [interpret_value(field) for field in card[3:]]
        assert len(list_fields) % 2 == 0, 'list_fields=%s' % list_fields
        for i in range(0, len(list_fields), 2):
            independent_label = list_fields[i]
            Ci = list_fields[i + 1]
            independent_labels.append(independent_label)
            Cis.append(Ci)
        return AELINK(id, label, independent_labels, Cis, comment=comment)

    #def uncross_reference(self):
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        list_fields : List[int/float/str]
            the fields that define the card
        """
        list_fields = ['AELINK', self.id, self.label]
        for (ivar, ival) in zip(self.independent_labels, self.Cis):
            list_fields += [ivar, ival]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class AELIST(BaseCard):
    """
    Defines a list of aerodynamic elements to undergo the motion prescribed
    with the AESURF Bulk Data entry for static aeroelasticity.

    +---------+------+------+------+------+------+------+------+------+
    |    1    |   2  |   3  |  4   |   5  |   6  |  7   |  8   |  9   |
    +=========+======+======+======+======+======+======+======+======+
    |  AELIST |  SID |  E1  |  E2  |  E3  |  E4  |  E5  |  E6  |  E7  |
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

    def __init__(self, sid, elements, comment=''):
        """
        Creates an AELIST card, which defines the aero boxes for
        an AESURF/SPLINEx.

        Parameters
        ----------
        sid : int
            unique id
        elements : List[int, ..., int]
            list of box ids
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        #: Set identification number. (Integer > 0)
        self.sid = sid
        #: List of aerodynamic boxes generated by CAERO1 entries to define a
        #: surface. (Integer > 0 or 'THRU')
        self.elements = expand_thru(elements)
        self.elements.sort()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AELIST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        elements = fields(integer_or_string, card, 'eid', i=2, j=len(card))
        return AELIST(sid, elements, comment=comment)

    def cross_reference(self, model):
        pass

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass

    def clean_ids(self):
        self.elements = list(set(self.elements))
        self.elements.sort()

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card
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
    |    1   | 2  |    3   |   4   |
    +========+====+========+=======+
    | AEPARM | ID | LABEL  | UNITS |
    +--------+----+--------+-------+
    | AEPARM | 5  | THRUST | LBS   |
    +--------+----+--------+-------+
    """
    type = 'AEPARM'
    _field_map = {
        1: 'id', 2:'label', 3:'units'
    }

    def __init__(self, aeparm_id, label, units, comment=''):
        """
        Creates an AEPARM card, which defines a new trim variable.

        Parameters
        ----------
        id : int
            the unique id
        label : str
            the variable name
        units : str
            unused by Nastran
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.id = aeparm_id
        self.label = label
        self.units = units

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AEPARM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        aeparm_id = integer(card, 1, 'aeparm_id')
        label = string(card, 2, 'label')
        units = card.field(3)
        units = '' if units is None else units

        assert len(card) <= 4, 'len(AEPARM card) = %i\ncard=%s' % (len(card), card)
        return AEPARM(aeparm_id, label, units, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds an AEPARM card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        id = data[0]
        label = data[1]
        units = data[2]
        assert len(data) == 3, 'data = %s' % data
        return AEPARM(id, label, units, comment=comment)

    def cross_reference(self, model):
        pass

    #def uncross_reference(self):
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card
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
    |    1   |   2  |    3   |
    +========+======+========+
    | AESTAT |  ID  | LABEL  |
    +--------+------+--------+
    | AESTAT | 5001 | ANGLEA |
    +--------+------+--------+
    """
    type = 'AESTAT'

    _field_map = {
        1: 'id', 2:'label',
    }
    def __init__(self, aestat_id, label, comment=''):
        """
        Creates an AESTAT card, which is a variable to be used in a TRIM analysis

        Parameters
        ----------
        id : int
            unique id
        label : str
            name for the id
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.id = aestat_id
        self.label = label

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AESTAT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        aestat_id = integer(card, 1, 'ID')
        label = string(card, 2, 'label')
        assert len(card) <= 3, 'len(AESTAT card) = %i\ncard=%s' % (len(card), card)
        return AESTAT(aestat_id, label, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        aestat_id = data[0]
        label = data[1]
        assert len(data) == 2, 'data = %s' % data
        return AESTAT(aestat_id, label, comment=comment)

    #def cross_reference(self, model):
        #pass

    #def uncross_reference(self):
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/str]
            the fields that define the card
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
    |    1   |   2    |   3   |   4   |   5   |   6    |    7   |   8    |   9    |
    +========+========+=======+=======+=======+========+========+========+========+
    | AESURF |   ID   | LABEL | CID1  | ALID1 |  CID2  | ALID2  |  EFF   |  LDW   |
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

    def __init__(self, aesid, label, cid1, alid1, cid2=None, alid2=None, eff=1.0, ldw='LDW',
                 crefc=1.0, crefs=1.0, pllim=-np.pi/2., pulim=np.pi/2.,
                 hmllim=None, hmulim=None, tqllim=None, tqulim=None,
                 comment=''):
        """
        Creates an AESURF card, which defines a control surface

        Parameters
        ----------
        aesid : int
            controller number
        label : str
            controller name
        cid1 / cid2 : int / None
            coordinate system id for primary/secondary control surface
        alid1 / alid2 : int / None
            AELIST id for primary/secondary control surface
        eff : float; default=1.0
            Control surface effectiveness
        ldw : str; default='LDW'
            Linear downwash flag;  ['LDW', 'NODLW']
        crefc : float; default=1.0
            reference chord for the control surface
        crefs : float; default=1.0
            reference area for the control surface
        pllim / pulim : float; default=-pi/2 / pi/2
            Lower/Upper deflection limits for the control surface in radians
        hmllim / hmulim : float; default=None
            Lower/Upper hinge moment limits for the control surface in
            force-length units
        tqllim / tqulim : int; default=None
            Set identification numbers of TABLEDi entries that provide the
            lower/upper deflection limits for the control surface as a
            function of the dynamic pressure
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        #: Controller identification number
        self.aesid = aesid

        #: Controller name.
        self.label = label

        #: Identification number of a rectangular coordinate system with a
        #: y-axis that defines the hinge line of the control surface
        #: component.
        self.cid1 = cid1
        #: Identification of an AELIST Bulk Data entry that identifies all
        #: aerodynamic elements that make up the control surface
        #: component. (Integer > 0)
        self.alid1 = alid1

        self.cid2 = cid2
        self.alid2 = alid2

        #: Control surface effectiveness. See Remark 4. (Real != 0.0;
        #: Default=1.0)
        self.eff = eff
        #: Linear downwash flag. See Remark 2.
        #: (Character, one of LDW or NOLDW; Default=LDW).
        self.ldw = ldw
        #: Reference chord length for the control surface. (Real>0.0;
        #: Default=1.0)
        self.crefc = crefc
        #: Reference surface area for the control surface. (Real>0.0;
        #: Default=1.0)
        self.crefs = crefs

        #: Lower and upper deflection limits for the control surface in
        #: radians. (Real, Default = +/- pi/2)
        self.pllim = pllim
        self.pulim = pulim

        #: Lower and upper hinge moment limits for the control surface in
        #: force-length units. (Real, Default = no limit) -> 1e8
        self.hmllim = hmllim
        self.hmulim = hmulim

        #: Set identification numbers of TABLEDi entries that provide the
        #: lower and upper deflection limits for the control surface as a
        #: function of the dynamic pressure. (Integer>0, Default = no limit)
        self.tqllim = tqllim
        self.tqulim = tqulim
        self.tqllim_ref = None
        self.tqulim_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AESURF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        aesid = integer(card, 1, 'aesid')
        label = string(card, 2, 'label')

        cid1 = integer(card, 3, 'cid1')
        alid1 = integer(card, 4, 'alid1')

        cid2 = integer_or_blank(card, 5, 'cid2')
        alid2 = integer_or_blank(card, 6, 'alid2')

        eff = double_or_blank(card, 7, 'eff', 1.0)
        ldw = string_or_blank(card, 8, 'ldw', 'LDW')
        crefc = double_or_blank(card, 9, 'crefc', 1.0)
        crefs = double_or_blank(card, 10, 'crefs', 1.0)

        pllim = double_or_blank(card, 11, 'pllim', -np.pi / 2.)
        pulim = double_or_blank(card, 12, 'pulim', np.pi / 2.)

        hmllim = double_or_blank(card, 13, 'hmllim')
        hmulim = double_or_blank(card, 14, 'hmulim')
        tqllim = integer_or_blank(card, 15, 'tqllim')
        tqulim = integer_or_blank(card, 16, 'tqulim')
        assert len(card) <= 17, 'len(AESURF card) = %i\ncard=%s' % (len(card), card)
        return AESURF(aesid, label, cid1, alid1, cid2, alid2, eff, ldw,
                      crefc, crefs, pllim, pulim, hmllim, hmulim,
                      tqllim, tqulim, comment=comment)

    def Cid1(self):
        if isinstance(self.cid1, integer_types):
            return self.cid1
        return self.cid1_ref.cid

    def Cid2(self):
        if isinstance(self.cid2, integer_types) or self.cid2 is None:
            return self.cid2
        return self.cid2_ref.cid

    def AELIST_id1(self):
        if isinstance(self.alid1, integer_types):
            return self.alid1
        return self.alid1_ref.sid

    def AELIST_id2(self):
        if isinstance(self.alid2, integer_types) or self.alid2 is None:
            return self.alid2
        return self.alid2_ref.sid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.cid1 = model.Coord(self.Cid1())
        self.cid1_ref = self.cid1
        if self.cid2 is not None:
            self.cid2 = model.Coord(self.Cid2())
            self.cid2_ref = self.cid2
        self.alid1 = model.AELIST(self.AELIST_id1())
        self.alid1_ref = self.alid1
        if self.alid2:
            self.alid2 = model.AELIST(self.AELIST_id2())
            self.alid2_ref = self.alid2
        if self.tqllim is not None:
            self.tqllim_ref = model.TableD(self.tqllim)
        if self.tqulim is not None:
            self.tqulim_ref = model.TableD(self.tqulim)

    def safe_cross_reference(self, model):
        try:
            self.cid1 = model.Coord(self.Cid1())
            self.cid1_ref = self.cid1
        except KeyError:
            pass
        if self.cid2 is not None:
            try:
                self.cid2 = model.Coord(self.Cid2())
                self.cid2_ref = self.cid2
            except KeyError:
                pass
        try:
            self.alid1 = model.AELIST(self.AELIST_id1())
            self.alid1_ref = self.alid1
        except KeyError:
            pass
        if self.alid2:
            try:
                self.alid2 = model.AELIST(self.AELIST_id2())
                self.alid2_ref = self.alid2
            except KeyError:
                pass
        if self.tqllim is not None:
            try:
                self.tqllim_ref = model.TableD(self.tqllim)
            except KeyError:
                pass
        if self.tqulim is not None:
            try:
                self.tqulim_ref = model.TableD(self.tqulim)
            except KeyError:
                pass

    def uncross_reference(self):
        self.cid1 = self.Cid1()
        self.cid2 = self.Cid2()
        del self.cid1_ref
        if self.cid2:
            del self.cid2_ref

        self.alid1 = self.AELIST_id1()
        del self.alid1_ref

        self.alid2 = self.AELIST_id2()
        if self.alid2:
            del self.alid2_ref

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fieldsreset_camera[int/float/str]
            the fields that define the card
        """
        list_fields = ['AESURF', self.aesid, self.label, self.Cid1(), self.AELIST_id1(),
                       self.Cid2(), self.AELIST_id2(), self.eff, self.ldw,
                       self.crefc, self.crefs, self.pllim, self.pulim, self.hmllim,
                       self.hmulim, self.tqllim, self.tqulim]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card
        """
        eff = set_blank_if_default(self.eff, 1.0)
        ldw = set_blank_if_default(self.ldw, 'LDW')
        crefc = set_blank_if_default(self.crefc, 1.0)
        crefs = set_blank_if_default(self.crefs, 1.0)

        pllim = set_blank_if_default(self.pllim, -np.pi / 2.)
        pulim = set_blank_if_default(self.pulim, np.pi / 2.)

        list_fields = ['AESURF', self.aesid, self.label, self.Cid1(), self.AELIST_id1(),
                       self.Cid2(), self.AELIST_id2(), eff, ldw, crefc, crefs,
                       pllim, pulim, self.hmllim, self.hmulim, self.tqllim,
                       self.tqulim]
        return list_fields

    def write_card(self, size=8, is_double=False):
        """
        Writes the card with the specified width and precision

        Parameters
        ----------
        size : int (default=8)
            size of the field; {8, 16}
        is_double : bool (default=False)
            is this card double precision

        Returns
        -------
        msg : str
            the string representation of the card
        """
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
    +=========+======+=======+===+=======+===+=======+
    | AESURFS |  ID  | LABEL |   | LIST1 |   | LIST2 |
    +---------+------+-------+---+-------+---+-------+
    | AESURFS | 6001 | ELEV  |   | 6002  |   | 6003  |
    +---------+------+-------+---+-------+---+-------+
    """
    type = 'AESURFS'

    def __init__(self, aesid, label, list1, list2, comment=''):
        """
        Creates an AESURFS card

        Parameters
        ----------
        aesid : int
            the unique id
        label : str
            the AESURF name
        list1 / list2 : int / None
            the list (AELIST) of node ids for the primary/secondary
            control surface(s) on the AESURF card
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.aesid = aesid
        self.label = label
        self.list1 = list1
        self.list2 = list2

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AESURFS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        aesid = integer(card, 1, 'ID')
        label = string(card, 2, 'label')
        list1 = integer(card, 4, 'list1')
        list2 = integer(card, 6, 'list2')
        assert len(card) <= 7, 'len(AESURFS card) = %i\ncard=%s' % (len(card), card)
        return AESURFS(aesid, label, list1, list2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        aesid = data[0]
        label = data[1]
        list1 = data[2]
        list2 = data[3]
        assert len(data) == 4, 'data = %s' % data
        return AESURFS(aesid, label, list1, list2, comment=comment)

    #def uncross_reference(self):
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
            the fields that define the card
        """
        list_fields = ['AESURFS', self.aesid, self.label, None, self.list1, None,
                       self.list2]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)


class Aero(BaseCard):
    """Base class for AERO and AEROS cards."""
    def __init__(self):
        """
        Common class for AERO, AEROS

        Attributes
        ----------
        acsid : int; default=0
            aerodyanmic coordinate system
            defines the direction of the wind
        sym_xz : int; default=0
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        sym_xy : int; default=0
            xy symmetry flag (+1=symmetry; -1=antisymmetric)
        """
        self.sym_xy = None
        self.sym_xz = None
        self.acsid = None
        self.acsid_ref = None

    def Acsid(self):
        try:
            return self.acsid_ref.cid
        except AttributeError:
            return self.acsid

    def is_symmetric_xy(self):
        if self.sym_xy == 1:
            return True
        return False

    def is_symmetric_xz(self):
        if self.sym_xz == 1:
            return True
        return False

    def is_anti_symmetric_xy(self):
        if self.sym_xy == -1:
            return True
        return False

    def is_anti_symmetric_xz(self):
        if self.sym_xy == -1:
            return True
        return False

    def set_ground_effect(self, enable):
        if enable:
            self.sym_xy = -1
        else:
            self.sym_xy = 1


class AERO(Aero):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.

    +------+-------+----------+------+--------+-------+-------+
    | 1    | 2     | 3        | 4    | 5      | 6     | 7     |
    +======+=======+==========+======+========+=======+=======+
    | AERO | ACSID | VELOCITY | REFC | RHOREF | SYMXZ | SYMXY |
    +------+-------+----------+------+--------+-------+-------+
    | AERO |   3   |   1.3+   | 100. |  1.-5  |   1   |  -1   |
    +------+-------+----------+------+--------+-------+-------+
    """
    type = 'AERO'
    _field_map = {
        1: 'acsid', 2:'velocity', 3:'cRef', 4:'rhoRef', 5:'symXZ',
        6:'symXY',
    }

    def __init__(self, velocity, cref, rho_ref, acsid=0, sym_xz=0, sym_xy=0, comment=''):
        """
        Creates an AERO card

        Parameters
        ----------
        velocity : float
            the airspeed
        cref : float
            the aerodynamic chord
        rho_ref : float
            FLFACT density scaling factor
        acsid : int; default=0
            aerodyanmic coordinate system
            defines the direction of the wind
        sym_xz : int; default=0
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        sym_xy : int; default=0
            xy symmetry flag (+1=symmetry; -1=antisymmetric)
        comment : str; default=''
            a comment for the card
        """
        Aero.__init__(self)
        if comment:
            self.comment = comment

        #: Aerodynamic coordinate system identification
        if acsid is None:
            acsid = 0
        self.acsid = acsid

        #: Velocity for aerodynamic force data recovery and to calculate the BOV
        #: parameter
        self.velocity = velocity

        #: Reference length for reduced frequency
        self.cref = cref

        #: Reference density
        self.rho_ref = rho_ref

        #: Symmetry key for the aero coordinate x-z plane. See Remark 6.
        #: (Integer = +1 for symmetry, 0 for no symmetry, and -1 for antisymmetry;
        #: Default = 0)
        self.sym_xz = sym_xz

        #: The symmetry key for the aero coordinate x-y plane can be used to
        #: simulate ground effect. (Integer = -1 for symmetry, 0 for no symmetry,
        #: and +1 for antisymmetry; Default = 0)
        self.sym_xy = sym_xy

    def validate(self):
        if not isinstance(self.acsid, integer_types):
            msg = 'AERO acsid=%s expected int, got %s' % (
                self.acsid, type(self.acsid))
            raise TypeError(msg)
        if not isinstance(self.sym_xz, integer_types):
            msg = 'AERO acsid=%s sym_xz=%s; expected int, got %s' % (
                self.acsid, self.sym_xz, type(self.sym_xz))
            raise TypeError(msg)
        if not isinstance(self.sym_xy, integer_types):
            msg = 'AERO acsid=%s sym_xy=%s; expected int, got %s' % (
                self.acsid, self.sym_xy, type(self.sym_xy))
            raise TypeError(msg)

    def cross_reference(self, model):
        """
        Cross refernece aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.
        """
        msg = ' which is required by AERO'
        self.acsid_ref = model.Coord(self.acsid, msg=msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AERO card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        acsid = integer_or_blank(card, 1, 'acsid', 0)
        velocity = double_or_blank(card, 2, 'velocity')
        cref = double(card, 3, 'cRef')
        rho_ref = double(card, 4, 'rho_ref')
        sym_xz = integer_or_blank(card, 5, 'symXZ', 0)
        sym_xy = integer_or_blank(card, 6, 'symXY', 0)
        assert len(card) <= 7, 'len(AERO card) = %i\ncard=%s' % (len(card), card)
        return AERO(velocity, cref, rho_ref, acsid=acsid, sym_xz=sym_xz, sym_xy=sym_xy,
                    comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        acsid = data[0]
        velocity = data[1]
        cref = data[2]
        rho_ref = data[3]
        sym_xz = data[4]
        sym_xy = data[5]
        assert len(data) == 6, 'data = %s' % data
        return AERO(acsid, velocity, cref, rho_ref, sym_xz, sym_xy,
                    comment=comment)

        # T is the tabular function
        #angle = self.wg*self.t*(t-(x-self.x0)/self.V)

    def uncross_reference(self):
        self.acsid_ref = None

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : List[int/float/str]
           the fields that define the card
        """
        list_fields = ['AERO', self.Acsid(), self.velocity, self.cref,
                       self.rho_ref, self.sym_xz, self.sym_xy]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[varies]
          the fields that define the card
        """
        sym_xz = set_blank_if_default(self.sym_xz, 0)
        sym_xy = set_blank_if_default(self.sym_xy, 0)
        list_fields = ['AERO', self.Acsid(), self.velocity, self.cref,
                       self.rho_ref, sym_xz, sym_xy]
        return list_fields

    def write_card(self, size=8, is_double=False):
        """
        Writes the card with the specified width and precision

        Parameters
        ----------
        size : int (default=8)
            size of the field; {8, 16}
        is_double : bool (default=False)
            is this card double precision

        Returns
        -------
        msg : str
            the string representation of the card
        """
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class AEROS(Aero):
    """
    Gives basic aerodynamic parameters for unsteady aerodynamics.

    +-------+-------+-------+------+------+-------+-------+-------+
    | 1     | 2     | 3     | 4    | 5    | 6     |   7   |   8   |
    +=======+=======+=======+======+======+=======+=======+=======+
    | AEROS | ACSID | RCSID | REFC | REFB | REFS  | SYMXZ | SYMXY |
    +-------+-------+-------+------+------+-------+-------+-------+

    +-------+-------+-------+------+------+-------+-------+-------+
    | AEROS |   10  |   20  | 10.  | 100. | 1000. |   1   |       |
    +-------+-------+-------+------+------+-------+-------+-------+
    """
    type = 'AEROS'
    _field_map = {
        1: 'acsid', 2:'rcsid', 3:'cRef', 4:'bRef', 5:'Sref',
        6:'symXZ', 7:'symXY',
    }

    def __init__(self, cref, bref, sref, acsid=0, rcsid=0, sym_xz=0, sym_xy=0, comment=''):
        """
        Creates an AEROS card

        Parameters
        ----------
        cref : float
            the aerodynamic chord
        bref : float
            the wing span
            for a half model, this should be the full span
            for a full model, this should be the full span
        sref : float
            the wing area
            for a half model, this should be the half area
            for a full model, this should be the full area
        acsid : int; default=0
            aerodyanmic coordinate system
            defines the direction of the wind
        rcsid : int; default=0
            coordinate system for rigid body motions
        sym_xz : int; default=0
            xz symmetry flag (+1=symmetry; -1=antisymmetric)
        sym_xy : int; default=0
            xy symmetry flag (+1=symmetry; -1=antisymmetric)
        comment : str; default=''
            a comment for the card
        """
        Aero.__init__(self)
        if comment:
            self.comment = comment

        #: Aerodynamic coordinate system identification.
        self.acsid = acsid

        #: Reference coordinate system identification for rigid body motions.
        self.rcsid = rcsid

        #: Reference chord length
        self.cref = cref

        #: Reference span
        self.bref = bref

        #: Reference wing area
        self.sref = sref

        #: Symmetry key for the aero coordinate x-z plane. See Remark 6.
        #: (Integer = +1 for symmetry, 0 for no symmetry, and -1 for antisymmetry;
        #: Default = 0)
        self.sym_xz = sym_xz

        #: The symmetry key for the aero coordinate x-y plane can be used to
        #: simulate ground effects. (Integer = +1 for antisymmetry, 0 for no
        #: symmetry, and -1 for symmetry; Default = 0)
        self.sym_xy = sym_xy
        if self.acsid is None:
            self.acsid = 0
        if self.rcsid is None:
            self.rcsid = 0
        if self.sym_xz is None:
            self.sym_xz = 0
        if self.sym_xy is None:
            self.sym_xy = 0

    def Acsid(self):
        try:
            return self.acsid_ref.cid
        except AttributeError:
            return self.acsid

    def Rcsid(self):
        try:
            return self.rcsid_ref.cid
        except AttributeError:
            return self.rcsid

    def validate(self):
        msg = ''
        if not isinstance(self.acsid, integer_types):
            msg += 'acsid=%s must be an integer; type=%s\n' % (self.acsid, type(self.acsid))
        if not isinstance(self.rcsid, integer_types):
            msg += 'rcsid=%s must be an integer; type=%s\n' % (self.rcsid, type(self.rcsid))
        if not isinstance(self.cref, float):
            msg += 'cref=%s must be an float; type=%s\n' % (self.cref, type(self.cref))
        if not isinstance(self.bref, float):
            msg += 'bref=%s must be an float; type=%s\n' % (self.bref, type(self.bref))
        if not isinstance(self.sref, float):
            msg += 'sref=%s must be an float; type=%s\n' % (self.sref, type(self.sref))
        if not isinstance(self.sym_xz, integer_types):
            msg += 'sym_xz=%s must be an integer; type=%s\n' % (self.sym_xz, type(self.sym_xz))
        if not isinstance(self.sym_xy, integer_types):
            msg += 'sym_xy=%s must be an integer; type=%s\n' % (self.sym_xy, type(self.sym_xy))
        if msg:
            raise RuntimeError('There are errors on the AEROS card:\n%s%s' % (msg, self))

    def cross_reference(self, model):
        """
        Cross refernece aerodynamic coordinate system.

        Parameters
        ----------
        model : BDF
            The BDF object.
        """
        msg = ' which is required by AEROS'
        self.acsid_ref = model.Coord(self.acsid, msg=msg)
        self.rcsid_ref = model.Coord(self.rcsid, msg=msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an AEROS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        acsid = integer_or_blank(card, 1, 'acsid', 0)
        rcsid = integer_or_blank(card, 2, 'rcsid', 0)
        cref = double(card, 3, 'cRef')
        bref = double(card, 4, 'bRef')
        sref = double(card, 5, 'Sref')
        sym_xz = integer_or_blank(card, 6, 'sym_xz', 0)
        sym_xy = integer_or_blank(card, 7, 'sym_xy', 0)
        assert len(card) <= 8, 'len(AEROS card) = %i\ncard=%s' % (len(card), card)
        return AEROS(cref, bref, sref, acsid, rcsid, sym_xz, sym_xy,
                     comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        acsid = data[0]
        rcsid = data[1]
        cref = data[2]
        bref = data[3]
        sref = data[4]
        sym_xz = data[5]
        sym_xy = data[6]
        assert len(data) == 7, 'data = %s' % data
        return AEROS(cref, bref, sref, acsid, rcsid, sym_xz, sym_xy,
                     comment=comment)

    def uncross_reference(self):
        self.acsid_ref = None
        self.rcsid_ref = None

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
        """
        list_fields = ['AEROS', self.Acsid(), self.Rcsid(), self.cref,
                       self.bref, self.sref, self.sym_xz, self.sym_xy]
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : List[varies]
          the fields that define the card
        """
        sym_xz = set_blank_if_default(self.sym_xz, 0)
        sym_xy = set_blank_if_default(self.sym_xy, 0)
        list_fields = ['AEROS', self.Acsid(), self.Rcsid(), self.cref,
                       self.bref, self.sref, sym_xz, sym_xy]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CSSCHD(Aero):
    """
    Defines a scheduled control surface deflection as a function of
    Mach number and angle of attack.

    +--------+-----+-------+--------+-------+-------+
    |    1   |  2  |   3   |   4    |   5   |   6   |
    +========+=====+=======+========+=======+=======+
    | CSSCHD | SlD | AESID | LALPHA | LMACH | LSCHD |
    +--------+-----+-------+--------+-------+-------+

    +--------+-----+-------+--------+-------+-------+
    | CSSCHD |  5  |  50   |   12   |   15  |   25  |
    +--------+-----+-------+--------+-------+-------+
    """
    type = 'CSSCHD'
    _field_map = {
        1: 'sid', 2:'aesid', 3:'lalpha', 4:'lmach', 5:'lschd',
    }

    def __init__(self, sid, aesid, lschd, lalpha=None, lmach=None, comment=''):
        """
        Creates an CSSCHD card, which defines a specified control surface
        deflection as a function of Mach and alpha (used in SOL 144/146).

        Parameters
        ----------
        sid : int
            the unique id
        aesid : int
            the control surface (AESURF) id
        lalpha : int; default=None
            the angle of attack profile (AEFACT) id
        lmach : int; default=None
            the mach profile (AEFACT) id
        lschd : int; default=None
            the control surface deflection profile (AEFACT) id
        comment : str; default=''
            a comment for the card
        """
        Aero.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.aesid = aesid
        self.lalpha = lalpha
        self.lmach = lmach
        self.lschd = lschd

    def validate(self):
        if not(self.lalpha is None or isinstance(self.lalpha, integer_types)):
            raise TypeError('lalpha=%r must be an int or None' % self.lalpha)

        if not(self.lmach is None or isinstance(self.lmach, integer_types)):
            raise TypeError('lmach=%r must be an int or None' % self.lmach)

        if self.lalpha is None and self.lmach is None:
            msg = ('CSSCHD sid=%s; lalpha and lmach are both None'
                   ' (one must be an integer (AEFACT)\n%s' % (self.sid, str(self)))
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CSSCHD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        aesid = integer(card, 2, 'aesid')             # AESURF
        lalpha = integer_or_blank(card, 3, 'lAlpha')  # AEFACT
        lmach = integer_or_blank(card, 4, 'lMach')    # AEFACT
        lschd = integer(card, 5, 'lSchd')             # AEFACT
        assert len(card) <= 6, 'len(CSSCHD card) = %i\ncard=%s' % (len(card), card)
        return CSSCHD(sid, aesid, lalpha, lmach, lschd, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        aesid = data[1]   # AESURF
        lalpha = data[2]  # AEFACT
        lmach = data[3]   # AEFACT
        lschd = data[4]   # AEFACT
        return CSSCHD(sid, aesid, lalpha, lmach, lschd, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CSSCHD sid=%s' % self.sid
        self.aesid = model.AESurf(self.aesid, msg=msg)
        self.aesid_ref = self.aesid
        self.lalpha = model.AEFact(self.lalpha, msg=msg)
        self.lalpha_ref = self.lalpha
        self.lmach = model.AEFact(self.lmach, msg=msg)
        self.lmach_ref = self.lmach
        self.lschd = model.AEFact(self.lschd, msg=msg)
        self.lschd_ref = self.lschd

    def safe_cross_reference(self, model):
        msg = ' which is required by CSSCHD sid=%s' % self.sid
        try:
            self.aesid = model.AESurf(self.aesid, msg=msg)
            self.aesid_ref = self.aesid
        except KeyError:
            pass

        try:
            self.lalpha = model.AEFact(self.lalpha, msg=msg)
            self.lalpha_ref = self.lalpha
        except KeyError:
            pass

        try:
            self.lmach = model.AEFact(self.lmach, msg=msg)
            self.lmach_ref = self.lmach
        except KeyError:
            pass

        try:
            self.lschd = model.AEFact(self.lschd, msg=msg)
            self.lschd_ref = self.lschd
        except KeyError:
            pass

    def uncross_reference(self):
        self.aesid = self.AESid()
        self.lalpha = self.LAlpha()
        self.lmach = self.LMach()
        self.lschd = self.LSchd()
        del self.aesid_ref
        del self.lalpha_ref
        del self.lmach_ref
        del self.lschd_ref

    def AESid(self):
        if isinstance(self.aesid, integer_types):
            return self.aesid
        return self.aesid_ref.aesid

    def LAlpha(self):
        if self.lalpha is None or isinstance(self.lalpha, integer_types):
            return self.lalpha
        return self.lalpha_ref.sid

    def LMach(self):
        if self.lmach is None or isinstance(self.lmach, integer_types):
            return self.lmach
        return self.lmach_ref.sid

    def LSchd(self):
        if self.lschd is None or isinstance(self.lschd, integer_types):
            return self.lschd
        return self.lschd_ref.sid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
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
    +========+=====+=====+====+=======+========+========+========+======+
    | CAERO1 | EID | PID | CP | NSPAN | NCHORD |  LSPAN | LCHORD | IGID |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |        |  X1 | Y1  | Z1 |  X12  |   X4   |   Y4   |   Z4   | X43  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+

        ::

          1
          | \
          |   \
          |     \
          |      4
          |      |
          |      |
          2------3

    Attributes
    ----------
    eid : int
        element id
    pid : int, PAERO1
        int : PAERO1 ID
        PAERO1 : PAERO1 object (xref)
    igid : int
        Group number
    p1 : (1, 3) ndarray float
        xyz location of point 1 (leading edge; inboard)
    p4 : (1, 3) ndarray float
        xyz location of point 4 (leading edge; outboard)
    x12 : float
        distance along the flow direction from node 1 to node 2; (typically x, root chord)
    x43 : float
        distance along the flow direction from node 4 to node 3; (typically x, tip chord)

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
    comment : str; default=''
         a comment for the card
    """
    type = 'CAERO1'
    _field_map = {
        1: 'sid', 2:'pid', 3:'cp', 4:'nspan', 5:'nchord',
        6:'lspan', 7:'lchord', 8:'igid', 12:'x12', 16:'x43',
    }
    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the CAERO1 card

        Parameters
        ----------
        n : int
            the field number to update
        value : int/float
            the value for the appropriate field
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
            raise KeyError('Field %r is an invalid CAERO1 entry.' % n)

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the CAERO1 card

        Parameters
        ----------
        n : int
            the field number to update
        value : int/float
            the value for the appropriate field
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
            raise KeyError('Field %r=%r is an invalid CAERO1 entry.' % (n, value))

    def __init__(self, eid, pid, igid, p1, x12, p4, x43,
                 cp=0, nspan=0, lspan=0, nchord=0, lchord=0, comment=''):
        """
        Defines a CAERO1 card, which defines a simplified lifting surface
        (e.g., wing/tail).

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO1
            int : PAERO1 ID
            PAERO1 : PAERO1 object (xref)
        igid : int
            Group number
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2; (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3; (typically x, tip chord)
        cp : int, CORDx; default=0
            int : coordinate system
            CORDx : Coordinate object (xref)
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        nchord : int; default=0
            int > 0 : N chordwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
            AEFACT : AEFACT object  (xref)
        lchord : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nchord
            int = 0 : use nchord
            AEFACT : AEFACT object  (xref)
        comment : str; default=''
             a comment for the card
        """
        if cp is None:
            cp = 0
        if lspan is None:
            lspan = 0
        if nspan is None:
            nspan = 0
        if nchord is None:
            nchord = 0
        if lchord is None:
            lchord = 0

        if comment:
            self.comment = comment
        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO2 entry.
        self.pid = pid

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.nchord = nchord
        self.lchord = lchord
        self.igid = igid
        self.p1 = p1
        self.x12 = x12
        self.p4 = p4
        self.x43 = x43

    def validate(self):
        msg = ''
        is_failed = False
        if self.nspan == 0 and self.lspan == 0:
            msg += 'NSPAN or LSPAN must be greater than 0; nspan=%r nlspan=%s\n' % (
                self.nspan, self.lspan)
            is_failed = True
        if self.nspan != 0 and self.lspan != 0:
            msg += 'Either NSPAN or LSPAN must 0; nspan=%r nlspan=%s\n' % (
                self.nspan, self.lspan)
            is_failed = True

        if self.nchord == 0 and self.lchord == 0:
            msg += 'NCHORD or LCHORD must be greater than 0; nchord=%r lchord=%s\n' % (
                self.nchord, self.lchord)
            is_failed = True
        if self.nchord != 0 and self.lchord != 0:
            msg += 'Either NCHORD or LCHORD must 0; nchord=%r lchord=%s\n' % (
                self.nchord, self.lchord)
            is_failed = True
        if is_failed:
            msg += str(self)
            raise ValueError(msg)
        assert len(self.p1) == 3, 'p1=%s' % self.p1
        assert len(self.p4) == 3, 'p4=%s' % self.p4

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', 0)
        nspan = integer_or_blank(card, 4, 'nspan', 0)
        nchord = integer_or_blank(card, 5, 'nchord', 0)
        lspan = integer_or_blank(card, 6, 'lspan', 0)
        lchord = integer_or_blank(card, 7, 'lchord', 0)
        igid = integer(card, 8, 'igid')

        p1 = np.array([
            double_or_blank(card, 9, 'x1', 0.0),
            double_or_blank(card, 10, 'y1', 0.0),
            double_or_blank(card, 11, 'z1', 0.0)])
        x12 = double_or_blank(card, 12, 'x12', 0.)

        p4 = np.array([
            double_or_blank(card, 13, 'x4', 0.0),
            double_or_blank(card, 14, 'y4', 0.0),
            double_or_blank(card, 15, 'z4', 0.0)])
        x43 = double_or_blank(card, 16, 'x43', 0.)

        assert len(card) <= 17, 'len(CAERO1 card) = %i\ncard=%s' % (len(card), card)
        return CAERO1(eid, pid, igid, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                      comment=comment)

    @classmethod
    def add_quad(cls, eid, pid, span, chord, igid,
                 p1, p2, p3, p4, cp=0, spanwise='y', comment=''):
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

        TODO: CP not handled correctly
        """
        x12 = p2[0] - p1[0]
        x43 = p3[0] - p4[0]
        nspan = 0
        lspan = 0
        nchord = 0
        lchord = 0

        if spanwise.lower() == 'y':
            y41 = p4[1] - p1[1]
            y32 = p3[1] - p2[1]
            dspan = max(y41, y32)
        elif spanwise.lower() == 'z':
            y41 = p4[2] - p1[2]
            y32 = p3[2] - p2[2]
            dspan = max(y41, y32)
        else:
            raise NotImplementedError('spanwise=%r; expected=[y, z]' % spanwise.lower())

        dx = max(x12, x43)
        if isinstance(span, integer_types):
            nspan = span
        elif isinstance(span, AEFACT):
            lspan = span.sid
        elif isinstance(span, float):
            nspan = int(math.ceil(dspan / span))
            if nspan <= 0:
                msg = 'y41=%s y32=%s; dspan=%s span=%s nspan=%s; nspan must be greater than 0' % (
                    y41, y32, dspan, span, nspan)
                raise ValueError(msg)
        else:
            raise TypeError(span)

        if isinstance(chord, integer_types):
            nchord = chord
        elif isinstance(chord, AEFACT):
            lchord = chord.sid
        elif isinstance(chord, float):
            nchord = int(math.ceil(dx / chord))
            if nchord <= 0:
                msg = 'x12=%s x43=%s; dx=%s chord=%s nchord=%s; nchord must be greater than 0' % (
                    x12, x43, dx, chord, nchord)
                raise ValueError(msg)
        else:
            raise TypeError(chord)

        return CAERO1(eid, pid, igid, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                      comment=comment)

    def _init_ids(self, dtype='int32'):
        """
        Fill `self.box_ids` with the sub-box ids. Shape is (nchord, nspan)
        """
        nchord, nspan = self.shape
        assert nchord >= 1, 'nchord=%s' % nchord
        assert nspan >= 1, 'nspan=%s' % nspan
        self.box_ids = np.zeros((nchord, nspan), dtype=dtype)

        try:
            for ichord in range(nchord):
                for ispan in range(nspan):
                    self.box_ids[ichord, ispan] = self.eid + ichord + ispan * nchord
        except OverflowError:
            if dtype == 'int64':
                msg = 'eid=%s ichord=%s ispan=%s nchord=%s' % (
                    self.eid, ichord, ispan, nchord)
                raise OverflowError(msg)
            self._init_ids(dtype='int64')

    def Cp(self):
        if isinstance(self.cp, integer_types):
            return self.cp
        return self.cp_ref.cid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid_ref.pid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CAERO1 eid=%s' % self.eid
        self.pid = model.PAero(self.pid, msg=msg)
        self.pid_ref = self.pid
        self.cp = model.Coord(self.cp, msg=msg)
        self.cp_ref = self.cp

        self.ascid_ref = model.Acsid(msg=msg)

        if self.nchord == 0:
            assert isinstance(self.lchord, integer_types), self.lchord
            self.lchord = model.AEFact(self.lchord, msg)
            self.lchord_ref = self.lchord
        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan = model.AEFact(self.lspan, msg)
            self.lspan_ref = self.lspan
        self._init_ids()

    def safe_cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CAERO1 eid=%s' % self.eid
        try:
            self.pid = model.PAero(self.pid, msg=msg)
            self.pid_ref = self.pid
        except KeyError:
            pass

        try:
            self.cp = model.Coord(self.cp, msg=msg)
            self.cp_ref = self.cp
        except KeyError:
            pass

        try:
            self.ascid_ref = model.Acsid(msg=msg)
        except KeyError:
            pass

        if self.nchord == 0:
            assert isinstance(self.lchord, integer_types), self.lchord
            try:
                self.lchord = model.AEFact(self.lchord, msg)
                self.lchord_ref = self.lchord
            except KeyError:
                pass

        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            try:
                self.lspan = model.AEFact(self.lspan, msg)
                self.lspan_ref = self.lspan
            except KeyError:
                pass

        self._init_ids()

    #def uncross_reference(self):
        #self.pid = self.Pid()
        #self.cp = self.Cp()
        #del self.pid_ref
        #del self.pid_ref, self.cp_ref

    def uncross_reference(self):
        self.pid = self.Pid()
        self.cp = self.Cp()
        self.lchord = self.get_LChord()
        self.lspan = self.get_LSpan()

    @property
    def min_max_eid(self):
        """
        Gets the min and max element ids of the CAERO card

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

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame
        """
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p4 = self.cp_ref.transform_node_to_global(self.p4)
        p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))
        p3 = p4 + self.ascid_ref.transform_vector_to_global(np.array([self.x43, 0., 0.]))
        return [p1, p2, p3, p4]

    @property
    def shape(self):
        """returns (nelements_nchord, nelements_span)"""
        if self.nchord == 0:
            x = self.lchord_ref.Di
            nchord = len(x) - 1
        else:
            nchord = self.nchord

        if self.nspan == 0:
            y = self.lspan_ref.Di
            nspan = len(y) - 1
        else:
            nspan = self.nspan
        if nchord < 1 or nspan < 1:
            msg = 'CAERO1 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                self.eid, self.nchord, self.nspan, self.lchord, self.lspan)
            raise RuntimeError(msg)
        return nchord, nspan

    def get_npanel_points_elements(self):
        """
        Gets the number of sub-points and sub-elements for the CAERO card

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
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel
        """
        if self.nchord == 0:
            x = self.lchord_ref.Di
            nchord = len(x) - 1
        else:
            nchord = self.nchord
            x = np.linspace(0., 1., nchord + 1)

        if self.nspan == 0:
            y = self.lspan_ref.Di
            nspan = len(y) - 1
        else:
            nspan = self.nspan
            y = np.linspace(0., 1., nspan + 1)

        if nchord < 1 or nspan < 1:
            msg = 'CAERO1 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                self.eid, self.nchord, self.nspan, self.lchord, self.lspan)
            raise RuntimeError(msg)
        return x, y

    def panel_points_elements(self):
        """
        Gets the sub-points and sub-elements for the CAERO card

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
        p2 = points[1]
        p3 = points[2]
        self.p4 = points[3]
        self.x12 = p2[0] - self.p1[0]
        self.x43 = p3[0] - self.p4[0]
        assert self.x12 >= 0., 'p1=%s p2=%s' % (self.p1, p2)
        assert self.x43 >= 0., 'p4=%s p3=%s' % (self.p4, p3)
        assert self.x12 > 0. or self.x43 > 0., 'points=%s' % (points)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

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
        return self.lchord_ref.sid

    def get_LSpan(self):
        if isinstance(self.lspan, integer_types):
            return self.lspan
        return self.lspan_ref.sid

    def repr_fields(self):
        """
        Gets the fields in their simplified form

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

    +--------+-----+-----+----+-----+------+-----+------+------+
    |    1   |  2  |  3  |  4 |  5  |   6  |   7 |   8  |  9   |
    +========+=====+=====+====+=====+======+=====+======+======+
    | CAERO2 | EID | PID | CP | NSB | NINT | LSB | LINT | IGID |
    +--------+-----+-----+----+-----+------+-----+------+------+
    |        | X1  |  Y1 | Z1 | X12 |      |     |      |      |
    +--------+-----+-----+----+-----+------+-----+------+------+
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
            raise KeyError('Field %r is an invalid CAERO2 entry.' % n)

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the CAERO2 card

        Parameters
        ----------
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
            raise KeyError('Field %r=%r is an invalid CAERO2 entry.' % (n, value))

    def __init__(self, eid, pid, igid, p1, x12,
                 cp=0, nsb=0, nint=0, lsb=0, lint=0, comment=''):
        """
        Defines a CAERO2 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO2
            int : PAERO2 ID
            PAERO2 : PAERO2 object (xref)
        igid : int
            Group number
        p1 : (1, 3) ndarray float
            xyz location of point 1 (forward position)
        x12 : float
            length of the CAERO2
        cp : int, CORDx; default=0
            int : coordinate system
            CORDx : Coordinate object (xref)
        nsb : int; default=0
            AEFACT id for defining the location of the slender body elements
        lsb : int; default=0
            AEFACT id for defining the location of interference elements
        nint : int; default=0
            Number of slender body elements
        lint : int; default=0
            Number of interference elements
        comment : str; default=''
            a comment for the card
        """
        if lsb is None:
            lsb = 0
        if lint is None:
            lint = 0
        if nint is None:
            nint = 0
        if nsb is None:
            nsb = 0

        if comment:
            self.comment = comment
        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO2 entry.
        self.pid = pid

        #: Coordinate system for locating point 1.
        self.cp = cp

        #: Number of slender body elements. If NSB > 0, then NSB equal
        #: divisions are assumed; if zero or blank, specify a list of
        #: divisions in LSB. (Integer >= 0)
        self.nsb = nsb

        #: Number of interference elements. If NINT > 0, then NINT equal
        #: divisions are assumed; if zero or blank, specify a list of
        #: divisions in LINT. (Integer >= 0)
        self.nint = nint

        #: ID of an AEFACT Bulk Data entry for slender body division
        #: points; used only if NSB is zero or blank. (Integer >= 0)
        self.lsb = lsb

        #: ID of an AEFACT data entry containing a list of division
        #: points for interference elements; used only if NINT is zero
        #: or blank. (Integer > 0)
        self.lint = lint

        #: Interference group identification. Aerodynamic elements with
        #: different IGIDs are uncoupled. (Integer >= 0)
        self.igid = igid

        #: Location of point 1 in coordinate system CP
        self.p1 = p1

        #: Length of body in the x-direction of the aerodynamic coordinate
        #: system.  (Real > 0)
        self.x12 = x12

    def validate(self):
        #print('nsb=%s lsb=%s' % (self.nsb, self.lsb))
        #print('nint=%s lint=%s' % (self.nint, self.lint))
        assert isinstance(self.lsb, integer_types), self.lsb
        assert isinstance(self.lint, integer_types), self.lint
        assert len(self.p1) == 3, 'p1=%s' % self.p1
        if self.nsb == 0 and self.lsb == 0:
            msg = 'nsb=%s lsb=%s; nsb or lsb must be > 0' % (self.nsb, self.lsb)
            raise ValueError(msg)
        if self.nint == 0 and self.lint == 0:
            msg = 'nint=%s lint=%s; nint or lint must be > 0' % (self.nint, self.lint)
            raise ValueError(msg)
        assert len(self.p1) == 3, 'p1=%s' % self.p1

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', 0)
        nsb = integer_or_blank(card, 4, 'nsb', 0)
        nint = integer_or_blank(card, 5, 'nint', 0)

        lsb = integer_or_blank(card, 6, 'nsb=%s lsb' % nsb, 0)
        lint = integer_or_blank(card, 7, 'nint=%s lint' % nint, 0)
        igid = integer(card, 8, 'igid')

        p1 = np.array([
            double_or_blank(card, 9, 'x1', 0.0),
            double_or_blank(card, 10, 'y1', 0.0),
            double_or_blank(card, 11, 'z1', 0.0)])
        x12 = double_or_blank(card, 12, 'x12', 0.)
        assert len(card) <= 13, 'len(CAERO2 card) = %i\ncard=%s' % (len(card), card)
        return CAERO2(eid, pid, igid, p1, x12,
                      cp=cp, nsb=nsb, nint=nint, lsb=lsb, lint=lint,
                      comment=comment)

    def Cp(self):
        if isinstance(self.cp, integer_types):
            return self.cp
        return self.cp_ref.cid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid_ref.pid

    def Lsb(self):  # AEFACT
        if self.lsb is None or isinstance(self.lsb, integer_types):
            return self.lsb
        return self.lsb_ref.sid

    def Lint(self):  # AEFACT
        if self.lint is None or isinstance(self.lint, integer_types):
            return self.lint
        return self.lint_ref.sid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CAERO2 eid=%s' % self.eid
        self.pid = model.PAero(self.pid, msg=msg)  # links to PAERO2
        self.pid_ref = self.pid
        self.cp = model.Coord(self.cp, msg=msg)
        self.cp_ref = self.cp
        if self.nsb == 0:
            self.lsb = model.AEFact(self.lsb, msg=msg)
            self.lsb_ref = self.lsb
        if self.nint == 0:
            self.lint = model.AEFact(self.lint, msg=msg)
            self.lint_ref = self.lint
        self.ascid_ref = model.Acsid(msg=msg)

    def safe_cross_reference(self, model, debug=False):
        msg = ' which is required by CAERO2 eid=%s' % self.eid
        try:
            self.pid = model.PAero(self.pid, msg=msg)  # links to PAERO2
            self.pid_ref = self.pid
        except KeyError:
            pass

        try:
            self.cp = model.Coord(self.cp, msg=msg)
            self.cp_ref = self.cp
        except KeyError:
            pass

        if self.nsb == 0:
            try:
                self.lsb = model.AEFact(self.lsb, msg=msg)
                self.lsb_ref = self.lsb
            except KeyError:
                pass
        if self.nint == 0:
            try:
                self.lint = model.AEFact(self.lint, msg=msg)
                self.lint_ref = self.lint
            except KeyError:
                pass
        try:
            self.ascid_ref = model.Acsid(msg=msg)
        except KeyError:
            pass

    def uncross_reference(self):
        self.pid = self.Pid()
        self.cp = self.Cp()
        if self.nsb == 0:
            self.lsb = self.Lsb()
            del self.lsb_ref
        if self.nint == 0:
            self.lint = self.Lint()
            del self.lint_ref
        del self.pid_ref, self.cp_ref

    def get_points(self):
        """
        creates a 1D representation of the CAERO2
        """
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))

        #print("x12 = %s" % self.x12)
        #print("pcaero[%s] = %s" % (self.eid, [p1,p2]))
        return [p1, p2]

    def get_points_elements_3d(self):
        """
        Gets the points/elements in 3d space as CQUAD4s
        The idea is that this is used by the GUI to display CAERO panels.

        TODO: doesn't support the aero coordinate system
        """
        paero2 = self.pid_ref

        if self.nsb == 0:
            xstation = self.lsb_ref.data
            nx = len(xstation) - 1
            #print('xstation = ', xstation)
        else:
            nx = self.nsb
            station = np.linspace(0., nx, num=nx+1) # *dx?
        assert nx > 0, 'nx=%s' % nx


        #print('paero2 - pid=%s lrsb=%s lrib=%s' % (paero2.pid, paero2.lrsb, paero2.lrib))
        if paero2.lrsb in [0, None]:
            radii_slender = np.ones(nx + 1) * paero2.width
        else:
            radii_slender = paero2.lrsb_ref.data

        # TODO: not suppported
        if paero2.lrib in [0, None]:
            radii_interference = np.ones(nx + 1) * paero2.width
        else:
            #print('lrib = ', paero2.lrib)
            radii_interference = paero2.lrib_ref.data
        radii = radii_slender

        # TODO: not suppported
        #theta_interference1 = paero2.theta1
        #theta_interference2 = paero2.theta2

        if self.nsb != 0:
            p1, p2 = self.get_points()
            L = p2 - p1
            #print('L=%s nx=%s' % (L, nx))
            dxyz = L / nx
            #print('dxyz\n%s' % (dxyz))
            dx, dy, dz = dxyz
            xstation = station * dx
            ystation = station * dy
            zstation = station * dz
        else:
            p1, p2 = self.get_points()
            L = p2 - p1
            dxi = xstation.max() - xstation.min()

            #print('L=%s nx=%s dxi=%s' % (L, nx, dxi))
            xratio = xstation / dxi
            #print('xstation/dxi=%s' % xratio)
            dxyz = np.zeros((nx+1, 3))
            for i, xr in enumerate(xratio):
                dxyz[i, :] = xr * L
            ystation = dxyz[:, 1]
            zstation = dxyz[:, 2]

        # I think this just lets you know what directions it can pivot in
        # and therefore doesn't affect visualization
        #assert paero2.orient == 'ZY', paero2.orient
        aspect_ratio = paero2.AR

        #Rs = []
        assert len(radii) == (nx + 1), 'len(radii)=%s nx=%s' % (len(radii), nx)
        if len(xstation) != (nx + 1):
            msg = 'len(xstation)=%s nx=%s\nxstation=%s\n%s' % (
                len(xstation), nx, xstation, str(self))
            raise RuntimeError(msg)

        xs = []
        ys = []
        zs = []
        yzs = []
        for i, xi, yi, zi, radius in zip(count(), xstation, ystation, zstation, radii):
            #print('  station=%s xi=%.4f radius=%s' % (i, xi, radius))
            yz = self.create_ellipse(aspect_ratio, radius)
            yzs.append(yz)
            try:
                y = yz[:, 0] + yi
                z = yz[:, 1] + zi
            except ValueError:
                print('yz = %s' % yz)
                print('yz.shape = %s' % str(yz.shape))
                print('dy = %s' % dy)
                print('dz = %s' % dz)
                raise
            ntheta = yz.shape[0]
            x = np.ones(ntheta) * xi
            xs.append(x)
            ys.append(y)
            zs.append(z)
            #Rs.append(np.sqrt(y**2 + z**2))
        #print('yz.shape=%s xs.shape=%s' % (str(np.array(yzs).shape), str(np.array(xs).shape)))
        #xyz = np.hstack([yzs, xs])
        xs = np.array(xs)
        ys = np.array(ys)
        zs = np.array(zs)
        try:
            xyz = np.vstack([
                np.hstack(xs),
                np.hstack(ys),
                np.hstack(zs),
            ]).T + p1
        except:
            print('xs =', xs.shape)
            print('ys =', ys.shape)
            print('zs =', zs.shape)
            raise

        #R = np.hstack(Rs)
        #print('xyz.shape =', xyz.shape)
        #print('xyz =', xyz)
        #print('R =', R)

        ny = ntheta
        elems = elements_from_quad(nx+1, ny)
        #print('elems =\n', elems)
        return xyz, elems

    @staticmethod
    def create_ellipse(aspect_ratio, radius, thetas=None):
        r"""
        a : major radius
        b : minor radius

        Parameters
        ----------
        aspect_ratio : float
            AR = height/width

        https://en.wikipedia.org/wiki/Ellipse#Polar_form_relative_to_center

        .. math::

            r(\theta )={\frac {ab}{\sqrt {(b\cos \theta )^{2}+(a\sin \theta )^{2}}}}

        R(theta) = a*b / ((b*cos(theta))**2 + (a*sin(theta))**2)

        TODO: doesn't support the aero coordinate system
        """
        if thetas is None: # 41
            thetas = np.radians(np.linspace(0., 360., 17)) # 4,8,12,16,... becomes 5,9,13,17,...
        ntheta = len(thetas)

        a = radius
        b = radius * aspect_ratio
        if a == 0.0 and b == 0.0:
            xy = np.zeros((ntheta, 2)) # this is just R
            return xy

        R = a * b / np.sqrt((b*np.cos(thetas))**2 + (a*np.sin(thetas))**2)
        x = R * np.cos(thetas)
        y = R * np.sin(thetas)

        xy = np.vstack([x, y]).T
        assert xy.shape == (ntheta, 2), xy.shape
        return xy

    def set_points(self, points):
        self.p1 = points[0]
        self.p2 = points[1]
        x12 = self.p2 - self.p1
        self.x12 = x12[0]

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card
        """
        list_fields = (['CAERO2', self.eid, self.Pid(), self.Cp(), self.nsb,
                        self.nint, self.Lsb(), self.Lint(), self.igid, ] + list(self.p1)
                       + [self.x12])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card
        """
        cp = set_blank_if_default(self.Cp(), 0)
        nint = set_blank_if_default(self.nint, 0)
        lsb = set_blank_if_default(self.Lsb(), 0)
        lint = set_blank_if_default(self.Lint(), 0)
        list_fields = (['CAERO2', self.eid, self.Pid(), cp, self.nsb, nint,
                        lsb, lint, self.igid, ] + list(self.p1) +
                       [self.x12])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO3(BaseCard):
    type = 'CAERO3'
    def __init__(self, eid, pid, list_w,
                 p1, x12, p4, x43,
                 cp=0, list_c1=None, list_c2=None,
                 comment=''):
        """
        eid : int
            element id
        pid : int
            PAERO3 property id
        cp : int; default=0
            coordinate system for locating point 1
        list_w : int
            ???
        list_c1 : int; default=None
            ???
        list_c2 : int; default=None
            ???
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment

        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO3 entry.
        self.pid = pid

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.list_w = list_w
        self.list_c1 = list_c1
        self.list_c2 = list_c2
        self.p1 = p1
        self.x12 = x12
        self.p4 = p4
        self.x43 = x43

    def validate(self):
        assert len(self.p1) == 3, 'p1=%s' % self.p1
        assert len(self.p4) == 3, 'p4=%s' % self.p4
        assert self.x12 > 0., 'x12=%s' % self.x12
        assert self.x43 >= 0., 'x43=%s' % self.x43

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', 0)
        list_w = integer(card, 4, 'list_w')
        list_c1 = integer_or_blank(card, 5, 'list_c1')
        list_c2 = integer_or_blank(card, 6, 'list_c2')
        p1 = np.array([
            double_or_blank(card, 9, 'x1', 0.0),
            double_or_blank(card, 10, 'y1', 0.0),
            double_or_blank(card, 11, 'z1', 0.0)])
        x12 = double(card, 12, 'x12')
        p4 = np.array([
            double_or_blank(card, 13, 'x4', 0.0),
            double_or_blank(card, 14, 'y4', 0.0),
            double_or_blank(card, 15, 'z4', 0.0)])
        x43 = double_or_blank(card, 16, 'x43', 0.0)
        assert len(card) <= 17, 'len(CAERO3 card) = %i\ncard=%s' % (len(card), card)
        return CAERO3(eid, pid, list_w, p1, x12, p4, x43,
                      cp=cp, list_c1=list_c1, list_c2=list_c2, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CAERO3 eid=%s' % self.eid
        self.pid = model.PAero(self.pid, msg=msg)  # links to PAERO3
        self.pid_ref = self.pid
        self.cp = model.Coord(self.cp, msg=msg)
        self.cp_ref = self.cp
        if self.list_w is not None:
            self.list_w = model.AEFact(self.list_w, msg=msg)
        if self.list_c1 is not None:
            self.list_c1 = model.AEFact(self.list_c1, msg=msg)
        if self.list_c2 is not None:
            self.list_c2 = model.AEFact(self.list_c2, msg=msg)
        self.ascid_ref = model.Acsid(msg=msg)

    def safe_cross_reference(self, model):
        msg = ', which is required by CAERO3 eid=%s' % self.eid
        try:
            self.pid = model.PAero(self.pid, msg=msg)  # links to PAERO3
            self.pid_ref = self.pid
        except KeyError:
            model.log.warning('cannot find PAERO3 pid=%s%s' % (self.pid, msg))

        try:
            self.cp = model.Coord(self.cp, msg=msg)
            self.cp_ref = self.cp
        except KeyError:
            model.log.warning('cannot find PAERO3 pid=%s%s' % (self.pid, msg))

        if self.list_w is not None:
            try:
                self.list_w = model.AEFact(self.list_w, msg=msg)
            except KeyError:
                model.log.warning('cannot find an AEFACT pid=%s%s' % (self.list_w, msg))

        if self.list_c1 is not None:
            try:
                self.list_c1 = model.AEFact(self.list_c1, msg=msg)
            except KeyError:
                model.log.warning('cannot find an AEFACT pid=%s%s' % (self.list_c1, msg))

        if self.list_c2 is not None:
            try:
                self.list_c2 = model.AEFact(self.list_c2, msg=msg)
            except KeyError:
                model.log.warning('cannot find an AEFACT list_c2=%s%s' % (self.list_c2, msg))
        try:
            self.ascid_ref = model.Acsid(msg=msg)
        except KeyError:
            model.log.warning('cannot find an aero coordinate system for %s' % msg)

    def uncross_reference(self):
        self.pid = self.Pid()
        self.cp = self.Cp()
        if self.list_w != self.List_w():
            self.list_w = self.List_w()
        if self.list_c1 != self.List_c1():
            self.list_c1 = self.List_c1()
        if self.list_c2 != self.List_c2():
            self.list_c2 = self.List_c2()
        del self.pid_ref, self.cp_ref, self.ascid_ref

    def get_points(self):
        """
        Get the 4 corner points for the CAERO card

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame
        """
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p4 = self.cp_ref.transform_node_to_global(self.p4)
        p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))
        p3 = p4 + self.ascid_ref.transform_vector_to_global(np.array([self.x43, 0., 0.]))
        return [p1, p2, p3, p4]

    def panel_points_elements(self):
        """
        Gets the sub-points and sub-elements for the CAERO card

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

    def get_npanel_points_elements(self):
        """
        Gets the number of sub-points and sub-elements for the CAERO card

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
    def shape(self):
        """returns (nelements_nchord, nelements_span)"""
        nchord = 2
        nspan = self.pid_ref.nbox
        return nchord, nspan

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel
        """
        nchord, nspan = self.shape
        x = np.linspace(0., 1., nchord + 1)
        y = np.linspace(0., 1., nspan + 1)

        if nchord < 1 or nspan < 1:
            msg = 'CAERO3 eid=%s nchord=%s nspan=%s' % (
                self.eid, nchord, nspan)
            raise RuntimeError(msg)
        return x, y

    #def get_points_elements_3d(self):
        #"""
        #Gets the points/elements in 3d space as CQUAD4s
        #The idea is that this is used by the GUI to display CAERO panels.

        #TODO: doesn't support the aero coordinate system
        #"""
        #paero2 = self.pid_ref

    def Cp(self):
        if isinstance(self.cp, integer_types):
            return self.cp
        return self.cp_ref.cid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid_ref.pid

    def List_w(self):
        if self.list_w is None or isinstance(self.list_w, integer_types):
            return self.list_w
        return self.list_w.sid

    def List_c1(self):
        if self.list_c1 is None or isinstance(self.list_c1, integer_types):
            return self.list_c1
        return self.list_c1.sid

    def List_c2(self):
        if self.list_c2 is None or isinstance(self.list_c2, integer_types):
            return self.list_c2
        return self.list_c2.sid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card
        """
        list_fields = (['CAERO3', self.eid, self.Pid(), self.Cp(), self.List_w(),
                        self.List_c1(), self.List_c2(), None, None] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card
        """
        cp = set_blank_if_default(self.Cp(), 0)
        list_fields = (['CAERO3', self.eid, self.Pid(), cp, self.List_w(),
                        self.List_c1(), self.List_c2(), None, None] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO4(BaseCard):
    """
    Aerodynamic Macro-Strip Element Connection
    Defines an aerodynamic macro element for Strip theory.

    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
    +========+=====+=====+====+=======+========+========+========+======+
    | CAERO4 | EID | PID | CP | NSPAN | NCHORD |        |        |      |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |        |  X1 | Y1  | Z1 |  X12  |   X4   |   Y4   |   Z4   | X43  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    """
    type = 'CAERO4'
    def __init__(self, eid, pid, p1, x12, p4, x43,
                 cp=0, nspan=0, lspan=0, comment=''):
        """
        Defines a CAERO4 card, which defines a strip theory surface.

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO4
            int : PAERO4 ID
            PAERO4 : PAERO4 object (xref)
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2; (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3; (typically x, tip chord)
        cp : int, CORDx; default=0
            int : coordinate system
            CORDx : Coordinate object (xref)
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
            AEFACT : AEFACT object  (xref)
        comment : str; default=''
             a comment for the card
        """
        if comment:
            self.comment = comment

        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO4 entry.
        self.pid = pid

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.p1 = p1
        self.x12 = x12
        self.p4 = p4
        self.x43 = x43

    def validate(self):
        if self.nspan == 0 and self.lspan == 0:
            msg = 'NSPAN or LSPAN must be greater than 0; nspan=%r nlspan=%s\n' % (
                self.nspan, self.lspan)
            raise RuntimeError(msg)
        assert len(self.p1) == 3, 'p1=%s' % self.p1
        assert len(self.p4) == 3, 'p4=%s' % self.p4
        assert self.x12 > 0., 'x12=%s' % self.x12
        assert self.x43 >= 0., 'x43=%s' % self.x43

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', 0)
        nspan = integer_or_blank(card, 4, 'nspan', 0)
        lspan = integer_or_blank(card, 5, 'lspan', 0)

        p1 = np.array([
            double_or_blank(card, 9, 'x1', 0.0),
            double_or_blank(card, 10, 'y1', 0.0),
            double_or_blank(card, 11, 'z1', 0.0)])
        x12 = double_or_blank(card, 12, 'x12', 0.)

        p4 = np.array([
            double_or_blank(card, 13, 'x4', 0.0),
            double_or_blank(card, 14, 'y4', 0.0),
            double_or_blank(card, 15, 'z4', 0.0)])
        x43 = double_or_blank(card, 16, 'x43', 0.)
        assert len(card) <= 17, 'len(CAERO4 card) = %i\ncard=%s' % (len(card), card)
        return CAERO4(eid, pid, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, lspan=lspan, comment=comment)

    def get_points(self):
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p4 = self.cp_ref.transform_node_to_global(self.p4)
        p2 = p1 + np.array([self.x12, 0., 0.])
        p3 = p4 + np.array([self.x43, 0., 0.])
        return [p1, p2, p3, p4]

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CAERO4 eid=%s' % self.eid

        self.pid = model.PAero(self.pid, msg=msg)  # links to PAERO4 (not added)
        self.cp = model.Coord(self.cp, msg=msg)

        self.cp_ref = self.cp
        self.pid_ref = self.pid

        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            self.lspan = model.AEFact(self.lspan, msg)
            self.lspan_ref = self.lspan
        self._init_ids()

    def safe_cross_reference(self, model):
        msg = ' which is required by CAERO4 eid=%s' % self.eid
        try:
            self.pid = model.PAero(self.pid, msg=msg)  # links to PAERO4 (not added)
            self.pid_ref = self.pid
        except KeyError:
            self.model.warning('cannot find PAERO4=%r' % self.pid)

        try:
            self.cp = model.Coord(self.cp, msg=msg)
            self.cp_ref = self.cp
        except KeyError:
            self.model.warning('cannot find CORDx=%r' % self.cp)

        if self.nspan == 0:
            assert isinstance(self.lspan, integer_types), self.lspan
            try:
                self.lspan = model.AEFact(self.lspan, msg)
                self.lspan_ref = self.lspan
            except KeyError:
                pass
        self._init_ids()

    def uncross_reference(self):
        self.pid = self.Pid()
        self.cp = self.Cp()
        del self.pid_ref, self.cp_ref

    def Cp(self):
        if isinstance(self.cp, integer_types):
            return self.cp
        return self.cp_ref.cid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid_ref.pid

    def _init_ids(self, dtype='int32'):
        """
        Fill `self.box_ids` with the sub-box ids. Shape is (nchord, nspan)
        """
        nchord, nspan = self.shape
        assert nchord >= 1, 'nchord=%s' % nchord
        assert nspan >= 1, 'nspan=%s' % nspan
        self.box_ids = np.zeros((nchord, nspan), dtype=dtype)

        try:
            for ichord in range(nchord):
                for ispan in range(nspan):
                    self.box_ids[ichord, ispan] = self.eid + ichord + ispan * nchord
        except OverflowError:
            if dtype == 'int64':
                msg = 'eid=%s ichord=%s ispan=%s nchord=%s' % (
                    self.eid, ichord, ispan, nchord)
                raise OverflowError(msg)
            self._init_ids(dtype='int64')

    @property
    def shape(self):
        """returns (nelements_nchord, nelements_span)"""
        nchord = 1

        if self.nspan == 0:
            y = self.lspan_ref.Di
            nspan = len(y) - 1
        else:
            nspan = self.nspan
        if nspan < 1:
            msg = 'CAERO4 eid=%s nspan=%s lspan=%s' % (
                self.eid, self.nspan, self.lspan)
            raise RuntimeError(msg)
        return nchord, nspan

    def get_npanel_points_elements(self):
        """
        Gets the number of sub-points and sub-elements for the CAERO card

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
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel
        """
        x = np.linspace(0., 1., num=2)  # nchord=1

        if self.nspan == 0:
            y = self.lspan_ref.Di
            nspan = len(y) - 1
        else:
            nspan = self.nspan
            y = np.linspace(0., 1., nspan + 1)

        if nspan < 1:
            msg = 'CAERO4 eid=%s nspan=%s lspan=%s' % (
                self.eid, self.nspan, self.lspan)
            raise RuntimeError(msg)
        return x, y

    def panel_points_elements(self):
        """
        Gets the sub-points and sub-elements for the CAERO card

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

    def get_LSpan(self):
        if isinstance(self.lspan, integer_types):
            return self.lspan
        return self.lspan_ref.sid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list
            The fields that define the card
        """
        list_fields = (['CAERO4', self.eid, self.Pid(), self.Cp(), self.nspan,
                        self.get_LSpan(), None, None, None,] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list
            The fields that define the card
        """
        cp = set_blank_if_default(self.Cp(), 0)

        nspan = set_blank_if_default(self.nspan, 0)
        lspan = set_blank_if_default(self.get_LSpan(), 0)

        list_fields = (['CAERO4', self.eid, self.Pid(), cp, nspan,
                        lspan, None, None, None,] + list(self.p1) + [self.x12] +
                       list(self.p4) + [self.x43])
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CAERO5(BaseCard):
    """
    Defines an aerodynamic macro element for Piston theory.

    +--------+------+------+-----+-------+-------+-------+--------+-------+
    |   1    |  2   |   3  |  4  |   5   |   6   |   7   |   8    |   9   |
    +========+======+======+=====+=======+=======+=======+========+=======+
    | CAERO5 | EID  | PID  | CP  | NSPAN | LSPAN | NTHRY | NTHICK |       |
    +--------+------+------+-----+-------+-------+-------+--------+-------+
    |        | X1   |  Y1  | Z1  |  X12  |  X4   |  Y4   |   Z4   |  X43  |
    +--------+------+------+-----+-------+-------+-------+--------+-------+

    +--------+------+------+-----+-------+-------+-------+--------+-------+
    | CAERO5 | 6000 | 6001 | 100 |       |  315  |   0   |   0    |       |
    +--------+------+------+-----+-------+-------+-------+--------+-------+
    |        | 0.0  |  0.0 | 0.0 |  1.0  |  0.2  |  1.0  |   0.   |  0.8  |
    +--------+------+------+-----+-------+-------+-------+--------+-------+
    """
    type = 'CAERO5'
    def __init__(self, eid, pid, p1, x12, p4, x43,
                 cp=0, nspan=0, lspan=0, ntheory=0, nthick=0,
                 comment=''):
        if comment:
            self.comment = comment

        #: Element identification number
        self.eid = eid

        #: Property identification number of a PAERO5 entry.
        self.pid = pid

        #: Coordinate system for locating point 1.
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.ntheory = ntheory
        self.nthick = nthick
        self.p1 = np.asarray(p1, dtype='float64')
        self.x12 = x12
        self.p4 = np.asarray(p4, dtype='float64')
        self.x43 = x43

        assert self.x12 > 0., 'x12=%s' % self.x12
        if not (self.nspan > 0 or self.lspan > 0):
            msg = 'nspan=%r or lspan=%r must be > 0' % (self.nspan, self.lspan)
            raise ValueError(msg)
        if not (self.x12 > 0.0 or self.x43 > 0.0):
            msg = 'x12=%r or x43=%r must be > 0.0' % (self.x12, self.x43)
            raise ValueError(msg)

    def validate(self):
        assert self.ntheory in [0, 1, 2], 'ntheory=%r' % self.ntheory

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CAERO5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', 0)
        nspan = integer_or_blank(card, 4, 'nspan', 0)
        lspan = integer_or_blank(card, 5, 'lspan', 0)
        ntheory = integer_or_blank(card, 6, 'ntheory', 0)
        nthick = integer_or_blank(card, 7, 'nthick')
        # 8 - blank
        p1 = np.array([
            double_or_blank(card, 9, 'x1', 0.0),
            double_or_blank(card, 10, 'y1', 0.0),
            double_or_blank(card, 11, 'z1', 0.0)])
        x12 = double(card, 12, 'x12')
        p4 = np.array([
            double_or_blank(card, 13, 'x4', 0.0),
            double_or_blank(card, 14, 'y4', 0.0),
            double_or_blank(card, 15, 'z4', 0.0)])
        x43 = double_or_blank(card, 16, 'x43', 0.0)
        assert len(card) <= 17, 'len(CAERO3 card) = %i\ncard=%s' % (len(card), card)
        return CAERO5(eid, pid, p1, x12, p4, x43,
                      cp=cp, nspan=nspan, lspan=lspan, ntheory=ntheory, nthick=nthick,
                      comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CAERO5 eid=%s' % self.eid
        self.pid = model.PAero(self.pid, msg=msg)
        self.pid_ref = self.pid
        self.cp = model.Coord(self.cp, msg=msg)
        self.cp_ref = self.cp
        if self.nspan == 0:
            self.lspan = model.AEFact(self.lspan, msg=msg)
            self.lspan_ref = self.lspan

    def safe_cross_reference(self, model):
        msg = ' which is required by CAERO5 eid=%s' % self.eid
        try:
            self.pid = model.PAero(self.pid, msg=msg)
            self.pid_ref = self.pid
        except KeyError:
            pass

        try:
            self.cp = model.Coord(self.cp, msg=msg)
            self.cp_ref = self.cp
        except KeyError:
            pass

        if self.nspan == 0:
            try:
                self.lspan = model.AEFact(self.lspan, msg=msg)
                self.lspan_ref = self.lspan
            except KeyError:
                pass

    def uncross_reference(self):
        self.pid = self.Pid()
        self.cp = self.Cp()
        if self.nspan == 0:
            self.lspan = self.LSpan()
            del self.lspan_ref
        del self.pid_ref, self.cp_ref

    def get_points(self):
        p1 = self.cp_ref.transform_node_to_global(self.p1)
        p4 = self.cp_ref.transform_node_to_global(self.p4)
        p2 = p1 + np.array([self.x12, 0., 0.])
        p3 = p4 + np.array([self.x43, 0., 0.])
        return [p1, p2, p3, p4]

    def get_npanel_points_elements(self):
        msg = 'CAERO5 eid=%s nspan=%s lspan=%s' % (
            self.eid, self.nspan, self.lspan)
        if self.nspan == 0:
            y = self.lspan_ref.Di
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

        msg = 'CAERO5 eid=%s nspan=%s lspan=%s' % (
            self.eid, self.nspan, self.lspan)
        if self.nspan == 0:
            y = self.lspan_ref.Di
            nspan = len(y) - 1
        else:
            nspan = self.nspan
            y = np.linspace(0., 1., nspan + 1)
        assert nspan >= 1, msg

        x = np.array([0., 1.], dtype='float64')
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
        L = np.linalg.norm(p4 - p1)

        # b = L * cos(Lambda)
        # ci = L * sin(Lambda)

        if self.ntheory == 0:
            # piston theory
            pass
        elif self.ntheory == 1:
            raise NotImplementedError('ntheory=%s' % self.ntheory)
            #gamma = 1.4
            #Lambda = 0.
            #c1 = 1.
            #secL = 1 / np.cos(Lambda)
            #secL2 = secL ** 2
            #ma2_secL2 = mach ** 2 - secL2
            #c1 = mach / ma2_secL2 ** 0.5
            #c2 = (mach ** 4 * (gamma + 1) - 4 * secL2 * ma2_secL2) / (4 * ma2_secL2 ** 2)
        else:
            gamma = 1.4

            # the advance in x
            ci = p4[0] - p1[0]

            # sweep angle
            Lambda = np.arcsin(ci / L)

            secL = 1 / np.cos(Lambda)
            secL2 = secL ** 2
            ma2_secL2 = mach ** 2 - secL2
            c1 = mach / ma2_secL2 ** 0.5
            c2 = (mach ** 4 * (gamma + 1) - 4 * secL2 * ma2_secL2) / (4 * ma2_secL2 ** 2)

    def Cp(self):
        if isinstance(self.cp, integer_types):
            return self.cp
        return self.cp_ref.cid

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid_ref.pid

    def LSpan(self):
        if isinstance(self.lspan, integer_types):
            return self.lspan
        return self.lspan_ref.sid

    def repr_fields(self):
        """
        Gets the fields in their simplified form

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

def elements_from_quad(nx, ny):
    assert nx > 1
    assert ny > 1

    nelements = (nx - 1) * (ny - 1)
    npoints = nx * ny

    # create a matrix with the point counter
    ipoints = np.arange(npoints, dtype='int32').reshape((nx, ny))

    # move around the CAERO quad and apply ipoints
    elements = np.zeros((nelements, 4), dtype='int32')
    elements[:, 0] = ipoints[:-1, :-1].ravel()  # (i,  j  )
    elements[:, 1] = ipoints[1:, :-1].ravel()   # (i+1,j  )
    elements[:, 2] = ipoints[1:, 1:].ravel()    # (i+1,j+1)
    elements[:, 3] = ipoints[:-1, 1:].ravel()   # (i,j+1  )
    return elements

def points_elements_from_quad_points(p1, p2, p3, p4, x, y):
    """
    Creates nodes and elements in a structured grid given 4 points.
    Used to make an CAERO1 panel.

    Parameters
    ----------
    p1 : (3, ) float ndarray
        leading edge root
    p2 : (3, ) float ndarray
        trailing edge root
    p3 : (3, ) float ndarray
        trailing edge tip
    p4 : (3, ) float ndarray
        leading edge tip
    x : (nchord, ) float ndarray
        points in the chordwise direction in percentage of the chord
    y : (nspan, ) float ndarray
        points in the spanwise direction in percentage of the span

    Returns
    -------
    points (nchord, nspan) float ndarray; might be backwards???
        the points
    elements (nquads, 4) int ndarray
        series of quad elements
        nquads = (nchord-1) * (nspan-1)
    """
    nx = x.shape[0]
    ny = y.shape[0]

    elements = elements_from_quad(nx, ny)
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
    xv = np.repeat(x, ny, axis=1).reshape(npoints, 1)
    yv = np.repeat(y, nx, axis=0).reshape(npoints, 1)

    # calculate the points a and b xv% along the chord
    a = xv * p2 + (1 - xv) * p1
    b = xv * p3 + (1 - xv) * p4

    # calculate the point yv% along the span
    points = yv * b + (1 - yv) * a
    assert points.shape == (npoints, 3), 'npoints=%s shape=%s' % (npoints, str(points.shape))

    # create a matrix with the point counter
    ipoints = np.arange(npoints, dtype='int32').reshape((nx, ny))

    return points, elements


class PAERO5(BaseCard):
    type = 'PAERO5'
    def __init__(self, pid, caoci,
                 nalpha=0, lalpha=0,
                 nxis=0, lxis=0,
                 ntaus=0, ltaus=0, comment=''):
        """
        +--------+-------+--------+--------+---------+-------+-------+-------+
        |   1    |   2   |    3   |   4    |    5    |   6   |   7   |   8   |
        +========+=======+========+========+=========+=======+=======+=======+
        | PAERO5 |  PID  | NALPHA | LALPHA |  NXIS   | LXIS  | NTAUS | LTAUS |
        +--------+-------+--------+--------+---------+-------+-------+-------+
        |        | CAOC1 | CAOC2  | CAOC3  |  CAOC4  | CAOC5 |       |       |
        +--------+-------+--------+--------+---------+-------+-------+-------+

        +--------+-------+--------+--------+---------+-------+-------+-------+
        | PAERO5 | 7001  |   1    |  702   |    1    | 701   |   1   |  700  |
        +--------+-------+--------+--------+---------+-------+-------+-------+
        |        |  0.0  |  0.0   |  5.25  | 3.99375 |  0.0  |       |       |
        +--------+-------+--------+--------+---------+-------+-------+-------+
        """
        if comment:
            self.comment = comment
        self.pid = pid
        self.nalpha = nalpha
        self.lalpha = lalpha

        # number of dimensionless chord coordinates in zeta ()
        self.nxis = nxis
        # ID of AEFACT that lists zeta
        self.lxis = lxis

        # number of dimensionless thickess coordinates in tau
        self.ntaus = ntaus
        # ID of AEFACT that lists thickness ratios (t/c)
        self.ltaus = ltaus

        # ca/c - control surface chord / strip chord
        self.caoci = np.array(caoci, dtype='float64')

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAERO5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'property_id')
        nalpha = integer_or_blank(card, 2, 'nalpha', default=0)
        lalpha = integer_or_blank(card, 3, 'lalpha', default=0)

        nxis = integer_or_blank(card, 4, 'nxis', default=0)
        lxis = integer_or_blank(card, 5, 'lxis', default=0)

        ntaus = integer_or_blank(card, 6, 'ntaus', default=0)
        ltaus = integer_or_blank(card, 7, 'ltaus', default=0)

        caoci = []
        for n, i in enumerate(range(9, len(card))):
            ca = double(card, i, 'ca/ci_%i' % (n+1))
            caoci.append(ca)
        return PAERO5(pid, caoci,
                      nalpha=nalpha, lalpha=lalpha, nxis=nxis, lxis=lxis,
                      ntaus=ntaus, ltaus=ltaus,
                      comment=comment)
    @property
    def lxis_id(self):
        return self.lxis if isinstance(self.lxis, integer_types) else self.lxis_ref.sid

    @property
    def ltaus_id(self):
        return self.ltaus if isinstance(self.ltaus, integer_types) else self.ltaus_ref.sid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.lxis = model.AEFact(self.lxis_id)
        self.ltaus = model.AEFact(self.ltaus_id)

        self.ltaus_ref = self.ltaus
        self.lxis_ref = self.lxis

    def safe_cross_reference(self, model):
        try:
            self.lxis = model.AEFact(self.lxis_id)
            self.lxis_ref = self.lxis
        except KeyError:
            pass

        try:
            self.ltaus = model.AEFact(self.ltaus_id)
            self.ltaus_ref = self.ltaus
        except KeyError:
            pass

    def uncross_reference(self):
        self.lxis = self.lxis_id
        self.ltaus = self.ltaus_id
        del self.ltaus_ref, self.lxis_ref

    def raw_fields(self):
        list_fields = ['PAERO5', self.pid, self.nalpha, self.lalpha, self.nxis,
                       self.lxis_id, self.ntaus, self.ltaus_id] + list(self.caoci)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

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

class DIVERG(BaseCard):
    """
    +--------+-----+--------+----+----+----+----+----+----+
    |   1    |  2  |   3    | 4  | 5  | 6  | 7  | 8  | 9  |
    +========+=====+========+====+====+====+====+====+====+
    | DIVERG | SID | NROOT  | M1 | M2 | M3 | M4 | M5 | M6 |
    +--------+-----+--------+----+----+----+----+----+----+
    |        |  M7 |  etc.  |    |    |    |    |    |    |
    +--------+-----+--------+----+----+----+----+----+----+

    Attributes
    ----------
    sid : int
        The name.
    nroots : int
        the number of roots
    machs : List[float, ..., float]
        list of Mach numbers
    """
    type = 'DIVERG'
    def __init__(self, sid, nroots, machs, comment=''):
        """
        Creates an DIVERG card, which is used in divergence
        analysis (SOL 144).

        Parameters
        ----------
        sid : int
            The name
        nroots : int
            the number of roots
        machs : List[float, ..., float]
            list of Mach numbers
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.sid = sid
        self.nroots = nroots
        self.machs = machs

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DIVERG card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        nroots = integer(card, 2, 'nroot')
        j = 1
        machs = []
        for i in range(3, len(card)):
            mach = double(card, i, 'Mach_%i' % j)
            machs.append(mach)
            j += 1
        return DIVERG(sid, nroots, machs, comment=comment)

    #def cross_reference(self, model):
        #pass

    #def uncross_reference(self):
        #pass

    def raw_fields(self):
        list_fields = ['DIVERG', self.sid, self.nroots] + list(self.machs)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class FLFACT(BaseCard):
    """
    +--------+-----+----+------+----+----+----+----+----+
    |   1    |  2  |  3 |   4  | 5  | 6  | 7  | 8  | 9  |
    +========+=====+====+======+====+====+====+====+====+
    | FLFACT | SID | F1 | F2   | F3 | F4 | F5 | F6 | F7 |
    +--------+-----+----+------+----+----+----+----+----+
    |        | F8  | F9 | etc. |    |    |    |    |    |
    +--------+-----+----+------+----+----+----+----+----+

    +--------+-----+----+------+-----+
    |   1    |  2  |  3 |   4  | 5   |
    +========+=====+====+======+=====+
    | FLFACT | 97  | .3 |  .7  | 3.5 |
    +--------+-----+----+------+-----+

    # delta quantity approach

    +--------+-----+-------+------+-------+----+--------+
    |   1    |  2  |  3    |   4  |   5   | 6  |     7  |
    +========+=====+=======+======+=======+====+========+
    | FLFACT | SID | F1    | THRU | FNF   | NF |  FMID  |
    +--------+-----+-------+------+-------+----+--------+
    | FLFACT | 201 | 0.200 | THRU | 0.100 | 11 | 0.1333 |
    +--------+-----+-------+------+-------+----+--------+
    """
    type = 'FLFACT'

    def __init__(self, sid, factors, comment=''):
        """
        Creates an FLFACT card

        Parameters
        ----------
        sid : int
            the id of a density, reduced_frequency, mach, or velocity table
            the FLUTTER card defines the meaning
        factors : varies
            values : List[float, ..., float]
                list of factors
            List[f1, THRU, fnf, nf, fmid]
                f1 : float
                    first value
                THRU : str
                    the word THRU
                nf : float
                    second value
                fmid : float; default=(f1 + fnf) / 2.
                    the mid point to bias the array
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.sid = sid
        #self.f1 = f1
        #self.fnf = fnf
        #self.nf = nf
        #self.fmid = fmid

        # the dumb string_types thing is because we also get floats
        if len(factors) > 1 and isinstance(factors[1], string_types) and factors[1] == 'THRU':
            #msg = 'embedded THRUs not supported yet on FLFACT card\n'
            nfactors = len(factors)
            if nfactors == 4:
                (f1, _thru, fnf, nf) = factors
                fmid = (f1 + fnf) / 2.
            elif nfactors == 5:
                (f1, thru, fnf, nf, fmid) = factors
            else:
                raise RuntimeError('factors must be length 4/5; factors=%s' % factors)
            i = np.linspace(0, nf, nf, endpoint=False) + 1
            factors = (
                (f1*(fnf - fmid) * (nf-i) + fnf * (fmid - f1) * (i-1)) /
                (   (fnf - fmid) * (nf-i) +       (fmid - f1) * (i-1))
            )
        self.factors = np.asarray(factors)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an FLFACT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        assert len(card) > 2, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
        field3 = double_string_or_blank(card, 3, 'THRU')
        if field3 is None:
            f1 = double(card, 2, 'f1')
            factors = [f1]
            assert len(card) == 3, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
        elif isinstance(field3, float):
            factors = fields(double, card, 'factors', i=2, j=len(card))
        elif isinstance(field3, string_types) and field3 == 'THRU':
            f1 = double(card, 2, 'f1')
            fnf = double(card, 4, 'fnf')
            nf = integer(card, 5, 'nf')
            fmid_default = (f1 + fnf) / 2.
            fmid = double_or_blank(card, 6, 'fmid', fmid_default)
            assert len(card) <= 7, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
            i = np.linspace(0, nf, nf, endpoint=False) + 1
            factors = (
                (f1*(fnf - fmid) * (nf-i) + fnf * (fmid - f1) * (i-1)) /
                (   (fnf - fmid) * (nf-i) +       (fmid - f1) * (i-1))
            )
        else:
            raise SyntaxError('expected a float or string for FLFACT field 3; value=%r' % field3)
        return FLFACT(sid, factors, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        factors = data[1:]
        return FLFACT(sid, factors, comment=comment)

    def max(self):
        return self.factors.max()

    def min(self):
        return self.factors.min()

    #def uncross_reference(self):
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
        """
        list_fields = ['FLFACT', self.sid] + list(self.factors)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        else:
            return self.comment + print_card_16(card)


class FLUTTER(BaseCard):
    """
    Defines data needed to perform flutter analysis.

    +---------+-----+--------+------+------+-------+-------+-------------+------+
    |    1    |  2  |   3    |  4   |  5   |   6   |   7   |      8      |  9   |
    +=========+=====+========+======+======+=======+=======+=============+======+
    | FLUTTER | SID | METHOD | DENS | MACH | RFREQ | IMETH | NVALUE/OMAX | EPS  |
    +---------+-----+--------+------+------+-------+-------+-------------+------+
    | FLUTTER | 19  |   K    | 119  | 219  | 319   |   S   |      5      | 1.-4 |
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

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : int/float/str
            the value for the appropriate field
        """
        if n == 7:
            if self.method in ['K', 'KE']:
                value = self.nvalue
            elif self.method in ['PKS', 'PKNLS']:
                value = self.omax
            else:
                value = self.nvalue
            return value
        else:
            raise KeyError('Field %r=%r is an invalid FLUTTER entry.' % (n, value))

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the FLUTTER card

        Parameters
        ----------
        n : int
            the field number to update
        value : int/float/str
            the value for the appropriate field
        """
        if n == 7:
            if self.method in ['K', 'KE']:
                self.nvalue = value
            elif self.method in ['PKS', 'PKNLS']:
                self.omax = value
            else:
                self.nvalue = value
        else:
            raise KeyError('Field %r=%r is an invalid FLUTTER entry.' % (n, value))

    def __init__(self, sid, method, density, mach, reduced_freq_velocity,
                 imethod='L', nvalue=None, omax=None, epsilon=1.0e-3, comment=''):
        """
        Creates a FLUTTER card, which is required for a flutter (SOL 145)
        analysis.

        Parameters
        ----------
        sid : int
            flutter id
        method : str
            valid methods = [K, KE,
                             PKS, PKNLS, PKNL, PKE]
        density : int
            defines a series of air densities in units of mass/volume
            PARAM,WTMASS does not affect this
            AERO affects this
            references an FLFACT id
        mach : int
            defines a series of the mach numbers
            references an FLFACT id
        reduced_freq_velocity : int
            Defines a series of either:
               1) reduced frequencies - K, KE
               2) velocities - PK, PKNL, PKS, PKNLS
            depending on the method chosen.
            references an FLFACT id
        imethod : str; default='L'
            Choice of interpolation method for aerodynamic matrix interpolation.
            imethods :
               1) L - linear
               2) S - surface
               3) TCUB - termwise cubic
        nvalue : int
            Number of eigenvalues beginning with the first eigenvalue for
            output and plots
        omax : float
            For the PKS and PKNLS methods, OMAX specifies the maximum frequency, in
            Hz., to be used in he flutter sweep.
            MSC only.
        epsilon : float; default=1.0e-3
            Convergence parameter for k. Used in the PK and PKNL methods only
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.sid = sid
        if method in ['PK', 'PKNL', 'PKNLS']:
            imethod = 'L'
        else:
            assert imethod in ['S', 'L', None], imethod
        self.method = method
        self.density = density
        self.mach = mach

        # KFREQ - K, KE
        # VEL - PK, PKNL, PKS, PKNLS
        self.reduced_freq_velocity = reduced_freq_velocity

        #
        self.imethod = imethod
        self.nvalue = nvalue
        self.omax = omax
        self.epsilon = epsilon

    def validate(self):
        if self.method not in ['K', 'KE', 'PK', 'PKNL', 'PKS', 'PKNLS']:
            msg = 'method = %r; allowed=[K, KE, PKS, PKNLS, PKNL, PK]' % self.method
            raise ValueError(msg)
        if self.imethod not in ['L', 'S', 'TCUB']:
            msg = 'imethod = %r; allowed=[L, S, TCUB]' % self.imethod
            raise ValueError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FLUTTER card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        method = string(card, 2, 'method (K, KE, PKS, PKNLS, PKNL, PK)')
        density_id = integer(card, 3, 'density')
        mach_id = integer(card, 4, 'mach')
        reduced_freq_velocity_id = integer(card, 5, 'reduced_freq_velocity')

        if method in ['K', 'KE']:
            imethod = string_or_blank(card, 6, 'imethod', 'L')
            nvalue = integer_or_blank(card, 7, 'nvalue')
            omax = None
            assert imethod in ['L', 'S', 'TCUB'], 'imethod = %s' % imethod  # linear-surface
        elif method in ['PKS', 'PKNLS']:
            imethod = None
            nvalue = None
            omax = double_or_blank(card, 7, 'omax')
        elif method == 'PKNL':
            nvalue = integer_or_blank(card, 7, 'nvalue')
            omax = None
            imethod = None
        elif method == 'PK':
            nvalue = integer_or_blank(card, 7, 'nvalue')
            omax = None
            imethod = None
        else:
            raise NotImplementedError('FLUTTER method=%r' % method)

        assert method in ['K', 'KE', 'PK', 'PKS', 'PKNL', 'PKNLS', None], method
        epsilon = double_or_blank(card, 8, 'epsilon', 1e-3)  # not defined in QRG
        assert len(card) <= 9, 'len(FLUTTER card) = %i\ncard=%s' % (len(card), card)
        return FLUTTER(sid, method, density_id, mach_id, reduced_freq_velocity_id,
                       imethod=imethod, nvalue=nvalue, omax=omax,
                       epsilon=epsilon, comment=comment)

    @property
    def headers(self):
        headers = ['density', 'mach']
        if self.method in ['PK', 'PKS', 'PKNL', 'PKNLS']:
            headers.append('velocity')
        elif self.method in ['K', 'KE']:
            headers.append('reduced_frequency')
        else:
            raise NotImplementedError('FLUTTER method=%r' % self.method)
        return headers

    @classmethod
    def add_op2_data(cls, data, comment=''):
        assert len(data) == 8, 'FLUTTER = %s' % data
        sid = data[0]
        method = data[1]
        density = data[2]
        mach = data[3]
        reduced_freq_velocity = data[4]
        method = data[5]
        imethod = data[6]
        nvalue = data[7]
        omax = data[8]
        epsilon = None
        return FLUTTER(sid, method, density, mach, reduced_freq_velocity,
                       imethod, nvalue, omax,
                       epsilon, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by FLUTTER sid=%s' % self.sid
        self.density = model.FLFACT(self.density, msg=msg)
        self.density_ref = self.density
        self.mach = model.FLFACT(self.mach, msg=msg)
        self.mach_ref = self.mach
        self.reduced_freq_velocity = model.FLFACT(self.reduced_freq_velocity, msg=msg)
        self.reduced_freq_velocity_ref = self.reduced_freq_velocity

    def safe_cross_reference(self, model):
        msg = ' which is required by FLUTTER sid=%s' % self.sid
        try:
            self.density = model.FLFACT(self.density, msg=msg)
            self.density_ref = self.density
        except KeyError:
            pass
        try:
            self.mach = model.FLFACT(self.mach, msg=msg)
            self.mach_ref = self.mach
        except KeyError:
            pass
        try:
            self.reduced_freq_velocity = model.FLFACT(self.reduced_freq_velocity, msg=msg)
            self.reduced_freq_velocity_ref = self.reduced_freq_velocity
        except KeyError:
            pass

    def uncross_reference(self):
        self.density = self.get_density()
        self.mach = self.get_mach()
        self.reduced_freq_velocity = self.get_rfreq_vel()
        del self.density_ref, self.mach_ref, self.reduced_freq_velocity_ref

    def get_density(self):
        if isinstance(self.density, integer_types):
            return self.density
        return self.density_ref.sid

    def get_mach(self):
        if isinstance(self.mach, integer_types):
            return self.mach
        return self.mach_ref.sid

    def get_rfreq_vel(self):
        if isinstance(self.reduced_freq_velocity, integer_types):
            return self.reduced_freq_velocity
        return self.reduced_freq_velocity_ref.sid

    def _get_raw_nvalue_omax(self):
        if self.method in ['K', 'KE']:
            #assert self.imethod in ['L', 'S'], 'imethod = %s' % self.imethod
            return(self.imethod, self.nvalue)
        elif self.method in ['PKS', 'PKNLS']:
            return(self.imethod, self.omax)
        else:
            return(self.imethod, self.nvalue)

    def _repr_nvalue_omax(self):
        if self.method in ['K', 'KE']:
            imethod = set_blank_if_default(self.imethod, 'L')
            #assert self.imethod in ['L', 'S'], 'imethod = %s' % self.imethods
            return (imethod, self.nvalue)
        elif self.method in ['PKS', 'PKNLS']:
            return(self.imethod, self.omax)
        else:
            return(self.imethod, self.nvalue)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
        """
        (imethod, nvalue) = self._get_raw_nvalue_omax()
        list_fields = ['FLUTTER', self.sid, self.method, self.get_density(),
                       self.get_mach(), self.get_rfreq_vel(), imethod, nvalue, self.epsilon]
        return list_fields

    #def repr_fields(self):
        #(imethod, nvalue) = self._get_raw_nvalue_omax()
        #list_fields = ['FLUTTER', self.sid, self.method, self.get_density(), self.get_mach(),
        #          self.get_rfreq_vel(), imethod, nvalue, self.epsilon]
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
    +======+=====+=======+=====+=====+======+
    | GUST | SID | DLOAD | WG  | X0  |  V   |
    +------+-----+-------+-----+-----+------+
    | GUST | 133 |   61  | 1.0 | 0.  | 1.+4 |
    +------+-----+-------+-----+-----+------+
    """
    type = 'GUST'
    _field_map = {
        1: 'sid', 2:'dload', 3:'wg', 4:'x0', 5:'V',
    }

    def __init__(self, sid, dload, wg, x0, V=None, comment=''):
        if comment:
            self.comment = comment
        self.sid = sid
        self.dload = dload
        self.wg = wg
        self.x0 = x0
        self.V = V

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a GUST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        dload = integer(card, 2, 'dload')
        wg = double(card, 3, 'wg')
        x0 = double(card, 4, 'x0')
        V = double_or_blank(card, 4, 'V')
        assert len(card) <= 6, 'len(GUST card) = %i\ncard=%s' % (len(card), card)
        return GUST(sid, dload, wg, x0, V=V, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        sid = data[0]
        dload = data[1]
        wg = data[2]
        x0 = data[3]
        V = data[4]
        assert len(data) == 5, 'data = %s' % data
        return GUST(sid, dload, wg, x0, V, comment=comment)

    #def Angle(self):
        #angle = self.wg*self.t*(t-(x-self.x0)/self.V) # T is the tabular
        #return angle

    #def uncross_reference(self):
        #pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
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
    +=========+====+====+====+====+====+====+====+====+
    | MKAERO1 | m1 | m2 | m3 | m4 | m5 | m6 | m7 | m8 |
    +---------+----+----+----+----+----+----+----+----+
    |         | k1 | k2 | k3 | k4 | k5 | k6 | k7 | k8 |
    +---------+----+----+----+----+----+----+----+----+
    """
    type = 'MKAERO1'

    def __init__(self, machs, reduced_freqs, comment=''):
        """
        Creates an MKAERO1 card, which defines a set of mach and
        reduced frequencies.

        Parameters
        ----------
        machs : List[float]
            series of Mach numbers
        reduced_freqs : List[float]
            series of reduced frequencies
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.machs = np.unique(machs)
        self.reduced_freqs = np.unique(reduced_freqs)

    def validate(self):
        if len(self.machs) == 0:
            msg = 'MKAERO1; nmachs=%s machs=%s' % (len(self.machs), self.machs)
            raise ValueError(msg)
        if len(self.reduced_freqs) == 0:
            msg = 'MKAERO1; nrfreqs=%s rfreqs=%s' % (len(self.reduced_freqs), self.reduced_freqs)
            raise ValueError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an MKAERO1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        list_fields = [interpret_value(field) for field in card[1:]]
        nfields = len(list_fields) - 8
        machs = []
        reduced_freqs = []
        for i in range(1, 1 + nfields):
            machs.append(double_or_blank(card, i, 'mach'))
            reduced_freqs.append(double_or_blank(card, i + 8, 'rFreq'))
        machs = wipe_empty_fields(machs)
        reduced_freqs = wipe_empty_fields(reduced_freqs)
        return MKAERO1(machs, reduced_freqs, comment=comment)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
        """
        #list_fields = ['MKAERO1']
        #for (i, mach, rfreq) in zip(count(), self.machs, self.reduced_freqs):
        #    list_fields += [mach, rfreq]

        # kind of a hack because there isn't a good way to do this for
        # duplicately-defined MKAERO1s
        machs = [None] * max(8, len(self.machs))
        freqs = [None] * max(8, len(self.reduced_freqs))
        for i, mach in enumerate(self.machs):
            machs[i] = mach
        for i, freq in enumerate(self.reduced_freqs):
            freqs[i] = freq
        list_fields = ['MKAERO1'] + machs + freqs
        return list_fields

    def write_card(self, size=8, is_double=False):
        #cards = []
        nmachs = len(self.machs)
        nreduced_freqs = len(self.reduced_freqs)
        if nmachs > 8 or nreduced_freqs > 8:
            #cards = []
            mach_sets = []
            rfreq_sets = []
            imach = 0
            ifreq = 0
            while imach < nmachs:
                mach_sets.append(self.machs[imach:imach+8])
                imach += 8
            while ifreq < nreduced_freqs:
                rfreq_sets.append(self.reduced_freqs[ifreq:ifreq+8])
                ifreq += 8
            msg = self.comment

            #print('mach_sets = %s' % mach_sets)
            #print('rfreq_sets = %s' % rfreq_sets)
            for mach_set in mach_sets:
                for rfreq_set in rfreq_sets:
                    msg += MKAERO1(mach_set, rfreq_set).write_card(
                        size=size, is_double=is_double)
            return msg

        machs = [None] * 8
        reduced_freqs = [None] * 8
        cards = []
        if not 0 < len(self.machs) <= 8:
            msg = 'MKAERO1; nmachs=%s machs=%s' % (len(self.machs), self.machs)
            raise ValueError(msg)
        if not 0 < len(self.reduced_freqs) <= 8:
            msg = 'MKAERO1; nrfreqs=%s rfreqs=%s' % (len(self.reduced_freqs), self.reduced_freqs)
            raise ValueError(msg)

        for i, mach in zip(count(), self.machs):
            machs[i] = mach
        for i, rfreq in zip(count(), self.reduced_freqs):
            reduced_freqs[i] = rfreq
        return self.comment + print_card_8(['MKAERO1'] + machs + reduced_freqs)

    def __repr__(self):
        return self.write_card()

class MKAERO2(BaseCard):
    """
    Provides a table of Mach numbers (m) and reduced frequencies (k) for
    aerodynamic matrix calculation.

    +---------+----+----+----+----+----+----+----+----+
    |    1    | 2  | 3  | 4  | 5  | 6  | 7  | 8  | 9  |
    +=========+====+====+====+====+====+====+====+====+
    | MKAERO2 | m1 | k1 | m2 | k2 | m3 | k3 | m4 | k4 |
    +---------+----+----+----+----+----+----+----+----+
    """
    type = 'MKAERO2'

    def __init__(self, machs, reduced_freqs, comment=''):
        """
        Creates an MKAERO2 card, which defines a set of mach and
        reduced frequency pairs.

        Parameters
        ----------
        machs : List[float]
            series of Mach numbers
        reduced_freqs : List[float]
            series of reduced frequencies
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.machs = machs
        self.reduced_freqs = reduced_freqs

    def validate(self):
        if len(self.machs) != len(self.reduced_freqs):
            msg = 'MKAERO2; len(machs)=%s len(rfreqs)=%s; should be the same' % (
                len(self.machs), len(self.reduced_freqs))
            raise ValueError(msg)

        if len(self.machs) == 0:
            msg = 'MKAERO2; nmachs=%s machs=%s' % (len(self.machs), self.machs)
            raise ValueError(msg)
        if len(self.reduced_freqs) == 0:
            msg = 'MKAERO2; nrfreqs=%s rfreqs=%s' % (len(self.reduced_freqs), self.reduced_freqs)
            raise ValueError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds an MKAERO2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        list_fields = card.fields(1)
        nfields = len(list_fields)
        machs = []
        reduced_freqs = []
        for i in range(1, 1 + nfields, 2):
            machs.append(double(card, i, 'mach'))
            reduced_freqs.append(double(card, i + 1, 'rFreq'))
        return MKAERO2(machs, reduced_freqs, comment=comment)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
        """
        list_fields = ['MKAERO2']
        for (i, mach, rfreq) in zip(count(), self.machs, self.reduced_freqs):
            list_fields += [mach, rfreq]
        return list_fields

    def write_card(self, size=8, is_double=False):
        cards = []
        list_fields = ['MKAERO2']
        nvalues = 0
        for mach, rfreq in zip(self.machs, self.reduced_freqs):
            list_fields += [mach, rfreq]
            nvalues += 1
            if nvalues == 4:
                cards.append(print_card_8(list_fields))
                list_fields = ['MKAERO2']
                nvalues = 0
        if nvalues:
            cards.append(print_card_8(list_fields))
        return self.comment + ''.join(cards)

    def __repr__(self):
        return self.write_card()


class MONPNT1(BaseCard):
    """
    +---------+---------+------+-----+-----+-------+------+----+----+
    |    1    |    2    |  3   |  4  |  5  |   6   |   7  | 8  | 9  |
    +=========+=========+======+=====+=====+=======+======+====+====+
    | MONPNT1 |  NAME   |                   LABEL                   |
    +---------+---------+------+-----+-----+-------+------+----+----+
    |         |  AXES   | COMP | CP  |  X  |   Y   |   Z  | CD |    |
    +---------+---------+------+-----+-----+-------+------+----+----+

    +---------+---------+------+-----+-----+-------+------+----+----+
    |    1    |    2    |  3   |  4  |  5  |   6   |   7  | 8  | 9  |
    +=========+=========+======+=====+=====+=======+======+====+====+
    | MONPNT1 | WING155 |    Wing Integrated Load to Butline 155    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    |         |    34   | WING |     | 0.0 | 155.0 | 15.0 |    |    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    """
    type = 'MONPNT1'
    def __init__(self, name, label, axes, comp, xyz, cp=0, cd=None, comment=''):
        """
        Creates a MONPNT1 card

        Parameters
        ----------
        name : str
            Character string of up to 8 characters identifying the
            monitor point
        label : str
            A string comprising no more than 56 characters
            that identifies and labels the monitor point.
        axes : str
            components {1,2,3,4,5,6}
        comp : str
            name of the AECOMP/AECOMPL entry
        xyz : List[float, float, float]; default=None
            The coordinates in the CP coordinate system about which the
            loads are to be monitored.
            None : [0., 0., 0.]
        cp : int, CORDx; default=0
           int : coordinate system
        cd : int; default=None -> cp
            the coordinate system for load outputs
        comment : str; default=''
            a comment for the card

        CD - MSC specific field
        """
        if comment:
            self.comment = comment
        if cd is None:
            cd = cp
        self.name = name
        self.label = label
        self.axes = axes
        self.comp = comp
        self.cp = cp
        self.xyz = xyz
        self.cd = cd
        assert len(xyz) == 3, xyz

    @classmethod
    def add_card(cls, card, comment=''):
        name = string(card, 1, 'name')

        label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        label = ''.join(label_fields).strip()
        assert len(label) <= 56, label

        axes = parse_components(card, 9, 'axes')
        comp = string(card, 10, 'comp')
        cp = integer_or_blank(card, 11, 'cp', 0)
        xyz = [
            double_or_blank(card, 12, 'x', default=0.0),
            double_or_blank(card, 13, 'y', default=0.0),
            double_or_blank(card, 14, 'z', default=0.0),
        ]
        cd = integer_or_blank(card, 15, 'cd', cp)
        return MONPNT1(name, label, axes, comp, xyz, cp=cp, cd=cd, comment=comment)

    #def cross_reference(self):
        #pass

    #def uncross_reference(self):
        #pass

    def raw_fields(self):
        list_fields = [
            'MONPNT1', self.name, self.label.strip(), self.axes, self.comp,
            self.cp,] + self.xyz + [self.cd]
        return list_fields

    def write_card(self, size=8, is_double=False):
        cp = self.cp
        x, y, z = self.xyz
        cd = self.cd
        if cd == cp:
            cd = ''
        msg = 'MONPNT1 %-8s%s\n' % (self.name, self.label)
        msg += '        %-8s%-8s%-8s%-8s%-8s%-8s%-8s\n' % (
            self.axes, self.comp, cp,
            print_float_8(x), print_float_8(y), print_float_8(z),
            cd)
        #card = self.repr_fields()
        return self.comment + msg

    def __repr__(self):
        return self.write_card()

class MONPNT2(BaseCard):
    """MSC Nastran specific card"""
    type = 'MONPNT2'
    def __init__(self, name, label, table, Type, nddl_item, eid, comment=''):
        if comment:
            self.comment = comment
        self.name = name
        self.label = label
        self.table = table
        self.Type = Type
        self.nddl_item = nddl_item
        self.eid = eid

    def validate(self):
        assert self.table in ['STRESS', 'FORCE', 'STRAIN'], self.table

    @classmethod
    def add_card(cls, card, comment=''):
        name = string(card, 1, 'name')

        label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        label = ''.join(label_fields).strip()
        assert len(label) <= 56, label

        table = string(card, 9, 'table')
        Type = string(card, 10, 'type')
        nddl_item = integer_or_blank(card, 11, 'nddl_item')
        eid = integer_or_blank(card, 12, 'eid')
        return MONPNT2(name, label, table, Type, nddl_item, eid, comment=comment)

    #def uncross_reference(self):
        #pass

    def raw_fields(self):
        list_fields = [
            'MONPNT2', self.name, self.label.strip(),
            self.table, self.Type, self.nddl_item, self.eid]
        return list_fields

    def write_card(self, size=8, is_double=False):
        msg = 'MONPNT2 %-8s%s\n' % (self.name, self.label)
        msg += ('        %-8s%-8s%-8s%-8s\n' % (
            self.table, self.Type, self.nddl_item, self.eid
        ))
        #card = self.repr_fields()
        return self.comment + msg.rstrip() + '\n'

    def __repr__(self):
        return self.write_card()

class MONPNT3(BaseCard):
    """MSC Nastran specific card"""
    type = 'MONPNT3'
    def __init__(self, name, label, axes, grid_set, elem_set, xyz,
                 cp=0, cd=None, xflag=None, comment=''):
        if comment:
            self.comment = comment
        if cd is None:
            cd = cp
        self.name = name
        self.label = label
        self.axes = axes
        #self.comp = comp
        self.grid_set = grid_set
        self.elem_set = elem_set
        self.xyz = xyz
        self.xflag = xflag
        self.cp = cp
        self.cd = cd

    @classmethod
    def add_card(cls, card, comment=''):
        name = string(card, 1, 'name')

        label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        label = ''.join(label_fields).strip()
        assert len(label) <= 56, label

        axes = parse_components(card, 9, 'axes')
        grid_set = integer(card, 10, 'grid_set')
        elem_set = integer_or_blank(card, 11, 'elem_set')
        cp = integer_or_blank(card, 12, 'cp', 0)
        xyz = [
            double_or_blank(card, 13, 'x', default=0.0),
            double_or_blank(card, 14, 'y', default=0.0),
            double_or_blank(card, 15, 'z', default=0.0),
        ]
        xflag = string_or_blank(card, 16, 'xflag')
        cd = integer_or_blank(card, 17, 'cd', cp)
        return MONPNT3(name, label, axes, grid_set, elem_set, xyz,
                       cp=cp, cd=cd, xflag=xflag, comment=comment)

    #def uncross_reference(self):
        #pass

    def raw_fields(self):
        list_fields = [
            'MONPNT3', self.name, self.label.strip(),
            self.axes, self.grid_set, self.elem_set, self.cp] + self.xyz + [self.xflag, self.cd]
        return list_fields

    def write_card(self, size=8, is_double=False):
        cp = self.cp
        cd = self.cd
        if cp == cd:
            cd = ''
        xflag = self.xflag
        if xflag is None:
            xflag = ''

        x, y, z = self.xyz
        msg = 'MONPNT3 %-8s%s\n' % (self.name, self.label)
        msg += ('        %-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s\n'
                '         %-8s' % (
                    self.axes, self.grid_set, self.elem_set, cp,
                    print_float_8(x), print_float_8(y), print_float_8(z),
                    xflag, cd
                ))
        #card = self.repr_fields()
        return self.comment + msg.rstrip() + '\n'

    def __repr__(self):
        return self.write_card()


class PAERO1(BaseCard):
    """
    Defines associated bodies for the panels in the Doublet-Lattice method.

    +--------+-----+----+----+----+----+----+----+
    |    1   |  2  |  3 |  4 |  5 |  6 |  7 |  8 |
    +========+=====+====+====+====+====+====+====+
    | PAERO1 | PID | B1 | B2 | B3 | B4 | B5 | B6 |
    +--------+-----+----+----+----+----+----+----+
    """
    type = 'PAERO1'
    _field_map = {1: 'pid'}

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the PAERO1 card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : int
            the value for the appropriate field
        """
        return self.Bi[n - 1]

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the PAERO1 card

        Parameters
        ----------
        n : int
            the field number to update
        value : varies
            the value for the appropriate field
        """
        self.Bi[n - 1] = value

    def __init__(self, pid, Bi=None, comment=''):
        """
        Creates a PAERO1 card, which defines associated bodies for the
        panels in the Doublet-Lattice method.

        Parameters
        ----------
        pid : int
            PAERO1 id
        Bi : List[int]; default=None
            CAERO2 ids that are within the same IGID group
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.pid = pid
        if Bi is None:
            Bi = []
        self.Bi = Bi

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAERO1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        Bi = [interpret_value(field) for field in card[2:]]
        Bi2 = []

        for bi in Bi:
            if isinstance(bi, integer_types) and bi >= 0:
                Bi2.append(bi)
            elif bi is not None:
                raise RuntimeError('invalid Bi value on PAERO1 bi=%r' % (bi))
            #else:
                #pass
        return PAERO1(pid, Bi, comment=comment)

    def cross_reference(self, model):
        pass

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass

    def Bodies(self):
        return self.Bi

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
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
    |    1   |  2   |   3    |    4  |   5  |   6  |   7  |   8  |   9  |
    +========+======+========+=======+======+======+======+======+======+
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

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : varies
            the value for the appropriate field
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

        Parameters
        ----------
        n : int
            the field number to update
        value : varies
            the value for the appropriate field
        """
        nnew = n - 8
        spot = nnew // 2
        i = nnew % 2
        if i == 0:
            self.thi[spot] = value
        else:
            self.thn[spot] = value

    def __init__(self, pid, orient, width, AR,
                 thi, thn, lrsb=None, lrib=None, lth1=None, lth2=None,
                 comment=''):
        """
        Creates a PAERO2 card, which defines additional cross-sectional
        properties for the CAERO2 geometry.

        Parameters
        ----------
        pid : int
            PAERO1 id
        orient : str
            Orientation flag. Type of motion allowed for bodies. Refers
            to the aerodynamic coordinate system of ACSID. See AERO entry.
            valid_orientations = {Z, Y, ZY}
        width : float
            Reference half-width of body and the width of the constant
            width interference tube
        AR : float
            Aspect ratio of the interference tube (height/width)
        thi / thn : List[int]
            The first (thi) and last (thn) interference element of a body
            to use the theta1/theta2 array
        lrsb : int; default=None
            int : AEFACT id containing a list of slender body half-widths
                  at the end points of the slender body elements
            None : use width
        lrib : int; default=None
            int : AEFACT id containing a list of interference body
                  half-widths at the end points of the interference elements
            None : use width
        lth1 / lth2 : int; default=None
            AEFACT id for defining theta arrays for interference calculations
            for theta1/theta2
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment

        #: Property identification number. (Integer > 0)
        self.pid = pid

        #: Orientation flag. Type of motion allowed for bodies. Refers to
        #: the aerodynamic coordinate system of ACSID. See AERO entry.
        #: (Character = 'Z', 'Y', or 'ZY')
        self.orient = orient

        #: Reference half-width of body and the width of the constant width
        #: interference tube. (Real > 0.0)
        self.width = width

        #: Aspect ratio of the interference tube (height/width). float>0.
        self.AR = AR

        #: Identification number of an AEFACT entry containing a list of
        #: slender body half-widths at the end points of the
        #: slender body elements. If blank, the value of WIDTH will be used.
        #: (Integer > 0 or blank)
        self.lrsb = lrsb

        #: Identification number of an AEFACT entry containing a list of
        #: slender body half-widths at the end points of the
        #: interference elements. If blank, the value of WIDTH will be used.
        #: (Integer > 0 or blank)
        self.lrib = lrib

        #: Identification number of AEFACT entries for defining theta arrays for
        #: interference calculations. (Integer >= 0)
        self.lth1 = lth1
        self.lth2 = lth2
        self.thi = thi
        self.thn = thn
        if self.lrsb == 0:
            self.lrsb = None
        if self.lrib is 0:
            self.lrib = None


    def validate(self):
        assert self.orient in ['Z', 'Y', 'ZY'], 'orient=%r' % self.orient
        assert isinstance(self.AR, float), 'AR=%r type=%s' % (self.AR, type(self.AR))
        assert isinstance(self.width, float), 'width=%r type=%s' % (self.width, type(self.width))
        assert isinstance(self.thi, list), 'thi=%s type=%s' % (self.thi, type(self.thi))
        assert isinstance(self.thn, list), 'thn=%s type=%s' % (self.thn, type(self.thn))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAERO2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        orient = string(card, 2, 'orient')
        width = double(card, 3, 'width')
        AR = double(card, 4, 'AR')
        lrsb = integer_or_blank(card, 5, 'lrsb')
        lrib = integer_or_blank(card, 6, 'lrib')
        lth1 = integer_or_blank(card, 7, 'lth1')
        lth2 = integer_or_blank(card, 8, 'lth2')
        thi = []
        thn = []
        list_fields = [interpret_value(field) for field in card[9:]]
        nfields = len(list_fields)
        for i in range(9, 9 + nfields, 2):
            thi.append(integer(card, i, 'lth'))
            thn.append(integer(card, i + 1, 'thn'))
        return PAERO2(pid, orient, width, AR, thi, thn,
                      lrsb=lrsb, lrib=lrib, lth1=lth1, lth2=lth2,
                      comment=comment)

    def cross_reference(self, model):
        msg = ' which is required by PAERO2 eid=%s' % self.pid
        if self.lrsb is not None and isinstance(self.lrsb, integer_types):
            self.lrsb_ref = model.AEFact(self.lrsb, msg=msg)
            self.lrsb = self.lrsb_ref
        if self.lrib is not None and isinstance(self.lrib, integer_types):
            self.lrib_ref = model.AEFact(self.lrib, msg=msg)
            self.lrib = self.lrib_ref

    def safe_cross_reference(self, model, debug=False):
        msg = ' which is required by PAERO2 eid=%s' % self.pid
        if self.lrsb is not None and isinstance(self.lrsb, integer_types):
            try:
                self.lrsb_ref = model.AEFact(self.lrsb, msg=msg)
                self.lrsb = self.lrsb_ref
            except KeyError:
                pass
        if self.lrib is not None and isinstance(self.lrib, integer_types):
            try:
                self.lrib_ref = model.AEFact(self.lrib, msg=msg)
                self.lrib = self.lrib_ref
            except KeyError:
                pass

    def uncross_reference(self):
        if self.lrsb is not None and isinstance(self.lrsb, integer_types):
            self.lrsb = self.lrsb_ref.sid # AEFACT id
            del self.lrsb_ref
        if self.lrib is not None and isinstance(self.lrib, integer_types):
            self.lrib = self.lrib_ref.sid # AEFACT id
            del self.lrib_ref

    def Lrsb(self):
        if self.lrsb is None or isinstance(self.lrsb, integer_types):
            return self.lrsb
        return self.lrsb_ref.sid

    def Lrib(self):
        if self.lrib is None or isinstance(self.lrib, integer_types):
            return self.lrib
        return self.lrib_ref.sid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
        """
        list_fields = ['PAERO2', self.pid, self.orient, self.width,
                       self.AR, self.Lrsb(), self.Lrib(), self.lth1, self.lth2]
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

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : varies
            the value for the appropriate field
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

        Parameters
        ----------
        n : int
            the field number to update
        value :varies
            the value for the appropriate field
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

    def __init__(self, pid, nbox, ncontrol_surfaces, x, y, comment=''):
        if comment:
            self.comment = comment
        #: Property identification number. (Integer > 0)
        self.pid = pid
        self.nbox = nbox
        self.ncontrol_surfaces = ncontrol_surfaces
        self.x = x
        self.y = y

    def validate(self):
        assert len(self.x) == len(self.y), 'nx=%s ny=%s' % (len(self.x), len(self.y))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAERO3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        nbox = integer(card, 2, 'nbox')
        ncontrol_surfaces = integer(card, 3, 'ncontrol_surfaces')
        x = []
        y = []
        nfields = card.nfields

        j = 0
        for i in range(6, nfields, 2):
            xi = double(card, i, 'x%i' % j)
            yi = double(card, i + 1, 'y%i' % j)
            x.append(xi)
            y.append(yi)
            j += 1
        return PAERO3(pid, nbox, ncontrol_surfaces, x, y, comment=comment)

    def cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass

    @property
    def npoints(self):
        return len(self.x)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
        """
        list_fields = ['PAERO3', self.pid, self.nbox, self.ncontrol_surfaces, None]
        for (x, y) in zip(self.x, self.y):
            list_fields += [x, y]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class PAERO4(BaseCard):
    """
    Defines properties of each strip element for Strip theory.
    """
    type = 'PAERO4'
    _field_map = {
        1: 'pid', #2:'orient', 3:'width', 4:'AR',
    }

    #def _get_field_helper(self, n):
        #"""
        #Gets complicated parameters on the PAERO3 card

        #Parameters
        #----------
        #n : int
            #the field number to update

        #Returns
        #-------
        #value : varies
            #the value for the appropriate field
        #"""
        #nnew = n - 6
        #if nnew < 0:
            #raise RuntimeError('field n=%i on PAERO3 is invalid' % n)
        #spot = nnew // 2
        #i = nnew % 2
        #if i == 0:
            #value = self.x[spot]
        #else:
            #value = self.y[spot]
        #return value

    #def _update_field_helper(self, n, value):
        #"""
        #Updates complicated parameters on the PAERO3 card

        #Parameters
        #----------
        #n : int
            #the field number to update
        #value :varies
            #the value for the appropriate field
        #"""
        #nnew = n - 6
        #if nnew < 0:
            #raise RuntimeError('field n=%i on PAERO3 is invalid' % n)
        #spot = nnew // 2
        #i = nnew % 2
        #if i == 0:
            #self.x[spot] = value
        #else:
            #self.y[spot] = value

    def __init__(self, pid, docs, caocs, gapocs,
                 cla=0, lcla=0, circ=0, lcirc=0, comment=''):
        if comment:
            self.comment = comment
        #: Property identification number. (Integer > 0)
        self.pid = pid
        self.cla = cla
        self.lcla = lcla
        self.circ = circ
        self.lcirc = lcirc
        self.docs = docs
        self.caocs = caocs
        self.gapocs = gapocs
        assert isinstance(docs, list), docs
        assert isinstance(caocs, list), caocs
        assert isinstance(gapocs, list), gapocs

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PAERO4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        cla = integer_or_blank(card, 2, 'cla', 0)
        lcla = integer_or_blank(card, 3, 'lcla', 0) # ???

        circ = integer_or_blank(card, 4, 'circ', 0)
        lcirc = integer_or_blank(card, 5, 'lcirc', 0) # ???
        nfields = card.nfields

        j = 0
        docs = []
        caocs = []
        gapocs = []
        for i in range(6, nfields, 3):
            doc = double(card, i, 'doc_%i' % j)
            caoc = double(card, i + 1, 'caoc_%i' % j)
            gapoc = double(card, i + 2, 'gapoc_%i' % j)
            docs.append(doc)
            caocs.append(caoc)
            gapocs.append(gapoc)
            j += 1
        return PAERO4(pid, docs, caocs, gapocs,
                      cla=cla, lcla=lcla, circ=circ, lcirc=lcirc, comment=comment)

    def cross_reference(self, model):
        pass

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
        """
        list_fields = ['PAERO4', self.pid, self.cla, self.lcla, self.circ, self.lcirc]
        for doc, caoc, gapoc in zip(self.docs, self.caocs, self.gapocs):
            list_fields += [doc, caoc, gapoc]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class Spline(BaseCard):
    def __init__(self):
        pass


class SPLINE1(Spline):
    """
    Surface Spline Methods
    Defines a surface spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by regular
    arrays of aerodynamic points.

      +---------+-------+-------+------+------+------+----+------+-------+
      |    1    |   2   |    3  |   4  |   5  |   6  |  7 |   8  |   9   |
      +=========+=======+=======+======+======+======+====+======+=======+
      | SPLINE1 | EID   | CAERO | BOX1 | BOX2 | SETG | DZ | METH | USAGE |
      +---------+-------+-------+------+------+------+----+------+-------+
      |         | NELEM | MELEM |      |      |      |    |      |       |
      +---------+-------+-------+------+------+------+----+------+-------+
      | SPLINE1 |   3   |  111  | 115  | 122  |  14  | 0. |      |       |
      +---------+-------+-------+------+------+------+----+------+-------+
    """
    type = 'SPLINE1'
    _field_map = {
        1: 'eid', 2:'caero', 3:'box1', 4:'box2', 5:'setg', 6:'dz',
        7: 'method', 8:'usage', 9:'nelements', 10:'melements',
    }

    def __init__(self, eid, caero, box1, box2, setg, dz=0., method='IPS',
                 usage='BOTH', nelements=10, melements=10, comment=''):
        """
        Creates a SPLINE1, which defines a surface spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the surface spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        method : str; default=IPS
            method for spline fit
            valid_methods = {IPS, TPS, FPS}
            IPS : Harder-Desmarais Infinite Plate Spline
            TPS : Thin Plate Spline
            FPS : Finite Plate Spline
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        nelements : int; default=10
            The number of FE elements along the local spline x-axis if
            using the FPS option
        melements : int; default=10
            The number of FE elements along the local spline y-axis if
            using the FPS option
        comment : str; default=''
            a comment for the card
        """
        Spline.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.caero = caero
        self.box1 = box1
        self.box2 = box2
        self.setg = setg
        self.dz = dz
        self.method = method
        self.usage = usage
        self.nelements = nelements
        self.melements = melements

    def validate(self):
        assert self.nelements > 0, 'nelements = %s' % self.nelements
        assert self.melements > 0, 'melements = %s' % self.melements
        assert self.box2 >= self.box1, 'box1=%s box2=%s' % (self.box1, self.box2)
        assert self.method in ['IPS', 'TPS', 'FPS'], 'method = %s' % self.method
        assert self.usage in ['FORCE', 'DISP', 'BOTH'], 'usage = %s' % self.usage

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        box1 = integer(card, 3, 'box1')
        box2 = integer(card, 4, 'box2')
        setg = integer(card, 5, 'setg')
        dz = double_or_blank(card, 6, 'dz', 0.0)
        method = string_or_blank(card, 7, 'method', 'IPS')
        usage = string_or_blank(card, 8, 'usage', 'BOTH')
        nelements = integer_or_blank(card, 9, 'nelements', 10)
        melements = integer_or_blank(card, 10, 'melements', 10)
        assert len(card) <= 11, 'len(SPLINE1 card) = %i\ncard=%s' % (len(card), card)
        return SPLINE1(eid, caero, box1, box2, setg, dz, method, usage,
                       nelements, melements, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        caero = data[1]
        box1 = data[2]
        box2 = data[3]
        setg = data[4]
        dz = data[5]
        method = data[6]
        usage = data[7]
        nelements = data[8]
        melements = data[9]
        assert len(data) == 10, 'data = %s' % data
        return SPLINE1(eid, caero, box1, box2, setg, dz, method, usage,
                       nelements, melements, comment=comment)

    @property
    def aero_element_ids(self):
        return np.arange(self.box1, self.box2 + 1)

    def CAero(self):
        if isinstance(self.caero, integer_types):
            return self.caero
        return self.caero_ref.eid

    def Set(self):
        if isinstance(self.setg, integer_types):
            return self.setg
        return self.setg_ref.sid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by SPLINE1 eid=%s' % self.eid
        self.caero = model.CAero(self.CAero(), msg=msg)
        self.caero_ref = self.caero
        self.setg = model.Set(self.Set(), msg=msg)
        self.setg.cross_reference(model, 'Node')
        self.setg_ref = self.setg

        nnodes = len(self.setg_ref.ids)
        if nnodes < 3:
            msg = 'SPLINE1 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
            msg += str(self)
            msg += str(self.setg_ref)
            raise RuntimeError(msg)

    def safe_cross_reference(self, model):
        msg = ' which is required by SPLINE1 eid=%s' % self.eid
        try:
            self.caero = model.CAero(self.CAero(), msg=msg)
            self.caero_ref = self.caero
        except KeyError:
            pass

        try:
            self.setg = model.Set(self.Set(), msg=msg)
            self.setg.cross_reference(model, 'Node')
            self.setg_ref = self.setg

            nnodes = len(self.setg_ref.ids)
            if nnodes < 3:
                msg = 'SPLINE1 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
                msg += str(self)
                msg += str(self.setg_ref)
                raise RuntimeError(msg)
        except KeyError:
            pass

    def uncross_reference(self):
        self.caero = self.CAero()
        self.setg = self.Set()
        del self.caero_ref, self.setg_ref

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
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
    Defines a beam spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by regular
    arrays of aerodynamic points.

      +---------+------+-------+-------+-------+------+----+------+-----+
      |    1    |   2  |   3   |   4   |   5   |  6   |  7 |   8  |  9  |
      +=========+======+=======+=======+=======+======+====+======+=====+
      | SPLINE2 | EID  | CAERO |  ID1  |  ID2  | SETG | DZ | DTOR | CID |
      +---------+------+-------+-------+-------+------+----+------+-----+
      |         | DTHX | DTHY  | None  | USAGE |      |    |      |     |
      +---------+------+-------+-------+-------+------+----+------+-----+

      +---------+------+-------+-------+-------+------+----+------+-----+
      |    1    |   2  |   3   |   4   |   5   |  6   |  7 |   8  |  9  |
      +=========+======+=======+=======+=======+======+====+======+=====+
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

    def __init__(self, eid, caero, id1, id2, setg, dz=0.0, dtor=1.0, cid=0,
                 dthx=None, dthy=None, usage='BOTH', comment=''):
        """
        Creates a SPLINE2 card, which defines a beam spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        id1 / id2 : int
            First/last box/body id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the beam spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        dtor : float; default=1.0
            Torsional flexibility ratio (EI/GJ).
            Use 1.0 for bodies (CAERO2).
        cid : int; default=0
            Rectangular coordinate system for which the y-axis defines the
            axis of the spline. Not used for bodies, CAERO2
        dthx : float; default=None
            Rotational attachment flexibility.
            DTHX : Used for rotation about the spline's x-axis (in-plane
                   bending rotations).  It is not used for bodies (CAERO2).
            DTHY : Used for rotation about the spline's y-axis (torsion).
                   It is used for slope of bodies.
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        comment : str; default=''
            a comment for the card
        """
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.caero = caero
        self.id1 = id1
        self.id2 = id2
        self.setg = setg
        self.dz = dz
        self.dtor = dtor
        self.cid = cid
        self.dthx = dthx
        self.dthy = dthy
        self.usage = usage
        assert self.id2 >= self.id1, 'id2=%s id1=%s' % (self.id2, self.id1)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        id1 = integer(card, 3, 'id1')
        id2 = integer(card, 4, 'id2')
        setg = integer(card, 5, 'setg')
        dz = double_or_blank(card, 6, 'dz', 0.0)
        dtor = double_or_blank(card, 7, 'dtor', 1.0)
        cid = integer_or_blank(card, 8, 'cid', 0)
        dthx = double_or_blank(card, 9, 'dthx')
        dthy = double_or_blank(card, 10, 'dthy')

        usage = string_or_blank(card, 12, 'usage', 'BOTH')
        assert len(card) <= 13, 'len(SPLINE2 card) = %i\ncard=%s' % (len(card), card)
        return SPLINE2(eid, caero, id1, id2, setg, dz, dtor, cid,
                       dthx, dthy, usage, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by SPLINE2 eid=%s' % self.eid
        self.cid = model.Coord(self.Cid(), msg=msg)
        self.cid_ref = self.cid
        self.caero = model.CAero(self.CAero(), msg=msg)
        self.caero_ref = self.caero
        self.setg = model.Set(self.Set(), msg=msg)
        self.setg_ref = self.setg
        self.setg.cross_reference(model, 'Node', msg=msg)

        nnodes = len(self.setg_ref.ids)
        if nnodes < 2:
            msg = 'SPLINE2 requires at least 2 nodes; nnodes=%s\n' % (nnodes)
            msg += str(self)
            msg += str(self.setg_ref)
            raise RuntimeError(msg)

    def safe_cross_reference(self, model):
        msg = ' which is required by SPLINE2 eid=%s' % self.eid
        try:
            self.cid = model.Coord(self.Cid(), msg=msg)
            self.cid_ref = self.cid
        except KeyError:
            pass

        try:
            self.caero = model.CAero(self.CAero(), msg=msg)
            self.caero_ref = self.caero
        except KeyError:
            pass

        try:
            self.setg = model.Set(self.Set(), msg=msg)
            self.setg_ref = self.setg
            self.setg.cross_reference(model, 'Node', msg=msg)

            nnodes = len(self.setg_ref.ids)
            if nnodes < 2:
                msg = 'SPLINE2 requires at least 2 nodes; nnodes=%s\n' % (nnodes)
                msg += str(self)
                msg += str(self.setg_ref)
                raise RuntimeError(msg)
        except KeyError:
            pass

    def uncross_reference(self):
        self.cid = self.Cid()
        self.caero = self.CAero()
        self.setg = self.Set()
        del self.cid_ref, self.caero_ref, self.setg_ref

    @property
    def aero_element_ids(self):
        return np.arange(self.id1, self.id2 + 1)

    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def CAero(self):
        if isinstance(self.caero, integer_types):
            return self.caero
        return self.caero_ref.eid

    def Set(self):
        if isinstance(self.setg, integer_types):
            return self.setg
        return self.setg_ref.sid

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
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
    """
    Defines a constraint equation for aeroelastic problems.
    Useful for control surface constraints.
    """
    type = 'SPLINE3'
    _field_map = {
        1: 'eid', 2:'caero', 3:'box_id',
        7: 'a1', 8:'usage',
    }
        #5:'g1', 6:'c1',
        #9:G2,C2,A2...

    def __init__(self, eid, caero, box_id, components,
                 nids, displacement_components, coeffs,
                 usage='BOTH', comment=''):
        Spline.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.caero = caero
        self.box_id = box_id
        self.components = components
        self.usage = usage
        self.nids = nids
        self.displacement_components = displacement_components
        self.coeffs = coeffs

    def validate(self):
        is_failed = False
        msg = ''
        if self.components not in [0, 1, 2, 3, 4, 5, 6]:
            msg += 'components=%s must be [0, 1, 2, 3, 4, 5, 6]\n' % (
                self.components)
            is_failed = True

        for i, disp_component, coeff  in zip(count(), self.displacement_components, self.coeffs):
            if disp_component not in [0, 1, 2, 3, 4, 5, 6]:
                msg += 'i=%s displacement_component=%s must be [0, 1, 2, 3, 4, 5, 6]\n' % (
                    i, disp_component)
                is_failed = True

        if self.usage not in ['FORCE', 'DISP', 'BOTH']:
            msg += 'usage=%r must be in [FORCE, DISP, BOTH]' % self.usage
            is_failed = True

        if is_failed:
            msg += str(self)
            raise RuntimeError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        box_id = integer(card, 3, 'box_id')
        components = integer(card, 4, 'comp')
        g1 = integer(card, 5, 'G1')
        c1 = integer(card, 5, 'C1')
        a1 = double(card, 5, 'A1')
        usage = string_or_blank(card, 8, 'usage', 'BOTH')

        nfields = len(card) - 1
        nrows = nfields // 8
        if nfields % 8:
            nrows += 1

        i = 1
        Gi = [g1]
        ci = [c1]
        ai = [a1]
        for irow in range(1, nrows):
            j = 1 + nrows * 8
            gii = integer(card, j, 'Gi_' % i)
            cii = integer(card, j + 1, 'Ci_' % i)
            aii = double(card, j + 2, 'Ai_' % i)
            Gi.append(gii)
            ci.append(cii)
            ai.append(aii)
            if card[j + 3] or card[j + 4] or card[j + 5]:
                i += 1
                gii = integer(card, j, 'Gi_' % i)
                cii = integer(card, j + 1, 'Ci_' % i)
                aii = double(card, j + 2, 'Ai_' % i)
                Gi.append(gii)
                ci.append(cii)
                ai.append(aii)
            i += 1
        return SPLINE3(eid, caero, box_id, components, Gi, ci, ai, usage,
                       comment=comment)

    def CAero(self):
        if isinstance(self.caero, integer_types):
            return self.caero
        return self.caero_ref.eid

    def Set(self):
        if isinstance(self.setg, integer_types):
            return self.setg
        return self.setg_ref.sid

    def cross_reference(self, model):
        msg = ' which is required by SPLINE3 eid=%s' % self.eid
        self.caero = model.CAero(self.CAero(), msg=msg)
        self.caero_ref = self.caero
        self.setg = model.Set(self.Set(), msg=msg)
        self.setg.cross_reference(model, 'Node')
        self.setg_ref = self.setg

        nnodes = len(self.setg_ref.ids)
        if nnodes < 3:
            msg = 'SPLINE3 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
            msg += str(self)
            msg += str(self.setg_ref)
            raise RuntimeError(msg)

    def uncross_reference(self):
        self.caero = self.CAero()
        self.setg = self.Set()
        del self.caero_ref, self.setg_ref

    def raw_fields(self):
        list_fields = [
            'SPLINE3', self.eid, self.caero, self.box_id,
            self.nids[0], self.displacement_components[0], self.coeffs[0], self.usage]
        for nid, disp_c, coeff in zip(self.nids[1:], self.displacement_components[1:],
                                      self.coeffs[1:]):
            list_fields += [nid, disp_c, coeff, None]
        return list_fields


class SPLINE4(Spline):
    """
    Surface Spline Methods
    Defines a curved surface spline for interpolating motion and/or forces for
    aeroelastic problems on general aerodynamic geometries using either the
    Infinite Plate, Thin Plate or Finite Plate splining method.

    +---------+-------+-------+--------+-----+------+----+------+-------+
    |    1    |   2   |   3   |    4   |  5  |   6  |  7 |   8  |   9   |
    +=========+=======+=======+========+=====+======+====+======+=======+
    | SPLINE4 |  EID  | CAERO | AELIST | --- | SETG | DZ | METH | USAGE |
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

    def __init__(self, eid, caero, aelist, setg, dz, method, usage,
                 nelements, melements, comment=''):
        Spline.__init__(self)
        if comment:
            self.comment = comment
        self.eid = eid
        self.caero = caero
        self.aelist = aelist
        self.setg = setg
        self.dz = dz
        self.method = method
        self.usage = usage
        self.nelements = nelements
        self.melements = melements

    def validate(self):
        assert self.method in ['IPS', 'TPS', 'FPS'], 'method = %s' % self.method
        assert self.usage in ['FORCE', 'DISP', 'BOTH'], 'uasge = %s' % self.usage

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        aelist = integer(card, 3, 'aelist')
        # None
        setg = integer(card, 5, 'setg')
        dz = double_or_blank(card, 6, 'dz', 0.0)
        method = string_or_blank(card, 7, 'method', 'IPS')
        usage = string_or_blank(card, 8, 'usage', 'BOTH')
        nelements = integer_or_blank(card, 9, 'nelements', 10)
        melements = integer_or_blank(card, 10, 'melements', 10)
        assert len(card) <= 11, 'len(SPLINE4 card) = %i\ncard=%s' % (len(card), card)
        return SPLINE4(eid, caero, aelist, setg, dz, method, usage,
                       nelements, melements, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        eid = data[0]
        caero = data[1]
        aelist = data[2]
        setg = data[3]
        dz = data[4]
        method = data[5]
        usage = data[6]
        nelements = data[7]
        melements = data[8]
        assert len(data) == 9, 'data = %s' % (data)
        return SPLINE4(eid, caero, aelist, setg, dz, method, usage,
                       nelements, melements, comment=comment)

    def CAero(self):
        if isinstance(self.caero, integer_types):
            return self.caero
        return self.caero_ref.eid

    def AEList(self):
        if isinstance(self.aelist, integer_types):
            return self.aelist
        return self.aelist_ref.sid

    def Set(self):
        if isinstance(self.setg, integer_types):
            return self.setg
        return self.setg_ref.sid

    @property
    def aero_element_ids(self):
        return self.aelist_ref.elements

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by SPLINE4 eid=%s' % self.eid
        self.caero = model.CAero(self.CAero(), msg=msg)
        self.setg = model.Set(self.Set(), msg=msg)
        self.aelist = model.AEList(self.aelist, msg=msg)
        self.setg.cross_reference(model, 'Node')

        self.caero_ref = self.caero
        self.setg_ref = self.setg
        self.aelist_ref = self.aelist

        nnodes = len(self.setg_ref.ids)
        if nnodes < 3:
            msg = 'SPLINE4 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
            msg += str(self)
            msg += str(self.setg_ref)
            raise RuntimeError(msg)

    def uncross_reference(self):
        self.caero = self.CAero()
        self.setg = self.Set()
        self.aelist = self.AEList()
        del  self.caero_ref, self.setg_ref, self.aelist_ref

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
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

    +=========+======+=======+========+=======+======+====+=======+=======+
    |    1    |  2   |    3  |    4   |   5   |   6  |  7 |   8   |   9   |
    +=========+======+=======+========+=======+======+====+=======+=======+
    | SPLINE5 | EID  | CAERO | AELIST |  ---  | SETG | DZ | DTOR  |  CID  |
    +---------+------+-------+--------+-------+------+----+-------+-------+
    |         | DTHX | DTHY  |  ---   | USAGE | METH |    | FTYPE | RCORE |
    +---------+------+-------+--------+-------+------+----+-------+-------+

    METH, FTYPE, RCORE are in 2012+ (not MSC.2005r2 or NX.10)
    """
    type = 'SPLINE5'
    _field_map = {
        1: 'eid', 2:'caero', 3:'aelist', 5:'setg', 6:'dz',
        7: 'dtor', 8:'cid', 9:'thx', 10:'thy', 12:'usage',
        13 : 'meth', 15 : 'ftype', 16 : 'rcore',
    }

    def __init__(self, eid, caero, aelist, setg, dz, dtor, cid,
                 thx, thy, usage, method, ftype, rcore, comment=''):
        Spline.__init__(self)
        if comment:
            self.comment = comment

        self.eid = eid
        self.caero = caero
        self.aelist = aelist
        self.setg = setg
        self.dz = dz
        self.dtor = dtor
        self.cid = cid
        self.thx = thx
        self.thy = thy
        self.usage = usage
        self.meth = method
        self.ftype = ftype
        self.rcore = rcore
        assert method in ['BEAM', 'RIS'], method
        assert ftype in ['WF0', 'WF2'], ftype

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a SPLINE5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        aelist = integer(card, 3, 'aelist')
        # None
        setg = integer(card, 5, 'setq')
        dz = double_or_blank(card, 6, 'dz', 0.0)
        dtor = double_or_blank(card, 7, 'dtor', 1.0)
        cid = integer_or_blank(card, 8, 'cid', 0)
        thx = double(card, 9, 'thx')
        thy = double(card, 10, 'thy')
        usage = string_or_blank(card, 12, 'usage', 'BOTH')
        # per nast/tpl/fmondsp.dat, METH can be a double(0.0) ???
        method = string_or_blank(card, 13, 'meth', 'BEAM')
        ftype = string_or_blank(card, 15, 'ftype', 'WF2')
        rcore = double_or_blank(card, 16, 'rcore')

        usage = string_or_blank(card, 12, 'usage', 'BOTH')
        assert len(card) <= 16, 'len(SPLINE5 card) = %i\n%s' % (len(card), card)
        return SPLINE5(eid, caero, aelist, setg, dz, dtor, cid, thx, thy,
                       usage, method, ftype, rcore,
                       comment=comment)
    @property
    def aero_element_ids(self):
        return self.aelist_ref.elements


    def Cid(self):
        if isinstance(self.cid, integer_types):
            return self.cid
        return self.cid_ref.cid

    def CAero(self):
        if isinstance(self.caero, integer_types):
            return self.caero
        return self.caero_ref.eid

    def AEList(self):
        if isinstance(self.aelist, integer_types):
            return self.aelist
        return self.aelist_ref.sid

    def Set(self):
        if isinstance(self.setg, integer_types):
            return self.setg
        return self.setg_ref.sid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by SPLINE5 eid=%s' % self.eid
        self.cid = model.Coord(self.Cid(), msg=msg)
        self.cid_ref = self.cid
        self.caero = model.CAero(self.CAero(), msg=msg)
        self.caero_ref = self.caero
        self.setg = model.Set(self.Set(), msg=msg)
        self.setg_ref = self.setg
        self.setg.cross_reference(model, 'Node')
        self.aelist = model.AEList(self.AEList(), msg=msg)
        self.aelist_ref = self.aelist

        nnodes = len(self.setg_ref.ids)
        if nnodes < 3:
            msg = 'SPLINE5 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
            msg += str(self)
            msg += str(self.setg_ref)
            raise RuntimeError(msg)

    def safe_cross_reference(self, model):
        msg = ' which is required by SPLINE5 eid=%s' % self.eid
        try:
            self.cid = model.Coord(self.Cid(), msg=msg)
            self.cid_ref = self.cid
        except KeyError:
            pass
        try:
            self.caero = model.CAero(self.CAero(), msg=msg)
            self.caero_ref = self.caero
        except KeyError:
            pass

        try:
            self.setg = model.Set(self.Set(), msg=msg)
            self.setg_ref = self.setg
            nnodes = len(self.setg_ref.ids)
            if nnodes < 3:
                msg = 'SPLINE5 requires at least 3 nodes; nnodes=%s\n' % (nnodes)
                msg += str(self)
                msg += str(self.setg_ref)
                raise RuntimeError(msg)
        except KeyError:
            pass

        try:
            self.setg.cross_reference(model, 'Node')
            self.aelist = model.AEList(self.AEList(), msg=msg)
            self.aelist_ref = self.aelist
        except KeyError:
            pass

    def uncross_reference(self):
        self.cid = self.Cid()
        self.caero = self.CAero()
        self.setg = self.Set()
        self.aelist = self.AEList()
        del  self.cid, self.caero_ref, self.setg_ref, self.aelist_ref

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card
        """
        list_fields = ['SPLINE5', self.eid, self.CAero(), self.AEList(), None,
                       self.Set(), self.dz, self.dtor, self.Cid(), self.thx,
                       self.thy, None, self.usage, None, self.ftype, self.rcore]
        return list_fields

    def repr_fields(self):
        dz = set_blank_if_default(self.dz, 0.)
        usage = set_blank_if_default(self.usage, 'BOTH')
        list_fields = ['SPLINE5', self.eid, self.CAero(), self.AEList(), None,
                       self.Set(), dz, self.dtor, self.Cid(), self.thx, self.thy,
                       None, usage, self.meth, None, self.ftype, self.rcore]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class TRIM(BaseCard):
    """
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    |   1  |   2    |   3  |    4   |    5   |  6  |    7   |  8  |     9    |
    +======+========+======+========+========+=====+========+=====+==========+
    | TRIM |   ID   | MACH |    Q   | LABEL1 | UX1 | LABEL2 | UX2 | IS_RIGID |
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    |      | LABEL3 |  UX3 | LABEL4 |   UX4  | ... |        |     |          |
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    """
    type = 'TRIM'
    _field_map = {
        1: 'sid', 2:'mach', 3:'q', 8:'aeqr',
    }

    def _get_field_helper(self, n):
        """
        Gets complicated parameters on the TRIM card

        Parameters
        ----------
        n : int
            the field number to update

        Returns
        -------
        value : varies
            the value for the appropriate field
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
        raise KeyError('Field %r=%r is an invalid TRIM entry.' % (n, value))

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the TRIM card

        Parameters
        ----------
        n : int
            the field number to update
        value : varies
            the value for the appropriate field
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
        raise KeyError('Field %r=%r is an invalid TRIM entry.' % (n, value))

    def __init__(self, sid, mach, q, labels, uxs, aeqr=0.0, comment=''):
        """
        Creates a TRIM card for a static aero (144) analysis.

        Parameters
        ----------
        sid : int
            the trim id; referenced by the Case Control TRIM field
        mach : float
            the mach number
        q : float
            dynamic pressure
        labels : List[str]
            names of the fixed variables
        uxs : List[float]
            values corresponding to labels
        aeqr : float
            0.0 : rigid trim analysis
            1.0 : elastic trim analysis
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        #: Trim set identification number. (Integer > 0)
        self.sid = sid
        #: Mach number. (Real > 0.0 and != 1.0)
        self.mach = mach
        #: Dynamic pressure. (Real > 0.0)
        self.q = q

        #: The label identifying aerodynamic trim variables defined on an
        #: AESTAT or AESURF entry.
        self.labels = labels

        #: The magnitude of the aerodynamic extra point degree-of-freedom.
        #: (Real)
        self.uxs = uxs

        #: Flag to request a rigid trim analysis (Real > 0.0 and < 1.0;
        #: Default = 1.0. A value of 0.0 provides a rigid trim analysis,
        #: not supported
        self.aeqr = aeqr

    def validate(self):
        assert self.mach >= 0.0, 'mach = %r' % self.mach
        assert self.mach != 1.0, 'mach = %r' % self.mach
        assert self.q > 0.0, 'q=%s' % self.q
        if len(set(self.labels)) != len(self.labels):
            msg = 'not all labels are unique; labels=%s' % str(self.labels)
            raise RuntimeError(msg)
        if len(self.labels) != len(self.uxs):
            msg = 'nlabels=%s != nux=%s; labels=%s uxs=%s' % (
                len(self.labels), len(self.uxs), str(self.labels), str(self.uxs))
            raise RuntimeError(msg)

    def _verify(self, suport, suport1, aestats, aeparms, aelinks, aesurf, xref=True):
        """
        Magic function that makes TRIM cards not frustrating.

        Warning
        -------
        TODO: This probably gets AELINKs/AEPARMs/AESURFSs wrong.

        The TRIM equality
        -----------------
        ndelta = (naestat + naesurf + naeparm) - (
               - (ntrim + ntrim_aesurf? + naelink + nsuport_dofs + nsuport1_dofs)
        ndelta = 0
        ntrim_aesurf is not included, but it might exist...

        Steps to a TRIM analysis
        ------------------------
        1.  Define the number of independent control surfaces (naesurf)
            Make an AESURF for each.  Dual link the AESURFs if you can
            to avoid needing an AELINK (e.g., +roll is left aileron down,
            right aileron up).
            Horizontal Tail : name it DPITCH
            Vertical Tail   : name it DYAW
            Aileron         : name it DROLL
        2.  Create AELINKs if necessary.
        3.  Add the AESTAT variables.  Include one for each DOF the
            aircraft can move in the frame of the model
            (e.g., half/full model).
                Half model (2.5g pullup, abrupt pitch):
                  - 2d pitch/plunge, 1 control : URDD3, URDD5, PITCH, ANGLEA
                Full model (2.5g pullup, abrupt pitch):
                  - 3d pitch/plunge, 3 control : URDD3, URDD5, PITCH, ANGLEA, YAW (???)
        4.  Add the TRIM card to lock the variables that could theoretically move
            in the plane of the analysis that are known.
                Half model:
                   2.5g pullup   : lock URDD3=2.5, URDD5=0, PITCH=0
                                   solve for ANGLEA, DPITCH
                                   use DPITCH
                   abrupt pitch  : lock URDD3=1.0, URDD5=0, ANGLEA=5
                                   solve for PITCH, DPITCH
                                   use DPITCH
                Full model:
                   2.5g pullup   : lock URDD3=2.5, URDD4=0, URDD5=0,  PITCH=0, YAW=0,
                                   lock SIDES=0,  ROLL=0
                                   solve for ANGLEA, DPITCH
                                   use DPITCH, DYAW, DROLL
                                   TODO: probably wrong
                   30 degree yaw : lock URDD3=1.0, URDD4=0, ANGLEA=5, PITCH=0, YAW=30,
                                   lock DPITCH=0, ROLL=0
                                   solve for SIDES, URDD5
                                   use DPITCH, DYAW, DROLL
                                   TODO: probably wrong

        5.  Note that we could have simplified our full model AESTAT/TRIM
            cards (they can be the same as for a half model), but we'd
            like to be able to do multiple load cases in the same deck.

        6.  Add some SUPORT/SUPORT1 DOFs to ignore non-relevant motion in
            certain DOFs (e.g., z-motion).  Add enough to satisfy the TRIM
            equality.

        Doesn't Consider
        ----------------
         - AELINK
         - AEPARM
         - AESURFS

        +------------------------------------------------+
        |                 Default AESTATs                |
        +--------+---------+-----------------------------+
        | ANGLEA | ur (R2) | Angle of Attack             |
        | YAW    | ur (R3) | Yaw Rate                    |
        | SIDES  | ur (R3) | Angle of Sideslip           |
        +--------+---------+-----------------------------+
        | ROLL   | ůr (R1) | Roll Rate                   |
        | PITCH  | ůr (R2) | Pitch Rate                  |
        +--------+---------+-----------------------------+
        | URDD1  | ür (T1) | Longitudinal (See Remark 3) |
        | URDD2  | ür (T2) | Lateral                     |
        | URDD3  | ür (T3) | Vertical                    |
        | URDD4  | ür (R1) | Roll                        |
        | URDD5  | ür (R2) | Pitch                       |
        | URDD6  | ür (R3) | Yaw                         |
        +--------+---------+-----------------------------+
        """
        if xref:
            nsuport_dofs = 0
            nsuport1_dofs = 0
            suport_dofs = set()
            assert isinstance(suport, list), type(suport)
            for suporti in suport:
                for nid, cs in zip(suporti.node_ids, suporti.Cs):
                    for ci in cs:
                        #print('  SUPORT: nid=%r C=%r' % (nid, ci))
                        dof = (nid, ci)
                        if dof in suport_dofs:
                            msg = 'Duplicate DOF\n  dof=%s suport_dofs=%s' % (
                                str(dof), str(suport_dofs))
                            raise RuntimeError(msg)
                        suport_dofs.add(dof)
                        nsuport_dofs += 1

            if suport1:
                conid = suport1.conid
                nids = suport1.node_ids
                for nid, cs in zip(nids, suport1.Cs):
                    for ci in cs:
                        #print('  SUPORT1: id=%r nid=%r C=%r' % (conid, nid, ci))
                        dof = (nid, ci)
                        if dof in suport_dofs:
                            msg = 'dof=%s suport_dofs=%s' % (str(dof), str(suport_dofs))
                            raise RuntimeError(msg)
                        suport_dofs.add(dof)
                        nsuport1_dofs += 1

            aesurf_names = [aesurfi.label for aesurfi in aesurf.values()]
            aestat_labels = [aestat.label for aestat in aestats.values()]
            aeparm_labels = [aeparm.label for aeparm in aeparms.values()]
            naestat = len(aestat_labels)
            ntrim = len(self.labels)
            naesurf = len(aesurf_names)
            naeparm = len(aeparm_labels)
            naelink = 0
            if self.sid in aelinks:
                naelink = len(aelinks[self.sid])
            if 0 in aelinks:
                #  TODO: what is this...is 0 the global subcase?
                naelink += len(aelinks[0])

            ntrim_aesurf = 0
            labels = aestat_labels + aesurf_names + aeparm_labels
            for label in self.labels:
                if label not in labels:
                    msg = 'label=%r\n aestat_labels=%s\n aeparm_labels=%s\n aesurf_names=%s' % (
                        label, aestat_labels, aeparm_labels, aesurf_names)
                    raise RuntimeError(msg)

                if label in aesurf_names:
                    #print('AESTAT/AESURF label = %r' % label)
                    ntrim_aesurf += 1

            # TODO: this doesn't work for multiple subcases
            #ntotal_suport_dofs = nsuport_dofs, nsuport1_dofs
            #ndelta = ntrim - nsuport_dofs - nsuport1_dofs - naesurf
            #if ndelta != 0:
                #msg = 'ntrim - nsuport_dofs - nsuport1_dofs - naesurf = ndelta = %s; ndelta != 0\n' % ndelta
                #msg += 'ntrim=%s nsuport_dofs=%s nsuport1_dofs=%s naesurfs=%s' % (
                    #ntrim, nsuport_dofs, nsuport1_dofs, naesurf)
                #raise RuntimeError(msg)

           #ndelta = (naestat + naesurf + naeparm + ntrim_aesurf) - (ntrim + naelink + nsuport_dofs + nsuport1_dofs)
           #if ndelta != 0:
               #msg = '(naestat + naesurf + naeparm + ntrim_aesurf) - (ntrim + naelink + nsuport_dofs + nsuport1_dofs) = ndelta = %s; ndelta != 0\n' % ndelta
               #msg += ('naestat=%s naesurf=%s naeparm=%s ntrim_aesurfs=%s\n'
                       #'ntrim=%s naelink=%s nsuport_dofs=%s nsuport1_dofs=%s' % (
                           #naestat, naesurf, naeparms, ntrim_aesurf,
                           #ntrim, naelink, nsuport_dofs, nsuport1_dofs))

            ndelta = (naestat + naesurf + naeparm) - (ntrim + naelink + nsuport_dofs + nsuport1_dofs) #+ ntrim_aesurfs
            if ndelta != 0:
                msg = '(naestat + naesurf + naeparm) - (ntrim + ntrim_aesurf? + naelink + nsuport_dofs + nsuport1_dofs) = ndelta = %s; ndelta != 0\n' % ndelta
                msg += 'naestat=%s naesurf=%s naeparm=%s ntrim=%s ntrim_aesurf=%s naelink=%s nsuport_dofs=%s nsuport1_dofs=%s' % (
                    naestat, naesurf, naeparm, ntrim, ntrim_aesurf, naelink, nsuport_dofs, nsuport1_dofs)
                raise RuntimeError(msg)

    def cross_reference(self, model):
        pass
        #self.suport = model.suport
        #self.suport1 = model.suport1
        #self.aestats = model.aestats
        #self.aelinks = model.aelinks
        #self.aesurf = model.aesurf

    def safe_cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        mach = double(card, 2, 'mach')
        q = double(card, 3, 'q')
        labels = []
        uxs = []

        label = string_or_blank(card, 4, 'label1')
        if label:
            ux = double(card, 5, 'ux1')
            uxs.append(ux)
            labels.append(label)

        label = string_or_blank(card, 6, 'label2')
        if label:
            ux = double(card, 7, 'ux1')
            uxs.append(ux)
            labels.append(label)
        aeqr = double_or_blank(card, 8, 'aeqr', 0.0)

        i = 9
        n = 3
        while i < len(card):
            label = string(card, i, 'label%i' % n)
            ux = double(card, i + 1, 'ux%i' % n)
            labels.append(label)
            uxs.append(ux)
            i += 2
        return TRIM(sid, mach, q, labels, uxs, aeqr, comment=comment)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['TRIM', self.sid, self.mach, self.q]
        nlabels = len(self.labels)
        assert nlabels > 0, self.labels
        for (i, label, ux) in zip(count(), self.labels, self.uxs):
            list_fields += [label, ux]
            if i == 1:
                list_fields += [self.aeqr]
        if nlabels == 1:
            list_fields += [None, None, self.aeqr]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)
