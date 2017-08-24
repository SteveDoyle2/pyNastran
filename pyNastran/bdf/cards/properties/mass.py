# pylint: disable=C0103,R0902,R0904,R0914
"""
All mass properties are defined in this file.  This includes:

 * NSM
 * PMASS

All mass properties are PointProperty and Property objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import integer_types
from pyNastran.bdf.cards.base_card import _node_ids, expand_thru
from pyNastran.bdf.cards.base_card import BaseCard, Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, double_or_blank, string)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class NSM1(Property):
    """
    Defines a set of non structural mass.

    +------+-----+------+-------+-----+----+----+----+----+
    |  1   |  2  |  3   |   4   |  5  | 6  | 7  | 8  | 9  |
    +======+=====+======+=======+=====+====+====+====+====+
    | NSM1 | SID | TYPE | VALUE | ID  | ID | ID | ID | ID |
    +------+-----+------+-------+-----+----+----+----+----+
    |      |  ID |  ID  |  ID   | etc |    |    |    |    |
    +------+-----+------+-------+-----+----+----+----+----+
    """
    type = 'NSM1'
    valid_properties = [
        'PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP',
        'PROD', 'CONROD', 'PBEND', 'PSHEAR', 'PTUBE', 'PCONEAX', 'PRAC2D',
        'ELEMENT',
    ]

    def __init__(self, sid, nsm_type, value, ids, comment=''):
        """
        Creates an NSM1 card

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        value : float
            the non-structural pass per unit length/area
        ids : List[int]
            property ids or element ids depending on nsm_type
        comment : str; default=''
            a comment for the card
        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        if isinstance(ids, integer_types):
            ids = [ids]
        self.sid = sid
        self.nsm_type = nsm_type
        self.ids = ids
        self.value = value
        if self.nsm_type not in self.valid_properties:
            raise TypeError('nsm_type=%r must be in [%s]' % (
                self.nsm_type, ', '.join(self.valid_properties)))
        assert isinstance(ids, list), 'ids=%r is not a list' % (ids)

    @property
    def Type(self):
        """gets the nsm_type"""
        return self.nsm_type
    @Type.setter
    def Type(self, nsm_type):
        """sets the nsm_type"""
        self.nsm_type = nsm_type

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a NSM1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        nsm_type = string(card, 2, 'Type')
        value = double(card, 3, 'value')
        ids = card[4:]
        return NSM1(sid, nsm_type, value, ids, comment=comment)

    def cross_reference(self, model):
        pass

    def raw_fields(self):
        #nodes = self.node_ids
        list_fields = ['NSM1', self.sid, self.nsm_type, self.value] + self.ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class NSM(Property):
    """
    Defines a set of non structural mass.

    +-----+-----+------+----+-------+----+-------+----+-------+
    |  1  |  2  |  3   |  4 |   5   | 6  |   7   | 8  |   9   |
    +=====+=====+======+====+=======+====+=======+====+=======+
    | NSM | SID | TYPE | ID | VALUE | ID | VALUE | ID | VALUE |
    +-----+-----+------+----+-------+----+-------+----+-------+
    """
    type = 'NSM'
    _field_map = {
        1: 'sid', 2:'Type', 3:'id', 4:'value'
    }

    #: Set points to either Property entries or Element entries.
    #: Properties are:
    valid_properties = [
        'PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP',
        'PROD', 'CONROD', 'PBEND', 'PSHEAR', 'PTUBE', 'PCONEAX', 'PRAC2D',
        'ELEMENT',
    ]

    def __init__(self, sid, nsm_type, id, value, comment=''):
        """
        Creates an NSM card

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        id : int
            property id or element id depending on nsm_type
        value : float
            the non-structural pass per unit length/area
        comment : str; default=''
            a comment for the card
        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.nsm_type = nsm_type
        self.id = id
        self.value = value
        if self.nsm_type not in self.valid_properties:
            raise TypeError('nsm_type=%r must be in [%s]' % (
                self.nsm_type, ', '.join(self.valid_properties)))

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds an NSM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            the index of the card that's being parsed
        comment : str; default=''
            a comment for the card
        """
        noffset = 2 * icard
        sid = integer(card, 1, 'sid')
        nsm_type = string(card, 2, 'Type')
        id = integer(card, 3 + noffset, 'id')
        value = double(card, 4 + noffset, 'value')
        return NSM(sid, nsm_type, id, value, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds an NSM card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid = data[0]
        #sid=9  prop_set=PBEA ID=538976333 value=0.0
        #sid=10 prop_set=PDUM ID=538976312 value=2.80259692865e-45
        #sid=10 prop_set=ELEM ID=542395973 value=0.0
        nsm_type = data[1]
        id = data[2]
        value = data[3]
        return NSM(sid, nsm_type, id, value, comment=comment)

    def cross_reference(self, model):
        pass

    @property
    def Type(self):
        """gets the nsm_type"""
        return self.nsm_type
    @Type.setter
    def Type(self, nsm_type):
        """sets the nsm_type"""
        self.nsm_type = nsm_type

    def raw_fields(self):
        #nodes = self.node_ids
        list_fields = [self.type, self.sid, self.nsm_type, self.id, self.value]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class NSML(NSM):
    """
    Defines a set of lumped non structural mass.

    +------+-----+------+----+-------+----+-------+----+-------+
    |  1   |  2  |  3   |  4 |   5   | 6  |   7   | 8  |   9   |
    +======+=====+======+====+=======+====+=======+====+=======+
    | NSML | SID | TYPE | ID | VALUE | ID | VALUE | ID | VALUE |
    +------+-----+------+----+-------+----+-------+----+-------+
    """
    type = 'NSML'
    def __init__(self, sid, nsm_type, id, value, comment=''):
        """
        Creates an NSML card, which defines lumped non-structural mass

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        id : int
            property id or element id depending on nsm_type
        value : float
            the non-structural pass per unit length/area
        comment : str; default=''
            a comment for the card
        """
        NSM.__init__(self, sid, nsm_type, id, value, comment=comment)

class NSML1(NSM1):
    """
    Defines lumped non structural mass entries by VALUE,ID list.

    +-------+-----+------+-------+-----+----+----+----+----+
    |   1   |  2  |  3   |   4   |  5  | 6  | 7  | 8  | 9  |
    +=======+=====+======+=======+=====+====+====+====+====+
    | NSML1 | SID | TYPE | VALUE | ID  | ID | ID | ID | ID |
    +-------+-----+------+-------+-----+----+----+----+----+
    |       |  ID |  ID  |  ID   | etc |    |    |    |    |
    +-------+-----+------+-------+-----+----+----+----+----+
    """
    type = 'NSML1'
    def __init__(self, sid, nsm_type, ids, value, comment=''):
        """
        Creates an NSML card, which defines lumped non-structural mass

        Parameters
        ----------
        sid : int
            Case control NSM id
        nsm_type : str
            Type of card the NSM is applied to
            valid_properties = {
                PSHELL, PCOMP, PBAR, PBARL, PBEAM, PBEAML, PBCOMP,
                PROD, CONROD, PBEND, PSHEAR, PTUBE, PCONEAX, PRAC2D,
                ELEMENT
            }
        id : int
            property id or element id depending on nsm_type
        value : float
            the non-structural pass per unit length/area
        comment : str; default=''
            a comment for the card
        """
        NSM1.__init__(self, sid, nsm_type, value, ids, comment=comment)


class NSMADD(BaseCard):
    """
    Defines non structural mass as the sum of the sets listed.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | NSMADD | 2  |  1 |  3  |
    +--------+----+----+-----+
    """
    type = 'NSMADD'

    def __init__(self, sid, sets, comment=''):
        """
        Creates an NSMADD card, which sum NSM sets

        Parameters
        ----------
        sid : int
            the NSM Case Control value
        sets : List[int]
            the NSM, NSM1, NSML, NSML1 values
        comment : str; default=''
            a comment for the card
        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.sets = expand_thru(sets)
        self.sets.sort()
        self.sets_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a NSMADD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        sets = card.fields(2)
        return NSMADD(sid, sets, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds an NSMADD card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        sid = data[0]
        sets = list(data[1:-1])
        return NSMADD(sid, sets, comment=comment)

    @property
    def nsm_ids(self):
        """gets the nonstructural-mass ids"""
        if self.sets_ref is None:
            return self.sets

        nsm_ids = []
        for nsm in self.sets_ref:
            if isinstance(nsm, integer_types):
                nsm_ids.append(nsm)
            elif isinstance(nsm, list):
                nsm_ids.append(nsm[0].sid)
            else:
                raise TypeError('type=%s; nsm=\n%s' % (type(nsm), nsm))
        return nsm_ids

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by NSMADD=%s' % self.sid
        nsms = []
        for i, nsm in enumerate(self.sets):
            nsms.append(model.NSM(nsm, msg=msg))
        self.sets_ref = nsms

    def safe_cross_reference(self, model, debug=True):
        nsms = []
        msg = ' which is required by NSMADD=%s' % self.sid
        for nsm_id in self.sets:
            try:
                nsm = model.NSM(nsm_id, msg=msg)
            except KeyError:
                if debug:
                    msg = 'Couldnt find NSM=%i, which is required by NSMADD=%s' % (
                        nsm_id, self.sid)
                    print(msg)
                continue
            nsms.append(nsm)
        self.sets_ref = nsms

    def uncross_reference(self):
        self.sets = self.nsm_ids
        self.sets_ref = None

    def raw_fields(self):
        fields = ['NSMADD', self.sid] + self.nsm_ids
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)


class PMASS(Property):
    type = 'PMASS'
    _field_map = {
        1: 'pid', 2:'mass',
    }

    def __init__(self, pid, mass, comment=''):
        """
        Creates an PMASS card, which defines a mass applied to a single DOF

        Parameters
        ----------
        pid : int
            Property id used by a CMASS1/CMASS3 card
        mass : float
            the mass to apply
        comment : str; default=''
            a comment for the card
        """
        Property.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.mass = mass

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds a PMASS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            the index of the card that's being parsed
        comment : str; default=''
            a comment for the card
        """
        icard *= 2
        #: Property ID
        pid = integer(card, 1 + icard, 'pid')
        mass = double_or_blank(card, 2 + icard, 'mass', 0.)
        return PMASS(pid, mass, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PMASS card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        pid = data[0]
        mass = data[1]
        return PMASS(pid, mass, comment=comment)

    def cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass

    def _verify(self, xref=False):
        pid = self.Pid()
        mass = self.Mass()
        assert isinstance(pid, integer_types), 'pid=%r' % pid
        assert isinstance(mass, float), 'mass=%r' % mass

    def Mass(self):
        return self.mass

    def raw_fields(self):
        list_fields = ['PMASS', self.pid, self.mass]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
