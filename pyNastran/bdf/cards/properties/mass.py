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


class PointProperty(Property):
    def __init__(self):
        Property.__init__(self)

    def cross_reference(self, model):
        pass

class NSM1(PointProperty):
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

    def __init__(self, sid, Type, value, ids, comment=''):
        PointProperty.__init__(self)
        if comment:
            self.comment = comment
        if isinstance(ids, integer_types):
            ids = [ids]
        self.sid = sid
        self.Type = Type
        self.ids = ids
        self.value = value
        if self.Type not in self.valid_properties:
            msg = 'Type=%r must be in [%s]' % (self.Type, ', '.join(self.valid_properties))
            raise TypeError(msg)
        assert isinstance(ids, list), 'ids=%r is not a list' % (ids)


    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        Type = string(card, 2, 'Type')
        value = double(card, 3, 'value')
        ids = card[4:]
        return NSM1(sid, Type, value, ids, comment=comment)

    def raw_fields(self):
        #nodes = self.node_ids
        list_fields = ['NSM1', self.sid, self.Type, self.value] + self.ids
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class NSM(PointProperty):
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

    def __init__(self, sid, Type, id, value, comment=''):
        PointProperty.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.Type = Type
        self.id = id
        self.value = value
        if self.Type not in self.valid_properties:
            msg = 'Type=%r must be in [%s]' % (self.Type, ', '.join(self.valid_properties))
            raise TypeError(msg)

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        noffset = 2 * icard
        sid = integer(card, 1, 'sid')
        Type = string(card, 2, 'Type')
        id = integer(card, 3 + noffset, 'id')
        value = double(card, 4 + noffset, 'value')
        return NSM(sid, Type, id, value, comment=comment)

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
        Type = data[1]
        id = data[2]
        value = data[3]
        return NSM(sid, Type, id, value, comment=comment)

    def raw_fields(self):
        #nodes = self.node_ids
        list_fields = [self.type, self.sid, self.Type, self.id, self.value]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class NSML(NSM):
    type = 'NSML'
    def __init__(self, sid, Type, id, value, comment=''):
        NSM.__init__(self, sid, Type, id, value, comment=comment)

class NSML1(NSM1):
    type = 'NSML1'
    def __init__(self, sid, Type, ids, value, comment=''):
        NSM1.__init__(self, sid, Type, value, ids, comment=comment)

class NSMADD(BaseCard):
    """
    Defines a single-point constraint set as a union of single-point constraint
    sets defined on SPC or SPC1 entries.

    +--------+----+----+-----+
    |    1   | 2  |  3 |  4  |
    +========+====+====+=====+
    | NSMADD | 2  |  1 |  3  |
    +--------+----+----+-----+
    """
    type = 'SPCADD'

    def __init__(self, sid, sets, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.sid = sid
        self.sets = expand_thru(sets)
        self.sets.sort()

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
        nsm_ids = []
        for nsm in self.sets:
            if isinstance(nsm, integer_types):
                nsm_ids.append(nsm)
            elif isinstance(nsm, list):
                nsm_ids.append(nsm[0].sid)
            else:
                raise TypeError('type=%s; nsm=\n%s' % (type(nsm), nsm))
        return nsm_ids

    @property
    def ids(self):
        return self.nsm_ids

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by NSMADD=%s' % self.sid
        for i, nsm in enumerate(self.sets):
            self.sets[i] = model.NSM(nsm, msg=msg)
        self.sets_ref = self.sets

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
        self.sets = nsms
        self.sets_ref = self.sets

    def uncross_reference(self):
        self.sets = self.nsm_ids
        del self.sets_ref

    def raw_fields(self):
        fields = ['NSMADD', self.sid] + self.nsm_ids
        return fields

    def write_card(self, size=8, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.raw_fields()
        return self.comment + print_card_16(card)


class PMASS(PointProperty):
    type = 'PMASS'
    _field_map = {
        1: 'pid', 2:'mass',
    }

    def __init__(self, pid, mass, comment=''):
        PointProperty.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.mass = mass

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
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
