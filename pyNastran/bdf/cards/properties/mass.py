# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All mass properties are defined in this file.  This includes:

 * NSM
 * PMASS

All mass properties are PointProperty and Property objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, double_or_blank, string)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class PointProperty(Property):
    def __init__(self):
        Property.__init__(self)

    def cross_reference(self, model):
        pass


class NSM(PointProperty):
    """
    Defines a set of non structural mass.
    """
    type = 'NSM'
    _field_map = {
        1: 'sid', 2:'Type', 3:'id', 4:'value'
    }

    #: Set points to either Property entries or Element entries.
    #: Properties are:
    validProperties = [
        'PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP',
        'PROD', 'CONROD', 'PBEND', 'PSHEAR', 'PTUBE', 'PCONEAX', 'PRAC2D']

    def __init__(self, sid, Type, id, value, comment=''):
        PointProperty.__init__(self, card, data)
        if comment:
            self.comment = comment
        self.sid = sid
        self.Type = Type
        self.id = id
        self.value = value
        assert self.Type in self.validProperties

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
        list_fields = ['NSM', self.sid, self.Type, self.id, self.value]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PMASS(PointProperty):
    type = 'PMASS'
    _field_map = {
        1: 'pid', 2:'mass',
    }

    def __init__(self, pid, mass, comment=''):
        PointProperty.__init__(self)
        self.pid = pid
        self.mass = mass

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        if comment:
            self.comment = comment
        icard *= 2
        #: Property ID
        pid = integer(card, 1 + icard, 'pid')
        mass = double_or_blank(card, 2 + icard, 'mass', 0.)
        return PMASS(pid, mass, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        pid = data[0]
        mass = data[1]
        return PMASS(pid, mass, comment=comment)

    def _verify(self, xref=False):
        pid = self.Pid()
        mass = self.Mass()
        assert isinstance(pid, int), 'pid=%r' % pid
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
