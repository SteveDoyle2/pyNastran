# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All mass properties are defined in this file.  This includes:
 * NSM
 * PMASS

All mass properties are PointProperty and Property objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.baseCard import Property
from pyNastran.bdf.assign_type import (integer,
    double, double_or_blank, string)


class PointProperty(Property):
    def __init__(self, card, data):
        Property.__init__(self, card, data)

    def cross_reference(self, model):
        pass


class NSM(PointProperty):
    """
    Defines a set of non structural mass.
    """
    type = 'NSM'
    ## Set points to either Property entries or Element entries.
    ## Properties are:
    validProperties = [
        'PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP',
        'PROD', 'CONROD', 'PBEND', 'PSHEAR', 'PTUBE', 'PCONEAX', 'PRAC2D']

    def __init__(self, card=None, nOffset=0, data=None, comment=''):
        PointProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            nOffset *= 2
            self.sid = integer(card, 1, 'sid')
            self.Type = string(card, 2, 'Type')
            self.id = integer(card, 3 + nOffset, 'id')
            self.value = double(card, 4 + nOffset, 'value')
        else:
            self.sid = data[0]
            #sid=9  propSet=PBEA ID=538976333 value=0.0
            #sid=10 propSet=PDUM ID=538976312 value=2.80259692865e-45
            #sid=10 propSet=ELEM ID=542395973 value=0.0
            self.Type = data[1]
            self.id = data[2]
            self.value = data[3]
        assert self.Type in self.validProperties

    def rawFields(self):
        #nodes = self.nodeIDs()
        list_fields = ['NSM', self.sid, self.Type, self.id, self.value]
        return list_fields

    def reprFields(self):
        return self.rawFields()


class PMASS(PointProperty):
    type = 'PMASS'
    def __init__(self, card=None, nOffset=0, data=None, comment=''):
        PointProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            nOffset *= 2
            ## Property ID
            self.pid = integer(card, 1 + nOffset, 'pid')
            self.mass = double_or_blank(card, 2 + nOffset, 'mass', 0.)
        else:
            self.pid = data[0]
            self.mass = data[1]

    def _verify(self, isxref=False):
        pid = self.Pid()
        mass = self.Mass()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(mass, float), 'mass=%r' % mass

    def Mass(self):
        return self.mass

    def rawFields(self):
        list_fields = ['PMASS', self.pid, self.mass]
        return list_fields

    def reprFields(self):
        return self.rawFields()