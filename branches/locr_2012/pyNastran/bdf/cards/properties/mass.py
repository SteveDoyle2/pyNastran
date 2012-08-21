# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.baseCard import Property


class PointProperty(Property):
    type = 'PointProperty'

    def __init__(self, card, data):
        Property.__init__(self, card, data)

    def cross_reference(self, model):
        pass


class NSM(PointProperty):
    """
    Defines a set of non structural mass.
    """
    ## Set points to either Property entries or Element entries.
    ## Properties are:
    validProperties = [
        'PSHELL', 'PCOMP', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP',
        'PROD', 'CONROD', 'PBEND', 'PSHEAR', 'PTUBE', 'PCONEAX', 'PRAC2D']

    def __init__(self, card=None, nOffset=0, data=None):
        PointProperty.__init__(self, card, data)
        if card:
            nOffset *= 2
            self.sid = card.field(1)
            self.Type = card.field(2)
            self.id = card.field(3 + nOffset)
            self.value = card.field(4 + nOffset)
        else:
            self.sid = data[0]
            #sid=9  propSet=PBEA ID=538976333 value=0.0
            #sid=10 propSet=PDUM ID=538976312 value=2.80259692865e-45
            #sid=10 propSet=ELEM ID=542395973 value=0.0
            self.Type = data[1]
            self.id = data[2]
            self.value = data[3]
        ###

    def rawFields(self):
        #nodes = self.nodeIDs()
        fields = ['NSM', self.sid, self.Type, self.id, self.value]
        return fields

    def reprFields(self):
        return self.rawFields()


class PMASS(PointProperty):
    def __init__(self, card=None, nOffset=0, data=None):
        PointProperty.__init__(self, card, data)
        if card:
            nOffset *= 2
            ## Property ID
            self.pid = card.field(1 + nOffset)
            self.mass = card.field(2 + nOffset, 0.)
        else:
            self.pid = data[0]
            self.mass = data[1]
        ###

    def Mass(self):
        return self.mass

    def rawFields(self):
        fields = ['PMASS', self.pid, self.mass]
        return fields

    def reprFields(self):
        return self.rawFields()
