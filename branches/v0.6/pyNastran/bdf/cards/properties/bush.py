# pylint: disable=C0103,R0902,R0904,R0915
"""
All bush properties are defined in this file.  This includes:

 *   PBUSH
 *   PBUSH1D
 *   PBUSH2D (not implemented)
 *   PBUSHT (not implemented)

All bush properties are BushingProperty and Property objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.baseCard import Property
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank,
    string, string_or_blank, blank, fields)


class BushingProperty(Property):
    type = 'BushingProperty'

    def __init__(self, card, data):
        Property.__init__(self, card, data)

    def cross_reference(self, model):
        pass


class PBUSH(BushingProperty):
    type = 'PBUSH'
    _field_map = {
        1: 'pid',
    }

    def __init__(self, card=None, data=None, comment=''):
        BushingProperty.__init__(self, card, data)
        if comment:
            self._comment = comment

        # K parameter
        self.Ki = []

        # B parameter
        self.Bi = []

        # GE parameter
        self.GEi = []

        # RCV parameters
        self.sa = None
        self.st = None
        self.ea = None
        self.et = None

        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')

            nFields = card.nFields()
            self.vars = []
            iStart = 2
            while iStart < nFields:
                pname = string(card, iStart, 'pname')
                if   pname == 'K':
                    self.getK(card, iStart)
                elif pname == 'B':
                    self.getB(card, iStart)
                elif pname == 'GE':
                    self.getGE(card, iStart)
                elif pname == 'RCV':
                    self.getRCV(card, iStart)
                else:
                    break
                iStart += 8
        else:
            self.pid = data[0]
            self.b = data[1]
            raise NotImplementedError('PBUSH data...')
        #print self

    def _verify(self, xref=False):
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r' % pid

    def getK(self, card, iStart):
        # Flag indicating that the next 1 to 6 fields are stiffness values in
        # the element coordinate system.
        #self.k = string(card, iStart, 'k')

        #: Nominal stiffness values in directions 1 through 6.
        #: See Remarks 2 and 3. (Real; Default = 0.0)
        self.Ki = fields(double_or_blank, card, 'Ki', iStart + 1, iStart + 7)
        #print "Ki = ",self.Ki
        self.vars.append('K')

    def getB(self, card, iStart):
        # Flag indicating that the next 1 to 6 fields are force-per-velocity
        # damping.
        #self.b = string(card, iStart, 'b')

        #: Force per unit velocity (Real)
        #: Nominal damping coefficients in direction 1 through 6 in units of
        #: force per unit velocity. See Remarks 2, 3, and 9. (Real; Default=0.)
        self.Bi = fields(double_or_blank, card, 'Bi', iStart + 1, iStart + 7)
        self.vars.append('B')

    def getGE(self, card, iStart):
        # Flag indicating that the next fields, 1 through 6 are structural
        # damping constants. See Remark 7. (Character)
        #self.ge = string(card, iStart, 'ge')

        #: Nominal structural damping constant in directions 1 through 6. See
        #: Remarks 2. and 3. (Real; Default = 0.0)
        self.GEi = fields(double_or_blank, card, 'GEi', iStart + 1, iStart + 7)
        self.vars.append('GE')

    def getRCV(self, card, iStart):
        # Flag indicating that the next 1 to 4 fields are stress or strain
        # coefficients. (Character)
        #self.rcv = string(card, iStart, 'rcv')

        self.sa = double_or_blank(card, iStart + 1, 'sa', 1.)
        self.st = double_or_blank(card, iStart + 2, 'st', 1.)
        self.ea = double_or_blank(card, iStart + 3, 'ea', 1.)
        self.et = double_or_blank(card, iStart + 4, 'et', 1.)
        self.vars.append('RCV')

    def rawFields(self):
        list_fields = ['PBUSH', self.pid]
        for var in self.vars:
            if var == 'K':
                list_fields += ['K'] + self.Ki
            elif var == 'B':
                list_fields += ['B'] + self.Bi
            elif var == 'GE':
                list_fields += ['GE'] + self.GEi
            elif var == 'RCV':
                list_fields += ['RCV', self.sa, self.st, self.ea, self.et]
            else:
                raise RuntimeError('not supported PBUSH field...')
            nSpaces = 8 - (len(list_fields) - 1) % 8

            #print("nSpaces = ",nSpaces)
            if nSpaces == 8:
                list_fields += [None]
            elif nSpaces < 8:
                list_fields += [None] * (nSpaces + 1)
        return list_fields

    def reprFields(self):
        return self.rawFields()


class PBUSH1D(BushingProperty):
    type = 'PBUSH1D'

    def __init__(self, card=None, data=None, comment=''):
        BushingProperty.__init__(self, card, data)
        if comment:
            self._comment = comment

        # SPRING parameters
        self.springType = None
        self.springIDT = None
        self.springIDC = None
        self.springIDTDU = None
        self.springIDCDU = None

        # DAMPER parameters
        self.damperType = None
        self.damperIDT = None
        self.damperIDC = None
        self.damperIDTDV = None
        self.damperIDCDV = None

        # GENER parameters
        self.generIDT = None
        self.generIDC = None
        self.generIDTDU = None
        self.generIDCDU = None
        self.generIDTDV = None
        self.generIDCDV = None

        # SHOCK parameters
        self.shockType = None
        self.shockCVT = None
        self.shockCVC = None
        self.shockExpVT = None
        self.shockExpVC = None
        self.shockIDTS = None

        self.shockIDETS = None
        self.shockIDECS = None
        self.shockIDETSD = None
        self.shockIDECSD = None

        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            self.k = double_or_blank(card, 2, 'k', 0.0)
            self.c = double_or_blank(card, 3, 'c', 0.0)
            self.m = double_or_blank(card, 4, 'm', 0.0)

            self.sa = double_or_blank(card, 6, 'sa', 0.0)
            self.se = double_or_blank(card, 7, 'se', 0.0)

            nFields = card.nFields()
            self.vars = []
            iStart = 9
            while iStart < nFields:
                pname = string(card, iStart, 'pname')
                if   pname == 'SHOCKA':
                    iStart = self.getShockA(card, iStart)
                elif pname == 'SPRING':
                    self.getSpring(card, iStart)
                elif pname == 'DAMPER':
                    self.getDamper(card, iStart)
                elif pname == 'GENER':
                    self.getGener(card, iStart)
                else:
                    break
                iStart += 8
        else:
            self.pid = data[0]
            self.b = data[1]
            raise NotImplementedError('PBUSH1D data...')

    def _verify(self, xref=False):
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r' % pid

    def getShockA(self, card, iStart):
        self.shockType = string_or_blank(card, iStart + 1, 'shockType')
        self.shockCVT = double(card, iStart + 2, 'shockCVT')
        self.shockCVC = double_or_blank(card, iStart + 3, 'shockCVC')
        self.shockExpVT = double_or_blank(card, iStart + 4, 'shockExpVT', 1.0)
        self.shockExpVC = double_or_blank(card, iStart + 5,
                                          'shockExpVC', self.shockExpVT)

        if self.shockType == 'TABLE':
            pass
            # self.shockIDTS = integer(card, iStart + 6, 'shockIDTS')
            # self.shockIDETS = blank(card, iStart + 9, 'shockIDETS')
            # self.shockIDECS = blank(card, iStart + 10, 'shockIDECS')
            # self.shockIDETSD = blank(card, iStart + 11, 'shockIDETSD')
            # self.shockIDECSD = blank(card, iStart + 12, 'shockIDECSD')
        elif self.shockType == 'EQUAT':
            self.shockIDTS = blank(card, iStart + 6, 'shockIDTS')
            self.shockIDETS = integer(card, iStart + 9, 'shockIDETS')
            self.shockIDECS = integer_or_blank(card, iStart + 10,
                                               'shockIDECS', self.shockIDETS)
            self.shockIDETSD = integer(card, iStart + 11, 'shockIDETSD')
            self.shockIDECSD = integer_or_blank(card, iStart + 11,
                                                'shockIDECSD', self.shockIDETSD)
        else:
            raise RuntimeError('Invalid shockType=|%s| on card\n%s' %(self.shockType, card))

        iStart += 8
        return iStart

    def getSpring(self, card, iStart):
        self.springType = string_or_blank(card, iStart + 1, 'springType')
        self.springIDT = integer(card, iStart + 2, 'springIDT')

        if self.springType == 'TABLE':
            self.springIDC = blank(card, iStart + 3, 'springIDC')
            self.springIDTDU = blank(card, iStart + 4, 'springIDTDU')
            self.springIDCDU = blank(card, iStart + 5, 'springIDCDU')
        elif self.springType == 'EQUAT':
            self.springIDC = integer_or_blank(card, iStart + 3,
                                              'springIDC', self.springIDT)
            self.springIDTDU = integer(card, iStart + 4, 'springIDTDU')
            self.springIDCDU = integer_or_blank(card, iStart + 5,
                                                'springIDCDU', self.springIDTDU)
        else:
            raise RuntimeError('Invalid springType=|%s| on card\n%s' %(self.springType, card))

        self.vars.append('SPRING')

    def getDamper(self, card, iStart):
        self.damperType = string_or_blank(card, iStart + 1, 'damperType')
        self.damperIDT = integer(card, iStart + 2, 'damperIDT')
        if self.damperType == 'TABLE':
            self.damperIDC = blank(card, iStart + 3, 'damperIDC')
            self.damperIDTDV = blank(card, iStart + 4, 'damperIDTDV')
            self.damperIDCDV = blank(card, iStart + 5, 'damperIDCDV')
        elif self.damperType == 'EQUAT':
            self.damperIDC = integer_or_blank(card, iStart + 3, 'damperIDC')
            self.damperIDTDV = integer(card, iStart + 4, 'damperIDTDV')
            self.damperIDCDV = integer_or_blank(card, iStart + 5, 'damperIDCDV', self.damperIDTDV)
        else:
            raise RuntimeError('Invalid springType=|%s| on card\n%s' %(self.springType, card))

        self.vars.append('DAMPER')

    def getGener(self, card, iStart):
        self.generIDT = integer(card, iStart + 2, 'generIDT')
        self.generIDC = integer_or_blank(card, iStart + 3,
                                         'generIDC', self.generIDT)
        self.generIDTDU = integer(card, iStart + 4, 'generIDTDU')
        self.generIDCDU = integer_or_blank(card, iStart + 5,
                                           'generIDCDU', self.generIDTDU)
        self.generIDTDV = integer(card, iStart + 6, 'generIDTDV')
        self.generIDCDV = integer_or_blank(card, iStart + 7,
                                           'generIDCDV', self.generIDTDV)
        self.vars.append('GENER')

    def _shockFields(self):
        list_fields = ['SHOCKA', self.shockType, self.shockCVT, self.shockCVC,
                  self.shockExpVT, self.shockExpVC, self.shockIDTS, None, None,
                  self.shockIDETS, self.shockIDECS, self.shockIDETSD,
                  self.shockIDECSD]
        return list_fields

    def _springFields(self):
        list_fields = ['SPRING', self.springType, self.springIDT,
                       self.springIDC, self.springIDTDU, self.springIDCDU]
        return list_fields

    def _damperFields(self):
        list_fields = ['DAMPER', self.damperType, self.damperIDT,
                       self.damperIDC, self.damperIDTDV, self.damperIDCDV]
        return list_fields

    def _generFields(self):
        list_fields = ['GENER', None, self.generIDT, self.generIDC,
                  self.generIDTDU, self.generIDCDU, self.generIDTDV,
                  self.generIDCDV]
        return list_fields

    def rawFields(self):
        list_fields = ['PBUSH1D', self.pid, self.k, self.c, self.m, None,
                  self.sa, self.se, None]
        for var in self.vars:
            if var == 'SHOCKA':
                list_fields += self._shockFields()
            elif var == 'SPRING':
                list_fields += self._springFields()
            elif var == 'DAMPER':
                list_fields += self._damperFields()
            elif var == 'GENER':
                list_fields += self._generFields()
            else:
                msg = 'var=%s not supported PBUSH1D field...' % var
                raise RuntimeError(msg)
            nSpaces = 8 - (len(list_fields) - 1) % 8

            if nSpaces < 8:
                list_fields += [None] * (nSpaces)
        return list_fields

    def reprFields(self):
        return self.rawFields()


class PBUSH2D(BushingProperty):
    type = 'PBUSH2D'

    def __init__(self, card=None, data=None, comment=''):
        BushingProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            raise NotImplementedError()
        else:
            raise NotImplementedError()


class PBUSHT(BushingProperty):
    type = 'PBUSHT'

    def __init__(self, card=None, data=None, comment=''):
        BushingProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            raise NotImplementedError()
        else:
            raise NotImplementedError()
