# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.bdf.cards.baseCard import Property


class BushingProperty(Property):
    type = 'BushingProperty'

    def __init__(self, card, data):
        Property.__init__(self, card, data)

    def cross_reference(self, model):
        pass


class PBUSH(BushingProperty):
    type = 'PBUSH'

    def __init__(self, card=None, data=None):
        BushingProperty.__init__(self, card, data)

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
            ## Property ID
            self.pid = card.field(1)

            nFields = card.nFields()
            self.vars = []
            iStart = 2
            while iStart < nFields:
                pname = card.field(iStart)
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
            ###
        else:
            self.pid = data[0]
            self.b = data[1]
            raise NotImplementedError('PBUSH data...')
        ###
        #print self

    def getK(self, card, iStart):
        ## Flag indicating that the next 1 to 6 fields are stiffness values in
        ## the element coordinate system.
        #self.k = card.field(iStart)
        ## Nominal stiffness values in directions 1 through 6.
        ## See Remarks 2 and 3. (Real; Default = 0.0)
        self.Ki = card.fields(i=iStart + 1, j=iStart + 6)
        #print "Ki = ",self.Ki
        self.vars.append('K')

    def getB(self, card, iStart):
        ## Flag indicating that the next 1 to 6 fields are force-per-velocity
        ## damping.
        #self.b = card.field(iStart)
        ## Force per unit velocity (Real)
        ## Nominal damping coefficients in direction 1 through 6 in units of
        ## force per unit velocity. See Remarks 2, 3, and 9. (Real; Default=0.)
        self.Bi = card.fields(i=iStart + 1, j=iStart + 6)
        self.vars.append('B')

    def getGE(self, card, iStart):
        ## Flag indicating that the next fields, 1 through 6 are structural
        ## damping constants. See Remark 7. (Character)
        #self.ge = card.field(iStart)
        ## Nominal structural damping constant in directions 1 through 6. See
        ## Remarks 2. and 3. (Real; Default = 0.0)
        self.GEi = card.fields(i=iStart + 1, j=iStart + 6)
        self.vars.append('GE')

    def getRCV(self, card, iStart):
        ## Flag indicating that the next 1 to 4 fields are stress or strain
        ## coefficients. (Character)
        #self.ge = card.field(iStart)
        self.sa = card.field(iStart + 1, 1.)
        self.st = card.field(iStart + 2, 1.)
        self.ea = card.field(iStart + 3, 1.)
        self.et = card.field(iStart + 4, 1.)
        self.vars.append('RCV')

    def rawFields(self):
        fields = ['PBUSH', self.pid]
        for var in self.vars:
            if var == 'K':
                fields += ['K'] + self.Ki
            elif var == 'B':
                fields += ['B'] + self.Bi
            elif var == 'GE':
                fields += ['GE'] + self.GEi
            elif var == 'RCV':
                fields += ['RCV', self.sa, self.st, self.ea, self.et]
            else:
                raise RuntimeError('not supported PBUSH field...')
            nSpaces = 8 - (len(fields) - 1) % 8

            #print "nSpaces = ",nSpaces
            if nSpaces < 8:
                fields += [None] * (nSpaces + 1)
            ###
        return fields

    def reprFields(self):
        return self.rawFields()


class PBUSH1D(BushingProperty):
    type = 'PBUSH1D'

    def __init__(self, card=None, data=None):
        BushingProperty.__init__(self, card, data)

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
            ## Property ID
            self.pid = card.field(1)
            self.k = card.field(2)
            self.c = card.field(3)
            self.m = card.field(4)
            #
            self.sa = card.field(6)
            self.se = card.field(7)

            nFields = card.nFields()
            self.vars = []
            iStart = 9
            while iStart < nFields:
                pname = card.field(iStart)
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
            ###
        else:
            self.pid = data[0]
            self.b = data[1]
            raise NotImplementedError('PBUSH1D data...')
        ###

    def getShockA(self, card, iStart):
        self.shockType = card.field(iStart + 1)
        self.shockCVT = card.field(iStart + 2)
        self.shockCVC = card.field(iStart + 3)
        self.shockExpVT = card.field(iStart + 4)
        self.shockExpVC = card.field(iStart + 5)
        self.shockIDTS = card.field(iStart + 6)

        self.shockIDETS = card.field(iStart + 9)
        self.shockIDECS = card.field(iStart + 10)
        self.shockIDETSD = card.field(iStart + 11)
        self.shockIDECSD = card.field(iStart + 12)
        iStart += 8
        return iStart

    def getSpring(self, card, iStart):
        self.springType = card.field(iStart + 1)
        self.springIDT = card.field(iStart + 2)
        self.springIDC = card.field(iStart + 3)
        self.springIDTDU = card.field(iStart + 4)
        self.springIDCDU = card.field(iStart + 5)
        self.vars.append('SPRING')

    def getDamper(self, card, iStart):
        self.damperType = card.field(iStart + 1)
        self.damperIDT = card.field(iStart + 2)
        self.damperIDC = card.field(iStart + 3)
        self.damperIDTDV = card.field(iStart + 4)
        self.damperIDCDV = card.field(iStart + 5)
        self.vars.append('DAMPER')

    def getGener(self, card, iStart):
        self.generIDT = card.field(iStart + 2)
        self.generIDC = card.field(iStart + 3)
        self.generIDTDU = card.field(iStart + 4)
        self.generIDCDU = card.field(iStart + 5)
        self.generIDTDV = card.field(iStart + 6)
        self.generIDCDV = card.field(iStart + 7)
        self.vars.append('GENER')

    def _shockFields(self):
        fields = ['SHOCKA', self.shockType, self.shockCVT, self.shockCVC, self.shockExpVT, self.shockExpVC, self.shockIDTS,
                  None, None, self.shockIDETS, self.shockIDECS, self.shockIDETSD, self.shockIDECSD]
        return fields

    def _springFields(self):
        fields = ['SPRING', self.springType, self.springIDT,
                  self.springIDC, self.springIDTDU, self.springIDCDU]
        return fields

    def _damperFields(self):
        fields = ['DAMPER', self.damperType, self.damperIDT,
                  self.damperIDC, self.damperIDTDV, self.damperIDCDV]
        return fields

    def _generFields(self):
        fields = ['GENER', None, self.generIDT, self.generIDC, self.generIDTDU,
                  self.generIDCDU, self.generIDTDV, self.generIDCDV]
        return fields

    def rawFields(self):
        fields = ['PBUSH1D', self.pid, self.k, self.c, self.m,
                  None, self.sa, self.se, None]
        for var in self.vars:
            if var == 'SHOCKA':
                fields += self._shockFields()
            elif var == 'SPRING':
                fields += self._springFields()
            elif var == 'DAMPER':
                fields += self._damperFields()
            elif var == 'GENER':
                fields += self._generFields()
            else:
                raise RuntimeError('var=%s not supported PBUSH1D field...' %
                                   (var))
            nSpaces = 8 - (len(fields) - 1) % 8

            #print "nSpaces = ",nSpaces
            if nSpaces < 8:
                fields += [None] * (nSpaces)
            ###
        return fields

    def reprFields(self):
        return self.rawFields()

#class PBUSH2D
#class PBUSHT
