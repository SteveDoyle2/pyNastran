# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.general.general import ListPrint


class Table(BaseCard):
    type = 'TABLE??'

    def __init__(self, card, data):
        pass

    def mapAxis(self, axis):
        if axis == 0:
            axisType = 'LINEAR'
        else:
            raise ValueError('axis=|%s|' % (axis))
        return axisType

    def parseFields(self, fields, nRepeated, isData=False):
        self.table = TableObj(fields, nRepeated, isData)


class TableObj(object):
    def __init__(self, fields, nRepeated, isData=False):
        self.table = []
        fields = self.cleanupFields(fields, isData)

        nFields = len(fields)

        if not isData:
            if nFields % nRepeated != 0:
                self.crashFields(fields, nRepeated, nFields)
            if nFields % nRepeated != 0:
                msg = 'invalid table length nRepeat=%s fields=%s' % (
                    nRepeated, ListPrint(fields))
                raise RuntimeError(msg)

        i = 0
        while i < nFields:
            pack = []
            for j in xrange(nRepeated):
                pack.append(fields[i + j])
            i += nRepeated
            self.table.append(pack)

    def crashFields(self, fields, nRepeated, nFields):
        try:
            msg = ''
            #k = 0
            for i in xrange(nFields):
                for j in xrange(nRepeated):
                    try:
                        msg += '%-8g ' % (fields[i * nRepeated + j])
                    except TypeError:
                        msg += '*%-8s ' % (fields[i * nRepeated + j])
                    except IndexError:
                        msg += 'IndexError'
                    #k+=1
                msg += '\n'
        except:
            print(fields)
            print(msg)
            #assert nFields%nRepeated==0,msg
            raise

    def cleanupFields(self, fields, isData=False):
        fields2 = []  # remove extra ENDTs
        foundENDT = False
        for field in fields:
            if isinstance(field, unicode) and 'ENDT' in field.upper():
                foundENDT = True
            else:
                fields2.append(field)

        if not isData:
            assert foundENDT == True, fields
        return fields2

    def fields(self):
        fields = []
        for pack in self.table:
            fields += pack
        return fields


class TABLED1(Table):
    type = 'TABLED1'

    def __init__(self, card=None, data=None):
        Table.__init__(self, card, data)
        if card:
            self.tid = card.field(1)
            self.xaxis = card.field(2, 'LINEAR')
            self.yaxis = card.field(3, 'LINEAR')
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.xaxis = self.mapAxis(data[1])
            self.yaxis = self.mapAxis(data[2])
            fields = data[3:]
            isData = True
        assert self.xaxis in ['LINEAR', 'LOG'], 'xaxis=|%s|' % (self.xaxis)
        assert self.yaxis in ['LINEAR', 'LOG'], 'yaxis=|%s|' % (self.yaxis)
        self.parseFields(fields, nRepeated=2, isData=isData)

    def rawFields(self):
        fields = ['TABLED1', self.tid, self.xaxis, self.yaxis, None,
                  None, None, None, None] + self.table.fields() + ['ENDT']
        return fields

    def reprFields(self):
        #xaxis = set_blank_if_default(self.xaxis,'LINEAR')
        #yaxis = set_blank_if_default(self.yaxis,'LINEAR')
        return self.rawFields()
        #fields = ['TABLED1',self.tid,self.xaxis,self.yaxis,None,None,None,None,None]+self.table.fields()+['ENDT']
        #return fields


class TABLED2(Table):
    type = 'TABLED2'

    def __init__(self, card=None, data=None):
        Table.__init__(self, card, data)
        if card:
            self.tid = card.field(1)
            self.x1 = card.field(2)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            fields = data[2:]
            isData = True
        self.parseFields(fields, nRepeated=2, isData=isData)

    def rawFields(self):
        fields = ['TABLED2', self.tid, self.x1, None, None, None,
                  None, None, None] + self.table.fields() + ['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()


class TABLED3(Table):
    type = 'TABLED3'

    def __init__(self, card=None, data=None):
        Table.__init__(self, card, data)
        if card:
            self.tid = card.field(1)
            self.x1 = card.field(2)
            self.x2 = card.field(3)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            self.x2 = data[2]
            fields = data[3:]
            isData = True
        self.parseFields(fields, nRepeated=2, isData=isData)

    def rawFields(self):
        fields = ['TABLED3', self.tid, self.x1, self.x2, None,
                  None, None, None, None] + self.table.fields() + ['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()


class TABLEM1(Table):
    type = 'TABLEM1'

    def __init__(self, card=None, data=None):
        Table.__init__(self, card, data)
        if card:
            self.tid = card.field(1)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            fields = data[1:]
            isData = True
        self.parseFields(fields, nRepeated=2, isData=isData)

    def rawFields(self):
        fields = ['TABLEM1', self.tid, None, None, None, None,
                  None, None, None] + self.table.fields() + ['ENDT']
        return fields


class TABLEM2(Table):
    type = 'TABLEM2'

    def __init__(self, card=None, data=None):
        Table.__init__(self, card, data)
        if card:
            self.tid = card.field(1)
            self.x1 = card.field(2)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            fields = data[2:]
            isData = True
        self.parseFields(fields, nRepeated=2, isData=isData)

    def rawFields(self):
        fields = ['TABLEM2', self.tid, self.x1, None, None, None,
                  None, None, None] + self.table.fields() + ['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()


class TABLEM3(Table):
    type = 'TABLEM3'

    def __init__(self, card=None, data=None):
        Table.__init__(self, card, data)
        if card:
            self.tid = card.field(1)
            self.x1 = card.field(2)
            self.x2 = card.field(3)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            self.x2 = data[2]
            fields = data[3:]
            isData = True
        self.parseFields(fields, nRepeated=2, isData=isData)

    def rawFields(self):
        fields = ['TABLEM3', self.tid, self.x1, self.x2, None,
                  None, None, None, None] + self.table.fields() + ['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()


class TABLEM4(Table):
    type = 'TABLEM4'

    def __init__(self, card=None, data=None):
        Table.__init__(self, card, data)
        if card:
            self.tid = card.field(1)
            self.x1 = card.field(2)
            self.x2 = card.field(3)
            self.x3 = card.field(4)
            self.x4 = card.field(5)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            self.x2 = data[2]
            self.x3 = data[3]
            self.x4 = data[4]
            fields = data[3:]
            isData = True
        self.parseFields(fields, nRepeated=1, isData=isData)

    def rawFields(self):
        fields = ['TABLEM4', self.tid, self.x1, self.x2, self.x3, self.x4,
                  None, None, None] + self.table.fields() + ['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()


class TABLES1(Table):
    type = 'TABLES1'

    def __init__(self, card=None, data=None):
        Table.__init__(self, card, data)
        if card:
            self.tid = card.field(1)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            fields = data[1:]
            isData = True
        self.parseFields(fields, nRepeated=2, isData=isData)

    def rawFields(self):
        fields = ['TABLES1', self.tid, None, None, None, None,
                  None, None, None] + self.table.fields() + ['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()


class TABLEST(Table):
    type = 'TABLEST'

    def __init__(self, card=None, data=None):
        Table.__init__(self, card, data)
        if card:
            self.tid = card.field(1)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            fields = data[1:]
            isData = True
        self.parseFields(fields, nRepeated=2, isData=isData)

    def rawFields(self):
        fields = ['TABLEST', self.tid, None, None, None, None,
                  None, None, None] + self.table.fields() + ['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()


class RandomTable(BaseCard):
    type = 'TABLE??'

    def __init__(self, card, data):
        pass


class TABRND1(RandomTable):
    type = 'TABRND1'

    def __init__(self, card=None, data=None):
        RandomTable.__init__(self, card, data)
        if card:
            self.tid = card.field(1)
            self.xaxis = card.field(2, 'LINEAR')
            self.yaxis = card.field(3, 'LINEAR')
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.xaxis = self.mapAxis(data[1])
            self.yaxis = self.mapAxis(data[2])
            fields = data[3:]
            isData = True
        assert self.xaxis in ['LINEAR', 'LOG'], 'xaxis=|%s|' % (self.xaxis)
        assert self.yaxis in ['LINEAR', 'LOG'], 'yaxis=|%s|' % (self.yaxis)
        self.parseFields(fields, nRepeated=2, isData=isData)

    def parseFields(self, fields, nRepeated, isData=False):
        self.table = TableObj(fields, nRepeated, isData)

    def mapAxis(self, axis):
        if axis == 0:
            axisType = 'LINEAR'
        else:
            raise ValueError('axis=|%s|' % (axis))
        return axisType

    def rawFields(self):
        fields = ['TABRND1', self.tid, self.xaxis, self.yaxis, None, None,
                  None, None, None] + self.table.fields() + ['ENDT']
        return fields

    def reprFields(self):
        xaxis = set_blank_if_default(self.xaxis, 'LINEAR')
        yaxis = set_blank_if_default(self.yaxis, 'LINEAR')
        fields = ['TABRND1', self.tid, xaxis, yaxis, None, None,
                  None, None, None] + self.table.fields() + ['ENDT']
        return fields


class TABRNDG(RandomTable):
    """
    Gust Power Spectral Density
    Defines the power spectral density (PSD) of a gust for aeroelastic response analysis.
    """
    type = 'TABRNDG'

    def __init__(self, card=None, data=None):
        RandomTable.__init__(self, card, data)
        if card:
            ## Table identification number. (Integer >0)
            self.tid = card.field(1)
            ## PSD Type: 1. von Karman; 2. Dryden
            self.Type = card.field(2)
            ## Scale of turbulence divided by velocity (units of time; Real)
            self.LU = card.field(3)
            ## Root-mean-square gust velocity. (Real)
            self.WG = card.field(4)
            assert self.Type in [1,
                                 2], 'Type must be 1 or 2.  Type=%s' % (self.Type)
        else:
            raise NotImplementedError()

    def rawFields(self):
        fields = ['TABRNDG', self.tid, self.Type, self.LU, self.WG]
        return fields

    def reprFields(self):
        return self.rawFields()


class TIC(Table):
    """Transient Initial Condition"""
    type = 'TIC'

    def __init__(self, card=None, data=None):
        """
        Defines values for the initial conditions of variables used in
        structural transient analysis. Both displacement and velocity values may
        be specified at independent degrees-of-freedom. This entry may not be
        used for heat transfer analysis.
        """
        Table.__init__(self, card, data)
        if card:
            self.sid = card.field(1)
            self.G = card.field(2)
            self.C = card.field(3)
            self.U0 = card.field(4)
            self.V0 = card.field(5)
        else:
            self.sid = data[0]
            self.G = data[1]
            self.C = data[2]
            self.U0 = data[3]
            self.V0 = data[4]

    def rawFields(self):
        fields = ['TIC', self.sid, self.G, self.C, self.U0, self.V0]
        return fields

    def reprFields(self):
        return self.rawFields()
