## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
#import sys

# my code
from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.general.general import ListPrint

class TIC(BaseCard): # Transient Initial Condition
    type = 'TIC'
    def __init__(self,card=None,data=None):
        """
        Defines values for the initial conditions of variables used in structural transient
        analysis. Both displacement and velocity values may be specified at independent
        degrees-of-freedom. This entry may not be used for heat transfer analysis.
        """
        Table.__init__(self,card,data)
        if card:
            self.sid = card.field(1)
            self.G   = card.field(2)
            self.C   = card.field(3)
            self.U0  = card.field(4)
            self.V0  = card.field(5)
        else:
            self.sid = data[0]
            self.G   = data[1]
            self.C   = data[2]
            self.U0  = data[3]
            self.V0  = data[4]
        ###

    def rawFields(self):
        fields = ['TIC',self.sid,self.G,self.C,self.U0,self.V0]
        return fields

    def reprFields(self):
        return self.rawFields()
   
class Table(BaseCard):
    type = 'TABLE??'
    def __init__(self,card,data):
        pass

    def mapAxis(self,axis):
        if axis==0:
            axisType = 'LINEAR'
        else:
            raise Exception('axis=|%s|' %(axis))
        return axisType

    def parseFields(self,fields,nRepeated,isData=False):
        self.table = TableObj(fields,nRepeated,isData)

class TableObj(object):
    def __init__(self,fields,nRepeated,isData=False):
        self.table = []
        fields = self.cleanupFields(fields,isData)
        
        nFields = len(fields)

        if not isData:
            if nFields%nRepeated!=0:
                self.crashFields(fields,nRepeated,nFields)
            assert nFields%nRepeated==0,'invalid table length nRepeat=%s fields=%s' %(nRepeated,ListPrint(fields))
            assert nFields%nRepeated==0,msg
        
        i=0
        while i<nFields:
            pack = []
            for j in range(nRepeated):
                pack.append(fields[i+j])
            ###
            i+=nRepeated
            self.table.append(pack)
        ###

    def crashFields(self,fields,nRepeated,nFields):
        try:
            msg = ''
            k = 0
            for i in range(nFields):
                for j in range(nRepeated):
                    try:
                        msg += '%-8g ' %(fields[i*nRepeated+j])
                    except TypeError:
                        msg += '*%-8s ' %(fields[i*nRepeated+j])
                    except IndexError:
                        msg += 'IndexError'
                    k+=1
                ###
                msg += '\n'
            ###
        except:
            print fields
            print msg
            #assert nFields%nRepeated==0,msg
            raise
        ###
    ###

    def cleanupFields(self,fields,isData=False):
        fields2 = [] # remove extra ENDTs
        foundENDT = False
        for field in fields:
            if isinstance(field,str) and 'ENDT' in field.upper():
                foundENDT = True
            else:
                fields2.append(field)
            ###
        ###
        if not isData:
            assert foundENDT==True,fields
        return fields2

    def fields(self):
        fields = []
        for pack in self.table:
            fields += pack
        return fields

class TABLED1(Table):
    type = 'TABLED1'
    def __init__(self,card=None,data=None):
        Table.__init__(self,card,data)
        if card:
            self.tid   = card.field(1)
            self.xaxis = card.field(2,'LINEAR')
            self.yaxis = card.field(3,'LINEAR')
            fields = card.fields(9)
            isData = False
        else:
            self.tid   = data[0]
            self.xaxis = self.mapAxis(data[1])
            self.yaxis = self.mapAxis(data[2])
            fields = data[3:]
            isData = True
        assert self.xaxis in ['LINEAR','LOG'],'xaxis=|%s|' %(self.xaxis)
        assert self.yaxis in ['LINEAR','LOG'],'yaxis=|%s|' %(self.yaxis)
        self.parseFields(fields,nRepeated=2,isData=isData)

    def rawFields(self):
        fields = ['TABLED1',self.tid,self.xaxis,self.yaxis,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

    def reprFields(self):
        xaxis = self.setBlankIfDefault(self.xaxis,'LINEAR')
        yaxis = self.setBlankIfDefault(self.yaxis,'LINEAR')
        fields = ['TABLED1',self.tid,xaxis,yaxis,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

class TABLED2(Table):
    type = 'TABLED2'
    def __init__(self,card=None,data=None):
        (footer) = Table.__init__(self,card,data)
        if card:
            self.tid = card.field(1)
            self.x1  = card.field(2)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.x1  = data[1]
            fields   = data[2:]
            isData = True
        self.parseFields(fields,nRepeated=2,isData=isData)

    def rawFields(self):
        fields = ['TABLED2',self.tid,self.x1,None,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()

class TABLED3(Table):
    type = 'TABLED3'
    def __init__(self,card=None,data=None):
        (footer) = Table.__init__(self,card,data)
        if card:
            self.tid = card.field(1)
            self.x1  = card.field(2)
            self.x2  = card.field(3)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.x1  = data[1]
            self.x2  = data[2]
            fields   = data[3:]
            isData = True
        self.parseFields(fields,nRepeated=2,isData=isData)

    def rawFields(self):
        fields = ['TABLED3',self.tid,self.x1,self.x2,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()

class TABLEM1(Table):
    type = 'TABLEM1'
    def __init__(self,card=None,data=None):
        (footer) = Table.__init__(self,card,data)
        if card:
            self.tid = card.field(1)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            fields   = data[1:]
            isDadta = True
        self.parseFields(fields,nRepeated=2,isData=isData)

    def rawFields(self):
        fields = ['TABLEM1',self.tid,None,None,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

class TABLEM2(Table):
    type = 'TABLEM2'
    def __init__(self,card=None,data=None):
        (footer) = Table.__init__(self,card,data)
        if card:
            self.tid = card.field(1)
            self.x1  = card.field(2)
            fields   = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.x1  = data[1]
            fields   = data[2:]
            isData = True
        self.parseFields(fields,nRepeated=2)

    def rawFields(self):
        fields = ['TABLEM2',self.tid,self.x1,None,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields
    
    def reprFields(self):
        return self.rawFields()


class TABLEM3(Table):
    type = 'TABLEM3'
    def __init__(self,card=None,data=None):
        (footer) = Table.__init__(self,card,data)
        if card:
            self.tid = card.field(1)
            self.x1  = card.field(2)
            self.x2  = card.field(3)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.x1  = data[1]
            self.x2  = data[2]
            fields   = data[3:]
            isData = True
        self.parseFields(fields,nRepeated=2,isData=isData)

    def rawFields(self):
        fields = ['TABLEM3',self.tid,self.x1,self.x2,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()

class TABLEM4(Table):
    type = 'TABLEM4'
    def __init__(self,card=None,data=None):
        (footer) = Table.__init__(self,card,data)
        if card:
            self.tid = card.field(1)
            self.x1  = card.field(2)
            self.x2  = card.field(3)
            self.x3  = card.field(4)
            self.x4  = card.field(5)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            self.x1  = data[1]
            self.x2  = data[2]
            self.x3  = data[3]
            self.x4  = data[4]
            fields   = data[3:]
            isData = True
        self.parseFields(fields,nRepeated=1,isData=isData)

    def rawFields(self):
        fields = ['TABLEM4',self.tid,self.x1,self.x2,self.x3,self.x4,None,None,None]+self.table.fields()+['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()

class TABLES1(Table):
    type = 'TABLES1'
    def __init__(self,card=None,data=None):
        (footer) = Table.__init__(self,card,data)
        if card:
            self.tid = card.field(1)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            fields   = data[1:]
            isData = True
        self.parseFields(fields,nRepeated=2,isData=isData)

    def rawFields(self):
        fields = ['TABLES1',self.tid,None,None,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()

class TABLEST(Table):
    type = 'TABLEST'
    def __init__(self,card=None,data=None):
        (footer) = Table.__init__(self,card,data)
        if card:
            self.tid = card.field(1)
            fields = card.fields(9)
            isData = False
        else:
            self.tid = data[0]
            fields   = data[1:]
            isData = True
        self.parseFields(fields,nRepeated=2,isData=isData)

    def rawFields(self):
        fields = ['TABLEST',self.tid,None,None,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

    def reprFields(self):
        return self.rawFields()

class TABRND1(Table):
    type = 'TABRND1'
    def __init__(self,card=None,data=None):
        (footer) = Table.__init__(self,card,data)
        if card:
            self.tid   = card.field(1)
            self.xaxis = card.field(2,'LINEAR')
            self.yaxis = card.field(3,'LINEAR')
            fields = card.fields(9)
            isData = False
        else:
            self.tid   = data[0]
            self.xaxis = self.mapAxis(data[1])
            self.yaxis = self.mapAxis(data[2])
            fields = data[3:]
            isData = True
        assert self.xaxis in ['LINEAR','LOG'],'xaxis=|%s|' %(self.xaxis)
        assert self.yaxis in ['LINEAR','LOG'],'yaxis=|%s|' %(self.yaxis)
        self.parseFields(fields,nRepeated=2,isData=isData)

    def rawFields(self):
        fields = ['TABRND1',self.tid,self.xaxis,self.yaxis,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

    def reprFields(self):
        xaxis = self.setBlankIfDefault(self.xaxis,'LINEAR')
        yaxis = self.setBlankIfDefault(self.yaxis,'LINEAR')
        fields = ['TABRND1',self.tid,xaxis,yaxis,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields
