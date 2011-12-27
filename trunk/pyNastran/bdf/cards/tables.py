#import sys

# my code
from baseCard import BaseCard
from pyNastran.general.general import ListPrint

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
                self.crashFields(fields,nFields)
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

    def crashFields(self,fields,nFields):
        try:
            msg = ''
            k = 0
            for i in range(nFields):
                for j in range(nRepeated):
                    msg += '%-8g ' %(fields[i*nRepeated+j])
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
        else:
            self.tid = data[0]
            self.x1  = data[1]
            fields   = data[2:]
        self.parseFields(fields,nRepeated=2,isData=isData)

    def reprFields(self):
        fields = ['TABLED2',self.tid,self.x1,None,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

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

    def reprFields(self):
        fields = ['TABLED3',self.tid,self.x1,self.x2,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

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

    def reprFields(self):
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

    def reprFields(self):
        fields = ['TABLEM2',self.tid,self.x1,None,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

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

    def reprFields(self):
        fields = ['TABLEM3',self.tid,self.x1,self.x2,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

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

    def reprFields(self):
        fields = ['TABLEM4',self.tid,self.x1,self.x2,self.x3,self.x4,None,None,None]+self.table.fields()+['ENDT']
        return fields

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

    def reprFields(self):
        fields = ['TABLES1',self.tid,None,None,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

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

    def reprFields(self):
        fields = ['TABLEST',self.tid,None,None,None,None,None,None,None]+self.table.fields()+['ENDT']
        return fields

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
