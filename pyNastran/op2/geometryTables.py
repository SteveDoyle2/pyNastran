import os
import sys
import struct
from struct import unpack
from op2Errors import *

class GeomObj(object):
    def __init__(self):
        pass
    def geomFunc(self,data):
        pass

class GeometryTables(object):

    def checkForNextTable(self):
        foundTable = False
        #print "---checking---"
        word = self.readTableName(rewind=True,debug=False,stopOnFailure=False)
        if word != None:
            foundTable = True
        #print '---checked---'
        #print "geomWord = ",word
        return foundTable

    def checkForNextSubTable(self,n):
        foundSubTable = False
        #print "tell = ",self.op2.tell()
        
        try:
            nOld = self.op2.tell()
            self.readMarkers([n,1,0])
            foundSubTable=True
            #print "subtable :) = ",foundSubTable
        except InvalidMarkersError:
            foundSubTable = False
        ###
        self.n = nOld
        self.op2.seek(self.n)

        return foundSubTable

    def readTable_Geom1(self):
        tableName = self.readTableName(rewind=False) # GEOM1
        self.tableInit(tableName)
        print "*tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        fields = self.readIntBlock()
        #print "fields = ",fields

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()

        iTable = -3
        while 1:  ## @todo could this cause an infinite loop...i dont this so...
            (tableName,isNextTable,isNextSubTable) = self.readGeomSubTable(iTable)
        
            if self.checkForNextTable():
                return
            iTable -= 1
        ###

    def readGeomSubTable(self,iTable):
        i=0
        isNextTable=False
        isNextSubTable=False
        self.readMarkers([iTable,1,0])

        tableName = self.readTableName(rewind=True,stopOnFailure=False)
        if tableName:
            print "**tableName = |%r|" %(tableName)
            return tableName,isNextTable,isNextSubTable

        while isNextSubTable==False and isNextTable==False:
            #print self.printSection(100)
            marker = self.getMarker()
            #print "marker = ",marker
            if marker<0:
                msg = 'marker is less than 0...'
                raise Exception(msg)
            data = self.readBlock()
            #self.op2Debug.write('ints = %s\n' %(str(ints)))

            isNextTable = self.checkForNextTable()
            isNextSubTable = self.checkForNextSubTable(iTable-1)
            print "i=%s tell=%s isNextTable=%s isNextSubTable=%s" %(i,self.op2.tell(),isNextTable,isNextSubTable)
            #if i==13:
            #    sys.exit('stopA')
            i+=1
        ### while
        
        #print "exiting the geom sub table"
        return (tableName,isNextTable,isNextSubTable)

    def readTable_Geom2(self):
        self.readTable_Geom1()
        return

        tableName = self.readTableName(rewind=False) # GEOM2/GEOM3/GEOM4/EPT/MPTS/DYNAMICS
        self.tableInit(tableName)
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-2,1,0]) #2
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()
        #print "word = |%r|" %(word)

        iTable = -3
        while 1:  ## @todo could this cause an infinite loop...i dont this so...
            (tableName,isNextTable,isNextSubTable) = self.readGeomSubTable(iTable)

            if self.checkForNextTable():
                return
            print "----end of SubTable=%s----" %(iTable)
            iTable -= 1
        ###

        
        self.printSection(100)
        sys.exit('end block...this should never happen...')


    def readTable_Geom3(self):
        self.readTable_Geom1()
        return

        ## GEOM3
        tableName = self.readTableName(rewind=False) # GEOM3
        self.tableInit(tableName)
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0])
        bufferWords = self.getMarker() # 2
        print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()
        print "word = |%r|" %(word)

        self.readMarkers([-3,1,0]) # 24
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords,bufferWords*4
        #self.printSection(4*187)
        ints = self.readIntBlock()
        print "ints = ",ints
        #self.skip(4*26)
        

        #self.printSection(4*30)
        #self.readMarkers([-4,1,0]) # 9
        #bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords,bufferWords*4
        #ints = self.readIntBlock()
        #print "4,-1,0,ints = ",ints
        #self.skip(4*11)


        data = self.readTableData([-4,1,0],'GEOM3')
        data = self.readTableData([-5,1,0],'GEOM3')
        print "time for block section 6..."
        
        self.readMarkers([-6,1,0])
        #assert self.op2.tell()==1488,self.op2.tell()
        

    def readTableData(self,markers,tableName):
        self.readMarkers(markers,tableName) # 3
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        self.printSection(80)
        data = self.readBlock()
        return data

    def readTable_Geom4(self):
        self.readTable_Geom1()
        return

        # GEOM4
        tableName = self.readTableName(rewind=False) # GEOM4
        self.tableInit(tableName)
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()
        print "word = |%r|" %(word)

        iTable = -3
        while 1:  ## @todo could this cause an infinite loop...i dont this so...
            (tableName,isNextTable,isNextSubTable) = self.readGeomSubTable(iTable)

            if self.checkForNextTable():
                return
            iTable -= 1
        ###
        sys.exit('end block of geom4...this should never happen...')

        print "------------"

    def readTable_EPT(self):
        self.readTable_Geom1()
        return

        tableName = self.readTableName(rewind=False) # EPT
        self.tableInit(tableName)
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        print "------------"
        word = self.readStringBlock()
        print "word = |%r|" %(word)

        iTable = -3
        while 1:  ## @todo could this cause an infinite loop...i dont this so...
            (tableName,isNextTable,isNextSubTable) = self.readGeomSubTable(iTable)

            if self.checkForNextTable():
                return
            iTable -= 1
        ###
        sys.exit('end block of EPT...this should never happen...')

    def readTable_MPTS(self):
        self.readTable_Geom1()
        return

        tableName = self.readTableName(rewind=False) # MPTS
        self.tableInit(tableName)
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords,bufferWords*4

        word = self.readStringBlock()
        print "word = |%r|" %(word)

        iTable = -3
        while 1:  ## @todo could this cause an infinite loop...i dont this so...
            (tableName,isNextTable,isNextSubTable) = self.readGeomSubTable(iTable)

            if self.checkForNextTable():
                return
            iTable -= 1
        ###
        sys.exit('end block of MPTS...this should never happen...')

    def readTable_DYNAMICS(self):
        self.readTable_Geom1()
        return

        tableName = self.readTableName(rewind=False) # DYNAMICS
        self.tableInit(tableName)
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords,bufferWords*4

        word = self.readStringBlock()
        print "word = |%r|" %(word)

        iTable = -3
        while 1:  ## @todo could this cause an infinite loop...i dont this so...
            (tableName,isNextTable,isNextSubTable) = self.readGeomSubTable(iTable)

            if self.checkForNextTable():
                return
            iTable -= 1
        ###
        sys.exit('end block of DYNAMICS...this should never happen...')
