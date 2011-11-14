from fortranFile import FortranFile
from op2Codes import Op2Codes
import os
import sys
import struct
from struct import unpack

from geometryTables import GeometryTables
from elementsStressStrain import ElementsStressStrain
from ougv1 import OUGV1
from oqg1  import OQG1
from oes   import OES


class Op2(FortranFile,Op2Codes,GeometryTables,ElementsStressStrain,OQG1,OUGV1,OES):
    def __init__(self,infileName): 
        self.infilename = infileName
        #self.tablesToRead = ['GEOM1','GEOM2','GEOM3','GEOM4','OQG1','OUGV1','OES1X1']  # 'OUGV1','GEOM1','GEOM2'
        #self.tablesToRead = ['GEOM1','GEOM2','OQG1','OUGV1','OES1X1']  # 'OUGV1','GEOM1','GEOM2'
        self.tablesToRead = ['GEOM1','GEOM2','GEOM3','OQG1',]  # 'OUGV1','GEOM1','GEOM2'
        self.tablesToRead = ['OUGV1',]  # 'OUGV1','GEOM1','GEOM2'
        ## GEOM1 & GEOM2 are skippable on simple problems...hmmm
    
    def readTapeCode(self):
        self.printSection(500)
        #sys.exit('op2-readTapeCode')
        self.readMarkers([3])
        ints = self.readIntBlock()
        #print "*ints = ",ints
        self.readMarkers([7])

        word = self.readStringBlock() # Nastran Fort Tape ID Code - 
        #print "word = |%r|" %(word)

        self.readMarkers([2])
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-1])

        #data = self.getData(60)
        #self.printBlock(data)

    def read(self):
        self.op2 = open(self.infilename,'rb')
        
        self.n = self.op2.tell()
        print "self.n = ",self.n
        self.readTapeCode()
        #sys.exit('end of tape code')

        isAnotherTable = True
        while isAnotherTable:
            tableName = self.readTableName(rewind=True)
            print "tableName = |%r|" %(tableName)
 
            if tableName in self.tablesToRead:
                if tableName=='GEOM1': # nodes,coords,etc.
                    self.readTable_Geom1()
                elif tableName=='GEOM2': # elements
                    self.readTable_Geom2()
                    print self.printSection(80)
                elif tableName=='GEOM3': # static/thermal loads
                    #print "**************"
                    self.readTable_Geom3()
                    print "**************"
                    #sys.exit('end of geom3')
                elif tableName=='GEOM4': # constraints
                    self.readTable_Geom4()
                    print "**************"

                elif tableName=='EPT':   # element properties
                    self.readTable_EPT()
                    print "**************"
                elif tableName=='MPTS':  # material properties
                    self.readTable_MPTS()
                    print "**************"


                elif tableName=='OEF1X':  # ???
                    self.readTable_OEF1X()
                    print "**************"
                elif tableName=='OQG1':  # spc forces
                    self.readTable_OQG1()
                    print "**************"
                elif tableName=='OUGV1': # displacements/velocity/acceleration
                    self.readTable_OUGV1()
                    print "**************"
                elif tableName=='OES1X1': # stress
                    self.readTable_OES1X1()
                    print "**************"
                    sys.exit('stopping after oes1x1')
                else:
                    raise Exception('unhandled tableName=|%s|' %(tableName))
                #print "---isAnotherTable---"
                #(isAnotherTable) = self.hasMoreTables()
                isAnotherTable = True
                self.printSection(100)
            else:
                (isAnotherTable) = self.skipNextTable()
                continue
            print "*** finished tableName = |%r|" %(tableName)
            ###
        ###
        #

        #tableName = self.readTableName(rewind=True)
        #self.readTable_Geom2()
        #self.readTable_Geom3()
        #self.readTable_Geom4()
        #self.readTable_OQG1()
        #self.readTable_OES1X1()
        print "---end of all tables---"

        #self.printSection(4*51+12)
        
    def parseAnalysisCode(self,data):
        #self.printBlock(data)
        (aCode,tCode,elementType,iSubcase) = unpack('iiii',data[:16])
        deviceCode   = aCode%10
        approachCode = (aCode-deviceCode)/10
        print "aCode=%s analysisCode=%s deviceCode=%s tCode=%s elementType=%s iSubcase=%s" %(aCode,approachCode,deviceCode,tCode,elementType,iSubcase)
        return (approachCode,deviceCode,tCode,elementType,iSubcase)
