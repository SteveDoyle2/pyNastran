from fortranFile import FortranFile
from op2Codes import Op2Codes
from op2Errors import *
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
        #self.tablesToRead = ['GEOM1','GEOM2','GEOM3','OQG1','OUGV1','OES1X1']  # 'OUGV1','GEOM1','GEOM2'
        self.tablesToRead = ['OQG1','OUGV1','OES1X1','OSTR1X']  # 'OUGV1','GEOM1','GEOM2'
        #self.tablesToRead = ['OUGV1',]  # 'OUGV1','GEOM1','GEOM2'
        ## GEOM1 & GEOM2 are skippable on simple problems...hmmm

        self.displacements = {}
        self.temperatures  = {}

        self.forces = {}
        self.fluxes = {}

        self.rodStress   = {}
        self.rodStrain   = {}
        self.barStress   = {}
        self.barStrain   = {}
        self.plateStress = {}
        self.plateStrain = {}
        self.solidStress = {}
        self.solidStrain = {}

    def printResults(self):
        results = [self.displacements,self.temperatures,
                   self.forces,self.fluxes,
                   self.rodStress,self.rodStrain,
                   self.barStress,self.barStrain,
                   self.plateStress,self.plateStrain,
                   self.solidStress,self.solidStrain]
        
        msg = '---ALL RESULTS---\n'
        for result in results:
            for iSubcase,res in sorted(result.items()):
                msg += 'iSubcase = %s\n' %(iSubcase)
                msg += str(res) + '\n'
            ###
        ###
        return msg
        
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
            print '-'*80
            try:
                tableName = self.readTableName(rewind=True)
            except EndOfFileError:  # the isAnotherTable method sucks...
                isAnotherTable = False
                print "***ok exit, but it could be better..."
                break
            except AssertionError:  # the isAnotherTable method sucks...
                isAnotherTable = False
                print "***poor exit, but it worked..."
                break
            print "tableName = |%r|" %(tableName)
 
            if tableName in self.tablesToRead:
                if tableName=='GEOM1': # nodes,coords,etc.
                    self.readTable_Geom1()
                elif tableName=='GEOM2': # elements
                    self.readTable_Geom2()
                    print self.printSection(80)
                elif tableName=='GEOM3': # static/thermal loads
                    self.readTable_Geom3()
                    #sys.exit('end of geom3')
                elif tableName=='GEOM4': # constraints
                    self.readTable_Geom4()

                elif tableName=='EPT':   # element properties
                    self.readTable_EPT()
                elif tableName=='MPTS':  # material properties
                    self.readTable_MPTS()


                elif tableName=='OEF1X':  # ???
                    self.readTable_OEF1X()
                elif tableName=='OQG1':  # spc forces
                    self.readTable_OQG1()
                elif tableName=='OUGV1': # displacements/velocity/acceleration
                    self.readTable_OUGV1()
                elif tableName=='OES1X1': # stress
                    self.readTable_OES1X1()
                elif tableName=='OSTR1X': # strain
                    self.readTable_OES1X1()
                else:
                    raise Exception('unhandled tableName=|%s|' %(tableName))
                #print "---isAnotherTable---"
                #(isAnotherTable) = self.hasMoreTables()
                isAnotherTable = True
                #self.printSection(100)
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
        
    def parseApproachCode(self,data):
        #self.printBlock(data)
        (aCode,tCode,elementType,iSubcase) = unpack('iiii',data[:16])
        self.deviceCode   = aCode%10
        self.approachCode = (aCode-self.deviceCode)/10
        print "aCode(1)=%s analysisCode=%s deviceCode=%s tCode(2)=%s elementType(3)=%s iSubcase(4)=%s" %(aCode,self.approachCode,self.deviceCode,tCode,elementType,iSubcase)
        print "tableType = ",self.printTableCode(tCode)
        return (tCode,elementType,iSubcase)
