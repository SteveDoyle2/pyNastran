import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
#from pyNastran.bdf.cards.nodes import GRID
#from pyNastran.bdf.cards.coordinateSystems import CORD1R,CORD2R,CORD2C,CORD3G #CORD1C,CORD1S,CORD2S


class LAMA(object):

    def readTable_LAMA(self):
        #self.tableName = 'DESTAB'
        tableName = self.readTableName(rewind=False) # DESTAB
        self.tableInit(tableName)
        print "tableName1 = |%r|" %(tableName)
        print "tableName2 = |%r|" %(self.tableName)

        self.readMarkers([-1,7],'LAMA')
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0],'LAMA')
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords
        
        
        word = self.readStringBlock() # DESTAB
        print "word = |%s|" %(word)
        self.readMarkers([-3,1,0],'LAMA')

        #data = self.getData(4*50)
        #print self.printBlock(data)

        self.readTable_LAMA_3(-3)

        self.readMarkers([-4,1,0],'LAMA')
        self.readTable_LAMA_4(-4)

        self.readMarkers([-5,1,0],'LAMA')
        #data = self.getData(4*30)
        #print self.printBlock(data)
        return
        sys.exit('stopping in LAMA')
        if 0:
            iTable = -3
            #imax   = -244

            while bufferWords: # read until bufferWords=0
                self.readMarkers([iTable,1,0],'LAMA')
                nOld = self.n
                bufferWords = self.getMarker()
                #print "bufferWords = ",bufferWords
                if bufferWords==0: # maybe read new buffer...
                    self.goto(nOld)
                    break
                data = self.readBlock()
                #print "len(data) = ",len(data)
                self.readDesvar(data)
                iTable -= 1
            self.printSection(80)
        ###
        #self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "1-bufferWords = ",bufferWords,bufferWords*4

        #print self.printSection(300)
        #sys.exit('asdf')

    def readTable_LAMA_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        #print self.printBlock(data)
        
        (three) = self.parseApproachCode(data)

        self.addDataParameter(data,'seven',  'i',10,False)   ## seven
        self.addDataParameter(data,'resFlag','i',11,False)   ## residual vector augmentation flag
        self.addDataParameter(data,'fldFlag','i',12,False)   ## fluid modes Flag
        
        print self.dataCode
        #self.addDataParameter(data,'formatCode',  'i',9,False)   ## format code
        #self.addDataParameter(data,'numWide',     'i',10,False)  ## number of words per entry in record; @note is this needed for this table ???
        
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)
        
        #print "*iSubcase=%s"%(self.iSubcase)
        #print "analysisCode=%s tableCode=%s thermal=%s" %(self.analysisCode,self.tableCode,self.thermal)

        #self.printBlock(data)
        self.readTitle()

    def readTable_LAMA_4(self,iTable): # iTable=-4
        bufferWords = self.getMarker() # 70*4=280
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4) # dummy - 70*4=280
        print self.printBlock(data)
        #print "280/3 = ",280/4
        nModes = bufferWords//7

        msg = ['                                              R E A L   E I G E N V A L U E S\n',
               '   MODE    EXTRACTION      EIGENVALUE            RADIANS             CYCLES            GENERALIZED         GENERALIZED\n',
               '    NO.       ORDER                                                                       MASS              STIFFNESS\n']
        
        for i in range(nModes):
            data = self.getData(28) # 4*7
            (mode,order,eigen,omega,freq,mass,stiff) = unpack('iifffff',data)
            print "mode=%s order=%s eigen=%s omega=%s freq=%s mass=%s stiff=%s" %(mode,order,eigen,omega,freq,mass,stiff)


            msg.append(' %8s  %8s       %13s       %13s       %13s       %13s       %13s\n'%(mode,order,eigen,omega,freq,mass,stiff))
        #print ''.join(msg)
        print "self.iSubcase = ",self.iSubcase
#                       '        1         1        8.232776E+06        2.869281E+03        4.566603E+02        8.719168E-03        7.178296E+04
#                       '        2         2        8.232776E+06        2.869281E+03        4.566603E+02        8.719168E-03        7.178296E+04

        data = self.getData(4)
