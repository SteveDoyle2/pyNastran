import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.elements      import CELAS1,CELAS2,CELAS3,CELAS4,CSHEAR,CONM2
from pyNastran.bdf.cards.elementsShell import CTRIA3,CQUAD4
from pyNastran.bdf.cards.elementsBars  import CROD,CBAR,CTUBE
from pyNastran.bdf.cards.elementsSolid import CTETRA4,CTETRA10,CPENTA6,CPENTA15,CHEXA8,CHEXA20
from pyNastran.bdf.cards.thermal       import CHBDYG,CHBDYP

class Geometry2(object):
    def readTable_Geom2(self):
        self.iTableMap = {
                           (2408, 24, 180):  self.readCBAR,   # record 8
                           (201,2,69):       self.readCDAMP1, # record 16 - not tested
                           (301,3,70):       self.readCDAMP2, # record 17 - not tested
                           (401,4,71):       self.readCDAMP3, # record 18 - not tested
                           (501,5,72):       self.readCDAMP4, # record 19 - not tested
                           (10608,106,404):  self.readCDAMP5, # record 20 - not tested
                           (601,6,73):       self.readCELAS1, # record 29
                           #(701,7,74):      self.readCELAS2, # record 30
                           (801,8,75):       self.readCELAS3, # record 31
                           (901,9,76):       self.readCELAS4, # record 32
                           (1501,15,64):     self.readCONM2,  # record 57
                           #(2958,51,177):    self.readCQUAD4, # record 69
                           #(13900,139,9989): self.readCQUAD4, # record 70
                           (3001,30,48):     self.readCROD,   # record 80
                           (12201,122,9013): self.readCTETP,  # record 86 - not done
                           #(5508,55,217):   self.readCTETRA, # record 87
                           #(5959,59,282):    self.readCTRIA3, # record 93
                           (3701,37,49):     self.readCTUBE,  # record 103
                           (5551,49,105):    self.readSPOINT, # record 118 - not done
                           (7308,73,253):    self.readCHEXA,  # record 45
                           (4108,41,280):    self.readCPENTA, # record 62
                           (10808,108,406):  self.readCHBDYG, # record 43
                           #(10908,109,407):  self.readCHBDYP, # record 44 - not done
                         }
        self.readRecordTable('GEOM2')

    def addOp2Element(self,elem):
        self.addElement(elem,allowOverwrites=True)


# 1-AEROQ4
# AEROT3
# BEAMAERO
# CAABSF
# CAXIF2
# CAXIF3
# CAXIF4
# CBAR
# CBARAO
# CBEAM
# CBEAMP
# CBEND
# CBUSH
# CBUSH1D
# CCONE

    def readCBAR(self,data):
        """
        CBAR(2408,24,180) - the marker for Record 8
        """
        #print "reading CBAR"
        while len(data)>=64: # 16*4
            eData = data[:64]
            data  = data[64:]
            f, = unpack('i',eData[28:32])
            #print "len(eData) = %s" %(len(eData))
            if   f==0:
                out = unpack('iiiifffiiiffffff',eData)
                (eid,pid,ga,gb,x1,x2,  x3,  f,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b) = out
                dataIn = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,x1,x2,x3]]
            elif f==1:
                out = unpack('iiiifffiiiffffff',eData)
                (eid,pid,ga,gb,x1,x2,  x3,  f,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b) = out
                dataIn = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,x1,x2,x3]]
            elif f==2:
                out = unpack('iiiiiiifiiffffff',eData)
                (eid,pid,ga,gb,g0,junk,junk,f,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b) = out
                dataIn = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,g0]]
            else:
                raise Exception('invalid f value...f=%s' %(f))
            ###
            elem = CBAR(None,dataIn)
            self.addOp2Element(elem)
        ###

    def readCDAMP1(self,data):
        """
        CDAMP1(201,2,69) - the marker for Record 16
        """
        print "reading CDAMP1"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            out = unpack('iiiiii',eData)
            (eid,pid,g1,g2,c1,c2) = out
            elem = CDAMP1(None,out)
            self.addOp2Element(elem)
        ###

    def readCDAMP2(self,data):
        """
        CDAMP2(301,3,70) - the marker for Record 17
        """
        print "reading CDAMP2"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            out = unpack('ifiiii',eData)
            (eid,b,g1,g2,c1,c2) = out
            elem = CDAMP2(None,out)
            self.addOp2Element(elem)
        ###

    def readCDAMP3(self,data):
        """
        CDAMP3(401,4,71) - the marker for Record 18
        """
        print "reading CDAMP3"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            out = unpack('iiii',eData)
            (eid,pid,s1,s2) = out
            elem = CDAMP3(None,out)
            self.addOp2Element(elem)
        ###

    def readCDAMP4(self,data):
        """
        CDAMP4(501,5,72) - the marker for Record 19
        """
        print "reading CDAMP4"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            out = unpack('ifii',eData)
            (eid,b,s1,s2) = out
            elem = CDAMP4(None,out)
            self.addOp2Element(elem)
        ###

    def readCDAMP5(self,data):
        """
        CDAMP5(10608,106,404) - the marker for Record 20
        """
        print "reading CDAMP5"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            out = unpack('iiii',eData)
            (eid,pid,s1,s2) = out
            elem = CDAMP5(None,out)
            self.addOp2Element(elem)
        ###

# CDUM2
# CDUM3
# CDUM4
# CDUM5
# CDUM6
# CDUM7
# CDUM8
# CDUM9

    def readCELAS1(self,data):
        """
        CELAS1(601,6,73) - the marker for Record 29
        """
        print "reading CELAS1"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            out = unpack('iiiiii',eData)
            (eid,pid,g1,g2,c1,c2) = out
            elem = CELAS1(None,out)
            self.addOp2Element(elem)
        ###

    def readCELAS2(self,data):
        """
        CELAS2(701,7,74) - the marker for Record 30
        """
        print "reading CELAS2"
        while len(data)>=32: # 8*4
            eData = data[:32]
            data  = data[32:]
            out = unpack('ifiiiiff',eData)
            (eid,k,g1,g2,c1,c2,ge,s) = out
            elem = CELAS2(None,out)
            self.addOp2Element(elem)
        ###

    def readCELAS3(self,data):
        """
        CELAS3(801,8,75) - the marker for Record 31
        """
        print "reading CELAS3"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            out = unpack('iiii',eData)
            (eid,pid,s1,s2) = out
            elem = CELAS3(None,out)
            self.addOp2Element(elem)
        ###

    def readCELAS4(self,data):
        """
        CELAS4(901,9,76) - the marker for Record 32
        """
        print "reading CELAS4"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            out = unpack('ifii',eData)
            (eid,k,s1,s2) = out
            elem = CELAS4(None,out)
            self.addOp2Element(elem)
        ###

    def readCONM2(self,data):
        """
        CONM2(1501,15,64) - the marker for Record 32
        """
        print "reading CONM2"
        while len(data)>=52: # 13*4
            eData = data[:52]
            data  = data[52:]
            out = unpack('iiiffffffffff',eData)
            (eid,g,cid,m,x1,x2,x3,i1,i2a,i2b,i3a,i3b,i3c) = out
            elem = CONM2(None,out)
            self.addOp2Element(elem)
        ###

# CFAST
# CFASTP
# CFLUID2
# CFLUID3
# CFLUID4
# CINT
# CGAP
# CHACAB
# CHACBR
# CHBDYE
# CHBDYG

    def readCHBDYG(self,data):
        """
        CHBDYG(10808,108,406) - the marker for Record 43
        """
        print "reading CHBDYG"
        while len(data)>=64: # 16*4
            eData = data[:64]
            data  = data[64:]
            (eid,blank,Type,iviewf,iviewb,radmidf,radmidb,blank2,
                      g1,g2,g3,g4,g5,g6,g7,g8) = unpack('iiiiiiiiiiiiiiii',eData)
            dataIn = [eid,Type,iviewf,iviewb,radmidf,radmidb,
                      g1,g2,g3,g4,g5,g6,g7,g8]
            elem = CHBDYG(None,dataIn)
            self.addOp2Element(elem)
        ###

# CHBDYP

    #def readCHBDYP(self,data):
    #    pass

    def readCHEXA(self,data):
        """
        CHEXA(7308,73,253) - the marker for Record 45
        """
        print "reading CHEXA"
        while len(data)>=88: # 22*4
            eData = data[:88]
            data  = data[88:]
            out = unpack('iiiiiiiiiiiiiiiiiiiiii',eData)
            (eid,pid,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,
             g11,g12,g13,g14,g15,g16,g17,g18,g19,g20) = out

            dataIn   = [eid,pid,g1,g2,g3,g4,g5,g6,g7,g8,]
            bigNodes = [g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20]
            if sum(bigNodes)>0:
                elem = CHEXA20(None,dataIn+bigNodes)
            else:
                elem = CHEXA8(None,dataIn)
            self.addOp2Element(elem)
        ###

# CHEXA20F
# CHEXAFD
# CHEXAL
# CHEXP
# CHEXPR
# CMASS1
# CMASS2
# CMASS3
# CMASS4
# CMFREE
# CONM1
# CONM2
# CONROD
# CONV
# CONVM
# CPENP

    def readCPENTA(self,data):
        """
        CPENTA(4108,41,280) - the marker for Record 62
        """
        print "reading CPENTA"
        while len(data)>=68: # 17*4
            eData = data[:68]
            data  = data[68:]
            out = unpack('iiiiiiiiiiiiiiiii',eData)
            (eid,pid,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,
             g11,g12,g13,g14,g15) = out
            
            dataIn   = [eid,pid,g1,g2,g3,g4,g5,g6]
            bigNodes = [g7,g8,g9,g10,g11,g12,g13,g14,g15]
            if sum(bigNodes)>0:
                elem = CPENTA15(None,dataIn+bigNodes)
            else:
                elem = CPENTA6(None,dataIn)
            self.addOp2Element(elem)
        ###

# CPENPR
# CPENT15F
# CPENT6FD
# CQDX4FD
# CQDX9FD
# CQUAD

    def readCQUAD4(self,data):
        """
        CQUAD4(2958,51,177)    - the marker for Record 69
        CQUAD4(13900,139,9989) - the marker for Record 70
        """
        #print "reading CQUAD4"
        while len(data)>=56: # 14*4
            eData = data[:56]
            data  = data[56:]
            (eid,pid,n1,n2,n3,n4,theta,zoffs,blank,tflag,t1,t2,t3,t4) = unpack('iiiiiiffiiffff',eData)
            #print "eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s blank=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" %(eid,pid,n1,n2,n3,n4,theta,zoffs,blank,tflag,t1,t2,t3,t4)
            dataInit = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            elem = CQUAD4(None,dataInit)
            self.addOp2Element(elem)
        ###

# CQUAD4FD
# CQUAD8
# CQUAD9FD
# CQUADP
# CQUADR
# CQUADX
# CRBAR
# CRBE1
# CRBE3
# CRJOINT

    def readCROD(self,data):
        """
        CROD(3001,30,48)    - the marker for Record 80
        """
        #print "reading CROD"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            out = unpack('iiii',eData)
            (eid,pid,n1,n2) = out
            elem = CROD(None,out)
            self.addOp2Element(elem)
        ###

# CRROD
# CSEAM
# CSHEAR
# CSLOT3
# CSLOT4

    def readCTETP(self,data):
        """
        CTETP(12201,122,9013)    - the marker for Record 86
        @todo create object
        """
        print "reading CTETP"
        raise Exception('needs work...')
        while len(data)>=108: # 27*4
            eData = data[:108]
            data  = data[108:]
            out = unpack('iiiiiiiiiiiiiiiiiiiiiiiiiii',eData)
            (eid,pid,n1,n2,n3,n4,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,
            f1,f2,f3,f4,b1,ee1,ee2,ee3,ee4) = out
            print "out = ",out
            e = [e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12]
            f = [f1,f2,f3,f4]
            ee = [ee1,ee2,ee3,ee4]
            
            print "e  = ",e
            print "f  = ",f
            print "b1  = ",b1
            print "ee = ",ee
            dataIn = [eid,pid,n1,n2,n2,n3,n4]
            elem = CTETRA4(None,dataIn)
            self.addOp2Element(elem)
        ###

    def readCTETRA(self,data):
        """
        CTETRA(5508,55,217)    - the marker for Record 87
        """
        print "reading CTETRA"
        while len(data)>=48: # 12*4
            eData = data[:48]
            data  = data[48:]
            out = unpack('iiiiiiiiiiii',eData)
            (eid,pid,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10) = out
            #print "out = ",out

            dataIn   = [eid,pid,n1,n2,n3,n4]
            bigNodes = [g5,g6,g7,g8,g9,g10]
            if sum(bigNodes)>0:
                elem = CTETRA10(None,dataIn+bigNodes)
            else:
                elem = CTETRA4(None,dataIn)
            self.addOp2Element(elem)
        ###

# CTETPR
# CTETR10F
# CTETR4FD
# CTQUAD
# CTTRIA

    def readCTRIA3(self,data):
        """
        (5959,59,282)    - the marker for Record 93
        """
        #print "reading CTRIA3"
        while len(data)>=52: # 13*4
            eData = data[:52]
            data  = data[52:]
            out = unpack('iiiiiffiiifff',eData)
            #print "eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s" %(eid,pid,n1,n2,n3,theta,zoffs,blank1,blank2,tflag,t1,t2,t3)
            (eid,pid,n1,n2,n3,theta,zoffs,blank1,blank2,tflag,t1,t2,t3) = out
            dataIn = [eid,pid,n1,n2,n3,theta,zoffs,tflag,t1,t2,t3]
            elem = CTRIA3(None,dataIn)
            self.addOp2Element(elem)
        ###

# CTRIAFD
# CTRIA6
# CTRIA6FD
# CTRIAP
# CTRIAR
# CTRIAX
# CTRIAX6
# CTRIX3FD
# CTRIX6FD

    def readCTUBE(self,data):
        """
        CTUBE(3701,37,49)    - the marker for Record 103
        """
        print "reading CTUBE"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            out = unpack('iiii',eData)
            (eid,pid,n1,n2) = out
            elem = CTUBE(None,out)
            self.addOp2Element(elem)
        ###

# CVISC
# CWELD
# CWELDC
# CWELDG
# CWSEAM
# GENEL
# GMDNDC
# GMBNDS
# GMINTC
# GMINTS
# PLOTEL
# RADBC
# RADINT
# SINT

    def readSPOINT(self,data):
        """
        (5551,49,105)    - the marker for Record 118
        @todo create object
        """
        print "reading SPOINT"
        while len(data)>=4: # 4*4
            eData = data[:4]
            data  = data[4:]
            (nid) = unpack('i',eData)
            #node = SPOINT(None,[nid])
            #self.addNode(node)
        ###

# VUBEAM
# VUHEXA
# VUQUAD4
# VUPENTA
# VUTETRA
# VUTRIA
# VUBEAM
# VUHEXA
# VUQUAD4
# WELDP
