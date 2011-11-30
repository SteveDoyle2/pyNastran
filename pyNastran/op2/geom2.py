import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.elementsShell import CTRIA3,CQUAD4
from pyNastran.bdf.cards.elementsBars  import CROD
from pyNastran.bdf.cards.elementsSolid import CTETRA4,CTETRA10,CPENTA15,CHEXA8,CHEXA20

class Geometry2(object):

    def readTable_Geom2(self):
        self.iTableMap = {
                           (201,2,69):       self.readCDamp1, # 16
                           (301,3,70):       self.readCDamp2, # 17
                           (401,4,71):       self.readCDamp3, # 18
                           (501,5,72):       self.readCDamp4, # 19
                           (10608,106,404):  self.readCDamp5, # 20
                           (601,6,73):       self.readCElas1, # 29
                           (701,7,74):       self.readCElas2, # 30
                           (801,8,75):       self.readCElas3, # 31
                           (901,9,76):       self.readCElas4, # 32
                           (2958,51,177):    self.readCQUAD4, # record 69
                           (13900,139,9989): self.readCQUAD4, # record 70
                           (3001,30,48):     self.readCROD,   # record 80
                          #(12201,122,9013): self.readCTETP,  # 86
                           (5508,55,217):    self.readCTETRA, # 87
                           (5959,59,282):    self.readCTRIA3, # record 93
                           (3701,37,49):     self.readCTUBE,  # record 103
                           (5551,49,105):    self.readSPOINT, # record 118
                         }
        self.readRecordTable('GEOM2')

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

    def readCDamp1(self,data):
        """
        (201,2,69) - the marker for Record 16
        @todo create object
        """
        print "reading CDAMP1"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            (eid,pid,g1,g2,c1,c2) = unpack('iiiiii',eData)
        ###

    def readCDamp2(self,data):
        """
        (301,3,70) - the marker for Record 17
        @todo create object
        """
        print "reading CDAMP2"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            (eid,b,g1,g2,c1,c2) = unpack('ifiiii',eData)
        ###

    def readCDamp3(self,data):
        """
        (401,4,71) - the marker for Record 18
        @todo create object
        """
        print "reading CDAMP3"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (eid,pid,s1,s2) = unpack('iiii',eData)
        ###

    def readCDamp4(self,data):
        """
        (501,5,72) - the marker for Record 19
        @todo create object
        """
        print "reading CDAMP4"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (eid,b,s1,s2) = unpack('ifii',eData)
        ###

    def readCDamp5(self,data):
        """
        (10608,106,404) - the marker for Record 20
        @todo create object
        """
        print "reading CDAMP5"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (eid,pid,s1,s2) = unpack('iiii',eData)
        ###

# CDUM2
# CDUM3
# CDUM4
# CDUM5
# CDUM6
# CDUM7
# CDUM8
# CDUM9

    def readCElas1(self,data):
        """
        (601,6,73) - the marker for Record 29
        @todo create object
        """
        print "reading CELAS1"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            (eid,pid,g1,g2,c1,c2) = unpack('iiiiii',eData)
        ###

    def readCElas2(self,data):
        """
        (701,7,74) - the marker for Record 30
        @todo create object
        """
        print "reading CELAS2"
        while len(data)>=32: # 8*4
            eData = data[:32]
            data  = data[32:]
            (eid,k,g1,g2,c1,c2,ge,s) = unpack('ifiiiiff',eData)
        ###

    def readCElas3(self,data):
        """
        (801,8,75) - the marker for Record 31
        @todo create object
        """
        print "reading CELAS3"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (eid,pid,s1,s2) = unpack('iiii',eData)
        ###

    def readCElas4(self,data):
        """
        (901,9,76) - the marker for Record 32
        @todo create object
        """
        print "reading CELAS4"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (eid,k,s1,s2) = unpack('ifii',eData)
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
# CHBDYP
# CHEXA
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
# CPENTA
# CPENPR
# CPENT15F
# CPENT6FD
# CQDX4FD
# CQDX9FD
# CQUAD

    def readCQUAD4(self,data):
        """
        (2958,51,177)    - the marker for Record 69
        (13900,139,9989) - the marker for Record 70
        """
        print "reading CQUAD4"
        while len(data)>=56: # 14*4
            eData = data[:56]
            data  = data[56:]
            (eid,pid,n1,n2,n3,n4,theta,zoffs,blank,tflag,t1,t2,t3,t4) = unpack('iiiiiiffiiffff',eData)
            #print "eid=%s pid=%s n1=%s n2=%s n3=%s n4=%s theta=%s zoffs=%s blank=%s tflag=%s t1=%s t2=%s t3=%s t4=%s" %(eid,pid,n1,n2,n3,n4,theta,zoffs,blank,tflag,t1,t2,t3,t4)
            if t1==-1.0: t1=1.0
            if t2==-1.0: t2=1.0
            if t3==-1.0: t3=1.0
            if t4==-1.0: t4=1.0
            dataInit = [eid,pid,n1,n2,n3,n4,theta,zoffs,tflag,t1,t2,t3,t4]
            elem = CQUAD4(None,dataInit)
            self.addElement(elem)
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
        (3001,30,48)    - the marker for Record 80
        @todo create object
        """
        print "reading CROD"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (eid,pid,n1,n2) = unpack('iiii',eData)
            dataInit = [eid,pid,n1,n2]
            elem = CROD(None,dataInit)
            self.addElement(elem)
        ###


# CRROD
# CSEAM
# CSHEAR
# CSLOT3
# CSLOT4

    def readCTETP(self,data):
        """
        (12201,122,9013)    - the marker for Record 86
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
            dataInit = [eid,pid,n1,n2,n2,n3,n4]
            elem = CTETRA4(None,dataInit)
            self.addElement(elem)
        ###

    def readCTETRA(self,data):
        """
        (5508,55,217)    - the marker for Record 87
        @todo create object
        """
        print "reading CTETRA"
        while len(data)>=48: # 12*4
            eData = data[:48]
            data  = data[48:]
            out = unpack('iiiiiiiiiiii',eData)
            (eid,pid,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10) = out
            #print "out = ",out
            n = [n1,n2,n3,n4,n5,n6,n7,n8,n9,n10]
            
            #print "n  = ",n
            dataInit = [eid,pid,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10]
            elem = CTETRA10(None,dataInit)
            self.addElement(elem)
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
        print "reading CTRIA3"
        while len(data)>=52: # 13*4
            eData = data[:52]
            data  = data[52:]
            (eid,pid,n1,n2,n3,theta,zoffs,blank1,blank2,tflag,t1,t2,t3) = unpack('iiiiiffiiifff',eData)
            #print "eid=%s pid=%s n1=%s n2=%s n3=%s theta=%s zoffs=%s blank1=%s blank2=%s tflag=%s t1=%s t2=%s t3=%s" %(eid,pid,n1,n2,n3,theta,zoffs,blank1,blank2,tflag,t1,t2,t3)
            if t1==-1.0: t1=1.0
            if t2==-1.0: t2=1.0
            if t3==-1.0: t3=1.0
            dataInit = [eid,pid,n1,n2,n3,theta,zoffs,tflag,t1,t2,t3]
            elem = CTRIA3(None,dataInit)
            self.addElement(elem)
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
        (3701,37,49)    - the marker for Record 103
        @todo create object
        """
        print "reading CTUBE"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (eid,pid,n1,n2) = unpack('iiii',eData)
            dataInit = [eid,pid,n1,n2]
            #elem = CTUBE(None,dataInit)
            self.addElement(elem)
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
