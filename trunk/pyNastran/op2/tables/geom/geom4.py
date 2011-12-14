import os
import sys
from struct import unpack

#from pyNastran.op2.op2Errors import *
#from pyNastran.bdf.cards.constraints import SPC,SPCADD

class Geometry4(object):

    def readTable_Geom4(self):
        self.iTableMap = {
                         (5561,76,215):   self.readASET,    # record 1  - not done
                         (5571,77,216):   self.readASET1,   # record 2  - not done
                         (10200,102,473): self.readBNDGRID, # record 3  - not done
                         (1510,15,328):   self.readCYAX,    # record 8  - not done
                         (5210,52,257):   self.readCYJOIN,  # record 9  - not done
                         (1710,17,330):   self.readCYSYM,   # record 11 - not done
                         (4901,49,17):    self.readMPC,     # record 16 - not done
                         (4891,60,83):    self.readMPCADD,  # record 17 - not done
                         (4951,63,92):    self.readOMIT1,   # record 19 - not done
                         (610, 6, 316):   self.readQSET1,   # record 21 - not done
                         (6601,66,292):   self.readRBAR,    # record 22 - not done
                         (6801,68,294):   self.readRBE1,    # record 23 - not done
                         (6901,69,295):   self.readRBE2,    # record 24 - not done
                         (7101,71,187):   self.readRBE3,    # record 25 - not done
                         (6501,65,291):   self.readRROD,    # record 30 - not done
                         (7001,70,186):   self.readRSPLINE, # record 31 - not done
                         (7201,72,398):   self.readRSSCON,  # record 32 - not done
                         (1210,12,322):   self.readSEQSET1, # record 40 - not done
                         (5501,55,16):    self.readSPC,     # record 44 - not done
                         (5481,58,12):    self.readSPC1,    # record 45 - not done
                         (5491,59,13):    self.readSPCADD,  # record 46 - not done
                         (5110,51,256):   self.readSPCD,    # record 47 - not done
                         (5601,56, 14):   self.readSUPORT,  # record 59 - not done
                         (10100,101,472): self.readSUPORT1, # record 60 - not done
                        #(4901,49,420017):self.readFake,    # record 
                        #(5561,76,0):     self.readFake,    # record 
                        #(5110,51,256):   self.readFake,    # record 
                        #(610,6,0):       self.readFake,    # record 
                        #(5501,55,620016):self.readFake,    # record 

                         }
        self.readRecordTable('GEOM4')

    def readASET(self,data):
        pass

    def readASET1(self,data):
        pass

    def readBNDGRID(self,data):
        pass

    def readCYAX(self,data):
        pass

    def readCYJOIN(self,data):
        pass

    def readCYSYM(self,data):
        pass

    def readMPC(self,data):
        pass

    def readMPCADD(self,data):
        pass

    def readOMIT1(self,data):
        pass

    def readQSET1(self,data):
        pass

    def readRBAR(self,data):
        pass

    def readRBE1(self,data):
        pass

    def readRBE2(self,data):
        pass

    def readRBE3(self,data):
        pass

    def readRROD(self,data):
        pass

    def readRSPLINE(self,data):
        pass

    def readRSSCON(self,data):
        pass

    def readSEQSET1(self,data):
        pass

    def readSPC(self,data):
        pass

    def readSPC1(self,data):
        pass

    def readSPCADD(self,data):
        pass

    def readSPCD(self,data):
        pass

    def readSUPORT(self,data):
        pass

    def readSUPORT1(self,data):
        pass

# ASET
# ASET1
# BNDGRID
# BSET
# BSET1
# CSET
# CSET1
# CYAX
# CYJOIN
# CYSUP
# CYSYM
# EGENDT
# GMBC
# GMSPC
# MPC
# MPCADD
# OMIT
# OMIT1
# QSET
# QSET1
# RBAR
# RBE1
# RBE2
# RBE3
# RBJOINT
# RBJSTIF
# RELEASE
# RPNOM
# RROD
# RSPLINE
# RSSCON
# RTRPLT
# RWELD
# SEBSET
# SEBSET1
# SECSET
# SECSET1
# SEQSET
# SEQSET1
# SESUP
# SEUSET
# SEUSET1
# SPC
# SPC1
# SPCADD
# SPCD
# SPCDE
# SPCDF
# SPCDG
# SPCE
# SPCEB
# SPCF
# SPCFB
# SPCGB
# SPCGRID
# SPCOFF
# SPCOFF1
# SUPORT
# SUPORT1
# TEMPBC
# USET
# USET1
