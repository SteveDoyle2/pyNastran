import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.materials import CREEP,MAT1,MAT2,MAT3,MAT4,MAT5,MAT8,MAT9,MAT10

class MPT(object):

    def readTable_MPTS(self):
        self.bigMaterials = {}
        self.iTableMap = {
                         (1003,10,245): self.readCreep, # record 1
                         (103,1,77):    self.readMat1,  # record 2
                         (203,2,78):    self.readMat2,  # record 3
                         (1403,14,122): self.readMat3,  # record 4
                         (2103,21,234): self.readMat4,  # record 5
                         (2203,22,235): self.readMat5,  # record 6
                         (2503,25,288): self.readMat8,  # record 7
                         (2801,28,365): self.readMat10, # record 9
                         #(3003,30,286): self.read
                         #(3103,31,337): self.read
                         }
        self.readRecordTable('MPTS')

    def readCreep(self,data):
        """
        CREEP(1003,10,245) - record 1
        """
        print "reading CREEP"
        while len(data)>=64: # 16*4
            eData = data[:64]
            data  = data[64:]
            out = unpack('iffiiiififffffff',eData)
            (mid,T0,exp,form,tidkp,tidcp,tidcs,thresh,Type,ag1,ag2,ag3,ag4,ag5,ag6,ag7) = out
            mat = CREEP(None,out)
            self.addMaterial(mat)
        ###

    def readMat1(self,data):
        """
        MAT1(103,1,77) - record 2
        """
        print "reading MAT1"
        while len(data)>=48: # 12*4
            eData = data[:48]
            data  = data[48:]
            out = unpack('iffffffffffi',eData)
            (mid,E,G,nu,rho,A,TRef,ge,St,Sc,Ss,mcsid) = out
            mat = MAT1(None,out)
            self.addMaterial(mat)
        ###

    def readMat2(self,data):
        """
        MAT2(203,2,78) - record 3
        """
        print "reading MAT2"
        while len(data)>=68: # 17*4
            eData = data[:68]
            data  = data[68:]
            out = unpack('ifffffffffffffffi',eData)
            (mid,g1,g2,g3,g4,g5,g6,rho,aj1,aj2,aj3,TRef,ge,St,Sc,Ss,mcsid) = out
            #print "MAT2 = ",out
            mat = MAT2(None,out)

            if mid>1e8:
                self.bigMaterials[mid] = mat
            else:
                self.addMaterial(mat)
            ###
        ###

    def readMat3(self,data):
        """
        MAT3(1403,14,122) - record 4
        """
        print "reading MAT3"
        while len(data)>=64: # 16*4
            eData = data[:64]
            data  = data[64:]
            out = unpack('iffffffffifffffi',eData)
            (mid,ex,eth,ez,nuxth,nuthz,nuzx,rho,gzx,blank,ax,ath,az,TRef,ge,blank) = out
            mat = MAT3(None,[mid,ex,eth,ez,nuxth,nuthz,nuzx,rho,gzx,ax,ath,az,TRef,ge])
            self.addMaterial(mat)
        ###

    def readMat4(self,data):
        """
        MAT4(2103,21,234) - record 5
        """
        print "reading MAT4"
        while len(data)>=44: # 11*4
            eData = data[:44]
            data  = data[44:]
            out = unpack('iffffffffff',eData)
            (mid,k,cp,rho,h,mu,hgen,refenth,tch,tdelta,qlat) = out
            mat = MAT4(None,out)
            self.addMaterial(mat)
        ###

    def readMat5(self,data):
        """
        MAT5(2203,22,235) - record 6
        """
        print "reading MAT5"
        while len(data)>=40: # 10*4
            eData = data[:40]
            data  = data[40:]
            out = unpack('ifffffffff',eData)
            (mid,k1,k2,k3,k4,k5,k6,cp,rho,hgen) = out
            mat = MAT5(None,out)
            self.addMaterial(mat)
        ###

    def readMat8(self,data):
        """
        MAT8(2503,25,288) - record 7
        """
        print "reading MAT8"
        while len(data)>=76: # 19*4
            eData = data[:76]
            data  = data[76:]
            out = unpack('iffffffffffffffffff',eData)
            (mid,E1,E2,nu12,G12,G1z,G2z,rho,a1,a2,TRef,Xt,Xc,Yt,Yc,S,ge,f12,strn) = out
            mat = MAT8(None,out)
            self.addMaterial(mat)
        ###

# MAT9

    def readMat10(self,data):
        """
        MAT10(2801,28,365) - record 9
        """
        print "reading MAT10"
        while len(data)>=20: # 5*4
            eData = data[:20]
            data  = data[20:]
            out = unpack('iffff',eData)
            (mid,bulk,rho,c,ge) = out
            mat = MAT10(None,out)
            self.addMaterial(mat)
        ###

# MAT11
# MATHP
# MATS1
# MATT1
# MATT2
# MATT3
# MATT4
# MATT5
# MATT8
# MATT9
# MBOLT
# MBOLTUS
# MSTACK
# NLAUTO
# RADBND
# RADM
# RADMT
# NLPARM
# NLPCI
# TSTEPNL
