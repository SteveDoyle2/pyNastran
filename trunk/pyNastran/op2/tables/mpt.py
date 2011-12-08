import os
import sys
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.materials import CREEP,MAT1,MAT2,MAT3,MAT4,MAT5,MAT8,MAT9,MAT10

class MPT(object):

    def readTable_MPTS(self):
        self.bigMaterials = {}
        self.iTableMap = {
                         (1003,10,245): self.readCREEP, # record 1
                         (103,1,77):    self.readMAT1,  # record 2
                         (203,2,78):    self.readMAT2,  # record 3
                         (1403,14,122): self.readMAT3,  # record 4
                         (2103,21,234): self.readMAT4,  # record 5
                         (2203,22,235): self.readMAT5,  # record 6
                         (2503,25,288): self.readMAT8,  # record 7
                         (2603,26,300): self.readMAT9,  # record 8 - not tested
                         (2801,28,365): self.readMAT10, # record 9
                         (503,5,90):    self.readMATS1, # record 12 - not done

                         (3003,30,286): self.readNLPARM,  # record 27 - not done
                         (3103,31,337): self.readTSTEPNL, # record 29 - not done
                         (503,  5, 90): self.readMATS1,   # record 12 - not done
                         (703,  7, 91): self.readMATT1,   # record 13 - not done
                         (803,  8,102): self.readMATT2,   # record 14 - not done
                         (1503,14,189): self.readMATT3,   # record 15 - not done
                         (2303,23,237): self.readMATT4,   # record 16 - not done
                         (2403,24,238): self.readMATT5,   # record 17 - not done
                         (2703,27,301): self.readMATT9,   # record 19 - not done
                         (8802,88,413): self.readRADM,    # record 25 - not done

                         }

        self.readRecordTable('MPTS')

    def addOp2Material(self,mat):
        self.addMaterial(mat,allowOverwrites=True)

    def readCREEP(self,data):
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
            self.addOp2Material(mat)
        ###

    def readMAT1(self,data):
        """
        MAT1(103,1,77) - record 2
        """
        #print "reading MAT1"
        while len(data)>=48: # 12*4
            eData = data[:48]
            data  = data[48:]
            out = unpack('iffffffffffi',eData)
            (mid,E,G,nu,rho,A,TRef,ge,St,Sc,Ss,mcsid) = out
            mat = MAT1(None,out)
            self.addOp2Material(mat)
        ###

    def readMAT2(self,data):
        """
        MAT2(203,2,78) - record 3
        """
        #print "reading MAT2"
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
                self.addOp2Material(mat)
            ###
        ###

    def readMAT3(self,data):
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
            self.addOp2Material(mat)
        ###

    def readMAT4(self,data):
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
            self.addOp2Material(mat)
        ###

    def readMAT5(self,data):
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
            self.addOp2Material(mat)
        ###

    def readMAT8(self,data):
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
            self.addOp2Material(mat)
        ###

    def readMAT9(self,data):
        """
        MAT9(2603,26,300) - record 9
        @todo buggy
        """
        print "reading MAT9"
        while len(data)>=140: # 35*4
            eData = data[:140]
            data  = data[140:]
            out = unpack('iiiiiiiiiiiiiiiiiiiiiifffffffffiiii',eData)
            
            (mid,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,
            rho,a1,a2,a3,a4,a5,a6,TRef,ge,blank1,blank2,blank3,blank4) = out
            dataIn = [mid,[g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21],
                      rho,[a1,a2,a3,a4,a5,a6],
                      TRef,ge]
            print "dataIn = ",dataIn
            #mat = MAT9(None,dataIn)
            #self.addOp2Material(mat)
        ###

    def readMAT10(self,data):
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
            self.addOp2Material(mat)
        ###

# MAT11 - unused
# MATHP

    def readMATS1(self,data):
        """
        MATS1(503,5,90) - record 12
        @todo add object
        """
        print "reading MATS1"
        while len(data)>=44: # 11*4
            eData = data[:44]
            data  = data[44:]
            out = unpack('iiifiiffiii',eData)
            (mid,tid,Type,h,yf,hr,limit1,limit2,a,b,c) = out
            dataIn = [mid,tid,Type,h,yf,hr,limit1,limit2]
            #mat = MATS1(None,dataIn)
            #self.addOp2Material(mat)
        ###

    def readMATT1(self,data):
        pass

    def readMATT2(self,data):
        pass

    def readMATT3(self,data):
        pass

    def readMATT4(self,data):
        pass

    def readMATT5(self,data):
        pass

# MATT8 - unused

    def readMATT9(self,data):
        pass

# MBOLT
# MBOLTUS
# MSTACK
# NLAUTO
# RADBND

    def readRADM(self,data):
        """
        RADM(8802,88,413) - record 25
        @todo add object
        """
        print "reading RADM"
        return
        while len(data)>=4: # 1*4
            eData = data[:4]
            data  = data[4:]
            number, = unpack('i',eData)
            
            strings = 'if'+'f'*number
            eDataLen = len(strings)*4
            
            eData = data[:eDataLen]
            data = data[eDataLen:]
            pack = list(unpack(strings,eData))
            packs = []
            
            while data:
                eData = data[:eDataLen]
                data = data[eDataLen:]
                pack = list(unpack(strings,eData))
                packs.append(pack)

            #mat = RADM(None,packs)
            #self.addOp2Material(mat)
        ###

# RADMT

    def readNLPARM(self,data):
        """
        NLPARM(3003,30,286) - record 27
        @todo add object
        """
        print "reading NLPARM"
        while len(data)>=76: # 19*4
            eData = data[:76]
            data  = data[76:]
            out = unpack('iifiiiiifffiiiffiff',eData)
            (sid,ninc,dt,kmethod,kstep,maxiter,conv,intout,epsu,epsp,epsw,
             maxdiv,maxqn,maxls,fstress,lstol,maxbis,maxr,rtolb) = out
            #mat = NLPARM(None,out)
            #self.addOp2Material(mat)
        ###

# NLPCI

    def readTSTEPNL(self,data):
        """
        TSTEPNL(3103,31,337) - record 29
        @todo add object
        """
        print "reading TSTEPNL"
        while len(data)>=88: # 19*4
            eData = data[:88]
            data  = data[88:]
            out = unpack('iifiiiiifffiiifiiiffff',eData)
            (sid,ndt,dt,no,kmethod,kstep,maxiter,conv,epsu,epsp,epsw,
             maxdiv,maxqn,maxls,fstress,lstol,maxbis,adjust,rb,maxr,utol,rtolb) = out
            #mat = TSTEPNL(None,out)
            #self.addOp2Material(mat)
        ###
