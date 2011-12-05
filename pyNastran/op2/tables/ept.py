import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.properties import PROD,PBAR,PBARL,PBEAM,PSHELL,PSOLID,PCOMP,PTUBE

class EPT(object):

    def readTable_EPT(self):
        self.bigProperties = {}
        self.iTableMap = {
                         (3201,32,55):    self.readNSM,     # record 2
                         (52,20,181):     self.readPBAR,    # record 11 - buggy
                         (9102,91,52):    self.readPBARL,   # record 12 - no PBARL object
                         #(5402,54,262):  self.readPBEAM,   # record 14 - not done
                         (2706,27,287):   self.readPCOMP,   # record 22 - buggy
                         (902,9,29):      self.readPROD,    # record 49
                         #(1002,10,42):    self.readPSHEAR, # record 50 - no PSHEAR object
                         (2402,24,281):   self.readPSOLID,  # record 51
                         (2302,23,283):   self.readPSHELL,  # record 52
                         (1602,16,30):    self.readPTUBE,   # record 56
                         
                         (5402, 54, 262): self.readFake,
                         (11001,110,411): self.readFake,
                         (2802, 28, 236): self.readFake,
                         (2606, 26, 289): self.readFake,
                         (152,  19, 147): self.readFake,
                         (2102, 21, 121): self.readFake,
                         

                         }
        self.readRecordTable('EPT')

    def addOp2Property(self,prop):
        self.addProperty(prop,allowOverwrites=True)

# HGSUPPR

    def readNSM(self,data):
        """
        NSM(3201,32,55) - the marker for Record 2
        """
        print "reading NSM"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            out = unpack('iccccif',eData)
            (sid,A,B,C,D,ID,value) = out
            propSet = A+B+C+D
            print "sid=%s propSet=%s ID=%s value=%s" %(sid,propSet,ID,value)
            prop = NSM(None,[sid,propSet,ID,value])
            self.addOp2Property(prop)
        ###


# NSM1
# NSML1
# NSMADD
# NSML
# NSML1
# PAABSF
# PACABS
# PACBAR

    def readPBAR(self,data):
        """
        PBAR(52,20,181) - the marker for Record 11
        @warning this makes a funny property...
        """
        #print "reading PBAR"
        while len(data)>=76: # 19*4
            eData = data[:76]
            data  = data[76:]
            out = unpack('iifffffffffffffffff',eData)
            #print "len(out) = ",len(out)
            #print out
            (pid,mid,a,I1,I2,J,nsm,fe,c1,c2,d1,d2,e1,e2,f1,f2,k1,k2,I12) = out
            prop = PBAR(None,out)
            self.addOp2Property(prop)
        ###


    def readPBARL(self,data):
        """
        PBARL(9102,91,52) - the marker for Record 12
        @todo create object
        """
        print "reading PBARL"
        #while len(data)>=28: # 7*4
        if 1:
            eData = data[:28]
            data  = data[28:]
            out = unpack('iiccccccccccccccccf',eData)
            (pid,mid,a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p,value) = out
            group1 = a+b+c+d
            group2 = e+f+g+h
            type1  = i+j+k+l
            type2  = m+n+o+p
            dataIn = [pid,mid,group1+group2,type1+type2,value]
            print "pid=%s mid=%s group1=%s group2=%s type1=%s type2=%s value=%s" %(pid,mid,group1,group2,type1,type2,value)
            
            while len(data)>4:
                value, = unpack('f',data[:4])
                #print "valueNew = %s" %(value)
                data = data[4:]
                dataIn.append(value)
            
            #print "len(out) = ",len(out)
            #print "PBARL = ",dataIn
            prop = PBARL(None,dataIn)
            self.addProperty(prop)
            #print self.printSection(20)
            #sys.exit()
        ###


# PBCOMP

    def readPBEAM(self,data):
        """
        PBEAM(5402,54,262) - the marker for Record 14
        @todo add object
        """
        print "reading PBEAM"
        while len(data)>=1072: # 44+12*84+20
            eData = data[:20]
            data  = data[20:]
            dataIn = list(unpack('iiiif',eData))
            #print "len(out) = ",len(out)
            #print out
            (pid,mid,nsegs,ccf,x) = dataIn

            for i in range(12):
                eData = data[84:]
                data  = data[:84]
                pack = unpack('ffffffffffffffff',eData)
                (so,xxb,a,i1,i2,i12,j,nsm,c1,c2,d1,d2,e1,e2,f1,f2) = pack
                dataIn.append(pack)
            
            eData = data[:44]
            data  = data[44:]

            dataIn = list(unpack('fffffffffff',eData))
            (k1,k2,s1,s2,nsia,nsib,cwa,cwb,m1a,m2a,m1b,m2b,n1a,n2a,n1b,n2b) = pack
            
                
            prop = PBEAM(None,out)
            self.addOp2Property(prop)
            sys.exit('ept-PBEAM')
        ###


# PBEAML
# PBEND
# PBMSECT
# PBRSECT
# PBUSH
# PBUSH1D
# PBUSHT

    def readPCOMP(self,data):
        """
        PCOMP(2706,27,287) - the marker for Record 22
        """
        #print "reading PCOMP"
        while len(data)>=32: # 8*4 - dynamic
            eData = data[:32]
            data  = data[32:]
            out = unpack('iifffiff',eData)
            (pid,nLayers,z0,nsm,sb,ft,Tref,ge,) = out
            
            eData = data[:16*(nLayers)]
            data  = data[16*(nLayers):]
            
            Mid=[]; T=[]; Theta=[]; Sout=[]
            for n in range(nLayers):
                #print "len(eData) = ",len(eData)
                (mid,t,theta,sout) = unpack('iffi',eData[0:16])
                Mid.append(mid)
                T.append(t)
                Theta.append(theta)
                Sout.append(sout)
                eData = eData[16:]
            ###
            
            dataIn = [pid,z0,nsm,sb,ft,Tref,ge,Mid,T,Theta,Sout]
            #print "PCOMP = %s" %(dataIn)
            prop = PCOMP(None,dataIn)
            self.addOp2Property(prop)
        ###

# PCOMPA
# PCONEAX
# PCONV
# PCONVM
# PDAMP
# PDAMPT
# PDAMP5
# PDUM1
# PDUM2
# PDUM3
# PDUM4
# PDUM5
# PDUM6
# PDUM7
# PDUM8
# PDUM9
# PELAS
# PFAST
# PELAST
# PGAP
# PHBDY
# PINTC
# PINTS
# PLPLANE
# PLSOLID
# PMASS

    def readPROD(self,data):
        """
        PROD(902,9,29) - the marker for Record 49
        """
        print "reading PSHEAR"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            out = unpack('iiffff',eData)
            (pid,mid,a,j,c,nsm) = out
            prop = PROD(None,out)
            self.addOp2Property(prop)
        ###

    def readPSHEAR(self,data):
        """
        PSHEAR(1002,10,42) - the marker for Record 50
        """
        print "reading PSHEAR"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            out = unpack('iiffff',eData)
            (pid,mid,t,nsm,f1,f2) = out
            prop = PSHEAR(None,out)
            self.addOp2Property(prop)
        ###

    def readPSHELL(self,data):
        """
        PSHELL(2302,23,283) - the marker for Record 51
        """
        #print "reading PSHELL"
        while len(data)>=44: # 11*4
            eData = data[:44]
            data  = data[44:]
            out = unpack('iifififfffi',eData)
            (pid,mid1,t,mid2,bk,mid3,ts,nsm,z1,z2,mid4) = out
            prop = PSHELL(None,out)

            if max(pid,mid1,mid2,mid3,mid4)>1e8:
                #print "PSHELL = ",out
                self.bigProperties[pid] = prop
            else:
                self.addOp2Property(prop)
            ###
        ###


    def readPSOLID(self,data):
        """
        PSOLID(2402,24,281) - the marker for Record 52
        """
        print "reading PSOLID"
        while len(data)>=28: # 7*4
            eData = data[:28]
            data  = data[28:]
            (pid,mid,cid,inp,stress,isop,fctnA,fctnB,fctnC,fctnD) = unpack('iiiiiicccc',eData)
            fctn = fctnA+fctnB+fctnC+fctnD
            dataIn = [pid,mid,cid,inp,stress,isop,fctn]
            prop = PSOLID(None,dataIn)
            self.addOp2Property(prop)
        ###

# PSOLIDL
# PTRIA6
# PTSHELL

    def readPTUBE(self,data):
        """
        PTUBE(1602,16,30) - the marker for Record 56
        @todo OD2 only exists for heat transfer...how do i know if there's heat transfer...
        @warning assuming OD2 is not written (only done for thermal)
        """
        print "reading PTUBE"
        while len(data)>=20: # 5*4
            eData = data[:20]
            data  = data[20:] # or 24
            (pid,mid,OD,t,nsm) = unpack('iifff',eData)
            dataIn = [pid,mid,OD,t,nsm]
            prop = PTUBE(None,dataIn)
            self.addOp2Property(prop)
        ###

# PSET
# PVAL
# PVISC
# PWELD
# PWSEAM
# PVIEW
# PVIEW3D
