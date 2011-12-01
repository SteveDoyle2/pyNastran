import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.properties import PROD,PBAR,PBEAM,PSHELL,PSOLID,PCOMP,PTUBE

class EPT(object):

    def readTable_EPT(self):
        self.bigProperties = {}
        self.iTableMap = {
                         (3201,32,55):    self.readNSM,     # record 2
                         (52,20,181):     self.readPBar,   # record 11 - buggy
                         (2706,27,287):   self.readPComp,   # record 22
                         (902,9,29):      self.readPRod,    # record 49
                         #(1002,10,42):    self.readPShear, # record 50
                         (2402, 24, 281): self.readPSolid,  # record 51
                         (2302,23,283):   self.readPShell,  # record 52
                         (1602,16,30):    self.readPTube,   # record 56
                         }
        self.readRecordTable('EPT')

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
            self.addProperty(prop)
        ###


# NSM1
# NSML1
# NSMADD
# NSML
# NSML1
# PAABSF
# PACABS
# PACBAR

    def readPBar(self,data):
        """
        PBAR(52,20,181) - the marker for Record 11
        @warning this makes a funny property...
        """
        print "reading PBAR"
        while len(data)>=76: # 19*4
            eData = data[:76]
            data  = data[76:]
            out = unpack('iifffffffffffffffff',eData)
            #print "len(out) = ",len(out)
            print out
            (pid,mid,a,I1,I2,J,nsm,fe,c1,c2,d1,d2,e1,e2,f1,f2,k1,k2,I12) = out
            prop = PBAR(None,out)
            self.addProperty(prop)
            #sys.exit()
        ###

# PBARL
# PBCOMP
# PBEAM
# PBEAML
# PBEND
# PBMSECT
# PBRSECT
# PBUSH
# PBUSH1D
# PBUSHT

    def readPComp(self,data):
        """
        PCOMP(2706,27,287) - the marker for Record 22
        """
        print "reading PCOMP"
        while len(data)>=32: # 8*4 - dynamic
            eData = data[:32]
            data  = data[32:]
            out = unpack('iifffiff',eData)
            (pid,nLayers,z0,nsm,sb,ft,Tref,ge,) = out
            
            eData = data[:16*(nLayers+1)]
            data  = data[16*(nLayers+1):]
            
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
            print "PCOMP = %s" %(dataIn)
            prop = PCOMP(None,dataIn)
            self.addProperty(prop)
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

    def readPRod(self,data):
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
            self.addProperty(prop)
        ###

    def readPShear(self,data):
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
            self.addProperty(prop)
        ###

    def readPShell(self,data):
        """
        PSHELL(2302,23,283) - the marker for Record 51
        """
        print "reading PSHELL"
        while len(data)>=44: # 11*4
            eData = data[:44]
            data  = data[44:]
            out = unpack('iifififfffi',eData)
            print "PSHELL = ",out
            (pid,mid1,t,mid2,bk,mid3,ts,nsm,z1,z2,mid4) = out
            prop = PSHELL(None,out)

            if max(pid,mid1,mid2,mid3,mid4)>1e8:
                self.bigProperties[pid] = prop
            else:
                self.addProperty(prop)
            ###
        ###


    def readPSolid(self,data):
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
            self.addProperty(prop)
        ###

# PSOLIDL
# PTRIA6
# PTSHELL

    def readPTube(self,data):
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
            self.addProperty(prop)
        ###

# PSET
# PVAL
# PVISC
# PWELD
# PWSEAM
# PVIEW
# PVIEW3D
