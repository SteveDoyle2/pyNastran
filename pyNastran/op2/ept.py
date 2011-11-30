import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.properties import PROD,PSHELL,PSOLID

class EPT(object):

    def readTable_EPT(self):
        self.iTableMap = {
                         (902,9,29):      self.readPRod,
                         #(1002,10,42):    self.readPShear,
                         (2402, 24, 281): self.readPSolid,
                         (2302,23,283):   self.readPShell,
                         
                         }
        self.readRecordTable('EPT')

# HGSUPPR
# NSM
# NSM1
# NSML1
# NSMADD
# NSML
# NSML1
# PAABSF
# PACABS
# PACBAR
# PBAR
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
# PCOMP
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
            (pid,mid1,t,mid2,bk,mid3,ts,nsm,z1,z2,mid4) = out
            prop = PSHELL(None,out)
            self.addProperty(prop)
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
# PTUBE
# PSET
# PVAL
# PVISC
# PWELD
# PWSEAM
# PVIEW
# PVIEW3D
