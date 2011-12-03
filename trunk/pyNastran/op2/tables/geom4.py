import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
#from pyNastran.bdf.cards.constraints import SPC,SPCADD

class Geometry4(object):

    def readTable_Geom4(self):
        self.iTableMap = {
                         #(4901, 49, 17):  self.readMPC,     # record 16
                         #(4891, 60, 83):  self.readMPCADD,  # record 17
                         #(6901,  69,295): self.readRBE2,    # record 24
                         #(5501,  55, 16): self.readSPC,     # record 44
                         #(5481,  58, 12): self.readSPC1,    # record 45
                         #(5491,  59, 13): self.readSPCADD,  # record 46
                         #(10100,101,472): self.readSUPORT1, # record 60
                         }
        self.readRecordTable('GEOM4')

    #def readSPC1
    #def readSPCADD

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
