import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
#from pyNastran.bdf.cards.constraints import SPC,SPCADD

class Geometry4(object):

    def readTable_Geom4(self):
        self.iTableMap = {
                         #(5481, 58, 12): self.
                         #(5491, 59, 13): self.
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
