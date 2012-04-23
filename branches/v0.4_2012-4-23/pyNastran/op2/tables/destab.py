## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
#from pyNastran.bdf.cards.nodes import GRID
#from pyNastran.bdf.cards.coordinateSystems import CORD1R,CORD2R,CORD2C,CORD3G #CORD1C,CORD1S,CORD2S


class DESTAB(object):

    def readDesvar(self,data):
        out = unpack('iiccccccccfff',data[0:28])
        (idvid,dvid,alabel1,bLabel1,cLabel1,dLabel1,alabel2,bLabel2,cLabel2,dLabel2,vmin,vmax,delx) = out
        label1 = alabel1+bLabel1+cLabel1+dLabel1
        label2 = alabel2+bLabel2+cLabel2+dLabel2
        #print "ivid=%s dvid=%s label1=%s label2=%s vmax=%g vmin=%g delx=%g" %(idvid,dvid,label1,label2,vmin,vmax,delx)
        #self.readData(8)

    def readTable_DesTab(self):
        self.tableName = 'DESTAB'
        tableName = self.readTableName(rewind=False) # DESTAB
        self.tableInit(tableName)
        #print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7],'DESTAB')
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-2,1,0],'DESTAB')
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords
        word = self.readStringBlock() # DESTAB
        #print "word = |%s|" %(word)

        iTable = -3
        #imax   = -244
        
        while bufferWords: # read until bufferWords=0
            self.readMarkers([iTable,1,0],'DESTAB')
            nOld = self.n
            bufferWords = self.getMarker()
            #print "bufferWords = ",bufferWords
            if bufferWords==0: # maybe read new buffer...
                self.goto(nOld)
                break
            data = self.readBlock()
            #print "len(data) = ",len(data)
            self.readDesvar(data)
            iTable -= 1

        self.printSection(80)

        #self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "1-bufferWords = ",bufferWords,bufferWords*4

        #print self.printSection(300)
        #sys.exit('asdf')

