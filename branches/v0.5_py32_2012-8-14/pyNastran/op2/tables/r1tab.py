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


#import os
#import sys
#from struct import unpack

#from pyNastran.bdf.cards.nodes import GRID
#from pyNastran.bdf.cards.coordinateSystems import CORD1R,CORD2R,CORD2C,CORD3G #CORD1C,CORD1S,CORD2S


class R1TAB(object):

    def readTable_R1TAB(self):
        tableName = self.readTableName(rewind=False) # R1TAB
        self.tableInit(tableName)
        #print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7],'R1TAB')
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-2,1,0],'R1TAB')
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords
        word = self.readStringBlock() # DESTAB

        iTable = -3
        while bufferWords: # read until bufferWords=0
            #print "iTable = ",iTable
            self.readMarkers([iTable,1,0],'R1TAB')
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

        #print self.printSection(200)
        
        #sys.exit('R1TAB in r1tab.py')

