from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import os
#import sys
#from struct import unpack

#from pyNastran.bdf.cards.nodes import GRID
#from pyNastran.bdf.cards.coordinateSystems import CORD1R,CORD2R,CORD2C,CORD3G #CORD1C,CORD1S,CORD2S


class R1TAB(object):

    def readTable_R1TAB(self):
        tableName = self.readTableName(rewind=False)  # R1TAB
        self.tableInit(tableName)
        #print "tableName = |%r|" %(tableName)

        self.readMarkers([-1, 7], 'R1TAB')
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-2, 1, 0], 'R1TAB')
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords
        word = self.readStringBlock()  # DESTAB

        iTable = -3
        while bufferWords:  # read until bufferWords=0
            #print "iTable = ",iTable
            self.readMarkers([iTable, 1, 0], 'R1TAB')
            nOld = self.n
            bufferWords = self.getMarker()
            #print "bufferWords = ",bufferWords
            if bufferWords == 0:  # maybe read new buffer...
                self.goto(nOld)
                break
            data = self.readBlock()
            #print "len(data) = ",len(data)
            self.readDesvar(data)
            iTable -= 1

        #print self.printSection(200)

