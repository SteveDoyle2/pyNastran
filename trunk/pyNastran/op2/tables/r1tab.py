from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import os
#import sys
#from struct import unpack

#from pyNastran.bdf.cards.nodes import GRID
#from pyNastran.bdf.cards.coordinateSystems import CORD1R,CORD2R,CORD2C,CORD3G #CORD1C,CORD1S,CORD2S


class R1TAB(object):

    def readTable_R1TAB(self):
        table_name = self.read_table_name(rewind=False)  # R1TAB
        self._table_init(table_name)
        #print "table_name = |%r|" %(table_name)

        self.read_markers([-1, 7], 'R1TAB')
        ints = self.read_int_block()
        #print "*ints = ",ints

        self.read_markers([-2, 1, 0], 'R1TAB')
        buffer_words = self.get_marker()
        #print "buffer_words = ",buffer_words
        word = self.read_string_block()  # DESTAB

        iTable = -3
        while buffer_words:  # read until buffer_words=0
            #print "iTable = ",iTable
            self.read_markers([iTable, 1, 0], 'R1TAB')
            nOld = self.n
            buffer_words = self.get_marker()
            #print "buffer_words = ",buffer_words
            if buffer_words == 0:  # maybe read new buffer...
                self.goto(nOld)
                break
            data = self.read_block()
            #print "len(data) = ",len(data)
            self.readDesvar(data)
            iTable -= 1

        #print self.print_section(200)

