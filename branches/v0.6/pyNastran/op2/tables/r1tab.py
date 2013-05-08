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
        self.table_init(table_name)
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

