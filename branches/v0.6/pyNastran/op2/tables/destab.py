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
from struct import unpack

#from pyNastran.bdf.cards.nodes import GRID
#from pyNastran.bdf.cards.coordinateSystems import CORD1R,CORD2R,CORD2C,CORD3G #CORD1C,CORD1S,CORD2S


class DESTAB(object):

    def readDesvar(self, data):
        out = unpack(b'ii4s4sfff', data[0:28])
        (idvid, dvid, label1, label2, vmin, vmax, delx) = out
        #print "ivid=%s dvid=%s label1=%s label2=%s vmax=%g vmin=%g delx=%g" %(idvid,dvid,label1,label2,vmin,vmax,delx)
        #self.readData(8)

    def readTable_DesTab(self):
        self.table_name = 'DESTAB'
        table_name = self.read_table_name(rewind=False)  # DESTAB
        self.table_init(table_name)
        #print "table_name = |%r|" % (table_name)

        self.read_markers([-1, 7], 'DESTAB')
        ints = self.read_int_block()
        #print "*ints = ",ints

        self.read_markers([-2, 1, 0], 'DESTAB')
        buffer_words = self.get_marker()
        #print "buffer_words = ",buffer_words
        word = self.read_string_block()  # DESTAB
        #print "word = |%s|" %(word)

        iTable = -3
        #imax   = -244

        while buffer_words:  # read until buffer_words=0
            self.read_markers([iTable, 1, 0], 'DESTAB')
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

        self.print_section(80)

        #self.op2Debug.write('buffer_words=%s\n' %(str(buffer_words)))
        #print "1-buffer_words = ",buffer_words,buffer_words*4

        #print self.print_section(300)
