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
class MaxDisplacement(object):
    def __init__(self, data):
        self.translations = {}
        self.rotations = {}
        for line in data:
            sid = line[0]
            self.translations[sid] = line[1:4]
            self.rotations[sid] = line[4:]

    def write_f06(self, pageStamp='', pageNum=1):
        msg = ['0                                                  MAXIMUM  DISPLACEMENTS',
               '  SUBCASE/',
               '  DAREA ID        T1             T2             T3             R1             R2             R3']
        for sid, trans in sorted(self.translations.iteritems()):
            rot = self.rotations[sid]
            msg.append('0 %8i  %13.8E %13.8E %13.8E %13.8E %13.8E %13.8E' %
                       (tuple([sid] + trans + rot)))
        msg.append(pageStamp + str(pageNum))
        return '\n'.join(msg)
