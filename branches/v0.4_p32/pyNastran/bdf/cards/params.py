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
# my code
from pyNastran.bdf.cards.baseCard import BaseCard

class PARAM(BaseCard):
    type = 'PARAM'
    def __init__(self,card):
        self.key   = card.field(1)
        self.value = card.field(2)

    #def isSameCard(self,param):
        #fields1 = [self.key, self.value ]
        #fields2 = [param.key,param.value]
        #for (field1,field2) in zip(fields1,fields2):
        #    if not self.isSame(field1,field2):
        #        return False
        #    ###
        ####
        #return True

    def rawFields(self):
        fields = ['PARAM',self.key,self.value]
        return fields
