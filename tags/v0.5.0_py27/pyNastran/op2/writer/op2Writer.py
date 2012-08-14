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
from struct import pack
from pyNastran.op2.writer.oesWriter import  Oes1Writer
from pyNastran.op2.writer.ougWriter import Ougv1Writer

class Op2Writer(Ougv1Writer,Oes1Writer):
    self.hollerith = pack('i',584) # assumes a 128 character string
    def pack(self,format,*vals):
        return [str(val) for val in vals]

    def writeStart(self):
        self.writeMarkers([3,7])
        self.writeStringBlock('NASTRAN FORT TAPE ID CODE - ',28)
        self.writeMarkers([2,-1])

    def writeOp2(self,op2Name):
        op2 = open(op2Name,'wb')
        self.writeStart()
        #self.writeGEOM1()
        #self.writeGEOM2()
        #self.writeGEOM3()
        #self.writeGEOM4()
        
        #self.writeEPT()
        #self.writeMPTS()

        #self.writeOQG1() # spc forces

        #self.writeOGP()
        self.writeOUGV1() # displacements/temps/heatFlux
        #self.writeOGF1()
        self.writeOES1()  # stress
        
    def printHeader(word,nChars):
        self.deviceCode = 1 # print the OP2...

        self.writeStringBlock(word,nChars)

        msg += writeMarkers([-1,7])
        out = [101, 0, 4136, 0, 0, 0, 1]  ## @todo what this is - DMAP -> "no def or month,year,one,one"...huh???
        msg += pack('iiiiiii',*out)

        msg += writeMarkers([-2,1,0])

        # approachCode=1, tableCode=1
        self.iTable = -3

    def writeStringBlock(self,word,nChars):
        """
        'OUG'      - 1 word  = 4 characters
        'OUGV'     - 1 word  = 4 characters
        'OUGV1'    - 2 words = 8 characters
        '12345678' - 2 words = 8 characters
        nWords = round(ceil(len(word)/4.))
        nChars = nWords*4 != len(word)
        just set nChars and dont overthink it too much
        """
        m = self.writeMarkers([1])
        #word = self.readStringBlock()
        ## creating a %Xs - where x is the number of 
        ## words*4,left justifying it, and setting the value
        value = '%%-%ss' %(nChars) %(word)

        out = pack('c'*nChars,value)
        out = ''
        return m+out+m
        
    def packTitle(iSubcase):
        Title = self.Title
        subtitle,label = self.iSubcaseNameMap[iSubcase]
        titleSubtitleLabel = '%128s%128s%128s' %(Title,subtitle,label)
        msg = pack('c'*384,list(titleSubtitleLabel)) + self.hollerith
        return msg

    def writeMarkers(markers):
        """
        takes -5,1,0  -> [4,5,4,  4,1,4,  4,0,4]
        and puts it into binary
        """
        out = []
        for marker in markers:
            out += [4,marker,4]
        n = len(out)
        return pack('i'*n,*out)
        
    def combineApproachDeviceCodes(self,approachCode):
        aCode = approachCode*10 + self.deviceCode
        return aCode

    def combineTableSortCodes(self,tableCode,sortCode):
        tCode = sortCode*1000+tableCode
        return tCode

    def aCode_tCode(approachCode,tableCode,sortCode):
        aCode = self.combineApproachDeviceCodes(approachCode)
        tCode = self.combineTableDeviceCodes(tableCode,sortCode)
        return (aCode,tCode)
