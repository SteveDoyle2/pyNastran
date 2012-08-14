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

class Oes1Writer(object):
    def writeOES1(self):
        """
        writes isotropic/composite stress/strain
        @todo assumes sCode=0 (stress) or 10 (strain)
        """
        msg += self.printHeader('OES1X1',8)
        # OES1X1
        stress = [
              self.rodStress,
              self.barStress,  
              self.beamStress, 
              self.plateStress,
              self.solidStress,
        ]
        cstress = [self.compositePlateStress,] # OES1C

        # approachCode=1, tableCode=1
        for data in stress:
            for iSubcase in data:
                msg += self.writeMarkers([self.iTable,1,0])
                msg += self.writeOES(iSubcase,data)
                self.iTable-=1
        msg += self.writeMarkers([self.iTable,1,0])
        
        # ----------------
        msg += self.printHeader('OSTR1X',8)
        # OSTR1X
        strain = [
              self.rodStrain,
              self.barStrain,
              self.beamStrain,
              self.plateStrain,
              self.solidStrain,
              self.compositePlateStrain, 
        ]
        for data in strain:
            for iSubcase in data:
                msg += self.writeMarkers([self.iTable,1,0])
                msg += self.writeOES(iSubcase,data)
                self.iTable-=1
        msg += self.writeMarkers([self.iTable,1,0])

        # ----------------
        msg = self.printHeader('OSTRIC',8)
        cstrain = [self.compositePlateStrain,] # OSTR1C
        for data in cstrain:
            for iSubcase in data:
                msg += self.writeMarkers([self.iTable,1,0])
                msg += self.writeOES(iSubcase,data)
                self.iTable-=1
        msg += self.writeMarkers([self.iTable,1,0])


