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
from pyNastran.op2.tables.oqg_constraintForces.oqg_spcForces import spcForcesObject #,complexSpcForcesObject
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpcForces import mpcForcesObject #,complexMpcForcesObject


class OQG(object):
    def __init__(self):
        self.spcForces = {}
        self.mpcForces = {}

    def getSpcForces(self):
        (subcaseName,iSubcase,transient,dt,analysisCode,isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        data = self.readTable([int,str,float,float,float,float,float,float])

        dataCode = {'log':self.log,'analysisCode':analysisCode,'deviceCode':1,'tableCode':3,
                    'sortCode':0,'sortBits':[0,0,0],'numWide':8,'tableName':'OQG',}

        if iSubcase in self.spcForces:
            self.spcForces[iSubcase].addF06Data(data,transient)
        else:
            spc = spcForcesObject(dataCode,iSubcase)
            spc.addF06Data(data,transient)
            self.spcForces[iSubcase] = spc
        self.iSubcases.append(iSubcase)

    def getMpcForces(self):
        (subcaseName,iSubcase,transient,dt,analysisCode,isSort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        data = self.readTable([int,str,float,float,float,float,float,float])

        dataCode = {'log':self.log,'analysisCode':analysisCode,'deviceCode':1,'tableCode':39,
                    'sortCode':0,'sortBits':[0,0,0],'numWide':8,'tableName':'OQG',}

        if iSubcase in self.mpcForces:
            self.mpcForces[iSubcase].addF06Data(data,transient)
        else:
            mpc = mpcForcesObject(dataCode,iSubcase)
            mpc.addF06Data(data,transient)
            self.mpcForces[iSubcase] = mpc
        self.iSubcases.append(iSubcase)

