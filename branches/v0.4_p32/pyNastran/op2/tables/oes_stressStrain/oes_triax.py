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
import sys
from pyNastran.op2.tables.oes_stressStrain.oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

class ctriaxStressObject(stressObject):
    """
    # formatCode=1 sortCode=0 stressCode=0
                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S
    ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  
       ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR
       5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
                4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CTRIAX6'
        
        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.radial  = {}
        self.azimuthal = {}
        self.axial = {}
        self.shear = {}
        self.omax = {}
        self.oms = {}
        self.ovm = {}

        if dt is not None:
            self.isTransient = True
            self.dt = self.nonlinearFactor
            self.addNewTransient()
            self.add = self.addTransient
            self.addNewEid = self.addNewEidTransient
        ###

    def addF06Data(self,data,transient):
        raise Exception('Not Implemented')
        if transient is None:
            for line in data:
                (eid,axial,MSa,torsion,MSt) = line
                if MSa==None: MSa = 0.
                if MSt==None: MSt = 0.
                self.axial[eid]      = axial
                self.MS_axial[eid]   = MSa
                self.torsion[eid]    = torsion
                self.MS_torsion[eid] = MSt
            ###
            return

        (dtName,dt) = transient
        self.dataCode['name'] = dtName
        if dt not in self.s1:
            self.updateDt(self.dataCode,dt)
            self.isTransient = True

        for line in data:
            (eid,axial,MSa,torsion,MSt) = line
            if MSa==None: MSa = 0.
            if MSt==None: MSt = 0.
            self.axial[dt][eid]      = axial
            self.MS_axial[dt][eid]   = MSa
            self.torsion[dt][eid]    = torsion
            self.MS_torsion[dt][eid] = MSt
        ###

    #def getLength_format1_sort0(self):
        #return (20,'iffff')

    def deleteTransient(self,dt):
        del self.radial[dt]
        del self.azimuthal[dt]
        del self.axial[dt]
        del self.shear[dt]
        del self.omax[dt]
        del self.oms[dt]
        del self.ovm[dt]

    def getTransients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        if self.dt not in self.axial:
            self.radial[self.dt]  = {}
            self.azimuthal[self.dt] = {}
            self.axial[self.dt] = {}
            self.shear[self.dt] = {}
            self.omax[self.dt] = {}
            self.oms[self.dt] = {}
            self.ovm[self.dt] = {}

    def addNewEid(self,eid,nid,rs,azs,As,ss,maxp,tmax,octs):
        #print "**?eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid]    = {nid: rs}
        self.azimuthal[eid] = {nid: azs}
        self.axial[eid]     = {nid: As}
        self.shear[eid]     = {nid: ss}
        self.omax[eid]      = {nid: maxp}
        self.oms[eid]       = {nid: tmax}
        self.ovm[eid]       = {nid: octs}

    def add(self,eid,nid,rs,azs,As,ss,maxp,tmax,octs):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid][nid]    = rs
        self.azimuthal[eid][nid] = azs
        self.axial[eid][nid]     = As
        self.shear[eid][nid]     = ss
        self.omax[eid][nid]      = maxp
        self.oms[eid][nid]       = tmax
        self.ovm[eid][nid]       = octs

    def addNewEidTransient(self,eid,nid,rs,azs,As,ss,maxp,tmax,octs):
        dt = self.dt
        #assert isinstance(eid,int)
        #assert eid >= 0
        #print "*  eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[dt][eid]    = {nid: rs}
        self.azimuthal[dt][eid] = {nid: azs}
        self.axial[dt][eid]     = {nid: As}
        self.shear[dt][eid]     = {nid: ss}
        self.omax[dt][eid]      = {nid: maxp}
        self.oms[dt][eid]       = {nid: tmax}
        self.ovm[dt][eid]       = {nid: octs}

    def addTransient(self,eid,nid,rs,azs,As,ss,maxp,tmax,octs):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        dt = self.dt
        self.radial[dt][eid][nid]    = rs
        self.azimuthal[dt][eid][nid] = azs
        self.axial[dt][eid][nid]     = As
        self.shear[dt][eid][nid]     = ss
        self.omax[dt][eid][nid]      = maxp
        self.oms[dt][eid][nid]       = tmax
        self.ovm[dt][eid][nid]       = octs

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)

        msg = header+['                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S\n',
               '   ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
               '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n',]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        #out = []
        for eid,radial in sorted(self.radial.items()):
            for nid in sorted(radial):
                rad   = self.radial[eid][nid]
                azimuth  = self.azimuthal[eid][nid]
                axial = self.axial[eid][nid] 
                shear = self.shear[eid][nid]
                omax  = self.omax[eid][nid]
                oms   = self.oms[eid][nid]
                ovm   = self.ovm[eid][nid]
                if nid==0:
                    Eid=eid
                else:
                    Eid=''
                ([rad,azimuth,axial,shear,omax,oms,ovm],isAllZeros) = self.writeF06Floats13E([rad,azimuth,axial,shear,omax,oms,ovm])
                msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' %(Eid,nid,radial,azimuth,axial,shear,omax,oms,ovm.rstrip()))
            ###
            msg.append('\n')
        ###

        msg.append(pageStamp+str(pageNum)+'\n')
        return(''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        words = ['                                      S T R E S S E S   I N   T R I A X 6   E L E M E N T S\n',
               '   ELEMENT  GRID ID       STRESSES  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
               '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n',]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        msg = []
        for dt,Radial in sorted(self.radial.items()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header+words
            for eid,radial in sorted(Radial.items()):
                for nid in sorted(radial):
                    rad   = self.radial[dt][eid][nid]
                    azimuth  = self.azimuthal[dt][eid][nid]
                    axial = self.axial[dt][eid][nid]
                    shear = self.shear[dt][eid][nid]
                    omax  = self.omax[dt][eid][nid]
                    oms   = self.oms[dt][eid][nid]
                    ovm   = self.ovm[dt][eid][nid]
                    if nid==0:
                        Eid=eid
                    else:
                        Eid=''
                    ([rad,azimuth,axial,shear,omax,oms,ovm],isAllZeros) = self.writeF06Floats13E([rad,azimuth,axial,shear,omax,oms,ovm])
                    msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' %(Eid,nid,rad,azimuth,axial,shear,omax,oms,ovm.rstrip()))
                ###
                msg.append('\n')
            ###

            msg.append(pageStamp+str(pageNum)+'\n')
        return(''.join(msg),pageNum-1)

    def __repr__(self):
        return self.writeF06(['',''],'PAGE ',1)[0]

class ctriaxStrainObject(strainObject):
    """
    # formatCode=1 sortCode=0 stressCode=0
                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S
    ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  
       ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR
       5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
                4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CTRIAX6'
        #raise NotImplementedError('CTRIAX6 strain...')
        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.radial  = {}
        self.azimuthal = {}
        self.axial = {}
        self.shear = {}
        self.emax = {}
        self.ems = {}
        self.evm = {}

        if dt is not None:
            self.isTransient = True
            self.dt = self.nonlinearFactor
            self.addNewTransient()
            self.add = self.addTransient
            self.addNewEid = self.addNewEidTransient
        ###

    def addF06Data(self,data,transient):
        raise Exception('Not Implemented')

    #def getLength_format1_sort0(self):
        #return (20,'iffff')

    def deleteTransient(self,dt):
        del self.radial[dt]
        del self.azimuthal[dt]
        del self.axial[dt]
        del self.shear[dt]
        del self.emax[dt]
        del self.ems[dt]
        del self.evm[dt]

    def getTransients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        if self.dt not in self.axial:
            self.radial[self.dt]  = {}
            self.azimuthal[self.dt] = {}
            self.axial[self.dt] = {}
            self.shear[self.dt] = {}
            self.emax[self.dt] = {}
            self.ems[self.dt] = {}
            self.evm[self.dt] = {}

    def addNewEid(self,eid,nid,rs,azs,As,ss,maxp,tmax,octs):
        #print "**?eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid]    = {nid: rs}
        self.azimuthal[eid] = {nid: azs}
        self.axial[eid]     = {nid: As}
        self.shear[eid]     = {nid: ss}
        self.emax[eid]      = {nid: maxp}
        self.ems[eid]       = {nid: emax}
        self.evm[eid]       = {nid: ects}

    def add(self,eid,nid,rs,azs,As,ss,maxp,emax,ects):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[eid][nid]    = rs
        self.azimuthal[eid][nid] = azs
        self.axial[eid][nid]     = As
        self.shear[eid][nid]     = ss
        self.emax[eid][nid]      = maxp
        self.ems[eid][nid]       = emax
        self.evm[eid][nid]       = ects

    def addNewEidTransient(self,eid,nid,rs,azs,As,ss,maxp,emax,ects):
        dt = self.dt
        #assert isinstance(eid,int)
        #assert eid >= 0
        #print "*  eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        self.radial[dt][eid]    = {nid: rs}
        self.azimuthal[dt][eid] = {nid: azs}
        self.axial[dt][eid]     = {nid: As}
        self.shear[dt][eid]     = {nid: ss}
        self.emax[dt][eid]      = {nid: maxp}
        self.ems[dt][eid]       = {nid: emax}
        self.evm[dt][eid]       = {nid: ects}

    def addTransient(self,eid,nid,rs,azs,As,ss,maxp,emax,ects):
        #print "***eid=%s loc=%s rs=%s azs=%s as=%s ss=%s maxp=%s tmx=%s octs=%s" %(eid,nid,rs,azs,As,ss,maxp,tmax,octs)
        dt = self.dt
        self.radial[dt][eid][nid]    = rs
        self.azimuthal[dt][eid][nid] = azs
        self.axial[dt][eid][nid]     = As
        self.shear[dt][eid][nid]     = ss
        self.emax[dt][eid][nid]      = maxp
        self.ems[dt][eid][nid]       = emax
        self.evm[dt][eid][nid]       = ects

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)

        msg = header+['                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S\n',
               '   ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
               '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n',]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        #out = []
        for eid,radial in sorted(self.radial.items()):
            for nid in sorted(radial):
                rad   = self.radial[eid][nid]
                azimuth  = self.azimuthal[eid][nid]
                axial = self.axial[eid][nid] 
                shear = self.shear[eid][nid]
                emax  = self.emax[eid][nid]
                ems   = self.ems[eid][nid]
                evm   = self.evm[eid][nid]
                if nid==0:
                    Eid=eid
                else:
                    Eid=''
                ([rad,azimuth,axial,shear,emax,ems,evm],isAllZeros) = self.writeF06Floats13E([rad,azimuth,axial,shear,emax,ems,evm])
                msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' %(Eid,nid,radial,azimuth,axial,shear,emax,ems,evm.rstrip()))
            ###
            msg.append('\n')
        ###

        msg.append(pageStamp+str(pageNum)+'\n')
        return(''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        words = ['                                      S T R A I N S   I N   T R I A X 6   E L E M E N T S\n',
               '   ELEMENT  GRID ID       STRAINS  IN  MATERIAL  COORD  SYSTEM                 MAX  MAG        MAX        VON MISES  \n',
               '      ID               RADIAL        AZIMUTHAL     AXIAL         SHEAR         PRINCIPAL       SHEAR\n',]
              #'      5351        0 -9.726205E+02 -1.678908E+03 -1.452340E+03 -1.325111E+02  -1.678908E+03  3.702285E+02  6.654553E+02
              #'               4389 -9.867789E+02 -1.624276E+03 -1.388424E+03 -9.212539E+01  -1.624276E+03  3.288099E+02  5.806334E+02

        msg = []
        for dt,Radial in sorted(self.radial.items()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header+words
            for eid,radial in sorted(Radial.items()):
                for nid in sorted(radial):
                    rad   = self.radial[dt][eid][nid]
                    azimuth  = self.azimuthal[dt][eid][nid]
                    axial = self.axial[dt][eid][nid]
                    shear = self.shear[dt][eid][nid]
                    emax  = self.emax[dt][eid][nid]
                    ems   = self.ems[dt][eid][nid]
                    evm   = self.evm[dt][eid][nid]
                    if nid==0:
                        Eid=eid
                    else:
                        Eid=''
                    ([rad,azimuth,axial,shear,emax,ems,evm],isAllZeros) = self.writeF06Floats13E([rad,azimuth,axial,shear,emax,ems,evm])
                    msg.append('  %8s %8s %s %s %s %s  %s %s %-s\n' %(Eid,nid,rad,azimuth,axial,shear,emax,ems,evm.rstrip()))
                ###
                msg.append('\n')
            ###

            msg.append(pageStamp+str(pageNum)+'\n')
        return(''.join(msg),pageNum-1)

    def __repr__(self):
        return self.writeF06(['',''],'PAGE ',1)[0]
