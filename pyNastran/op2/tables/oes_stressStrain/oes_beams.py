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
from struct import unpack
from oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

class beamStressObject(stressObject):
    """
    [1,0,0]
                 S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )
                      STAT DIST/
     ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C
            1       1   0.000   -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04          
                    2   1.000   -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04          

    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CBEAM'
        
        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.xxb = {}
        self.grids = {}
        self.smax = {}
        self.smin = {}
        self.MS_tension = {}
        self.MS_compression = {}
        
        if self.code in [[1,0,0]]: # ,[1,0,1]
            #self.MS_axial   = {}
            #self.MS_torsion = {}
            self.sxc = {}
            self.sxd = {}
            self.sxe = {}
            self.sxf = {}
            self.getLength1     = self.getLength1_format1_sort0
            self.getLength2     = self.getLength2_format1_sort0
            self.getLengthTotal = self.getLengthTotal_format1_sort0

            #self.isImaginary = False
            if dt is not None:
                self.addNewTransient = self.addNewTransient_format1_sort0
                self.addNewEid       = self.addNewEidTransient_format1_sort0
                self.add             = self.addTransient_format1_sort0
            else:
                self.addNewEid = self.addNewEid_format1_sort0
                self.add       = self.add_format1_sort0
            ###
        #elif self.code==[2,1,0]:
        #    self.getLength       = self.getLength_format1_sort0
        #    self.addNewTransient = self.addNewTransient_format2_sort1
        #    self.addNewEid       = self.addNewEidTransient_format2_sort1
        #    #self.isImaginary = True
        else:
            raise InvalidCodeError('beamStress - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt is not None:
            self.isTransient = True
            self.dt = self.nonlinearFactor
            self.addNewTransient()
        ###

    def getLengthTotal_format1_sort0(self):
        return 444  # 44+10*40   (11 nodes)

    def getLength1_format1_sort0(self):
        return (44,'iifffffffff')

    def getLength2_format1_sort0(self):
        return (40,'ifffffffff')

    def deleteTransient(self,dt):
        del self.sxc[dt]
        del self.sxd[dt]
        del self.sxe[dt]
        del self.sxf[dt]
        del self.smax[dt]
        del self.smin[dt]
        del self.MS_tension[dt]
        del self.MS_compression[dt]

    def getTransients(self):
        k = self.smax.keys()
        k.sort()
        return k

    def addNewTransient_format1_sort0(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        #print "addNewTransient_beam+1+0"
        if self.dt not in self.smax:
            self.sxc[self.dt] = {}
            self.sxd[self.dt] = {}
            self.sxe[self.dt] = {}
            self.sxf[self.dt] = {}
            self.smax[self.dt] = {}
            self.smin[self.dt] = {}
            self.MS_tension[self.dt]     = {}
            self.MS_compression[self.dt] = {}

    def addNewTransient_format2_sort1(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        raise Exception('not supported')
        #print self.dataCode
        self.axial[self.dt]     = {}
        self.torsion[self.dt]   = {}

    def addNewEid_format1_sort0(self,out):
        #print "Beam Stress addNewEid..."
        (eid,grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
        eid = (eid-self.deviceCode)/10
        #print "eid=%s grid=%s" %(eid,grid)
        assert eid >= 0
        #assert isinstance(eid,int)
        #assert isinstance(grid,int)
        self.grids[eid] = [grid]
        self.xxb[eid]  = [sd]
        self.sxc[eid] = [sxc]
        self.sxd[eid] = [sxd]
        self.sxe[eid] = [sxe]
        self.sxf[eid] = [sxf]
        self.smax[eid] = [smax]
        self.smin[eid] = [smin]
        self.MS_tension[eid] = [mst]
        self.MS_compression[eid] = [msc]
        return eid

    def addNewEid_format2_sort1(self,out):
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        assert eid >= 0
        self.axial[eid]      = [axialReal,axialImag]
        self.torsion[eid]    = [torsionReal,torsionImag]
        return eid

    def addNewEidTransient_format1_sort0(self,out):
        #print "Beam Transient Stress addNewEid..."
        (eid,grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out

        eid = (eid-self.deviceCode)/10
        dt = self.dt
        assert eid  >= 0
        self.grids[eid] = [grid]
        self.xxb[eid] = [sd]
        self.sxc[dt][eid] = [sxc]
        self.sxd[dt][eid] = [sxd]
        self.sxe[dt][eid] = [sxe]
        self.sxf[dt][eid] = [sxf]
        self.smax[dt][eid] = [smax]
        self.smin[dt][eid] = [smin]
        self.MS_tension[dt][eid] = [mst]
        self.MS_compression[dt][eid] = [msc]
        return eid

    def addNewEidTransient_format2_sort1(self,out):
        raise Exception('not supported')
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        dt = self.dt
        assert eid >= 0
        self.axial[dt][eid]      = [axialReal,axialImag]
        self.torsion[dt][eid]    = [torsionReal,torsionImag]
        return eid

    def add_format1_sort0(self,eid,out):
        #print "Beam Stress add..."
        (grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
        if grid:
            self.grids[eid].append(grid)
            self.xxb[eid].append(sd)
            self.sxc[eid].append(sxc)
            self.sxd[eid].append(sxd)
            self.sxe[eid].append(sxe)
            self.sxf[eid].append(sxf)
            self.smax[eid].append(smax)
            self.smin[eid].append(smin)
            self.MS_tension[eid].append(mst)
            self.MS_compression[eid].append(msc)
        ###

    def addTransient_format1_sort0(self,eid,out):
        #print "Beam Transient Stress add..."
        (grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
        dt = self.dt
        if grid:
            self.grids[eid].append(grid)
            self.xxb[eid].append(sd)
            #self.sd[dt][eid].append(sd)
            self.sxc[dt][eid].append(sxc)
            self.sxd[dt][eid].append(sxd)
            self.sxe[dt][eid].append(sxe)
            self.sxf[dt][eid].append(sxf)
            self.smax[dt][eid].append(smax)
            self.smin[dt][eid].append(smin)
            self.MS_tension[dt][eid].append(mst)
            self.MS_compression[dt][eid].append(msc)
        ###

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)

        msg = header + ['                                  S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                        '                    STAT DIST/\n',
                        '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']

        for eid in sorted(self.smax):
            msg.append('0  %8i\n' %(eid))
            #print self.xxb[eid]
            for i,nid in enumerate(self.grids[eid]):
                #print i,nid
                xxb  = self.xxb[eid][i]
                sxc = self.sxc[eid][i]
                sxd = self.sxd[eid][i]
                sxe = self.sxe[eid][i]
                sxf = self.sxf[eid][i]
                sMax = self.smax[eid][i]
                sMin = self.smin[eid][i]
                SMt  = self.MS_tension[eid][i]
                SMc  = self.MS_compression[eid][i]
                (vals2,isAllZeros) = self.writeF06Floats13E([sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc])
                (sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc) = vals2
                msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' %(nid,xxb,sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc.strip()))
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        words = ['                                  S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                 '                    STAT DIST/\n',
                 '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        msg = []
        for dt,SMaxs in sorted(self.smax.items()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header + words
            for eid,Smax in sorted(SMaxs.items()):
                msg.append('0  %8i\n' %(eid))
                for i,nid in enumerate(self.grids[eid]):
                    xxb  = self.xxb[eid][i]
                    #sd  = self.sd[eid][i]
                    sxc = self.sxc[dt][eid][i]
                    sxd = self.sxd[dt][eid][i]
                    sxe = self.sxe[dt][eid][i]
                    sxf = self.sxf[dt][eid][i]
                    sMax = self.smax[dt][eid][i]
                    sMin = self.smin[dt][eid][i]
                    SMt  = self.MS_tension[dt][eid][i]
                    SMc  = self.MS_compression[dt][eid][i]
                    (vals2,isAllZeros) = self.writeF06Floats13E([sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc])
                    (sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc) = vals2
                    msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' %(nid,xxb,sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc.strip()))
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        return (''.join(msg),pageNum-1)

    def __reprTransient_format1_sort0__(self):
        msg = '---BEAM STRESSES---\n'
        msg += '%-6s %6s %6s %7s' %('EID','eType','NID','xxb')
        headers = ['sMax','sMin','MS_tension','MS_compression']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,smax in sorted(self.smax.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid in sorted(smax):
                for i,nid in enumerate(self.grids[eid]):
                    xxb  = self.xxb[eid][i]
                    sMax = self.smax[dt][eid][i]
                    sMin = self.smin[dt][eid][i]
                    SMt  = self.MS_tension[dt][eid][i]
                    SMc  = self.MS_compression[dt][eid][i]
                    xxb = round(xxb,2)

                    msg += '%-6i %6s %6i %7.2f ' %(eid,self.eType,nid,xxb)
                    vals = [sMax,sMin,SMt,SMc]
                    for val in vals:
                        if abs(val)<1e-6:
                            msg += '%10s ' %('0')
                        else:
                            msg += '%10g ' %(val)
                        ###
                    msg += '\n'
            ###
        #print msg
        #sys.exit('beamT')
        return msg

    def __reprTransient_format2_sort1__(self):
        raise Exception('not implemented...')
        msg = '---COMPLEX BEAM STRESSES---\n'
        msg += '%-10s %10s ' %('EID','eType')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,axial in sorted(self.axial.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid in sorted(axial):
                axial   = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                msg += '%-6i %6s ' %(eid,self.eType)
                vals = axial + torsion # concatination
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10i ' %(val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
            ###
        return msg

    def __repr__(self):
        if   self.isTransient and self.code in [[1,0,0],[1,0,1]]:
            return self.__reprTransient_format1_sort0__()
        elif self.code==[2,1,0]:
            return self.__reprTransient_format2_sort1__()
        #else:
        #    raise Exception('code=%s' %(self.code))
        msg = '---BEAM STRESSES---\n'
        msg += '%-6s %6s %6s %6s' %('EID','eType','NID','xxb')
        headers = ['sMax','sMin','MS_tension','MS_compression']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'
        #print "self.code = ",self.code
        for eid in sorted(self.smax):
            #print self.xxb[eid]
            for i,nid in enumerate(self.grids[eid]):
                #print i,nid
                xxb  = self.xxb[eid][i]
                sMax = self.smax[eid][i]
                sMin = self.smin[eid][i]
                SMt  = self.MS_tension[eid][i]
                SMc  = self.MS_compression[eid][i]

                xxb = round(xxb,2)
                msg += '%-6i %6s %6i %4.2f ' %(eid,self.eType,nid,xxb)
                
                vals = [sMax,sMin,SMt,SMc]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10g ' %(val)
                    ###
                msg += '\n'
        #print msg
        return msg

class beamStrainObject(strainObject):
    """
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CBEAM' #{} # 'CBEAM/CONBEAM'

        self.code = [self.formatCode,self.sortCode,self.sCode]
        
        self.xxb = {}
        self.grids = {}
        self.sxc = {}
        self.sxd = {}
        self.sxe = {}
        self.sxf = {}
        self.smax = {}
        self.smin = {}
        self.MS_tension = {}
        self.MS_compression = {}
        
        if self.code in [[1,0,10]]:
            self.getLength1     = self.getLength1_format1_sort0
            self.getLength2     = self.getLength2_format1_sort0
            self.getLengthTotal = self.getLengthTotal_format1_sort0

            #self.isImaginary = False
            if dt is not None:
                self.addNewTransient = self.addNewTransient_format1_sort0
                self.addNewEid       = self.addNewEidTransient_format1_sort0
                self.add             = self.addTransient_format1_sort0
            else:
                self.addNewEid = self.addNewEid_format1_sort0
                self.add       = self.add_format1_sort0
            ###
        #elif self.code==[2,1,0]:
        #    self.getLength       = self.getLength_format1_sort0
        #    self.addNewTransient = self.addNewTransient_format2_sort1
        #    self.addNewEid       = self.addNewEidTransient_format2_sort1
        #    #self.isImaginary = True
        else:
            raise InvalidCodeError('beamStress - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt is not None:
            self.isTransient = True
            self.dt = self.nonlinearFactor
            self.addNewTransient()
        ###

    def getLengthTotal_format1_sort0(self):
        return 444  # 44+10*40   (11 nodes)

    def getLength1_format1_sort0(self):
        return (44,'iifffffffff')

    def getLength2_format1_sort0(self):
        return (40,'ifffffffff')

    def deleteTransient(self,dt):
        del self.sxc[dt]
        del self.sxd[dt]
        del self.sxe[dt]
        del self.sxf[dt]
        del self.smax[dt]
        del self.smin[dt]
        del self.MS_tension[dt]
        del self.MS_compression[dt]

    def getTransients(self):
        k = self.smax.keys()
        k.sort()
        return k

    def addNewTransient_format1_sort0(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        #print "addNewTransient_beam+1+0"
        if self.dt not in self.smax:
            self.sxc[self.dt] = {}
            self.sxd[self.dt] = {}
            self.sxe[self.dt] = {}
            self.sxf[self.dt] = {}
            self.smax[self.dt] = {}
            self.smin[self.dt] = {}
            self.MS_tension[self.dt]     = {}
            self.MS_compression[self.dt] = {}

    def addNewTransient_format2_sort1(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        raise Exception('not supported')
        #print self.dataCode
        if self.dt not in self.smax:
            self.axial[self.dt]     = {}
            self.torsion[self.dt]   = {}

    def addNewEid_format1_sort0(self,out):
        #print "Beam Stress addNewEid..."
        (eid,grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
        eid = (eid-self.deviceCode)/10
        #print "eid=%s grid=%s" %(eid,grid)
        assert eid >= 0
        #assert isinstance(eid,int)
        #assert isinstance(grid,int)
        self.grids[eid] = [grid]
        self.xxb[eid]  = [sd]
        self.sxc[eid] = [sxc]
        self.sxd[eid] = [sxd]
        self.sxe[eid] = [sxe]
        self.sxf[eid] = [sxf]
        self.smax[eid] = [smax]
        self.smin[eid] = [smin]
        self.MS_tension[eid] = [mst]
        self.MS_compression[eid] = [msc]
        return eid

    def addNewEid_format2_sort1(self,out):
        raise Exception('not supported')
        assert eid >= 0
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        self.axial[eid]      = [axialReal,axialImag]
        self.torsion[eid]    = [torsionReal,torsionImag]

    def addNewEidTransient_format1_sort0(self,out):
        #print "Beam Transient Stress addNewEid..."
        (eid,grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out

        eid = (eid-self.deviceCode)/10
        dt = self.dt
        assert eid  >= 0
        self.grids[eid] = [grid]
        self.xxb[eid] = [sd]
        self.sxc[dt][eid] = [sxc]
        self.sxd[dt][eid] = [sxd]
        self.sxe[dt][eid] = [sxe]
        self.sxf[dt][eid] = [sxf]
        self.smax[dt][eid] = [smax]
        self.smin[dt][eid] = [smin]
        self.MS_tension[dt][eid] = [mst]
        self.MS_compression[dt][eid] = [msc]
        return eid

    def addNewEidTransient_format2_sort1(self,out):
        raise Exception('not supported')
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        assert eid >= 0
        dt = self.dt
        self.axial[dt][eid]      = [axialReal,axialImag]
        self.torsion[dt][eid]    = [torsionReal,torsionImag]

    def add_format1_sort0(self,eid,out):
        #print "Beam Stress add..."
        (grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
        if grid:
            self.grids[eid].append(grid)
            self.xxb[eid].append(sd)
            self.sxc[eid].append(sxc)
            self.sxd[eid].append(sxd)
            self.sxe[eid].append(sxe)
            self.sxf[eid].append(sxf)
            self.smax[eid].append(smax)
            self.smin[eid].append(smin)
            self.MS_tension[eid].append(mst)
            self.MS_compression[eid].append(msc)
        ###

    def addTransient_format1_sort0(self,eid,out):
        #print "Beam Transient Stress add..."
        (grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
        dt = self.dt
        if grid:
            self.grids[eid].append(grid)
            self.xxb[eid].append(sd)
            self.smax[dt][eid].append(smax)
            self.smin[dt][eid].append(smin)
            self.MS_tension[dt][eid].append(mst)
            self.MS_compression[dt][eid].append(msc)
        ###

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)

        msg = header + ['                                  S T R A I N S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                        '                    STAT DIST/\n',
                        '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']

        for eid in sorted(self.smax):
            msg.append('0  %8i\n' %(eid))
            #print self.xxb[eid]
            for i,nid in enumerate(self.grids[eid]):
                #print i,nid
                xxb  = self.xxb[eid][i]
                sxc = self.sxc[eid][i]
                sxd = self.sxd[eid][i]
                sxe = self.sxe[eid][i]
                sxf = self.sxf[eid][i]
                sMax = self.smax[eid][i]
                sMin = self.smin[eid][i]
                SMt  = self.MS_tension[eid][i]
                SMc  = self.MS_compression[eid][i]
                (vals2,isAllZeros) = self.writeF06Floats13E([sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc])
                (sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc) = vals2
                msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' %(nid,xxb,sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc.strip()))
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        words = ['                                  S T R A I N S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                 '                    STAT DIST/\n',
                 '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        msg = []
        for dt,SMaxs in sorted(self.smax.items()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header + words
            for eid,Smax in sorted(SMaxs.items()):
                msg.append('0  %8i\n' %(eid))
                for i,nid in enumerate(self.grids[eid]):
                    xxb  = self.xxb[eid][i]
                    sxc = self.sxc[dt][eid][i]
                    sxd = self.sxd[dt][eid][i]
                    sxe = self.sxe[dt][eid][i]
                    sxf = self.sxf[dt][eid][i]
                    sMax = self.smax[dt][eid][i]
                    sMin = self.smin[dt][eid][i]
                    SMt  = self.MS_tension[dt][eid][i]
                    SMc  = self.MS_compression[dt][eid][i]
                    (vals2,isAllZeros) = self.writeF06Floats13E([sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc])
                    (sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc) = vals2
                    msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' %(nid,xxb,sxc,sxd,sxe,sxf,sMax,sMin,SMt,SMc.strip()))
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        return (''.join(msg),pageNum-1)

    def __reprTransient_format2_sort1__(self):
        raise Exception('not supported')
        msg = '---COMPLEX BEAM STRAINS---\n'
        msg += '%-10s %10s ' %('EID','eType')
        headers = ['axialReal','axialImag','torsionReal','torsionImag']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,axial in sorted(self.axial.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid in sorted(axial):
                axial   = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                msg += '%-6i %6s ' %(eid,self.eType)
                vals = axial + torsion # concatination
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10g ' %(val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
            ###
        return msg

    def __repr__(self):
        if   self.isTransient and self.code==[1,0,10]:
            return self.__reprTransient_format1_sort0__()
        elif self.code==[2,1,10]:
            return self.__reprTransient_format2_sort1__()

        msg = '---BEAM STRAINS---\n'
        msg += '%-6s %6s %6s %6s' %('EID','eType','NID','xxb')
        headers = ['sMax','sMin','MS_tension','MS_compression']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'
        #print "self.code = ",self.code
        for eid in sorted(self.smax):
            #print self.xxb[eid]
            for i,nid in enumerate(self.grids[eid]):
                #print i,nid
                xxb  = self.xxb[eid][i]
                sMax = self.smax[eid][i]
                sMin = self.smin[eid][i]
                SMt  = self.MS_tension[eid][i]
                SMc  = self.MS_compression[eid][i]

                xxb = round(xxb,2)
                msg += '%-6i %6s %6i %4.2f ' %(eid,self.eType,nid,xxb)
                
                vals = [sMax,sMin,SMt,SMc]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10.3e ' %(val)
                    ###
                msg += '\n'
        #print msg
        return msg
