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
from oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

class compositePlateStressObject(stressObject):
    """
    # sCode = 0
                    S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
    ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
      ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR

    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType  = {}

        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.o11    = {}
        self.o22    = {}
        self.t12    = {}
        self.t1z    = {}
        self.t2z    = {}
        self.angle  = {}
        self.majorP = {}
        self.minorP = {}
        #self.fiberCurvature = {}
        self.ovmShear       = {}
        #print "self.dataCode = ",self.dataCode
        #if self.isVonMisesStress():
        if self.code == [1,0,0]:
            assert self.isVonMises()==False,'isVonMises=%s' %(self.isVonMises())
            #self.isVonMises      = True
        else:
            raise InvalidCodeError('compositePlateStress - get the format/sort/stressCode=%s' %(self.code))
        ###

        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.add       = self.addTransient
            self.addNewEid = self.addNewEidTransient
        ###

    def addF06Data(self,data,transient,eType):
        """
                       S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
        ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
          ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
            151    1  -1.02406E+04  4.18348E+05  4.14359E+02   -8.62021E+00  1.86352E+04   89.94  4.18348E+05 -1.02410E+04  2.14295E+05
        """
        if transient is None:
            #for line in data:
                #print line
            #sys.exit()
            for line in data:
                #print line
                if eType=='CQUAD4':
                    (eid,iLayer,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovmShear) = line
                    #if nid=='CEN/4': nid='C'
                    if eid not in self.eType:
                        self.eType[eid] = 'CQUAD4'
                        self.o11[eid]    = [o11]
                        self.o22[eid]    = [o22]
                        self.t12[eid]    = [t12]
                        self.t1z[eid]    = [t1z]
                        self.t2z[eid]    = [t2z]
                        self.angle[eid]  = [angle]
                        self.majorP[eid] = [majorP]
                        self.minorP[eid] = [minorP]
                        self.ovmShear[eid] = [ovmShear]
                    else:
                        self.o11[eid].append(o11)
                        self.o22[eid].append(o22)
                        self.t12[eid].append(t12)
                        self.t1z[eid].append(t1z)
                        self.t2z[eid].append(t2z)
                        self.angle[eid].append(angle)
                        self.majorP[eid].append(majorP)
                        self.minorP[eid].append(minorP)
                        self.ovmShear[eid].append(ovmShear)
                else:
                    raise NotImplementedError('line=%s not supported...' %(line))
            return
        #for line in data:
            #print line
        raise NotImplementedError('transient results not supported')

    def deleteTransient(self,dt):
        del self.o11[dt]
        del self.o22[dt]
        del self.t12[dt]
        del self.t1z[dt]
        del self.t2z[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.ovmShear[dt]

    def getTransients(self):
        k = self.o11.keys()
        k.sort()
        return k

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        #self.fiberDistance[self.dt] = {}
        if self.dt not in self.o11:
            self.o11[self.dt]    = {}
            self.o22[self.dt]    = {}
            self.t12[self.dt]    = {}
            self.t1z[self.dt]    = {}
            self.t2z[self.dt]    = {}
            self.angle[self.dt]  = {}
            self.majorP[self.dt] = {}
            self.minorP[self.dt] = {}
            self.ovmShear[self.dt]    = {}

    def addNewEid(self,eType,eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        msg = "eid=%s eType=%s o11=%g o22=%g t12=%g t1z=%g t2z=%g \nangle=%g major=%g minor=%g vm=%g" %(eid,eType,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm)
        if eid in self.o11:
            return self.add(eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm)
        assert eid not in self.o11,msg+'\n  o11=%s eType=%s code=%s' %(self.o11[eid],self.eType[eid],self.dataCode)

        self.eType[eid]  = eType
        self.o11[eid]    = [o11]
        self.o22[eid]    = [o22]
        self.t12[eid]    = [t12]
        self.t1z[eid]    = [t1z]
        self.t2z[eid]    = [t2z]
        self.angle[eid]  = [angle]
        self.majorP[eid] = [majorP]
        self.minorP[eid] = [minorP]
        self.ovmShear[eid]    = [ovm]
        #print msg
        #if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        dt = self.dt
        if eid in self.o11[dt]:
            return self.add(eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm)
        assert eid not in self.o11[dt]
        self.eType[eid]  = eType
        self.o11[dt][eid]    = [o11]
        self.o22[dt][eid]    = [o22]
        self.t12[dt][eid]    = [t12]
        self.t1z[dt][eid]    = [t1z]
        self.t2z[dt][eid]    = [t2z]
        self.angle[dt][eid]  = [angle]
        self.majorP[dt][eid] = [majorP]
        self.minorP[dt][eid] = [minorP]
        self.ovmShear[dt][eid]    = [ovm]
        msg = "eid=%s o11=%g o22=%g t12=%g t1z=%g t2z=%g \nangle=%g major=%g minor=%g vm=%g" %(eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add(self,eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm):
        #print "***add"
        msg = "eid=%s o11=%g o22=%g t12=%g t1z=%g t2z=%g \nangle=%g major=%g minor=%g vm=%g" %(eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm)
        #print msg
        #print self.o11
        self.o11[eid].append(o11)
        self.o22[eid].append(o22)
        self.t12[eid].append(t12)
        self.t1z[eid].append(t1z)
        self.t2z[eid].append(t2z)
        self.angle[eid].append(angle)
        self.majorP[eid].append(majorP)
        self.minorP[eid].append(minorP)
        self.ovmShear[eid].append(ovm)
        #if nodeID==0: raise Exception(msg)

    def addTransient(self,eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm):
        #print "***add"
        msg = "eid=%s o11=%g o22=%g t12=%g t1z=%g t2z=%g \nangle=%g major=%g minor=%g vm=%g" %(eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm)
        #print msg
        #print self.o11
        dt = self.dt
        self.o11[dt][eid].append(o11)
        self.o22[dt][eid].append(o22)
        self.t12[dt][eid].append(t12)
        self.t1z[dt][eid].append(t1z)
        self.t2z[dt][eid].append(t2z)
        self.angle[dt][eid].append(angle)
        self.majorP[dt][eid].append(majorP)
        self.minorP[dt][eid].append(minorP)
        self.ovmShear[dt][eid].append(ovm)
        #if nodeID==0: raise Exception(msg)

    def getHeaders(self):
        headers = ['o11','o22','t12','t1z','t2z']
        if self.isVonMises:
            headers.append('oVonMises')
        else:
            headers.append('maxShear')
        return headers

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)

        if self.isVonMises():
            von   = 'VON'
            mises = 'MISES'
        else:
            von   = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' %(von),
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n'  %(mises)]

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes or 'QUAD4LC' in eTypes:
            quadMsg = header+['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n']+words
            isQuad = True
        else:
            quadMsg = []
            isQuad = False

        if 'CTRIA3' in eTypes or 'TRIA3LC' in eTypes:
            isTri = True
            triMsg = header+['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n']+words
        else:
            isTri = False
            triMsg = []
        ###
        for eid,o11s in sorted(self.o11.items()):
            out = ''
            eType = self.eType[eid]
            for iLayer in range(len(o11s)):
                o11 = self.o11[eid][iLayer]
                o22 = self.o22[eid][iLayer]
                t12 = self.t12[eid][iLayer]
                t1z = self.t1z[eid][iLayer]
                t2z = self.t2z[eid][iLayer]

                angle = self.angle[eid][iLayer]
                major = self.majorP[eid][iLayer]
                minor = self.minorP[eid][iLayer]
                ovm   = self.ovmShear[eid][iLayer]
                (vals2,isAllZeros) = self.writeF06Floats12E([o11,o22,t12,t1z,t2z,major,minor,ovm])
                [o11,o22,t12,t1z,t2z,major,minor,ovm] = vals2
                out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' %(eid,iLayer+1,o11,o22,t12,t1z,t2z,angle,major,minor,ovm)

            if eType in ['CQUAD4','QUAD4LC']:
                quadMsg.append(out)
            elif eType in ['CTRIA3','TRIA3LC']:
                triMsg.append(out)
            #else:
            #    raise NotImplementedError('eType = |%r|' %(eType)) # CQUAD8LC
            ###
        ###
        if isQuad:
            quadMsg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        if isTri:
            triMsg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1

        msg = ''.join(quadMsg+triMsg)
        return (msg,pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        if self.isVonMises():
            von   = 'VON'
            mises = 'MISES'
        else:
            von   = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' %(von),
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n'  %(mises)]

        eTypes = self.eType.values()        
        if 'CQUAD4' in eTypes or 'QUAD4LC' in eTypes:
            quadWords = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n']+words
            isQuad = True
        else:
            quadWords = []
            isQuad = False

        if 'CTRIA3' in eTypes or 'TRIA3LC' in eTypes:
            isTri = True
            triWords = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n']+words
        else:
            isTri = False
            triWords = []
        ###
        
        msg = []
        for dt,O11s in sorted(self.o11.items()):
            quadMsg = []
            triMsg = []
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            if isQuad:
                quadMsg = header + quadWords
            if isTri:
                triMsg = header + triWords

            for eid,o11s in sorted(O11s.items()):
                out = ''
                eType = self.eType[eid]
                for iLayer in range(len(o11s)):
                    o11 = self.o11[dt][eid][iLayer]
                    o22 = self.o22[dt][eid][iLayer]
                    t12 = self.t12[dt][eid][iLayer]
                    t1z = self.t1z[dt][eid][iLayer]
                    t2z = self.t2z[dt][eid][iLayer]

                    angle = self.angle[dt][eid][iLayer]
                    major = self.majorP[dt][eid][iLayer]
                    minor = self.minorP[dt][eid][iLayer]
                    ovm   = self.ovmShear[dt][eid][iLayer]
                    (vals2,isAllZeros) = self.writeF06Floats12E([o11,o22,t12,t1z,t2z,major,minor,ovm])
                    [o11,o22,t12,t1z,t2z,major,minor,ovm] = vals2
                    out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' %(eid,iLayer+1,o11,o22,t12,t1z,t2z,angle,major,minor,ovm)

                if eType in ['CQUAD4','QUAD4LC']:
                    quadMsg.append(out)
                elif eType in ['CTRIA3','TRIA3LC']:
                    triMsg.append(out)
                #else:
                #    raise NotImplementedError('eType = |%r|' %(eType)) # CQUAD8LC
                ###
            ###
            if isQuad:
                quadMsg.append(pageStamp+str(pageNum)+'\n')
                pageNum+=1
            if isTri:
                triMsg.append(pageStamp+str(pageNum)+'\n')
                pageNum+=1
            ###
            msg += quadMsg+triMsg
        ###
        return (''.join(msg),pageNum-1)

    def __repr__(self):
        return self.writeF06(['','',''],'PAGE ',1)[0]
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---COMPOSITE PLATE STRESS---\n'
        msg += '%-6s %8s %8s ' %('EID','eType','iLayer')
        headers = self.getHeaders()
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for eid,o11s in sorted(self.o11.items()):
            eType = self.eType[eid]
            for iLayer in range(len(o11s)):
                o11 = self.o11[eid][iLayer]
                o22 = self.o22[eid][iLayer]
                t12 = self.t12[eid][iLayer]
                t1z = self.t1z[eid][iLayer]
                t2z = self.t2z[eid][iLayer]

                angle = self.angle[eid][iLayer]
                major = self.majorP[eid][iLayer]
                minor = self.minorP[eid][iLayer]
                ovm   = self.ovmShear[eid][iLayer]

                msg += '%-6i %8s %8s ' %(eid,eType,iLayer+1,)
                vals = [o11,o22,t12,t1z,t2z,ovm]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10i ' %(val)
                    ###
                msg += '\n'
            ###
        ###
        return msg

    def __reprTransient__(self):
        msg = '---COMPOSITE PLATE STRESS---\n'
        msg += '%-6s %8s %8s ' %('EID','eType','iLayer')
        headers = self.getHeaders()
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,O11s in sorted(self.o11.items()):
            msg += "dt = %s\n" %(dt)
            for eid,o11s in sorted(O11s.items()):
                eType = self.eType[eid]
                for iLayer in range(len(o11s)):
                    o11 = self.o11[dt][eid][iLayer]
                    o22 = self.o22[dt][eid][iLayer]
                    t12 = self.t12[dt][eid][iLayer]
                    t1z = self.t1z[dt][eid][iLayer]
                    t2z = self.t2z[dt][eid][iLayer]

                    angle = self.angle[dt][eid][iLayer]
                    major = self.majorP[dt][eid][iLayer]
                    minor = self.minorP[dt][eid][iLayer]
                    ovm   = self.ovmShear[dt][eid][iLayer]

                    msg += '%-6i %8s %8s ' %(eid,eType,iLayer+1,)
                    vals = [o11,o22,t12,t1z,t2z,ovm]
                    for val in vals:
                        if abs(val)<1e-6:
                            msg += '%10s ' %('0')
                        else:
                            msg += '%10i ' %(val)
                        ###
                    msg += '\n'
                ###
        ###
        return msg

class compositePlateStrainObject(strainObject):
    """
    ???
    ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
      ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)

        self.eType  = {}
        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.e11    = {}
        self.e22    = {}
        self.e12    = {}
        self.e1z    = {}
        self.e2z    = {}
        self.angle  = {}
        self.majorP = {}
        self.minorP = {}
        
        if self.code == [1,0,14]:
            self.evmShear   = {}
            assert self.isVonMises()==False
        else:
            raise InvalidCodeError('compositePlateStrain - get the format/sort/stressCode=%s' %(self.code))
        ###

        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.add       = self.addTransient
            self.addNewEid = self.addNewEidTransient
        ###

    def deleteTransient(self,dt):
        del self.e11[dt]
        del self.e22[dt]
        del self.e12[dt]
        del self.e1z[dt]
        del self.e2z[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.evmShear[dt]

    def getTransients(self):
        k = self.e11.keys()
        k.sort()
        return k

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        #self.fiberDistance[self.dt] = {}
        if self.dt not in self.e11:
            self.e11[self.dt]    = {}
            self.e22[self.dt]    = {}
            self.e12[self.dt]    = {}
            self.e1z[self.dt]    = {}
            self.e2z[self.dt]    = {}
            self.angle[self.dt]  = {}
            self.majorP[self.dt] = {}
            self.minorP[self.dt] = {}
            self.evmShear[self.dt]    = {}

    def addNewEid(self,eType,eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        if eid in self.e11:
            return self.add(eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm)
        assert eid not in self.e11
        assert isinstance(eid,int)
        self.eType[eid]  = eType
        self.e11[eid]    = [e11]
        self.e22[eid]    = [e22]
        self.e12[eid]    = [e12]
        self.e1z[eid]    = [e1z]
        self.e2z[eid]    = [e2z]
        self.angle[eid]  = [angle]
        self.majorP[eid] = [majorP]
        self.minorP[eid] = [minorP]
        self.evmShear[eid] = [evm]
        msg = "eid=%s e11=%g e22=%g e12=%g e1z=%g e2z=%g \nangle=%g major=%g minor=%g vm=%g" %(eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        dt = self.dt
        assert eid not in self.e11[dt]
        assert isinstance(eid,int)
        self.eType[eid]  = eType
        self.e11[dt][eid]    = [e11]
        self.e22[dt][eid]    = [e22]
        self.e12[dt][eid]    = [e12]
        self.e1z[dt][eid]    = [e1z]
        self.e2z[dt][eid]    = [e2z]
        self.angle[dt][eid]  = [angle]
        self.majorP[dt][eid] = [majorP]
        self.minorP[dt][eid] = [minorP]
        self.evmShear[dt][eid] = [evm]
        msg = "eid=%s e11=%g e22=%g e12=%g e1z=%g e2z=%g \nangle=%g major=%g minor=%g vm=%g" %(eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add(self,eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm):
        #print "***add"
        msg = "eid=%s e11=%g e22=%g e12=%g e1z=%g e2z=%g \nangle=%g major=%g minor=%g vm=%g" %(eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm)
        #print msg
        #print self.o11
        self.e11[eid].append(e11)
        self.e22[eid].append(e22)
        self.e12[eid].append(e12)
        self.e1z[eid].append(e1z)
        self.e2z[eid].append(e2z)
        self.angle[eid].append(angle)
        self.majorP[eid].append(majorP)
        self.minorP[eid].append(minorP)
        self.evmShear[eid].append(evm)
        #if nodeID==0: raise Exception(msg)

    def addTransient(self,eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm):
        #print "***add"
        msg = "eid=%s e11=%g e22=%g e12=%g e1z=%g e2z=%g \nangle=%g major=%g minor=%g vm=%g" %(eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm)
        #print msg
        #print self.o11
        dt = self.dt
        self.e11[dt][eid].append(e11)
        self.e22[dt][eid].append(e22)
        self.e12[dt][eid].append(e12)
        self.e1z[dt][eid].append(e1z)
        self.e2z[dt][eid].append(e2z)
        self.angle[dt][eid].append(angle)
        self.majorP[dt][eid].append(majorP)
        self.minorP[dt][eid].append(minorP)
        self.evmShear[dt][eid].append(evm)
        #if nodeID==0: raise Exception(msg)

    def getHeaders(self):
        headers = ['e11','e22','e12','e1z','e2z']
        if self.isVonMises:
            headers.append('eVonMises')
        else:
            headers.append('maxShear')
        return headers

    def writeF06(self,header,pageStamp,pageNum=1):  # @todo no transient
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)

        if self.isVonMises():
            von   = 'VON'
            mises = 'MISES'
        else:
            von   = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY   STRAINS IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR   STRAINS  PRINCIPAL  STRAINS (ZERO SHEAR)      %s\n' %(von),
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n'  %(mises)]

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes or 'QUAD4LC' in eTypes:
            quadMsg = header+['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n']+words
            isQuad = True
        else:
            quadMsg = []
            isQuad = False

        if 'CTRIA3' in eTypes or 'TRIA3LC' in eTypes:
            isTri = True
            triMsg = header+['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n']+words
        else:
            isTri = False
            triMsg = []
        ###
        for eid,e11s in sorted(self.e11.items()):
            out = ''
            eType = self.eType[eid]
            for iLayer in range(len(e11s)):
                e11 = self.e11[eid][iLayer]
                e22 = self.e22[eid][iLayer]
                e12 = self.e12[eid][iLayer]
                e1z = self.e1z[eid][iLayer]
                e2z = self.e2z[eid][iLayer]

                angle = self.angle[eid][iLayer]
                major = self.majorP[eid][iLayer]
                minor = self.minorP[eid][iLayer]
                evm   = self.evmShear[eid][iLayer]
                
                (vals2,isAllZeros) = self.writeF06Floats12E([e11,e22,e12,e1z,e2z,major,minor,evm])
                [e11,e22,e12,e1z,e2z,major,minor,evm] = vals2
                out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' %(eid,iLayer+1,e11,e22,e12,e1z,e2z,angle,major,minor,evm)

            if eType in ['CQUAD4','QUAD4LC']:
                quadMsg.append(out)
            elif eType in ['CTRIA3','TRIA3LC']:
                triMsg.append(out)
            #else:
            #    raise NotImplementedError('eType = |%r|' %(eType)) # CQUAD8LC
            ###
        ###
        if isQuad:
            quadMsg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        if isTri:
            triMsg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1

        msg = ''.join(quadMsg+triMsg)
        return (msg,pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        if self.isVonMises():
            von   = 'VON'
            mises = 'MISES'
        else:
            von   = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' %(von),
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n'  %(mises)]

        eTypes = self.eType.values()        
        if 'CQUAD4' in eTypes or 'QUAD4LC' in eTypes:
            quadWords = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n']+words
            isQuad = True
        else:
            quadWords = []
            isQuad = False

        if 'CTRIA3' in eTypes or 'TRIA3LC' in eTypes:
            isTri = True
            triWords = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n']+words
        else:
            isTri = False
            triWords = []
        ###
        
        msg = []
        for dt,e11s in sorted(self.e11.items()):
            quadMsg = []
            triMsg = []
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            if isQuad:
                quadMsg = header + quadWords
            if isTri:
                triMsg = header + triWords

            for eid,e11s in sorted(e11s.items()):
                out = ''
                eType = self.eType[eid]
                for iLayer in range(len(e11s)):
                    e11 = self.e11[dt][eid][iLayer]
                    e22 = self.e22[dt][eid][iLayer]
                    e12 = self.e12[dt][eid][iLayer]
                    e1z = self.e1z[dt][eid][iLayer]
                    e2z = self.e2z[dt][eid][iLayer]

                    angle = self.angle[dt][eid][iLayer]
                    major = self.majorP[dt][eid][iLayer]
                    minor = self.minorP[dt][eid][iLayer]
                    evm   = self.evmShear[dt][eid][iLayer]
                    (vals2,isAllZeros) = self.writeF06Floats12E([e11,e22,e12,e1z,e2z,major,minor,evm])
                    [e11,e22,e12,e1z,e2z,major,minor,evm] = vals2
                    out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' %(eid,iLayer+1,e11,e22,e12,e1z,e2z,angle,major,minor,evm)

                if eType in ['CQUAD4','QUAD4LC']:
                    quadMsg.append(out)
                elif eType in ['CTRIA3','TRIA3LC']:
                    triMsg.append(out)
                else:
                    raise NotImplementedError('eType = |%r|' %(eType))
                ###
            ###
            if isQuad:
                quadMsg.append(pageStamp+str(pageNum)+'\n')
                pageNum+=1
            if isTri:
                triMsg.append(pageStamp+str(pageNum)+'\n')
                pageNum+=1
            ###
            msg += quadMsg+triMsg
        ###
        return (''.join(msg),pageNum-1)

    def __repr__(self):
        return self.writeF06(['','',''],'PAGE ',1)[0]
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---COMPOSITE PLATE STAIN---\n'
        msg += '%-6s %8s %8s ' %('EID','eType','iLayer')
        headers = self.getHeaders()
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for eid,e11s in sorted(self.e11.items()):
            eType = self.eType[eid]
            for iLayer in range(len(e11s)):
                e11 = self.e11[eid][iLayer]
                e22 = self.e22[eid][iLayer]
                e12 = self.e12[eid][iLayer]
                e1z = self.e1z[eid][iLayer]
                e2z = self.e2z[eid][iLayer]

                angle = self.angle[eid][iLayer]
                major = self.majorP[eid][iLayer]
                minor = self.minorP[eid][iLayer]
                evm   = self.evmShear[eid][iLayer]

                msg += '%-6i %8s %8s ' %(eid,eType,iLayer+1,)
                vals = [e11,e22,e12,e1z,e2z,evm]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10.3g ' %(val)
                    ###
                msg += '\n'
            ###
        ###
        return msg

    def __reprTransient__(self):
        msg = '---COMPOSITE PLATE STAIN---\n'
        headers = self.getHeaders()
        msg += '%-6s %8s %8s ' %('EID','eType','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,E11s in sorted(self.e11.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid,e11s in sorted(E11s.items()):
                eType = self.eType[eid]
                for iLayer in range(len(e11s)):
                    e11 = self.e11[dt][eid][iLayer]
                    e22 = self.e22[dt][eid][iLayer]
                    e12 = self.e12[dt][eid][iLayer]
                    e1z = self.e1z[dt][eid][iLayer]
                    e2z = self.e2z[dt][eid][iLayer]

                    angle = self.angle[dt][eid][iLayer]
                    major = self.majorP[dt][eid][iLayer]
                    minor = self.minorP[dt][eid][iLayer]
                    evm   = self.evmShear[dt][eid][iLayer]

                    msg += '%-6i %8s %8s ' %(eid,eType,iLayer+1,)
                    vals = [e11,e22,e12,e1z,e2z,evm]
                    for val in vals:
                        if abs(val)<1e-6:
                            msg += '%10s ' %('0')
                        else:
                            msg += '%10.3g ' %(val)
                        ###
                    msg += '\n'
                ###
            ###
        ###
        return msg
