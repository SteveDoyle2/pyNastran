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
from .oes_objects import stressObject,strainObject
from pyNastran.op2.op2Errors import *


class plateStressObject(stressObject):
    """
    ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 
      ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
          6    CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  4.227345E+03
                       1.250000E-01   5.406062E+02  1.201854E+04 -4.174177E+01   -89.7916   1.201869E+04  5.404544E+02  5.739119E+03


                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  
    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)          MAX  
      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR         SHEAR   
          6    CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  4.227345E+03
                       1.250000E-01   5.406062E+02  1.201854E+04 -4.174177E+01   -89.7916   1.201869E+04  5.404544E+02  5.739119E+03
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = {}

        #self.appendDataMember('sCodes','sCode')
        #print "self.sCode = ",self.sCode
        self.code = [self.formatCode,self.sortCode,self.sCode]

        self.fiberCurvature = {}
        self.oxx    = {}
        self.oyy    = {}
        self.txy    = {}
        self.angle  = {}
        self.majorP = {}
        self.minorP = {}
        self.ovmShear = {}

        #print "%%%%%dt = ",dt
        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.add        = self.addTransient
            self.addNewEid  = self.addNewEidTransient
            self.addNewNode = self.addNewNodeTransient
        #else:
        #    self.dt = None
        ###

    def addF06Data(self,data,transient):
        if transient is None:
            eType = data[0][0]
            for line in data:
                if eType=='CTRIA3':
                    (eType,eid,f1,ox1,oy1,txy1,angle1,o11,o21,ovm1,
                               f2,ox2,oy2,txy2,angle2,o12,o22,ovm2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {'C':[f1,f2]}
                    self.oxx[eid]      = {'C':[ox1,ox2]}
                    self.oyy[eid]      = {'C':[oy1,oy2]}
                    self.txy[eid]      = {'C':[txy1,txy2]}
                    self.angle[eid]    = {'C':[angle1,angle2]}
                    self.majorP[eid]   = {'C':[o11,o12]}
                    self.minorP[eid]   = {'C':[o21,o22]}
                    self.ovmShear[eid] = {'C':[ovm1,ovm2]}
                elif eType=='CQUAD4':
                    #assert len(line)==19,'len(line)=%s' %(len(line))
                    if len(line)==19: # Centroid - bilinear
                        (eType,eid,nid,f1,ox1,oy1,txy1,angle1,o11,o21,ovm1,
                                       f2,ox2,oy2,txy2,angle2,o12,o22,ovm2) = line
                        if nid=='CEN/4': nid='C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid:[f1,f2]}
                        self.oxx[eid]      = {nid:[ox1,ox2]}
                        self.oyy[eid]      = {nid:[oy1,oy2]}
                        self.txy[eid]      = {nid:[txy1,txy2]}
                        self.angle[eid]    = {nid:[angle1,angle2]}
                        self.majorP[eid]   = {nid:[o11,o12]}
                        self.minorP[eid]   = {nid:[o21,o22]}
                        self.ovmShear[eid] = {nid:[ovm1,ovm2]}
                    elif len(line)==18: # Centroid
                        (eType,eid,f1,ox1,oy1,txy1,angle1,o11,o21,ovm1,
                                   f2,ox2,oy2,txy2,angle2,o12,o22,ovm2) = line
                        nid='C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid:[f1,f2]}
                        self.oxx[eid]      = {nid:[ox1,ox2]}
                        self.oyy[eid]      = {nid:[oy1,oy2]}
                        self.txy[eid]      = {nid:[txy1,txy2]}
                        self.angle[eid]    = {nid:[angle1,angle2]}
                        self.majorP[eid]   = {nid:[o11,o12]}
                        self.minorP[eid]   = {nid:[o21,o22]}
                        self.ovmShear[eid] = {nid:[ovm1,ovm2]}
                    elif len(line)==17: ## Bilinear
                        #print line
                        (nid,f1,ox1,oy1,txy1,angle1,o11,o21,ovm1,
                             f2,ox2,oy2,txy2,angle2,o12,o22,ovm2) = line
                        self.fiberCurvature[eid][nid] = [f1,f2]
                        self.oxx[eid][nid]      = [ox1,ox2]
                        self.oyy[eid][nid]      = [oy1,oy2]
                        self.txy[eid][nid]      = [txy1,txy2]
                        self.angle[eid][nid]    = [angle1,angle2]
                        self.majorP[eid][nid]   = [o11,o12]
                        self.minorP[eid][nid]   = [o21,o22]
                        self.ovmShear[eid][nid] = [ovm1,ovm2]
                    else:
                        assert len(line)==19,'len(line)=%s' %(len(line))
                        #raise NotImplementedError()
                    ###
                else:
                    raise NotImplementedError('line=%s not supported...' %(line))
            return
        #for line in data:
            #print line
        raise NotImplementedError('transient results not supported')

    def deleteTransient(self,dt):
        del self.fiberCurvature[dt]
        del self.oxx[dt]
        del self.oyy[dt]
        del self.txy[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.ovmShear[dt]

    def getTransients(self):
        k = self.oxx.keys()
        k.sort()
        return k

    def addNewTransient(self):
        """
        initializes the transient variables
        """
        if self.dt not in self.oxx:
            self.fiberCurvature[self.dt] = {}
            self.oxx[self.dt]      = {}
            self.oyy[self.dt]      = {}
            self.txy[self.dt]      = {}
            self.angle[self.dt]    = {}
            self.majorP[self.dt]   = {}
            self.minorP[self.dt]   = {}
            self.ovmShear[self.dt] = {}

    def addNewEid(self,eType,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***addNewEid..."
        if eid in self.oxx:
            return self.add(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)

        assert eid not in self.oxx
        #print self.oxx
        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [fd]}
        assert isinstance(eid,int)
        self.oxx[eid]    = {nodeID: [oxx]}
        self.oyy[eid]    = {nodeID: [oyy]}
        self.txy[eid]    = {nodeID: [txy]}
        self.angle[eid]  = {nodeID: [angle]}
        self.majorP[eid] = {nodeID: [majorP]}
        self.minorP[eid] = {nodeID: [minorP]}
        self.ovmShear[eid] = {nodeID: [ovm]}
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***addNewEidTransient..."
        dt = self.dt
        #msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(dt,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g major=%g vm=%g" %(dt,eid,nodeID,fd,oxx,majorP,ovm)
        #print msg
        #if eid in self.ovmShear[dt]:
        #    return self.add(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)

        if eid in self.ovmShear[dt]:  # SOL200, erase the old result
            nid = nodeID
            #msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s major=%s vm=%s" %(dt,eid,nodeID,str(self.fiberCurvature[dt][eid][nid]),str(self.oxx[dt][eid][nid]),str(self.majorP[dt][eid][nid]),str(self.ovmShear[dt][eid][nid]))
            self.deleteTransient(dt)
            self.addNewTransient()
        
        assert isinstance(eid,int)
        self.eType[eid] = eType
        self.fiberCurvature[dt][eid] = {nodeID: [fd]}
        self.oxx[dt][eid]      = {nodeID: [oxx]}
        self.oyy[dt][eid]      = {nodeID: [oyy]}
        self.txy[dt][eid]      = {nodeID: [txy]}
        self.angle[dt][eid]    = {nodeID: [angle]}
        self.majorP[dt][eid]   = {nodeID: [majorP]}
        self.minorP[dt][eid]   = {nodeID: [minorP]}
        self.ovmShear[dt][eid] = {nodeID: [ovm]}
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***add"
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        assert isinstance(eid,int)
        self.fiberCurvature[eid][nodeID].append(fd)
        self.oxx[eid][nodeID].append(oxx)
        self.oyy[eid][nodeID].append(oyy)
        self.txy[eid][nodeID].append(txy)
        self.angle[eid][nodeID].append(angle)
        self.majorP[eid][nodeID].append(majorP)
        self.minorP[eid][nodeID].append(minorP)
        self.ovmShear[eid][nodeID].append(ovm)
        if nodeID==0: raise Exception(msg)

    def addTransient(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***addTransient"
        dt = self.dt
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(dt,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #print self.oxx
        #print self.fiberCurvatrure
        assert eid is not None
        self.fiberCurvature[dt][eid][nodeID].append(fd)
        self.oxx[dt][eid][nodeID].append(oxx)
        self.oyy[dt][eid][nodeID].append(oyy)
        self.txy[dt][eid][nodeID].append(txy)
        self.angle[dt][eid][nodeID].append(angle)
        self.majorP[dt][eid][nodeID].append(majorP)
        self.minorP[dt][eid][nodeID].append(minorP)
        self.ovmShear[dt][eid][nodeID].append(ovm)
        if nodeID==0: raise Exception(msg)

    def addNewNode(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***addNewNode"
        #print self.oxx
        assert eid is not None
        #assert nodeID not in self.oxx[eid]
        self.fiberCurvature[eid][nodeID] = [fd]
        self.oxx[eid][nodeID]    = [oxx]
        self.oyy[eid][nodeID]    = [oyy]
        self.txy[eid][nodeID]    = [txy]
        self.angle[eid][nodeID]  = [angle]
        self.majorP[eid][nodeID] = [majorP]
        self.minorP[eid][nodeID] = [minorP]
        self.ovmShear[eid][nodeID]    = [ovm]
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def addNewNodeTransient(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***addNewNodeTransient"
        #print self.oxx
        assert eid is not None
        dt = self.dt
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #assert nodeID not in self.oxx[dt][eid]
        self.fiberCurvature[dt][eid][nodeID] = [fd]
        self.oxx[dt][eid][nodeID]    = [oxx]
        self.oyy[dt][eid][nodeID]    = [oyy]
        self.txy[dt][eid][nodeID]    = [txy]
        self.angle[dt][eid][nodeID]  = [angle]
        self.majorP[dt][eid][nodeID] = [majorP]
        self.minorP[dt][eid][nodeID] = [minorP]
        self.ovmShear[dt][eid][nodeID]    = [ovm]
        if nodeID==0: raise Exception(msg)

    def getHeaders(self):
        if self.isFiberDistance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']
        headers += ['oxx','oyy','txy','majorP','minorP']
        if self.isVonMises():
            headers.append('oVonMises')
        else:
            headers.append('maxShear')
        return headers

    def __reprTransient__(self):
        msg = '---ISOTROPIC PLATE STRESS---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' %('EID','eType','nodeID','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        #print self.oxx.keys()
        for dt,oxxs in sorted(self.oxx.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid,oxxNodes in sorted(oxxs.items()):
                eType = self.eType[eid]
                for nid in sorted(oxxNodes):
                    for iLayer in range(len(self.oxx[dt][eid][nid])):
                        fd    = self.fiberCurvature[dt][eid][nid][iLayer]
                        oxx   = self.oxx[   dt][eid][nid][iLayer]
                        oyy   = self.oyy[   dt][eid][nid][iLayer]
                        txy   = self.txy[   dt][eid][nid][iLayer]
                        angle = self.angle[ dt][eid][nid][iLayer]
                        major = self.majorP[dt][eid][nid][iLayer]
                        minor = self.minorP[dt][eid][nid][iLayer]
                        ovm = self.ovmShear[dt][eid][nid][iLayer]

                        msg += '%-6i %6s %8s %7s %10g ' %(eid,eType,nid,iLayer+1,fd)
                        vals = [oxx,oyy,txy,major,minor,ovm]
                        for val in vals:
                            if abs(val)<1e-6:
                                msg += '%10s ' %('0')
                            else:
                                try:
                                    msg += '%10i ' %(int(val))
                                except:
                                    print("bad val = ",val)
                                    raise
                            ###
                        msg += '\n'
                    ###
                ###
            ###
        ###
        return msg

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' %(vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                          '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' %(vonMises)]

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            qkey = eTypes.index('CQUAD4')
            kkey = self.eType.keys()[qkey]
            ekey = self.oxx[kkey].keys()
            isBilinear=True
            quadMsg = header+['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n']+quadMsgTemp
            if len(ekey)==1:
                isBilinear=False
                quadMsg = header+['                           S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n \n']+triMsgTemp
            isQuad = True
        else:
            quadMsg = []
            isQuad = False

        if 'CTRIA3' in eTypes:
            isTri = True
            triMsg = header+['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']+triMsgTemp
        else:
            isTri = False
            triMsg = []

        for eid,oxxNodes in sorted(self.oxx.items()):
            eType = self.eType[eid]
            if eType in 'CQUAD4':
                if isBilinear:
                    out = self.writeF06_Quad4_Bilinear(eid,oxxNodes)
                else:
                    out = self.writeF06_Tri3(eid,oxxNodes)
                quadMsg.append(out)
            elif eType in 'CTRIA3':
                out = self.writeF06_Tri3(eid,oxxNodes)
                triMsg.append(out)
            #else:
            #    raise NotImplementedError('eType = |%r|' %(eType)) # CQUAD8, CTRIA6
            ###
        ###
        if isQuad:
            quadMsg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        if isTri:
            triMsg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1

        #print "quadMsg = ",quadMsg
        #print "triMsg = ",triMsg
        return (''.join(quadMsg+triMsg),pageNum-1)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance():
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' %(vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                           '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)\n',
                          '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' %(vonMises)]

        #print self.oxx
        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            ElemKey = eTypes.index('CQUAD4')
            #print qkey
            eid = self.eType.keys()[ElemKey]
            #print "self.oxx = ",self.oxx
            #print "eid=%s" %(eid)
            dt = self.oxx.keys()[0]
            #print "dt=%s" %(dt)
            nLayers = len(self.oxx[dt][eid])
            #print "elementKeys = ",elementKeys
            isBilinear=True
            quadMsg2 = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n']+quadMsgTemp
            if nLayers==1:
                isBilinear=False
                quadMsg2 = ['                           S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n \n']+triMsgTemp
            isQuad = True
        else:
            quadMsg2 = []
            isQuad = False

        if 'CTRIA3' in eTypes:
            isTri = True
            triMsg2 = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']+triMsgTemp
        else:
            isTri = False
            triMsg2 = []
        
        msg = []
        for dt,OxxNodes in sorted(self.oxx.items()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            #msg += header+words
            quadMsg = []
            triMsg = []
            for eid,oxxNodes in sorted(OxxNodes.items()):
                eType = self.eType[eid]
                if eType in 'CQUAD4':
                    if isBilinear:
                        out = self.writeF06_Quad4_BilinearTransient(dt,eid,oxxNodes)
                    else:
                        out = self.writeF06_Tri3Transient(dt,eid,oxxNodes)
                    quadMsg.append(out)
                elif eType in 'CTRIA3':
                    out = self.writeF06_Tri3Transient(dt,eid,oxxNodes)
                    triMsg.append(out)
                #else:
                #    raise NotImplementedError('eType = |%r|' %(eType)) # CQUAD8, CTRIA6
                ###
            ###
            if isQuad:
                quadMsg.append(pageStamp+str(pageNum)+'\n')
                msg += header+quadMsg2+quadMsg
                pageNum+=1
            if isTri:
                triMsg.append(pageStamp+str(pageNum)+'\n')
                msg += header+triMsg2+triMsg
                pageNum+=1
            ###
        ###
        return (''.join(msg),pageNum-1)

    def writeF06_Quad4_Bilinear(self,eid,oxxNodes):
        msg = ''
        k = oxxNodes.keys()
        k.sort()
        k.pop(-1)
        for nid in ['C']+k:
            for iLayer in range(len(self.oxx[eid][nid])):
                fd    = self.fiberCurvature[eid][nid][iLayer]
                oxx   = self.oxx[eid][nid][iLayer]
                oyy   = self.oyy[eid][nid][iLayer]
                txy   = self.txy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                ovm = self.ovmShear[eid][nid][iLayer]
                ([fd,oxx,oyy,txy,major,minor,ovm],isAllZeros) = self.writeF06Floats13E([fd,oxx,oyy,txy,major,minor,ovm])
                ([angle],isAllZeros) = self.writeF06Floats8p4F([angle])

                if nid=='C' and iLayer==0:
                    msg += '0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' %(eid,'CEN/4',fd,oxx,oyy,txy,angle,major,minor,ovm)
                elif iLayer==0:
                    msg += '   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' %('',nid,     fd,oxx,oyy,txy,angle,major,minor,ovm)
                elif iLayer==1:
                    msg += '   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n\n' %('','',    fd,oxx,oyy,txy,angle,major,minor,ovm)
                else:
                    raise Exception('Invalid option for cquad4')
                ###
            ###
        ###
        return msg

    def writeF06_Quad4_BilinearTransient(self,dt,eid,oxxNodes):
        msg = ''
        k = oxxNodes.keys()
        k.sort()
        k.pop(-1)
        for nid in ['C']+k:
            for iLayer in range(len(self.oxx[dt][eid][nid])):
                fd    = self.fiberCurvature[dt][eid][nid][iLayer]
                oxx   = self.oxx[dt][eid][nid][iLayer]
                oyy   = self.oyy[dt][eid][nid][iLayer]
                txy   = self.txy[dt][eid][nid][iLayer]
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                ovm = self.ovmShear[dt][eid][nid][iLayer]
                ([fd,oxx,oyy,txy,major,minor,ovm],isAllZeros) = self.writeF06Floats13E([fd,oxx,oyy,txy,major,minor,ovm])
                ([angle],isAllZeros) = self.writeF06Floats8p4F([angle])

                if nid=='C' and iLayer==0:
                    msg += '0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' %(eid,'CEN/4',fd,oxx,oyy,txy,angle,major,minor,ovm.rstrip())
                elif iLayer==0:
                    msg += '   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' %('',nid,     fd,oxx,oyy,txy,angle,major,minor,ovm.rstrip())
                elif iLayer==1:
                    msg += '   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n\n' %('','',    fd,oxx,oyy,txy,angle,major,minor,ovm.rstrip())
                else:
                    #msg += '   %8s %8s  %13E  %13E %13E %13E   %8.4F  %13E %13E %13E\n' %('','',  fd,oxx,oyy,txy,angle,major,minor,ovm)
                    raise Exception('Invalid option for cquad4')
                ###
            ###
        ###
        return msg

    def writeF06_Tri3(self,eid,oxxNodes):
        msg = ''
        for nid in sorted(oxxNodes):
            for iLayer in range(len(self.oxx[eid][nid])):
                fd    = self.fiberCurvature[eid][nid][iLayer]
                oxx   = self.oxx[eid][nid][iLayer]
                oyy   = self.oyy[eid][nid][iLayer]
                txy   = self.txy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                ovm = self.ovmShear[eid][nid][iLayer]
                ([fd,oxx,oyy,txy,major,minor,ovm],isAllZeros) = self.writeF06Floats13E([fd,oxx,oyy,txy,major,minor,ovm])
                ([angle],isAllZeros) = self.writeF06Floats8p4F([angle])

                if iLayer==0:
                    msg += '0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' %(eid,fd,oxx,oyy,txy,angle,major,minor,ovm.rstrip())
                else:
                    msg += '   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' %('', fd,oxx,oyy,txy,angle,major,minor,ovm.rstrip())
                ###
            ###
        ###
        return msg

    def writeF06_Tri3Transient(self,dt,eid,oxxNodes):
        msg = ''
        for nid in sorted(oxxNodes):
            for iLayer in range(len(self.oxx[dt][eid][nid])):
                fd    = self.fiberCurvature[dt][eid][nid][iLayer]
                oxx   = self.oxx[dt][eid][nid][iLayer]
                oyy   = self.oyy[dt][eid][nid][iLayer]
                txy   = self.txy[dt][eid][nid][iLayer]
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                ovm = self.ovmShear[dt][eid][nid][iLayer]
                ([fd,oxx,oyy,txy,major,minor,ovm],isAllZeros) = self.writeF06Floats13E([fd,oxx,oyy,txy,major,minor,ovm])
                ([angle],isAllZeros) = self.writeF06Floats8p4F([angle])

                if iLayer==0:
                    msg += '0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' %(eid,fd,oxx,oyy,txy,angle,major,minor,ovm)
                else:
                    msg += '   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' %('',  fd,oxx,oyy,txy,angle,major,minor,ovm)
                ###
            ###
        ###
        return msg

    def __repr__(self):
        #print "sCodes = ",self.sCodes
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---ISOTROPIC PLATE STRESS---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' %('EID','eType','nodeID','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        #print self.oxx.keys()
        for eid,oxxNodes in sorted(self.oxx.items()):
            eType = self.eType[eid]
            for nid in sorted(oxxNodes):
                for iLayer in range(len(self.oxx[eid][nid])):
                    fd    = self.fiberCurvature[eid][nid][iLayer]
                    oxx   = self.oxx[eid][nid][iLayer]
                    oyy   = self.oyy[eid][nid][iLayer]
                    txy   = self.txy[eid][nid][iLayer]
                    angle = self.angle[eid][nid][iLayer]
                    major = self.majorP[eid][nid][iLayer]
                    minor = self.minorP[eid][nid][iLayer]
                    ovm = self.ovmShear[eid][nid][iLayer]

                    msg += '%-6i %6s %8s %7s %10g ' %(eid,eType,nid,iLayer+1,fd)
                    vals = [oxx,oyy,txy,major,minor,ovm]
                    for val in vals:
                        if abs(val)<1e-6:
                            msg += '%10s ' %('0')
                        else:
                            try:
                                msg += '%10i ' %(val)
                            except:
                                print("bad val = ",val)
                                raise
                        ###
                    msg += '\n'
                ###
            ###
        ###
        return msg

class plateStrainObject(strainObject):
    """
    # ??? - is this just 11
    ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 
      ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
    
    # sCode=11
                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  
    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               
      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       VON MISES 

    # sCode=15
                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )
    ELEMENT      FIBER                STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 
      ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES

    # sCode=10
                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  
    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)          MAX  
      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR         SHEAR   
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)
        self.eType     = {}

        self.code = [self.formatCode,self.sortCode,self.sCode]
        
        #print self.dataCode
        self.fiberCurvature = {}
        self.exx    = {}
        self.eyy    = {}
        self.exy    = {}
        self.angle  = {}
        self.majorP = {}
        self.minorP = {}
        self.evmShear = {}

        if dt is not None:
            self.dt = dt
            self.addNewTransient()
            self.add       = self.addTransient
            self.addNewEid = self.addNewEidTransient
            self.isTransient = True
        ###

    def addF06Data(self,data,transient):
        if transient is None:
            eType = data[0][0]
            for line in data:
                if eType=='CTRIA3':
                    (eType,eid,f1,ex1,ey1,exy1,angle1,e11,e21,evm1,
                               f2,ex2,ey2,exy2,angle2,e12,e22,evm2) = line
                    self.eType[eid] = eType
                    self.fiberCurvature[eid] = {'C':[f1,f2]}
                    self.exx[eid]      = {'C':[ex1,ex2]}
                    self.eyy[eid]      = {'C':[ey1,ey2]}
                    self.exy[eid]      = {'C':[exy1,exy2]}
                    self.angle[eid]    = {'C':[angle1,angle2]}
                    self.majorP[eid]   = {'C':[e11,e12]}
                    self.minorP[eid]   = {'C':[e21,e22]}
                    self.evmShear[eid] = {'C':[evm1,evm2]}
                elif eType=='CQUAD4':
                    #assert len(line)==19,'len(line)=%s' %(len(line))
                    #print line
                    if len(line)==19: # Centroid - bilinear
                        (eType,eid,nid,f1,ex1,ey1,exy1,angle1,e11,e21,evm1,
                                       f2,ex2,ey2,exy2,angle2,e12,e22,evm2) = line
                        if nid=='CEN/4': nid='C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid:[f1,f2]}
                        self.exx[eid]      = {nid:[ex1,ex2]}
                        self.eyy[eid]      = {nid:[ey1,ey2]}
                        self.exy[eid]      = {nid:[exy1,exy2]}
                        self.angle[eid]    = {nid:[angle1,angle2]}
                        self.majorP[eid]   = {nid:[e11,e12]}
                        self.minorP[eid]   = {nid:[e21,e22]}
                        self.evmShear[eid] = {nid:[evm1,evm2]}
                    elif len(line)==18: # Centroid
                        (eType,eid,f1,ex1,ey1,exy1,angle1,e11,e21,evm1,
                                   f2,ex2,ey2,exy2,angle2,e12,e22,evm2) = line
                        nid='C'
                        self.eType[eid] = eType
                        self.fiberCurvature[eid] = {nid:[f1,f2]}
                        self.exx[eid]      = {nid:[ex1,ex2]}
                        self.eyy[eid]      = {nid:[ey1,ey2]}
                        self.exy[eid]      = {nid:[exy1,exy2]}
                        self.angle[eid]    = {nid:[angle1,angle2]}
                        self.majorP[eid]   = {nid:[e11,e12]}
                        self.minorP[eid]   = {nid:[e21,e22]}
                        self.evmShear[eid] = {nid:[evm1,evm2]}
                    elif len(line)==17: ## Bilinear node
                        #print line
                        (nid,f1,ex1,ey1,exy1,angle1,e11,e21,evm1,
                             f2,ex2,ey2,exy2,angle2,e12,e22,evm2) = line
                        self.fiberCurvature[eid][nid] = [f1,f2]
                        self.exx[eid][nid]      = [ex1,ex2]
                        self.eyy[eid][nid]      = [ey1,ey2]
                        self.txy[eid][nid]      = [exy1,exy2]
                        self.angle[eid][nid]    = [angle1,angle2]
                        self.majorP[eid][nid]   = [e11,e12]
                        self.minorP[eid][nid]   = [e21,e22]
                        self.evmShear[eid][nid] = [evm1,evm2]
                    else:
                        assert len(line)==19,'len(line)=%s' %(len(line))
                        #raise NotImplementedError()
                    ###
                else:
                    raise NotImplementedError('line=%s not supported...' %(line))
            return
        raise NotImplementedError('transient results not supported')

    def deleteTransient(self,dt):
        del self.fiberCurvature[dt]
        del self.exx[dt]
        del self.eyy[dt]
        del self.exy[dt]
        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]
        del self.evmShear[dt]

    def getTransients(self):
        k = self.exx.keys()
        k.sort()
        return k

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        if self.dt not in self.exx:
            self.fiberCurvature[self.dt] = {}
            self.exx[self.dt]   = {}
            self.eyy[self.dt]   = {}
            self.exy[self.dt]   = {}
            self.angle[self.dt] = {}
            self.majorP[self.dt] = {}
            self.minorP[self.dt] = {}
            self.evmShear[self.dt]   = {}

    def addNewEid(self,eType,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        
        if nodeID is not 'C': # centroid
            assert 0<nodeID<1000000000, 'nodeID=%s %s' %(nodeID,msg)

        if eid in self.exx:
            return self.add(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        assert eid not in self.exx
        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [curvature]}
        self.exx[eid]    = {nodeID: [exx]}
        self.eyy[eid]    = {nodeID: [eyy]}
        self.exy[eid]    = {nodeID: [exy]}
        self.angle[eid]  = {nodeID: [angle]}
        self.majorP[eid] = {nodeID: [majorP]}
        self.minorP[eid] = {nodeID: [minorP]}
        self.evmShear[eid]    = {nodeID: [evm]}
        #print msg
        if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        #print msg
        dt = self.dt
        if nodeID is not 'C': # centroid
            assert 0<nodeID<1000000000, 'nodeID=%s %s' %(nodeID,msg)

        if eid in self.evmShear[dt]:  # SOL200, erase the old result
            nid = nodeID
            #msg = "dt=%s eid=%s nodeID=%s fd=%s oxx=%s major=%s vm=%s" %(dt,eid,nodeID,str(self.fiberCurvature[dt][eid][nid]),str(self.oxx[dt][eid][nid]),str(self.majorP[dt][eid][nid]),str(self.ovmShear[dt][eid][nid]))
            self.deleteTransient(dt)
            self.addNewTransient()
        
        self.eType[eid] = eType
        self.fiberCurvature[dt][eid] = {nodeID: [curvature]}
        self.exx[dt][eid]    = {nodeID: [exx]}
        self.eyy[dt][eid]    = {nodeID: [eyy]}
        self.exy[dt][eid]    = {nodeID: [exy]}
        self.angle[dt][eid]  = {nodeID: [angle]}
        self.majorP[dt][eid] = {nodeID: [majorP]}
        self.minorP[dt][eid] = {nodeID: [minorP]}
        self.evmShear[dt][eid]    = {nodeID: [evm]}
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        #print "***add"
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        if nodeID is not 'C': # centroid
            assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        self.fiberCurvature[eid][nodeID].append(curvature)
        self.exx[eid][nodeID].append(exx)
        self.eyy[eid][nodeID].append(eyy)
        self.exy[eid][nodeID].append(exy)
        self.angle[eid][nodeID].append(angle)
        self.majorP[eid][nodeID].append(majorP)
        self.minorP[eid][nodeID].append(minorP)
        self.evmShear[eid][nodeID].append(evm)
        if nodeID==0: raise Exception(msg)

    def addTransient(self,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        #print "***add"
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        #print msg
        #print self.oxx
        #print self.fiberCurvature
        if nodeID is not 'C': # centroid
            assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        dt = self.dt
        self.fiberCurvature[dt][eid][nodeID].append(curvature)
        self.exx[dt][eid][nodeID].append(exx)
        self.eyy[dt][eid][nodeID].append(eyy)
        self.exy[dt][eid][nodeID].append(exy)
        self.angle[dt][eid][nodeID].append(angle)
        self.majorP[dt][eid][nodeID].append(majorP)
        self.minorP[dt][eid][nodeID].append(minorP)
        self.evmShear[dt][eid][nodeID].append(evm)
        if nodeID==0: raise Exception(msg)

    def addNewNode(self,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        #print "***addNewNode"
        #print self.oxx
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        assert nodeID not in self.exx[eid],msg
        self.fiberCurvature[eid][nodeID] = [curvature]
        self.exx[eid][nodeID]    = [exx]
        self.eyy[eid][nodeID]    = [eyy]
        self.exy[eid][nodeID]    = [exy]
        self.angle[eid][nodeID]  = [angle]
        self.majorP[eid][nodeID] = [majorP]
        self.minorP[eid][nodeID] = [minorP]
        self.evmShear[eid][nodeID]    = [evm]
        #print msg
        if nodeID==0: raise Exception(msg)

    def getHeaders(self):
        if self.isFiberDistance():
            headers = ['fiberDist']
        else:
            headers = ['curvature']
    
        headers += ['exx','eyy','exy','eMajor','eMinor']
        if self.isVonMises():
            headers.append('eVonMises')
        else:
            headers.append('maxShear')
        return headers

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance(): 
            quadMsgTemp = ['    ELEMENT              FIBER                STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' %(vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                           '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            triMsgTemp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                          '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' %(vonMises)]
        ###

        triMsg = []
        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            qkey = eTypes.index('CQUAD4')
            kkey = self.eType.keys()[qkey]
            ekey = self.exx[kkey].keys()
            isBilinear=True
            quadMsg     = header+['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n']+quadMsgTemp
            if len(ekey)==1:
                isBilinear=False
                quadMsg = header+['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']+triMsgTemp
            isQuad = True
        else:
            quadMsg = []
            isQuad = False

        if 'CTRIA3' in eTypes:
            isTri = True
            triMsg = header+['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']+triMsgTemp
        else:
            isTri = False
            triMsg = []

        for eid,exxNodes in sorted(self.exx.items()):
            eType = self.eType[eid]
            if eType in 'CQUAD4':
                if isBilinear:
                    out = self.writeF06_Quad4_Bilinear(eid,exxNodes)
                else:
                    out = self.writeF06_Tri3(eid,exxNodes)
                quadMsg.append(out)
            elif eType in 'CTRIA3':
                out = self.writeF06_Tri3(eid,exxNodes)
                triMsg.append(out)
            #else:
            #    raise NotImplementedError('eType = |%r|' %(eType)) # CQUAD8, CTRIA6
            ###
        ###
        if isQuad:
            quadMsg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        if isTri:
            triMsg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1

        #print "quadMsg = ",quadMsg
        #print "triMsg = ",triMsg
        msg = ''.join(quadMsg+triMsg)
        return (msg,pageNum-1)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        if self.isFiberDistance(): 
            quadMsgTemp = ['    ELEMENT              FIBER                STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                           '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            triMsgTemp = ['  ELEMENT      FIBER               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                          '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' %(vonMises)]
        else:
            quadMsgTemp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                           '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' %(vonMises)]
            triMsgTemp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)\n',
                          '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' %(vonMises)]
        ###

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes:
            ElemKey = eTypes.index('CQUAD4')
            #print qkey
            eid = self.eType.keys()[ElemKey]
            #print "self.oxx = ",self.oxx
            #print "eid=%s" %(eid)
            dt = self.exx.keys()[0]
            #print "dt=%s" %(dt)
            nLayers = len(self.exx[dt][eid])
            #print "elementKeys = ",elementKeys
            isBilinear=True
            quadMsg2     = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n']+quadMsgTemp
            if nLayers==1:
                isBilinear=False
                quadMsg2 = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']+triMsgTemp
            isQuad = True
        else:
            quadMsg2 = []
            isQuad = False

        if 'CTRIA3' in eTypes:
            isTri = True
            triMsg2 = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']+triMsgTemp
        else:
            isTri = False
            triMsg2 = []

        msg = []
        for dt,ExxNodes in sorted(self.exx.items()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            quadMsg = []
            triMsg = []
            if isQuad: quadMsg = header+quadMsg2
            if isTri:  triMsg  = header+triMsg2
            for eid,exxNodes in sorted(ExxNodes.items()):
                eType = self.eType[eid]
                if eType in 'CQUAD4':
                    if isBilinear:
                        out = self.writeF06_Quad4_BilinearTransient(dt,eid,exxNodes)
                    else:
                        out = self.writeF06_Tri3Transient(dt,eid,exxNodes)
                    quadMsg.append(out)
                elif eType in 'CTRIA3':
                    out = self.writeF06_Tri3Transient(dt,eid,exxNodes)
                    triMsg.append(out)
                #else:
                #    raise NotImplementedError('eType = |%r|' %(eType)) # CQUAD8, CTRIA6
                ###
            ###
            if isQuad:
                quadMsg.append(pageStamp+str(pageNum)+'\n')
                msg += header+quadMsg2+quadMsg
                pageNum+=1
            if isTri:
                triMsg.append(pageStamp+str(pageNum)+'\n')
                msg += header+triMsg2+triMsg
                pageNum+=1
            ###
        ###
        return (''.join(msg),pageNum-1)

    def writeF06_Quad4_Bilinear(self,eid,exxNodes):
        msg = ''
        k = exxNodes.keys()
        k.sort()
        k.pop(-1)
        for nid in ['C']+k:
            for iLayer in range(len(self.exx[eid][nid])):
                fd    = self.fiberCurvature[eid][nid][iLayer]
                exx   = self.exx[eid][nid][iLayer]
                eyy   = self.eyy[eid][nid][iLayer]
                exy   = self.exy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                evm   = self.evmShear[eid][nid][iLayer]
                ([fd,exx,eyy,exy,major,minor,evm],isAllZeros) = self.writeF06Floats13E([fd,exx,eyy,exy,major,minor,evm])
                ([angle],isAllZeros) = self.writeF06Floats8p4F([angle])

                if nid=='C' and iLayer==0:
                    msg += '0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' %(eid,'CEN/4',fd,exx,eyy,exy,angle,major,minor,evm.rstrip())
                elif iLayer==0:
                    msg += '   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' %('',nid,     fd,exx,eyy,exy,angle,major,minor,evm.rstrip())
                elif iLayer==1:
                    msg += '   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n\n' %('','',    fd,exx,eyy,exy,angle,major,minor,evm.rstrip())
                else:
                    raise Exception('Invalid option for cquad4')
                ###
            ###
        ###
        return msg

    def writeF06_Quad4_BilinearTransient(self,dt,eid,exxNodes):
        msg = ''
        k = exxNodes.keys()
        k.sort()
        k.pop(-1)
        for nid in ['C']+k:
            for iLayer in range(len(self.exx[dt][eid][nid])):
                fd    = self.fiberCurvature[dt][eid][nid][iLayer]
                exx   = self.exx[dt][eid][nid][iLayer]
                eyy   = self.eyy[dt][eid][nid][iLayer]
                exy   = self.exy[dt][eid][nid][iLayer]
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                evm   = self.evmShear[dt][eid][nid][iLayer]

                ([fd,exx,eyy,exy,major,minor,evm],isAllZeros) = self.writeF06Floats13E([fd,exx,eyy,exy,major,minor,evm])
                ([angle],isAllZeros) = self.writeF06Floats8p4F([angle])

                if nid=='C' and iLayer==0:
                    msg += '0  %8i %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' %(eid,'CEN/4',fd,exx,eyy,exy,angle,major,minor,evm.rstrip())
                elif iLayer==0:
                    msg += '   %8s %8i  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n' %('',nid,     fd,exx,eyy,exy,angle,major,minor,evm.rstrip())
                elif iLayer==1:
                    msg += '   %8s %8s  %13s  %13s %13s %13s   %8s  %13s %13s %-s\n\n' %('','',    fd,exx,eyy,exy,angle,major,minor,evm.rstrip())
                else:
                    raise Exception('Invalid option for cquad4')
                ###
            ###
        ###
        return msg

    def writeF06_Tri3(self,eid,exxNodes):
        msg = ''
        for nid in sorted(exxNodes):
            for iLayer in range(len(self.exx[eid][nid])):
                fd    = self.fiberCurvature[eid][nid][iLayer]
                exx   = self.exx[eid][nid][iLayer]
                eyy   = self.eyy[eid][nid][iLayer]
                exy   = self.exy[eid][nid][iLayer]
                angle = self.angle[eid][nid][iLayer]
                major = self.majorP[eid][nid][iLayer]
                minor = self.minorP[eid][nid][iLayer]
                evm   = self.evmShear[eid][nid][iLayer]

                ([fd,exx,eyy,exy,major,minor,evm],isAllZeros) = self.writeF06Floats13E([fd,exx,eyy,exy,major,minor,evm])
                ([angle],isAllZeros) = self.writeF06Floats8p4F([angle])
                if iLayer==0:
                    msg += '0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' %(eid,fd,exx,eyy,exy,angle,major,minor,evm.rstrip())
                else:
                    msg += '   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' %('', fd,exx,eyy,exy,angle,major,minor,evm.rstrip())
                ###
            ###
        ###
        return msg

    def writeF06_Tri3Transient(self,dt,eid,exxNodes):
        msg = ''
        for nid in sorted(exxNodes):
            for iLayer in range(len(self.exx[dt][eid][nid])):
                fd    = self.fiberCurvature[dt][eid][nid][iLayer]
                exx   = self.exx[dt][eid][nid][iLayer]
                eyy   = self.eyy[dt][eid][nid][iLayer]
                exy   = self.exy[dt][eid][nid][iLayer]
                angle = self.angle[dt][eid][nid][iLayer]
                major = self.majorP[dt][eid][nid][iLayer]
                minor = self.minorP[dt][eid][nid][iLayer]
                evm   = self.evmShear[dt][eid][nid][iLayer]

                ([fd,exx,eyy,exy,major,minor,evm],isAllZeros) = self.writeF06Floats13E([fd,exx,eyy,exy,major,minor,evm])
                ([angle],isAllZeros) = self.writeF06Floats8p4F([angle])
                if iLayer==0:
                    msg += '0  %6i   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' %(eid,fd,exx,eyy,exy,angle,major,minor,evm)
                else:
                    msg += '   %6s   %13s     %13s  %13s  %13s   %8s   %13s   %13s  %-s\n' %('', fd,exx,eyy,exy,angle,major,minor,evm)
                ###
            ###
        ###
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---ISOTROPIC PLATE STRAIN---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' %('EID','eType','nodeID','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for eid,exxNodes in sorted(self.exx.items()):
            eType = self.eType[eid]
            for nid in sorted(exxNodes):
                for iLayer in range(len(self.exx[eid][nid])):
                    fd    = self.fiberCurvature[eid][nid][iLayer]
                    exx   = self.exx[eid][nid][iLayer]
                    eyy   = self.eyy[eid][nid][iLayer]
                    exy   = self.exy[eid][nid][iLayer]
                    angle = self.angle[eid][nid][iLayer]
                    major = self.majorP[eid][nid][iLayer]
                    minor = self.minorP[eid][nid][iLayer]
                    evm   = self.evmShear[eid][nid][iLayer]
                    
                    msg += '%-6i %6s %8s %7s %10g ' %(eid,eType,nid,iLayer+1,fd)
                    vals = [exx,eyy,exy,major,minor,evm]
                    for val in vals:
                        if abs(val)<1e-6:
                            msg += '%10s ' %('0.')
                        else:
                            msg += '%10.3g ' %(val)
                        ###
                    msg += '\n'

                    #msg += "eid=%s eType=%s nid=%s iLayer=%s exx=%-9.3g eyy=%-9.3g exy=%-9.3g evm=%-9.3g\n" %(eid,eType,nid,iLayer,exx,eyy,exy,evm)
                ###
            ###
        ###
        return msg

    def __reprTransient__(self):
        msg = '---ISOTROPIC PLATE STRAIN---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s %7s ' %('EID','eType','nodeID','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,exx in sorted(self.exx.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid,exxNodes in sorted(exx.items()):
                eType = self.eType[eid]
                for nid in sorted(exxNodes):
                    for iLayer in range(len(self.exx[dt][eid][nid])):
                        fd    = self.fiberCurvature[dt][eid][nid][iLayer]
                        exx   = self.exx[   dt][eid][nid][iLayer]
                        eyy   = self.eyy[   dt][eid][nid][iLayer]
                        exy   = self.exy[   dt][eid][nid][iLayer]
                        angle = self.angle[ dt][eid][nid][iLayer]
                        major = self.majorP[dt][eid][nid][iLayer]
                        minor = self.minorP[dt][eid][nid][iLayer]
                        evm   = self.evmShear[   dt][eid][nid][iLayer]

                        msg += '%-6i %6s %8s %7s %10g ' %(eid,eType,nid,iLayer+1,fd)
                        vals = [exx,eyy,exy,major,minor,evm]
                        for val in vals:
                            if abs(val)<1e-6:
                                msg += '%10s ' %('0.')
                            else:
                                msg += '%10.3g ' %(val)
                            ###
                        msg += '\n'

                        #msg += "eid=%s eType=%s nid=%s iLayer=%s exx=%-9.3g eyy=%-9.3g exy=%-9.3g evm=%-9.3g\n" %(eid,eType,nid,iLayer,exx,eyy,exy,evm)
                    ###
                ###
            ###
        ###
        return msg
