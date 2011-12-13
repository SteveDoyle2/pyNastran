import sys
from oes_objects import stressObject,strainObject
from pyNastran.op2.op2Errors import *

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
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)
        self.eType     = {}

        self.appendDataMember('sCodes','sCode')
        self.code = [self.formatCode,self.sortCode,self.sCode]
        
        print self.dataCode
        if self.code == [1,0,11]:
            self.fiberCurvature = {}
            self.exx    = {}
            self.eyy    = {}
            self.exy    = {}
            self.angle  = {}
            self.majorP = {}
            self.minorP = {}
            self.evm    = {}
            self.isFiberDistance = False
            self.isBilinear = True
        
        elif self.code == [1,0,15]:
            self.fiberCurvature = {}
            self.exx    = {}
            self.eyy    = {}
            self.exy    = {}
            self.angle  = {}
            self.majorP = {}
            self.minorP = {}
            self.evm    = {}
            self.isFiberDistance = True
            self.isBilinear = False
        else:
            raise InvalidCodeError('plateStrain - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt:
            self.dt = dt
            self.addNewTransient()
            self.add       = self.addTransient
            self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.fiberCurvature[self.dt] = {}
        self.exx[self.dt]   = {}
        self.eyy[self.dt]   = {}
        self.exy[self.dt]   = {}
        self.angle[self.dt] = {}
        self.majorP[self.dt] = {}
        self.minorP[self.dt] = {}
        self.evm[self.dt]   = {}

    def addNewEid(self,eType,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        
        if nodeID is not 'C': # centroid
            assert 0<nodeID<1000000000, 'nodeID=%s %s' %(nodeID,msg)
        assert eid not in self.exx
        self.eType[eid] = eType
        self.fiberCurvature[eid] = {nodeID: [curvature]}
        self.exx[eid]    = {nodeID: [exx]}
        self.eyy[eid]    = {nodeID: [eyy]}
        self.exy[eid]    = {nodeID: [exy]}
        self.angle[eid]  = {nodeID: [angle]}
        self.majorP[eid] = {nodeID: [majorP]}
        self.minorP[eid] = {nodeID: [minorP]}
        self.evm[eid]    = {nodeID: [evm]}
        #print msg
        if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        dt = self.dt
        if nodeID is not 'C': # centroid
            assert 0<nodeID<1000000000, 'nodeID=%s %s' %(nodeID,msg)
        assert eid not in self.exx[dt]
        self.eType[eid] = eType
        self.fiberCurvature[dt][eid] = {nodeID: [curvature]}
        self.exx[dt][eid]    = {nodeID: [exx]}
        self.eyy[dt][eid]    = {nodeID: [eyy]}
        self.exy[dt][eid]    = {nodeID: [exy]}
        self.angle[dt][eid]  = {nodeID: [angle]}
        self.majorP[dt][eid] = {nodeID: [majorP]}
        self.minorP[dt][eid] = {nodeID: [minorP]}
        self.evm[dt][eid]    = {nodeID: [evm]}
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        #print "***add"
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        #print msg
        #print self.oxx
        #print self.fiberDistance
        if nodeID is not 'C': # centroid
            assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        self.fiberCurvature[eid][nodeID].append(curvature)
        self.exx[eid][nodeID].append(exx)
        self.eyy[eid][nodeID].append(eyy)
        self.exy[eid][nodeID].append(exy)
        self.angle[eid][nodeID].append(angle)
        self.majorP[eid][nodeID].append(majorP)
        self.minorP[eid][nodeID].append(minorP)
        self.evm[eid][nodeID].append(evm)
        if nodeID==0: raise Exception(msg)

    def addTransient(self,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        #print "***add"
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        #print msg
        #print self.oxx
        #print self.fiberDistance
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
        self.evm[dt][eid][nodeID].append(evm)
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
        self.evm[eid][nodeID]    = [evm]
        #print msg
        if nodeID==0: raise Exception(msg)

    def __reprTransient__(self):
        msg = '---ISOTROPIC PLATE STRAIN---\n'
        headers = ['exx','eyy','exy','eMajor','eMinor','evm']
        msg += '%-6s %6s %8s %7s ' %('EID','eType','nodeID','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,exx in sorted(self.exx.items()):
            msg += 'dt = %s\n' %(dt)
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
                        evm   = self.evm[   dt][eid][nid][iLayer]

                        msg += '%-6i %6s %8s %7s ' %(eid,eType,nid,iLayer+1)
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

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---ISOTROPIC PLATE STRAIN---\n'
        headers = ['exx','eyy','exy','eMajor','eMinor','evm']
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
                    evm   = self.evm[eid][nid][iLayer]
                    
                    msg += '%-6i %6s %8s %7s ' %(eid,eType,nid,iLayer+1)
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

class plateStressObject(stressObject):
    """
    ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 
      ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = {}

        #self.appendDataMember('sCodes','sCode')
        self.code = [self.formatCode,self.sortCode,self.sCode]
        if self.code == [1,0,1]:
            self.fiberDistance = {}
            self.oxx    = {}
            self.oyy    = {}
            self.txy    = {}
            self.angle  = {}
            self.majorP = {}
            self.minorP = {}
            self.ovm    = {}
        else:
            raise InvalidCodeError('plateStress - get the format/sort/stressCode=%s' %(self.code))
        ###


        #print "%%%%%dt = ",dt
        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.add       = self.addTransient
            self.addNewEid = self.addNewEidTransient
        else:
            self.dt = None
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.fiberDistance[self.dt] = {}
        self.oxx[self.dt]    = {}
        self.oyy[self.dt]    = {}
        self.txy[self.dt]    = {}
        self.angle[self.dt]  = {}
        self.majorP[self.dt] = {}
        self.minorP[self.dt] = {}
        self.ovm[self.dt]    = {}

    def addNewEid(self,eType,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "Plate Stress add..."
        #assert eid not in self.oxx
        #print self.oxx

        self.eType[eid] = eType
        self.fiberDistance[eid] = {nodeID: [fd]}
        assert eid is not None
        self.oxx[eid]    = {nodeID: [oxx]}
        self.oyy[eid]    = {nodeID: [oyy]}
        self.txy[eid]    = {nodeID: [txy]}
        self.angle[eid]  = {nodeID: [angle]}
        self.majorP[eid] = {nodeID: [majorP]}
        self.minorP[eid] = {nodeID: [minorP]}
        self.ovm[eid]    = {nodeID: [ovm]}
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "Plate Stress Transient add..."
        dt = self.dt
        #msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(dt,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g major=%g vm=%g" %(dt,eid,nodeID,fd,oxx,majorP,ovm)
        #print msg
        assert eid not in self.oxx[self.dt],msg
        assert eid is not None
        self.eType[eid] = eType
        self.fiberDistance[dt][eid] = {nodeID: [fd]}
        self.oxx[dt][eid]    = {nodeID: [oxx]}
        self.oyy[dt][eid]    = {nodeID: [oyy]}
        self.txy[dt][eid]    = {nodeID: [txy]}
        self.angle[dt][eid]  = {nodeID: [angle]}
        self.majorP[dt][eid] = {nodeID: [majorP]}
        self.minorP[dt][eid] = {nodeID: [minorP]}
        self.ovm[dt][eid]    = {nodeID: [ovm]}
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***add"
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #print self.oxx
        #print self.fiberDistance
        assert eid is not None
        self.fiberDistance[eid][nodeID].append(fd)
        self.oxx[eid][nodeID].append(oxx)
        self.oyy[eid][nodeID].append(oyy)
        self.txy[eid][nodeID].append(txy)
        self.angle[eid][nodeID].append(angle)
        self.majorP[eid][nodeID].append(majorP)
        self.minorP[eid][nodeID].append(minorP)
        self.ovm[eid][nodeID].append(ovm)
        if nodeID==0: raise Exception(msg)

    def addTransient(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***add"
        dt = self.dt
        msg = "dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(dt,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #print self.oxx
        #print self.fiberDistance
        assert eid is not None
        self.fiberDistance[dt][eid][nodeID].append(fd)
        self.oxx[dt][eid][nodeID].append(oxx)
        self.oyy[dt][eid][nodeID].append(oyy)
        self.txy[dt][eid][nodeID].append(txy)
        self.angle[dt][eid][nodeID].append(angle)
        self.majorP[dt][eid][nodeID].append(majorP)
        self.minorP[dt][eid][nodeID].append(minorP)
        self.ovm[dt][eid][nodeID].append(ovm)
        if nodeID==0: raise Exception(msg)

    def addNewNode(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***addNewNode"
        #print self.oxx
        assert eid is not None
        #assert nodeID not in self.oxx[eid]
        self.fiberDistance[eid][nodeID] = [fd]
        self.oxx[eid][nodeID]    = [oxx]
        self.oyy[eid][nodeID]    = [oyy]
        self.txy[eid][nodeID]    = [txy]
        self.angle[eid][nodeID]  = [angle]
        self.majorP[eid][nodeID] = [majorP]
        self.minorP[eid][nodeID] = [minorP]
        self.ovm[eid][nodeID]    = [ovm]
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def __reprTransient__(self):
        msg = '---ISOTROPIC PLATE STRESS---\n'
        headers = ['oxx','oyy','txy','ovm']
        msg += '%-6s %6s %8s %7s ' %('EID','eType','nodeID','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        #print self.oxx.keys()
        for dt,oxxs in sorted(self.oxx.items()):
            msg += 'dt = %s\n' %(dt)
            for eid,oxxNodes in sorted(oxxs.items()):
                eType = self.eType[eid]
                for nid in sorted(oxxNodes):
                    for iLayer in range(len(self.oxx[dt][eid][nid])):
                        fd    = self.fiberDistance[dt][eid][nid][iLayer]
                        oxx   = self.oxx[   dt][eid][nid][iLayer]
                        oyy   = self.oyy[   dt][eid][nid][iLayer]
                        txy   = self.txy[   dt][eid][nid][iLayer]
                        angle = self.angle[ dt][eid][nid][iLayer]
                        major = self.majorP[dt][eid][nid][iLayer]
                        minor = self.minorP[dt][eid][nid][iLayer]
                        ovm   = self.ovm[   dt][eid][nid][iLayer]

                        msg += '%-6i %6s %8s %7s ' %(eid,eType,nid,iLayer+1)
                        vals = [oxx,oyy,txy,ovm]
                        for val in vals:
                            if abs(val)<1e-6:
                                msg += '%10s ' %('0')
                            else:
                                msg += '%10i ' %(val)
                            ###
                        msg += '\n'
                    ###
                ###
            ###
        ###
        return msg

    def __repr__(self):
        #print "sCodes = ",self.sCodes
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---ISOTROPIC PLATE STRESS---\n'
        headers = ['oxx','oyy','txy','ovm']
        msg += '%-6s %6s %8s %7s ' %('EID','eType','nodeID','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        #print self.oxx.keys()
        for eid,oxxNodes in sorted(self.oxx.items()):
            eType = self.eType[eid]
            for nid in sorted(oxxNodes):
                for iLayer in range(len(self.oxx[eid][nid])):
                    fd    = self.fiberDistance[eid][nid][iLayer]
                    oxx   = self.oxx[eid][nid][iLayer]
                    oyy   = self.oyy[eid][nid][iLayer]
                    txy   = self.txy[eid][nid][iLayer]
                    angle = self.angle[eid][nid][iLayer]
                    major = self.majorP[eid][nid][iLayer]
                    minor = self.minorP[eid][nid][iLayer]
                    ovm   = self.ovm[eid][nid][iLayer]

                    msg += '%-6i %6s %8s %7s ' %(eid,eType,nid,iLayer+1)
                    vals = [oxx,oyy,txy,ovm]
                    for val in vals:
                        if abs(val)<1e-6:
                            msg += '%10s ' %('0')
                        else:
                            msg += '%10i ' %(val)
                        ###
                    msg += '\n'
                ###
            ###
        ###
        return msg
