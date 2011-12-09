import sys
from op2_Objects import scalarObject #,array

class plateStrainObject(scalarObject):
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
        scalarObject.__init__(self,dataCode,iSubcase)
        self.eType     = {}

        self.appendDataMember('sCodes','sCode')
        print "self.sCode = ",self.sCode
        if self.sCode in [11]:
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
        
        elif self.sCode in [15]:
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
            raise Exception('get the sCode')
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
        self.ecc[self.dt]   = {}
        self.eyy[self.dt]   = {}
        self.exy[self.dt]   = {}
        self.angle[self.dt] = {}
        self.major[self.dt] = {}
        self.minor[self.dt] = {}
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

class plateStressObject(scalarObject):
    """
    ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 
      ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.eType = {}


        #self.appendDataMember('sCodes','sCode')
        #print "self.sCodes = ",self.sCodes
        print "self.sCode = ",self.sCode
        if self.sCode in [1]:
            self.fiberDistance = {}
            self.oxx    = {}
            self.oyy    = {}
            self.txy    = {}
            self.angle  = {}
            self.majorP = {}
            self.minorP = {}
            self.ovm    = {}
        else:
            raise Exception('get the sCode')
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
        assert eid not in self.oxx
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
        assert nodeID not in self.oxx[eid]
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
                        oxx   = self.oxx[dt][eid][nid][iLayer]
                        oyy   = self.oyy[dt][eid][nid][iLayer]
                        txy   = self.txy[dt][eid][nid][iLayer]
                        angle = self.angle[dt][eid][nid][iLayer]
                        major = self.majorP[dt][eid][nid][iLayer]
                        minor = self.minorP[dt][eid][nid][iLayer]
                        ovm   = self.ovm[dt][eid][nid][iLayer]

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


class compositePlateStressObject(scalarObject):
    """
    # sCode = 0
                    S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
    ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
      ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR

    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.eType  = {}

        print "self.sCode = ",self.sCode
        if self.sCode in [0]:
            self.fiberDistance = {}
            self.o11    = {}
            self.o22    = {}
            self.t12    = {}
            self.t1z    = {}
            self.t2z    = {}
            self.angle  = {}
            self.majorP = {}
            self.minorP = {}
            self.ovm    = {}
        else:
            raise Exception('get the sCode')
        ###

        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.add       = self.addTransient
            self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        #self.fiberDistance[self.dt] = {}
        self.o11[self.dt]    = {}
        self.o22[self.dt]    = {}
        self.t12[self.dt]    = {}
        self.t1z[self.dt]    = {}
        self.t2z[self.dt]    = {}
        self.angle[self.dt]  = {}
        self.majorP[self.dt] = {}
        self.minorP[self.dt] = {}
        self.ovm[self.dt]    = {}

    def addNewEid(self,eType,eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        #assert eid not in self.o11
        self.eType[eid]  = eType
        self.o11[eid]    = [o11]
        self.o22[eid]    = [o22]
        self.t12[eid]    = [t12]
        self.t1z[eid]    = [t1z]
        self.t2z[eid]    = [t2z]
        self.angle[eid]  = [angle]
        self.majorP[eid] = [majorP]
        self.minorP[eid] = [minorP]
        self.ovm[eid]    = [ovm]
        msg = "eid=%s o11=%g o22=%g t12=%g t1z=%g t2z=%g \nangle=%g major=%g minor=%g vm=%g" %(eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        dt = self.dt
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
        self.ovm[dt][eid]    = [ovm]
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
        self.ovm[eid].append(ovm)
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
        self.ovm[dt][eid].append(ovm)
        #if nodeID==0: raise Exception(msg)

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---COMPOSITE PLATE STRESS---\n'
        msg += '%-6s %8s %8s ' %('EID','eType','iLayer')
        headers = ['o11','o22','t12','t1z','t2z','ovm']
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
                ovm   = self.ovm[eid][iLayer]

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

class compositePlateStrainObject(scalarObject):
    """
    ???
    ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
      ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)

        self.eType  = {}
        print "self.sCode = ",self.sCode
        if self.sCode in [14]:
            self.e11    = {}
            self.e22    = {}
            self.e12    = {}
            self.e1z    = {}
            self.e2z    = {}
            self.angle  = {}
            self.majorP = {}
            self.minorP = {}
            self.evm    = {}
        else:
            raise Exception('get the sCode')
        ###

        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.add       = self.addTransient
            self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.fiberDistance[self.dt] = {}
        self.e11[self.dt]    = {}
        self.e22[self.dt]    = {}
        self.e12[self.dt]    = {}
        self.e1z[self.dt]    = {}
        self.e2z[self.dt]    = {}
        self.angle[self.dt]  = {}
        self.majorP[self.dt] = {}
        self.minorP[self.dt] = {}
        self.evm[self.dt]    = {}

    def addNewEid(self,eType,eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        assert eid not in self.e11
        self.eType[eid]  = eType
        self.e11[eid]    = [e11]
        self.e22[eid]    = [e22]
        self.e12[eid]    = [e12]
        self.e1z[eid]    = [e1z]
        self.e2z[eid]    = [e2z]
        self.angle[eid]  = [angle]
        self.majorP[eid] = [majorP]
        self.minorP[eid] = [minorP]
        self.evm[eid]    = [evm]
        msg = "eid=%s e11=%g e22=%g e12=%g e1z=%g e2z=%g \nangle=%g major=%g minor=%g vm=%g" %(eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        dt = self.dt
        assert eid not in self.e11[dt]
        self.eType[eid]  = eType
        self.e11[dt][eid]    = [e11]
        self.e22[dt][eid]    = [e22]
        self.e12[dt][eid]    = [e12]
        self.e1z[dt][eid]    = [e1z]
        self.e2z[dt][eid]    = [e2z]
        self.angle[dt][eid]  = [angle]
        self.majorP[dt][eid] = [majorP]
        self.minorP[dt][eid] = [minorP]
        self.evm[dt][eid]    = [evm]
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
        self.evm[eid].append(evm)
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
        self.evm[dt][eid].append(evm)
        #if nodeID==0: raise Exception(msg)

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---COMPOSITE PLATE STAIN---\n'
        headers = ['e11','e22','e12','e1z','e2z','evm']
        msg += '%-6s %8s %8s ' %('EID','eType','iLayer')
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
                evm   = self.evm[eid][iLayer]

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


class rodStressObject(scalarObject):
    """
    # sCode=0
                                  S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
    ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
      ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CROD'
        print "self.sCode = ",self.sCode
        if self.sCode in [0]:
            self.axial      = {}
            self.MS_axial   = {}
            self.torsion    = {}
            self.MS_torsion = {}
        else:
            raise Exception('get the sCode')
        ###
        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.axial[self.dt]      = {}
        self.MS_axial[self.dt]   = {}
        self.torsion[self.dt]    = {}
        self.MS_torsion[self.dt] = {}

    def addNewEid(self,eid,axial,SMa,torsion,SMt):
        #print "Rod Stress add..."
        self.axial[eid]      = axial
        self.MS_axial[eid]   = SMa
        self.torsion[eid]    = torsion
        self.MS_torsion[eid] = SMt

    def addNewEidTransient(self,eid,axial,SMa,torsion,SMt):
        #print "Rod Stress add..."
        dt = self.dt
        self.axial[dt][eid]      = axial
        self.MS_axial[dt][eid]   = SMa
        self.torsion[dt][eid]    = torsion
        self.MS_torsion[dt][eid] = SMt

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---ROD STRESSES---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['axial','torsion','MS_axial','MS_torsion']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for eid in sorted(self.axial):
            axial   = self.axial[eid]
            torsion = self.torsion[eid]
            SMa     = self.MS_axial[eid]
            SMt     = self.MS_torsion[eid]
            msg += '%-6i %6s ' %(eid,self.eType)
            vals = [axial,torsion,SMa,SMt]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %('0')
                else:
                    msg += '%10i ' %(val)
                ###
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

class rodStrainObject(scalarObject):
    """
                                     S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )
    ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
      ID.        STRAIN       MARGIN        STRAIN      MARGIN
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CROD'
        print "self.sCode = ",self.sCode
        if self.sCode in [1]:
            self.axial      = {}
            self.MS_axial   = {}
            self.torsion    = {}
            self.MS_torsion = {}
        else:
            raise Exception('get the sCode')
        ###
        raise Exception('get the sCode')
        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.axial[self.dt]      = {}
        self.MS_axial[self.dt]   = {}
        self.torsion[self.dt]    = {}
        self.MS_torsion[self.dt] = {}

    def addNewEid(self,axial,SMa,torsion,SMt):
        #print "Rod Strain add..."
        self.eType[eid] = self.eType
        self.axial[eid]      = axial
        self.MS_axial[eid]   = SMa
        self.torsion[eid]    = torsion
        self.MS_torsion[eid] = SMt

    def addNewEidTransient(self,axial,SMa,torsion,SMt):
        #print "Rod Strain add..."
        dt = self.dt
        self.eType[dt][eid] = self.eType
        self.axial[dt][eid]      = axial
        self.MS_axial[dt][eid]   = SMa
        self.torsion[dt][eid]    = torsion
        self.MS_torsion[dt][eid] = SMt

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---ROD STRAINS---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['axial','torsion','MS_tension','MS_compression']
        for header in headers:
            msg += '%8s ' %(header)
        msg += '\n'

        for eid in sorted(self.axial):
            axial   = self.axial[eid]
            torsion = self.torsion[eid]
            SMa     = self.MS_axial[eid]
            SMt     = self.MS_tension[eid]
            msg += '%-6i %6s ' %(eid,self.eType)
            vals = [axial,torsion,SMa,SMt]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%8s ' %('0')
                else:
                    msg += '%8i ' %(val)
                ###
            msg += '\n'
        return msg

class barStressObject(scalarObject):
    """
    # sCode=0
                               S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )
    ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
      ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.eType = {}
        print "self.sCode = ",self.sCode
        if self.sCode in [0]:
            self.s1    = {}
            self.s2    = {}
            self.s3    = {}
            self.s4    = {}
            self.axial = {}
            self.smax  = {}
            self.smin  = {}
            self.MS_tension = {}
            self.MS_compression = {}
        else:
            raise Exception('get the sCode')
        ###
        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.s1[self.dt]    = {}
        self.s2[self.dt]    = {}
        self.s3[self.dt]    = {}
        self.s4[self.dt]    = {}
        self.axial[self.dt] = {}
        self.smax[self.dt]  = {}
        self.smin[self.dt]  = {}
        #self.MS_tension[self.dt]     = {}
        #self.MS_compression[self.dt] = {}

    def addNewEid(self,eType,eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                                 s1b,s2b,s3b,s4b,      smaxb,sminb,MSc):
        #print "Bar Stress add..."
        self.eType[eid] = eType
        self.s1[eid]    = [s1a,s1b]
        self.s2[eid]    = [s2a,s2b]
        self.s3[eid]    = [s3a,s3b]
        self.s4[eid]    = [s4a,s4b]
        self.axial[eid] = axial
        self.smax[eid]  = [smaxa,smaxb]
        self.smin[eid]  = [smina,sminb]
        #self.MS_tension[eid]     = MSt
        #self.MS_compression[eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                                          s1b,s2b,s3b,s4b,      smaxb,sminb,MSc):
        #print "Bar Stress add..."
        dt = self.dt
        self.eType[eid] = eType
        self.s1[dt][eid]    = [s1a,s1b]
        self.s2[dt][eid]    = [s2a,s2b]
        self.s3[dt][eid]    = [s3a,s3b]
        self.s4[dt][eid]    = [s4a,s4b]
        self.axial[dt][eid] = axial
        self.smax[dt][eid]  = [smaxa,smaxb]
        self.smin[dt][eid]  = [smina,sminb]
        #self.MS_tension[dt][eid]     = MSt
        #self.MS_compression[dt][eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def __reprTransient__(self):
        msg = '---BAR STRESS---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['s1','s2','s3','s4','Axial','sMax','sMin']
        for header in headers:
            msg += '%6s ' %(header)
        msg += '\n'

        for dt,S1ss in sorted(self.s1.items()):
            msg += 'dt = %s' %(self.dt)
            for eid,S1s in sorted(S1ss.items()):
                eType = self.eType[eid]
                axial = self.axial[dt][eid]
                #MSt = self.MSt[dt][eid]
                #MSc = self.MSc[dt][eid]

                s1   = self.s1[dt][eid]
                s2   = self.s2[dt][eid]
                s3   = self.s3[dt][eid]
                s4   = self.s4[dt][eid]
                smax = self.smax[dt][eid]
                smin = self.smin[dt][eid]
                msg += '%-6i %6s ' %(eid,eType)
                vals = [s1[0],s2[0],s3[0],s4[0],axial,smax[0],smin[0]]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%6s ' %('0')
                    else:
                        msg += '%6i ' %(val)
                    ###
                msg += '\n'

                msg += '%s ' %(' '*13)
                vals = [s1[1],s2[1],s3[1],s4[1],'',smax[1],smin[1]]
                for val in vals:
                    if isinstance(val,str):
                        msg += '%6s ' %(val)
                    elif abs(val)<1e-6:
                        msg += '%6s ' %('0')
                    else:
                        msg += '%6i ' %(val)
                    ###
                msg += '\n'


                #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i smax=%-5i smax=%-4i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial, smax[0],smin[0])
                #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s         smax=%-5i smax=%-4i\n" %(' '*4,    s1[1],s2[1],s3[1],s4[1],'    ',smax[1],smin[1])
            ###
        ###
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---BAR STRESS---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['s1','s2','s3','s4','Axial','sMax','sMin']
        for header in headers:
            msg += '%6s ' %(header)
        msg += '\n'

        for eid,S1s in sorted(self.s1.items()):
            eType = self.eType[eid]
            axial = self.axial[eid]
            #MSt = self.MSt[eid]
            #MSc = self.MSc[eid]

            s1   = self.s1[eid]
            s2   = self.s2[eid]
            s3   = self.s3[eid]
            s4   = self.s4[eid]
            smax = self.smax[eid]
            smin = self.smin[eid]
            msg += '%-6i %6s ' %(eid,eType)
            vals = [s1[0],s2[0],s3[0],s4[0],axial,smax[0],smin[0]]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%6s ' %('0')
                else:
                    msg += '%6i ' %(val)
                ###
            msg += '\n'

            msg += '%s ' %(' '*13)
            vals = [s1[1],s2[1],s3[1],s4[1],'',smax[1],smin[1]]
            for val in vals:
                if isinstance(val,str):
                    msg += '%6s ' %(val)
                elif abs(val)<1e-6:
                    msg += '%6s ' %('0')
                else:
                    msg += '%6i ' %(val)
                ###
            msg += '\n'


            #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i smax=%-5i smax=%-4i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial, smax[0],smin[0])
            #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s         smax=%-5i smax=%-4i\n" %(' '*4,    s1[1],s2[1],s3[1],s4[1],'    ',smax[1],smin[1])
        ###
        return msg


class barStrainObject(scalarObject):
    """
    # sCode=10
                                    S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )
    ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
      ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C

    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.eType = {}
        print "self.sCode = ",self.sCode
        if self.sCode in [1]:
            raise Exception('verify...')
            self.e1    = {}
            self.e2    = {}
            self.e3    = {}
            self.e4    = {}
            self.axial = {}
            self.emax  = {}
            self.emin  = {}
            #self.MS_tension = {}
            #self.MS_compression = {}
        elif self.sCode in [10]:
            self.e1    = {}
            self.e2    = {}
            self.e3    = {}
            self.e4    = {}
            self.axial = {}
            self.emax  = {}
            self.emin  = {}
            self.MS_tension = {}
            self.MS_compression = {}
        else:
            raise Exception('get the sCode')
        ###
        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.e1[self.dt]    = {}
        self.e2[self.dt]    = {}
        self.e3[self.dt]    = {}
        self.e4[self.dt]    = {}
        self.axial[self.dt] = {}
        self.emax[self.dt]  = {}
        self.emin[self.dt]  = {}
        #self.MS_tension[self.dt]     = {}
        #self.MS_compression[self.dt] = {}

    def addNewEid(self,eType,eid,e1a,e2a,e3a,e4a,axial,emaxa,emina,MSt,
                                 e1b,e2b,e3b,e4b,      emaxb,eminb,MSc):
        #print "Bar Stress add..."
        self.eType[eid] = eType
        self.e1[eid]    = [e1a,e1b]
        self.e2[eid]    = [e2a,e2b]
        self.e3[eid]    = [e3a,e3b]
        self.e4[eid]    = [e4a,e4b]
        self.axial[eid] = axial
        self.emax[eid]  = [emaxa,emaxb]
        self.emin[eid]  = [emina,eminb]
        #self.MS_tension[eid]     = MSt
        #self.MS_compression[eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,e1a,e2a,e3a,e4a,axial,emaxa,emina,MSt,
                                      e1b,e2b,e3b,e4b,      emaxb,eminb,MSc):
        #print "Bar Stress add..."
        dt = self.dt
        self.eType[dt][eid] = eType
        self.e1[dt][eid]    = [e1a,e1b]
        self.e2[dt][eid]    = [e2a,e2b]
        self.e3[dt][eid]    = [e3a,e3b]
        self.e4[dt][eid]    = [e4a,e4b]
        self.axial[dt][eid] = axial
        self.emax[dt][eid]  = [emaxa,emaxb]
        self.emin[dt][eid]  = [emina,eminb]
        #self.MS_tension[dt][eid]     = MSt
        #self.MS_compression[dt][eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---BAR STRAIN---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['e1','e2','e3','e4','Axial','eMax','eMin']
        for header in headers:
            msg += '%6s ' %(header)
        msg += '\n'

        for eid,E1s in sorted(self.e1.items()):
            eType = self.eType[eid]
            axial = self.axial[eid]
            #MSt  = self.MS_tension[eid]
            #MSc  = self.MS_compression[eid]

            e1   = self.e1[eid]
            e2   = self.e2[eid]
            e3   = self.e3[eid]
            e4   = self.e4[eid]
            emax = self.emax[eid]
            emin = self.emin[eid]
            msg += '%-6i %6s ' %(eid,eType)
            vals = [e1[0],e2[0],e3[0],e4[0],axial,emax[0],emin[0]]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%6s ' %('0')
                else:
                    msg += '%6i ' %(val)
                ###
            msg += '\n'

            msg += '%s ' %(' '*13)
            vals = [e1[1],e2[1],e3[1],e4[1],'',emax[1],emin[1]]
            for val in vals:
                if isinstance(val,str):
                    msg += '%6s ' %(val)
                elif abs(val)<1e-6:
                    msg += '%6s ' %('0')
                else:
                    msg += '%6i ' %(val)
                ###
            msg += '\n'


            #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i smax=%-5i smax=%-4i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial, smax[0],smin[0])
            #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s         smax=%-5i smax=%-4i\n" %(' '*4,    s1[1],s2[1],s3[1],s4[1],'    ',smax[1],smin[1])
        ###
        return msg


class solidStressObject(scalarObject):
    """
                          S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )
                   CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   
    ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES 
            1           0GRID CS  8 GP
                   CENTER  X   4.499200E+02  XY  -5.544791E+02   A   1.000000E+04  LX 0.00 0.69-0.72  -3.619779E+03    9.618462E+03
                           Y   4.094179E+02  YZ   5.456968E-12   B  -1.251798E+02  LY 0.00 0.72 0.69
                           Z   1.000000E+04  ZX  -4.547474E-13   C   9.845177E+02  LZ 1.00 0.00 0.00

    # sCode=1
                     S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )
                   CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   
    ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES 

    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)

        self.eType = {}
        print "self.sCode = ",self.sCode
        if self.sCode in [1]:
            self.cid = {}
            self.oxx = {}
            self.oyy = {}
            self.ozz = {}
            self.txy = {}
            self.tyz = {}
            self.txz = {}
            #self.aCos = {}
            #self.bCos = {}
            #self.cCos = {}
            #self.pressure = {}
            self.ovm = {}
        else:
            raise Exception('get the sCode')

        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.add       = self.addTransient
            self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.oxx[self.dt] = {}
        self.oyy[self.dt] = {}
        self.ozz[self.dt] = {}
        self.txy[self.dt] = {}
        self.tyz[self.dt] = {}
        self.txz[self.dt] = {}
        
        #self.aCos[self.dt] = {}
        #self.bCos[self.dt] = {}
        #self.cCos[self.dt] = {}
        #self.pressure[self.dt] = {}
        self.ovm[self.dt]      = {}

    def addNewEid(self,eType,cid,eid,nodeID,oxx,oyy,ozz,txy,tyz,txz,aCos,bCos,cCos,pressure,ovm):
        #print "Solid Stress add..."
        assert cid >= 0
        assert eid >= 0
        self.eType[eid] = eType
        self.cid[eid]  = cid
        self.oxx[eid]  = {nodeID: oxx}
        self.oyy[eid]  = {nodeID: oyy}
        self.ozz[eid]  = {nodeID: ozz}
        self.txy[eid]  = {nodeID: txy}
        self.tyz[eid]  = {nodeID: tyz}
        self.txz[eid]  = {nodeID: txz}
        #self.aCos[eid] = {nodeID: aCos}
        #self.bCos[eid] = {nodeID: bCos}
        #self.cCos[eid] = {nodeID: cCos}
        #self.pressure[eid] = {nodeID: pressure}
        self.ovm[eid]      = {nodeID: ovm}
        msg = "*eid=%s nodeID=%s vm=%g" %(eid,nodeID,ovm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,cid,eid,nodeID,oxx,oyy,ozz,txy,tyz,txz,aCos,bCos,cCos,pressure,ovm):
        #print "Solid Stress add transient..."
        assert cid >= 0
        assert eid >= 0
        dt = self.dt
        self.eType[eid] = eType
        self.cid[eid]   = cid
        self.oxx[dt][eid]  = {nodeID: oxx}
        self.oyy[dt][eid]  = {nodeID: oyy}
        self.ozz[dt][eid]  = {nodeID: ozz}
        self.txy[dt][eid]  = {nodeID: txy}
        self.tyz[dt][eid]  = {nodeID: tyz}
        self.txz[dt][eid]  = {nodeID: txz}
        #self.aCos[dt][eid] = {nodeID: aCos}
        #self.bCos[dt][eid] = {nodeID: bCos}
        #self.cCos[dt][eid] = {nodeID: cCos}
        #self.pressure[dt][eid] = {nodeID: pressure}
        self.ovm[dt][eid]      = {nodeID: ovm}
        msg = "*eid=%s nodeID=%s vm=%g" %(eid,nodeID,ovm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,oxx,oyy,ozz,txy,tyz,txz,aCos,bCos,cCos,pressure,ovm):
        #print "***add"
        msg = "eid=%s nodeID=%s vm=%g" %(eid,nodeID,ovm)
        #print msg
        #print self.oxx
        #print self.fiberDistance
        self.oxx[eid][nodeID] = oxx
        self.oyy[eid][nodeID] = oyy
        self.ozz[eid][nodeID] = ozz

        self.txy[eid][nodeID] = txy
        self.tyz[eid][nodeID] = tyz
        self.txz[eid][nodeID] = txz

        #self.aCos[eid][nodeID] = aCos
        #self.bCos[eid][nodeID] = bCos
        #self.cCos[eid][nodeID] = cCos
        #self.pressure[eid][nodeID] = pressure
        self.ovm[eid][nodeID] = ovm
        if nodeID==0: raise Exception(msg)

    def addTransient(self,eid,nodeID,oxx,oyy,ozz,txy,tyz,txz,aCos,bCos,cCos,pressure,ovm):
        #print "***add"
        msg = "eid=%s nodeID=%s vm=%g" %(eid,nodeID,ovm)
        #print msg
        #print self.oxx
        #print self.fiberDistance
        dt = self.dt
        self.oxx[dt][eid][nodeID] = oxx
        self.oyy[dt][eid][nodeID] = oyy
        self.ozz[dt][eid][nodeID] = ozz

        self.txy[dt][eid][nodeID] = txy
        self.tyz[dt][eid][nodeID] = tyz
        self.txz[dt][eid][nodeID] = txz

        #self.aCos[dt][eid][nodeID] = aCos
        #self.bCos[dt][eid][nodeID] = bCos
        #self.cCos[dt][eid][nodeID] = cCos
        #self.pressure[dt][eid][nodeID] = pressure
        self.ovm[dt][eid][nodeID] = ovm
        if nodeID==0: raise Exception(msg)

    def __reprTransient__(self):
        msg = '---SOLID STRESS---\n'
        headers = ['oxx','oyy','ozz','txy','tyz','txz','ovm']
        msg += '%-6s %6s %8s ' %('EID','eType','nodeID')
        for header in headers:
            msg += '%9s ' %(header)
        msg += '\n'

        for dt,oxxs in sorted(self.oxx.items()):
            msg += 'dt = %g\n' %(dt)
            for eid,oxxNodes in sorted(oxxs.items()):
                eType = self.eType[eid]
                for nid in sorted(oxxNodes):
                    oxx = self.oxx[dt][eid][nid]
                    oyy = self.oyy[dt][eid][nid]
                    ozz = self.ozz[dt][eid][nid]
                    txy = self.txy[dt][eid][nid]
                    tyz = self.tyz[dt][eid][nid]
                    txz = self.txz[dt][eid][nid]
                    ovm = self.ovm[dt][eid][nid]
                    msg += '%-6i %6s %8s ' %(eid,eType,nid)
                    vals = [oxx,oyy,ozz,txy,tyz,txz,ovm]
                    for val in vals:
                        if abs(val)<1e-6:
                            msg += '%9s ' %('0')
                        else:
                            msg += '%9i ' %(val)
                        ###
                    msg += '\n'
                    #msg += "eid=%-4s eType=%-6s nid=%-2i oxx=%-5i oyy=%-5i ozz=%-5i txy=%-5i tyz=%-5i txz=%-5i ovm=%-5i\n" %(eid,eType,nid,oxx,oyy,ozz,txy,tyz,txz,ovm)
                ###
            ###
        ###
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---SOLID STRESS---\n'
        headers = ['oxx','oyy','ozz','txy','tyz','txz','ovm']
        msg += '%-6s %6s %8s ' %('EID','eType','nodeID')
        for header in headers:
            msg += '%9s ' %(header)
        msg += '\n'
        for eid,oxxNodes in sorted(self.oxx.items()):
            eType = self.eType[eid]
            for nid in sorted(oxxNodes):
                oxx = self.oxx[eid][nid]
                oyy = self.oyy[eid][nid]
                ozz = self.ozz[eid][nid]
                txy = self.txy[eid][nid]
                tyz = self.tyz[eid][nid]
                txz = self.txz[eid][nid]
                ovm = self.ovm[eid][nid]
                msg += '%-6i %6s %8s ' %(eid,eType,nid)
                vals = [oxx,oyy,ozz,txy,tyz,txz,ovm]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%9s ' %('0')
                    else:
                        msg += '%9i ' %(val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%-6s nid=%-2i oxx=%-5i oyy=%-5i ozz=%-5i txy=%-5i tyz=%-5i txz=%-5i ovm=%-5i\n" %(eid,eType,nid,oxx,oyy,ozz,txy,tyz,txz,ovm)
            ###
        ###
        return msg

class solidStrainObject(scalarObject):
    """
                          S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )
                   CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   
    ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES 
            1           0GRID CS  8 GP
                   CENTER  X   4.499200E+02  XY  -5.544791E+02   A   1.000000E+04  LX 0.00 0.69-0.72  -3.619779E+03    9.618462E+03
                           Y   4.094179E+02  YZ   5.456968E-12   B  -1.251798E+02  LY 0.00 0.72 0.69
                           Z   1.000000E+04  ZX  -4.547474E-13   C   9.845177E+02  LZ 1.00 0.00 0.00

    """
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.eType = {}
        self.cid = {}
        self.exx = {}
        self.eyy = {}
        self.ezz = {}
        self.exy = {}
        self.eyz = {}
        self.exz = {}
        #self.aCos = {}
        #self.bCos = {}
        #self.cCos = {}
        #self.pressure = {}
        self.evm = {}
        print "self.sCode = ",self.sCode
        raise Exception('get the sCode')
        if dt is not None:
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
        self.exx[self.dt] = {}
        self.eyy[self.dt] = {}
        self.ezz[self.dt] = {}
        self.exy[self.dt] = {}
        self.eyz[self.dt] = {}
        self.exz[self.dt] = {}
        
        #self.aCos[self.dt] = {}
        #self.bCos[self.dt] = {}
        #self.cCos[self.dt] = {}
        #self.pressure[self.dt] = {}
        self.evm[self.dt]      = {}

    def addNewEid(self,eType,cid,eid,nodeID,exx,eyy,ezz,exy,eyz,exz,aCos,bCos,cCos,pressure,evm):
        #print "Solid Strain add..."
        assert cid >= 0
        assert eid >= 0
        self.eType[eid] = eType
        self.cid[eid]  = cid
        self.exx[eid]  = {nodeID: exx}
        self.eyy[eid]  = {nodeID: eyy}
        self.ezz[eid]  = {nodeID: ezz}
        self.exy[eid]  = {nodeID: exy}
        self.eyz[eid]  = {nodeID: eyz}
        self.exz[eid]  = {nodeID: exz}
        #self.aCos[eid] = {nodeID: aCos}
        #self.bCos[eid] = {nodeID: bCos}
        #self.cCos[eid] = {nodeID: cCos}
        #self.pressure[eid] = {nodeID: pressure}
        self.evm[eid]      = {nodeID: evm}
        msg = "*eid=%s nodeID=%s vm=%g" %(eid,nodeID,evm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,cid,eid,nodeID,exx,eyy,ezz,exy,eyz,exz,aCos,bCos,cCos,pressure,evm):
        #print "Solid Strain add..."
        assert cid >= 0
        assert eid >= 0
        self.eType[eid] = eType
        self.cid[eid]  = cid
        dt = self.dt
        self.exx[dt][eid]  = {nodeID: exx}
        self.eyy[dt][eid]  = {nodeID: eyy}
        self.ezz[dt][eid]  = {nodeID: ezz}
        self.exy[dt][eid]  = {nodeID: exy}
        self.eyz[dt][eid]  = {nodeID: eyz}
        self.exz[dt][eid]  = {nodeID: exz}
        #self.aCos[dt][eid] = {nodeID: aCos}
        #self.bCos[dt][eid] = {nodeID: bCos}
        #self.cCos[dt][eid] = {nodeID: cCos}
        #self.pressure[dt][eid] = {nodeID: pressure}
        self.evm[dt][eid]      = {nodeID: evm}
        msg = "*eid=%s nodeID=%s vm=%g" %(eid,nodeID,evm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,exx,eyy,ezz,exy,eyz,exz,aCos,bCos,cCos,pressure,evm):
        #print "***add"
        msg = "eid=%s nodeID=%s vm=%g" %(eid,nodeID,evm)
        #print msg
        #print self.exx
        #print self.fiberDistance
        self.exx[eid][nodeID] = exx
        self.eyy[eid][nodeID] = eyy
        self.ezz[eid][nodeID] = ezz

        self.exy[eid][nodeID] = exy
        self.eyz[eid][nodeID] = eyz
        self.exz[eid][nodeID] = exz

        #self.aCos[eid][nodeID] = aCos
        #self.bCos[eid][nodeID] = bCos
        #self.cCos[eid][nodeID] = cCos
        #self.pressure[eid][nodeID] = pressure
        self.evm[eid][nodeID] = evm

        if nodeID==0: raise Exception(msg)

    def addTransient(self,eid,nodeID,exx,eyy,ezz,exy,eyz,exz,aCos,bCos,cCos,pressure,evm):
        #print "***add"
        msg = "eid=%s nodeID=%s vm=%g" %(eid,nodeID,evm)
        #print msg
        #print self.exx
        #print self.fiberDistance
        dt = self.dt
        self.exx[dt][eid][nodeID] = exx
        self.eyy[dt][eid][nodeID] = eyy
        self.ezz[dt][eid][nodeID] = ezz

        self.exy[dt][eid][nodeID] = exy
        self.eyz[dt][eid][nodeID] = eyz
        self.exz[dt][eid][nodeID] = exz

        #self.aCos[dt][eid][nodeID] = aCos
        #self.bCos[dt][eid][nodeID] = bCos
        #self.cCos[dt][eid][nodeID] = cCos
        #self.pressure[dt][eid][nodeID] = pressure
        self.evm[dt][eid][nodeID] = evm

        if nodeID==0: raise Exception(msg)

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---SOLID STRAIN---\n'
        headers = ['exx','eyy','ezz','exy','eyz','exz','evm']
        msg += '%-6s %6s %8s ' %('EID','eType','nodeID')
        for header in headers:
            msg += '%9s ' %(header)
        msg += '\n'
        for eid,exxNodes in sorted(self.exx.items()):
            eType = self.eType[eid]
            for nid in sorted(exxNodes):
                exx = self.exx[eid][nid]
                eyy = self.eyy[eid][nid]
                ezz = self.ezz[eid][nid]
                exy = self.exy[eid][nid]
                eyz = self.eyz[eid][nid]
                exz = self.exz[eid][nid]
                evm = self.evm[eid][nid]
                msg += '%-6i %6s %8s ' %(eid,eType,nid)
                vals = [exx,eyy,ezz,exy,eyz,exz,evm]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%9s ' %('0')
                    else:
                        msg += '%9e ' %(val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%-6s nid=%-2i exx=%-5i eyy=%-5i ezz=%-5i exy=%-5i eyz=%-5i exz=%-5i evm=%-5i\n" %(eid,eType,nid,exx,eyy,ezz,exy,eyz,exz,evm)
            ###
        ###
        return msg

