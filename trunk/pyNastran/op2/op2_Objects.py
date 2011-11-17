from numpy import array

class scalarObject(object):
    def __init__(self,iSubcase):
        self.iSubcase = iSubcase
    def add6(self,nodeID,v1,v2,v3,v4,v5,v6):
        pass
    def add3(self,nodeID,v1,v2,v3):
        pass

class displacementObject(scalarObject):
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.displacements = {}
        self.rotations     = {}

    def add(self,nodeID,v1,v2,v3,v4,v5,v6):
        self.displacements[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]     = array([v4,v5,v6]) # rx,ry,rz


    def __repr__(self):
        msg = '---DISPLACEMENTS---\n'
        headers = ['Dx','Dy','Dz','Rx','Ry','Rz']
        msg += '%9s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,displacement in sorted(self.displacements.items()):
            rotation = self.rotations[nodeID]
            (dx,dy,dz) = displacement
            (rx,ry,rz) = rotation

            msg += '%9i ' %(nodeID)
            vals = [dx,dy,dz,rx,ry,rz]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.3e ' %(val)
                ###
            msg += '\n'
        return msg


class spcForcesObject(scalarObject):
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.forces  = {}
        self.moments = {}

    def add(self,nodeID,v1,v2,v3,v4,v5,v6):
        self.forces[nodeID]  = array([v1,v2,v3]) # Fx,Fy,Fz
        self.moments[nodeID] = array([v4,v5,v6]) # Mx,My,Mz

    def __repr__(self):
        msg = '---SPC FORCES---\n'

        headers = ['Fx','Fy','Fz','Mx','My','Mz']
        msg += '%9s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,force in sorted(self.forces.items()):
            moment = self.moments[nodeID]
            (Fx,Fy,Fz) = force
            (Mx,My,Mz) = moment

            msg += '%9i ' %(nodeID)
            vals = [Fx,Fy,Fz,Mx,My,Mx]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.3e ' %(val)
                ###
            msg += '\n'
        return msg

class temperatureObject(scalarObject):
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.temperatures = {}

    def add(self,nodeID,v1,v2=None,v3=None,v4=None,v5=None,v6=None):
        self.temperatures[nodeID] = v1

    def __repr__(self):
        msg = '---TEMPERATURE---\n'
        msg += '%-8s %-8s\n' %('GRID','TEMPERATURE')
        for nodeID,T in sorted(self.temperatures.items()):
            msg += '%9i ' %(nodeID)

            if abs(val)<1e-6:
                msg += '%10s' %(0)
            else:
                msg += '%10i\n' %(val)
            ###
        return msg

class fluxObject(scalarObject):
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.fluxs = {}

    def add(self,nodeID,v1,v2,v3,v4=None,v5=None,v6=None):
        fluxs[nodeID] = array([v1,v2,v3])

    def __repr__(self):
        msg = '---HEAT FLUX---\n'
        msg += '%-8s %-8s %-8s %-8s\n' %('GRID','Fx','Fy','Fz')
        for nodeID,flux in sorted(self.fluxs.items()):
            msg += '%9i ' %(nodeID)

            for val in flux:
                if abs(val)<1e-6:
                    msg += '%10s' %(0)
                else:
                    msg += '%10.3e ' %(val)
                ###
            msg += '\n'
        return msg

class strainObject(scalarObject):
    """
    ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 
      ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
    """
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.eType = {}
        self.curvature = {}
        self.exx = {}
        self.eyy = {}
        self.exy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        self.evm = {}

    def addNewEid(self,eType,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        self.eType[eid] = eType
        self.curvature[eid] = {nodeID: [curvature]}
        self.exx[eid] = {nodeID: [exx]}
        self.eyy[eid] = {nodeID: [eyy]}
        self.exy[eid] = {nodeID: [exy]}
        self.angle[eid] = {nodeID: [angle]}
        self.majorP[eid] = {nodeID: [majorP]}
        self.minorP[eid] = {nodeID: [minorP]}
        self.evm[eid]    = {nodeID: [evm]}
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm):
        #print "***add"
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        #print msg
        #print self.oxx
        #print self.fiberDistance
        self.curvature[eid][nodeID].append(curvature)
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
        self.curvature[eid][nodeID] = [curvature]
        self.exx[eid][nodeID] = [exx]
        self.eyy[eid][nodeID] = [eyy]
        self.exy[eid][nodeID] = [exy]
        self.angle[eid][nodeID] = [angle]
        self.majorP[eid][nodeID] = [majorP]
        self.minorP[eid][nodeID] = [minorP]
        self.evm[eid][nodeID] = [evm]
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def __repr__(self):
        msg = ''
        for eid,exxNodes in sorted(self.exx.items()):
            eType = self.eType[eid]
            for nid in sorted(exxNodes):
                for iLayer in range(len(self.exx[eid][nid])):
                    fd    = self.curvature[eid][nid][iLayer]
                    exx = self.exx[eid][nid][iLayer]
                    eyy = self.eyy[eid][nid][iLayer]
                    exy = self.exy[eid][nid][iLayer]
                    angle = self.angle[eid][nid][iLayer]
                    major = self.majorP[eid][nid][iLayer]
                    minor = self.minorP[eid][nid][iLayer]
                    evm = self.evm[eid][nid][iLayer]
                    print "eid=%s eType=%s nid=%s iLayer=%s exx=%-9.3g eyy=%-9.3g exy=%-9.3g evm=%-9.3g" %(eid,eType,nid,iLayer,exx,eyy,exy,evm)
            ###
        ###
        return msg

class stressObject(scalarObject):
    """
    ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 
      ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
    """
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.eType = {}
        self.fiberDistance = {}
        self.oxx = {}
        self.oyy = {}
        self.txy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        self.ovm = {}

    def addNewEid(self,eType,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        self.eType[eid] = eType
        self.fiberDistance[eid] = {nodeID: [fd]}
        self.oxx[eid] = {nodeID: [oxx]}
        self.oyy[eid] = {nodeID: [oyy]}
        self.txy[eid] = {nodeID: [txy]}
        self.angle[eid] = {nodeID: [angle]}
        self.majorP[eid] = {nodeID: [majorP]}
        self.minorP[eid] = {nodeID: [minorP]}
        self.ovm[eid]    = {nodeID: [ovm]}
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***add"
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #print self.oxx
        #print self.fiberDistance
        self.fiberDistance[eid][nodeID].append(fd)
        self.oxx[eid][nodeID].append(oxx)
        self.oyy[eid][nodeID].append(oyy)
        self.txy[eid][nodeID].append(txy)
        self.angle[eid][nodeID].append(angle)
        self.majorP[eid][nodeID].append(majorP)
        self.minorP[eid][nodeID].append(minorP)
        self.ovm[eid][nodeID].append(ovm)
        if nodeID==0: raise Exception(msg)

    def addNewNode(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print "***addNewNode"
        #print self.oxx
        self.fiberDistance[eid][nodeID] = [fd]
        self.oxx[eid][nodeID] = [oxx]
        self.oyy[eid][nodeID] = [oyy]
        self.txy[eid][nodeID] = [txy]
        self.angle[eid][nodeID] = [angle]
        self.majorP[eid][nodeID] = [majorP]
        self.minorP[eid][nodeID] = [minorP]
        self.ovm[eid][nodeID] = [ovm]
        msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def __repr__(self):
        msg = ''
        for eid,oxxNodes in sorted(self.oxx.items()):
            eType = self.eType[eid]
            for nid in sorted(oxxNodes):
                for iLayer in range(len(self.oxx[eid][nid])):
                    fd    = self.fiberDistance[eid][nid][iLayer]
                    oxx = self.oxx[eid][nid][iLayer]
                    oyy = self.oyy[eid][nid][iLayer]
                    txy = self.txy[eid][nid][iLayer]
                    angle = self.angle[eid][nid][iLayer]
                    major = self.majorP[eid][nid][iLayer]
                    minor = self.minorP[eid][nid][iLayer]
                    ovm = self.ovm[eid][nid][iLayer]
                    print "eid=%s eType=%s nid=%s iLayer=%s oxx=%-4i oyy=%-4i txy=%-4i ovm=%-4i" %(eid,eType,nid,iLayer,oxx,oyy,txy,ovm)
            ###
        ###
        return msg
