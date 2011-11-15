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

#class strainObject(scalarObject):
#    """
#    ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 
#      ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
#    """
#    def __init__(self,iSubcase):
#        scalarObject.__init__(self,iSubcase)
#        self.curvature = {}
#        self.exx = {}
#        self.eyy = {}
#        self.exy = {}
#        self.angle = {}
#        self.majorP = {}
#        self.minorP = {}
#        self.evm = {}

class stressObject(scalarObject):
    """
    ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 
      ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
    """
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.fiberDistance = {}
        self.oxx = {}
        self.oyy = {}
        self.txy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        self.ovm = {}

    def addNewEid(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        self.fiberDistance[eid] = {nodeID: fd}
        self.oxx[eid] = {nodeID: oxx}
        self.oyy[eid] = {nodeID: oyy}
        self.txy[eid] = {nodeID: txy}
        self.angle[eid] = {nodeID: angle}
        self.majorP[eid] = {nodeID: majorP}
        self.minorP[eid] = {nodeID: minorP}
        self.ovm[eid]    = {nodeID: ovm}

    def add(self,eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm):
        #print self.oxx
        self.fiberDistance[eid][nodeID] = fd
        self.oxx[eid][nodeID] = oxx
        self.oyy[eid][nodeID] = oyy
        self.txy[eid][nodeID] = txy
        self.angle[eid][nodeID] = angle
        self.majorP[eid][nodeID] = majorP
        self.minorP[eid][nodeID] = minorP
        self.ovm[eid][nodeID] = ovm

    def __repr__(self):
        msg = ''
        for eid,oxxNodes in self.oxx.items():
            for nid in oxxNodes:
                fd    = self.fiberDistance[eid][nid]
                oxx = self.oxx[eid][nid]
                oyy = self.oyy[eid][nid]
                txy = self.txy[eid][nid]
                angle = self.angle[eid][nid]
                major = self.majorP[eid][nid]
                minor = self.minorP[eid][nid]
                ovm = self.ovm[eid][nid]
                print "eid=%s nid=%s oxx=%-4i oyy=%-4i txy=%-4i ovm=%-4i" %(eid,nid,oxx,oyy,txy,ovm)
            ###
        ###
        return msg
