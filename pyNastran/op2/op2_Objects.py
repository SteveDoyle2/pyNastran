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
        msg += '%9s  %-9s %-9s %-9s %-9s %-9s %-9s\n' %('GRID','Dx','Dy','Dz','Rx','Ry','Rz')
        for nodeID,displacement in sorted(self.displacements.items()):
            rotation = self.rotations[nodeID]
            (dx,dy,dz) = displacement
            (rx,ry,rz) = rotation
            if abs(dx)<1e-5:  dx=0
            if abs(dy)<1e-5:  dy=0
            if abs(dz)<1e-5:  dz=0

            if abs(rx)<1e-5:  rx=0.
            if abs(ry)<1e-5:  ry=0.
            if abs(rz)<1e-5:  rz=0.
            msg += '%9i %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n' %(nodeID,dx,dy,dz,rx,ry,rz)
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
        msg += '%-8s %-8s %-8s %-8s %-8s %-8s %-8s\n' %('GRID','Fx','Fy','Fz','Mx','My','Mz')
        for nodeID,force in sorted(self.forces.items()):
            moment = self.moments[nodeID]
            (Fx,Fy,Fz) = force
            (Mx,My,Mz) = moment
            if abs(Fx)<1e-5:  Fx=0
            if abs(Fy)<1e-5:  Fy=0
            if abs(Fz)<1e-5:  Fz=0

            if abs(Mx)<1e-5:  Mx=0.
            if abs(My)<1e-5:  My=0.
            if abs(Mz)<1e-5:  Mz=0.
            msg += '%-8g %-8g %-8g %-8g %-8g %-8g %-8g\n' %(nodeID,Fx,Fy,Fz,Mx,My,Mz)
        return msg

class temperatureObject(scalarObject):
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.temperatures = {}

    def add(self,nodeID,v1,v2=None,v3=None,v4=None,v5=None,v6=None):
        self.temperatures[nodeID] = v1

class fluxObject(scalarObject):
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.fluxs = {}

    def add(self,nodeID,v1,v2,v3,v4=None,v5=None,v6=None):
        fluxs[nodeID] = array([v1,v2,v3])

