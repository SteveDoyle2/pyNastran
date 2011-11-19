from numpy import array

class scalarObject(object):
    def __init__(self,iSubcase):
        self.iSubcase = iSubcase

class displacementObject(scalarObject):
    def __init__(self,iSubcase,dt=None):
        scalarObject.__init__(self,iSubcase)
        self.dt = dt
        
        ## this could get very bad very quick, but it could be great!
        ## basically it's a way to handle transients without making
        ## a whole new class
        if self.dt is None:
            self.displacements = {}
            self.rotations     = {}
        else:
            assert dt>=0.
            self.displacements = {dt: {}}
            self.rotations     = {dt: {}}
            self.add = self.addTransient
            self.addBinary = self.addBinaryTransient
            self.__repr__ = self.__reprTransient__  # why cant i do this...
        ###

    def updateDt(self,dt=None):
        """@todo not really sure how to use this yet..."""
        assert dt>=0.
        self.dt = dt
        self.displacements[dt] = {}

    def addBinary(self,deviceCode,data):
        (nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.displacements

        self.displacements[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addBinaryTransient(self,deviceCode,data):
        raise Exception('not implemented...')
        (nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.displacements

        self.displacements[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def add(self,nodeID,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.displacements

        self.displacements[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addTransient(self,nodeID,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.displacements
        
        self.displacements[self.dt][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[self.dt][nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def __reprTransient__(self):
        print "Transient..."
        raise Exception('this could be cool...')
        return self.__repr__()

    def __repr__(self):
        msg = '---DISPLACEMENTS---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)
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
    def __init__(self,iSubcase,dt=None):
        scalarObject.__init__(self,iSubcase)
        self.dt = dt

        if self.dt is None:
            self.forces  = {}
            self.moments = {}
        else:
            assert dt>=0.
            self.forces  = {dt: {}}
            self.moments = {dt: {}}
            self.add = self.addTransient
            self.__repr__ = self.__reprTransient__  # why cant i do this...
        ###

    def updateDt(self,dt=None):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        assert dt>=0.
        self.dt = dt
        self.forces[dt]    = {}
        self.momenents[dt] = {}

    #def addBinary(self,deviceCode,data):
    #    print "*******add********"
    #    (nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)

    def add(self,nodeID,v1,v2,v3,v4,v5,v6):
        msg = 'nodeID=%s' %(nodeID)
        assert 0<nodeID<1000000000,msg
        assert nodeID not in self.forces
        self.forces[ nodeID] = array([v1,v2,v3]) # Fx,Fy,Fz
        self.moments[nodeID] = array([v4,v5,v6]) # Mx,My,Mz

    def addTransient(self,nodeID,v1,v2,v3,v4,v5,v6):
        msg = 'nodeID=%s' %(nodeID)
        assert 0<nodeID<1000000000,msg
        assert nodeID not in self.forces
        self.forces[ self.dt][nodeID] = array([v1,v2,v3]) # Fx,Fy,Fz
        self.moments[self.dt][nodeID] = array([v4,v5,v6]) # Mx,My,Mz

    def __repr__(self):
        msg = '---SPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

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
    def __init__(self,iSubcase,dt=None):
        scalarObject.__init__(self,iSubcase)
        self.dt = dt
        self.temperatures = {}
        
        print "dt = ",self.dt
        if self.dt is None:
            self.temperatures = {}
        else:
            self.temperatures = {self.dt:{}}
            self.add = self.addTransient
            self.addBinary = self.addBinaryTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...            
        ###

    def updateDt(self,dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        assert dt>=0.
        self.dt = dt
        self.temperatures[dt] = {}

    def addBinary(self,deviceCode,data):
        (nodeID,v1) = unpack('if',data[0:8])
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.temperatures
        self.temperatures[nodeID] = v1

    def addBinaryTransient(self,deviceCode,data):
        raise Exception('not implemented...')
        (nodeID,v1) = unpack('if',data[0:8])
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.temperatures
        self.temperatures[nodeID] = v1

    def add(self,nodeID,v1,v2=None,v3=None,v4=None,v5=None,v6=None):
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.temperatures
        self.temperatures[nodeID] = v1

    def addTransient(self,nodeID,v1,v2=None,v3=None,v4=None,v5=None,v6=None):
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.temperatures
        self.temperatures[self.dt][nodeID] = v1

    def __reprTransient__(self):
        msg = '---TEMPERATURE---\n'
        msg += '%-8s %-8s\n' %('GRID','TEMPERATURE')

        for dt,temperatures in sorted(self.temperatures.items()):
            msg += 'dt = %g\n' %(dt)
            for nodeID,T in sorted(temperatures.items()):
                msg += '%9i ' %(nodeID)

                if abs(T)<1e-6:
                    msg += '%10s\n' %(0)
                else:
                    msg += '%10g\n' %(T)
                ###
            ###
        return msg

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---TEMPERATURE---\n'
        msg += '%-8s %-8s\n' %('GRID','TEMPERATURE')
        for nodeID,T in sorted(self.temperatures.items()):
            msg += '%9i ' %(nodeID)

            if abs(T)<1e-6:
                msg += '%10s\n' %(0)
            else:
                msg += '%10g\n' %(T)
            ###
        return msg

class fluxObject(scalarObject):
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.fluxs = {}

    def add(self,nodeID,v1,v2,v3,v4=None,v5=None,v6=None):
        assert 0>nodeID>1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.fluxs
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

class plateStrainObject(scalarObject):
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
        #print "Plate add..."
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        assert 0<nodeID<1000000000, 'nodeID=%s %s' %(nodeID,msg)
        assert eid not in self.exx
        self.eType[eid] = eType
        self.curvature[eid] = {nodeID: [curvature]}
        self.exx[eid] = {nodeID: [exx]}
        self.eyy[eid] = {nodeID: [eyy]}
        self.exy[eid] = {nodeID: [exy]}
        self.angle[eid] = {nodeID: [angle]}
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
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
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
        msg = "eid=%s nodeID=%s curvature=%g exx=%g eyy=%g \nexy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,curvature,exx,eyy,exy,angle,majorP,minorP,evm)
        assert nodeID not in self.exx[eid],msg
        self.curvature[eid][nodeID] = [curvature]
        self.exx[eid][nodeID] = [exx]
        self.eyy[eid][nodeID] = [eyy]
        self.exy[eid][nodeID] = [exy]
        self.angle[eid][nodeID] = [angle]
        self.majorP[eid][nodeID] = [majorP]
        self.minorP[eid][nodeID] = [minorP]
        self.evm[eid][nodeID] = [evm]
        #print msg
        if nodeID==0: raise Exception(msg)

    def __repr__(self):
        msg = '---PLATE STRAIN---\n'
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
                    msg += "eid=%s eType=%s nid=%s iLayer=%s exx=%-9.3g eyy=%-9.3g exy=%-9.3g evm=%-9.3g\n" %(eid,eType,nid,iLayer,exx,eyy,exy,evm)
                ###
            ###
        ###
        return msg

class plateStressObject(scalarObject):
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
        #print "Plate Strain add..."
        assert eid not in self.oxx
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
        assert nodeID not in self.oxx[eid]
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
        msg = '---PLATE STRESS---\n'
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
                    msg += "eid=%-4s eType=%s nid=%-4s iLayer=%s oxx=%-4i oyy=%-4i txy=%-4i ovm=%-4i\n" %(eid,eType,nid,iLayer,oxx,oyy,txy,ovm)
                ###
            ###
        ###
        return msg



class rodStressObject(scalarObject):
    """
                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
    ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
      ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
    """
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.eType = 'CROD'
        self.axial = {}
        self.SMa = {}
        self.torsion = {}
        self.SMt = {}

    def addNewEid(self,eid,axial,SMa,torsion,SMt):
        #print "Rod Stress add..."
        self.axial[eid] = axial
        self.SMa[eid] = SMa
        self.torsion[eid] = torsion
        self.SMt[eid] = SMt

    def __repr__(self):
        msg = '---ROD STRESSES---\n'
        for eid in sorted(self.axial):
            axial   = self.axial[eid]
            torsion = self.torsion[eid]
            SMa = self.SMa[eid]
            SMt = self.SMt[eid]
            msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

class rodStrainObject(scalarObject):
    """
                                     S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )
    ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
      ID.        STRAIN       MARGIN        STRAIN      MARGIN
    """
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.eType = 'CROD'
        self.axial = {}
        #self.SMa = {}
        self.torsion = {}
        #self.SMt = {}

    def addNewEid(self,axial,SMa,torsion,SMt):
        #print "Rod Strain add..."
        self.eType[eid] = self.eType
        self.axial[eid] = axial
        #self.SMa[eid] = SMa
        self.torsion[eid] = torsion
        #self.SMt[eid] = SMt

    def __repr__(self):
        msg = '---ROD STRAINS---\n'
        for eid in sorted(self.axial):
            axial   = self.axial[eid]
            torsion = self.torsion[eid]
            SMa = self.SMa[eid]
            SMt = self.SMt[eid]
            msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,eType,axial,torsion)
        return msg

class barStressObject(scalarObject):
    """
                               S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )
    ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
      ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C
    """
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.eType = {}
        self.s1 = {}
        self.s2 = {}
        self.s3 = {}
        self.s4 = {}
        self.axial = {}
        self.smax = {}
        self.smin = {}
        #self.MS_tension = {}
        #self.MS_compression = {}

    def addNewEid(self,eType,eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                                 s1b,s2b,s3b,s4b,      smaxb,sminb,MSc):
        #print "Bar Stress add..."
        self.eType[eid] = eType
        self.s1[eid] = [s1a,s1b]
        self.s2[eid] = [s2a,s2b]
        self.s3[eid] = [s3a,s3b]
        self.s4[eid] = [s4a,s4b]
        self.axial[eid] = axial
        self.smax[eid] = [smaxa,smaxb]
        self.smin[eid] = [smina,sminb]
        #self.MS_tension[eid]     = MSt
        #self.MS_compression[eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def __repr__(self):
        msg = '---BAR STRESS---\n'

        for eid,S1s in sorted(self.s1.items()):
            eType = self.eType[eid]
            axial = self.axial[eid]
            #MSt = self.MSt[eid]
            #MSc = self.MSc[eid]

            s1 = self.s1[eid]
            s2 = self.s2[eid]
            s3 = self.s3[eid]
            s4 = self.s4[eid]
            smax  = self.smax[eid]
            smin  = self.smin[eid]
            msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%4i smax=%-4i smax=%-4i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial, smax[0],smin[0])
            msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s         smax=%-4i smax=%-4i\n" %(' '*4,    s1[1],s2[1],s3[1],s4[1],'    ',smax[1],smin[1])
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

    """
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.eType = {}
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

    def addNewEid(self,eType,cid,eid,nodeID,oxx,oyy,ozz,txy,tyz,txz,aCos,bCos,cCos,pressure,ovm):
        #print "Solid Stress add..."
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

    def __repr__(self):
        msg = '---SOLID STRESS---\n'
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
                msg += "eid=%-4s eType=%-6s nid=%-2i oxx=%-5i oyy=%-5i ozz=%-5i txy=%-5i tyz=%-5i txz=%-5i ovm=%-5i\n" %(eid,eType,nid,oxx,oyy,ozz,txy,tyz,txz,ovm)
            ###
        ###
        return msg

