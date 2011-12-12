import sys
from op2_Objects import stressObject,strainObject #,array
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
        if self.code == [1,0,0]:
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
            raise InvalidCodeError('compositePlateStress - get the format/sort/stressCode=%s' %(self.code))
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
        #assert eid not in self.o11[dt]
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

    def __reprTransient__(self):
        msg = '---COMPOSITE PLATE STRESS---\n'
        msg += '%-6s %8s %8s ' %('EID','eType','iLayer')
        headers = ['o11','o22','t12','t1z','t2z','ovm']
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
                    ovm   = self.ovm[dt][eid][iLayer]

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
        if self.code == [1,0,14]:
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
            raise InvalidCodeError('compositePlateStrain - get the format/sort/stressCode=%s' %(self.code))
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
        #assert eid not in self.e11
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

    def __reprTransient__(self):
        msg = '---COMPOSITE PLATE STAIN---\n'
        headers = ['e11','e22','e12','e1z','e2z','evm']
        msg += '%-6s %8s %8s ' %('EID','eType','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,E11s in sorted(self.e11.items()):
            msg += 'dt = %s\n' %(dt)
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
                    evm   = self.evm[dt][eid][iLayer]

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
