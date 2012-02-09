import sys
from oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

class solidStressObject(stressObject):
    """
    # sCode=0
            N O N L I N E A R   S T R E S S E S   I N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )
    ELEMENT GRID/   POINT                         STRESSES/ TOTAL STRAINS                          EQUIVALENT EFF. STRAIN  EFF. CREEP
       ID   GAUSS     ID       X           Y           Z           XY          YZ          ZX        STRESS   PLAS/NLELAS   STRAIN
        3    GRID   CENTER  6.6667E+02  2.6667E+03  2.6667E+03 -1.3333E+03  2.6667E+03 -1.3333E+03  6.0000E+03  1.5000E-04   0.0
                            1.6667E-05  6.6667E-05  6.6667E-05 -6.6667E-05  1.3333E-04 -6.6667E-05

    # sCode=1
                          S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )
                   CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   
    ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES 
        1               0GRID CS  8 GP
                   CENTER  X   4.499200E+02  XY  -5.544791E+02   A   1.000000E+04  LX 0.00 0.69-0.72  -3.619779E+03    9.618462E+03
                           Y   4.094179E+02  YZ   5.456968E-12   B  -1.251798E+02  LY 0.00 0.72 0.69
                           Z   1.000000E+04  ZX  -4.547474E-13   C   9.845177E+02  LZ 1.00 0.00 0.00
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)

        self.eType = {}

        self.code = [self.formatCode,self.sortCode,self.sCode]

        self.cid = {} # gridGauss
        self.oxx = {}
        self.oyy = {}
        self.ozz = {}
        self.txy = {}
        self.tyz = {}
        self.txz = {}
        self.o1 = {}
        self.o2 = {}
        self.o3 = {}
        self.ovmShear = {}
        if self.code == [1,0,0]:  # not done...
            pass
            #EQUIVALENT EFF. STRAIN  EFF. CREEP
            #  STRESS   PLAS/NLELAS   STRAIN
            #self.eqStress       = {}
            #self.effStrain      = {}
            #self.effCreepStrain = {}
        elif self.code == [1,0,1]:
            pass
            #self.aCos = {}
            #self.bCos = {}
            #self.cCos = {}
            #self.pressure = {}
        else:
            raise InvalidCodeError('solidStress - get the format/sort/stressCode=%s' %(self.code))

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
        self.o1[self.dt] = {}
        self.o2[self.dt] = {}
        self.o3[self.dt] = {}
        
        #self.aCos[self.dt] = {}
        #self.bCos[self.dt] = {}
        #self.cCos[self.dt] = {}
        #self.pressure[self.dt] = {}
        self.ovmShear[self.dt]  = {}

    def addNewEid(self,eType,cid,eid,nodeID,oxx,oyy,ozz,txy,tyz,txz,o1,o2,o3,aCos,bCos,cCos,pressure,ovm):
        #print "Solid Stress add..."
        #assert eid not in self.oxx
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
        self.o1[eid]   = {nodeID: o1}
        self.o2[eid]   = {nodeID: o2}
        self.o3[eid]   = {nodeID: o3}
        #self.aCos[eid] = {nodeID: aCos}
        #self.bCos[eid] = {nodeID: bCos}
        #self.cCos[eid] = {nodeID: cCos}
        #self.pressure[eid] = {nodeID: pressure}
        self.ovmShear[eid]  = {nodeID: ovm}
        msg = "*eid=%s nodeID=%s vm=%g" %(eid,nodeID,ovm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,cid,eid,nodeID,oxx,oyy,ozz,txy,tyz,txz,o1,o2,o3,aCos,bCos,cCos,pressure,ovm):
        #print "Solid Stress add transient..."
        assert cid >= 0
        assert eid >= 0
        dt = self.dt
        assert eid not in self.oxx[dt]
        self.eType[eid] = eType
        self.cid[eid]   = cid
        self.oxx[dt][eid]  = {nodeID: oxx}
        self.oyy[dt][eid]  = {nodeID: oyy}
        self.ozz[dt][eid]  = {nodeID: ozz}
        self.txy[dt][eid]  = {nodeID: txy}
        self.tyz[dt][eid]  = {nodeID: tyz}
        self.txz[dt][eid]  = {nodeID: txz}
        
        self.o1[dt][eid]  = {nodeID: o1}
        self.o2[dt][eid]  = {nodeID: o2}
        self.o3[dt][eid]  = {nodeID: o3}

        #self.aCos[dt][eid] = {nodeID: aCos}
        #self.bCos[dt][eid] = {nodeID: bCos}
        #self.cCos[dt][eid] = {nodeID: cCos}
        #self.pressure[dt][eid] = {nodeID: pressure}
        self.ovmShear[dt][eid]  = {nodeID: ovm}
        msg = "*eid=%s nodeID=%s vm=%g" %(eid,nodeID,ovm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,oxx,oyy,ozz,txy,tyz,txz,o1,o2,o3,aCos,bCos,cCos,pressure,ovm):
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

        self.o1[eid][nodeID] = o1
        self.o2[eid][nodeID] = o2
        self.o3[eid][nodeID] = o3

        #self.aCos[eid][nodeID] = aCos
        #self.bCos[eid][nodeID] = bCos
        #self.cCos[eid][nodeID] = cCos
        #self.pressure[eid][nodeID] = pressure
        self.ovmShear[eid][nodeID] = ovm
        if nodeID==0: raise Exception(msg)

    def addTransient(self,eid,nodeID,oxx,oyy,ozz,txy,tyz,txz,o1,o2,o3,aCos,bCos,cCos,pressure,ovm):
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

        self.o1[dt][eid][nodeID] = o1
        self.o2[dt][eid][nodeID] = o2
        self.o3[dt][eid][nodeID] = o3

        #self.aCos[dt][eid][nodeID] = aCos
        #self.bCos[dt][eid][nodeID] = bCos
        #self.cCos[dt][eid][nodeID] = cCos
        #self.pressure[dt][eid][nodeID] = pressure
        self.ovmShear[dt][eid][nodeID] = ovm
        if nodeID==0: raise Exception(msg)

    def getHeaders(self):
        headers = ['oxx','oyy','ozz','txy','tyz','txz']
        if self.isVonMises():
            headers.append('oVonMises')
        else:
            headers.append('oMaxShear')
        return headers

    def __reprTransient__(self):
        msg = '---SOLID STRESS---\n'
        msg += '%-6s %6s %8s ' %('EID','eType','nodeID')
        headers = self.getHeaders()
        for header in headers:
            msg += '%9s ' %(header)
        msg += '\n'

        for dt,oxxs in sorted(self.oxx.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid,oxxNodes in sorted(oxxs.items()):
                eType = self.eType[eid]
                for nid in sorted(oxxNodes):
                    oxx = self.oxx[dt][eid][nid]
                    oyy = self.oyy[dt][eid][nid]
                    ozz = self.ozz[dt][eid][nid]
                    txy = self.txy[dt][eid][nid]
                    tyz = self.tyz[dt][eid][nid]
                    txz = self.txz[dt][eid][nid]

                    #o1 = self.o1[dt][eid][nid]
                    #o2 = self.o2[dt][eid][nid]
                    #o3 = self.o3[dt][eid][nid]

                    ovm = self.ovmShear[dt][eid][nid]
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
        headers = self.getHeaders()
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

                #o1 = self.o1[eid][nid]
                #o2 = self.o2[eid][nid]
                #o3 = self.o3[eid][nid]

                ovm = self.ovmShear[eid][nid]
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

class solidStrainObject(strainObject):
    """
    # code=[1,0,11]
                          S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )
                   CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   
    ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES 
            1           0GRID CS  8 GP
                   CENTER  X   4.499200E+02  XY  -5.544791E+02   A   1.000000E+04  LX 0.00 0.69-0.72  -3.619779E+03    9.618462E+03
                           Y   4.094179E+02  YZ   5.456968E-12   B  -1.251798E+02  LY 0.00 0.72 0.69
                           Z   1.000000E+04  ZX  -4.547474E-13   C   9.845177E+02  LZ 1.00 0.00 0.00

    # code=[1,0,10]
                       S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )
                   CORNER        ------CENTER AND CORNER POINT  STRAINS---------       DIR.  COSINES       MEAN         OCTAHEDRAL
    ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE         SHEAR   
            4           0GRID CS  4 GP
                   CENTER  X  -2.288232E-04  XY   1.240506E-04   A   9.631978E-04  LX-0.10-0.71-0.70  -1.601805E-04    5.692614E-04
                           Y  -2.289814E-04  YZ  -2.369997E-04   B  -2.909276E-04  LY-0.10 0.71-0.70
                           Z   9.383460E-04  ZX  -2.369997E-04   C  -1.917288E-04  LZ 0.99 0.00-0.15
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)
        self.eType = {}

        self.code = [self.formatCode,self.sortCode,self.sCode]

        self.cid = {}
        self.exx = {}
        self.eyy = {}
        self.ezz = {}
        self.exy = {}
        self.eyz = {}
        self.exz = {}
        self.e1 = {}
        self.e2 = {}
        self.e3 = {}
        #self.aCos = {}
        #self.bCos = {}
        #self.cCos = {}
        #self.pressure = {}
        self.evmShear = {}
        
        #if self.isVonMises():
        #if   self.code == [1,0,10]:
        #    self.isVonMises = False
        #elif self.code == [1,0,11]:
        #    self.isVonMises = True
        #else:
        #    raise InvalidCodeError('solidStrain - get the format/sort/stressCode=%s' %(self.code))
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
        self.exx[self.dt] = {}
        self.eyy[self.dt] = {}
        self.ezz[self.dt] = {}
        self.exy[self.dt] = {}
        self.eyz[self.dt] = {}
        self.exz[self.dt] = {}
        self.e1[self.dt] = {}
        self.e2[self.dt] = {}
        self.e3[self.dt] = {}
        
        #self.aCos[self.dt] = {}
        #self.bCos[self.dt] = {}
        #self.cCos[self.dt] = {}
        #self.pressure[self.dt] = {}
        self.evmShear[self.dt]      = {}

    def addNewEid(self,eType,cid,eid,nodeID,exx,eyy,ezz,exy,eyz,exz,e1,e2,e3,aCos,bCos,cCos,pressure,evm):
        #print "Solid Strain add..."
        assert eid not in self.exx
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
        self.e1[eid]   = {nodeID: e1}
        self.e2[eid]   = {nodeID: e2}
        self.e3[eid]   = {nodeID: e3}
        #self.aCos[eid] = {nodeID: aCos}
        #self.bCos[eid] = {nodeID: bCos}
        #self.cCos[eid] = {nodeID: cCos}
        #self.pressure[eid] = {nodeID: pressure}
        self.evmShear[eid]      = {nodeID: evm}
        msg = "*eid=%s nodeID=%s evmShear=%g" %(eid,nodeID,evm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,cid,eid,nodeID,exx,eyy,ezz,exy,eyz,exz,e1,e2,e3,aCos,bCos,cCos,pressure,evm):
        #print "Solid Strain add..."
        assert cid >= 0
        assert eid >= 0
        dt = self.dt
        assert eid not in self.exx[dt]
        self.eType[eid] = eType
        self.cid[eid]  = cid
        self.exx[dt][eid]  = {nodeID: exx}
        self.eyy[dt][eid]  = {nodeID: eyy}
        self.ezz[dt][eid]  = {nodeID: ezz}
        self.exy[dt][eid]  = {nodeID: exy}
        self.eyz[dt][eid]  = {nodeID: eyz}
        self.exz[dt][eid]  = {nodeID: exz}
        self.e1[dt][eid]   = {nodeID: e1}
        self.e2[dt][eid]   = {nodeID: e2}
        self.e3[dt][eid]   = {nodeID: e3}
        #self.aCos[dt][eid] = {nodeID: aCos}
        #self.bCos[dt][eid] = {nodeID: bCos}
        #self.cCos[dt][eid] = {nodeID: cCos}
        #self.pressure[dt][eid] = {nodeID: pressure}
        self.evmShear[dt][eid]      = {nodeID: evm}
        msg = "*eid=%s nodeID=%s evmShear=%g" %(eid,nodeID,evm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,exx,eyy,ezz,exy,eyz,exz,e1,e2,e3,aCos,bCos,cCos,pressure,evm):
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
        self.e1[eid][nodeID] = e1
        self.e2[eid][nodeID] = e2
        self.e3[eid][nodeID] = e3

        #self.aCos[eid][nodeID] = aCos
        #self.bCos[eid][nodeID] = bCos
        #self.cCos[eid][nodeID] = cCos
        #self.pressure[eid][nodeID] = pressure
        self.evmShear[eid][nodeID] = evm

        if nodeID==0: raise Exception(msg)

    def addTransient(self,eid,nodeID,exx,eyy,ezz,exy,eyz,exz,e1,e2,e3,aCos,bCos,cCos,pressure,evm):
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

        self.e1[dt][eid][nodeID] = e1
        self.e2[dt][eid][nodeID] = e2
        self.e3[dt][eid][nodeID] = e3

        #self.aCos[dt][eid][nodeID] = aCos
        #self.bCos[dt][eid][nodeID] = bCos
        #self.cCos[dt][eid][nodeID] = cCos
        #self.pressure[dt][eid][nodeID] = pressure
        self.evmShear[dt][eid][nodeID] = evm

        if nodeID==0: raise Exception(msg)

    def getHeaders(self):
        headers = ['exx','eyy','ezz','exy','eyz','exz']
        if self.isVonMises():
            headers.append('evm')
        else:
            headers.append('maxShear')
        return headers

    def __reprTransient__(self):
        msg = '---SOLID STRAIN---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s ' %('EID','eType','nodeID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'
        for dt,exxs in sorted(self.exx.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],self.dt)
            for eid,exxNodes in sorted(exxs.items()):
                eType = self.eType[eid]
                for nid in sorted(exxNodes):
                    exx = self.exx[dt][eid][nid]
                    eyy = self.eyy[dt][eid][nid]
                    ezz = self.ezz[dt][eid][nid]
                    exy = self.exy[dt][eid][nid]
                    eyz = self.eyz[dt][eid][nid]
                    exz = self.exz[dt][eid][nid]
                    evm = self.evmShear[dt][eid][nid]
                    msg += '%-6i %6s %8s ' %(eid,eType,nid)
                    vals = [exx,eyy,ezz,exy,eyz,exz,evm]
                    for val in vals:
                        if abs(val)<1e-6:
                            msg += '%10s ' %('0')
                        else:
                            msg += '%10.3e ' %(val)
                        ###
                    msg += '\n'
                    #msg += "eid=%-4s eType=%-6s nid=%-2i exx=%-5i eyy=%-5i ezz=%-5i exy=%-5i eyz=%-5i exz=%-5i evm=%-5i\n" %(eid,eType,nid,exx,eyy,ezz,exy,eyz,exz,evm)
                ###
            ###
        ###
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---SOLID STRAIN---\n'
        headers = self.getHeaders()
        msg += '%-6s %6s %8s ' %('EID','eType','nodeID')
        for header in headers:
            msg += '%10s ' %(header)
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
                evm = self.evmShear[eid][nid]
                msg += '%-6i %6s %8s ' %(eid,eType,nid)
                vals = [exx,eyy,ezz,exy,eyz,exz,evm]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10.3e ' %(val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%-6s nid=%-2i exx=%-5i eyy=%-5i ezz=%-5i exy=%-5i eyz=%-5i exz=%-5i evm=%-5i\n" %(eid,eType,nid,exx,eyy,ezz,exy,eyz,exz,evm)
            ###
        ###
        return msg

