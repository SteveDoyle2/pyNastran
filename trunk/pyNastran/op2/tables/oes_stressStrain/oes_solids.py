import sys
from numpy.linalg import eigh

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

    def addF06Data(self,data,transient):
        if transient is None:
            if not hasattr(self,'data'):
                self.data = []
            self.data += data
        else:
            if not hasattr(self,'data'):
                self.data = {}
            if dt not in self.data:
                self.data[dt] = []
            for line in data:
                self.data[dt] += data
            ###
        ###

    def processF06Data(self):
        """
                          S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )
                       CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   
        ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES 
                2           0GRID CS  6 GP
                       CENTER  X  -1.829319E+03  XY   7.883865E+01   A   1.033182E+04  LX-0.12 0.71 0.70  -2.115135E+03    1.232595E+04
                               Y  -1.825509E+03  YZ  -1.415218E+03   B  -2.080181E+03  LY-0.12 0.69-0.71
                               Z   1.000023E+04  ZX  -1.415218E+03   C  -1.906232E+03  LZ 0.99 0.16 0.00
        """
        eMap = {'CTETRA':5,'CPENTA':7,'CHEXA':9,'HEXA':9,'PENTA':7,'TETRA':5,}   # +1 for the centroid
        if self.dt is None:
            pack = []
            i=0
            n=0
            while n<len(self.data):
                line = self.data[n]
                #print n,line

                eType = line[0]
                eid   = int(line[1])
                #print "eType = ",eType
                nNodes = eMap[eType]
                #nodeID = None
                self.eType[eid] = eType
                self.oxx[eid] = {}
                self.oyy[eid] = {}
                self.ozz[eid] = {}
                self.txy[eid] = {}
                self.tyz[eid] = {}
                self.txz[eid] = {}
                self.o1[eid]  = {}
                self.o2[eid]  = {}
                self.o3[eid]  = {}
                self.ovmShear[eid] = {}
                n+=1
                for j in range(nNodes):
                    #print self.data[n]
                    (blank,nodeID,x,oxx,xy,txy,a,o1,lx,d1,d2,d3,pressure,ovmShear) = self.data[n]
                    (blank,blank, y,oyy,yz,tyz,b,o2,ly,d1,d2,d3,blank,blank) = self.data[n+1]
                    (blank,blank, z,ozz,zx,txz,c,o3,lz,d1,d2,d3,blank,blank) = self.data[n+2]
                    if nodeID.strip()=='CENTER': nodeID='C'
                    self.oxx[eid][nodeID] = float(oxx)
                    self.oyy[eid][nodeID] = float(oyy)
                    self.ozz[eid][nodeID] = float(ozz)
                    self.txy[eid][nodeID] = float(txy)
                    self.tyz[eid][nodeID] = float(tyz)
                    self.txz[eid][nodeID] = float(txz)
                    self.o1[eid][nodeID] = float(o1)
                    self.o2[eid][nodeID] = float(o2)
                    self.o3[eid][nodeID] = float(o3)
                    self.ovmShear[eid][nodeID] = float(ovmShear)
                    n+=3
                ###
            ###
            return
        raise NotImplementedError()

        (dtName,dt) = transient
        self.dataCode['name'] = dtName
        if dt not in self.gridTypes:
            self.addNewTransient()

        for line in data:
            (nodeID,gridType,t1,t2,t3,r1,r2,r3) = line
            self.gridTypes[dt][nodeID]    = array([t1,t2,t3])
            self.translations[dt][nodeID] = array([t1,t2,t3])
            self.rotations[dt][nodeID]    = array([r1,r2,r3])
        ###
        del self.data

    def deleteTransient(self,dt):
        del self.oxx[dt]
        del self.oyy[dt]
        del self.ozz[dt]
        del self.txy[dt]
        del self.tyz[dt]
        del self.txz[dt]
        del self.o1[dt]
        del self.o2[dt]
        del self.o3[dt]

    def getTransients(self):
        k = self.oxx.keys()
        k.sort()
        return k

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

    def pressure(self,o1,o2,o3):
        """
        returns the hydrostatic pressure
        (o1+o2+o3)/-3.
        http://en.wikipedia.org/wiki/Stress_%28mechanics%29
        """
        return (o1+o2+o3)/-3.

    def directionalVectors(self,oxx,oyy,ozz,txy,tyz,txz):
        A = [[oxx,txy,txz],
             [txy,oyy,tyz],
             [txz,tyz,ozz]]
        (Lambda,v) = eigh(A) # a hermitian matrix is a symmetric-real matrix
        return v
    
    def writeF06(self,header,pageStamp,pageNum=1):
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)
        msg = []

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        tetraMsg = ['                   S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n',
                    '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN\n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        pentaMsg = ['                    S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n',
                    '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN\n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        hexaMsg = ['                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n',
                   '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN\n',
                   '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        eTypes = self.eType.values()
        isTetra=False; isPenta=False; isHexa=False
        if 'CTETRA' in eTypes or 'TETRA' in eTypes:
            isTetra = True
        if 'CPENTA' in eTypes or 'PENTA' in eTypes:
            isPenta = True
        if 'CHEXA' in eTypes or 'HEXA' in eTypes:
            isHexa = True
        
        nNodes = {'CTETRA':4,'CPENTA':6,'CHEXA':8,'HEXA':8,'PENTA':6,'TETRA':4,}
        for eid,oxxNodes in sorted(self.oxx.items()):
            eType = self.eType[eid]

            k = oxxNodes.keys()
            #kc = k.index('C')
            #k.pop(kc)
            k.sort()
            k.pop(-1)
            #print k
            msgA  = '0  %8s           0GRID CS  %i GP\n' %(eid,nNodes[eType])
            for nid in ['C']+k:
                oxx = self.oxx[eid][nid]
                oyy = self.oyy[eid][nid]
                ozz = self.ozz[eid][nid]
                txy = self.txy[eid][nid]
                tyz = self.tyz[eid][nid]
                txz = self.txz[eid][nid]

                o1 = self.o1[eid][nid]
                o2 = self.o2[eid][nid]
                o3 = self.o3[eid][nid]
                ovm = self.ovmShear[eid][nid]
                p = (o1+o2+o3)/-3.

                if nid=='C': nid='CENTER'
                print "nid = |%r|" %(nid)
                A = [[oxx,txy,txz],
                     [txy,oyy,tyz],
                     [txz,tyz,ozz]]
                (Lambda,v) = eigh(A) # a hermitian matrix is a symmetric-real matrix

                msgA += '0              %8s  X  %13.6E  XY  %13.6E   A  %13.6E  LX%5.2f%5.2f%5.2f  %13.6E   %13.6E\n' %(nid,oxx,txy,o1,v[0,0],v[0,1],v[0,2],p,ovm)
                msgA += '               %8s  Y  %13.6E  YZ  %13.6E   B  %13.6E  LY%5.2f%5.2f%5.2f\n'                  %('', oyy,tyz,o2,v[1,0],v[1,1],v[1,2])
                msgA += '               %8s  Z  %13.6E  ZX  %13.6E   C  %13.6E  LZ%5.2f%5.2f%5.2f\n'                  %('', ozz,txz,o3,v[2,0],v[2,1],v[2,2])
            ###
            if eType=='CTETRA' or eType=='TETRA':
                tetraMsg.append(msgA)
            elif eType=='CPENTA' or eType=='PENTA':
                pentaMsg.append(msgA)
            elif eType=='CHEXA' or eType=='HEXA':
                hexaMsg.append(msgA)
            else:
                raise NotImplementedError('eType=|%r|' %(eType))
            ###                
        ###
        if isTetra:
            msg += header+tetraMsg
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        if isPenta:
            msg += header+pentaMsg
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        if isHexa:
            msg += header+hexaMsg
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        return (''.join(msg),pageNum-1)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        msg = []

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        tetraMsg = ['                   S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n',
                    '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN\n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        pentaMsg = ['                    S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n',
                    '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN\n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        hexaMsg = ['                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n',
                   '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN\n',
                   '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        eTypes = self.eType.values()
        isTetra=False; isPenta=False; isHexa=False
        if 'CTETRA' in eTypes or 'TETRA' in eTypes:
            isTetra = True
        if 'CPENTA' in eTypes or 'PENTA' in eTypes:
            isPenta = True
        if 'CHEXA' in eTypes or 'HEXA' in eTypes:
            isHexa = True
        
        nNodes = {'CTETRA':4,'CPENTA':6,'CHEXA':8,'HEXA':8,'PENTA':6,'TETRA':4,}
        
        for dt,oxxs in sorted(self.oxx.items()):            
            dtLine = '%14s = %12.5E\n'%(self.dataCode['name'],dt)
            header[2] = dtLine
            for eid,oxxNodes in sorted(oxxs.items()):
                eType = self.eType[eid]

                k = oxxNodes.keys()
                #kc = k.index('C')
                #k.pop(kc)
                k.sort()
                k.pop(-1)
                #print k
                msgA  = '0  %8s           0GRID CS  %i GP\n' %(eid,nNodes[eType])
                for nid in ['C']+k:
                    oxx = self.oxx[dt][eid][nid]
                    oyy = self.oyy[dt][eid][nid]
                    ozz = self.ozz[dt][eid][nid]
                    txy = self.txy[dt][eid][nid]
                    tyz = self.tyz[dt][eid][nid]
                    txz = self.txz[dt][eid][nid]

                    o1 = self.o1[dt][eid][nid]
                    o2 = self.o2[dt][eid][nid]
                    o3 = self.o3[dt][eid][nid]
                    ovm = self.ovmShear[dt][eid][nid]
                    p = (o1+o2+o3)/-3.

                    if nid=='C': nid='CENTER'
                    #print "nid = |%r|" %(nid)
                    A = [[oxx,txy,txz],
                         [txy,oyy,tyz],
                         [txz,tyz,ozz]]
                    (Lambda,v) = eigh(A) # a hermitian matrix is a symmetric-real matrix

                    msgA += '0              %8s  X  %13.6E  XY  %13.6E   A  %13.6E  LX%5.2f%5.2f%5.2f  %13.6E   %13.6E\n' %(nid,oxx,txy,o1,v[0,0],v[0,1],v[0,2],p,ovm)
                    msgA += '               %8s  Y  %13.6E  YZ  %13.6E   B  %13.6E  LY%5.2f%5.2f%5.2f\n'                  %('', oyy,tyz,o2,v[1,0],v[1,1],v[1,2])
                    msgA += '               %8s  Z  %13.6E  ZX  %13.6E   C  %13.6E  LZ%5.2f%5.2f%5.2f\n'                  %('', ozz,txz,o3,v[2,0],v[2,1],v[2,2])
                ###
                if eType=='CTETRA' or eType=='TETRA':
                    tetraMsg.append(msgA)
                elif eType=='CPENTA' or eType=='PENTA':
                    pentaMsg.append(msgA)
                elif eType=='CHEXA' or eType=='HEXA':
                    hexaMsg.append(msgA)
                else:
                    raise NotImplementedError('eType=|%r|' %(eType))
                ###                
            ###
            if isTetra:
                msg += header+tetraMsg
                msg.append(pageStamp+str(pageNum)+'\n')
                pageNum+=1
            if isPenta:
                msg += header+pentaMsg
                msg.append(pageStamp+str(pageNum)+'\n')
                pageNum+=1
            if isHexa:
                msg += header+hexaMsg
                msg.append(pageStamp+str(pageNum)+'\n')
                pageNum+=1
            ###
        ###
        return (''.join(msg),pageNum-1)

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

    def addF06Data(self,data,transient):
        if transient is None:
            if not hasattr(self,'data'):
                self.data = []
            self.data += data
        else:
            if not hasattr(self,'data'):
                self.data = {}
            if dt not in self.data:
                self.data[dt] = []
            for line in data:
                self.data[dt] += data
            ###
        ###

    def processF06Data(self):
        """
                          S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )
                       CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   
        ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES 
                2           0GRID CS  6 GP
                       CENTER  X  -1.829319E+03  XY   7.883865E+01   A   1.033182E+04  LX-0.12 0.71 0.70  -2.115135E+03    1.232595E+04
                               Y  -1.825509E+03  YZ  -1.415218E+03   B  -2.080181E+03  LY-0.12 0.69-0.71
                               Z   1.000023E+04  ZX  -1.415218E+03   C  -1.906232E+03  LZ 0.99 0.16 0.00
        """
        eMap = {'CTETRA':5,'CPENTA':7,'CHEXA':9,'HEXA':9,'PENTA':7,'TETRA':5,}   # +1 for the centroid
        if self.dt is None:
            pack = []
            i=0
            n=0
            while n<len(self.data):
                line = self.data[n]
                #print n,line

                eType = line[0]
                eid   = int(line[1])
                #print "eType = ",eType
                nNodes = eMap[eType]
                #nodeID = None
                self.eType[eid] = eType
                self.exx[eid] = {}
                self.eyy[eid] = {}
                self.ezz[eid] = {}
                self.exy[eid] = {}
                self.eyz[eid] = {}
                self.exz[eid] = {}
                self.e1[eid]  = {}
                self.e2[eid]  = {}
                self.e3[eid]  = {}
                self.evmShear[eid] = {}
                n+=1
                for j in range(nNodes):
                    #print self.data[n]
                    (blank,nodeID,x,exx,xy,exy,a,e1,lx,d1,d2,d3,pressure,evmShear) = self.data[n]
                    (blank,blank, y,eyy,yz,eyz,b,e2,ly,d1,d2,d3,blank,blank) = self.data[n+1]
                    (blank,blank, z,ezz,zx,exz,c,e3,lz,d1,d2,d3,blank,blank) = self.data[n+2]
                    if nodeID.strip()=='CENTER': nodeID='C'
                    self.exx[eid][nodeID] = float(exx)
                    self.eyy[eid][nodeID] = float(eyy)
                    self.ezz[eid][nodeID] = float(ezz)
                    self.exy[eid][nodeID] = float(exy)
                    self.eyz[eid][nodeID] = float(eyz)
                    self.exz[eid][nodeID] = float(exz)
                    self.e1[eid][nodeID] = float(e1)
                    self.e2[eid][nodeID] = float(e2)
                    self.e3[eid][nodeID] = float(e3)
                    self.evmShear[eid][nodeID] = float(evmShear)
                    n+=3
                ###
            ###
            return
        raise NotImplementedError()
        (dtName,dt) = transient
        self.dataCode['name'] = dtName
        if dt not in self.gridTypes:
            self.addNewTransient()

        for line in data:
            (nodeID,gridType,t1,t2,t3,r1,r2,r3) = line
            self.gridTypes[dt][nodeID]    = array([t1,t2,t3])
            self.translations[dt][nodeID] = array([t1,t2,t3])
            self.rotations[dt][nodeID]    = array([r1,r2,r3])
        ###
        del self.data

    def deleteTransient(self,dt):
        del self.exx[dt]
        del self.eyy[dt]
        del self.ezz[dt]
        del self.exy[dt]
        del self.eyz[dt]
        del self.exz[dt]
        del self.e1[dt]
        del self.e2[dt]
        del self.e3[dt]

    def getTransients(self):
        k = self.exx.keys()
        k.sort()
        return k

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

    def ovm(self,o11,o22,o33,o12,o13,o23):
        """http://en.wikipedia.org/wiki/Von_Mises_yield_criterion"""
        ovm = 0.5 * ( (o11-o22)**2+(o22-o33)**2+(o11-o33)**2+6*(o23**2+o31**2+o12**2))
        return ovm

    #def ovmPlane(self,o11,o22,o12):
        #"""http://en.wikipedia.org/wiki/Von_Mises_yield_criterion"""
        #ovm = (o11**2+o22**2-o1*o2+3*o12**2)**0.5
        #return ovm

    def octahedral(self,o11,o22,o33,o12,o13,o23):
        """http://en.wikipedia.org/wiki/Von_Mises_yield_criterion"""
        ovm = self.ovm(o11,o22,o33,o12,o13,o23)
        return ovm/3*2**0.5

    def pressure(self,e1,e2,e3):
        """
        returns the hydrostatic pressure
        (e1+e2+e3)/-3.
        http://en.wikipedia.org/wiki/Stress_%28mechanics%29
        """
        return (e1+e2+e3)/-3.

    def directionalVectors(self,exx,eyy,ezz,exy,eyz,exz):
        A = [[exx,exy,exz],
             [exy,eyy,eyz],
             [exz,eyz,ezz]]
        (Lambda,v) = eigh(A) # a hermitian matrix is a symmetric-real matrix
        return v

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)
        msg = []

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        tetraMsg = ['                     S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n',
                   '0                CORNER        ------CENTER AND CORNER POINT STRAINS---------       DIR.  COSINES       MEAN\n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        pentaMsg = ['                      S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n',
                   '0                CORNER        ------CENTER AND CORNER POINT STRAINS---------       DIR.  COSINES       MEAN\n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        hexaMsg  = ['                      S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n',
                   '0                CORNER        ------CENTER AND CORNER POINT STRAINS---------       DIR.  COSINES       MEAN\n',
                   '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        eTypes = self.eType.values()
        isTetra=False; isPenta=False; isHexa=False
        if 'CTETRA' in eTypes or 'TETRA' in eTypes:
            isTetra = True
        if 'CPENTA' in eTypes or 'PENTA' in eTypes:
            isPenta = True
        if 'CHEXA' in eTypes or 'HEXA' in eTypes:
            isHexa = True
        
        nNodes = {'CTETRA':4,'CPENTA':6,'CHEXA':8,'HEXA':8,'PENTA':6,'TETRA':4,}
        for eid,exxNodes in sorted(self.exx.items()):
            eType = self.eType[eid]

            k = exxNodes.keys()
            #kc = k.index('C')
            #k.pop(kc)
            k.sort()
            k.pop(-1)
            #print k
            msgA  = '0  %8s           0GRID CS  %i GP\n' %(eid,nNodes[eType])
            for nid in ['C']+k:
                exx = self.exx[eid][nid]
                eyy = self.eyy[eid][nid]
                ezz = self.ezz[eid][nid]
                exy = self.exy[eid][nid]
                eyz = self.eyz[eid][nid]
                exz = self.exz[eid][nid]

                e1 = self.e1[eid][nid]
                e2 = self.e2[eid][nid]
                e3 = self.e3[eid][nid]
                evm = self.evmShear[eid][nid]
                p = (e1+e2+e3)/-3.

                if nid=='C': nodeID='CENTER'
                A = [[exx,exy,exz],
                     [exy,eyy,eyz],
                     [exz,eyz,ezz]]
                (Lambda,v) = eigh(A) # a hermitian matrix is a symmetric-real matrix

                msgA += '0              %8s  X  %13.6E  XY  %13.6E   A  %13.6E  LX%5.2f%5.2f%5.2f  %13.6E   %13.6E\n' %(nid,exx,exy,e1,v[0,0],v[0,1],v[0,2],p,evm)
                msgA += '0              %8s  Y  %13.6E  YZ  %13.6E   B  %13.6E  LY%5.2f%5.2f%5.2f\n'                  %('', eyy,eyz,e2,v[1,0],v[1,1],v[1,2])
                msgA += '0              %8s  Z  %13.6E  ZX  %13.6E   C  %13.6E  LZ%5.2f%5.2f%5.2f\n'                  %('', ezz,exz,e3,v[2,0],v[2,1],v[2,2])
            ###
            if eType=='CTETRA' or eType=='TETRA':
                tetraMsg.append(msgA)
            elif eType=='CPENTA' or eType=='PENTA':
                pentaMsg.append(msgA)
            elif eType=='CHEXA' or eType=='HEXA':
                hexaMsg.append(msgA)
            else:
                raise NotImplementedError('eType=|%r|' %(eType))
            ###
        ###
        msg = []
        if isTetra:
            msg += header+tetraMsg
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        if isPenta:
            msg += header+pentaMsg
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        if isHexa:
            msg += header+hexaMsg
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        return (''.join(msg),pageNum-1)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        msg = []

        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        tetraMsg = ['                     S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n',
                   '0                CORNER        ------CENTER AND CORNER POINT STRAINS---------       DIR.  COSINES       MEAN\n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        pentaMsg = ['                      S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n',
                   '0                CORNER        ------CENTER AND CORNER POINT STRAINS---------       DIR.  COSINES       MEAN\n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        hexaMsg  = ['                      S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n',
                   '0                CORNER        ------CENTER AND CORNER POINT STRAINS---------       DIR.  COSINES       MEAN\n',
                   '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s\n' %(vonMises)]

        eTypes = self.eType.values()
        isTetra=False; isPenta=False; isHexa=False
        if 'CTETRA' in eTypes or 'TETRA' in eTypes:
            isTetra = True
        if 'CPENTA' in eTypes or 'PENTA' in eTypes:
            isPenta = True
        if 'CHEXA' in eTypes or 'HEXA' in eTypes:
            isHexa = True
        
        nNodes = {'CTETRA':4,'CPENTA':6,'CHEXA':8,'HEXA':8,'PENTA':6,'TETRA':4,}
        
        for dt,exxs in sorted(self.exx.items()):            
            dtLine = '%14s = %12.5E\n'%(self.dataCode['name'],dt)
            header[2] = dtLine
            for eid,oxxNodes in sorted(exxs.items()):
                eType = self.eType[eid]

                k = exxNodes.keys()
                #kc = k.index('C')
                #k.pop(kc)
                k.sort()
                k.pop(-1)
                #print k
                msgA  = '0  %8s           0GRID CS  %i GP\n' %(eid,nNodes[eType])
                for nid in ['C']+k:
                    exx = self.exx[dt][eid][nid]
                    eyy = self.eyy[dt][eid][nid]
                    ezz = self.ezz[dt][eid][nid]
                    exy = self.exy[dt][eid][nid]
                    eyz = self.eyz[dt][eid][nid]
                    exz = self.exz[dt][eid][nid]

                    e1 = self.e1[dt][eid][nid]
                    e2 = self.e2[dt][eid][nid]
                    e3 = self.e3[dt][eid][nid]
                    evm = self.evmShear[dt][eid][nid]
                    p = (o1+o2+o3)/-3.

                    if nid=='C': nid='CENTER'
                    #print "nid = |%r|" %(nid)
                    A = [[exx,exy,exz],
                         [exy,eyy,eyz],
                         [exz,eyz,ezz]]
                    (Lambda,v) = eigh(A) # a hermitian matrix is a symmetric-real matrix

                    msgA += '0              %8s  X  %13.6E  XY  %13.6E   A  %13.6E  LX%5.2f%5.2f%5.2f  %13.6E   %13.6E\n' %(nid,exx,exy,o1,v[0,0],v[0,1],v[0,2],p,evm)
                    msgA += '               %8s  Y  %13.6E  YZ  %13.6E   B  %13.6E  LY%5.2f%5.2f%5.2f\n'                  %('', eyy,eyz,o2,v[1,0],v[1,1],v[1,2])
                    msgA += '               %8s  Z  %13.6E  ZX  %13.6E   C  %13.6E  LZ%5.2f%5.2f%5.2f\n'                  %('', ezz,exz,o3,v[2,0],v[2,1],v[2,2])
                ###
                if eType=='CTETRA' or eType=='TETRA':
                    tetraMsg.append(msgA)
                elif eType=='CPENTA' or eType=='PENTA':
                    pentaMsg.append(msgA)
                elif eType=='CHEXA' or eType=='HEXA':
                    hexaMsg.append(msgA)
                else:
                    raise NotImplementedError('eType=|%r|' %(eType))
                ###                
            ###
            if isTetra:
                msg += header+tetraMsg
                msg.append(pageStamp+str(pageNum)+'\n')
                pageNum+=1
            if isPenta:
                msg += header+pentaMsg
                msg.append(pageStamp+str(pageNum)+'\n')
                pageNum+=1
            if isHexa:
                msg += header+hexaMsg
                msg.append(pageStamp+str(pageNum)+'\n')
                pageNum+=1
            ###
        ###
        return (''.join(msg),pageNum-1)

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
