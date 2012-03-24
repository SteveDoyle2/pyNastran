import sys
from oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

class nonlinearQuadObject(stressObject):
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = 'QUAD4FD'
        
        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.Type = {}
        self.IDs = {}
        self.oxx = {}
        self.oyy = {}
        self.txy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}

    def addNewTransient(self):
        dt = self.dt
        self.oxx[dt] = {}
        self.oyy[dt] = {}
        self.txy[dt] = {}
        self.angle[dt]  = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}

    def addNewEid(self,data):
        dt = self.dt
        if dt not in self.oxx:
            self.addNewTransient()
        (eid,Type,oxx,oyy,txy,angle,majorP,minorP) = data
        self.Type[eid] = Type
        self.oxx[dt] = {eid:[oxx]}
        self.oyy[dt] = {eid:[oyy]}
        self.txy[dt] = {eid:[txy]}
        self.angle[dt]  = {eid:[angle]}
        self.majorP[dt] = {eid:[majorP]}
        self.minorP[dt] = {eid:[minorP]}

    def add(self,eid,data):
        dt = self.dt
        (ID,oxx,oyy,txy,angle,majorP,minorP) = data
        self.oxx[dt][eid].append(oxx)
        self.oyy[dt][eid].append(oyy)
        self.txy[dt][eid].append(txy)
        self.angle[dt][eid].append(angle)
        self.majorP[dt][eid].append(majorP)
        self.minorP[dt][eid].append(minorP)

    def writeF06(self,header,pageStamp,pageNum=1):
        print "hyper quad...."
        msg = ['           S T R E S S E S   I N   H Y P E R E L A S T I C   Q U A D R I L A T E R A L   E L E M E N T S  ( QUAD4FD )',
               '  ELEMENT     GRID/    POINT       ---------CAUCHY STRESSES--------             PRINCIPAL STRESSES (ZERO SHEAR)',
               '     ID       GAUSS      ID      NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR',]
               #0       1     GAUS         1   7.318995E+00   6.367099E-01  -6.551054E+00   -31.4888    1.133173E+01   -3.376026E+00
               #                           2   1.097933E+01   4.149028E+00   6.278160E+00    30.7275    1.471111E+01    4.172537E-01
        
        for dt,Oxxs in sorted(self.oxx.items()):
            #header[-1] = '     LOAD STEP = %12.5E' %(dt)
            msg += header
            for eid,oxxs in sorted(Oxxs.items()):
                gauss = self.Type[eid]
                oxx = self.oxx[dt][eid]
                oyy = self.oyy[dt][eid]
                txy = self.txy[dt][eid]
                angle = self.angle[dt][eid]
                majorP = self.majorP[dt][eid]
                minorP = self.minorP[dt][eid]
                
                for i in range(4): # 1,2,3,4
                    if i==0:
                        msg.append('0%8i %8s  %8i  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6' %(eid,gauss,i+1,oxx[i],oyy[i],txy[i],angle[i],majorP[i],minorP[i]))
                    else:
                        msg.append(' %8s %8s  %8i  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6' %('','',    i+1,oxx[i],oyy[i],txy[i],angle[i],majorP[i],minorP[i]))
                    ###
                ###
            ###
        ###                    
        return ('\n'.join(msg),pageNum)    

    def __repr__(self):
        return self.writeF06([],'PAGE ',1)[0]

class nonlinearRodObject(stressObject):
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CROD'
        
        self.code = [self.formatCode,self.sortCode,self.sCode]

        self.axialStress = {}
        self.equivStress = {}
        self.totalStrain = {}
        self.effectivePlasticCreepStrain = {}
        self.effectiveCreepStrain  = {}
        self.linearTorsionalStress = {}
    
    def addNewTransient(self):
        dt = self.dt
        self.axialStress[dt] = {}
        self.equivStress[dt] = {}
        self.totalStrain[dt] = {}
        self.effectivePlasticCreepStrain[dt] = {}
        self.effectiveCreepStrain[dt]  = {}
        self.linearTorsionalStress[dt] = {}

    def add(self,data):
        dt = self.dt
        if dt not in self.axialStress:
            self.addNewTransient()
        eid = data[0]
        self.axialStress[dt][eid] = data[1]
        self.equivStress[dt][eid] = data[2]
        self.totalStrain[dt][eid] = data[3]
        self.effectivePlasticCreepStrain[dt][eid] = data[4]
        self.effectiveCreepStrain[dt][eid]  = data[5]
        self.linearTorsionalStress[dt][eid] = data[6]
        #print data

    def writeF06(self,header,pageStamp,pageNum=1):
        """
        ELEMENT-ID =     102
                                 N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
          TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL
                                               STRESS                             PLASTIC/NLELAST          STRAIN              STRESS
        2.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
        3.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
        """
        msg = []
        msgStart = ['                         N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )',
               ' ',
               '    TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL',
               '                                         STRESS                             PLASTIC/NLELAST          STRAIN              STRESS']
        msgE = {}
        msgT = {}
        for dt,axials in sorted(self.axialStress.items()):
            for eid,axial in sorted(axials.items()):
                eqs  = self.equivStress[dt][eid]
                ts   = self.totalStrain[dt][eid]
                epcs = self.effectivePlasticCreepStrain[dt][eid]
                ecs  = self.effectiveCreepStrain[dt][eid]
                lts  = self.linearTorsionalStress[dt][eid]
                msgE[eid] = '      ELEMENT-ID = %8i' %(eid)
                if eid not in msgT:
                    msgT[eid] = []
                msgT[eid].append('  %9.3E       %13.6E       %13.6E       %13.6E       %13.6E       %13.6E       %13.6E'%(dt,axial,eqs,ts,epcs,ecs,lts))
            ###
        ###
        for eid,e in sorted(msgE.items()):
            msg += header+[e]+msgStart+msgT[eid]
            msg.append(pageStamp+str(pageNum))

        return ('\n'.join(msg),pageNum)    

    def __repr__(self):
        return self.writeF06([],'PAGE ',1)[0]

