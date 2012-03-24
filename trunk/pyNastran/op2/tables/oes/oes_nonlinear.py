import sys
from oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

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

