import sys
from oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

class celasStressObject(stressObject):
    """
                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )
        TIME         STRESS              TIME         STRESS              TIME         STRESS              TIME         STRESS
    0.0            0.0               5.000000E-02   0.0               1.000000E-01   0.0               1.500000E-01   0.0
    2.000000E-01   0.0               2.500000E-01   0.0               3.000000E-01   0.0               3.500000E-01   0.0
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = self.dataCode['ElementName']
        
        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.stress = {}
        if self.code in [[1,0,0]]:
            self.getLength = self.getLength_format1_sort0
            #self.isImaginary = False
            if dt is not None:
                self.addNewTransient = self.addNewTransient_format1_sort0
                self.addNewEid       = self.addNewEidTransient_format1_sort0
                self.isTransient = True
            else:
                self.addNewEid = self.addNewEid_format1_sort0
            ###
        elif 0: #self.code==[2,1,0]:
            self.getLength       = self.getLength_format1_sort0
            self.addNewTransient = self.addNewTransient_format2_sort1
            self.addNewEid       = self.addNewEidTransient_format2_sort1
            self.isImaginary = True
        else:
            raise InvalidCodeError('rodStress - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt is not None:
            #self.isTransient = True
            self.dt = self.nonlinearFactor
            self.addNewTransient()
        ###

    def getLength_format1_sort0(self):
        return (8,'if')

    def addNewTransient_format1_sort0(self):
        """
        initializes the transient variables
        """
        self.stress[self.dt] = {}

    def addNewTransient_format2_sort1(self):
        """
        initializes the transient variables
        """
        print self.dataCode
        self.axial[self.dt]     = {}
        self.torsion[self.dt]   = {}

    def addNewEid_format1_sort0(self,out):
        #print "Rod Stress add..."
        (eid,stress) = out
        eid = (eid-self.deviceCode)/10
        #assert isinstance(eid,int)
        self.stress[eid] = stress

    def addNewEid_format2_sort1(self,out):
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        assert eid >= 0
        self.axial[eid]      = [axialReal,axialImag]
        self.torsion[eid]    = [torsionReal,torsionImag]

    def addNewEidTransient_format1_sort0(self,out):
        (eid,stress) = out
        eid = (eid-self.deviceCode)/10
        dt = self.dt
        #assert isinstance(eid,int)
        #assert eid >= 0
        self.stress[dt][eid]      = stress

    def addNewEidTransient_format2_sort1(self,out):
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        dt = self.dt
        assert eid >= 0
        self.axial[dt][eid]      = [axialReal,axialImag]
        self.torsion[dt][eid]    = [torsionReal,torsionImag]

    def __reprTransient_format1_sort0__(self):
        msg = '---CELASx STRESSES---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['stress']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,stress in sorted(self.stress.items()):
            msg += 'dt = %g\n' %(dt)
            for eid,istress in sorted(stress.items()):
                msg += '%-6g %6s ' %(eid,self.eType)
                if abs(istress)<1e-6:
                    msg += '%10s ' %('0')
                else:
                    msg += '%10g ' %(val)
                ###
                msg += '\n'
            ###
        return msg

    def __reprTransient_format2_sort1__(self):
        msg = '---COMPLEX CELASx STRESSES---\n'
        msg += '%-10s %10s ' %('EID','eType')
        headers = ['axialReal','axialImag','torsionReal','torsionImag']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,axial in sorted(self.axial.items()):
            msg += 'dt = %g\n' %(dt)
            for eid in sorted(axial):
                axial   = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                msg += '%-6g %6s ' %(eid,self.eType)
                vals = axial + torsion # concatination
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10i ' %(val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
            ###
        return msg

    def __repr__(self):
        if   self.isTransient and self.code==[1,0,0]:
            return self.__reprTransient_format1_sort0__()
        elif self.code==[2,1,0]:
            return self.__reprTransient_format2_sort1__()

        msg = '---CELASx STRESSES---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['stress']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'
        print "self.code = ",self.code
        for eid,istress in sorted(self.stress.items()):

            msg += '%-6i %6s ' %(eid,self.eType)
            if abs(istress)<1e-6:
                msg += '%10s ' %('0')
            else:
                msg += '%10i ' %(val)
            ###
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

class celasStrainObject(strainObject):
    """
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CELASx'

        self.code = [self.formatCode,self.sortCode,self.sCode]
        
        self.isTransient = False
        if self.code == [1,0,10]:
            self.getLength = self.getLength_format1_sort0
            self.strain = {}
            #self.isImaginary = False
            if dt is not None:
                self.addNewTransient = self.addNewTransient_format1_sort0
                self.addNewEid       = self.addNewEidTransient_format1_sort0
                self.isTransient = True
            ###
            else:
                self.addNewEid = self.addNewEid_format1_sort0
            ###

        elif self.code==[2,1,10]:
            self.getLength = self.getLength_format1_sort0
            self.axial      = {}
            self.torsion    = {}
            self.addNewTransient = self.addNewTransient_format2_sort1
            self.addNewEid       = self.addNewEidTransient_format2_sort1
            #self.isImaginary = True
        else:
            raise InvalidCodeError('celasStrain - get the format/sort/stressCode=%s' %(self.code))
        ###
        self.dt = self.nonlinearFactor
        if dt is not None:
            self.addNewTransient()
        ###

    def getLength_format1_sort0(self):
        return (8,'if')

    def addNewTransient_format1_sort0(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.strain[self.dt] = {}

    def addNewTransient_format2_sort1(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        print self.dataCode
        self.axial[self.dt]     = {}
        self.torsion[self.dt]   = {}

    def addNewEid_format1_sort0(self,out):
        (eid,strain) = out
        eid = (eid-self.deviceCode)/10
        #print "Rod Strain add..."
        assert eid >= 0
        #self.eType = self.eType
        self.strainl[eid] = strain

    def addNewEid_format2_sort1(self,out):
        assert eid >= 0
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        self.axial[eid]      = [axialReal,axialImag]
        self.torsion[eid]    = [torsionReal,torsionImag]

    def addNewEidTransient_format1_sort0(self,out):
        #print "Rod Strain add..."
        print out
        (eid,strain) = out
        eid = (eid-self.deviceCode)/10
        assert eid >= 0
        dt = self.dt

        self.eType[eid] = self.elementType
        self.strain[dt][eid] = strain

    def addNewEidTransient_format2_sort1(self,out):
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        assert eid >= 0
        dt = self.dt
        self.axial[dt][eid]      = [axialReal,axialImag]
        self.torsion[dt][eid]    = [torsionReal,torsionImag]

    def __reprTransient_format2_sort1__(self):
        msg = '---COMPLEX CELASx STRAINS---\n'
        msg += '%-10s %10s ' %('EID','eType')
        headers = ['axialReal','axialImag','torsionReal','torsionImag']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,axial in sorted(self.axial.items()):
            msg += 'dt = %s\n' %(dt)
            for eid in sorted(axial):
                axial   = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                msg += '%-6i %6s ' %(eid,self.eType)
                vals = axial + torsion # concatination
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10g ' %(val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
            ###
        return msg

    def __repr__(self):
        if   self.isTransient and self.code==[1,0,10]:
            return self.__reprTransient_format1_sort0__()
        elif self.code==[2,1,10]:
            return self.__reprTransient_format2_sort1__()

        msg = '---CELASx STRAINS---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['strain']
        for header in headers:
            msg += '%8s ' %(header)
        msg += '\n'

        for eid in sorted(self.axial):
            strain = self.strain[eid]
            msg += '%-6i %6s ' %(eid,self.eType)

            if abs(strain)<1e-7:
                msg += '%8s ' %('0')
            else:
                msg += '%1.3g ' %(val)
            ###
            msg += '\n'
        return msg
