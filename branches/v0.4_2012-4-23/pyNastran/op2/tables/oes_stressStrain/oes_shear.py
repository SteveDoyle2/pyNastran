import sys
from oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

class shearStressObject(stressObject):
    """
    # formatCode=1 sortCode=0 stressCode=0
                                   S T R E S S E S   I N   S H E A R   P A N E L S      ( C S H E A R )
    ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY
      ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN
        328        1.721350E+03   1.570314E+03   7.2E+01
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CSHEAR'
        
        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.maxShear = {}
        self.avgShear = {}
        self.margin   = {}
        
        if self.code in [[1,0,0],]:
            self.getLength = self.getLength_format1_sort0
            self.isImaginary = False
            if dt is not None:
                self.addNewTransient = self.addNewTransient_format1_sort0
                self.addNewEid       = self.addNewEidTransient_format1_sort0
            else:
                self.addNewEid = self.addNewEid_format1_sort0
            ###
        #elif self.code==[2,1,0]:
        #    self.getLength       = self.getLength_format1_sort0
        #    self.addNewTransient = self.addNewTransient_format2_sort1
        #    self.addNewEid       = self.addNewEidTransient_format2_sort1
        #    self.isImaginary = True
        else:
            raise InvalidCodeError('rodStress - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt is not None:
            self.isTransient = True
            self.dt = self.nonlinearFactor
            self.addNewTransient()
        ###

    def getLength_format1_sort0(self):
        return (16,'ifff')

    def deleteTransient(self,dt):
        del self.maxShear[dt]
        del self.avgShear[dt]
        del self.margin[dt]

    def getTransients(self):
        k = self.maxShear.keys()
        k.sort()
        return k

    def addNewTransient_format1_sort0(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        if self.dt not in self.maxShear:
            self.maxShear[self.dt] = {}
            self.avgShear[self.dt] = {}
            self.margin[self.dt]   = {}

    def addNewTransient_format2_sort1(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        #print self.dataCode
        raise Exeption('not supported')
        if self.dt not in self.axial:
            self.axial[self.dt]     = {}
            self.torsion[self.dt]   = {}

    def addNewEid_format1_sort0(self,out):
        #print "Rod Stress add..."
        (eid,maxShear,avgShear,margin) = out
        eid = (eid-self.deviceCode) // 10
        assert isinstance(eid,int)
        self.maxShear = {}
        self.avgShear = {}
        self.margin   = {}
        self.maxShear[eid] = maxShear
        self.avgShear[eid] = avgShear
        self.margin[eid]   = margin

    def addNewEid_format2_sort1(self,out):
        raise Exception('not supported')
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode) // 10
        assert eid >= 0
        self.axial[eid]      = [axialReal,axialImag]
        self.torsion[eid]    = [torsionReal,torsionImag]

    def addNewEidTransient_format1_sort0(self,out):
        (eid,maxShear,avgShear,margin) = out
        eid = (eid-self.deviceCode) // 10
        dt = self.dt
        assert isinstance(eid,int)
        assert eid >= 0
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid]   = margin

    def addNewEidTransient_format2_sort1(self,out):
        raise Exception('not supported')
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode) // 10
        dt = self.dt
        assert eid >= 0
        self.axial[dt][eid]      = [axialReal,axialImag]
        self.torsion[dt][eid]    = [torsionReal,torsionImag]

    def __reprTransient_format1_sort0__(self):
        msg = '---TRANSIENT CSHEAR STRESSES---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['maxShear','avgShear','Margin']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,maxShears in sorted(self.maxShear.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid in sorted(maxShears):
                maxShear = self.maxShear[dt][eid]
                avgShear = self.avgShear[dt][eid]
                margin   = self.margin[dt][eid]
                msg += '%-6i %6s ' %(eid,self.eType)
                vals = [maxShear,avgShear,margin]
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

    def __reprTransient_format2_sort1__(self):
        raise Exeption('not supported')
        msg = '---COMPLEX ROD STRESSES---\n'
        msg += '%-10s %10s ' %('EID','eType')
        headers = ['axialReal','axialImag','torsionReal','torsionImag']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,axial in sorted(self.axial.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid in sorted(axial):
                axial   = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                msg += '%-6i %6s ' %(eid,self.eType)
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
        if   self.isTransient and self.code in [[1,0,0],[1,0,1]]:
            return self.__reprTransient_format1_sort0__()
        elif self.code==[2,1,0]:
            return self.__reprTransient_format2_sort1__()
        #else:
        #    raise Exception('code=%s' %(self.code))
        msg = '---CSHEAR STRESSES---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['maxShear','avgShear','margin']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'
        #print "self.code = ",self.code
        for eid in sorted(self.maxShear):
            #print self.__dict__.keys()
            maxShear = self.maxShear[eid]
            avgShear = self.avgShear[eid]
            margin   = self.margin[eid]
            msg += '%-6i %6s ' %(eid,self.eType)
            vals = [maxShear,avgShear,margin]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %('0')
                else:
                    msg += '%10i ' %(val)
                ###
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

class shearStrainObject(strainObject):
    """
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CSHEAR'
        raise Exception('not supported...CSHEAR strain')
        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.maxShear = {}
        self.avgShear = {}
        self.margin   = {}

        self.isTransient = False
        if self.code == [1,0,10]:
            self.getLength = self.getLength_format1_sort0
            self.isImaginary = False
            if dt is not None:
                self.addNewTransient = self.addNewTransient_format1_sort0
                self.addNewEid       = self.addNewEidTransient_format1_sort0
            else:
                self.addNewEid = self.addNewEid_format1_sort0
            ###

        #elif self.code==[2,1,10]:
        #    self.getLength = self.getLength_format1_sort0
        #    self.addNewTransient = self.addNewTransient_format2_sort1
        #    self.addNewEid       = self.addNewEidTransient_format2_sort1
        #    self.isImaginary = True
        else:
            raise InvalidCodeError('rodStrain - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt is not None:
            self.dt = self.nonlinearFactor
            self.isTransient = True
            self.addNewTransient()
        ###

    def getLength_format1_sort0(self):
        return (16,'ifff')

    def deleteTransient(self,dt):
        del self.maxShear[dt]
        del self.avgShear[dt]
        del self.margin[dt]

    def getTransients(self):
        k = self.maxShear.keys()
        k.sort()
        return k

    def addNewTransient_format1_sort0(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        if self.dt not in self.maxShear:
            self.maxShear[self.dt] = {}
            self.avgShear[self.dt] = {}
            self.margin[self.dt]   = {}

    def addNewTransient_format2_sort1(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        raise Exeption('not supported')
        #print self.dataCode
        if self.dt not in self.axial:
            self.axial[self.dt]     = {}
            self.torsion[self.dt]   = {}

    def addNewEid_format1_sort0(self,out):
        raise Exception('not supported')
        (eid,axial,SMa,torsion,SMt) = out
        eid = (eid-self.deviceCode) // 10
        #print "Rod Strain add..."
        assert eid >= 0
        #self.eType = self.eType
        self.maxShearl[eid] = axial
        self.avgShear[eid]  = SMa
        self.margin[eid]    = torsion

    def addNewEid_format2_sort1(self,out):
        raise Exeption('not supported')
        assert eid >= 0
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode) // 10
        self.axial[eid]      = [axialReal,axialImag]
        self.torsion[eid]    = [torsionReal,torsionImag]

    def addNewEidTransient_format1_sort0(self,out):
        #print "Rod Strain add..."
        #print out
        (eid,maxShear,avgShear,margin) = out
        eid = (eid-self.deviceCode) // 10
        assert eid >= 0
        dt = self.dt

        #self.eType[eid] = self.elementType
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid]   = margin

    def addNewEidTransient_format2_sort1(self,out):
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode) // 10
        assert eid >= 0
        dt = self.dt
        self.axial[dt][eid]      = [axialReal,axialImag]
        self.torsion[dt][eid]    = [torsionReal,torsionImag]

    def __reprTransient_format1_sort0__(self):
        raise Exception('not supported')
        msg = '---TRANSIENT CSHEAR STRAINS---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['maxShear','avgShear','Margin']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,maxShears in sorted(self.maxShear.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid in sorted(maxShears):
                maxShear = self.maxShear[dt][eid]
                avgShear = self.avgShear[dt][eid]
                margin   = self.margin[dt][eid]
                msg += '%-6i %6s ' %(eid,self.eType)
                vals = [maxShear,avgShear,margin]
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

    def __reprTransient_format2_sort1__(self):
        raise Exeption('not supported')
        msg = '---COMPLEX ROD STRAINS---\n'
        msg += '%-10s %10s ' %('EID','eType')
        headers = ['axialReal','axialImag','torsionReal','torsionImag']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,axial in sorted(self.axial.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
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

        msg = '---CSHEAR STRAINS---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['maxShear','avgShear','margin']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        #print "self.code = ",self.code
        for eid in sorted(self.maxShear):
            #print self.__dict__.keys()
            maxShear = self.maxShear[eid]
            avgShear = self.avgShear[eid]
            margin   = self.margin[eid]
            msg += '%-6i %6s ' %(eid,self.eType)
            vals = [maxShear,avgShear,margin]

            for val in vals:
                if abs(val)<1e-7:
                    msg += '%8s ' %('0')
                else:
                    msg += '%8.3g ' %(val)
                ###
            msg += '\n'
        return msg
