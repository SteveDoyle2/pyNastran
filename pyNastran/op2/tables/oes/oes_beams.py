import sys
from struct import unpack
from oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

class beamStressObject(stressObject):
    """
    [1,0,0]
                 S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )
                      STAT DIST/
     ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C
            1       1   0.000   -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04          
                    2   1.000   -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04          

    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CBEAM'
        
        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.xxb = {}
        self.grids = {}
        self.sxc = {}
        self.sxd = {}
        self.sxe = {}
        self.sxf = {}
        self.smax = {}
        self.smin = {}
        self.MS_tension = {}
        self.MS_compression = {}
        
        if self.code in [[1,0,0]]: # ,[1,0,1]
            #self.MS_axial   = {}
            #self.MS_torsion = {}
            self.getLength1     = self.getLength1_format1_sort0
            self.getLength2     = self.getLength2_format1_sort0
            self.getLengthTotal = self.getLengthTotal_format1_sort0

            #self.isImaginary = False
            if dt is not None:
                self.addNewTransient = self.addNewTransient_format1_sort0
                self.addNewEid       = self.addNewEidTransient_format1_sort0
                self.add             = self.addTransient_format1_sort0
                self.isTransient = True
            else:
                self.addNewEid = self.addNewEid_format1_sort0
                self.add       = self.add_format1_sort0
            ###
        #elif self.code==[2,1,0]:
        elif 0:
            self.getLength       = self.getLength_format1_sort0
            self.addNewTransient = self.addNewTransient_format2_sort1
            self.addNewEid       = self.addNewEidTransient_format2_sort1
            #self.isImaginary = True
        else:
            raise InvalidCodeError('beamStress - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt is not None:
            #self.isTransient = True
            self.dt = self.nonlinearFactor
            self.addNewTransient()
        ###

    def getLengthTotal_format1_sort0(self):
        return 444  # 44+10*40   (11 nodes)

    def getLength1_format1_sort0(self):
        return (44,'iffffffffff')

    def getLength2_format1_sort0(self):
        return (40,'ifffffffff')

    def addNewTransient_format1_sort0(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        print "addNewTransient_beam+1+0"
        #self.axial[self.dt]      = {}
        #self.MS_axial[self.dt]   = {}
        #self.torsion[self.dt]    = {}
        #self.MS_torsion[self.dt] = {}

    def addNewTransient_format2_sort1(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        print self.dataCode
        self.axial[self.dt]     = {}
        self.torsion[self.dt]   = {}

    def addNewEid_format1_sort0(self,out):
        #print "Rod Stress add..."
        (eid,grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
        eid = (eid-self.deviceCode)/10
        assert isinstance(eid,int)
        #self.axial[eid]      = axial
        #self.MS_axial[eid]   = SMa
        #self.torsion[eid]    = torsion
        #self.MS_torsion[eid] = SMt
        return eid

    def addNewEid_format2_sort1(self,out):
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        assert eid >= 0
        self.axial[eid]      = [axialReal,axialImag]
        self.torsion[eid]    = [torsionReal,torsionImag]
        return eid

    def addNewEidTransient_format1_sort0(self,out):
        (eid,grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
        eid = (eid-self.deviceCode)/10
        dt = self.dt
        assert isinstance(eid,int)
        assert eid >= 0
        #self.axial[dt][eid]      = axial
        #self.MS_axial[dt][eid]   = SMa
        #self.torsion[dt][eid]    = torsion
        #self.MS_torsion[dt][eid] = SMt
        return edi

    def addNewEidTransient_format2_sort1(self,out):
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        dt = self.dt
        assert eid >= 0
        self.axial[dt][eid]      = [axialReal,axialImag]
        self.torsion[dt][eid]    = [torsionReal,torsionImag]
        return edi

    def add_format1_sort0(self,eid,out):
        #print "Rod Stress add..."
        (grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
        #self.axial[eid]      = axial
        #self.MS_axial[eid]   = SMa
        #self.torsion[eid]    = torsion
        #self.MS_torsion[eid] = SMt

    def addTransient_format1_sort0(self,eid,out):
        #print "Rod Stress add..."
        (grid,sd,sxc,sxd,sxe,sxf,smax,smin,mst,msc) = out
        dt = self.dt
        #self.axial[dt][eid]      = axial
        #self.MS_axial[dt][eid]   = SMa
        #self.torsion[dt][eid]    = torsion
        #self.MS_torsion[dt][eid] = SMt

    def __reprTransient_format1_sort0__(self):
        msg = '---BEAM STRESSES---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['axial','torsion','MS_axial','MS_torsion']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,axial in sorted(self.axial.items()):
            msg += 'dt = %s' %(dt)
            for eid in sorted(axial):
                axial   = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                SMa     = self.MS_axial[dt][eid]
                SMt     = self.MS_torsion[dt][eid]
                msg += '%-6i %6s ' %(eid,self.eType)
                vals = [axial,torsion,SMa,SMt]
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
        msg = '---COMPLEX BEAM STRESSES---\n'
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
        msg = '---BEAM STRESSES---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['axial','torsion','MS_axial','MS_torsion']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'
        #print "self.code = ",self.code
        for eid in sorted(self.axial):
            #print self.__dict__.keys()
            axial   = self.axial[eid]
            torsion = self.torsion[eid]
            SMa     = self.MS_axial[eid]
            SMt     = self.MS_torsion[eid]
            msg += '%-6i %6s ' %(eid,self.eType)
            vals = [axial,torsion,SMa,SMt]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %('0')
                else:
                    msg += '%10i ' %(val)
                ###
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

class beamStrainObject(strainObject):
    """
    # sCode=1???
                                     S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )
    ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
      ID.        STRAIN       MARGIN        STRAIN      MARGIN

    # sCode=10
                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C O N R O D )
    ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
      ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN
       1001    1.000000E+00   1.0E+00    1.250000E+00   3.0E+00         1007    1.000000E+00   1.0E+00    1.250000E+00   3.0E+00
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)
        self.eType = 'CBEAM' #{} # 'CBEAM/CONBEAM'

        self.code = [self.formatCode,self.sortCode,self.sCode]
        
        self.isTransient = False
        if self.code == [1,0,10]:
            self.getLength = self.getLength_format1_sort0
            self.axial      = {}
            self.MS_axial   = {}
            self.torsion    = {}
            self.MS_torsion = {}
            self.isImaginary = False
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
            self.isImaginary = True
        else:
            raise InvalidCodeError('beamStrain - get the format/sort/stressCode=%s' %(self.code))
        ###
        self.dt = self.nonlinearFactor
        if dt is not None:
            self.addNewTransient()
        ###

    def getLength_format1_sort0(self):
        return (20,'iffff')

    def addNewTransient_format1_sort0(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.axial[self.dt]      = {}
        self.MS_axial[self.dt]   = {}
        self.torsion[self.dt]    = {}
        self.MS_torsion[self.dt] = {}

    def addNewTransient_format2_sort1(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        print self.dataCode
        self.axial[self.dt]     = {}
        self.torsion[self.dt]   = {}

    def addNewEid_format1_sort0(self,out):
    
        (eid,axial,SMa,torsion,SMt) = out
        eid = (eid-self.deviceCode)/10
        #print "Rod Strain add..."
        assert eid >= 0
        #self.eType = self.eType
        self.axial[eid]      = axial
        self.MS_axial[eid]   = SMa
        self.torsion[eid]    = torsion
        self.MS_torsion[eid] = SMt

    def addNewEid_format2_sort1(self,out):
        assert eid >= 0
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        self.axial[eid]      = [axialReal,axialImag]
        self.torsion[eid]    = [torsionReal,torsionImag]

    def addNewEidTransient_format1_sort0(self,out):
        #print "Rod Strain add..."
        print out
        (eid,axial,SMa,torsion,SMt) = out
        eid = (eid-self.deviceCode)/10
        assert eid >= 0
        dt = self.dt

        self.eType[eid] = self.elementType
        self.axial[dt][eid]      = axial
        self.MS_axial[dt][eid]   = SMa
        self.torsion[dt][eid]    = torsion
        self.MS_torsion[dt][eid] = SMt

    def addNewEidTransient_format2_sort1(self,out):
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode)/10
        assert eid >= 0
        dt = self.dt
        self.axial[dt][eid]      = [axialReal,axialImag]
        self.torsion[dt][eid]    = [torsionReal,torsionImag]

    def __reprTransient_format2_sort1__(self):
        msg = '---COMPLEX BEAM STRAINS---\n'
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

        msg = '---BEAM STRAINS---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['axial','torsion','MS_tension','MS_compression']
        for header in headers:
            msg += '%8s ' %(header)
        msg += '\n'

        for eid in sorted(self.axial):
            axial   = self.axial[eid]
            torsion = self.torsion[eid]
            SMa     = self.MS_axial[eid]
            SMt     = self.MS_torsion[eid]
            msg += '%-6i %6s ' %(eid,self.eType)
            vals = [axial,torsion,SMa,SMt]
            for val in vals:
                if abs(val)<1e-7:
                    msg += '%8s ' %('0')
                else:
                    msg += '%1.3g ' %(val)
                ###
            msg += '\n'
        return msg
