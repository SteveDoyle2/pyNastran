import sys
from oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

class barStressObject(stressObject):
    """
    # sCode=0
                               S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )
    ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
      ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = {}

        self.code = [self.formatCode,self.sortCode,self.sCode]
        if self.code==[1,0,0]:
            self.s1    = {}
            self.s2    = {}
            self.s3    = {}
            self.s4    = {}
            self.axial = {}
            self.smax  = {}
            self.smin  = {}
            self.MS_tension = {}
            self.MS_compression = {}
        else:
            raise InvalidCodeError('barStress - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        """
        if self.dt not in self.s1:
            self.s1[self.dt]    = {}
            self.s2[self.dt]    = {}
            self.s3[self.dt]    = {}
            self.s4[self.dt]    = {}
            self.axial[self.dt] = {}
            self.smax[self.dt]  = {}
            self.smin[self.dt]  = {}
            #self.MS_tension[self.dt]     = {}
            #self.MS_compression[self.dt] = {}

    def addNewEid(self,eType,eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                                 s1b,s2b,s3b,s4b,      smaxb,sminb,MSc):
        #print "Bar Stress add..."
        self.eType[eid] = eType
        if self.dt not in self.s1:
            self.s1[eid]    = [s1a,s1b]
            self.s2[eid]    = [s2a,s2b]
            self.s3[eid]    = [s3a,s3b]
            self.s4[eid]    = [s4a,s4b]
            self.axial[eid] = axial
            self.smax[eid]  = [smaxa,smaxb]
            self.smin[eid]  = [smina,sminb]
            #self.MS_tension[eid]     = MSt
            #self.MS_compression[eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,s1a,s2a,s3a,s4a,axial,smaxa,smina,MSt,
                                          s1b,s2b,s3b,s4b,      smaxb,sminb,MSc):
        #print "Bar Stress add..."
        dt = self.dt
        self.eType[eid] = eType
        self.s1[dt][eid]    = [s1a,s1b]
        self.s2[dt][eid]    = [s2a,s2b]
        self.s3[dt][eid]    = [s3a,s3b]
        self.s4[dt][eid]    = [s4a,s4b]
        self.axial[dt][eid] = axial
        self.smax[dt][eid]  = [smaxa,smaxb]
        self.smin[dt][eid]  = [smina,sminb]
        #self.MS_tension[dt][eid]     = MSt
        #self.MS_compression[dt][eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def __reprTransient__(self):
        msg = '---BAR STRESS---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['s1','s2','s3','s4','Axial','sMax','sMin']
        for header in headers:
            msg += '%6s ' %(header)
        msg += '\n'

        for dt,S1ss in sorted(self.s1.items()):
            msg += '%s = %s' %(self.dataCode['name'],self.dt)
            for eid,S1s in sorted(S1ss.items()):
                eType = self.eType[eid]
                axial = self.axial[dt][eid]
                #MSt = self.MSt[dt][eid]
                #MSc = self.MSc[dt][eid]

                s1   = self.s1[dt][eid]
                s2   = self.s2[dt][eid]
                s3   = self.s3[dt][eid]
                s4   = self.s4[dt][eid]
                smax = self.smax[dt][eid]
                smin = self.smin[dt][eid]
                msg += '%-6i %6s ' %(eid,eType)
                vals = [s1[0],s2[0],s3[0],s4[0],axial,smax[0],smin[0]]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%6s ' %('0')
                    else:
                        msg += '%6i ' %(val)
                    ###
                msg += '\n'

                msg += '%s ' %(' '*13)
                vals = [s1[1],s2[1],s3[1],s4[1],'',smax[1],smin[1]]
                for val in vals:
                    if isinstance(val,str):
                        msg += '%6s ' %(val)
                    elif abs(val)<1e-6:
                        msg += '%6s ' %('0')
                    else:
                        msg += '%6i ' %(val)
                    ###
                msg += '\n'


                #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i smax=%-5i smax=%-4i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial, smax[0],smin[0])
                #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s         smax=%-5i smax=%-4i\n" %(' '*4,    s1[1],s2[1],s3[1],s4[1],'    ',smax[1],smin[1])
            ###
        ###
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---BAR STRESS---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['s1','s2','s3','s4','Axial','sMax','sMin']
        for header in headers:
            msg += '%6s ' %(header)
        msg += '\n'

        for eid,S1s in sorted(self.s1.items()):
            eType = self.eType[eid]
            axial = self.axial[eid]
            #MSt = self.MSt[eid]
            #MSc = self.MSc[eid]

            s1   = self.s1[eid]
            s2   = self.s2[eid]
            s3   = self.s3[eid]
            s4   = self.s4[eid]
            smax = self.smax[eid]
            smin = self.smin[eid]
            msg += '%-6i %6s ' %(eid,eType)
            vals = [s1[0],s2[0],s3[0],s4[0],axial,smax[0],smin[0]]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%6s ' %('0')
                else:
                    msg += '%6i ' %(val)
                ###
            msg += '\n'

            msg += '%s ' %(' '*13)
            vals = [s1[1],s2[1],s3[1],s4[1],'',smax[1],smin[1]]
            for val in vals:
                if isinstance(val,str):
                    msg += '%6s ' %(val)
                elif abs(val)<1e-6:
                    msg += '%6s ' %('0')
                else:
                    msg += '%6i ' %(val)
                ###
            msg += '\n'


            #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i smax=%-5i smax=%-4i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial, smax[0],smin[0])
            #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s         smax=%-5i smax=%-4i\n" %(' '*4,    s1[1],s2[1],s3[1],s4[1],'    ',smax[1],smin[1])
        ###
        return msg


class barStrainObject(strainObject):
    """
    # sCode=10
                                    S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )
    ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
      ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C

    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)
        self.eType = {}

        self.code = [self.formatCode,self.sortCode,self.sCode]
        if self.code == [1,0,1]:
            raise InvalidCodeError('barStrain - get the format/sort/stressCode=%s' %(self.code))
            self.e1    = {}
            self.e2    = {}
            self.e3    = {}
            self.e4    = {}
            self.axial = {}
            self.emax  = {}
            self.emin  = {}
            #self.MS_tension = {}
            #self.MS_compression = {}
        elif self.code == [1,0,10]:
            self.e1    = {}
            self.e2    = {}
            self.e3    = {}
            self.e4    = {}
            self.axial = {}
            self.emax  = {}
            self.emin  = {}
            self.MS_tension = {}
            self.MS_compression = {}
        else:
            raise InvalidCodeError('barStrain - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        """
        if self.dt not in self.e1:
            self.e1[self.dt]    = {}
            self.e2[self.dt]    = {}
            self.e3[self.dt]    = {}
            self.e4[self.dt]    = {}
            self.axial[self.dt] = {}
            self.emax[self.dt]  = {}
            self.emin[self.dt]  = {}
            #self.MS_tension[self.dt]     = {}
            #self.MS_compression[self.dt] = {}

    def addNewEid(self,eType,eid,e1a,e2a,e3a,e4a,axial,emaxa,emina,MSt,
                                 e1b,e2b,e3b,e4b,      emaxb,eminb,MSc):
        #print "Bar Stress add..."
        self.eType[eid] = eType
        if self.dt not in self.e1:
            self.e1[eid]    = [e1a,e1b]
            self.e2[eid]    = [e2a,e2b]
            self.e3[eid]    = [e3a,e3b]
            self.e4[eid]    = [e4a,e4b]
            self.axial[eid] = axial
            self.emax[eid]  = [emaxa,emaxb]
            self.emin[eid]  = [emina,eminb]
            #self.MS_tension[eid]     = MSt
            #self.MS_compression[eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def addNewEidTransient(self,eType,eid,e1a,e2a,e3a,e4a,axial,emaxa,emina,MSt,
                                      e1b,e2b,e3b,e4b,      emaxb,eminb,MSc):
        #print "Bar Stress add..."
        dt = self.dt
        self.eType[eid] = eType
        self.e1[dt][eid]    = [e1a,e1b]
        self.e2[dt][eid]    = [e2a,e2b]
        self.e3[dt][eid]    = [e3a,e3b]
        self.e4[dt][eid]    = [e4a,e4b]
        self.axial[dt][eid] = axial
        self.emax[dt][eid]  = [emaxa,emaxb]
        self.emin[dt][eid]  = [emina,eminb]
        #self.MS_tension[dt][eid]     = MSt
        #self.MS_compression[dt][eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def __reprTransient__(self):
        """
        @warning untested
        """
        msg = '---BAR STRAIN---\n'
        msg += '%-8s %6s ' %('EID','eType')
        headers = ['e1','e2','e3','e4','Axial','eMax','eMin']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,E1s in sorted(self.e1.items()):
            print "%s = %s\n" %(self.dataCode['name'],self.dt)
            for eid,e1s in sorted(Els.items()):
                eType = self.eType[eid]
                axial = self.axial[dt][eid]
                #MSt  = self.MS_tension[dt][eid]
                #MSc  = self.MS_compression[dt][eid]

                e1   = self.e1[dt][eid]
                e2   = self.e2[dt][eid]
                e3   = self.e3[dt][eid]
                e4   = self.e4[dt][eid]
                emax = self.emax[dt][eid]
                emin = self.emin[dt][eid]
                msg += '%-8i %6s ' %(eid,eType)
                vals = [e1[0],e2[0],e3[0],e4[0],axial,emax[0],emin[0]]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10.3g ' %(val)
                    ###
                msg += '\n'

                msg += '%s ' %(' '*17)
                vals = [e1[1],e2[1],e3[1],e4[1],'',emax[1],emin[1]]
                for val in vals:
                    if isinstance(val,str):
                        msg += '%10s ' %(val)
                    elif abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10.3g ' %(val)
                    ###
                msg += '\n'

                #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i smax=%-5i smax=%-4i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial, smax[0],smin[0])
                #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s         smax=%-5i smax=%-4i\n" %(' '*4,    s1[1],s2[1],s3[1],s4[1],'    ',smax[1],smin[1])
            ###
        ###
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---BAR STRAIN---\n'
        msg += '%-8s %6s ' %('EID','eType')
        headers = ['e1','e2','e3','e4','Axial','eMax','eMin']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for eid,E1s in sorted(self.e1.items()):
            eType = self.eType[eid]
            axial = self.axial[eid]
            #MSt  = self.MS_tension[eid]
            #MSc  = self.MS_compression[eid]

            e1   = self.e1[eid]
            e2   = self.e2[eid]
            e3   = self.e3[eid]
            e4   = self.e4[eid]
            emax = self.emax[eid]
            emin = self.emin[eid]
            msg += '%-8i %6s ' %(eid,eType)
            vals = [e1[0],e2[0],e3[0],e4[0],axial,emax[0],emin[0]]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %('0')
                else:
                    msg += '%10.3g ' %(val)
                ###
            msg += '\n'

            msg += '%s ' %(' '*17)
            vals = [e1[1],e2[1],e3[1],e4[1],'',emax[1],emin[1]]
            for val in vals:
                if isinstance(val,str):
                    msg += '%10s ' %(val)
                elif abs(val)<1e-6:
                    msg += '%10s ' %('0')
                else:
                    msg += '%10.3g ' %(val)
                ###
            msg += '\n'

            #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i smax=%-5i smax=%-4i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial, smax[0],smin[0])
            #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s         smax=%-5i smax=%-4i\n" %(' '*4,    s1[1],s2[1],s3[1],s4[1],'    ',smax[1],smin[1])
        ###
        return msg
