from op2_Objects import scalarObject #,array

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
        
        if nodeID is not 'C': # centroid
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
        if nodeID is not 'C': # centroid
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
        msg = '---ISOTROPIC PLATE STRAIN---\n'
        headers = ['exx','eyy','exy','eMajor','eMinor','evm']
        msg += '%-6s %6s %8s %7s ' %('EID','eType','nodeID','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

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
                    
                    msg += '%-6i %6s %8s %7s ' %(eid,eType,nid,iLayer+1)
                    vals = [exx,eyy,exy,major,minor,evm]
                    for val in vals:
                        if abs(val)<1e-6:
                            msg += '%10s ' %('0.')
                        else:
                            msg += '%10.3g ' %(val)
                        ###
                    msg += '\n'

                    #msg += "eid=%s eType=%s nid=%s iLayer=%s exx=%-9.3g eyy=%-9.3g exy=%-9.3g evm=%-9.3g\n" %(eid,eType,nid,iLayer,exx,eyy,exy,evm)
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
        msg = '---ISOTROPIC PLATE STRESS---\n'
        headers = ['oxx','oyy','txy','ovm']
        msg += '%-6s %6s %8s %7s ' %('EID','eType','nodeID','iLayer')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

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

                    msg += '%-6i %6s %8s %7s ' %(eid,eType,nid,iLayer+1)
                    vals = [oxx,oyy,txy,ovm]
                    for val in vals:
                        if abs(val)<1e-6:
                            msg += '%10s ' %('0')
                        else:
                            msg += '%10i ' %(val)
                        ###
                    msg += '\n'
                ###
            ###
        ###
        return msg


class compositePlateStressObject(scalarObject):
    """
    ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
      ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
    """
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.eType = {}
        self.o11 = {}
        self.o22 = {}
        self.t12 = {}
        self.t1z = {}
        self.t2z = {}
        self.angle  = {}
        self.majorP = {}
        self.minorP = {}
        self.ovm = {}

    def addNewEid(self,eType,eid,o11,o22,t12,t1z,t2z,angle,majorP,minorP,ovm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        #assert eid not in self.o11
        self.eType[eid] = eType
        self.o11[eid] = [o11]
        self.o22[eid] = [o22]
        self.t12[eid] = [t12]
        self.t1z[eid] = [t1z]
        self.t2z[eid] = [t2z]
        self.angle[eid]  = [angle]
        self.majorP[eid] = [majorP]
        self.minorP[eid] = [minorP]
        self.ovm[eid]    = [ovm]
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

    def __repr__(self):
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

class compositePlateStrainObject(scalarObject):
    """
    ???
    ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
      ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
    """
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.eType = {}
        self.e11 = {}
        self.e22 = {}
        self.e12 = {}
        self.e1z = {}
        self.e2z = {}
        self.angle  = {}
        self.majorP = {}
        self.minorP = {}
        self.evm = {}

    def addNewEid(self,eType,eid,e11,e22,e12,e1z,e2z,angle,majorP,minorP,evm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        assert eid not in self.e11
        self.eType[eid] = eType
        self.e11[eid] = [e11]
        self.e22[eid] = [e22]
        self.e12[eid] = [e12]
        self.e1z[eid] = [e1z]
        self.e2z[eid] = [e2z]
        self.angle[eid]  = [angle]
        self.majorP[eid] = [majorP]
        self.minorP[eid] = [minorP]
        self.evm[eid]    = [evm]
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

    def __repr__(self):
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
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['axial','torsion','MS_tension','MS_comp']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for eid in sorted(self.axial):
            axial   = self.axial[eid]
            torsion = self.torsion[eid]
            SMa = self.SMa[eid]
            SMt = self.SMt[eid]
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
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['axial','torsion','MS_tension','MS_compression']
        for header in headers:
            msg += '%8s ' %(header)
        msg += '\n'

        for eid in sorted(self.axial):
            axial   = self.axial[eid]
            torsion = self.torsion[eid]
            SMa = self.SMa[eid]
            SMt = self.SMt[eid]
            msg += '%-6i %6s ' %(eid,self.eType)
            vals = [axial,torsion,SMa,SMt]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%8s ' %('0')
                else:
                    msg += '%8i ' %(val)
                ###
            msg += '\n'
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

            s1 = self.s1[eid]
            s2 = self.s2[eid]
            s3 = self.s3[eid]
            s4 = self.s4[eid]
            smax  = self.smax[eid]
            smin  = self.smin[eid]
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


class barStrainObject(scalarObject):
    """
                               S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )
    ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
      ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C
    """
    def __init__(self,iSubcase):
        scalarObject.__init__(self,iSubcase)
        self.eType = {}
        self.e1 = {}
        self.e2 = {}
        self.e3 = {}
        self.e4 = {}
        self.axial = {}
        self.emax = {}
        self.emin = {}
        #self.MS_tension = {}
        #self.MS_compression = {}

    def addNewEid(self,eType,eid,e1a,e2a,e3a,e4a,axial,emaxa,emina,MSt,
                                 e1b,e2b,e3b,e4b,      emaxb,eminb,MSc):
        #print "Bar Stress add..."
        self.eType[eid] = eType
        self.e1[eid] = [e1a,e1b]
        self.e2[eid] = [e2a,e2b]
        self.e3[eid] = [e3a,e3b]
        self.e4[eid] = [e4a,e4b]
        self.axial[eid] = axial
        self.emax[eid] = [emaxa,emaxb]
        self.emin[eid] = [emina,eminb]
        #self.MS_tension[eid]     = MSt
        #self.MS_compression[eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def __repr__(self):
        msg = '---BAR STRAIN---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['e1','e2','e3','e4','Axial','eMax','eMin']
        for header in headers:
            msg += '%6s ' %(header)
        msg += '\n'

        for eid,E1s in sorted(self.e1.items()):
            eType = self.eType[eid]
            axial = self.axial[eid]
            #MSt = self.MSt[eid]
            #MSc = self.MSc[eid]

            e1 = self.e1[eid]
            e2 = self.e2[eid]
            e3 = self.e3[eid]
            e4 = self.e4[eid]
            emax  = self.emax[eid]
            emin  = self.emin[eid]
            msg += '%-6i %6s ' %(eid,eType)
            vals = [e1[0],e2[0],e3[0],e4[0],axial,emax[0],emin[0]]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%6s ' %('0')
                else:
                    msg += '%6i ' %(val)
                ###
            msg += '\n'

            msg += '%s ' %(' '*13)
            vals = [e1[1],e2[1],e3[1],e4[1],'',emax[1],emin[1]]
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
        headers = ['oxx','oyy','ozz','txy','tyz','txz','ovm']
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
                ovm = self.ovm[eid][nid]
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

class solidStrainObject(scalarObject):
    """
                          S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )
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
        self.exx = {}
        self.eyy = {}
        self.ezz = {}
        self.exy = {}
        self.eyz = {}
        self.exz = {}
        #self.aCos = {}
        #self.bCos = {}
        #self.cCos = {}
        #self.pressure = {}
        self.evm = {}

    def addNewEid(self,eType,cid,eid,nodeID,exx,eyy,ezz,exy,eyz,exz,aCos,bCos,cCos,pressure,evm):
        #print "Solid Strain add..."
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
        #self.aCos[eid] = {nodeID: aCos}
        #self.bCos[eid] = {nodeID: bCos}
        #self.cCos[eid] = {nodeID: cCos}
        #self.pressure[eid] = {nodeID: pressure}
        self.evm[eid]      = {nodeID: evm}
        msg = "*eid=%s nodeID=%s vm=%g" %(eid,nodeID,evm)
        #print msg
        if nodeID==0: raise Exception(msg)

    def add(self,eid,nodeID,exx,eyy,ezz,exy,eyz,exz,aCos,bCos,cCos,pressure,evm):
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

        #self.aCos[eid][nodeID] = aCos
        #self.bCos[eid][nodeID] = bCos
        #self.cCos[eid][nodeID] = cCos
        #self.pressure[eid][nodeID] = pressure
        self.evm[eid][nodeID] = evm

        if nodeID==0: raise Exception(msg)

    def __repr__(self):
        msg = '---SOLID STRAIN---\n'
        headers = ['exx','eyy','ezz','exy','eyz','exz','evm']
        msg += '%-6s %6s %8s ' %('EID','eType','nodeID')
        for header in headers:
            msg += '%9s ' %(header)
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
                evm = self.evm[eid][nid]
                msg += '%-6i %6s %8s ' %(eid,eType,nid)
                vals = [exx,eyy,ezz,exy,eyz,exz,evm]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%9s ' %('0')
                    else:
                        msg += '%9e ' %(val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%-6s nid=%-2i exx=%-5i eyy=%-5i ezz=%-5i exy=%-5i eyz=%-5i exz=%-5i evm=%-5i\n" %(eid,eType,nid,exx,eyy,ezz,exy,eyz,exz,evm)
            ###
        ###
        return msg

