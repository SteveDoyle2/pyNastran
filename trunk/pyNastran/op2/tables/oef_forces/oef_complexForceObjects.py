class ComplexRodForce(object): # 1-ROD, 3-TUBE, 10-CONROD
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.axialForceReal = {}
        self.axialForceImag = {}
        self.torqueReal = {}
        self.torqueImag = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.axialForceReal[dt] = {}
        self.axialForceImag[dt] = {}
        self.torqueReal[dt] = {}
        self.torqueImag[dt] = {}

    def add(self,dt,data):
        [eid,axialForceReal,torqueReal,axialForceImag,torqueImag] = data

        #self.eType[eid] = eType
        self.axialForceReal[eid] = axialForceReal
        self.axialForceImag[eid] = axialForceImag
        self.torqueReal[eid] = torqueReal
        self.torqueImag[eid] = torqueImag

    def addSort1(self,dt,data):
        [eid,axialForceReal,torqueReal,axialForceImag,torqueImag] = data
        if dt not in self.axialForceReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForceReal[dt][eid] = axialForceReal
        self.axialForceImag[dt][eid] = axialForceImag
        self.torqueReal[dt][eid] = torqueReal
        self.torqueImag[dt][eid] = torqueImag

    def addSort2(self,eid,data):
        [dt,axialForceReal,torqueReal,axialForceImag,torqueImag] = data
        if dt not in self.axialForceReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForceReal[dt][eid] = axialForceReal
        self.axialForceImag[dt][eid] = axialForceImag
        self.torqueReal[dt][eid] = torqueReal
        self.torqueImag[dt][eid] = torqueImag

    def __repr__(self):
        return str(self.axialForceReal)

class ComplexCBEAMForce(object): # 2-CBEAM
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.bendingMomentReal = {}
        self.shearReal = {}
        self.axialReal = {}
        self.totalTorqueReal = {}
        self.warpingTorqueReal = {}

        self.bendingMomentImag = {}
        self.shearImag = {}
        self.axialImag = {}
        self.totalTorqueImag = {}
        self.warpingTorqueImag = {}

        if isSort1:
            if dt is not None:
                self.addNewElement = self.addNewElementSort1
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.addNewElement = self.addNewElementSort2
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.bendingMomentReal[dt] = {}
        self.shearReal[dt] = {}
        self.axialReal[dt] = {}
        self.totalTorqueReal[dt] = {}
        self.warpingTorqueReal[dt] = {}

        self.bendingMomentImag[dt] = {}
        self.shearImag[dt] = {}
        self.axialImag[dt] = {}
        self.totalTorqueImag[dt] = {}
        self.warpingTorqueImag[dt] = {}

    def addNewElement(self,dt,data):
        [eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                    bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi] = data
        print "CBEAM addnew",data
        #self.eType[eid] = eType
        self.bendingMomentReal[eid] = {sd:[bm1r,bm2r]}
        self.shearReal[eid] = {sd:[ts1r,ts2r]}
        self.axialReal[eid] = {sd:afr}
        self.totalTorqueReal[eid] = {sd:ttrqr}
        self.warpingTorqueReal[eid] = {sd:wtrqr}

        self.bendingMomentImag[eid] = {sd:[bm1i,bm2i]}
        self.shearImag[eid] = {sd:[ts1i,ts2i]}
        self.axialImag[eid] = {sd:afi}
        self.totalTorqueImag[eid] = {sd:ttrqi}
        self.warpingTorqueImag[eid] = {sd:wtrqi}

    def add(self,dt,data):
        [eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                    bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi] = data
        print "CBEAM add   ",data

        #self.eType[eid] = eType
        self.bendingMomentReal[eid][sd] = [bm1r,bm2r]
        self.shearReal[eid][sd] = [ts1r,ts2r]
        self.axialReal[eid][sd] = afr
        self.totalTorqueReal[eid][sd] = ttrqr
        self.warpingTorqueReal[eid][sd] = wtrqr

        self.bendingMomentImag[eid][sd] = [bm1i,bm2i]
        self.shearImag[eid][sd] = [ts1i,ts2i]
        self.axialImag[eid][sd] = afi
        self.totalTorqueImag[eid][sd] = ttrqi
        self.warpingTorqueImag[eid][sd] = wtrqi

    def addNewElementSort1(self,dt,data):
        [eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                    bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi] = data

        self._fillNewObject(dt,eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                          bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi)

    def addSort1(self,dt,data):
        [eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                    bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi] = data
        self._fillObject(dt,eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                       bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi)

    def addNewElementSort2(self,eid,data):
        [dt,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                   bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi] = data
        self._fillNewObject(dt,eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                          bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi)

    def addSort2(self,eid,data):
        [dt,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                   bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi] = data
        self._fillObject(dt,eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                       bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi)

    def _fillNewObject(self,dt,eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                          bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi):
        if dt not in self.axialReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMomentReal[dt][eid] = {sd:[bm1r,bm2r]}
        self.shearReal[dt][eid] = {sd:[ts1r,ts2r]}
        self.axialReal[dt][eid] = {sd:afr}
        self.totalTorqueReal[dt][eid] = {sd:ttrqr}
        self.warpingTorqueReal[dt][eid] = {sd:wtrqr}

        self.bendingMomentImag[dt][eid] = {sd:[bm1i,bm2i]}
        self.shearImag[dt][eid] = {sd:[ts1i,ts2i]}
        self.axialImag[dt][eid] = {sd:afi}
        self.totalTorqueImag[dt][eid] = {sd:ttrqi}
        self.warpingTorqueImag[dt][eid] = {sd:wtrqi}

    def _fillObject(self,dt,eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                       bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi):
        #if dt not in self.axialReal:
            #self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMomentReal[dt][eid][sd] = [bm1r,bm2r]
        self.shearReal[dt][eid][sd] = [ts1r,ts2r]
        self.axialReal[dt][eid][sd] = afr
        self.totalTorqueReal[dt][eid][sd] = ttrqr
        self.warpingTorqueReal[dt][eid][sd] = wtrqr

        self.bendingMomentImag[dt][eid][sd] = [bm1i,bm2i]
        self.shearImag[dt][eid][sd] = [ts1i,ts2i]
        self.axialImag[dt][eid][sd] = afi
        self.totalTorqueImag[dt][eid][sd] = ttrqi
        self.warpingTorqueImag[dt][eid][sd] = wtrqi

    def __repr__(self):
        return str(self.axialReal)

class ComplexCShearForce(object): # 4-CSHEAR
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.force41Real = {}
        self.force14Real = {}
        self.force21Real = {}
        self.force12Real = {}
        self.force32Real = {}
        self.force23Real = {}
        self.force43Real = {}
        self.force34Real = {}
        self.kickForce1Real = {}
        self.kickForce2Real = {}
        self.kickForce3Real = {}
        self.kickForce4Real = {}
        self.shear12Real = {}
        self.shear23Real = {}
        self.shear34Real = {}
        self.shear41Real = {}

        self.force41Imag = {}
        self.force14Imag = {}
        self.force21Imag = {}
        self.force12Imag = {}
        self.force32Imag = {}
        self.force23Imag = {}
        self.force43Imag = {}
        self.force34Imag = {}
        self.kickForce1Imag = {}
        self.kickForce2Imag = {}
        self.kickForce3Imag = {}
        self.kickForce4Imag = {}
        self.shear12Imag = {}
        self.shear23Imag = {}
        self.shear34Imag = {}
        self.shear41Imag = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.force41Real[dt] = {}
        self.force14Real[dt] = {}
        self.force21Real[dt] = {}
        self.force12Real[dt] = {}
        self.force32Real[dt] = {}
        self.force23Real[dt] = {}
        self.force43Real[dt] = {}
        self.force34Real[dt] = {}
        self.kickForce1Real[dt] = {}
        self.kickForce2Real[dt] = {}
        self.kickForce3Real[dt] = {}
        self.kickForce4Real[dt] = {}
        self.shear12Real[dt] = {}
        self.shear23Real[dt] = {}
        self.shear34Real[dt] = {}
        self.shear41Real[dt] = {}

        self.force41Imag[dt] = {}
        self.force14Imag[dt] = {}
        self.force21Imag[dt] = {}
        self.force12Imag[dt] = {}
        self.force32Imag[dt] = {}
        self.force23Imag[dt] = {}
        self.force43Imag[dt] = {}
        self.force34Imag[dt] = {}
        self.kickForce1Imag[dt] = {}
        self.kickForce2Imag[dt] = {}
        self.kickForce3Imag[dt] = {}
        self.kickForce4Imag[dt] = {}
        self.shear12Imag[dt] = {}
        self.shear23Imag[dt] = {}
        self.shear34Imag[dt] = {}
        self.shear41Imag[dt] = {}

    def add(self,dt,data):
        [eid,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,
             f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,
             kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
             kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i] = data
        #self.eType[eid] = eType
        self.force41Real[eid] = f41r
        self.force14Real[eid] = f14r
        self.force21Real[eid] = f21r
        self.force12Real[eid] = f12r
        self.force32Real[eid] = f32r
        self.force23Real[eid] = f23r
        self.force43Real[eid] = f43r
        self.force34Real[eid] = f34r
        self.kickForce1Real[eid] = kf1r
        self.kickForce2Real[eid] = kf2r
        self.kickForce3Real[eid] = kf3r
        self.kickForce4Real[eid] = kf4r
        self.shear12Real[eid] = s12r
        self.shear23Real[eid] = s23r
        self.shear34Real[eid] = s34r
        self.shear41Real[eid] = s41r

        self.force41Imag[eid] = f41i
        self.force14Imag[eid] = f14i
        self.force21Imag[eid] = f21i
        self.force12Imag[eid] = f12i
        self.force32Imag[eid] = f32i
        self.force23Imag[eid] = f23i
        self.force43Imag[eid] = f43i
        self.force34Imag[eid] = f34i
        self.kickForce1Imag[eid] = kf1i
        self.kickForce2Imag[eid] = kf2i
        self.kickForce3Imag[eid] = kf3i
        self.kickForce4Imag[eid] = kf4i
        self.shear12Imag[eid] = s12i
        self.shear23Imag[eid] = s23i
        self.shear34Imag[eid] = s34i
        self.shear41Imag[eid] = s41i
        
    def addSort1(self,dt,data):
        [eid,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,
             f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,
             kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
             kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i] = data
        self._fillObject(dt,eid,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,
                                f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,
                                kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
                                kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i)

    def addSort2(self,eid,data):
        [dt,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,
            f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,
            kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
            kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i] = data

        self._fillObject(dt,eid,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,
                                f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,
                                kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
                                kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i)

    def _fillObject(self,dt,eid,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,
                                f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,
                                kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
                                kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i):
        if dt not in self.force41Real:
            self.addNewTransient(dt)
        #self.eType[eid] = eType
        self.force41Real[dt][eid] = f41r
        self.force14Real[dt][eid] = f14r
        self.force21Real[dt][eid] = f21r
        self.force12Real[dt][eid] = f12r
        self.force32Real[dt][eid] = f32r
        self.force23Real[dt][eid] = f23r
        self.force43Real[dt][eid] = f43r
        self.force34Real[dt][eid] = f34r
        self.kickForce1Real[dt][eid] = kf1r
        self.kickForce2Real[dt][eid] = kf2r
        self.kickForce3Real[dt][eid] = kf3r
        self.kickForce4Real[dt][eid] = kf4r
        self.shear12Real[dt][eid] = s12r
        self.shear23Real[dt][eid] = s23r
        self.shear34Real[dt][eid] = s34r
        self.shear41Real[dt][eid] = s41r

        self.force41Imag[dt][eid] = f41i
        self.force14Imag[dt][eid] = f14i
        self.force21Imag[dt][eid] = f21i
        self.force12Imag[dt][eid] = f12i
        self.force32Imag[dt][eid] = f32i
        self.force23Imag[dt][eid] = f23i
        self.force43Imag[dt][eid] = f43i
        self.force34Imag[dt][eid] = f34i
        self.kickForce1Imag[dt][eid] = kf1i
        self.kickForce2Imag[dt][eid] = kf2i
        self.kickForce3Imag[dt][eid] = kf3i
        self.kickForce4Imag[dt][eid] = kf4i
        self.shear12Imag[dt][eid] = s12i
        self.shear23Imag[dt][eid] = s23i
        self.shear34Imag[dt][eid] = s34i
        self.shear41Imag[dt][eid] = s41i

    def __repr__(self):
        return str(self.force41)

class ComplexSpringForce(object): # 11-CELAS1,12-CELAS2,13-CELAS3, 14-CELAS4
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.force = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.force[dt] = {}

    def add(self,dt,data):
        [eid,forceReal,forceImag] = data

        #self.eType[eid] = eType
        self.force[eid] = [forceReal,forceImag]

    def addSort1(self,dt,data):
        [eid,forceReal,forceImag] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = [forceReal,forceImag]

    def addSort2(self,eid,data):
        [dt,forceReal,forceImag] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = [forceReal,forceImag]

    def __repr__(self):
        return str(self.force)

class ComplexViscForce(object): # 24-CVISC
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.axialForceReal = {}
        self.axialForceImag = {}
        self.torqueReal = {}
        self.torqueImag = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.axialForceReal[dt] = {}
        self.axialForceImag[dt] = {}
        self.torqueReal[dt] = {}
        self.torqueImag[dt] = {}

    def add(self,dt,data):
        [eid,axialForceReal,torqueReal,axialForceImag,torqueImag] = data

        #self.eType[eid] = eType
        self.axialForceReal[eid] = axialForceReal
        self.axialForceImag[eid] = axialForceImag
        self.torqueReal[eid] = torqueReal
        self.torqueImag[eid] = torqueImag

    def addSort1(self,dt,data):
        [eid,axialForceReal,torqueReal,axialForceImag,torqueImag] = data
        if dt not in self.axialForceReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForceReal[dt][eid] = axialForceReal
        self.axialForceImag[dt][eid] = axialForceImag
        self.torqueReal[dt][eid] = torqueReal
        self.torqueImag[dt][eid] = torqueImag

    def addSort2(self,eid,data):
        [dt,axialForceReal,torqueReal,axialForceImag,torqueImag] = data
        if dt not in self.axialForceReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForceReal[dt][eid] = axialForceReal
        self.axialForceImag[dt][eid] = axialForceImag
        self.torqueReal[dt][eid] = torqueReal
        self.torqueImag[dt][eid] = torqueImag

    def __repr__(self):
        return str(self.axialForceReal)

class ComplexPlateForce(object): # 33-CQUAD4, 74-CTRIA3
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.mx = {}
        self.my = {}
        self.mxy = {}
        self.bmx = {}
        self.bmy = {}
        self.bmxy = {}
        self.tx = {}
        self.ty = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.mx[dt] = {}
        self.my[dt] = {}
        self.mxy[dt] = {}
        self.bmx[dt] = {}
        self.bmy[dt] = {}
        self.bmxy[dt] = {}
        self.tx[dt] = {}
        self.ty[dt] = {}

    def add(self,dt,data):
        [eid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
             mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data

        #self.eType[eid] = eType
        self.mx[eid] = [mxr,mxi]
        self.my[eid] = [myr,myi]
        self.mxy[eid] = [mxyr,mxyi]
        self.bmx[eid] = [bmxr,bmxi]
        self.bmy[eid] = [bmyr,bmyi]
        self.bmxy[eid] = [bmxyr,bmxyi]
        self.tx[eid] = [txr,txi]
        self.ty[eid] = [tyr,tyi]

    def addSort1(self,dt,data):
        [eid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
             mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid] = [mxr,mxi]
        self.my[dt][eid] = [myr,myi]
        self.mxy[dt][eid] = [mxyr,mxyi]
        self.bmx[dt][eid] = [bmxr,bmxi]
        self.bmy[dt][eid] = [bmyr,bmyi]
        self.bmxy[dt][eid] = [bmxyr,bmxyi]
        self.tx[dt][eid] = [txr,txi]
        self.ty[dt][eid] = [tyr,tyi]

    def addSort2(self,eid,data):
        [dt,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
            mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid] = [mxr,mxi]
        self.my[dt][eid] = [myr,myi]
        self.mxy[dt][eid] = [mxyr,mxyi]
        self.bmx[dt][eid] = [bmxr,bmxi]
        self.bmy[dt][eid] = [bmyr,bmyi]
        self.bmxy[dt][eid] = [bmxyr,bmxyi]
        self.tx[dt][eid] = [txr,txi]
        self.ty[dt][eid] = [tyr,tyi]

    def __repr__(self):
        return str(self.mx)

class ComplexPLATE2Force(object): # 64-CQUAD8, 75-CTRIA6, 82-CQUADR
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.term = {}
        self.ngrids = {}
        self.mxReal = {}
        self.myReal = {}
        self.mxyReal = {}
        self.bmxReal = {}
        self.bmyReal = {}
        self.bmxyReal = {}
        self.txReal = {}
        self.tyReal = {}

        self.mxImag = {}
        self.myImag = {}
        self.mxyImag = {}
        self.bmxImag = {}
        self.bmyImag = {}
        self.bmxyImag = {}
        self.txImag = {}
        self.tyImag = {}

        if isSort1:
            if dt is not None:
                self.addNewElement = self.addNewElementSort1
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.addNewElement = self.addNewElementSort2
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.term[dt] = {}
        self.ngrids[dt] = {}

        self.mxReal[dt] = {}
        self.myReal[dt] = {}
        self.mxyReal[dt] = {}
        self.bmxReal[dt] = {}
        self.bmyReal[dt] = {}
        self.bmxyReal[dt] = {}
        self.txReal[dt] = {}
        self.tyReal[dt] = {}

        self.mxImag[dt] = {}
        self.myImag[dt] = {}
        self.mxyImag[dt] = {}
        self.bmxImag[dt] = {}
        self.bmyImag[dt] = {}
        self.bmxyImag[dt] = {}
        self.txImag[dt] = {}
        self.tyImag[dt] = {}

    def addNewElement(self,eid,dt,data):
        #print "eid = ",eid
        [term,nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                  mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data

        #self.eType[eid] = eType
        self.term[eid] = term
        self.ngrids[eid] = nid

        self.mxReal[eid] = [mxr]
        self.myReal[eid] = [myr]
        self.mxyReal[eid] = [mxyr]
        self.bmxReal[eid] = [bmxr]
        self.bmyReal[eid] = [bmyr]
        self.bmxyReal[eid] = [bmxyr]
        self.txReal[eid] = [txr]
        self.tyReal[eid] = [tyr]

        self.mxImag[eid] = [mxi]
        self.myImag[eid] = [myi]
        self.mxyImag[eid] = [mxyi]
        self.bmxImag[eid] = [bmxi]
        self.bmyImag[eid] = [bmyi]
        self.bmxyImag[eid] = [bmxyi]
        self.txImag[eid] = [txi]
        self.tyImag[eid] = [tyi]


    def add(self,eid,dt,data):
        [nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
             mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data

        #self.eType[eid] = eType
        #print "mx = ",self.mx,mx
        self.mxReal[eid].append(mxr)
        self.myReal[eid].append(myr)
        self.mxyReal[eid].append(mxyr)
        self.bmxReal[eid].append(bmxr)
        self.bmyReal[eid].append(bmyr)
        self.bmxyReal[eid].append(bmxyr)
        self.txReal[eid].append(txr)
        self.tyReal[eid].append(tyr)

        self.mxImag[eid].append(mxi)
        self.myImag[eid].append(myi)
        self.mxyImag[eid].append(mxyi)
        self.bmxImag[eid].append(bmxi)
        self.bmyImag[eid].append(bmyi)
        self.bmxyImag[eid].append(bmxyi)
        self.txImag[eid].append(txi)
        self.tyImag[eid].append(tyi)

    def addNewElementSort1(self,eid,dt,data):
        [term,nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                  mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mxReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.term[dt][eid] = term
        self.ngrids[dt][eid] = nid
        self.mxReal[dt][eid] = [mxr]
        self.myReal[dt][eid] = [myr]
        self.mxyReal[dt][eid] = [mxyr]
        self.bmxReal[dt][eid] = [bmxr]
        self.bmyReal[dt][eid] = [bmyr]
        self.bmxyReal[dt][eid] = [bmxyr]
        self.txReal[dt][eid] = [txr]
        self.tyReal[dt][eid] = [tyr]

        self.mxImag[dt][eid] = [mxi]
        self.myImag[dt][eid] = [myi]
        self.mxyImag[dt][eid] = [mxyi]
        self.bmxImag[dt][eid] = [bmxi]
        self.bmyImag[dt][eid] = [bmyi]
        self.bmxyImag[dt][eid] = [bmxyi]
        self.txImag[dt][eid] = [txi]
        self.tyImag[dt][eid] = [tyi]

    def addSort1(self,eid,dt,data):
        [nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
             mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mxReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mxReal[dt][eid].append(mxr)
        self.myReal[dt][eid].append(myr)
        self.mxyReal[dt][eid].append(mxyr)
        self.bmxReal[dt][eid].append(bmxr)
        self.bmyReal[dt][eid].append(bmyr)
        self.bmxyReal[dt][eid].append(bmxyr)
        self.txReal[dt][eid].append(txr)
        self.tyReal[dt][eid].append(tyr)

        self.mxImag[dt][eid].append(mxi)
        self.myImag[dt][eid].append(myi)
        self.mxyImag[dt][eid].append(mxyi)
        self.bmxImag[dt][eid].append(bmxi)
        self.bmyImag[dt][eid].append(bmyi)
        self.bmxyImag[dt][eid].append(bmxyi)
        self.txImag[dt][eid].append(txi)
        self.tyImag[dt][eid].append(tyi)

    def addNewElementSort2(self,dt,eid,data):
        [term,nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                  mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mxReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.term[dt][eid] = term
        self.ngrids[dt][eid] = nid

        self.mxReal[dt][eid] = [mxr]
        self.myReal[dt][eid] = [myr]
        self.mxyReal[dt][eid] = [mxyr]
        self.bmxReal[dt][eid] = [bmxr]
        self.bmyReal[dt][eid] = [bmyr]
        self.bmxyReal[dt][eid] = [bmxyr]
        self.txReal[dt][eid] = [txr]
        self.tyReal[dt][eid] = [tyr]

        self.mxImag[dt][eid] = [mxi]
        self.myImag[dt][eid] = [myi]
        self.mxyImag[dt][eid] = [mxyi]
        self.bmxImag[dt][eid] = [bmxi]
        self.bmyImag[dt][eid] = [bmyi]
        self.bmxyImag[dt][eid] = [bmxyi]
        self.txImag[dt][eid] = [txi]
        self.tyImag[dt][eid] = [tyi]

    def addSort2(self,dt,eid,data):
        [nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
             mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mxReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mxReal[dt][eid].append(mxr)
        self.myReal[dt][eid].append(myr)
        self.mxyReal[dt][eid].append(mxyr)
        self.bmxReal[dt][eid].append(bmxr)
        self.bmyReal[dt][eid].append(bmyr)
        self.bmxyReal[dt][eid].append(bmxyr)
        self.txReal[dt][eid].append(txr)
        self.tyReal[dt][eid].append(tyr)

        self.mxImag[dt][eid].append(mxi)
        self.myImag[dt][eid].append(myi)
        self.mxyImag[dt][eid].append(mxyi)
        self.bmxImag[dt][eid].append(bmxi)
        self.bmyImag[dt][eid].append(bmyi)
        self.bmxyImag[dt][eid].append(bmxyi)
        self.txImag[dt][eid].append(txi)
        self.tyImag[dt][eid].append(tyi)

    def __repr__(self):
        return str(self.mxReal)


class ComplexCBARForce(object): # 34-CBAR
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.bendingMomentAReal = {}
        self.bendingMomentBReal = {}
        self.shearReal = {}
        self.axialReal = {}
        self.torqueReal = {}

        self.bendingMomentAImag = {}
        self.bendingMomentBImag = {}
        self.shearImag = {}
        self.axialImag = {}
        self.torqueImag = {}


        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.bendingMomentAReal[dt] = {}
        self.bendingMomentBReal[dt] = {}
        self.shearReal[dt] = {}
        self.axialReal[dt] = {}
        self.torqueReal[dt] = {}

        self.bendingMomentAImag[dt] = {}
        self.bendingMomentBImag[dt] = {}
        self.shearImag[dt] = {}
        self.axialImag[dt] = {}
        self.torqueImag[dt] = {}

    def add(self,dt,data):
        [eid,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
             bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi] = data

        #self.eType[eid] = eType
        self.bendingMomentAReal[eid] = [bm1ar,bm2ar]
        self.bendingMomentBReal[eid] = [bm1br,bm2br]
        self.shearReal[eid] = [ts1r,ts2r]
        self.axialReal[eid] = afr
        self.torqueReal[eid] = trqr

        self.bendingMomentAImag[eid] = [bm1ai,bm2ai]
        self.bendingMomentBImag[eid] = [bm1bi,bm2bi]
        self.shearImag[eid] = [ts1i,ts2i]
        self.axialImag[eid] = afi
        self.torqueImag[eid] = trqi

    def addSort1(self,dt,data):
        [eid,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
             bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi] = data
        if dt not in self.axialReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMomentAReal[dt][eid] = [bm1ar,bm2ar]
        self.bendingMomentBReal[dt][eid] = [bm1br,bm2br]
        self.shearReal[dt][eid] = [ts1r,ts2r]
        self.axialReal[dt][eid] = afr
        self.torqueReal[dt][eid] = trqr

        self.bendingMomentAImag[dt][eid] = [bm1ai,bm2ai]
        self.bendingMomentBImag[dt][eid] = [bm1bi,bm2bi]
        self.shearImag[dt][eid] = [ts1i,ts2i]
        self.axialImag[dt][eid] = afi
        self.torqueImag[dt][eid] = trqi

    def addSort2(self,eid,data):
        [dt,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
            bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi] = data
        if dt not in self.axialReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMomentAReal[dt][eid] = [bm1ar,bm2ar]
        self.bendingMomentBReal[dt][eid] = [bm1br,bm2br]
        self.shearReal[dt][eid] = [ts1r,ts2r]
        self.axialReal[dt][eid] = afr
        self.torqueReal[dt][eid] = trqr

        self.bendingMomentAImag[dt][eid] = [bm1ai,bm2ai]
        self.bendingMomentBImag[dt][eid] = [bm1bi,bm2bi]
        self.shearImag[dt][eid] = [ts1i,ts2i]
        self.axialImag[dt][eid] = afi
        self.torqueImag[dt][eid] = trqi

    def __repr__(self):
        return str(self.axialReal)

class ComplexBendForce(object): # 69-CBEND
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.nodeIDs = {}
        self.bendingMoment1Real = {}
        self.bendingMoment2Real = {}
        self.shearPlane1Real = {}
        self.shearPlane2Real = {}
        self.axialReal  = {}
        self.torqueReal = {}

        self.bendingMoment1Imag = {}
        self.bendingMoment2Imag = {}
        self.shearPlane1Imag = {}
        self.shearPlane2Imag = {}
        self.axialImag  = {}
        self.torqueImag = {}
        
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.bendingMoment1Real[dt] = {}
        self.bendingMoment2Real[dt] = {}
        self.shearPlane1Real[dt] = {}
        self.shearPlane2Real[dt] = {}
        self.axialReal[dt]  = {}
        self.torqueReal[dt] = {}

        self.bendingMoment1Imag[dt] = {}
        self.bendingMoment2Imag[dt] = {}
        self.shearPlane1Imag[dt] = {}
        self.shearPlane2Imag[dt] = {}
        self.axialImag[dt]  = {}
        self.torqueImag[dt] = {}

    def add(self,dt,data):
        [eid,nidA,bm1Ar,bm2Ar,sp1Ar,sp2Ar,axialAr,torqueAr,
                  bm1Ai,bm2Ai,sp1Ai,sp2Ai,axialAi,torqueAi,
             nidB,bm1Br,bm2Br,sp1Br,sp2Br,axialBr,torqueBr,
                  bm1Bi,bm2Bi,sp1Bi,sp2Bi,axialBi,torqueBi] = data

        #self.eType[eid] = eType
        self.nodeIDs[eid] = [nidA,nidB]
        self.bendingMoment1Real[eid] = [bm1Ar,bm1Br]
        self.bendingMoment2Real[eid] = [bm2Ar,bm2Br]
        self.shearPlane1Real[eid] = [sp1Ar,sp1Br]
        self.shearPlane2Real[eid] = [sp2Ar,sp2Br]
        self.axialReal[eid]  = [axialAr,axialBr]
        self.torqueReal[eid] = [torqueAr,torqueBr]

        self.bendingMoment1Imag[eid] = [bm1Ai,bm1Bi]
        self.bendingMoment2Imag[eid] = [bm2Ai,bm2Bi]
        self.shearPlane1Imag[eid] = [sp1Ai,sp1Bi]
        self.shearPlane2Imag[eid] = [sp2Ai,sp2Bi]
        self.axialImag[eid]  = [axialAi,axialBi]
        self.torqueImag[eid] = [torqueAi,torqueBi]

    def addSort1(self,dt,data):
        [eid,nidA,bm1Ar,bm2Ar,sp1Ar,sp2Ar,axialAr,torqueAr,
                  bm1Ai,bm2Ai,sp1Ai,sp2Ai,axialAi,torqueAi,
             nidB,bm1Br,bm2Br,sp1Br,sp2Br,axialBr,torqueBr,
                  bm1Bi,bm2Bi,sp1Bi,sp2Bi,axialBi,torqueBi] = data
        self._fillObject(dt,eid,nidA,bm1Ar,bm2Ar,sp1Ar,sp2Ar,axialAr,torqueAr,
                                     bm1Ai,bm2Ai,sp1Ai,sp2Ai,axialAi,torqueAi,
                                nidB,bm1Br,bm2Br,sp1Br,sp2Br,axialBr,torqueBr,
                                     bm1Bi,bm2Bi,sp1Bi,sp2Bi,axialBi,torqueBi)

    def addSort2(self,eid,data):
        [dt,nidA,bm1Ar,bm2Ar,sp1Ar,sp2Ar,axialAr,torqueAr,
                 bm1Ai,bm2Ai,sp1Ai,sp2Ai,axialAi,torqueAi,
            nidB,bm1Br,bm2Br,sp1Br,sp2Br,axialBr,torqueBr,
                 bm1Bi,bm2Bi,sp1Bi,sp2Bi,axialBi,torqueBi] = data
        self._fillObject(dt,eid,nidA,bm1Ar,bm2Ar,sp1Ar,sp2Ar,axialAr,torqueAr,
                                     bm1Ai,bm2Ai,sp1Ai,sp2Ai,axialAi,torqueAi,
                                nidB,bm1Br,bm2Br,sp1Br,sp2Br,axialBr,torqueBr,
                                     bm1Bi,bm2Bi,sp1Bi,sp2Bi,axialBi,torqueBi)

    def _fillObject(self,dt,eid,nidA,bm1Ar,bm2Ar,sp1Ar,sp2Ar,axialAr,torqueAr,
                                     bm1Ai,bm2Ai,sp1Ai,sp2Ai,axialAi,torqueAi,
                                nidB,bm1Br,bm2Br,sp1Br,sp2Br,axialBr,torqueBr,
                                     bm1Bi,bm2Bi,sp1Bi,sp2Bi,axialBi,torqueBi):
        if dt not in self.axialReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.nodeIDs[eid] = [nidA,nidB]
        self.bendingMoment1Real[dt][eid] = [bm1Ar,bm1Br]
        self.bendingMoment2Real[dt][eid] = [bm2Ar,bm2Br]
        self.shearPlane1Real[dt][eid] = [sp1Ar,sp1Br]
        self.shearPlane2Real[dt][eid] = [sp2Ar,sp2Br]
        self.axialReal[dt][eid]  = [axialAr,axialBr]
        self.torqueReal[dt][eid] = [torqueAr,torqueBr]

        self.bendingMoment1Imag[dt][eid] = [bm1Ai,bm1Bi]
        self.bendingMoment2Imag[dt][eid] = [bm2Ai,bm2Bi]
        self.shearPlane1Imag[dt][eid] = [sp1Ai,sp1Bi]
        self.shearPlane2Imag[dt][eid] = [sp2Ai,sp2Bi]
        self.axialImag[dt][eid]  = [axialAi,axialBi]
        self.torqueImag[dt][eid] = [torqueAi,torqueBi]

    def __repr__(self):
        return str(self.axialReal)

class ComplexPentaPressureForce(object): # 76-CHEXA_PR,77-PENTA_PR,78-TETRA_PR
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.accelerationReal = {}
        self.accelerationImag = {}
        self.velocityReal = {}
        self.velocityImag = {}
        self.pressure = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.accelerationReal[dt] = {}
        self.accelerationImag[dt] = {}
        self.velocityReal[dt] = {}
        self.velocityImag[dt] = {}
        self.pressure[dt] = {}

    def add(self,dt,data):
        [eid,eName,axr,ayr,azr,vxr,vyr,vzr,pressure,
                   axi,ayi,azi,vxi,vyi,vzi] = data

        #self.eType[eid] = eType
        self.accelerationReal[eid] = [axr,ayr,azr]
        self.accelerationImag[eid] = [axi,ayi,azi]
        self.velocityReal[eid] = [vxr,vyr,vzr]
        self.velocityImag[eid] = [vxi,vyi,vzi]
        self.pressure[eid] = pressure

    def addSort1(self,dt,data):
        [eid,eName,axr,ayr,azr,vxr,vyr,vzr,pressure,
                  axi,ayi,azi,vxi,vyi,vzi] = data
        if dt not in self.accelerationReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.accelerationReal[dt][eid] = [axr,ayr,azr]
        self.accelerationImag[dt][eid] = [axi,ayi,azi]
        self.velocityReal[dt][eid] = [vxr,vyr,vzr]
        self.velocityImag[dt][eid] = [vxi,vyi,vzi]
        self.pressure[dt][eid] = pressure

    def addSort2(self,eid,data):
        [dt,eName,axr,ayr,azr,vxr,vyr,vzr,pressure,
                  axi,ayi,azi,vxi,vyi,vzi] = data
        if dt not in self.accelerationReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.accelerationReal[dt][eid] = [axr,ayr,azr]
        self.accelerationImag[dt][eid] = [axi,ayi,azi]
        self.velocityReal[dt][eid] = [vxr,vyr,vzr]
        self.velocityImag[dt][eid] = [vxi,vyi,vzi]
        self.pressure[dt][eid] = pressure


    def __repr__(self):
        return str(self.accelerationReal)

class ComplexCBUSHForce(object): # 102-CBUSH
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.forceReal = {}
        self.forceImag = {}
        self.momentReal = {}
        self.momentImag = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.forceReal[dt] = {}
        self.forceImag[dt] = {}
        self.momentReal[dt] = {}
        self.momentImag[dt] = {}

    def add(self,dt,data):
        [eid,fxr,fyr,fzr,mxr,myr,mzr,
             fxi,fyi,fzi,mxi,myi,mzi] = data

        #self.eType[eid] = eType
        self.forceReal[eid] = [fxr,fyr,fzr]
        self.forceImag[eid] = [fxi,fyi,fzi]
        self.momentReal[eid] = [mxr,myr,mzr]
        self.momentImag[eid] = [mxi,myi,mzi]

    def addSort1(self,dt,data):
        [eid,fxr,fyr,fzr,mxr,myr,mzr,
             fxi,fyi,fzi,mxi,myi,mzi] = data
        if dt not in self.forceReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.forceReal[dt][eid] = [fxr,fyr,fzr]
        self.forceImag[dt][eid] = [fxi,fyi,fzi]
        self.momentReal[dt][eid] = [mxr,myr,mzr]
        self.momentImag[dt][eid] = [mxi,myi,mzi]

    def addSort2(self,eid,data):
        [dt,fxr,fyr,fzr,mxr,myr,mzr,
            fxi,fyi,fzi,mxi,myi,mzi] = data
        if dt not in self.forceReal:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.forceReal[dt][eid] = [fxr,fyr,fzr]
        self.forceImag[dt][eid] = [fxi,fyi,fzi]
        self.momentReal[dt][eid] = [mxr,myr,mzr]
        self.momentImag[dt][eid] = [mxi,myi,mzi]

    def __repr__(self):
        return str(self.forceReal)

class ComplexForce_VU(object): # 191-VUBEAM
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.parent = {}
        self.coord = {}
        self.icord = {}
        
        self.forceXReal = {}
        self.shearYReal = {}
        self.shearZReal = {}
        self.torsionReal  = {}
        self.bendingYReal = {}
        self.bendingZReal = {}

        self.forceXImag  = {}
        self.shearYImag  = {}
        self.shearZImag  = {}
        self.torsionImag  = {}
        self.bendingYImag = {}
        self.bendingZImag = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.forceXReal[dt] = {}
        self.shearYReal[dt] = {}
        self.shearZReal[dt] = {}
        self.torsionReal[dt]  = {}
        self.bendingYReal[dt] = {}
        self.bendingZReal[dt] = {}

        self.forceXImag[dt]  = {}
        self.shearYImag[dt]  = {}
        self.shearZImag[dt]  = {}
        self.torsionImag[dt]  = {}
        self.bendingYImag[dt] = {}
        self.bendingZImag[dt] = {}

    def add(self,nNodes,dt,data):
        [eid,parent,coord,icord,forces] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType
        
        self.forceXReal[eid]  = {}
        self.shearYReal[eid]  = {}
        self.shearZReal[eid]  = {}
        self.torsionReal[eid]  = {}
        self.bendingYReal[eid] = {}
        self.bendingZReal[eid] = {}

        self.forceXImag[eid]  = {}
        self.shearYImag[eid]  = {}
        self.shearZImag[eid]  = {}
        self.torsionImag[eid]  = {}
        self.bendingYImag[eid] = {}
        self.bendingZImag[eid] = {}

        for force in forces:
            [nid,posit,forceXr,shearYr,shearZr,torsionr,bendingYr,bendingZr,
                       forceXi,shearYi,shearZi,torsioni,bendingYi,bendingZi] = force
            self.forceXReal[eid][nid]  = forceXr
            self.shearYReal[eid][nid]  = shearYr
            self.shearZReal[eid][nid]  = shearZr
            self.torsionReal[eid][nid]  = torsionr
            self.bendingYReal[eid][nid] = bendingYr
            self.bendingZReal[eid][nid] = bendingZr

            self.forceXImag[eid][nid]  = forceXi
            self.shearYImag[eid][nid]  = shearYi
            self.shearZImag[eid][nid]  = shearZi
            self.torsionImag[eid][nid]  = torsioni
            self.bendingYImag[eid][nid] = bendingYi
            self.bendingZImag[eid][nid] = bendingZi

    def addSort1(self,nNodes,dt,data):
        [eid,parent,coord,icord,forces] = data
        if dt not in self.forceXReal:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType
        
        self.forceXReal[dt][eid]  = {}
        self.shearYReal[dt][eid]  = {}
        self.shearZReal[dt][eid]  = {}
        self.torsionReal[dt][eid]  = {}
        self.bendingYReal[dt][eid] = {}
        self.bendingZReal[dt][eid] = {}

        self.forceXImag[dt][eid]  = {}
        self.shearYImag[dt][eid]  = {}
        self.shearZImag[dt][eid]  = {}
        self.torsionImag[dt][eid]  = {}
        self.bendingYImag[dt][eid] = {}
        self.bendingZImag[dt][eid] = {}
        for force in forces:
            [nid,posit,forceXr,shearYr,shearZr,torsionr,bendingYr,bendingZr,
                       forceXi,shearYi,shearZi,torsioni,bendingYi,bendingZi] = force
            self.forceXReal[dt][eid][nid]  = forceXr
            self.shearYReal[dt][eid][nid]  = shearYr
            self.shearZReal[dt][eid][nid]  = shearZr
            self.torsionReal[dt][eid][nid]  = torsionr
            self.bendingYReal[dt][eid][nid] = bendingYr
            self.bendingZReal[dt][eid][nid] = bendingZr

            self.forceXImag[dt][eid][nid]  = forceXi
            self.shearYImag[dt][eid][nid]  = shearYi
            self.shearZImag[dt][eid][nid]  = shearZi
            self.torsionImag[dt][eid][nid]  = torsioni
            self.bendingYImag[dt][eid][nid] = bendingYi
            self.bendingZImag[dt][eid][nid] = bendingZi

    def addSort2(self,nNodes,eid,data):
        [dt,parent,coord,icord,forces] = data
        if dt not in self.forceXReal:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.forceXReal[dt][eid]  = {}
        self.shearYReal[dt][eid]  = {}
        self.shearZReal[dt][eid]  = {}
        self.torsionReal[dt][eid]  = {}
        self.bendingYReal[dt][eid] = {}
        self.bendingZReal[dt][eid] = {}

        self.forceXImag[dt][eid]  = {}
        self.shearYImag[dt][eid]  = {}
        self.shearZImag[dt][eid]  = {}
        self.torsionImag[dt][eid]  = {}
        self.bendingYImag[dt][eid] = {}
        self.bendingZImag[dt][eid] = {}
        for force in forces:
            [nid,posit,forceXr,shearYr,shearZr,torsionr,bendingYr,bendingZr,
                       forceXi,shearYi,shearZi,torsioni,bendingYi,bendingZi] = force
            self.forceXReal[dt][eid][nid]  = forceXr
            self.shearYReal[dt][eid][nid]  = shearYr
            self.shearZReal[dt][eid][nid]  = shearZr
            self.torsionReal[dt][eid][nid]  = torsionr
            self.bendingYReal[dt][eid][nid] = bendingYr
            self.bendingZReal[dt][eid][nid] = bendingZr

            self.forceXImag[dt][eid][nid]  = forceXi
            self.shearYImag[dt][eid][nid]  = shearYi
            self.shearZImag[dt][eid][nid]  = shearZi
            self.torsionImag[dt][eid][nid]  = torsioni
            self.bendingYImag[dt][eid][nid] = bendingYi
            self.bendingZImag[dt][eid][nid] = bendingZi

    def __repr__(self):
        return str(self.forceXReal)

class ComplexForce_VU_2D(object): # 189-VUQUAD,190-VUTRIA
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.parent = {}
        self.coord = {}
        self.icord = {}
        self.theta = {}
        
        self.membraneXReal  = {}
        self.membraneYReal  = {}
        self.membraneXYReal = {}
        self.bendingXReal  = {}
        self.bendingYReal  = {}
        self.bendingXYReal = {}
        self.shearYZReal = {}
        self.shearXZReal = {}

        self.membraneXImag  = {}
        self.membraneYImag  = {}
        self.membraneXYImag = {}
        self.bendingXImag  = {}
        self.bendingYImag  = {}
        self.bendingXYImag = {}
        self.shearYZImag = {}
        self.shearXZImag = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.membraneXReal[dt]  = {}
        self.membraneYReal[dt]  = {}
        self.membraneXYReal[dt] = {}
        self.bendingXReal[dt]  = {}
        self.bendingYReal[dt]  = {}
        self.bendingXYReal[dt] = {}
        self.shearYZReal[dt] = {}
        self.shearXZReal[dt] = {}

        self.membraneXImag[dt]  = {}
        self.membraneYImag[dt]  = {}
        self.membraneXYImag[dt] = {}
        self.bendingXImag[dt]  = {}
        self.bendingYImag[dt]  = {}
        self.bendingXYImag[dt] = {}
        self.shearYZImag[dt] = {}
        self.shearXZImag[dt] = {}

    def add(self,nNodes,dt,data):
        [eid,parent,coord,icord,theta,forces] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType
        
        self.membraneXReal[eid]  = {}
        self.membraneYReal[eid]  = {}
        self.membraneXYReal[eid] = {}
        self.bendingXReal[eid]  = {}
        self.bendingYReal[eid]  = {}
        self.bendingXYReal[eid] = {}
        self.shearYZReal[eid] = {}
        self.shearXZReal[eid] = {}

        self.membraneXImag[eid]  = {}
        self.membraneYImag[eid]  = {}
        self.membraneXYImag[eid] = {}
        self.bendingXImag[eid]  = {}
        self.bendingYImag[eid]  = {}
        self.bendingXYImag[eid] = {}
        self.shearYZImag[eid] = {}
        self.shearXZImag[eid] = {}

        for force in forces:
            [nid,membraneXr,membraneYr,membraneXYr,_,_,_,bendingXr,bendingYr,bendingXYr,shearYZr,shearXZr,_,
                 membraneXi,membraneYi,membraneXYi,_,_,_,bendingXi,bendingYi,bendingXYi,shearYZi,shearXZi,_] = force
            self.membraneXReal[dt][eid][nid]  = membraneXr
            self.membraneYReal[dt][eid][nid]  = membraneYr
            self.membraneXYReal[dt][eid][nid] = membraneXYr
            self.bendingXReal[dt][eid][nid]  = bendingXr
            self.bendingYReal[dt][eid][nid]  = bendingYr
            self.bendingXYReal[dt][eid][nid] = bendingXYr
            self.shearYZReal[dt][eid][nid]  = shearYZr
            self.shearXZReal[dt][eid][nid]  = shearXZr

            self.membraneXImag[eid][nid]  = membraneXi
            self.membraneYImag[eid][nid]  = membraneYi
            self.membraneXYImag[eid][nid] = membraneXYi
            self.bendingXImag[eid][nid]  = bendingXi
            self.bendingYImag[eid][nid]  = bendingYi
            self.bendingXYImag[eid][nid] = bendingXYi
            self.shearYZImag[eid][nid]  = shearYZi
            self.shearXZImag[eid][nid]  = shearXZi

    def addSort1(self,nNodes,dt,data):
        [eid,parent,coord,icord,theta,forces] = data
        self._fillObject(dt,eid,parent,coord,icord,theta,forces)

    def addSort2(self,nNodes,eid,data):
        [dt,parent,coord,icord,theta,forces] = data
        self._fillObject(dt,eid,parent,coord,icord,theta,forces)

    def _fillObject(self,dt,eid,parent,coord,icord,theta,forces):
        if dt not in self.membraneXReal:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.membraneXReal[dt][eid]  = {}
        self.membraneYReal[dt][eid]  = {}
        self.membraneXYReal[dt][eid] = {}
        self.bendingXReal[dt][eid]  = {}
        self.bendingYReal[dt][eid]  = {}
        self.bendingXYReal[dt][eid] = {}
        self.shearYZReal[dt][eid] = {}
        self.shearXZReal[dt][eid] = {}

        self.membraneXImag[dt][eid]  = {}
        self.membraneYImag[dt][eid]  = {}
        self.membraneXYImag[dt][eid] = {}
        self.bendingXImag[dt][eid]  = {}
        self.bendingYImag[dt][eid]  = {}
        self.bendingXYImag[dt][eid] = {}
        self.shearYZImag[dt][eid] = {}
        self.shearXZImag[dt][eid] = {}
        for force in forces:
            [nid,membraneXr,membraneYr,membraneXYr,_,_,_,bendingXr,bendingYr,bendingXYr,shearYZr,shearXZr,_,
                 membraneXi,membraneYi,membraneXYi,_,_,_,bendingXi,bendingYi,bendingXYi,shearYZi,shearXZi,_] = force
            self.membraneXReal[dt][eid][nid]  = membraneXr
            self.membraneYReal[dt][eid][nid]  = membraneYr
            self.membraneXYReal[dt][eid][nid] = membraneXYr
            self.bendingXReal[dt][eid][nid]  = bendingXr
            self.bendingYReal[dt][eid][nid]  = bendingYr
            self.bendingXYReal[dt][eid][nid] = bendingXYr
            self.shearYZReal[dt][eid][nid]  = shearYZr
            self.shearXZReal[dt][eid][nid]  = shearXZr

            self.membraneXImag[dt][eid][nid]  = membraneXi
            self.membraneYImag[dt][eid][nid]  = membraneYi
            self.membraneXYImag[dt][eid][nid] = membraneXYi
            self.bendingXImag[dt][eid][nid]  = bendingXi
            self.bendingYImag[dt][eid][nid]  = bendingYi
            self.bendingXYImag[dt][eid][nid] = bendingXYi
            self.shearYZImag[dt][eid][nid]  = shearYZi
            self.shearXZImag[dt][eid][nid]  = shearXZi

    def __repr__(self):
        return str(self.membraneXReal)
