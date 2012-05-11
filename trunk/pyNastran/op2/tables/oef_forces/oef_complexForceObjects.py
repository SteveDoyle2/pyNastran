from pyNastran.op2.resultObjects.op2_Objects import scalarObject

class ComplexRodForce(scalarObject): # 1-ROD, 3-TUBE, 10-CONROD
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.axialForce = {}
        self.torque = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.axialForce[dt] = {}
        self.torque[dt] = {}

    def add(self,dt,data):
        [eid,axialForceReal,torqueReal,axialForceImag,torqueImag] = data

        #self.eType[eid] = eType
        self.axialForce[eid] = complex(axialForceReal,axialForceImag)
        self.torque[eid]     = complex(torqueReal,torqueImag)

    def addSort1(self,dt,data):
        [eid,axialForceReal,torqueReal,axialForceImag,torqueImag] = data
        if dt not in self.axialForce:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = complex(axialForceReal,axialForceImag)
        self.torque[dt][eid]     = complex(torqueReal,torqueImag)

    def addSort2(self,eid,data):
        [dt,axialForceReal,torqueReal,axialForceImag,torqueImag] = data
        if dt not in self.axialForce:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = complex(axialForceReal,axialForceImag)
        self.torque[dt][eid]     = complex(torqueReal,torqueImag)

    def __repr__(self):
        return str(self.axialForce)

class ComplexCBEAMForce(scalarObject): # 2-CBEAM
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.bendingMoment = {}
        self.shear = {}
        self.axial = {}
        self.totalTorque = {}
        self.warpingTorque = {}

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
        self.bendingMoment[dt] = {}
        self.shear[dt] = {}
        self.axial[dt] = {}
        self.totalTorque[dt] = {}
        self.warpingTorque[dt] = {}

    def addNewElement(self,dt,data):
        [eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                    bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi] = data
        print "CBEAM addnew",data
        #self.eType[eid] = eType
        self.bendingMoment[eid] = {sd:[complex(bm1r,bm1i),complex(bm2r,bm2i)]}
        self.shear[eid] = {sd:[complex(ts1r,ts1i),complex(ts2r,ts2i)]}
        self.axial[eid] = {sd:complex(afr,afi)}
        self.totalTorque[eid] = {sd:complex(ttrqr,ttrqi)}
        self.warpingTorque[eid] = {sd:complex(wtrqr,wtrqi)}

    def add(self,dt,data):
        [eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                    bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi] = data
        print "CBEAM add   ",data

        #self.eType[eid] = eType
        self.bendingMoment[eid][sd] = [complex(bm1r,bm1i),complex(bm2r,bm2i)]
        self.shear[eid][sd] = [complex(ts1r,ts1i),complex(ts2r,ts2i)]
        self.axial[eid][sd] = complex(afr,afi)
        self.totalTorque[eid][sd]   = complex(ttrqr,ttrqi)
        self.warpingTorque[eid][sd] = complex(wtrqr,wtrqi)

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
        if dt not in self.axial:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMoment[dt][eid] = {sd:[complex(bm1r,bm1i),complex(bm2r,bm2i)]}
        self.shear[dt][eid] = {sd:[complex(ts1r,ts1i),complex(ts2r,ts2i)]}
        self.axial[dt][eid] = {sd:complex(afr,afi)}
        self.totalTorque[dt][eid] = {sd:complex(ttrqr,ttrqi)}
        self.warpingTorque[dt][eid] = {sd:complex(wtrqr,wtrqi)}

    def _fillObject(self,dt,eid,nid,sd,bm1r,bm2r,ts1r,ts2r,afr,ttrqr,wtrqr,
                                       bm1i,bm2i,ts1i,ts2i,afi,ttrqi,wtrqi):
        #if dt not in self.axial:
            #self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMoment[dt][eid][sd] = [complex(bm1r,bm1i),complex(bm2r,bm2i)]
        self.shear[dt][eid][sd] = [complex(ts1r,ts1i),complex(ts2r,ts2i)]
        self.axial[dt][eid][sd] = complex(afr,afi)
        self.totalTorque[dt][eid][sd]   = complex(ttrqr,ttrqi)
        self.warpingTorque[dt][eid][sd] = complex(wtrqr,wtrqi)

    def __repr__(self):
        return str(self.axial)

class ComplexCShearForce(scalarObject): # 4-CSHEAR
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.force41 = {}
        self.force14 = {}
        self.force21 = {}
        self.force12 = {}
        self.force32 = {}
        self.force23 = {}
        self.force43 = {}
        self.force34 = {}
        self.kickForce1 = {}
        self.kickForce2 = {}
        self.kickForce3 = {}
        self.kickForce4 = {}
        self.shear12 = {}
        self.shear23 = {}
        self.shear34 = {}
        self.shear41 = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.force41[dt] = {}
        self.force14[dt] = {}
        self.force21[dt] = {}
        self.force12[dt] = {}
        self.force32[dt] = {}
        self.force23[dt] = {}
        self.force43[dt] = {}
        self.force34[dt] = {}
        self.kickForce1[dt] = {}
        self.kickForce2[dt] = {}
        self.kickForce3[dt] = {}
        self.kickForce4[dt] = {}
        self.shear12[dt] = {}
        self.shear23[dt] = {}
        self.shear34[dt] = {}
        self.shear41[dt] = {}


    def add(self,dt,data):
        [eid,f41r,f21r,f12r,f32r,f23r,f43r,f34r,f14r,
             f41i,f21i,f12i,f32i,f23i,f43i,f34i,f14i,
             kf1r,s12r,kf2r,s23r,kf3r,s34r,kf4r,s41r,
             kf1i,s12i,kf2i,s23i,kf3i,s34i,kf4i,s41i] = data
        #self.eType[eid] = eType
        self.force41[eid] = complex(f41r,f41i)
        self.force14[eid] = complex(f14r,f14i)
        self.force21[eid] = complex(f21r,f21i)
        self.force12[eid] = complex(f12r,f12i)
        self.force32[eid] = complex(f32r,f32i)
        self.force23[eid] = complex(f23r,f23i)
        self.force43[eid] = complex(f43r,f43i)
        self.force34[eid] = complex(f34r,f34i)
        self.kickForce1[eid] = complex(kf1r,kf1i)
        self.kickForce2[eid] = complex(kf2r,kf2i)
        self.kickForce3[eid] = complex(kf3r,kf3i)
        self.kickForce4[eid] = complex(kf4r,kf4i)
        self.shear12[eid] = complex(s12r,s12i)
        self.shear23[eid] = complex(s23r,s23i)
        self.shear34[eid] = complex(s34r,s34i)
        self.shear41[eid] = complex(s41r,s41i)
        
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

        if dt not in self.force41:
            self.addNewTransient(dt)
        #self.eType[eid] = eType
        self.force41[dt][eid] = complex(f41r,f41i)
        self.force14[dt][eid] = complex(f14r,f14i)
        self.force21[dt][eid] = complex(f21r,f21i)
        self.force12[dt][eid] = complex(f12r,f12i)
        self.force32[dt][eid] = complex(f32r,f32i)
        self.force23[dt][eid] = complex(f23r,f23i)
        self.force43[dt][eid] = complex(f43r,f43i)
        self.force34[dt][eid] = complex(f34r,f34i)
        self.kickForce1[dt][eid] = complex(kf1r,kf1i)
        self.kickForce2[dt][eid] = complex(kf2r,kf2i)
        self.kickForce3[dt][eid] = complex(kf3r,kf3i)
        self.kickForce4[dt][eid] = complex(kf4r,kf4i)
        self.shear12[dt][eid] = complex(s12r,s12i)
        self.shear23[dt][eid] = complex(s23r,s23i)
        self.shear34[dt][eid] = complex(s34r,s34i)
        self.shear41[dt][eid] = complex(s41r,s41i)

    def __repr__(self):
        return str(self.force41)

class ComplexSpringForce(scalarObject): # 11-CELAS1,12-CELAS2,13-CELAS3, 14-CELAS4
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
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
        self.force[eid] = complex(forceReal,forceImag)

    def addSort1(self,dt,data):
        [eid,forceReal,forceImag] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = complex(forceReal,forceImag)
    def addSort2(self,eid,data):
        [dt,forceReal,forceImag] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = complex(forceReal,forceImag)

    def __repr__(self):
        return str(self.force)

class ComplexViscForce(scalarObject): # 24-CVISC
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.axialForce = {}
        self.torque = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.axialForce[dt] = {}
        self.torque[dt] = {}

    def add(self,dt,data):
        [eid,axialForceReal,torqueReal,axialForceImag,torqueImag] = data

        #self.eType[eid] = eType
        self.axialForce[eid] = complex(axialForceReal,axialForceImag)
        self.torque[eid]     = complex(torqueReal,torqueImag)

    def addSort1(self,dt,data):
        [eid,axialForceReal,torqueReal,axialForceImag,torqueImag] = data
        if dt not in self.axialForce:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = complex(axialForceReal,axialForceImag)
        self.torque[dt][eid]     = complex(torqueReal,torqueImag)

    def addSort2(self,eid,data):
        [dt,axialForceReal,torqueReal,axialForceImag,torqueImag] = data
        if dt not in self.axialForce:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = complex(axialForceReal,axialForceImag)
        self.torque[dt][eid]     = complex(torqueReal,torqueImag)

    def __repr__(self):
        return str(self.axialForce)

class ComplexPlateForce(scalarObject): # 33-CQUAD4, 74-CTRIA3
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
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
        self.mx[eid] = complex(mxr,mxi)
        self.my[eid] = complex(myr,myi)
        self.mxy[eid] = complex(mxyr,mxyi)
        self.bmx[eid] = complex(bmxr,bmxi)
        self.bmy[eid] = complex(bmyr,bmyi)
        self.bmxy[eid] = complex(bmxyr,bmxyi)
        self.tx[eid] = complex(txr,txi)
        self.ty[eid] = complex(tyr,tyi)

    def addSort1(self,dt,data):
        [eid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
             mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid] = complex(mxr,mxi)
        self.my[dt][eid] = complex(myr,myi)
        self.mxy[dt][eid] = complex(mxyr,mxyi)
        self.bmx[dt][eid] = complex(bmxr,bmxi)
        self.bmy[dt][eid] = complex(bmyr,bmyi)
        self.bmxy[dt][eid] = complex(bmxyr,bmxyi)
        self.tx[dt][eid] = complex(txr,txi)
        self.ty[dt][eid] = complex(tyr,tyi)

    def addSort2(self,eid,data):
        [dt,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
            mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid] = complex(mxr,mxi)
        self.my[dt][eid] = complex(myr,myi)
        self.mxy[dt][eid] = complex(mxyr,mxyi)
        self.bmx[dt][eid] = complex(bmxr,bmxi)
        self.bmy[dt][eid] = complex(bmyr,bmyi)
        self.bmxy[dt][eid] = complex(bmxyr,bmxyi)
        self.tx[dt][eid] = complex(txr,txi)
        self.ty[dt][eid] = complex(tyr,tyi)

    def __repr__(self):
        return str(self.mx)

class ComplexPLATE2Force(scalarObject): # 64-CQUAD8, 75-CTRIA6, 82-CQUADR
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.term = {}
        self.ngrids = {}
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

        self.mx[dt] = {}
        self.my[dt] = {}
        self.mxy[dt] = {}
        self.bmx[dt] = {}
        self.bmy[dt] = {}
        self.bmxy[dt] = {}
        self.tx[dt] = {}
        self.ty[dt] = {}

    def addNewElement(self,eid,dt,data):
        #print "eid = ",eid
        [term,nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                  mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data

        #self.eType[eid] = eType
        self.term[eid] = term
        self.ngrids[eid] = nid

        self.mx[eid] = [complex(mxr,mxi)]
        self.my[eid] = [complex(myr,myi)]
        self.mxy[eid] = [complex(mxyr,mxyi)]
        self.bmx[eid] = [complex(bmxr,bmxi)]
        self.bmy[eid] = [complex(bmyr,bmyi)]
        self.bmxy[eid] = [complex(bmxyr,bmxyi)]
        self.tx[eid] = [complex(txr,txi)]
        self.ty[eid] = [complex(tyr,tyi)]

    def add(self,eid,dt,data):
        [nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
             mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data

        #self.eType[eid] = eType
        #print "mx = ",self.mx,mx
        self.mx[eid].append(complex(mxr,mxi))
        self.my[eid].append(complex(myr,myi))
        self.mxy[eid].append(complex(mxyr,mxyi))
        self.bmx[eid].append(complex(bmxr,bmxi))
        self.bmy[eid].append(complex(bmyr,bmyi))
        self.bmxy[eid].append(complex(bmxyr,bmxyi))
        self.tx[eid].append(complex(txr,txi))
        self.ty[eid].append(complex(tyr,tyi))

    def addNewElementSort1(self,eid,dt,data):
        [term,nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                  mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.term[dt][eid] = term
        self.ngrids[dt][eid] = nid
        self.mx[dt][eid] = [complex(mxr,mxi)]
        self.my[dt][eid] = [complex(myr,myi)]
        self.mxy[dt][eid] = [complex(mxyr,mxyi)]
        self.bmx[dt][eid] = [complex(bmxr,bmxi)]
        self.bmy[dt][eid] = [complex(bmyr,bmyi)]
        self.bmxy[dt][eid] = [complex(bmxyr,bmxyi)]
        self.tx[dt][eid] = [complex(txr,txi)]
        self.ty[dt][eid] = [complex(tyr,tyi)]

    def addSort1(self,eid,dt,data):
        [nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
             mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid].append(complex(mxr,mxi))
        self.my[dt][eid].append(complex(myr,myi))
        self.mxy[dt][eid].append(complex(mxyr,mxyi))
        self.bmx[dt][eid].append(complex(bmxr,bmxi))
        self.bmy[dt][eid].append(complex(bmyr,bmyi))
        self.bmxy[dt][eid].append(complex(bmxyr,bmxyi))
        self.tx[dt][eid].append(complex(txr,txi))
        self.ty[dt][eid].append(complex(tyr,tyi))

    def addNewElementSort2(self,dt,eid,data):
        [term,nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
                  mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.term[dt][eid] = term
        self.ngrids[dt][eid] = nid

        self.mx[dt][eid] = [complex(mxr,mxi)]
        self.my[dt][eid] = [complex(myr,myi)]
        self.mxy[dt][eid] = [complex(mxyr,mxyi)]
        self.bmx[dt][eid] = [complex(bmxr,bmxi)]
        self.bmy[dt][eid] = [complex(bmyr,bmyi)]
        self.bmxy[dt][eid] = [complex(bmxyr,bmxyi)]
        self.tx[dt][eid] = [complex(txr,txi)]
        self.ty[dt][eid] = [complex(tyr,tyi)]

    def addSort2(self,dt,eid,data):
        [nid,mxr,myr,mxyr,bmxr,bmyr,bmxyr,txr,tyr,
             mxi,myi,mxyi,bmxi,bmyi,bmxyi,txi,tyi] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid].append(complex(mxr,mxi))
        self.my[dt][eid].append(complex(myr,myi))
        self.mxy[dt][eid].append(complex(mxyr,mxyi))
        self.bmx[dt][eid].append(complex(bmxr,bmxi))
        self.bmy[dt][eid].append(complex(bmyr,bmyi))
        self.bmxy[dt][eid].append(complex(bmxyr,bmxyi))
        self.tx[dt][eid].append(complex(txr,txi))
        self.ty[dt][eid].append(complex(tyr,tyi))

    def __repr__(self):
        return str(self.mx)


class ComplexCBARForce(scalarObject): # 34-CBAR
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.bendingMomentA = {}
        self.bendingMomentB = {}
        self.shear = {}
        self.axial = {}
        self.torque = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.bendingMomentA[dt] = {}
        self.bendingMomentB[dt] = {}
        self.shear[dt] = {}
        self.axial[dt] = {}
        self.torque[dt] = {}

    def add(self,dt,data):
        [eid,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
             bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi] = data

        #self.eType[eid] = eType
        self.bendingMomentA[eid] = [complex(bm1ar,bm1ai),complex(bm2ar,bm2ai)]
        self.bendingMomentB[eid] = [complex(bm1br,bm1bi),complex(bm2br,bm2bi)]
        self.shear[eid] = [complex(ts1r,ts1i),complex(ts2r,ts2i)]
        self.axial[eid] = complex(afr,afi)
        self.torque[eid] = complex(trqr,trqi)

    def addSort1(self,dt,data):
        [eid,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
             bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi] = data
        if dt not in self.axial:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMomentA[dt][eid] = [complex(bm1ar,bm1ai),complex(bm2ar,bm2ai)]
        self.bendingMomentB[dt][eid] = [complex(bm1br,bm1bi),complex(bm2br,bm2bi)]
        self.shear[dt][eid] = [complex(ts1r,ts1i),complex(ts2r,ts2i)]
        self.axial[dt][eid] = complex(afr,afi)
        self.torque[dt][eid] = complex(trqr,trqi)

    def addSort2(self,eid,data):
        [dt,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
            bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi] = data
        if dt not in self.axial:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMomentA[dt][eid] = [complex(bm1ar,bm1ai),complex(bm2ar,bm2ai)]
        self.bendingMomentB[dt][eid] = [complex(bm1br,bm1bi),complex(bm2br,bm2bi)]
        self.shear[dt][eid] = [complex(ts1r,ts1i),complex(ts2r,ts2i)]
        self.axial[dt][eid] = complex(afr,afi)
        self.torque[dt][eid] = complex(trqr,trqi)

    def __repr__(self):
        return str(self.axial)

class ComplexBendForce(scalarObject): # 69-CBEND
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.nodeIDs = {}
        self.bendingMoment1 = {}
        self.bendingMoment2 = {}
        self.shearPlane1 = {}
        self.shearPlane2 = {}
        self.axial  = {}
        self.torque = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.bendingMoment1[dt] = {}
        self.bendingMoment2[dt] = {}
        self.shearPlane1[dt] = {}
        self.shearPlane2[dt] = {}
        self.axial[dt]  = {}
        self.torque[dt] = {}

    def add(self,dt,data):
        [eid,nidA,bm1Ar,bm2Ar,sp1Ar,sp2Ar,axialAr,torqueAr,
                  bm1Ai,bm2Ai,sp1Ai,sp2Ai,axialAi,torqueAi,
             nidB,bm1Br,bm2Br,sp1Br,sp2Br,axialBr,torqueBr,
                  bm1Bi,bm2Bi,sp1Bi,sp2Bi,axialBi,torqueBi] = data

        #self.eType[eid] = eType
        self.nodeIDs[eid] = [nidA,nidB]
        self.bendingMoment1[eid] = [complex(bm1Ar,bm1Ai),complex(bm1Br,bm1Bi)]
        self.bendingMoment2[eid] = [complex(bm2Ar,bm2Ai),complex(bm2Br,bm2Bi)]
        self.shearPlane1[eid] = [complex(sp1Ar,sp1Ai),complex(sp1Br,sp1Bi)]
        self.shearPlane2[eid] = [complex(sp2Ar,sp2Ai),complex(sp2Br,sp2Bi)]
        self.axial[eid]  = [complex(axialAr,axialAi),complex(axialBr,axialBi)]
        self.torque[eid] = [complex(torqueAr,torqueAi),complex(torqueBr,torqueBi)]

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
        if dt not in self.axial:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.nodeIDs[eid] = [nidA,nidB]
        self.bendingMoment1[dt][eid] = [complex(bm1Ar,bm1Ai),complex(bm1Br,bm1Bi)]
        self.bendingMoment2[dt][eid] = [complex(bm2Ar,bm2Ai),complex(bm2Br,bm2Bi)]
        self.shearPlane1[dt][eid] = [complex(sp1Ar,sp1Ai),complex(sp1Br,sp1Bi)]
        self.shearPlane2[dt][eid] = [complex(sp2Ar,sp2Ai),complex(sp2Br,sp2Bi)]
        self.axial[dt][eid]  = [complex(axialAr,axialAi),complex(axialBr,axialBi)]
        self.torque[dt][eid] = [complex(torqueAr,torqueAi),complex(torqueBr,torqueBi)]

    def __repr__(self):
        return str(self.axial)

class ComplexPentaPressureForce(scalarObject): # 76-CHEXA_PR,77-PENTA_PR,78-TETRA_PR
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.acceleration = {}
        self.velocity = {}
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
        self.acceleration[dt] = {}
        self.velocity[dt] = {}
        self.pressure[dt] = {}

    def add(self,dt,data):
        [eid,eName,axr,ayr,azr,vxr,vyr,vzr,pressure,
                   axi,ayi,azi,vxi,vyi,vzi] = data

        #self.eType[eid] = eType
        self.acceleration[eid] = [complex(axr,axi),complex(ayr,ayi),complex(azr,azi)]
        self.velocity[eid]     = [complex(vxr,vxi),complex(vyr,vyi),complex(vzr,vzi)]
        self.pressure[eid] = pressure

    def addSort1(self,dt,data):
        [eid,eName,axr,ayr,azr,vxr,vyr,vzr,pressure,
                  axi,ayi,azi,vxi,vyi,vzi] = data
        if dt not in self.acceleration:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.acceleration[dt][eid] = [complex(axr,axi),complex(ayr,ayi),complex(azr,azi)]
        self.velocity[dt][eid]     = [complex(vxr,vxi),complex(vyr,vyi),complex(vzr,vzi)]
        self.pressure[dt][eid] = pressure

    def addSort2(self,eid,data):
        [dt,eName,axr,ayr,azr,vxr,vyr,vzr,pressure,
                  axi,ayi,azi,vxi,vyi,vzi] = data
        if dt not in self.acceleration:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.acceleration[dt][eid] = [complex(axr,axi),complex(ayr,ayi),complex(azr,azi)]
        self.velocity[dt][eid]     = [complex(vxr,vxi),complex(vyr,vyi),complex(vzr,vzi)]
        self.pressure[dt][eid] = pressure


    def __repr__(self):
        return str(self.acceleration)

class ComplexCBUSHForce(scalarObject): # 102-CBUSH
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.force = {}
        self.moment = {}

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
        self.moment[dt] = {}

    def add(self,dt,data):
        [eid,fxr,fyr,fzr,mxr,myr,mzr,
             fxi,fyi,fzi,mxi,myi,mzi] = data

        #self.eType[eid] = eType
        self.force[eid]  = [complex(fxr,fxi),complex(fyr,fyi),complex(fzr,fzi)]
        self.moment[eid] = [complex(mxr,mxi),complex(myr,myi),complex(mzr,mzi)]

    def addSort1(self,dt,data):
        [eid,fxr,fyr,fzr,mxr,myr,mzr,
             fxi,fyi,fzi,mxi,myi,mzi] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid]  = [complex(fxr,fxi),complex(fyr,fyi),complex(fzr,fzi)]
        self.moment[dt][eid] = [complex(mxr,mxi),complex(myr,myi),complex(mzr,mzi)]

    def addSort2(self,eid,data):
        [dt,fxr,fyr,fzr,mxr,myr,mzr,
            fxi,fyi,fzi,mxi,myi,mzi] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid]  = [complex(fxr,fxi),complex(fyr,fyi),complex(fzr,fzi)]
        self.moment[dt][eid] = [complex(mxr,mxi),complex(myr,myi),complex(mzr,mzi)]

    def __repr__(self):
        return str(self.force)

class ComplexForce_VU(scalarObject): # 191-VUBEAM
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.parent = {}
        self.coord = {}
        self.icord = {}
        
        self.forceX = {}
        self.shearY = {}
        self.shearZ = {}
        self.torsion  = {}
        self.bendingY = {}
        self.bendingZ = {}

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
        self.forceX[dt] = {}
        self.shearY[dt] = {}
        self.shearZ[dt] = {}
        self.torsion[dt]  = {}
        self.bendingY[dt] = {}
        self.bendingZ[dt] = {}

    def add(self,nNodes,dt,data):
        [eid,parent,coord,icord,forces] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType
        
        self.forceX[eid]  = {}
        self.shearY[eid]  = {}
        self.shearZ[eid]  = {}
        self.torsion[eid]  = {}
        self.bendingY[eid] = {}
        self.bendingZ[eid] = {}

        for force in forces:
            [nid,posit,forceXr,shearYr,shearZr,torsionr,bendingYr,bendingZr,
                       forceXi,shearYi,shearZi,torsioni,bendingYi,bendingZi] = force
            self.forceX[eid][nid]   = complex(forceXr,forceXi)
            self.shearY[eid][nid]   = complex(shearYr,shearYi)
            self.shearZ[eid][nid]   = complex(shearZr,shearZi)
            self.torsion[eid][nid]  = complex(torsionr, torsioni)
            self.bendingY[eid][nid] = complex(bendingYr,bendingYi)
            self.bendingZ[eid][nid] = complex(bendingZr,bendingZi)

    def addSort1(self,nNodes,dt,data):
        [eid,parent,coord,icord,forces] = data
        if dt not in self.forceX:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType
        
        self.forceX[dt][eid]  = {}
        self.shearY[dt][eid]  = {}
        self.shearZ[dt][eid]  = {}
        self.torsion[dt][eid]  = {}
        self.bendingY[dt][eid] = {}
        self.bendingZ[dt][eid] = {}

        for force in forces:
            [nid,posit,forceXr,shearYr,shearZr,torsionr,bendingYr,bendingZr,
                       forceXi,shearYi,shearZi,torsioni,bendingYi,bendingZi] = force
            self.forceX[dt][eid][nid]   = complex(forceXr,forceXi)
            self.shearY[dt][eid][nid]   = complex(shearYr,shearYi)
            self.shearZ[dt][eid][nid]   = complex(shearZr,shearZi)
            self.torsion[dt][eid][nid]  = complex(torsionr, torsioni)
            self.bendingY[dt][eid][nid] = complex(bendingYr,bendingYi)
            self.bendingZ[dt][eid][nid] = complex(bendingZr,bendingZi)

    def addSort2(self,nNodes,eid,data):
        [dt,parent,coord,icord,forces] = data
        if dt not in self.forceX:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.forceX[dt][eid]  = {}
        self.shearY[dt][eid]  = {}
        self.shearZ[dt][eid]  = {}
        self.torsion[dt][eid]  = {}
        self.bendingY[dt][eid] = {}
        self.bendingZ[dt][eid] = {}
        for force in forces:
            [nid,posit,forceXr,shearYr,shearZr,torsionr,bendingYr,bendingZr,
                       forceXi,shearYi,shearZi,torsioni,bendingYi,bendingZi] = force
            self.forceX[dt][eid][nid]   = complex(forceXr,forceXi)
            self.shearY[dt][eid][nid]   = complex(shearYr,shearYi)
            self.shearZ[dt][eid][nid]   = complex(shearZr,shearZi)
            self.torsion[dt][eid][nid]  = complex(torsionr, torsioni)
            self.bendingY[dt][eid][nid] = complex(bendingYr,bendingYi)
            self.bendingZ[dt][eid][nid] = complex(bendingZr,bendingZi)

    def __repr__(self):
        return str(self.forceX)

class ComplexForce_VU_2D(scalarObject): # 189-VUQUAD,190-VUTRIA
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.parent = {}
        self.coord = {}
        self.icord = {}
        self.theta = {}
        
        self.membraneX  = {}
        self.membraneY  = {}
        self.membraneXY = {}
        self.bendingX  = {}
        self.bendingY  = {}
        self.bendingXY = {}
        self.shearYZ = {}
        self.shearXZ = {}

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
        self.membraneX[dt]  = {}
        self.membraneY[dt]  = {}
        self.membraneXY[dt] = {}
        self.bendingX[dt]  = {}
        self.bendingY[dt]  = {}
        self.bendingXY[dt] = {}
        self.shearYZ[dt] = {}
        self.shearXZ[dt] = {}

    def add(self,nNodes,dt,data):
        [eid,parent,coord,icord,theta,forces] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType
        
        self.membraneX[eid]  = {}
        self.membraneY[eid]  = {}
        self.membraneXY[eid] = {}
        self.bendingX[eid]  = {}
        self.bendingY[eid]  = {}
        self.bendingXY[eid] = {}
        self.shearYZ[eid] = {}
        self.shearXZ[eid] = {}

        for force in forces:
            [nid,membraneXr,membraneYr,membraneXYr,bendingXr,bendingYr,bendingXYr,shearYZr,shearXZr,
                 membraneXi,membraneYi,membraneXYi,bendingXi,bendingYi,bendingXYi,shearYZi,shearXZi] = force
            self.membraneX[eid][nid]  = complex(membraneXr, membraneXi)
            self.membraneY[eid][nid]  = complex(membraneYr, membraneYi)
            self.membraneXY[eid][nid] = complex(membraneXYr,membraneXYi)
            self.bendingX[eid][nid]   = complex(bendingXr,  bendingXi)
            self.bendingY[eid][nid]   = complex(bendingYr,  bendingYi)
            self.bendingXY[eid][nid]  = complex(bendingXYr, bendingXYi)
            self.shearYZ[eid][nid]  = complex(shearYZr,shearYZi)
            self.shearXZ[eid][nid]  = complex(shearXZr,shearXZi)

    def addSort1(self,nNodes,dt,data):
        [eid,parent,coord,icord,theta,forces] = data
        self._fillObject(dt,eid,parent,coord,icord,theta,forces)

    def addSort2(self,nNodes,eid,data):
        [dt,parent,coord,icord,theta,forces] = data
        self._fillObject(dt,eid,parent,coord,icord,theta,forces)

    def _fillObject(self,dt,eid,parent,coord,icord,theta,forces):
        if dt not in self.membraneX:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.membraneX[dt][eid]  = {}
        self.membraneY[dt][eid]  = {}
        self.membraneXY[dt][eid] = {}
        self.bendingX[dt][eid]  = {}
        self.bendingY[dt][eid]  = {}
        self.bendingXY[dt][eid] = {}
        self.shearYZ[dt][eid] = {}
        self.shearXZ[dt][eid] = {}

        for force in forces:
            [nid,membraneXr,membraneYr,membraneXYr,bendingXr,bendingYr,bendingXYr,shearYZr,shearXZr,
                 membraneXi,membraneYi,membraneXYi,bendingXi,bendingYi,bendingXYi,shearYZi,shearXZi] = force
            self.membraneX[dt][eid][nid]  = complex(membraneXr, membraneXi)
            self.membraneY[dt][eid][nid]  = complex(membraneYr, membraneYi)
            self.membraneXY[dt][eid][nid] = complex(membraneXYr,membraneXYi)
            self.bendingX[dt][eid][nid]   = complex(bendingXr,  bendingXi)
            self.bendingY[dt][eid][nid]   = complex(bendingYr,  bendingYi)
            self.bendingXY[dt][eid][nid]  = complex(bendingXYr, bendingXYi)
            self.shearYZ[dt][eid][nid]  = complex(shearYZr,shearYZi)
            self.shearXZ[dt][eid][nid]  = complex(shearXZr,shearXZi)

    def __repr__(self):
        return str(self.membraneX)
