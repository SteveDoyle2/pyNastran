from pyNastran.op2.resultObjects.op2_Objects import scalarObject

class ComplexRodForce(scalarObject): # 1-ROD, 3-TUBE, 10-CONROD
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.axialForce = {}
        self.torque = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.dt = dt
        self.axialForce[dt] = {}
        self.torque[dt] = {}

    def add(self,dt,data):
        [eid,axialForce,torque] = data

        #self.eType[eid] = eType
        self.axialForce[eid] = axialForce
        self.torque[eid]     = torque

    def addSort1(self,dt,data):
        [eid,axialForce,torque] = data
        if dt not in self.axialForce:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid]     = torque

    def addSort2(self,eid,data):
        [dt,axialForce,torque] = data
        if dt not in self.axialForce:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid]     = torque

    def __repr__(self):
        return str(self.axialForce)

class ComplexCBeamForce(scalarObject): # 2-CBEAM
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.bendingMoment = {}
        self.shear = {}
        self.axial = {}
        self.totalTorque = {}
        self.warpingTorque = {}

        self.dt = dt
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
        self.dt = dt
        self.bendingMoment[dt] = {}
        self.shear[dt] = {}
        self.axial[dt] = {}
        self.totalTorque[dt] = {}
        self.warpingTorque[dt] = {}

    def addNewElement(self,dt,data):
        [eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        #print "CBEAM addnew",data
        #self.eType[eid] = eType
        self.bendingMoment[eid] = {sd:[bm1,bm2]}
        self.shear[eid] = {sd:[ts1,ts2]}
        self.axial[eid] = {sd:af}
        self.totalTorque[eid] = {sd:ttrq}
        self.warpingTorque[eid] = {sd:wtrq}

    def add(self,dt,data):
        [eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        #print "CBEAM add   ",data

        #self.eType[eid] = eType
        self.bendingMoment[eid][sd] = [bm1,bm2]
        self.shear[eid][sd] = [ts1,ts2]
        self.axial[eid][sd] = af
        self.totalTorque[eid][sd]   = ttrq
        self.warpingTorque[eid][sd] = wtrq

    def addNewElementSort1(self,dt,data):
        [eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data

        self._fillNewObject(dt,eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq)

    def addSort1(self,dt,data):
        [eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        self._fillObject(dt,eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq)

    def addNewElementSort2(self,eid,data):
        [dt,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        self._fillNewObject(dt,eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq)

    def addSort2(self,eid,data):
        [dt,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        self._fillObject(dt,eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq)

    def _fillObject(self,dt,eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq):
        #if dt not in self.axial:
            #self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMoment[dt][eid][sd] = [bm1,bm2]
        self.shear[dt][eid][sd] = [ts1,ts2]
        self.axial[dt][eid][sd] = af
        self.totalTorque[dt][eid][sd]   = ttrq
        self.warpingTorque[dt][eid][sd] = wtrq

    def _fillNewObject(self,dt,eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq):
        if dt not in self.axial:
            self.addNewTransient(dt)
        #self.eType[eid] = eType
        self.bendingMoment[dt][eid] = {sd:[bm1,bm2]}
        self.shear[dt][eid] = {sd:[ts1,ts2]}
        self.axial[dt][eid] = {sd:af}
        self.totalTorque[dt][eid] = {sd:ttrq}
        self.warpingTorque[dt][eid] = {sd:wtrq}

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

        self.dt = dt
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
        [eid,f41,f21,f12,f32,f23,f43,f34,f14,
             kf1,s12,kf2,s23,kf3,s34,kf4,s41] = data
        #self.eType[eid] = eType
        self.force41[eid] = f41
        self.force14[eid] = f14
        self.force21[eid] = f21
        self.force12[eid] = f12
        self.force32[eid] = f32
        self.force23[eid] = f23
        self.force43[eid] = f43
        self.force34[eid] = f34
        self.kickForce1[eid] = kf1
        self.kickForce2[eid] = kf2
        self.kickForce3[eid] = kf3
        self.kickForce4[eid] = kf4
        self.shear12[eid] = s12
        self.shear23[eid] = s23
        self.shear34[eid] = s34
        self.shear41[eid] = s41
        
    def addSort1(self,dt,data):
        [eid,f41,f21,f12,f32,f23,f43,f34,f14,
             kf1,s12,kf2,s23,kf3,s34,kf4,s41] = data
        self._fillObject(dt,eid,f41,f21,f12,f32,f23,f43,f34,f14,
                                kf1,s12,kf2,s23,kf3,s34,kf4,s41)

    def addSort2(self,eid,data):
        [dt,f41,f21,f12,f32,f23,f43,f34,f14,
            kf1,s12,kf2,s23,kf3,s34,kf4,s41] = data

        self._fillObject(dt,eid,f41,f21,f12,f32,f23,f43,f34,f14,
                                kf1,s12,kf2,s23,kf3,s34,kf4,s41)

    def _fillObject(self,dt,eid,f41,f21,f12,f32,f23,f43,f34,f14,
                                kf1,s12,kf2,s23,kf3,s34,kf4,s41):
        if dt not in self.force41:
            self.addNewTransient(dt)
        #self.eType[eid] = eType
        self.force41[dt][eid] = f41
        self.force14[dt][eid] = f14
        self.force21[dt][eid] = f21
        self.force12[dt][eid] = f12
        self.force32[dt][eid] = f32
        self.force23[dt][eid] = f23
        self.force43[dt][eid] = f43
        self.force34[dt][eid] = f34
        self.kickForce1[dt][eid] = kf1
        self.kickForce2[dt][eid] = kf2
        self.kickForce3[dt][eid] = kf3
        self.kickForce4[dt][eid] = kf4
        self.shear12[dt][eid] = s12
        self.shear23[dt][eid] = s23
        self.shear34[dt][eid] = s34
        self.shear41[dt][eid] = s41

    def __repr__(self):
        return str(self.force41)

class ComplexSpringForce(scalarObject): # 11-CELAS1,12-CELAS2,13-CELAS3, 14-CELAS4
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.force = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.dt = dt
        self.force[dt] = {}

    def add(self,dt,data):
        [eid,force] = data

        #self.eType[eid] = eType
        self.force[eid] = force

    def addSort1(self,dt,data):
        [eid,force] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = force

    def addSort2(self,eid,data):
        [dt,force] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = force

    def writeF06(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header,pageStamp,pageNum,f,isMagPhase)
        msg = header+['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 4 )\n',
                      '                                                          (REAL/IMAGINARY)\n',
                      ' \n',
                      '            FREQUENCY                    FORCE                        FREQUENCY                    FORCE\n']
        packs = []
        forces = []
        elements = []
        line = '   '
        for eid,force in sorted(self.force.items()):
            elements.append(eid)
            forces.append(force)
            #pack.append(eid)
            #pack.append(f)
            line += '%13s  %13s / %13s     ' %(eid,force.real,force.imag)
            if len(forces)==3:
                msg.append(line.rstrip()+'\n')
            ###
        ###
        if forces:
            msg.append(line.rstrip()+'\n')
        ###
        msg.append(pageStamp+str(pageNum)+'\n')

        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        words=['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 4 )\n',
               '                                                          (REAL/IMAGINARY)\n',
               ' \n',
               '                ELEMENT                                                   ELEMENT\n',
               '                  ID.                    FORCE                              ID.                    FORCE\n']
        msg = []
        for dt,Force in sorted(self.force.items()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header+words
            packs = []
            forces = []
            elements = []
            line = ''
            for eid,force in sorted(Force.items()):
                elements.append(eid)
                forces.append(force)
                #pack.append(eid)
                #pack.append(f)
                ([forceReal,forceImag],isAllZeros) = self.writeF06Floats13E([force.real,force.imag])

                line += '          %13s      %13s / %13s' %(eid,forceReal,forceImag)
                if len(forces)==2:
                    msg.append(line.rstrip()+'\n')
                    line = ''
                    forces = []
                ###
            ###
            if forces:
                msg.append(line.rstrip()+'\n')
            ###
            msg.append(pageStamp+str(pageNum)+'\n')

            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum+=1
        return (''.join(msg),pageNum-1)

    def __repr__(self):
        return str(self.force)

class ComplexDamperForce(scalarObject): # 20-CDAMP1,21-CDAMP2,22-CDAMP3,23-CDAMP4
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.force = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.dt = dt
        self.force[dt] = {}

    def add(self,dt,data):
        [eid,force] = data

        #self.eType[eid] = eType
        self.force[eid] = force

    def addSort1(self,dt,data):
        [eid,force] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = force

    def addSort2(self,eid,data):
        [dt,force] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = force

    def __repr__(self):
        return str(self.force)

class ComplexViscForce(scalarObject): # 24-CVISC
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.axialForce = {}
        self.torque = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.dt = dt
        self.axialForce[dt] = {}
        self.torque[dt] = {}

    def add(self,dt,data):
        [eid,axialForce,torque] = data

        #self.eType[eid] = eType
        self.axialForce[eid] = axialForce
        self.torque[eid] = torque

    def addSort1(self,dt,data):
        [eid,axialForce,torque] = data
        if dt not in self.axialForce:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid] = torque

    def addSort2(self,eid,data):
        [dt,axialForce,torque] = data
        if dt not in self.axialForce:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid] = torque

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

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.dt = dt
        self.mx[dt] = {}
        self.my[dt] = {}
        self.mxy[dt] = {}
        self.bmx[dt] = {}
        self.bmy[dt] = {}
        self.bmxy[dt] = {}
        self.tx[dt] = {}
        self.ty[dt] = {}

    def add(self,dt,data):
        [eid,mx,my,mxy,bmx,bmy,bmxy,tx,ty] = data

        #self.eType[eid] = eType
        self.mx[eid] = mx
        self.my[eid] = my
        self.mxy[eid] = mxy
        self.bmx[eid] = bmx
        self.bmy[eid] = bmy
        self.bmxy[eid] = bmxy
        self.tx[eid] = tx
        self.ty[eid] = ty

    def addSort1(self,dt,data):
        [eid,mx,my,mxy,bmx,bmy,bmxy,tx,ty] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid] = mx
        self.my[dt][eid] = my
        self.mxy[dt][eid] = mxy
        self.bmx[dt][eid] = bmx
        self.bmy[dt][eid] = bmy
        self.bmxy[dt][eid] = bmxy
        self.tx[dt][eid] = tx
        self.ty[dt][eid] = ty

    def addSort2(self,eid,data):
        [dt,mx,my,mxy,bmx,bmy,bmxy,tx,ty] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid] = mx
        self.my[dt][eid] = my
        self.mxy[dt][eid] = mxy
        self.bmx[dt][eid] = bmx
        self.bmy[dt][eid] = bmy
        self.bmxy[dt][eid] = bmxy
        self.tx[dt][eid] = tx
        self.ty[dt][eid] = ty

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

        self.dt = dt
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
        self.dt = dt
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
        [term,nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty] = data

        #self.eType[eid] = eType
        self.term[eid] = term
        self.ngrids[eid] = nid

        self.mx[eid] = [mx]
        self.my[eid] = [my]
        self.mxy[eid] = [mxy]
        self.bmx[eid] = [bmx]
        self.bmy[eid] = [bmy]
        self.bmxy[eid] = [bmxy]
        self.tx[eid] = [tx]
        self.ty[eid] = [ty]

    def add(self,eid,dt,data):
        [nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty] = data

        #self.eType[eid] = eType
        #print "mx = ",self.mx,mx
        self.mx[eid].append(mx)
        self.my[eid].append(my)
        self.mxy[eid].append(mxy)
        self.bmx[eid].append(bmx)
        self.bmy[eid].append(bmy)
        self.bmxy[eid].append(bmxy)
        self.tx[eid].append(tx)
        self.ty[eid].append(ty)

    def addNewElementSort1(self,eid,dt,data):
        [term,nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.term[eid] = term
        self.ngrids[eid] = nid
        self.mx[dt][eid] = [mx]
        self.my[dt][eid] = [my]
        self.mxy[dt][eid] = [mxy]
        self.bmx[dt][eid] = [bmx]
        self.bmy[dt][eid] = [bmy]
        self.bmxy[dt][eid] = [bmxy]
        self.tx[dt][eid] = [tx]
        self.ty[dt][eid] = [ty]

    def addSort1(self,eid,dt,data):
        [nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid].append(mx)
        self.my[dt][eid].append(my)
        self.mxy[dt][eid].append(mxy)
        self.bmx[dt][eid].append(bmx)
        self.bmy[dt][eid].append(bmy)
        self.bmxy[dt][eid].append(bmxy)
        self.tx[dt][eid].append(tx)
        self.ty[dt][eid].append(ty)

    def addNewElementSort2(self,dt,eid,data):
        [term,nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.term[eid] = term
        self.ngrids[eid] = nid

        self.mx[dt][eid] = [mx]
        self.my[dt][eid] = [my]
        self.mxy[dt][eid] = [mxy]
        self.bmx[dt][eid] = [bmx]
        self.bmy[dt][eid] = [bmy]
        self.bmxy[dt][eid] = [bmxy]
        self.tx[dt][eid] = [tx]
        self.ty[dt][eid] = [ty]

    def addSort2(self,dt,eid,data):
        [nid,mx,my,mxy,bmx,bmy,bmxy,tx,ty] = data
        if dt not in self.mx:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid].append(mx)
        self.my[dt][eid].append(my)
        self.mxy[dt][eid].append(mxy)
        self.bmx[dt][eid].append(bmx)
        self.bmy[dt][eid].append(bmy)
        self.bmxy[dt][eid].append(bmxy)
        self.tx[dt][eid].append(tx)
        self.ty[dt][eid].append(ty)

    def writeF06(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header,pageStamp,pageNum,f)
        raise NotImplementedError()
        #words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
        #       ' \n',
        #       '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()
        #return self._writeF06Block(words,header,pageStamp,pageNum,f)

    def writeF06Transient(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        raise NotImplementedError()

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

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.dt = dt
        self.bendingMomentA[dt] = {}
        self.bendingMomentB[dt] = {}
        self.shear[dt] = {}
        self.axial[dt] = {}
        self.torque[dt] = {}

    def add(self,dt,data):
        [eid,bm1a,bm2a,bm1b,bm2b,ts1,ts2,af,trq] = data

        #self.eType[eid] = eType
        self.bendingMomentA[eid] = [bm1a,bm2a]
        self.bendingMomentB[eid] = [bm1b,bm2b]
        self.shear[eid] = [ts1,ts2]
        self.axial[eid] = af
        self.torque[eid] = trq

    def addSort1(self,dt,data):
        [eid,bm1a,bm2a,bm1b,bm2b,ts1,ts2,af,trq] = data
        if dt not in self.axial:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMomentA[dt][eid] = [bm1a,bm2a]
        self.bendingMomentB[dt][eid] = [bm1b,bm2b]
        self.shear[dt][eid] = [ts1,ts2]
        self.axial[dt][eid] = af
        self.torque[dt][eid] = trq

    def addSort2(self,eid,data):
        [dt,bm1a,bm2a,bm1b,bm2b,ts1,ts2,af,trq] = data
        if dt not in self.axial:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMomentA[dt][eid] = [bm1a,bm2a]
        self.bendingMomentB[dt][eid] = [bm1b,bm2b]
        self.shear[dt][eid] = [ts1,ts2]
        self.axial[dt][eid] = af
        self.torque[dt][eid] = trq

    def writeF06Transient(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):

        words = ['                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B A R )\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']

        msg = []
        for dt,bendA in sorted(self.bendingMomentA.iteritems()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header+words
            for eid in sorted(bendA):
                bm1a,bm2a = self.bendingMomentA[dt][eid]
                bm1b,bm2b = self.bendingMomentB[dt][eid]
                ts1,ts2 = self.shear[dt][eid]
                af = self.axial[dt][eid]
                trq = self.torque[dt][eid]
                (vals2,isAllZeros) = self.writeF06ImagFloats13E([bm1a,bm2a,bm1b,bm2b,ts1,ts2,af,trq])
                [bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
                 bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi] = vals2
                msg.append('     %8i    %13s %13s  %13s %13s  %13s %13s  %13s  %-s\n' %(eid,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr))
                msg.append('     %8s    %13s %13s  %13s %13s  %13s %13s  %13s  %-s\n' %('', bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi))
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        ###
        return (''.join(msg),pageNum-1)

    def writeF06(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header,pageStamp,pageNum,f,isMagPhase)
        msg = header+['                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B A R )\n',
               '                                                          (REAL/IMAGINARY)\n',
               '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
               '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']

        for eid in sorted(self.bendingMomentA):
            bm1a,bm2a = self.bendingMomentA[eid]
            bm1b,bm2b = self.bendingMomentB[eid]
            ts1,ts2 = self.shear[eid]
            af = self.axial[eid]
            trq = self.torque[eid]
            (vals2,isAllZeros) = self.writeF06ImagFloats13E([bm1a,bm2a,bm1b,bm2b,ts1,ts2,af,trq])
            [bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr,
             bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi] = vals2
            msg.append('     %8i    %13s %13s  %13s %13s  %13s %13s  %13s  %-s\n' %(eid,bm1ar,bm2ar,bm1br,bm2br,ts1r,ts2r,afr,trqr))
            msg.append('     %8s    %13s %13s  %13s %13s  %13s %13s  %13s  %-s\n' %('', bm1ai,bm2ai,bm1bi,bm2bi,ts1i,ts2i,afi,trqi))
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg),pageNum)

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

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.dt = dt
        self.bendingMoment1[dt] = {}
        self.bendingMoment2[dt] = {}
        self.shearPlane1[dt] = {}
        self.shearPlane2[dt] = {}
        self.axial[dt]  = {}
        self.torque[dt] = {}

    def add(self,dt,data):
        [eid,nidA,bm1A,bm2A,sp1A,sp2A,axialA,torqueA,
             nidB,bm1B,bm2B,sp1B,sp2B,axialB,torqueB] = data

        #self.eType[eid] = eType
        self.nodeIDs[eid] = [nidA,nidB]
        self.bendingMoment1[eid] = [bm1A,bm1B]
        self.bendingMoment2[eid] = [bm2A,bm2B]
        self.shearPlane1[eid] = [sp1A,sp1B]
        self.shearPlane2[eid] = [sp2A,sp2B]
        self.axial[eid]  = [axialA,axialB]
        self.torque[eid] = [torqueA,torqueB]

    def addSort1(self,dt,data):
        [eid,nidA,bm1A,bm2A,sp1A,sp2A,axialA,torqueA,
             nidB,bm1B,bm2B,sp1B,sp2B,axialB,torqueB] = data
        self._fillObject(dt,eid,nidA,bm1A,bm2A,sp1A,sp2A,axialA,torqueA,
                                nidB,bm1B,bm2B,sp1B,sp2B,axialB,torqueB)

    def addSort2(self,eid,data):
        [dt,nidA,bm1A,bm2A,sp1A,sp2A,axialA,torqueA,
            nidB,bm1B,bm2B,sp1B,sp2B,axialB,torqueB] = data
        self._fillObject(dt,eid,nidA,bm1A,bm2A,sp1A,sp2A,axialA,torqueA,
                                nidB,bm1B,bm2B,sp1B,sp2B,axialB,torqueB)

    def _fillObject(self,dt,eid,nidA,bm1A,bm2A,sp1A,sp2A,axialA,torqueA,
                                nidB,bm1B,bm2B,sp1B,sp2B,axialB,torqueB):
        if dt not in self.axial:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.nodeIDs[eid] = [nidA,nidB]
        self.bendingMoment1[dt][eid] = [bm1A,bm1B]
        self.bendingMoment2[dt][eid] = [bm2A,bm2B]
        self.shearPlane1[dt][eid] = [sp1A,sp1B]
        self.shearPlane2[dt][eid] = [sp2A,sp2B]
        self.axial[dt][eid]  = [axialA,axialB]
        self.torque[dt][eid] = [torqueA,torqueB]

    def __repr__(self):
        return str(self.axial)

class ComplexPentaPressureForce(scalarObject): # 76-CHEXA_PR,77-PENTA_PR,78-TETRA_PR
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.acceleration = {}
        self.velocity = {}
        self.pressure = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.dt = dt
        self.acceleration[dt] = {}
        self.velocity[dt] = {}
        self.pressure[dt] = {}

    def add(self,dt,data):
        [eid,eName,ax,ay,az,vx,vy,vz,pressure] = data

        #self.eType[eid] = eType
        self.acceleration[eid] = [ax,ay,az]
        self.velocity[eid] = [vx,vy,vz]
        self.pressure[eid] = pressure

    def addSort1(self,dt,data):
        [eid,eName,ax,ay,az,vx,vy,vz,pressure] = data
        if dt not in self.acceleration:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.acceleration[dt][eid] = [ax,ay,az]
        self.velocity[dt][eid] = [vx,vy,vz]
        self.pressure[dt][eid] = pressure

    def addSort2(self,eid,data):
        [dt,eName,ax,ay,az,vx,vy,vz,pressure] = data
        if dt not in self.acceleration:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.acceleration[dt][eid] = [ax,ay,az]
        self.velocity[dt][eid] = [vx,vy,vz]
        self.pressure[dt][eid] = pressure


    def __repr__(self):
        return str(self.acceleration)

class ComplexCBUSHForce(scalarObject): # 102-CBUSH
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eType = {}
        self.force = {}
        self.moment = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.dt = dt
        self.force[dt] = {}
        self.moment[dt] = {}

    def add(self,dt,data):
        [eid,fx,fy,fz,mx,my,mz] = data

        #self.eType[eid] = eType
        self.force[eid] = [fx,fy,fz]
        self.moment[eid] = [mx,my,mz]

    def addSort1(self,dt,data):
        [eid,fx,fy,fz,mx,my,mz] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = [fx,fy,fz]
        self.moment[dt][eid] = [mx,my,mz]

    def addSort2(self,eid,data):
        [dt,fx,fy,fz,mx,my,mz] = data
        if dt not in self.force:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = [fx,fy,fz]
        self.moment[dt][eid] = [mx,my,mz]

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
        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.dt = dt
        self.forceX[dt]  = {}
        self.shearY[dt]  = {}
        self.shearZ[dt]  = {}
        self.torsion[dt]  = {}
        self.bendingY[dt]  = {}
        self.bendingZ[dt]  = {}

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
            [nid,posit,forceX,shearY,shearZ,torsion,bendingY,bendingZ] = force
            self.forceX[eid][nid]  = forceX
            self.shearY[eid][nid]  = shearY
            self.shearZ[eid][nid]  = shearZ
            self.torsion[eid][nid]  = torsion
            self.bendingY[eid][nid] = bendingY
            self.bendingZ[eid][nid] = bendingZ

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
            [nid,posit,forceX,shearY,shearZ,torsion,bendingY,bendingZ] = force
            self.forceX[dt][eid][nid]   = forceX
            self.shearY[dt][eid][nid]   = shearY
            self.shearZ[dt][eid][nid]   = shearZ
            self.torsion[dt][eid][nid]  = torsion
            self.bendingY[dt][eid][nid] = bendingY
            self.bendingZ[dt][eid][nid] = bendingZ

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
            [nid,posit,forceX,shearY,shearZ,torsion,bendingY,bendingZ] = force
            self.forceX[dt][eid][nid]   = forceX
            self.shearY[dt][eid][nid]   = shearY
            self.shearZ[dt][eid][nid]   = shearZ
            self.torsion[dt][eid][nid]  = torsion
            self.bendingY[dt][eid][nid] = bendingY
            self.bendingZ[dt][eid][nid] = bendingZ

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
        self.dt = dt
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
            [nid,membraneX,membraneY,membraneXY,bendingX,bendingY,bendingXY,shearYZ,shearXZ] = force
            self.membraneX[eid][nid]  = membraneX
            self.membraneY[eid][nid]  = membraneY
            self.membraneXY[eid][nid] = membraneXY
            self.bendingX[eid][nid]   = bendingX
            self.bendingY[eid][nid]   = bendingY
            self.bendingXY[eid][nid]  = bendingXY
            self.shearYZ[eid][nid]  = shearYZ
            self.shearXZ[eid][nid]  = shearXZ

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
            [nid,membraneX,membraneY,membraneXY,bendingX,bendingY,bendingXY,shearYZ,shearXZ] = force
            self.membraneX[dt][eid][nid]  = membraneX
            self.membraneY[dt][eid][nid]  = membraneY
            self.membraneXY[dt][eid][nid] = membraneXY
            self.bendingX[dt][eid][nid]   = bendingX
            self.bendingY[dt][eid][nid]   = bendingY
            self.bendingXY[dt][eid][nid]  = bendingXY
            self.shearYZ[dt][eid][nid]  = shearYZ
            self.shearXZ[dt][eid][nid]  = shearXZ

    def __repr__(self):
        return str(self.membraneX)
