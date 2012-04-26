class ROD(object): # 1-ROD, 3-TUBE, 10-CONROD
    def __init__(self,isSort1,dt):
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
        [eid,axialForce,torque] = data

        #self.eType[eid] = eType        
        self.axialForce[eid] = {}
        self.torque[eid] = {}

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

class CBEAM(object): # 2-CBEAM
    def __init__(self,isSort1,dt):
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
        [eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        print "CBEAM addnew",data
        #self.eType[eid] = eType
        self.bendingMoment[eid] = {nid:[bm1,bm2]}
        self.shear[eid] = {nid:[ts1,ts2]}
        self.axial[eid] = {nid:af}
        self.totalTorque[eid] = {nid:ttrq}
        self.warpingTorque[eid] = {nid:wtrq}

    def add(self,dt,data):
        [eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        print "CBEAM add   ",data

        #self.eType[eid] = eType
        self.bendingMoment[eid][nid] = [bm1,bm2]
        self.shear[eid][nid] = [ts1,ts2]
        self.axial[eid][nid] = af
        self.totalTorque[eid][nid] = ttrq
        self.warpingTorque[eid][nid] = wtrq

    def addNewElementSort1(self,dt,data):
        [eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        if dt not in self.axial:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMoment[dt][eid] = {nid:[bm1,bm2]}
        self.shear[dt][eid] = {nid:[ts1,ts2]}
        self.axial[dt][eid] = {nid:af}
        self.totalTorque[dt][eid] = {nid:ttrq}
        self.warpingTorque[dt][eid] = {nid:wtrq}

    def addSort1(self,dt,data):
        [eid,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        #if dt not in self.axial:
            #self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMoment[dt][eid][nid] = [bm1,bm2]
        self.shear[dt][eid][nid] = [ts1,ts2]
        self.axial[dt][eid][nid] = af
        self.totalTorque[dt][eid][nid] = ttrq
        self.warpingTorque[dt][eid][nid] = wtrq

    def addNewElementSort2(self,eid,data):
        [dt,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        if dt not in self.axial:
            self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMoment[dt][eid] = {nid:[bm1,bm2]}
        self.shear[dt][eid] = {nid:[ts1,ts2]}
        self.axial[dt][eid] = {nid:af}
        self.totalTorque[dt][eid] = {nid:ttrq}
        self.warpingTorque[dt][eid] = {nid:wtrq}

    def addSort2(self,eid,data):
        [dt,nid,sd,bm1,bm2,ts1,ts2,af,ttrq,wtrq] = data
        #if dt not in self.axial:
            #self.addNewTransient(dt)

        #self.eType[eid] = eType
        self.bendingMoment[dt][eid][nid] = [bm1,bm2]
        self.shear[dt][eid][nid] = [ts1,ts2]
        self.axial[dt][eid][nid] = af
        self.totalTorque[dt][eid][nid] = ttrq
        self.warpingTorque[dt][eid][nid] = wtrq

    def __repr__(self):
        return str(self.axial)

class SPRING(object): # 11-CELAS1,12-CELAS2,13-CELAS3, 14-CELAS4
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
        [eid,force] = data

        #self.eType[eid] = eType
        self.force[eid] = {}

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

class CBAR(object): # 34-CBAR
    def __init__(self,isSort1,dt):
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

    def __repr__(self):
        return str(self.axial)

class CBUSH(object): # 102-CBUSH
    def __init__(self,isSort1,dt):
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

