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
