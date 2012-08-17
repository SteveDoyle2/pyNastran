from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.op2.resultObjects.op2_Objects import scalarObject


class HeatFlux_VU_3D(scalarObject):  # 146-VUPENTA, 147-VUTETRA, 148-VUPENTA
    def __init__(self, dataCode, isSort1, iSubcase, dt):
        scalarObject.__init__(self, dataCode, iSubcase)
        #self.eType = {}
        self.parent = {}

        self.grad = {}
        self.flux = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self, dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def add(self, nNodes, dt, data):
        [eid, parent, gradFluxes] = data
        self.parent[eid] = parent
        #self.eType[eid]    = eType

        self.grad[eid] = {}
        self.flux[eid] = {}
        for gradFlux in gradFluxes:
            [nid, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = gradFlux
            self.grad[eid][nid] = [xGrad, yGrad, zGrad]
            self.flux[eid][nid] = [xFlux, yFlux, zFlux]

    def addSort1(self, nNodes, dt, data):
        [eid, parent, gradFluxes] = data
        if dt not in self.grad:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for gradFlux in gradFluxes:
            [nid, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = gradFlux
            self.grad[dt][eid][nid] = [xGrad, yGrad, zGrad]
            self.flux[dt][eid][nid] = [xFlux, yFlux, zFlux]

    def addSort2(self, nNodes, eid, data):
        [dt, parent, gradFluxes] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for gradFlux in gradFluxes:
            [nid, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = gradFlux
            self.grad[dt][eid][nid] = [xGrad, yGrad, zGrad]
            self.flux[dt][eid][nid] = [xFlux, yFlux, zFlux]

    def __repr__(self):
        return str(self.grad)


class HeatFlux_VU(scalarObject):  # 189-VUQUAD 190-VUTRIA,191-VUBEAM
    def __init__(self, dataCode, isSort1, iSubcase, dt):
        scalarObject.__init__(self, dataCode, iSubcase)
        #self.eType = {}
        self.parent = {}
        self.coord = {}
        self.icord = {}
        self.theta = {}

        self.grad = {}
        self.flux = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def addNewTransient(self, dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def add(self, nNodes, dt, data):
        [eid, parent, coord, icord, theta, gradFluxes] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.grad[eid] = {}
        self.flux[eid] = {}
        for gradFlux in gradFluxes:
            [nid, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = gradFlux
            self.grad[eid][nid] = [xGrad, yGrad, zGrad]
            self.flux[eid][nid] = [xFlux, yFlux, zFlux]

    def addSort1(self, nNodes, dt, data):
        [eid, parent, coord, icord, theta, gradFluxes] = data
        if dt not in self.grad:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for gradFlux in gradFluxes:
            [nid, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = gradFlux
            self.grad[dt][eid][nid] = [xGrad, yGrad, zGrad]
            self.flux[dt][eid][nid] = [xFlux, yFlux, zFlux]

    def addSort2(self, nNodes, eid, data):
        [dt, parent, coord, icord, theta, gradFluxes] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for gradFlux in gradFluxes:
            [nid, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = gradFlux
            self.grad[dt][eid][nid] = [xGrad, yGrad, zGrad]
            self.flux[dt][eid][nid] = [xFlux, yFlux, zFlux]

    def __repr__(self):
        return str(self.grad)


class HeatFlux_VUBEAM(scalarObject):  # 191-VUBEAM
    def __init__(self, dataCode, isSort1, iSubcase, dt):
        scalarObject.__init__(self, dataCode, iSubcase)
        #self.eType = {}
        self.parent = {}
        self.coord = {}
        self.icord = {}

        self.grad = {}
        self.flux = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self, dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def add(self, nNodes, dt, data):
        [eid, parent, coord, icord, gradFluxes] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.grad[eid] = {}
        self.flux[eid] = {}
        for gradFlux in gradFluxes:
            [nid, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = gradFlux
            self.grad[eid][nid] = [xGrad, yGrad, zGrad]
            self.flux[eid][nid] = [xFlux, yFlux, zFlux]

    def addSort1(self, nNodes, dt, data):
        [eid, parent, coord, icord, gradFluxes] = data
        if dt not in self.grad:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for gradFlux in gradFluxes:
            [nid, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = gradFlux
            self.grad[dt][eid][nid] = [xGrad, yGrad, zGrad]
            self.flux[dt][eid][nid] = [xFlux, yFlux, zFlux]

    def addSort2(self, nNodes, eid, data):
        [dt, parent, coord, icord, gradFluxes] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for gradFlux in gradFluxes:
            [nid, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = gradFlux
            self.grad[dt][eid][nid] = [xGrad, yGrad, zGrad]
            self.flux[dt][eid][nid] = [xFlux, yFlux, zFlux]

    def __repr__(self):
        return str(self.grad)


class HeatFlux_1D(scalarObject):  # 1-ROD, 2-BEAM, 3-TUBE, 10-CONROD, 34-BAR, 69-BEND
    def __init__(self, dataCode, isSort1, iSubcase, dt):
        scalarObject.__init__(self, dataCode, iSubcase)
        self.eType = {}
        self.grad = {}
        self.flux = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

    def addNewTransient(self, dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def add(self, dt, data):
        [eid, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = data
        self.eType[eid] = eType
        self.grad[eid] = [xGrad, yGrad, zGrad]
        self.flux[eid] = [xFlux, yFlux, zFlux]

    def addSort1(self, dt, data):
        [eid, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = data
        if dt not in self.grad:
            self.addNewTransient(dt)
        self.eType[eid] = eType
        self.grad[dt][eid] = [xGrad, yGrad, zGrad]
        self.flux[dt][eid] = [xFlux, yFlux, zFlux]

    def addSort2(self, eid, data):
        [dt, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.eType[eid] = eType
        self.grad[dt][eid] = [xGrad, yGrad, zGrad]
        self.flux[dt][eid] = [xFlux, yFlux, zFlux]

    def __repr__(self):
        return str(self.grad)


class HeatFlux_2D_3D(scalarObject):  # 33-QUAD4, 39-TETRA, 53-TRIAX6,64-QUAD8, 67-HEXA, 68-PENTA, 74-TRIA3, 75-TRIA6
    def __init__(self, dataCode, isSort1, iSubcase, dt):
        scalarObject.__init__(self, dataCode, iSubcase)
        self.eType = {}
        self.grad = {}
        self.flux = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            self.add = self.addSort1
        else:
            self.add = self.addSort2
        ###

    def addNewTransient(self, dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def addSort1(self, dt, data):
        [eid, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = data
        if dt not in self.grad:
            self.addNewTransient(dt)
        self.eType[eid] = eType
        self.grad[dt][eid] = [xGrad, yGrad, zGrad]
        self.flux[dt][eid] = [xFlux, yFlux, zFlux]

    def addSort2(self, eid, data):
        [dt, eType, xGrad, yGrad, zGrad, xFlux, yFlux, zFlux] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.eType[eid] = eType
        self.grad[dt][eid] = [xGrad, yGrad, zGrad]
        self.flux[dt][eid] = [xFlux, yFlux, zFlux]

    def __repr__(self):
        return str(self.grad)


class HeatFlux_CONV(scalarObject):  # 110-CONV
    def __init__(self, dataCode, isSort1, iSubcase, dt):
        scalarObject.__init__(self, dataCode, iSubcase)
        self.cntlNode = {}
        self.freeConv = {}
        self.freeConvK = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def addNewTransient(self, dt):
        self.cntlNode[dt] = {}
        self.freeConv[dt] = {}
        self.freeConvK[dt] = {}

    def add(self, dt, data):
        [eid, cntlNode, freeConv, freeConvK] = data
        #self.eType[eid]     = eType
        self.cntlNode[eid] = cntlNode
        self.freeConv[eid] = freeConv
        self.freeConvK[eid] = freeConvK

    def addSort1(self, dt, data):
        [eid, cntlNode, freeConv, freeConvK] = data
        if dt not in self.freeConv:
            self.addNewTransient(dt)
        #self.eType[eid]     = eType
        self.cntlNode[dt][eid] = cntlNode
        self.freeConv[dt][eid] = freeConv
        self.freeConvK[dt][eid] = freeConvK

    def addSort2(self, eid, data):
        [dt, eType, fApplied, freeConv, forceConv, fRad, fTotal] = data
        if dt not in self.freeConv:
            self.addNewTransient(dt)
        #self.eType[eid]     = eType
        self.fApplied[dt][eid] = fApplied
        self.freeConv[dt][eid] = freeConv
        self.forceConv[dt][eid] = forceConv
        self.fRad[dt][eid] = fRad
        self.fTotal[dt][eid] = fTotal

    def __repr__(self):
        return str(self.fApplied)


class HeatFlux_CHBDYx(scalarObject):  # 107-CHBDYE 108-CHBDYG 109-CHBDYP
    def __init__(self, dataCode, isSort1, iSubcase, dt):
        scalarObject.__init__(self, dataCode, iSubcase)
        self.eType = {}
        self.fApplied = {}
        self.freeConv = {}
        self.forceConv = {}
        self.fRad = {}
        self.fTotal = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def addNewTransient(self, dt):
        self.fApplied[dt] = {}
        self.freeConv[dt] = {}
        self.forceConv[dt] = {}
        self.fRad[dt] = {}
        self.fTotal[dt] = {}

    def add(self, dt, data):
        [eid, eType, fApplied, freeConv, forceConv, fRad, fTotal] = data

        self.eType[eid] = eType
        self.fApplied[eid] = fApplied
        self.freeConv[eid] = freeConv
        self.forceConv[eid] = forceConv
        self.fRad[eid] = fRad
        self.fTotal[eid] = fTotal

    def addSort1(self, dt, data):
        [eid, eType, fApplied, freeConv, forceConv, fRad, fTotal] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.eType[eid] = eType
        self.fApplied[dt][eid] = fApplied
        self.freeConv[dt][eid] = freeConv
        self.forceConv[dt][eid] = forceConv
        self.fRad[dt][eid] = fRad
        self.fTotal[dt][eid] = fTotal

    def addSort2(self, eid, data):
        [dt, eType, fApplied, freeConv, forceConv, fRad, fTotal] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.eType[eid] = eType
        self.fApplied[dt][eid] = fApplied
        self.freeConv[dt][eid] = freeConv
        self.forceConv[dt][eid] = forceConv
        self.fRad[dt][eid] = fRad
        self.fTotal[dt][eid] = fTotal

    def __repr__(self):
        return str(self.fApplied)
