from pyNastran.op2.resultObjects.op2_Objects import scalarObject

class stressObject(scalarObject):
    def __init__(self,dataCode,iSubcase):
        scalarObject.__init__(self,dataCode,iSubcase)

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        print "updating dt...dt=%s" %(dt)
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###

class strainObject(scalarObject):
    def __init__(self,dataCode,iSubcase):
        scalarObject.__init__(self,dataCode,iSubcase)
    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        print "updating dt...dt=%s" %(dt)
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###


