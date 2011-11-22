from struct import pack,unpack

from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array

class nonlinearFluxObject(scalarObject): # approachCode=10, sortCode=0
    def __init__(self,iSubcase,loadStep):
        scalarObject.__init__(self,iSubcase)

        self.loadStep = loadStep
        self.eTypes = {}
        if loadStep is None:
            self.fluxes = {}
            self.gradients = {}
        else:
            self.fluxes    = {loadStep: {}}
            self.gradients = {loadStep: {}}
            #raise Exception('transient not supported for flux yet...')
        ###

    def add(self,nodeID,eType,v1,v2,v3,v4=None,v5=None,v6=None):
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        #print "nodeID=%s eType=%s v1=%s v2=%s v3=%s v4=%s v5=%s v6=%s" %(nodeID,eType,v1,v2,v3,v4,v5,v6)
        assert nodeID not in self.fluxes[self.loadStep],'nodeID=%s' %(nodeID)
        self.gradients[self.loadStep][nodeID] = array([v1,v2,v3])
        self.fluxes[   self.loadStep][nodeID] = array([v4,v5,v6])
        self.eTypes[nodeID] = eType

    def updateDt(self,loadStep):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        assert loadStep>=0.
        self.loadStep = loadStep
        self.temperatures[loadStep] = {}

    def __repr__(self):
        msg = '---NONLINEAR GRADIENTS & HEAT FLUX---\n'
        msg += 'loadStep = %g\n' %(self.loadStep)

        for dt,fluxPack in sorted(self.fluxes.items()):
            msg += '%-10s %-8s %-10s %-10s %-10s %-10s %-10s %-10s\n' %('GRID','eType','xGrad','yGrad','zGrad','xFlux','yFlux','zFlux')
            
            for nodeID,flux in sorted(fluxPack.items()):
                eType = self.eTypes[nodeID]
                msg += '%-10i %-8s ' %(nodeID,eType)
                gradients = self.gradients[dt][nodeID]

                for val in list(gradients)+list(flux):
                    if abs(val)<1e-6:
                        msg += '%-10s ' %('0.')
                    else:
                        msg += '%-10i ' %(val)
                    ###
                msg += '\n'
            ###
        return msg

