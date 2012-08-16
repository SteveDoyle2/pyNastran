import sys
from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import scalarObject


#class staticFluxObj(scalarObject): # approachCode=1, tableCode=3 - whatever the static version of this is...


class fluxObject(scalarObject):  # approachCode=1, tableCode=3, thermal=1
    def __init__(self, dataCode, iSubcase, dt=None):
        scalarObject.__init__(self, dataCode, iSubcase)

        self.dt = dt
        self.fluxes = {}
        if dt is not None:
            self.fluxes = {}
            self.isTransient = True
            raise NotImplementedError('transient fluxObject is supported...')

    def deleteTransient(self, dt):
        del self.fluxes[dt]

    def getTransients(self):
        k = self.fluxes.keys()
        k.sort()
        return k

    def add(self, nodeID, gridType, v1, v2, v3, v4=None, v5=None, v6=None):
        assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)
        assert nodeID not in self.fluxes
        self.fluxes[nodeID] = array([v1, v2, v3])

    def writeOp2(self, block3, deviceCode=1):
        """
        creates the binary data for writing the table
        @warning hasnt been tested...
        """
        msg = block3
        for nodeID, flux in sorted(self.fluxes.iteritems()):
            grid = nodeID * 10 + deviceCode
            msg += pack('iffffff', grid, flux[0], flux[1], flux[2], 0, 0, 0)
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---HEAT FLUX---\n'
        msg += '%-10s %-8s %-8s %-8s\n' % ('NodeID', 'xFlux', 'yFlux', 'zFlux')
        for nodeID, flux in sorted(self.fluxes.iteritems()):
            msg += '%10i ' % (nodeID)

            for val in flux:
                if abs(val) < 1e-6:
                    msg += '%10s' % (0)
                else:
                    msg += '%10.3e ' % (val)
            msg += '\n'
        return msg
