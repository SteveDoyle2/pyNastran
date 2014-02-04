from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import scalarObject


#class staticFluxObj(scalarObject): # approach_code=1, table_code=3 - whatever the static version of this is...


class FluxObject(scalarObject):  # approach_code=1, table_code=3, thermal=1
    def __init__(self, data_code, isubcase, dt=None):
        scalarObject.__init__(self, data_code, isubcase)

        self.dt = dt
        self.fluxes = {}
        if dt is not None:
            self.fluxes = {}
            self.isTransient = True
            raise NotImplementedError('transient fluxObject is supported...')

    def delete_transient(self, dt):
        del self.fluxes[dt]

    def get_transients(self):
        k = self.fluxes.keys()
        k.sort()
        return k

    def add(self, nodeID, gridType, v1, v2, v3, v4=None, v5=None, v6=None):
        assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)
        assert nodeID not in self.fluxes
        self.fluxes[nodeID] = array([v1, v2, v3])

    def write_op2(self, block3, device_code=1):
        """
        Creates the binary data for writing the table
        .. warning:: hasnt been tested...
        """
        msg = block3
        for nodeID, flux in sorted(self.fluxes.iteritems()):
            grid = nodeID * 10 + device_code
            msg += pack('i6f', grid, flux[0], flux[1], flux[2], 0, 0, 0)
        return msg
