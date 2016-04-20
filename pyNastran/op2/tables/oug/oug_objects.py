from six import  iteritems
from struct import pack
from pyNastran.op2.result_objects.op2_objects import ScalarObject
from numpy import array

#class staticFlux(scalarObject): # approach_code=1, table_code=3 - whatever the static version of this is...


class Flux(ScalarObject):  # approach_code=1, table_code=3, thermal=1
    def __init__(self, data_code, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)

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

    def add(self, node_id, grid_type, v1, v2, v3, v4=None, v5=None, v6=None):
        assert 0 < node_id < 1000000000, 'node_id=%s' % (node_id)
        assert node_id not in self.fluxes
        self.fluxes[node_id] = array([v1, v2, v3])

    def write_op2(self, block3, device_code=1):
        """
        Creates the binary data for writing the table
        .. warning:: hasnt been tested...
        """
        msg = block3
        for node_id, flux in sorted(iteritems(self.fluxes)):
            grid = node_id * 10 + device_code
            msg += pack('i6f', grid, flux[0], flux[1], flux[2], 0, 0, 0)
        return msg
