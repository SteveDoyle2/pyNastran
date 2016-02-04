# pylint: disable=E1101,C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from numpy import array

from pyNastran.op2.result_objects.op2_objects import ScalarObject


class NonlinearFlux(ScalarObject):  # approach_code=10, sort_code=0
    def __init__(self, data_code, isubcase, load_step):
        ScalarObject.__init__(self, data_code, isubcase)

        self.load_step = load_step
        self.eTypes = {}
        self.fluxes = {}
        self.gradients = {}
        if load_step is not None:
            self.add_new_transient()
            #self.isTransient = True
            #raise Exception('transient not supported for flux yet...')

    def update_dt(self, data_code, load_step):
        self.data_code = data_code
        self.apply_data_code()
        assert load_step >= 0.
        self.load_step = load_step
        self.add_new_transient()

    def add_new_transient(self):
        """
        initializes the transient variables
        .. note:: make sure you set self.dt first
        """
        self.fluxes[self.load_step] = {}
        self.gradients[self.load_step] = {}

    def add(self, nodeID, eType, v1, v2, v3, v4=None, v5=None, v6=None):
        assert 0 < nodeID < 1000000000, 'nodeID=%s' % (nodeID)
        #print("nodeID=%s eType=%s v1=%s v2=%s v3=%s v4=%s v5=%s v6=%s"
        #    %(nodeID,eType,v1,v2,v3,v4,v5,v6))
        assert nodeID not in self.fluxes[self.load_step], 'nodeID=%s' % (nodeID)
        self.gradients[self.load_step][nodeID] = array([v1, v2, v3])
        self.fluxes[self.load_step][nodeID] = array([v4, v5, v6])
        self.eTypes[nodeID] = eType

