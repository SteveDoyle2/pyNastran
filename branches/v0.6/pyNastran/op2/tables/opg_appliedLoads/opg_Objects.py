from numpy import array
from pyNastran.op2.resultObjects.op2_Objects import scalarObject


class AppliedLoadsObject(scalarObject):  # approach_code=1, sort_code=0

    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        scalarObject.__init__(self, data_code, isubcase)
        self.dt = dt

        self.eids = {}
        self.sources = {}
        self.forces = {}
        self.moments = {}
        if self.dt is not None:
            assert dt >= 0.
            raise NotImplementedError(
                'transient appliedLoads not implemented...')
            self.eids = {dt: []}
            self.sources = {dt: []}
            self.forces = {dt: []}
            self.moments = {dt: []}
            self.add = self.addTransient

    def add_node(self, nodeID, eid, source, v1, v2, v3, v4, v5, v6):
        msg = "nodeID=%s eid=%s source=|%s| v1=%i v2=%i v3=%i v4=%i v5=%i v6=%i" % (nodeID, eid, source, v1, v2, v3, v4, v5, v6)
        assert 0 < nodeID < 1000000000, msg
        assert nodeID not in self.forces, msg

        self.eids[nodeID] = [eid]
        self.sources[nodeID] = [source]
        self.forces[nodeID] = [array([v1, v2, v3])]  # Fx,Fy,Fz
        self.moments[nodeID] = [array([v4, v5, v6])]  # Mx,My,Mz

    def add(self, nodeID, eid, source, v1, v2, v3, v4, v5, v6):
        msg = "nodeID=%s eid=%s source=|%s| v1=%i v2=%i v3=%i v4=%i v5=%i v6=%i" % (nodeID, eid, source, v1, v2, v3, v4, v5, v6)
        assert 0 < nodeID < 1000000000, msg
        if nodeID not in self.forces:
            self.add_node(nodeID, eid, source, v1, v2, v3, v4, v5, v6)
            return None
        #assert nodeID not in self.forces,msg

        self.eids[nodeID].append(eid)
        self.sources[nodeID].append(source)
        self.forces[nodeID].append(array([v1, v2, v3]))  # Fx,Fy,Fz
        self.moments[nodeID].append(array([v4, v5, v6]))  # Mx,My,Mz

    def addTransient(self, nodeID, eid, source, v1, v2, v3, v4, v5, v6):
        raise Exception('no implemented')
        msg = "nodeID=%s v1=%s v2=%s v3=%s" % (nodeID, v1, v2, v3)
        assert 0 < nodeID < 1000000000, msg
        assert nodeID not in self.forces

        self.forces[self.dt][nodeID] = array([v1, v2, v3])  # Fx,Fy,Fz
        self.moments[self.dt][nodeID] = array([v4, v5, v6])  # Mx,My,Mz
