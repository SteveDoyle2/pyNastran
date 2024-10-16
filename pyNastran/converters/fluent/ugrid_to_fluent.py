import numpy as np
from pyNastran.utils import PathLike
from pyNastran.converters.aflr.ugrid.ugrid_reader import read_ugrid
from pyNastran.converters.fluent.fluent import Fluent, read_fluent

def ugrid_to_fluent_filename(ugrid_filename: PathLike,
                             fluent_filename: PathLike) -> Fluent:
    model = ugrid_to_fluent(ugrid_filename)
    model.write_fluent(fluent_filename)
    return model

def ugrid_to_fluent(ugrid_filename: PathLike) -> Fluent:
    ugrid_model = read_ugrid(ugrid_filename)
    fluent_model = Fluent(auto_read_write_h5=False)

    # self.tets = np.array([], dtype='int32')
    # self.penta5s = np.array([], dtype='int32')
    # self.penta6s = np.array([], dtype='int32')
    # self.hexas = np.array([], dtype='int32')
    nnodes = len(ugrid_model.nodes)
    ntri = len(ugrid_model.tris)
    nquad = len(ugrid_model.quads)
    neids = ntri + nquad
    assert neids > 0, (ntri, nquad)
    element_id = np.arange(1, neids+1)
    property_id = ugrid_model.pids

    eid_tri = element_id[:ntri]
    eid_quad = element_id[ntri:]
    pid_tri = property_id[:ntri]
    pid_quad = property_id[ntri:]
    assert len(eid_quad) == nquad
    if ntri:
        assert len(eid_tri) == ntri
        tris = np.column_stack([eid_tri, pid_tri, ugrid_model.tris])
        print(ugrid_model.tris.shape, tris.shape)
        assert tris.shape[1] == 6, tris.shape
        fluent_model.tris = tris

    if nquad:
        quads = np.column_stack([eid_quad, pid_quad, ugrid_model.quads])
        assert quads.shape[1] == 6, quads.shape
        fluent_model.quads = quads

    fluent_model.node_id = np.arange(1, nnodes+1)
    fluent_model.xyz = ugrid_model.nodes
    fluent_model.element_id = element_id
    fluent_model.titles = ['ShellID', 'PropertyID']
    fluent_model.results = property_id.reshape((neids, 1))
    return fluent_model
