from __future__ import print_function
from six import iteritems, itervalues
from six.moves import zip
import scipy
from pyNastran.bdf.bdf import BDF
from numpy import array, zeros, unique, where, arange, hstack, vstack, searchsorted

def bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol):
    """
    Equivalences nodes; keeps the lower node id; creates two nodes with the same

    .. warning:: only handles CQUAD4, CTRIA3
    .. warning:: assumes cid=0
    """
    model = BDF()
    model.read_bdf(bdf_filename, xref=True)

    # quads / tris
    quads = []
    quadmap = []
    tris = []
    trimap = []

    inode = 1
    nid_map = {}
    for nid, node in sorted(iteritems(model.nodes)):
        node.nid = inode
        nid_map[inode - 1] = nid
        inode += 1
    #model.write_bdf('A_' + bdf_filename_out)

    nids = array([node.nid for nid, node in sorted(iteritems(model.nodes))], dtype='int32')
    #nids_used = set([])
    #for element in itervalues(model.elements):
    #    nids_used.update(element.node_ids)
    #nids_used = array(list(nids_used), dtype='int32')
    #nids_used.sort()
    #nids = nids_used

    nnodes = len(nids)
    i = arange(nnodes, dtype='int32')
    nids2 = vstack([i, nids]).T

    nodes_xyz = array([node.xyz for nid, node in sorted(iteritems(model.nodes))], dtype='float32')

    nids_new = set([])
    for eid, element in sorted(iteritems(model.elements)):
        emap = []

        if element.type == 'CQUAD4':
            quads.append(element.node_ids)
            quadmap.append(element.eid)
        elif element.type == 'CTRIA3':
            tris.append(element.node_ids)
            trimap.append(element.eid)
        else:
            raise NotImplementedError(element.type)

    quads = array(quads, dtype='int32') - 1
    quadmap = array(quadmap, dtype='int32')
    tris = array(tris, dtype='int32') - 1
    trimap = array(trimap, dtype='int32')

    # build the kdtree
    try:
        kdt = scipy.spatial.cKDTree(list(nodes_xyz))
    except RuntimeError:
        print(nodes_xyz)
        raise RuntimeError(nodes_xyz)

    # find the node ids of interest
    nids_new = hstack([unique(quads), unique(tris)])
    nids_new.sort()
    inew = searchsorted(nids, nids_new, side='left')

    # check the closest 10 nodes for equality
    deq, ieq = kdt.query(nodes_xyz[inew, :], k=10, distance_upper_bound=tol)

    # get the ids of the duplicate nodes
    slots = where(ieq[:, 1:] < nnodes)
    irows, icols = slots
    replacer = unique(ieq[slots])

    skip_nodes = []
    for (irow, icol) in zip(irows, icols):
        nid1 = nid_map[irow]
        nid2 = nid_map[icol]

        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]

        node2.nid = node1.nid
        node2.xyz = node1.xyz
        assert node2.cp == node1.cp
        assert node2.cd == node1.cd
        assert node2.ps == node1.ps
        assert node2.seid == node1.seid
        skip_nodes.append(nid2)

    model.remove_nodes = skip_nodes
    #model._write_nodes = _write_nodes
    model.write_bdf(bdf_filename_out)

def cut_model(model, axis='-y'):
    #model.read_bdf(bdf_filename, xref=False)
    #model.cross_reference(xref_loads=xref_loads)

    if axis == '-x':
        iaxis = 0
    elif axis == '-y':
        iaxis = 1
    elif axis == '-z':
        iaxis = 2
    else:
        raise NotImplementedError(axis)

    remove_nids = []
    for nid, node in iteritems(model.nodes):
        c = node.Position()
        if c[iaxis] < 0.0:
            remove_nids.append(nid)

    remove_eids = []
    for eid, element in iteritems(model.elements):
        c = element.Centroid()
        if c[iaxis] < 0.0:
            remove_eids.append(eid)

    for nid in remove_nids:
        del model.nodes[nid]
    for eid in remove_eids:
        del model.elements[eid]

    loads2 = {}
    for load_id, loadcase in iteritems(model.loads):
        loadcase2 = []
        for load in loadcase:
            if load.type == 'LOAD':
                loadcase2.append(load)
            elif load.type == 'PLOAD4':
                if load.eid not in remove_eids:
                    loadcase2.append(load)
            else:
                loadcase2.append(load)
        loads2[load_id] = loadcase2
    model.loads = loads2

def _write_nodes(self, outfile, size, is_double):
    """
    Writes the NODE-type cards

    :param self: the BDF object
    """
    if self.spoints:
        msg = []
        msg.append('$SPOINTS\n')
        msg.append(self.spoints.write_card(size, is_double))
        outfile.write(''.join(msg))

    if self.nodes:
        msg = []
        msg.append('$NODES\n')
        if self.gridSet:
            msg.append(self.gridSet.print_card(size))
        for (nid, node) in sorted(iteritems(self.nodes)):
            if nid not in self.remove_nodes:
                msg.append(node.write_card(size, is_double))
        outfile.write(''.join(msg))


def bdf_renumber(bdf_filename, bdf_filename_out):
    """
    ..warning:: only supports nodes/properties/elements
    """
    model = BDF()
    model.read_bdf(bdf_filename)

    nid = 1
    for nidi, node in iteritems(model.nodes):
        node.nid = nid
        nid += 1

    pid = 1
    for pidi, prop in iteritems(model.properties):
        prop.pid = pid
        pid += 1

    eid = 1
    for eidi, element in iteritems(model.elements):
        element.eid = eid
        eid += 1
    model.write_bdf(bdf_filename_out)


def main():
    msg = 'CEND\n'
    msg += 'BEGIN BULK\n'
    msg += 'GRID,1,,0.,0.,0.\n'
    msg += 'GRID,2,,0.,0.,0.5\n'
    msg += 'GRID,3,,0.,0.,0.51\n'
    msg += 'GRID,10,,0.,0.,1.\n'
    msg += 'GRID,11,,0.,0.,1.\n'
    msg += 'CTRIA3,1,1,1,2,11\n'
    msg += 'CTRIA3,2,1,1,2,11\n'
    msg += 'CTRIA3,3,1,2,3,11\n'
    msg += 'CTRIA3,4,1,1,2,10\n'
    msg += 'PSHELL,1,1,0.1\n'
    msg += 'MAT1,1,3.0,, 0.3\n'
    msg += 'ENDDATA'


    bdf_filename = 'nonunique.bdf'

    f = open(bdf_filename, 'wb')
    f.write(msg)
    f.close()

    bdf_filename_out ='unique.bdf'
    tol = 0.2
    bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol)

if __name__ == '__main__':
    main()
