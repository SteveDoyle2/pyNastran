from __future__ import print_function
from six import iteritems, itervalues
from six.moves import zip, range
import scipy
from pyNastran.bdf.bdf import BDF
from numpy import array, unique, where, arange, hstack, vstack, searchsorted, unique


def remove_unassociated_nodes(bdf_filename, bdf_filename_out, renumber=False):
    """
    Removes nodes from a model that are not referenced.

    .. warning only considers elements
    .. renumber=False is not supported
    """
    model = BDF()
    model.read_bdf(bdf_filename, xref=True)

    nids_used = set([])
    for element in itervalues(model.elements):
        nids_used.update(element.node_ids)
    all_nids = set(model.nodes.keys())

    nodes_to_remove = all_nids - nids_used
    #print('nodes_to_remove = %s' % nodes_to_remove)
    for nid in nodes_to_remove:
        del model.nodes[nid]
    model.write_bdf(bdf_filename_out)


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
        kdt = scipy.spatial.cKDTree(nodes_xyz)
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
    #replacer = unique(ieq[slots])

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
    """
    Removes the elements on one side of a model.

    Typically aircraft are defined as x-aft, y-right, z-up, so
    all node xyz locations are positive.  We then have a xz plane
    of symmetry with the axis of symmetry being y.

    Considers
    =========
      - nodes
      - elements
      - loads (LOAD/PLOAD4)

    ..note doesn't support +cuts (e.g. +x, +y, +z)
    """
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


def bdf_renumber(bdf_filename, bdf_filename_out, size=8, is_double=False):
    """
    Supports
    ========
     - GRIDs
       - superelements?
     - COORDx
     - elements
        - CELASx/CONROD/CBAR/CBEAM/CQUAD4/CTRIA3/CTETRA/CPENTA/CHEXA
     - properties
        - PSHELL/PCOMP/PCOMPG/PSOLID/PSHEAR/PBAR/PBARL
          PROD/PTUBE/PBEAM
     - mass
        - CMASSx/CONMx/PMASS

     - materials
     - loads should be updated...mabye not all...
     - RBAR/RBAR1/RBE1/RBE2/RBE3...not SUPORT/SUPORT1

    Not Done
    ========
     - SPOINT
     - any cards with SPOINTs
       - DMIG/DMI/DMIJ/DMIJI/DMIK/etc.
       - CELASx
       - CDAMPx
     - superelements
     - SPC/SPC1/SPCADD renumbering (nodes are supported)
     - SPCD/SPCAX
     - MPC/MPCADD renumbering (nodes are supported)
     - aero cards
       - FLFACT
       - CAEROx
       - PAEROx
       - SPLINEx
     - thermal cards?
     - optimization cards
     - case control
     - SETx

    ..warning:: spoints might be problematic...check
    ..warning:: still in development
    """
    model = BDF()
    model.read_bdf(bdf_filename)

    nid = 1
    if model.spoints is None:
        spoints = []
    else:
        spoints = list(model.spoints.spoints)
    nids = model.nodes.keys()

    spoints_nids = spoints + nids
    spoints_nids.sort()
    i = 1
    nnodes = len(spoints_nids)
    nid_map = {}
    reverse_nid_map = {}

    i = 1
    j = 1
    #print(spoints_nids)
    k = 0

    i = 1
    #banned_nodes = spoints
    for nid in spoints_nids:
        if nid in spoints:
            pass
            #print('sid=%s -> %s' % (nid, i))
            #i += 1
        else:
            while i in spoints:
                #print('*bump')
                i += 1
            #print('nid=%s -> %s' % (nid, i))
            nid_map[nid] = i
            reverse_nid_map[i] = nid
            i += 1
    #for nid in sorted(nids):
        #nid_map[nid] = i
        #reverse_nid_map[i] = nid
        #i += 1
    #print(nid_map)
    #print(reverse_nid_map)
    mids = (
        model.materials.keys() +
        model.creepMaterials.keys() +
        model.MATT1.keys() +
        model.MATT2.keys() +
        model.MATT3.keys() +
        model.MATT4.keys() +
        model.MATT5.keys() +
        #model.MATT6.keys() +
        #model.MATT7.keys() +
        model.MATT8.keys() +
        model.MATT9.keys() +
        model.MATS1.keys() +
        model.MATS3.keys() +
        model.MATS8.keys()
    )
    all_materials = (
        model.materials,
        model.creepMaterials,
        model.MATT1,
        model.MATT2,
        model.MATT3,
        model.MATT4,
        model.MATT5,
        #model.MATT6,
        #model.MATT7,
        model.MATT8,
        model.MATT9,
        model.MATS1,
        model.MATS3,
        model.MATS8,
    )
    mids = unique(mids)
    mids.sort()
    nmaterials = len(mids)

    mid_map = {}
    for i in range(nmaterials):
        mid = mids[i]
        mid_map[mid] = i + 1

    #spoints2 = arange(1, len(spoints) + 1)
    for nid, node in iteritems(model.nodes):
        nid_new = nid_map[nid]
        node.nid = nid_new

    cid = 0
    cid_map = {}
    for cidi, coord in iteritems(model.coords):
        coord.cid = cid
        cid_map[cidi] = cid
        cid += 1

    pid = 1
    for pidi, prop in iteritems(model.properties):
        prop.pid = pid
        pid += 1
    for pidi, prop in iteritems(model.properties_mass):
        # PMASS
        prop.pid = pid
        pid += 1

    eid = 1
    for eidi, element in iteritems(model.elements):
        element.eid = eid
        eid += 1
    for eidi, element in iteritems(model.masses):
        # CONM1, CONM2, CMASSx
        element.eid = eid
        eid += 1

    #mid = 1
    for materials in all_materials:
        for midi, material in iteritems(materials):
            mid = mid_map[midi]
            material.mid = mid

    # TODO: spc ids not updated - case control
    for spc_id, spc_group in iteritems(model.spcs):
        for i, spc in enumerate(spc_group):
            if spc.type in ['SPC', 'SPCD']:
                gids = []
                #print(spc.gids)
                for spci, gidi in enumerate(spc.gids):
                    gidii = nid_map[gidi]
                    gids.append(gidii)
                spc.gids = gids
            elif spc.type == 'SPC1':
                gids = []
                #print(spc.nodes)
                for spci, gidi in enumerate(spc.nodes):
                    gidii = nid_map[gidi]
                    gids.append(gidii)
                spc.nodes = gids
            #elif spc.type == 'SPCD':
                #raise NotImplementedError('SPCD')
            elif spc.type == 'SPCAX':
                raise NotImplementedError('SPCAX')
            else:
                raise NotImplementedError(spc.type)


    # TODO: mpc ids not updated - case control
    for mpc_id, mpc_group in iteritems(model.mpcs):
        for mpc in mpc_group:
            if mpc.type == 'MPC':
                gids = []
                for mpci, gidi in enumerate(mpc.gids):
                    gidii = nid_map[gidi]
                    gids.append(gidii)
                mpc.gids = gids
            else:
                raise NotImplementedError(mpci.type)

    for eidi, elem in iteritems(model.rigidElements):
        #$RBAR, EID, GA, GB, CNA
        #RBAR           5       1       2  123456
        elem.eid = eid
        if elem.type in ['RBAR', 'RBAR1']:
            #print('ga=%s gb=%s' % (elem.ga, elem.gb))
            elem.ga = nid_map[elem.ga]
            elem.gb = nid_map[elem.gb]
        elif elem.type == 'RBE1':
            gids = []
            for gidi in elem.Gmi:
                gidii = nid_map[gidi]
                gids.append(gidii)
            elem.Gmi = gids

            gids = []
            for gidi in elem.Gni:
                gidii = nid_map[gidi]
                gids.append(gidii)
            elem.Gni = gids
        elif elem.type == 'RBE3':
            elem.refgrid = nid_map[elem.refgrid]
            for i, (wt, ci, Gij) in enumerate(elem.WtCG_groups):
                gids = []
                for gidi in Gij:
                    gidii = nid_map[gidi]
                    gids.append(gidii)
                elem.WtCG_groups[i][2] = gids

            gids = []
            for gidi in elem.Gmi:
                gidii = nid_map[gidi]
                gids.append(gidii)
            elem.Gmi = gids
        elif elem.type == 'RBE2':
            #RBE2        3690    3389  123456    3387    3386    3385    3388       9+
            #+             10       2       1       6       5
            elem.gn = nid_map[elem.gn]
            gids = []
            for gidi in elem.Gmi:
                gidii = nid_map[gidi]
                gids.append(gidii)
            elem.Gmi = gids

            #raise NotImplementedError('RBE2')
        else:
            raise NotImplementedError(elem.type)
        eid += 1

    model.write_bdf(bdf_filename_out, size=8, is_double=False,
                    interspersed=False)

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

    bdf_filename_out = 'unique.bdf'
    tol = 0.2
    bdf_equivalence_nodes(bdf_filename, bdf_filename_out, tol)

if __name__ == '__main__':
    main()
