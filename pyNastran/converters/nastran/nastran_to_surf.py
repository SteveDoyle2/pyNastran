"""
defines:
 - clear_out_solids(bdf_filename, bdf_filename_out=None,
                    equivalence=True, renumber=True, equivalence_tol=0.01)
 - nastran_to_surf(bdf_filename, pid_to_element_flags, surf_filename,
                   renumber_pids=None, line_map=None,
                   scale=1.0, tol=1e-10,
                   xref=True)

"""
from collections import defaultdict
from numpy import array, allclose, unique, zeros
from pyNastran.bdf.bdf import read_bdf

from pyNastran.bdf.mesh_utils.bdf_equivalence import bdf_equivalence_nodes
from pyNastran.bdf.mesh_utils.bdf_renumber import bdf_renumber
from pyNastran.bdf.mesh_utils.remove_unused import remove_unused

def remove_unassociated_nodes(bdf_filename, unused_bdf_filename_out, renumber=False):
    """dummy function"""
    assert renumber is False, renumber
    remove_unused(bdf_filename, remove_nids=True, remove_cids=False,
                  remove_pids=False, remove_mids=False)

def clear_out_solids(bdf_filename, bdf_filename_out=None,
                     equivalence=True, renumber=True, equivalence_tol=0.01):
    """removes solid elements"""
    if bdf_filename_out is None:
        if renumber or equivalence:
            msg = ('bdf_filename_out=%s must be specified if renumber=%s '
                   'or equivalence=%s are True' % (
                       bdf_filename_out, renumber, equivalence))
            raise RuntimeError(msg)

    if isinstance(bdf_filename, str):
        print('clearing out solids from %s' % bdf_filename)
        model = read_bdf(bdf_filename, xref=False)
    else:
        model = bdf_filename
    #nodes2    = {nid, node for nid, node in model.nodes.items()}
    #elements2 = {eid, element for eid, element in model.elements.items()
                 #if element.type in ['CTRIA3', 'CQUAD4']}

    out_dict = model.get_card_ids_by_card_types(card_types=['CTRIA3', 'CQUAD4'])
    save_eids = set(out_dict['CTRIA3'] + out_dict['CQUAD4'])
    all_eids = set(model.element_ids)
    #print('all_eids =', all_eids)
    #print('save_eids =', save_eids)
    remove_eids = all_eids - save_eids
    #print('remove_eids =', remove_eids)

    for eid in remove_eids:
        #print('eid =', eid)
        del model.elements[eid]

    # TODO: seems like we could be more efficient...
    #nids = unique(hstack([model.elements[eid].node_ids for eid in save_eids]))

    # get nodes that are remaining in the model
    nids = set()
    unused_elements2 = {}
    #print(model.elements)
    for eid, element in model.elements.items():
        #if element.type not in ['CTRIA3', 'CQUAD4']:
            #continue
        #elements2[eid] = element
        nids.update(element.node_ids)
    nids = list(nids)
    nids.sort()

    # filter out old nodes & properties
    nodes2 = {nid : node for nid, node in model.nodes.items() if nid in nids}
    properties2 = {pid : prop for pid, prop in model.properties.items() if prop.type == 'PSHELL'}

    model.nodes = nodes2
    #model.elements = elements2
    model.properties = properties2

    # already equivalenced?
    #remove_unassociated_nodes(bdf_filename, bdf_filename_out, renumber=False)

    #bdf_filename_out = 'equivalence.bdf'
    starting_id_dict = {
        'cid' : 1,
        'nid' : 1,
        'eid' : 1,
        'pid' : 1,
        'mid' : 1,
    }
    if equivalence:
        if renumber:
            bdf_equivalenced_filename = 'equivalence.bdf'
        else:
            bdf_equivalenced_filename = bdf_filename_out

        model.write_bdf('remove_unused_nodes.bdf')
        bdf_equivalence_nodes(model, bdf_equivalenced_filename, equivalence_tol,
                              renumber_nodes=False, neq_max=4, xref=True)
        if renumber:
            bdf_renumber(bdf_equivalenced_filename, bdf_filename_out, size=8, is_double=False,
                         starting_id_dict=starting_id_dict)
    elif renumber:
        model.cross_reference()
        bdf_renumber(model, bdf_filename_out, size=8, is_double=False,
                     starting_id_dict=starting_id_dict)

    return model

def nastran_to_surf(bdf_filename, pid_to_element_flags, surf_filename,
                    renumber_pids=None, line_map=None,
                    scale=1.0, tol=1e-10,
                    xref=True):
    """
    Converts a BDF to an AFLR3 surf file

    Parameters
    ----------
    bdf_filename : str/BDF
        str : the input BDF filename
        BDF : a BDF model that has been cross-referenced
    surf_filename : str
        the output SURF filename
    pid_to_element_flags : dict[key] = value
        key=PSHELL value=[layer0 thickness, BL_thickness, grid_bc]
    renumber_pids : bool; default=None
        a mapping of pid to surface ID
        None = no remapping
    line_map : dict[key] = value
        same as pid_to_element_flags, but for the specific intersections
        where there are BC differences
        NOTE:  we only check [thickness, BL_thickness] because we're
        averaging this data for the nodes
    scale : float; default=1.0
        scales the mesh by scale for unit conversion
    tol : float; default=1e-16
        I hate 1e-16 values in my model
    xref : bool; default=True
        does the model need to be cross-referenced to calculate the
        node positions?

    # these pids correspond to the BDF
    pid_to_element_flags = {
        1 : [wall_initial_normal_spacing,           wall_bl_thickness, grid_bc], # top_wall
        2 : [side_wall_initial_normal_spacing, side_wall_bl_thickness, grid_bc], # right_wall
        3 : [side_wall_initial_normal_spacing, side_wall_bl_thickness, grid_bc], # left_wall
        4 : [side_wall_initial_normal_spacing, side_wall_bl_thickness, grid_bc], # outlet

        5 : [far_field_initial_normal_spacing, far_field_bl_thickness, grid_bc], # bottom_wall
        6 : [wall_initial_normal_spacing,           wall_bl_thickness, grid_bc], # bay

        11 : [far_field_initial_normal_spacing, far_field_bl_thickness, grid_bc],  # inlet_btm
        12 : [side_wall_initial_normal_spacing, side_wall_bl_thickness, grid_bc],  # inlet_front
        13 : [side_wall_initial_normal_spacing, side_wall_bl_thickness, grid_bc],  # inlet_left
        14 : [side_wall_initial_normal_spacing, side_wall_bl_thickness, grid_bc],  # inlet_right
        15 : [wall_initial_normal_spacing,           wall_bl_thickness, grid_bc],  # inlet_visc
    }

    # these pids correspond to the BDF
    # the pid_to_element_flag at the intersection between pids
    line_map = {
        (1, 2) : [wall_initial_normal_spacing, wall_bl_thickness, grid_bc],
        (1, 3) : [wall_initial_normal_spacing, wall_bl_thickness, grid_bc],
        (1, 4) : [wall_initial_normal_spacing, wall_bl_thickness, grid_bc],
    }

    # these are done at the last step to make the output "nice"
    renumber_pids = {
        11 : 7,
        12 : 8,
        13 : 9,
        14 : 10,
        15 : 11,
    }
    scale = 0.0254  # inches to meters;
    """
    if renumber_pids is None:
        renumber_pids = {}
    if line_map is None:
        line_map = {}

    if isinstance(bdf_filename, str):
        model = read_bdf(bdf_filename, xref=xref)
    else:
        model = bdf_filename

    unused_nnodes = len(model.nodes)
    nodes = []
    quads = []
    tris = []

    maxnode, nodes, node_flags, node_flags_temp = _get_nodes(model, scale, xref)

    node_remaps = {}
    #if 0:
        #xyz_array = array(nodes, dtype='float64')
        #for nid, xyz in enumerate(xyz_array):
            #for nidi, xyz2 in enumerate(xyz_array[nid+1:, :]):
                #nid2 = nid + nidi + 1
                #if not allclose(nid + 1, nid2 + 1):
                    #msg = 'nid=%s nid2=%s xyz=%s' % (nid+1, nid2+1, xyz)
                    #raise RuntimeError(msg)
                #if allclose(xyz, xyz2):
                    ##print(nid, nid2, nidi)
                    ##if nid + 1 in node_remaps:

                    #node_remaps[nid2 + 1] = nid + 1
                    #print('nid=%s nid2=%s xyz=%s xyz2=%s' % (nid+1, nid2+1, xyz, xyz2))
                #assert not(allclose(xyz, xyz2)), 'nid=%s nid2=%s xyz=%s' % (nid+1, nid2+1, xyz)
        #del xyz_array

    pid0 = 1
    for pid, prop in sorted(model.properties.items()):
        if pid != pid0:
            msg = 'properties must go from 1 to N, no gaps; pid=%s expected=%s' % (
                pid, pid0)
            raise RuntimeError(msg)
        #assert pid in pid_to_element_flags, pid
        if prop.type in ['PSOLID']:
            continue
        if prop.type not in ['PSHELL', 'PCOMP', 'PCOMPG']:
            raise NotImplementedError(prop)
        pid0 += 1

    nid_to_eid_map = get_nid_to_eid_map(
        model,
        node_flags_temp, pid_to_element_flags, node_remaps,
        tris, quads)

    initial_normal_spacing0 = 0
    bl_thickness0 = 0
    for nid, node_flagsi in node_flags_temp.items():
        nodes_flags_array = array(node_flagsi)  # (N, 2)
        nflags = nodes_flags_array.shape[0]
        if nflags == 0:
            #node_flags[nid] = [initial_normal_spacing0, bl_thickness0]
            continue
        try:
            avg_node_flagsi = nodes_flags_array.mean(axis=0)
            max_node_flagsi = nodes_flags_array.max(axis=0)
        except ValueError:
            print('nid=%s node_flagsi=%s' % (nid, node_flagsi))
            raise RuntimeError('node %i is duplicated (equivalence your nodes)'
                               ' or you have unused nodes' % nid)

        if not allclose(avg_node_flagsi, max_node_flagsi):
            eidsi = unique(nid_to_eid_map[nid])
            pidsi = unique([model.elements[eid].Pid() for eid in eidsi])
            pidsi.sort()
            pidsi = tuple(pidsi)
            if pidsi in line_map:
                element_flag = line_map[pidsi]
                unused_name, spacing, thickness, unused_grid_bc = element_flag
                avg_node_flagsi = [spacing, thickness]
            else:
                msg = ('\nERROR BL THICKNESS MISMATCH:\n  define a  line_map to resolve for nid=%s;'
                       ' map=%s; pids=%s\n' % (
                           nid, nid_to_eid_map[nid], pidsi))
                for pid in pidsi:
                    msg += '   pid=%s name=%s\n' % (pid, pid_to_element_flags[pid][0])
                msg += "  nodes_flags_array = \n%s\n" % nodes_flags_array
                avg_node_flagsi = nodes_flags_array.max(axis=0)
                msg += '  FOUND:    avg_node_flag =\n%s\n'  % (avg_node_flagsi)
                msg += '  ASSUMING: max_node_flags =\n%s'  % (max_node_flagsi)
                #raise NotImplementedError(msg)
                print(msg)
                avg_node_flagsi = max_node_flagsi
                #avg_node_flagsi = avg_node_flagsi.max()

        if len(avg_node_flagsi) != 2:
            msg = 'len([normal_spacing, bl_thickness])=2; actual=%s' % len(avg_node_flagsi)
            raise RuntimeError(msg)

        assert nid > 0, nid
        node_flags[nid] = avg_node_flagsi

    _write_surf(surf_filename, maxnode,
                nodes, tris, quads,
                node_flags, initial_normal_spacing0, pid_to_element_flags, bl_thickness0,
                renumber_pids, tol)

def _get_nodes(model, scale, xref):
    """helper method for ``nastran_to_surf``"""
    # assume nodes go from 1:#
    #nid0 = 1
    node_flags = {}
    node_flags_temp = {}

    maxnode = max(model.nodes.keys())
    nodes = zeros((maxnode, 3), dtype='float64')
    if xref:
        for nid, node in sorted(model.nodes.items()):
            #if nid != nid0:
                #msg = 'nodes must go from 1 to N, no gaps; nid=%s expected=%s' % (nid, nid0)
                #raise RuntimeError(msg)
            xyz = node.get_position()
            nodes[nid-1] = xyz * scale
            node_flags[nid] = []
            node_flags_temp[nid] = []
            #nid0 += 1
    else:
        for nid, node in sorted(model.nodes.items()):
            #if nid != nid0:
                #msg = 'nodes must go from 1 to N, no gaps; nid=%s expected=%s' % (nid, nid0)
                #raise RuntimeError(msg)
            xyz = node.xyz
            nodes[nid-1] = xyz * scale
            node_flags[nid] = []
            node_flags_temp[nid] = []
            #nid0 += 1
    return maxnode, nodes, node_flags, node_flags_temp

def get_nid_to_eid_map(model,
                       node_flags_temp, pid_to_element_flags, node_remaps,
                       tris, quads):
    """helper method for ``nastran_to_surf``"""
    nid_to_eid_map = defaultdict(list)
    for eid, element in sorted(model.elements.items()):
        #if element.type not in ['CQUAD4', 'CTRIA3']:
            #continue
        nids = element.node_ids
        nids2 = []
        for nid in nids:
            nid_to_eid_map[nid].append(eid)
            if nid in node_remaps:
                nid = node_remaps[nid]
            nids2.append(nid)
        pid = element.Pid()

        element_flag = pid_to_element_flags[pid]
        unused_name, spacing, thickness, unused_grid_bc = element_flag
        element_flag_node = [spacing, thickness]

        for nid in nids:
            try:
                node_flags_temp[nid].append(element_flag_node)
            except KeyError:
                print('max_nid=%s max_temp_nid=%s' % (max(nids), max(node_flags_temp)))
                print("nids = %s" % nids)
                print("node_flags_temp = %s" % node_flags_temp.keys())
                raise
            assert nid > 0, element

        if element.type == 'CTRIA3':
            tris.append([nids2, pid])
        elif element.type == 'CQUAD4':
            quads.append([nids2, pid])
        else:
            raise NotImplementedError(element)
    #del pid, nids
    return nid_to_eid_map

def _write_surf(surf_filename, maxnode,
                nodes, tris, quads,
                node_flags, initial_normal_spacing0, pid_to_element_flags, bl_thickness0,
                renumber_pids, tol):
    """writes the actual surf file"""
    ntris = len(tris)
    nquads = len(quads)
    assert ntris + nquads > 0, 'nelements=%s' % (ntris + nquads)

    with open(surf_filename, 'w', encoding='ascii') as surf_file:
        #surf_file.write('ntris nquads nnodes\n')
        surf_file.write('%i %i %i\n' % (ntris, nquads, maxnode))

        # writing nodes
        for nid, node in enumerate(nodes):
            x, y, z = node
            if abs(x) < tol:
                x = 0.0
            if abs(y) < tol:
                y = 0.0
            if abs(z) < tol:
                z = 0.0

            try:
                initial_normal_spacing, bl_thickness = node_flags[nid + 1]
            except KeyError:
                initial_normal_spacing = initial_normal_spacing0
                bl_thickness = bl_thickness0
            except ValueError:
                initial_normal_spacing = initial_normal_spacing0
                bl_thickness = bl_thickness0
            surf_file.write('%.10e %.10e %.10e %g %g\n' % (
                x, y, z, initial_normal_spacing, bl_thickness))

        # writing triangles
        rf = 0  # recon flag is used for tris and is super cobnfusing; quads are better
        for unused_eid, (element, pid) in enumerate(tris):
            unused_name, initial_normal_spacing, bl_thickness, grid_bc = pid_to_element_flags[pid]
            if pid in renumber_pids:
                pid = renumber_pids[pid]
            n1, n2, n3 = element
            surf_file.write('%i %i %i %i %s %s\n' % (n1, n2, n3, pid, rf, grid_bc))

        # writing quads
        rf = 0
        for unused_eid, (element, pid) in enumerate(quads):
            unused_name, initial_normal_spacing, bl_thickness, grid_bc = pid_to_element_flags[pid]
            if pid in renumber_pids:
                pid = renumber_pids[pid]
            n1, n2, n3, n4 = element
            surf_file.write('%i %i %i %i %i %s %s\n' % (n1, n2, n3, n4, pid, rf, grid_bc))
