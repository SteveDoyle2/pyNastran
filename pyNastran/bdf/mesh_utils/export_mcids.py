"""
Defines:
 - nodes, bars = export_mcids(bdf_filename, csv_filename=None)
"""
from __future__ import print_function
from six import integer_types
import numpy as np
from pyNastran.bdf.bdf import BDF, read_bdf


def export_mcids(bdf_filename, csv_filename=None,
                 eids=None, export_xaxis=True, export_yaxis=True,
                 iply=0, log=None, debug=False):
    """
    Exports the element material coordinates systems, so you can
    load it into pyNastranGUI.

    Parameters
    ----------
    bdf_filename : str/BDF
        a bdf filename or BDF model
    csv_filename : str; default=None
        str : the path to the output csv
        None : don't write a CSV
    eids : List[int]
        the element ids to consider
    export_xaxis : bool; default=True
        export the x-axis
    export_yaxis : bool; default=True
        export the x-axis
    iply : int; default=0
        TODO: not validated
        the ply to consider

        **PSHELL**

        iply   location
        ----   --------
         0      mid1 or mid2
         1      mid1
         2      mid2
         3      mid3
         4      mid4

        **PCOMP/PCOMPG**

        iply   location
        ----   --------
        0      layer1
        1      layer2

    Returns
    -------
    nodes : (nnodes, 3) float list
        the nodes
    bars : (nbars, 2) int list
        the "bars" that represent the x/y axes of the coordinate systems

    """
    if isinstance(eids, integer_types):
        eids = [eids]


    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename, xref=False, log=log, debug=debug)
        #print(model.get_bdf_stats())
        model.safe_cross_reference()

    skip_types = [
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'CELAS5',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
        'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CGAP',
        'CROD', 'CONROD',
        'CBAR', 'CBEAM',
        'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
        'CRAC2D', 'CRAC3D',
    ]

    eid = 1
    nid = 1
    nodes = []
    bars = []
    consider_property_rotation = True  # not tested
    export_both_axes = export_xaxis and export_yaxis
    assert export_xaxis or export_yaxis

    if eids is None:
        elements = model.elements
    else:
        elements = {eid : model.elements[eid] for eid in eids}

    pids_failed = set()
    for unused_eidi, elem in sorted(elements.items()):
        if elem.type in ['CQUAD4', 'CQUAD8', 'CQUAD']:
            nid, eid = _export_quad(elem, nodes,
                                    iply, nid, eid,
                                    pids_failed, bars,
                                    export_both_axes, export_xaxis,
                                    consider_property_rotation)
        elif elem.type in ['CTRIA3', 'CTRIA6']:
            nid, eid = _export_tria(elem, nodes,
                                    iply, nid, eid,
                                    pids_failed, bars,
                                    export_both_axes, export_xaxis,
                                    consider_property_rotation)
        elif elem.type in skip_types:
            pass
        else:
            msg = 'element type=%r is not supported\n%s' % (elem.type, elem)
            raise NotImplementedError(msg)

    if len(nodes) == 0 and pids_failed:
        msg = 'No material coordinate systems could be found for iply=%s\n' % iply
        pids_failed_list = list(pids_failed)
        pids_failed_list.sort()
        model.log.warning('pids_failed_list = %s' % str(pids_failed_list))
        pid_str = [str(pid) for pid in pids_failed_list]
        msg += 'iPly=%r; Property IDs failed: [%s]\n' % (iply, ', '.join(pid_str))

        for pid in pids_failed_list:
            prop = model.properties[pid]
            msg += 'Property %s:\n%s\n' % (pid, prop)
        raise RuntimeError(msg)
    _export_coord_axes(nodes, bars, csv_filename)
    return nodes, bars

def _export_quad(elem, nodes,
                 iply, nid, eid,
                 pids_failed, bars,
                 export_both_axes, export_xaxis,
                 consider_property_rotation):
    pid_ref = elem.pid_ref
    if pid_ref.type == 'PSHELL':
        mids = [mat.type for mat in pid_ref.materials()
                if mat is not None and mat.mid > 0]
        if 'MAT8' not in mids:
            return nid, eid
    elif pid_ref.type in ['PCOMP', 'PCOMPG']:
        pass
    else:
        raise NotImplementedError(pid_ref)
    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    xyz3 = elem.nodes_ref[2].get_position()
    xyz4 = elem.nodes_ref[3].get_position()

    dxyz21 = np.linalg.norm(xyz2 - xyz1)
    dxyz32 = np.linalg.norm(xyz3 - xyz2)
    dxyz43 = np.linalg.norm(xyz4 - xyz3)
    dxyz14 = np.linalg.norm(xyz1 - xyz4)
    dxyz = np.mean([dxyz21, dxyz32, dxyz43, dxyz14]) / 2.
    centroid, imat, jmat, normal = elem.material_coordinate_system()

    try:
        imat, jmat = _rotate_mcid(
            elem, pid_ref, iply, imat, jmat, normal,
            consider_property_rotation=consider_property_rotation)
    except IndexError:
        pids_failed.add(pid_ref.pid)
        return nid, eid

    iaxis = centroid + imat * dxyz
    jaxis = centroid + jmat * dxyz
    nid, eid = _add_elements(nid, eid, nodes, bars,
                             centroid, iaxis, jaxis,
                             export_both_axes, export_xaxis)

    #nid, eid = _export_quad(elem, bars, nid, eid, export_both_axes, export_xaxis)
    return nid, eid

def _export_tria(elem, nodes,
                 iply, nid, eid,
                 pids_failed, bars,
                 export_both_axes, export_xaxis,
                 consider_property_rotation):
    """helper method for ``export_mcids``"""
    pid_ref = elem.pid_ref
    if pid_ref.type == 'PSHELL':
        mids = [mat.type for mat in pid_ref.materials()
                if mat is not None and mat.mid > 0]
        if 'MAT8' not in mids:
            return nid, eid
    elif pid_ref.type in ['PCOMP', 'PCOMPG']:
        pass
    else:
        raise NotImplementedError(pid_ref)

    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    xyz3 = elem.nodes_ref[2].get_position()

    dxyz21 = np.linalg.norm(xyz2 - xyz1)
    dxyz32 = np.linalg.norm(xyz3 - xyz2)
    dxyz13 = np.linalg.norm(xyz1 - xyz3)
    dxyz = np.mean([dxyz21, dxyz32, dxyz13]) / 2.
    centroid, imat, jmat, normal = elem.material_coordinate_system()

    try:
        imat, jmat = _rotate_mcid(
            elem, pid_ref, iply, imat, jmat, normal,
            consider_property_rotation=consider_property_rotation)
    except IndexError:
        pids_failed.add(pid_ref.pid)
        return nid, eid

    iaxis = centroid + imat * dxyz
    jaxis = centroid + jmat * dxyz
    nid, eid = _add_elements(nid, eid, nodes, bars,
                             centroid, iaxis, jaxis,
                             export_both_axes, export_xaxis)
    return nid, eid

def _export_coord_axes(nodes, bars, csv_filename):
    """save the coordinate systems in a csv file"""
    if csv_filename:
        with open(csv_filename, 'w') as out_file:
            for node in nodes:
                out_file.write('GRID,%i,%s,%s,%s\n' % node)
            for bari in bars:
                out_file.write('BAR,%i,%i,%i\n' % bari)

def _rotate_mcid(elem, pid_ref, iply, imat, jmat, normal,
                 consider_property_rotation=True):
    """rotates a material coordinate system"""
    if not consider_property_rotation:
        return imat, jmat

    pid_ref = elem.pid_ref
    if pid_ref.type == 'PSHELL':
        theta_mcid = elem.theta_mcid
        if isinstance(theta_mcid, float):
            thetad = theta_mcid
        elif isinstance(theta_mcid, integer_types):
            return imat, jmat
        else:
            msg = 'theta/mcid=%r is not an int/float; type=%s\n%s' % (
                theta_mcid, type(theta_mcid), elem)
            raise TypeError(msg)
    elif pid_ref.type in ['PCOMP', 'PCOMPG']:
        #print(pid_ref.nplies)
        thetad = pid_ref.get_theta(iply)
    else:
        msg = 'property type=%r is not supported\n%s%s' % (
            elem.pid_ref.type, elem, elem.pid_ref)
        raise NotImplementedError(msg)
    if isinstance(thetad, float) and thetad == 0.0:
        return imat, jmat

    theta = np.radians(thetad)
    cos = np.cos(theta)
    sin = np.sin(theta)

    theta_rotation = np.array([
        [cos, -sin, 0.],
        [sin, cos, 0.],
        [0., 0., 1.]
        ], dtype='float64')

    element_axes = np.vstack([imat, jmat, normal])
    rotated_axes = np.dot(theta_rotation, element_axes)  ## TODO: validate
    imat2 = rotated_axes[0, :]
    jmat2 = rotated_axes[1, :]
    return imat2, jmat2

def _add_elements(nid, eid, nodes, bars,
                  centroid, iaxis, jaxis,
                  export_both_axes, export_xaxis):
    """adds the element data"""
    if export_both_axes:
        nodes.append((nid, centroid[0], centroid[1], centroid[2]))
        nodes.append((nid + 1, iaxis[0], iaxis[1], iaxis[2]))
        nodes.append((nid + 2, jaxis[0], jaxis[1], jaxis[2]))
        bars.append((eid, nid, nid + 1))  # x-axis
        bars.append((eid + 1, nid, nid + 2))  # y-axis
        nid += 3
        eid += 2
    elif export_xaxis:
        nodes.append((nid, centroid[0], centroid[1], centroid[2]))
        nodes.append((nid + 1, iaxis[0], iaxis[1], iaxis[2]))
        bars.append((eid, nid, nid + 1))  # x-axis
        nid += 2
        eid += 1
    else:
        # export_yaxis
        nodes.append((nid, centroid[0], centroid[1], centroid[2]))
        nodes.append((nid + 1, jaxis[0], jaxis[1], jaxis[2]))
        bars.append((eid, nid, nid + 1))  # y-axis
        nid += 2
        eid += 1
    return nid, eid
