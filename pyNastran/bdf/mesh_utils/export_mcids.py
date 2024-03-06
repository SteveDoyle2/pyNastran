"""
Defines:
 - nodes, bars = export_mcids(bdf_filename, csv_filename=None)

"""
from __future__ import annotations
from collections import defaultdict
from typing import Union, Optional, TYPE_CHECKING
import numpy as np
from pyNastran.bdf.bdf import BDF, read_bdf, PCOMP, PCOMPG, PSHELL
from pyNastran.bdf.cards.elements.shell import (
    CTRIA3, CTRIA6, CQUAD4, CQUAD8, CTRIAR, CQUADR, rotate_by_thetad)
ShellElement = Union[CTRIA3, CTRIA6, CQUAD4, CQUAD8, CTRIAR, CQUADR]
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger

SKIP_ETYPES = {
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'CELAS5',
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CGAP', 'CVISC', 'CFAST',
    'CROD', 'CONROD', 'CTUBE',
    'CBAR', 'CBEAM', 'CBEND', 'CBEAM3',
    'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
    'CRAC2D', 'CRAC3D',
    'CSHEAR',
    'CHBDYE', 'CHBDYP', 'CHBDYG', 'GENEL',
    'CQUADX', 'CTRIAX', 'CQUADX4', 'CQUADX8', 'CTRIAX6',
    'CTRAX3', 'CTRAX6',
    'CPLSTN3', 'CPLSTN4', 'CPLSTN8', 'CPLSTN6',
    'CPLSTS3', 'CPLSTS4', 'CPLSTS8', 'CPLSTS6',
    'CHACAB', 'CAABSF',
}

def export_mcids(bdf_filename: Union[BDF, str], csv_filename: Optional[str]=None,
                 eids: Optional[list[int]]=None,
                 export_xaxis: bool=True, export_yaxis: bool=True,
                 consider_property_rotation: bool=True,
                 iply: int=0, log=None, debug=False):
    """
    Exports the element material coordinates systems for non-isotropic
    materials.

    Parameters
    ----------
    bdf_filename : str/BDF
        a bdf filename or BDF model
    csv_filename : str; default=None
        str : the path to the output csv
        None : don't write a CSV
    eids : list[int]
        the element ids to consider
    export_xaxis : bool; default=True
        export the x-axis
    export_yaxis : bool; default=True
        export the x-axis
    consider_property_rotation : bool; default=True
        rotate the coordinate system
    iply : int; default=0
        TODO: not validated
        the ply to consider
    pid_to_nplies : dict[int pid, int nplies]; default=None -> auto
        optional dictionary to speed up analysis

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
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename, xref=False, log=log, debug=debug)
        #print(model.get_bdf_stats())
        model.safe_cross_reference()

    elements = _get_elements(model, eids)
    pid_to_nplies, nplies_max = get_pid_to_nplies(model)
    if nplies_max == 0:
        return {}, 0
    if iply >= nplies_max:
        raise RuntimeError(f'no ply {iply} found')

    eid = 1
    nid = 1
    nodes = []
    bars = []
    consider_property_rotation = True  # not tested
    export_both_axes = export_xaxis and export_yaxis
    assert export_xaxis or export_yaxis

    #pids_failed = set()
    for unused_eidi, elem in sorted(elements.items()):
        if elem.type in {'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR'}:
            pid = elem.pid
            nplies = pid_to_nplies[pid]
            if iply >= nplies:
                continue
            nid, eid = _export_quad_mcid(model, elem, nodes,
                                         iply, nid, eid,
                                         #pids_failed,
                                         bars,
                                         export_both_axes, export_xaxis,
                                         consider_property_rotation)
        elif elem.type in {'CTRIA3', 'CTRIA6', 'CTRIAR'}:
            pid = elem.pid
            nplies = pid_to_nplies[pid]
            if iply >= nplies:
                continue
            nid, eid = _export_tria_mcid(model, elem, nodes,
                                         iply, nid, eid,
                                         #pids_failed,
                                         bars,
                                         export_both_axes, export_xaxis,
                                         consider_property_rotation)
        elif elem.type in SKIP_ETYPES:
            continue
        else:
            raise NotImplementedError(f'element type={elem.type!r} is not supported\n{elem}')

    #if len(nodes) == 0 and pids_failed:
        #msg = 'No material coordinate systems could be found for iply=%s\n' % iply
        #pids_failed_list = list(pids_failed)
        #pids_failed_list.sort()
        #model.log.warning('pids_failed_list = %s' % str(pids_failed_list))
        #pid_str = [str(pid) for pid in pids_failed_list]
        #msg += 'iPly=%r; Property IDs failed: [%s]\n' % (iply, ', '.join(pid_str))

        #for pid in pids_failed_list:
            #prop = model.properties[pid]
            #msg += f'Property {pid}:\n{prop}\n'
        #raise RuntimeError(msg)
    _export_coord_axes(nodes, bars, csv_filename)
    return nodes, bars

def _get_elements(model, eids):
    if isinstance(eids, int):
        eids = [eids]

    if eids is None:
        elements = model.elements
    else:
        elements = {eid : model.elements[eid] for eid in eids}
    return elements

def export_element_cid(bdf_filename: Union[BDF, str],
                       eids: Optional[list[int]]=None,
                       log=None, debug=False):
    """
    Exports the element coordinates systems for non-isotropic materials.

    Note that for two quads identically oriented/numbered PSHELL quads
    with theta different between the two, the element cid will be same.
    """
    model = _load_bdf(bdf_filename, log=log, debug=debug)
    elements = _get_elements(model, eids)
    for unused_eidi, elem in sorted(elements.items()):
        if elem.type in {'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR'}:
            #pid = elem.pid
            x = 1
        elif elem.type in {'CTRIA3', 'CTRIA6', 'CTRIAR'}:
            pid = elem.pid
            #nplies = pid_to_nplies[pid]
            #if nplies == 0:
                #continue
            y = 1
            #_export_tri_all(
                #model, elem, nplies,
                #iply_to_nids,
                #iply_to_nodes,
                #iply_to_bars)
        elif elem.type in SKIP_ETYPES:
            pass
        else:
            raise NotImplementedError(f'element type={elem.type!r} is not supported\n{elem}')
    #print(iply_to_nids)
    #_export_coord_axes(nodes, bars, csv_filename)
    return nodes, bars


def _load_bdf(bdf_filename: Union[BDF, str],
              log: Optional[SimpleLogger]=None,
              debug: bool=True) -> BDF:
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        model = read_bdf(bdf_filename, xref=False, log=log, debug=debug)
        #print(model.get_bdf_stats())
        model.safe_cross_reference()
    return model

def export_mcids_all(bdf_filename: Union[BDF, str],
                     eids: Optional[list[int]]=None,
                     log: Optional[SimpleLogger]=None, debug: bool=False):
    """
    Exports the element material coordinates systems for non-isotropic
    materials.

    Note that for two quads identically oriented/numbered PSHELL quads
    with theta different between the two, the mcid will be different.

    Parameters
    ----------
    bdf_filename : str/BDF
        a bdf filename or BDF model
    csv_filename : str; default=None
        str : the path to the output csv
        None : don't write a CSV
    eids : list[int]
        the element ids to consider

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
    model = _load_bdf(bdf_filename, log=None, debug=True)

    pid_to_nplies, nplies_max = get_pid_to_nplies(model)
    if nplies_max == 0:
        return {}, 0

    elements = _get_elements(model, eids)

    # -2: elment coord
    #  -1: MCID
    #  0: ply 1
    #  1: ply 2
    # ...
    iply_to_nids = {iply : 0 for iply in range(-2, nplies_max)}
    iply_to_nodes = {iply : [] for iply in range(-2, nplies_max)}
    iply_to_bars = {iply : [] for iply in range(-2, nplies_max)}

    for unused_eidi, elem in sorted(elements.items()):
        if elem.type in {'CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR'}:
            pid = elem.pid
            nplies = pid_to_nplies[pid]
            if nplies == 0:
                continue
            _export_quad_mcid_all(
                model, elem, nplies,
                iply_to_nids,
                iply_to_nodes,
                iply_to_bars)
        elif elem.type in {'CTRIA3', 'CTRIA6', 'CTRIAR'}:
            pid = elem.pid
            nplies = pid_to_nplies[pid]
            if nplies == 0:
                continue
            _export_tri_mcid_all(
                model, elem, nplies,
                iply_to_nids,
                iply_to_nodes,
                iply_to_bars)
        elif elem.type in SKIP_ETYPES:
            continue
        else:
            raise NotImplementedError(f'element type={elem.type!r} is not supported\n{elem}')
        dxyz, centroid, ielement, jelement, normal = elem.element_coordinate_system()
        iaxis = centroid + ielement * dxyz
        #iaxis = centroid + imati * dxyz
        #jaxis = centroid + jmati * dxyz
        ilayer = -2
        iply_to_nids[ilayer] = _export_xaxis(
            iply_to_nids[ilayer],
            iply_to_nodes[ilayer],
            iply_to_bars[ilayer],
            centroid, iaxis)


    #print(iply_to_nids)
    #_export_coord_axes(nodes, bars, csv_filename)
    return iply_to_nodes, iply_to_bars

def get_pid_to_nplies(model: BDF) -> tuple[dict[int, int], int]:
    pid_to_nplies = defaultdict(int)
    for pid, prop in model.properties.items():
        if prop.type in ['PCOMP', 'PCOMPG']:
            pid_to_nplies[pid] = prop.nplies
        elif prop.type == 'PSHELL':
            pid_to_nplies[pid] = 1

    all_plies = list(pid_to_nplies.values())
    if len(all_plies) == 0:
        return {}, 0
    nplies_max = max(all_plies)
    assert isinstance(nplies_max, int), nplies_max
    return pid_to_nplies, nplies_max

def get_pid_ref_prop_type(model: BDF, elem) -> tuple[Union[PCOMP, PCOMPG, PSHELL], str]:
    """helper method for ``export_mcids``"""
    pid_ref = elem.pid_ref
    try:
        prop_type = pid_ref.type
    except AttributeError:
        print(elem.get_stats())
        print(model.properties)
        raise
    assert hasattr(elem, 'theta_mcid'), elem.get_stats()
    return pid_ref, prop_type

def _export_quad_mcid(model: BDF,
                      elem, nodes,
                      iply: int, nid: int, eid: int,
                      #pids_failed: set[int],
                      bars: list[list[int]],
                      export_both_axes: bool, export_xaxis: bool,
                      consider_property_rotation: bool) -> tuple[int, int]:
    """helper method for ``export_mcids``"""
    pid_ref, prop_type = get_pid_ref_prop_type(model, elem)

    if prop_type == 'PSHELL':
        mids = [mat.type for mat in pid_ref.materials()
                if mat is not None and mat.mid > 0]
        if 'MAT8' not in mids:
            return nid, eid
    elif prop_type in ['PCOMP', 'PCOMPG']:
        pass
    #elif prop_type in ['PLPLANE']:
        #return nid, eid
    else:
        raise NotImplementedError(pid_ref)

    dxyz, centroid, imat, jmat, normal = _get_quad_vectors_mcid(elem)

    nid, eid = _rotate_single_coord(
        elem, pid_ref, iply,
        nid, eid, nodes, bars,
        dxyz, centroid, imat, jmat, normal,
        export_both_axes=export_both_axes, export_xaxis=export_xaxis,
        consider_property_rotation=consider_property_rotation)
    return nid, eid

def _export_quad_mcid_all(model: BDF,
                          elem: ShellElement,
                          nplies: int,
                          nids: dict[int, int],
                          nodes: dict[int, list[np.ndarray]],
                          bars: dict[int, list[tuple[int, int]]]) -> None:
    """helper method for ``export_mcids``"""
    pid_ref, prop_type = get_pid_ref_prop_type(model, elem)

    if prop_type == 'PSHELL':
        mids = [mat.type for mat in pid_ref.materials()
                if mat is not None and mat.mid > 0]
        if 'MAT8' not in mids:
            _make_element_coord_quad(elem, pid_ref, nids, nodes, bars)
            return
    elif prop_type in ['PCOMP', 'PCOMPG']:
        pass
    #elif prop_type in ['PLPLANE']:
        #return nid, eid
    else:
        raise NotImplementedError(pid_ref)

    dxyz, centroid, imat, jmat, normal = _get_quad_vectors_mcid(elem)
    _rotate_coords(elem, pid_ref, nplies, nids, nodes, bars, dxyz, centroid, imat, jmat, normal)

def _make_element_coord_tri(elem: ShellElement, pid_ref, nids, nodes, bars):
    dxyz, centroid, imat, jmat, normal = _get_tri_vectors_mcid(elem)
    nplies = 0 # only make the element coordinate system
    _rotate_coords(elem, pid_ref, nplies, nids, nodes, bars, dxyz, centroid, imat, jmat, normal)

def _make_element_coord_quad(elem: ShellElement, pid_ref, nids, nodes, bars):
    dxyz, centroid, imat, jmat, normal = _get_quad_vectors_mcid(elem)
    nplies = 0 # only make the element coordinate system
    _rotate_coords(elem, pid_ref, nplies, nids, nodes, bars, dxyz, centroid, imat, jmat, normal)

def _rotate_coords(elem: ShellElement, pid_ref,
                   nplies: int,
                   nids: dict[int, int],
                   nodes: dict[int, list[Any]],
                   bars: dict[int, list[Any]],
                   dxyz: float,
                   centroid: np.ndarray, imat: np.ndarray, jmat: np.ndarray,
                   normal: np.ndarray) -> None:
    """
    iply   Final Label  Description
    ====   ===========  ===========
    -1     0            material coordinate system
    0      1            ply 1
    1      2            ply 2
    2      3            ply 3

    #nplies = 3
    nids - {-1: 0, 0: 0, 1: 0, 2: 0}
    nodes = {-1: [], 0: [], 1: [], 2: []}
    bars = {-1: [], 0: [], 1: [], 2: []}

    """
    # material coordinate system
    iaxis = centroid + imat * dxyz
    nids[-1] = _export_xaxis(nids[-1], nodes[-1], bars[-1], centroid, iaxis)

    for iply in range(nplies):
        imati, unused_jmat = _rotate_mcid(
            elem, pid_ref, iply, imat, jmat, normal,
            consider_property_rotation=True)

        iaxis = centroid + imati * dxyz
        #jaxis = centroid + jmati * dxyz
        nids[iply] = _export_xaxis(nids[iply], nodes[iply], bars[iply], centroid, iaxis)

def _export_tri_mcid_all(model: BDF,
                         elem: ShellElement,
                         nplies: int,
                         nids: dict[int, int],
                         nodes: dict[int, list[np.ndarray]],
                         bars: dict[int, list[tuple[int, int]]]) -> None:
    """helper method for ``export_mcids``"""
    pid_ref, prop_type = get_pid_ref_prop_type(model, elem)

    if prop_type == 'PSHELL':
        mids = [mat.type for mat in pid_ref.materials()
                if mat is not None and mat.mid > 0]
        if 'MAT8' not in mids:
            _make_element_coord_tri(elem, pid_ref, nids, nodes, bars)
            return
    elif prop_type in ['PCOMP', 'PCOMPG']:
        pass
    #elif prop_type in ['PLPLANE']:
        #return nid, eid
    else:
        raise NotImplementedError(pid_ref)

    dxyz, centroid, imat, jmat, normal = _get_tri_vectors_mcid(elem)
    _rotate_coords(elem, pid_ref, nplies, nids, nodes, bars, dxyz, centroid,
                   imat, jmat, normal)

def _get_tri_vectors_mcid(elem: CTRIA3):
    try:
        node1, node2, node3 = elem.nodes_ref[:3]
    except (IndexError, ValueError):
        print(elem.get_stats())
        raise
    xyz1 = node1.get_position()
    xyz2 = node2.get_position()
    xyz3 = node3.get_position()

    # take the mean edge length to size the vectors in the GUI
    dxyz21 = np.linalg.norm(xyz2 - xyz1)
    dxyz32 = np.linalg.norm(xyz3 - xyz2)
    dxyz13 = np.linalg.norm(xyz1 - xyz3)
    dxyz = np.mean([dxyz21, dxyz32, dxyz13]) / 2.

    # adjusted for element coordinate system
    # equivalent to PCOMP thetad=0.0
    dxyz2, centroid, imat, jmat, normal = elem.material_coordinate_system()
    return dxyz, centroid, imat, jmat, normal

def _get_quad_vectors_mcid(elem: CQUAD4):
    try:
        node1, node2, node3, node4 = elem.nodes_ref[:4]
    except (IndexError, ValueError):
        print(elem.get_stats())
        raise

    xyz1 = node1.get_position()
    xyz2 = node2.get_position()
    xyz3 = node3.get_position()
    xyz4 = node4.get_position()

    # take the mean length to size the vectors in the GUI
    dxyz21 = np.linalg.norm(xyz2 - xyz1)
    dxyz32 = np.linalg.norm(xyz3 - xyz2)
    dxyz43 = np.linalg.norm(xyz4 - xyz3)
    dxyz14 = np.linalg.norm(xyz1 - xyz4)
    dxyz = np.mean([dxyz21, dxyz32, dxyz43, dxyz14]) / 2.

    # adjusted for element coordinate system
    # equivalent to PCOMP thetad=0.0
    dxyz2, centroid, imat, jmat, normal = elem.material_coordinate_system()
    return dxyz, centroid, imat, jmat, normal

def _export_tria_mcid(model: BDF,
                      elem: ShellElement,
                      nodes,
                      iply: int, nid: int, eid: int,
                      #pids_failed: set[int],
                      bars,
                      export_both_axes: bool, export_xaxis: bool,
                      consider_property_rotation: bool) -> tuple[int, int]:
    """helper method for ``export_mcids``"""
    pid_ref, prop_type = get_pid_ref_prop_type(model, elem)
    if prop_type == 'PSHELL':
        mids = [mat.type for mat in pid_ref.materials()
                if mat is not None and mat.mid > 0]
        if 'MAT8' not in mids:
            return nid, eid
    elif prop_type in ['PCOMP', 'PCOMPG']:
        pass
    #elif prop_type in ['PLPLANE']:
        #return nid, eid
    else:
        raise NotImplementedError(pid_ref)

    # adjusted for element coordinate system
    # equivalent to PCOMP thetad=0.0
    dxyz, centroid, imat, jmat, normal = _get_tri_vectors_mcid(elem)

    nid, eid = _rotate_single_coord(
        elem, pid_ref, iply,
        nid, eid, nodes, bars,
        dxyz, centroid, imat, jmat, normal,
        export_both_axes=export_both_axes, export_xaxis=export_xaxis,
        consider_property_rotation=consider_property_rotation)
    return nid, eid

def _rotate_single_coord(elem: ShellElement,
                         pid_ref, iply: int,
                         nid, eid, nodes, bars,
                         dxyz: float,
                         centroid: np.ndarray,
                         imat: np.ndarray,
                         jmat: np.ndarray,
                         normal: np.ndarray,
                         export_both_axes: bool, export_xaxis: bool,
                         consider_property_rotation: bool) -> tuple[int, int]:
    # rotate the coord
    imat, jmat = _rotate_mcid(
        elem, pid_ref, iply, imat, jmat, normal,
        consider_property_rotation=consider_property_rotation)

    # size the axes based on the mean edge length
    iaxis = centroid + imat * dxyz
    jaxis = centroid + jmat * dxyz
    nid, eid = _add_elements(nid, eid, nodes, bars,
                             centroid, iaxis, jaxis,
                             export_both_axes, export_xaxis)
    return nid, eid

def _export_coord_axes(nodes, bars, csv_filename: str):
    """save the coordinate systems in a csv file"""
    if csv_filename:
        with open(csv_filename, 'w') as out_file:
            for node in nodes:
                out_file.write('GRID,%i,%s,%s,%s\n' % node)
            for bari in bars:
                out_file.write('BAR,%i,%i,%i\n' % bari)

def _rotate_mcid(elem: ShellElement,
                 pid_ref: Union[PCOMP, PCOMPG, PSHELL], iply: int,
                 imat: np.ndarray, jmat: np.ndarray, normal: np.ndarray,
                 consider_property_rotation: bool=True) -> tuple[np.ndarray, np.ndarray]:
    """
    Rotates a material coordinate system.  Assumes the element theta/mcid
    has already been acounted for.
    """
    if not consider_property_rotation:
        return imat, jmat

    if pid_ref.type == 'PSHELL':
        return imat, jmat
    elif pid_ref.type in ['PCOMP', 'PCOMPG']:
        thetad = pid_ref.get_theta(iply)
    else:
        raise NotImplementedError(f'property type={elem.pid_ref.type!r} is not supported\n'
                                  f'{elem}{elem.pid_ref}')

    assert isinstance(thetad, float), thetad
    if isinstance(thetad, float) and thetad == 0.0:
        return imat, jmat

    imat2, jmat2 = rotate_by_thetad(thetad, imat, jmat, normal)
    return imat2, jmat2

def _add_elements(nid: int, eid: int,
                  nodes: list[tuple[int, float, float, float]],
                  bars: list[tuple[int, int, int]],
                  centroid: np.ndarray,
                  iaxis: np.ndarray,
                  jaxis: np.ndarray,
                  export_both_axes: bool, export_xaxis: bool) -> tuple[int, int]:
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

def _export_xaxis(nid: int,
                  nodes: list[np.ndarray],
                  bars: list[np.ndarray],
                  centroid: np.ndarray,
                  iaxis: np.ndarray) -> int:
    nodes.append(centroid)
    nodes.append(iaxis)
    elemi = (nid, nid + 1)
    bars.append(elemi)  # x-axis
    nid += 2
    return nid
