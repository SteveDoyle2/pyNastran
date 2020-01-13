from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from cpylog import properties as log_properties
from pyNastran.bdf.cards.elements.shell import ShellElement
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

def get_shell_material_coord(element) -> Tuple[int, float]:
    """
    used by:
     - CQUAD4, CQUADR, CQUAD8, CQUAD, CQUADX
     - CTRIA3, CTRIAR, CTRIA6
    """
    if isinstance(element.theta_mcid, float):
        return -1, element.theta_mcid
    else:
        return element.theta_mcid, np.nan

def get_nastran_gui_layer_word(i: int, ilayer: int, is_pshell_pcomp: bool) -> str:
    """gets the PSHELL/PCOMP layer word"""
    ## TODO: this makes no sense...
    is_pshell, unused_is_pcomp = is_pshell_pcomp
    word = ''
    if i == 0:
        if ilayer == 0:
            if is_pshell:
                word += 'PSHELL: ilayer=1 & others'
            else:
                word += 'Other Properties'
        else:
            word += 'PSHELL: ilayer=%i' % (ilayer + 1)
    else:
        if ilayer == 0:
            word += 'PCOMP: Total'
        else:
            word += 'PCOMP: ilayer=%i' % (ilayer)
    return word

def check_for_missing_control_surface_boxes(name: str, cs_box_ids: List[int],
                                            box_id_to_caero_element_map: Dict[int, int],
                                            log, store_msg: bool=False) -> Tuple[List[int], str]:
    """helper method for creating control surface"""
    boxes_to_show = []
    missing_boxes = []
    for box_id in cs_box_ids:
        try:
            unused_ipoints = box_id_to_caero_element_map[box_id]
        except KeyError:
            missing_boxes.append(box_id)
            continue
        boxes_to_show.append(box_id)

    out_msg = ''
    if missing_boxes:
        #print('\nboxes_to_show =', boxes_to_show)
        msg = 'Missing CAERO AELIST/SPLINE control surface %r boxes: %s\n' % (
            name, str(missing_boxes))
        if not boxes_to_show:
            msg += 'boxes_to_show=%s\n' % boxes_to_show

        out_msg = store_error(log, store_msg, msg.rstrip())
    return boxes_to_show, out_msg

def store_error(log, store_msg: bool, msg: str) -> str:
    out_msg = ''
    if store_msg:
        n, filename = log_properties(nframe=2)
        out_msg = 'ERROR: %s:%s %s\n' % (filename, n, msg)
    else:
        log.warning(msg)
    return out_msg

def store_warning(log, store_msg: bool, msg: str) -> str:
    out_msg = ''
    if store_msg:
        n, filename = log_properties(nframe=2)
        out_msg = 'WARNING: %s:%s %s\n' % (filename, n, msg)
    else:
        log.warning(msg)
    return out_msg


def make_nid_map(nid_map: Dict[int, int], nids: List[int]) -> Dict[int, int]:
    """make the node map"""
    for i, nid in enumerate(nids):
        nid_map[nid] = i
    return nid_map

def get_elements_nelements_unvectorized(model: BDF) -> Tuple[Any, int, List[Dict[int, Any]]]:
    nelements = len(model.elements)
    #eid_map = self.gui.eid_map
    elements = model.elements
    superelements = None
    if model.superelement_models:
        superelements = []
        if nelements:
            superelements.append(np.zeros(nelements, dtype='int32'))
        for super_id, superelement in sorted(model.superelement_models.items()):
            nelements2 = len(superelement.elements)
            if nelements2:
                superelements.append(np.ones(nelements2, dtype='int32') * super_id)
            nelements += nelements2
            elements.update(superelement.elements)
            #eid_map.update(superelement.eid_map)
        superelements = np.hstack(superelements)
        assert len(superelements) == nelements, 'len(superelements)=%s nelements=%s' % (len(superelements), nelements)
    return elements, nelements, superelements


def build_offset_normals_dims(model: BDF, eid_map: Dict[int, int], nelements: int):
    normals = np.full((nelements, 3), np.nan, dtype='float32')
    offset = np.full(nelements, np.nan, dtype='float32')
    xoffset = np.full(nelements, np.nan, dtype='float32')
    yoffset = np.full(nelements, np.nan, dtype='float32')
    zoffset = np.full(nelements, np.nan, dtype='float32')
    element_dim = np.full(nelements, -1, dtype='int32')
    nnodes_array = np.full(nelements, -1, dtype='int32')
    log = model.log

    #eid_map = self.gui.eid_map
    assert eid_map is not None
    etype_to_nnodes_map = {
        'CTRIA3' : 3, 'CTRIAR' : 3, 'CTRAX3' : 3, 'CPLSTN3' : 3,
        # no a CTRIAX really has 6 nodes because reasons...
        'CTRIA6' : 6, 'CTRIAX' : 6, 'CTRIAX6' : 6, 'CPLSTN6' : 6, 'CTRAX6' : 6,
        'CQUAD4' : 4, 'CQUADR' : 4, 'CPLSTN4' : 4, 'CSHEAR' : 4, 'CQUADX4' : 4,
        'CQUAD8' : 8, 'CPLSTN8' : 8, 'CQUADX8' : 8,
        'CQUAD' : 9, 'CQUADX' : 9,
        'CPLSTN3': 3, 'CPLSTN4': 4, 'CPLSTN6': 6, 'CPLSTN8': 8,
        'CPLSTS3': 3, 'CPLSTS4': 4, 'CPLSTS6': 6, 'CPLSTS8': 8,
    }
    for eid, element in sorted(model.elements.items()):
        etype = element.type
        if isinstance(element, ShellElement):
            ieid = None
            element_dimi = 2
            #assert element.nodes_ref is not None, element.nodes_ref
            try:
                normali = element.Normal()
            except AttributeError:
                msg = str(element)
                msg += 'nodes_ref = %s\n' % element
                msg += 'nodes = %s' % element.node_ids
                raise AttributeError(msg)
            except RuntimeError:
                # this happens when you have a degenerate tri
                msg = (
                    'eid=%i normal=NaN...\n'
                    '%s'
                    'nodes = %s' % (eid, element, str(element.nodes)))
                log.error(msg)
                normali = np.ones(3) * np.nan
                #raise

            prop = element.pid_ref
            if prop is None:
                # F:\work\pyNastran\examples\Dropbox\move_tpl\ehbus69.op2
                # CTRIAX
                z0 = np.nan
            else:
                ptype = prop.type
                if ptype == 'PSHELL':
                    z0 = prop.z1
                elif ptype in ['PCOMP', 'PCOMPG']:
                    z0 = prop.z0
                elif ptype == 'PLPLANE':
                    z0 = 0.
                elif ptype in ['PSHEAR', 'PSOLID', 'PLSOLID', 'PPLANE']:
                    z0 = np.nan
                else:
                    raise NotImplementedError(ptype) # PSHEAR, PCOMPG

            if z0 is None:
                if etype in ['CTRIA3', 'CTRIAR']:
                    #node_ids = self.nodes[3:]
                    z0 = (element.T1 + element.T2 + element.T3) / 3.
                    nnodesi = 3
                elif etype == 'CTRIA6':
                    #node_ids = self.nodes[3:]
                    z0 = (element.T1 + element.T2 + element.T3) / 3.
                    nnodesi = 6
                elif etype in ['CQUAD4', 'CQUADR']:
                    #node_ids = self.nodes[4:]
                    z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                    nnodesi = 4
                elif etype == 'CQUAD8':
                    #node_ids = self.nodes[4:]
                    z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                    nnodesi = 8
                elif etype == 'CQUAD':
                    #node_ids = self.nodes[4:]
                    z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                    nnodesi = 9

                # axisymmetric
                elif etype == 'CTRAX3':
                    #node_ids = self.nodes[3:]
                    nnodesi = 3
                    z0 = np.nan
                elif etype == 'CTRAX6':
                    #node_ids = self.nodes[3:]
                    nnodesi = 6
                    z0 = np.nan
                elif etype in ['CTRIAX', 'CTRIAX6']:
                    # the CTRIAX6 uses a non-standard node orientation
                    #node_ids = self.nodes[3:]
                    z0 = np.nan
                    nnodesi = 6
                elif etype == 'CQUADX':
                    #node_ids = self.nodes[4:]
                    nnodesi = 9
                    z0 = np.nan
                elif etype == 'CQUADX4':
                    #node_ids = self.nodes[4:]
                    nnodesi = 4
                    z0 = np.nan
                elif etype == 'CQUADX8':
                    #node_ids = self.nodes[4:]
                    nnodesi = 8
                    z0 = np.nan
                else:
                    raise NotImplementedError(element)
            else:
                nnodesi = etype_to_nnodes_map[etype]

            ieid = eid_map[eid]
            normals[ieid, :] = normali
            if element.type in ['CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8']:
                element_dim[ieid] = element_dimi
                nnodes_array[ieid] = nnodesi
                log.debug('continue...element.type=%r' % element.type)
                continue

            offset[ieid] = z0
            xoffset[ieid] = z0 * normali[0]
            yoffset[ieid] = z0 * normali[1]
            zoffset[ieid] = z0 * normali[2]

        elif etype == 'CTETRA':
            ieid = eid_map[eid]
            element_dimi = 3
            nnodesi = 4
        elif etype == 'CPENTA':
            ieid = eid_map[eid]
            element_dimi = 3
            nnodesi = 6
        elif etype == 'CPYRAM':
            ieid = eid_map[eid]
            element_dimi = 3
            nnodesi = 5
        elif etype in ['CHEXA', 'CIHEX1', 'CIHEX2']:
            ieid = eid_map[eid]
            element_dimi = 3
            nnodesi = 8

        elif etype in ['CROD', 'CONROD', 'CBEND', 'CBAR', 'CBEAM', 'CGAP', 'CTUBE']:
            ieid = eid_map[eid]
            element_dimi = 1
            nnodesi = 2
        elif etype in ['CBUSH', 'CBUSH1D', 'CBUSH2D',
                       'CFAST', 'CVISC',
                       'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                       'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5']:
            ieid = eid_map[eid]
            element_dimi = 0
            nnodesi = 2
        elif etype == 'CHBDYG':
            ieid = eid_map[eid]
            surface_type = element.surface_type
            if surface_type == 'AREA3':
                nnodesi = 3
                element_dimi = 2
            elif surface_type == 'AREA4':
                nnodesi = 4
                element_dimi = 2
            elif surface_type == 'AREA6':
                nnodesi = 6
                element_dimi = 2
            elif surface_type == 'AREA8':
                nnodesi = 8
                element_dimi = 2
            #elif surface_type == 'REV':
                #nnodesi = 2 # ???
                #element_dimi = 1 # ???
            else:
                element_dimi = -1
                nnodesi = -1
                print('element.type=%s doesnt have a dimension' % element.type)
        elif etype == 'CHBDYP':
            continue
        else:
            try:
                ieid = eid_map[eid]
            except KeyError:
                print('didnt add element to eid_map')
                print(model.elements[eid])
                raise
            element_dimi = -1
            nnodesi = -1
            print('element.type=%s doesnt have a dimension' % element.type)
        assert ieid is not None
        element_dim[ieid] = element_dimi
        nnodes_array[ieid] = nnodesi
        #ielement += 1
    return normals, offset, xoffset, yoffset, zoffset, element_dim, nnodes_array

def build_map_centroidal_result(model: BDF, nid_map: Dict[int, int], stop_on_failure: bool=True) -> None:
    """
    Sets up map_centroidal_result.  Used for:
     - cutting plane
     - nodal stress/strain
    """
    try:
        # test_gui_superelement_1
        _build_map_centroidal_result(model, nid_map)
    except Exception as error:
        if len(model.superelement_models):
            model.log.error('superelements not supported in build_map_centroidal_result')
            model.log.error(str(error))
            return
        model.log.error('Cannot run build_map_centroidal_result')
        model.log.error(str(error))
        if stop_on_failure:
            raise
        raise

def _build_map_centroidal_result(model: BDF, nid_map: Dict[int, int]) -> None:
    """
    Sets up map_centroidal_result.  Used for:
     - cutting plane
     - nodal stress/strain
    """
    if hasattr(model, 'map_centroidal_result'):
        return
    mapped_node_ids = []
    nnodes = model.npoints
    node_count = np.zeros(nnodes, dtype='float32')

    eids = []
    etypes_all_nodes = {
        'CBEAM', 'CBAR', 'CROD', 'CONROD', 'CTUBE', 'CBEND',
        'CTRIA3', 'CQUAD4', 'CQUADR', 'CTRIAR', 'CSHEAR',
        'CGAP', 'CFAST', 'CVISC', 'CBUSH', 'CBUSH1D', 'CBUSH2D',
        'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8',
        'CPLSTS3', 'CPLSTS4', 'CPLSTS6', 'CPLSTS8',
    }

    springs_dampers = {'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                       'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4'}
    nnodes_map = {
        'CTRIA6' : (3, 6),
        'CQUAD8' : (4, 8),
        'CQUAD' : (9, 9),
        'CHEXA' : (8, 20),
        'CTETRA' : (4, 10),
        'CPENTA' : (6, 15),
        'CPYRAM' : (5, 13),
        'CQUADX' : (4, 9),
        'CTRIAX' : (3, 6),
        'CQUADX8' : (4, 8),
        'CTRIAX6' : (3, 6),
    }
    #['CTRIA6', 'CQUAD8', 'CHEXA', 'CTETRA', 'CPENTA', 'CPYRAM', 'CQUADX', 'CTRIAX']
    etypes_mixed_nodes = set(list(nnodes_map.keys()))

    for eid, elem in sorted(model.elements.items()):
        if elem.type in etypes_all_nodes:
            try:
                node_ids = [nid_map[nid] for nid in elem.nodes]
            except KeyError:  # pragma: no cover
                print(elem)
                raise
        elif elem.type in etypes_mixed_nodes:
            node_ids = [nid_map[nid] if nid is not None else 0
                        for nid in elem.nodes]
            nnodes_min, nnodes_max = nnodes_map[elem.type]
            nnodesi = len(node_ids)
            if nnodesi == nnodes_min:
                pass
            elif 0 in node_ids:
                node_ids = node_ids[:nnodes_min]
                assert len(node_ids) == nnodes_min, 'nnodes=%s min=%s\n%s' % (len(node_ids), nnodes_min, elem)
            else:
                assert nnodesi == nnodes_max, 'nnodes=%s max=%s\n%s' % (nnodesi, nnodes_max, elem)
        elif elem.type in springs_dampers:
            node_ids = [nid_map[nid] for nid in elem.nodes
                        if nid is not None]
        elif elem.type in {'CHBDYP', 'CHBDYG'}:
            node_ids = [nid_map[nid] for nid in elem.nodes
                        if nid not in [0, None]]
        else:
            raise NotImplementedError(elem)
        #print(elem.nodes, node_ids)
        mapped_node_ids.append(node_ids)
        eids.append(eid)

        # these are indicies
        for node_id in node_ids:
            node_count[node_id] += 1

    # calculate inv_node_count
    izero = np.where(node_count == 0)
    node_count[izero] = 1.
    inv_node_count = 1. / node_count

    # build the centroidal mapper
    def map_centroidal_result(centroidal_data):
        """maps centroidal data onto nodal data"""
        nodal_data = np.zeros(nnodes, dtype=centroidal_data.dtype)
        for unused_eid, datai, node_ids in zip(eids, centroidal_data, mapped_node_ids):
            for nid in node_ids:
                nodal_data[nid] += datai
        return nodal_data * inv_node_count
    model.map_centroidal_result = map_centroidal_result
