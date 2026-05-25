from __future__ import annotations
from typing import Any, TYPE_CHECKING
import numpy as np

from cpylog import properties as log_properties
from pyNastran.bdf.cards.elements.shell import ShellElement
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, CTRIA3, CTRIA6, CTRIAR, CQUAD4, CQUAD8, CQUADR, CQUAD


def get_shell_material_coord(element: CTRIA3 | CTRIA6 | CTRIAR |
                                      CQUAD4 | CQUAD8 | CQUADR | CQUAD) -> tuple[int, float]:
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
            word += f'PSHELL: ilayer={ilayer + 1:d}'
    else:
        if ilayer == 0:
            word += 'PCOMP: Total'
        else:
            word += f'PCOMP: ilayer={ilayer:d}'
    return word


def check_for_missing_control_surface_boxes(name: str, cs_box_ids: list[int],
                                            box_id_to_caero_element_map: dict[int, int],
                                            log, store_msg: bool=False) -> tuple[list[int], str]:
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
            msg += f'boxes_to_show={boxes_to_show}\n'

        out_msg = store_error(log, store_msg, msg.rstrip())
    return boxes_to_show, out_msg


def store_error(log, store_msg: bool, msg: str) -> str:
    out_msg = ''
    if store_msg:
        n, filename = log_properties(nframe=2)
        out_msg = f'ERROR: {filename}:{n:d} {msg}\n'
    else:
        log.warning(msg)
    return out_msg


def store_warning(log, store_msg: bool, msg: str) -> str:
    out_msg = ''
    if store_msg:
        n, filename = log_properties(nframe=2)
        out_msg = f'WARNING: {filename}:{n:d} {msg}\n'
    else:
        log.warning(msg)
    return out_msg


def make_nid_map(nid_map: dict[int, int], nids: list[int]) -> dict[int, int]:
    """make the node map"""
    for i, nid in enumerate(nids):
        nid_map[nid] = i
    return nid_map


def get_elements_nelements_unvectorized(model: BDF) -> tuple[Any, int, list[dict[int, Any]]]:
    nelements = len(model.elements)
    #eid_map = self.gui.eid_map
    elements = model.elements
    superelements = None
    if model.superelement_models:
        superelements = []
        if nelements:
            superelements.append(np.zeros(nelements, dtype='int32'))
        for super_tuple, superelement in sorted(model.superelement_models.items()):
            if isinstance(super_tuple, int):
                super_id = super_tuple
            else:
                super_id = super_tuple[1]
            nelements2 = len(superelement.elements)
            if nelements2:
                superelements.append(np.ones(nelements2, dtype='int32') * super_id)
            nelements += nelements2
            elements.update(superelement.elements)
            #eid_map.update(superelement.eid_map)
        superelements = np.hstack(superelements)
        assert len(superelements) == nelements, f'len(superelements)={len(superelements):d} nelements={nelements:d}'
    return elements, nelements, superelements


def _vectorized_shell_normals(model: BDF, eid_map: dict[int, int],
                              nid_map: dict[int, int],
                              xyz_cid0: np.ndarray,
                              normals: np.ndarray) -> None:
    """Compute shell element normals vectorized using cross products on bulk arrays."""
    tri_types = {'CTRIA3', 'CTRIAR', 'CTRAX3', 'CPLSTN3', 'CPLSTS3',
                 'CTRIA6', 'CTRIAX', 'CTRIAX6', 'CTRAX6', 'CPLSTN6', 'CPLSTS6'}
    quad_types = {'CQUAD4', 'CQUADR', 'CSHEAR', 'CQUADX4', 'CPLSTN4', 'CPLSTS4',
                  'CQUAD8', 'CQUADX8', 'CPLSTN8', 'CPLSTS8', 'CQUAD', 'CQUADX'}

    tri_ieids = []
    tri_n1 = []
    tri_n2 = []
    tri_n3 = []
    quad_ieids = []
    quad_n1 = []
    quad_n2 = []
    quad_n3 = []
    quad_n4 = []

    for eid, elem in model.elements.items():
        etype = elem.type
        ieid = eid_map.get(eid, -1)
        if ieid == -1:
            continue
        if etype in tri_types:
            nids = elem.node_ids
            tri_ieids.append(ieid)
            tri_n1.append(nid_map[nids[0]])
            tri_n2.append(nid_map[nids[1]])
            tri_n3.append(nid_map[nids[2]])
        elif etype in quad_types:
            nids = elem.node_ids
            quad_ieids.append(ieid)
            quad_n1.append(nid_map[nids[0]])
            quad_n2.append(nid_map[nids[1]])
            quad_n3.append(nid_map[nids[2]])
            quad_n4.append(nid_map[nids[3]])

    if tri_ieids:
        tri_ieids = np.array(tri_ieids, dtype='int32')
        p1 = xyz_cid0[tri_n1]
        p2 = xyz_cid0[tri_n2]
        p3 = xyz_cid0[tri_n3]
        v = np.cross(p2 - p1, p3 - p1)
        norms = np.linalg.norm(v, axis=1, keepdims=True)
        norms[norms == 0.0] = 1.0
        normals[tri_ieids] = (v / norms).astype('float32')

    if quad_ieids:
        quad_ieids = np.array(quad_ieids, dtype='int32')
        p1 = xyz_cid0[quad_n1]
        p2 = xyz_cid0[quad_n2]
        p3 = xyz_cid0[quad_n3]
        p4 = xyz_cid0[quad_n4]
        # average of two triangle normals (diagonal cross)
        v = np.cross(p3 - p1, p4 - p2)
        norms = np.linalg.norm(v, axis=1, keepdims=True)
        norms[norms == 0.0] = 1.0
        normals[quad_ieids] = (v / norms).astype('float32')


def _vectorized_element_dims_offsets(model: BDF, eid_map: dict[int, int],
                                     normals: np.ndarray,
                                     element_dim: np.ndarray,
                                     nnodes_array: np.ndarray,
                                     offset: np.ndarray,
                                     xoffset: np.ndarray,
                                     yoffset: np.ndarray,
                                     zoffset: np.ndarray,
                                     has_normals: bool=False) -> set[int]:
    """Vectorized computation of element_dim, nnodes, and offsets.

    Returns the set of eids that were handled (so the fallback loop can skip them).
    """
    etype_dim_nnodes = {
        'CTETRA': (3, 4), 'CPENTA': (3, 6), 'CPYRAM': (3, 5), 'CHEXA': (3, 8),
        'CROD': (1, 2), 'CONROD': (1, 2), 'CBEND': (1, 2),
        'CBAR': (1, 2), 'CBEAM': (1, 2), 'CGAP': (1, 2), 'CTUBE': (1, 2),
        'CBUSH': (0, 2), 'CBUSH1D': (0, 2), 'CBUSH2D': (0, 2),
        'CFAST': (0, 2), 'CVISC': (0, 2),
        'CELAS1': (0, 2), 'CELAS2': (0, 2), 'CELAS3': (0, 2), 'CELAS4': (0, 2),
        'CDAMP1': (0, 2), 'CDAMP2': (0, 2), 'CDAMP3': (0, 2), 'CDAMP4': (0, 2), 'CDAMP5': (0, 2),
    }
    shell_nnodes = {
        'CTRIA3': 3, 'CTRIAR': 3, 'CTRAX3': 3,
        'CTRIA6': 6, 'CTRIAX': 6, 'CTRIAX6': 6, 'CTRAX6': 6,
        'CQUAD4': 4, 'CQUADR': 4, 'CSHEAR': 4, 'CQUADX4': 4,
        'CQUAD8': 8, 'CQUADX8': 8,
        'CQUAD': 9, 'CQUADX': 9,
        'CPLSTN3': 3, 'CPLSTN4': 4, 'CPLSTN6': 6, 'CPLSTN8': 8,
        'CPLSTS3': 3, 'CPLSTS4': 4, 'CPLSTS6': 6, 'CPLSTS8': 8,
    }
    skip_types = {'CHBDYP', 'CAABSF'}
    plane_strain_types = {'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8'}

    # collect non-shell elements for bulk assignment
    dim_groups: dict[tuple[int, int], list[int]] = {}
    # collect shell elements for offset computation
    shell_ieids = []
    shell_z0s = []
    handled_eids: set[int] = set()

    for eid, element in model.elements.items():
        etype = element.type
        if etype in skip_types:
            handled_eids.add(eid)
            continue

        ieid = eid_map.get(eid, -1)
        if ieid == -1:
            handled_eids.add(eid)
            continue

        if etype in etype_dim_nnodes:
            dim_val, nnodes_val = etype_dim_nnodes[etype]
            key = (dim_val, nnodes_val)
            if key not in dim_groups:
                dim_groups[key] = []
            dim_groups[key].append(ieid)
            handled_eids.add(eid)
        elif has_normals and etype in shell_nnodes:
            nnodesi = shell_nnodes[etype]
            element_dim[ieid] = 2
            nnodes_array[ieid] = nnodesi

            if etype in plane_strain_types:
                handled_eids.add(eid)
                continue

            # extract z0
            prop = element.pid_ref
            if prop is None:
                z0 = np.nan
            else:
                ptype = prop.type
                if ptype == 'PSHELL':
                    z0 = prop.z1
                elif ptype in {'PCOMP', 'PCOMPG'}:
                    z0 = prop.z0
                elif ptype == 'PLPLANE':
                    z0 = 0.
                elif ptype in {'PSHEAR', 'PSOLID', 'PLSOLID', 'PPLANE', 'PMIC'}:
                    z0 = np.nan
                else:
                    raise NotImplementedError(f'prop={prop}\n ptype={ptype!r}')

            if z0 is None:
                if etype in {'CTRIA3', 'CTRIAR'}:
                    z0 = (element.T1 + element.T2 + element.T3) / 3.
                elif etype == 'CTRIA6':
                    z0 = (element.T1 + element.T2 + element.T3) / 3.
                elif etype in {'CQUAD4', 'CQUADR'}:
                    z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                elif etype in {'CQUAD8', 'CQUAD'}:
                    z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                elif etype in {'CTRAX3', 'CTRAX6', 'CTRIAX', 'CTRIAX6',
                               'CQUADX', 'CQUADX4', 'CQUADX8', 'CSHEAR'}:
                    z0 = np.nan
                else:  # pragma: no cover
                    raise NotImplementedError(element)

            shell_ieids.append(ieid)
            shell_z0s.append(z0)
            handled_eids.add(eid)

    # bulk assign dim/nnodes for non-shell elements
    for (dim_val, nnodes_val), ieids in dim_groups.items():
        idx = np.array(ieids, dtype='int32')
        element_dim[idx] = dim_val
        nnodes_array[idx] = nnodes_val

    # bulk compute offsets for shell elements
    if shell_ieids:
        idx = np.array(shell_ieids, dtype='int32')
        z0_arr = np.array(shell_z0s, dtype='float32')
        offset[idx] = z0_arr
        xoffset[idx] = z0_arr * normals[idx, 0]
        yoffset[idx] = z0_arr * normals[idx, 1]
        zoffset[idx] = z0_arr * normals[idx, 2]

    return handled_eids


def build_offset_normals_dims(model: BDF, eid_map: dict[int, int],
                              nelements: int,
                              xyz_cid0: np.ndarray=None,
                              nid_map: dict[int, int]=None):
    normals = np.full((nelements, 3), np.nan, dtype='float32')
    offset = np.full(nelements, np.nan, dtype='float32')
    xoffset = np.full(nelements, np.nan, dtype='float32')
    yoffset = np.full(nelements, np.nan, dtype='float32')
    zoffset = np.full(nelements, np.nan, dtype='float32')
    element_dim = np.full(nelements, -1, dtype='int32')
    nnodes_array = np.full(nelements, -1, dtype='int32')
    log = model.log

    assert eid_map is not None

    # vectorized shell normal computation
    if xyz_cid0 is not None and nid_map is not None:
        _vectorized_shell_normals(model, eid_map, nid_map, xyz_cid0, normals)

    # vectorized element_dim, nnodes, and offset computation
    has_normals = xyz_cid0 is not None and nid_map is not None
    handled_eids = _vectorized_element_dims_offsets(
        model, eid_map, normals, element_dim, nnodes_array,
        offset, xoffset, yoffset, zoffset, has_normals=has_normals)

    # fallback loop for unhandled elements (CHBDYG, unknown types)
    etype_to_nnodes_map = {
        'CTRIA3': 3, 'CTRIAR': 3, 'CTRAX3': 3,
        'CTRIA6': 6, 'CTRIAX': 6, 'CTRIAX6': 6, 'CTRAX6': 6,
        'CQUAD4': 4, 'CQUADR': 4, 'CSHEAR': 4, 'CQUADX4': 4,
        'CQUAD8': 8, 'CQUADX8': 8,
        'CQUAD': 9, 'CQUADX': 9,
        'CPLSTN3': 3, 'CPLSTN4': 4, 'CPLSTN6': 6, 'CPLSTN8': 8,
        'CPLSTS3': 3, 'CPLSTS4': 4, 'CPLSTS6': 6, 'CPLSTS8': 8,
    }
    for eid, element in model.elements.items():
        if eid in handled_eids:
            continue
        etype = element.type
        if etype == 'CHBDYG':
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
                print(f'element.type={element.type} doesnt have a dimension')
        elif etype in {'CHBDYP', 'CAABSF', 'CBEND'}:
            continue
        elif isinstance(element, ShellElement):
            if etype not in etype_to_nnodes_map:
                log.warning(f'unexpected element ({etype}) in etype_to_nnodes_map\n{str(element)}')
            element_dimi = 2
            ieid = eid_map.get(eid, -1)
            if ieid == -1:
                log.warning(f'failed to find shell element {element.type}={eid}\n{str(element)}')
                continue
            if not np.isnan(normals[ieid, 0]):
                normali = normals[ieid]
            else:
                try:
                    normali = element.Normal()
                except AttributeError:
                    msg += (
                        f'{element}'
                        f'nodes_ref = {element}\n'
                        f'nodes = {element.node_ids}'
                    )
                    raise AttributeError(msg)
                except RuntimeError:
                    # degenerate tri
                    msg = (
                        f'eid={eid:d} normal=NaN...\n'
                        f'{element}'
                        f'nodes = {element.nodes}')
                    log.error(msg)
                    normali = np.ones(3) * np.nan

            prop = element.pid_ref
            if prop is None:
                # F:\work\pyNastran\examples\Dropbox\move_tpl\ehbus69.op2
                # CTRIAX
                z0 = np.nan
            else:
                ptype = prop.type
                if ptype == 'PSHELL':
                    z0 = prop.z1
                elif ptype in {'PCOMP', 'PCOMPG'}:
                    z0 = prop.z0
                elif ptype in {'PLPLANE'}:  #, 'PTRSHL', 'PQUAD1'}:  # ? PTRSHL, PQUAD1
                    z0 = 0.
                elif ptype in {'PSHEAR', 'PSOLID', 'PLSOLID', 'PPLANE', 'PMIC'}:
                    z0 = np.nan
                else:
                    # PCOMPG
                    raise NotImplementedError(f'prop={prop}\n ptype={ptype!r}')

            if z0 is None:
                if etype in {'CTRIA3', 'CTRIAR'}:
                    z0 = (element.T1 + element.T2 + element.T3) / 3.
                    nnodesi = 3
                elif etype == 'CTRIA6':
                    z0 = (element.T1 + element.T2 + element.T3) / 3.
                    nnodesi = 6
                elif etype in {'CQUAD4', 'CQUADR'}:
                    z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                    nnodesi = 4
                elif etype == 'CQUAD8':
                    z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                    nnodesi = 8
                elif etype == 'CQUAD':
                    z0 = (element.T1 + element.T2 + element.T3 + element.T4) / 4.
                    nnodesi = 9
                elif etype in {'CTRAX3', 'CTRAX6', 'CTRIAX', 'CTRIAX6',
                               'CQUADX', 'CQUADX4', 'CQUADX8'}:
                    z0 = np.nan
                    nnodesi = etype_to_nnodes_map[etype]
                else:  # pragma: no cover
                    raise NotImplementedError(element)
            else:
                nnodesi = etype_to_nnodes_map[etype]

            normals[ieid, :] = normali
            if element.type in {'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8'}:
                element_dim[ieid] = element_dimi
                nnodes_array[ieid] = nnodesi
                continue

            offset[ieid] = z0
            xoffset[ieid] = z0 * normali[0]
            yoffset[ieid] = z0 * normali[1]
            zoffset[ieid] = z0 * normali[2]
            element_dim[ieid] = element_dimi
            nnodes_array[ieid] = nnodesi
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
            log.warning(f'shell element.type={etype} doesnt have a dimension')
        element_dim[ieid] = element_dimi
        nnodes_array[ieid] = nnodesi
    return normals, offset, xoffset, yoffset, zoffset, element_dim, nnodes_array


def build_map_centroidal_result(model: BDF, nid_map: dict[int, int],
                                stop_on_failure: bool=True) -> None:
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


def _build_map_centroidal_result(model: BDF, nid_map: dict[int, int]) -> None:
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
        #'CTRSHL',
    }

    springs_dampers = {'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                       'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4'}
    nnodes_map = {
        'CTRIA6': (3, 6),
        'CQUAD8': (4, 8),
        'CQUAD': (9, 9),
        'CHEXA': (8, 20),
        'CTETRA': (4, 10),
        'CPENTA': (6, 15),
        'CPYRAM': (5, 13),
        'CQUADX': (4, 9),
        'CTRIAX': (3, 6),
        'CQUADX8': (4, 8),
        'CTRIAX6': (3, 6),
        # nastran 95
        # 'CQUAD1': (4, 4),
        # 'CTRSHL': (6, 6),
        # 'CHEXA1': (8, 8),
        # 'CHEXA20': (20, 20),
    }
    skip_cards = ['CAABSF', 'GENEL']
    #['CTRIA6', 'CQUAD8', 'CHEXA', 'CTETRA', 'CPENTA', 'CPYRAM', 'CQUADX', 'CTRIAX']
    etypes_mixed_nodes = set(list(nnodes_map.keys()))

    for eid, elem in sorted(model.elements.items()):
        elem_type = elem.type
        if elem_type in etypes_all_nodes:
            try:
                node_ids = [nid_map[nid] for nid in elem.nodes]
            except KeyError:  # pragma: no cover
                print(elem)
                raise
        elif elem_type in etypes_mixed_nodes:
            node_ids = [nid_map[nid] if nid is not None else 0
                        for nid in elem.nodes]
            nnodes_min, nnodes_max = nnodes_map[elem_type]
            nnodesi = len(node_ids)
            if nnodesi == nnodes_min:
                pass
            elif 0 in node_ids:
                node_ids = node_ids[:nnodes_min]
                assert len(node_ids) == nnodes_min, f'nnodes={len(node_ids):d} min={nnodes_min:d}\n{elem}'
            else:
                assert nnodesi == nnodes_max, f'nnodes={nnodesi:d} max={nnodes_max:d}\n{elem}'
        elif elem_type in springs_dampers:
            node_ids = [nid_map[nid] for nid in elem.nodes
                        if nid is not None]
        elif elem_type in {'CHBDYP', 'CHBDYG'}:
            node_ids = [nid_map[nid] for nid in elem.nodes
                        if nid not in [0, None]]
        elif elem_type in skip_cards:
            continue
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

    # build flat arrays for vectorized scatter-add
    # elem_indices[k] = which element index this entry came from
    # node_indices[k] = which node index to accumulate into
    elem_indices_list = []
    node_indices_list = []
    for i, node_idsi in enumerate(mapped_node_ids):
        for nid in node_idsi:
            elem_indices_list.append(i)
            node_indices_list.append(nid)
    elem_indices_arr = np.array(elem_indices_list, dtype='int32')
    node_indices_arr = np.array(node_indices_list, dtype='int32')

    def map_centroidal_result(centroidal_data):
        """maps centroidal data onto nodal data via np.add.at"""
        nodal_data = np.zeros(nnodes, dtype=centroidal_data.dtype)
        np.add.at(nodal_data, node_indices_arr, centroidal_data[elem_indices_arr])
        return nodal_data * inv_node_count
    model.map_centroidal_result = map_centroidal_result
