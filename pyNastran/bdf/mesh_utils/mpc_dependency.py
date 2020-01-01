"""
defines methods to access MPC/rigid element data:
  - get_mpc_node_ids( mpc_id, stop_on_failure=True)
  - get_mpc_node_ids_c1( mpc_id, stop_on_failure=True)
  - get_rigid_elements_with_node_ids(self, node_ids)
  - get_dependent_nid_to_components(self, mpc_id=None, stop_on_failure=True)
  - get_lines_rigid(model: BDF)
  - get_mpcs(model, mpc_id, mpc_id, consider_mpcadd=True,
             stop_on_failure=True)

"""
from __future__ import annotations
from collections import defaultdict
from typing import Tuple, List, Dict, Any, TYPE_CHECKING

import numpy as np
from pyNastran.utils.numpy_utils import integer_types

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

def get_mpc_node_ids(model: BDF, mpc_id: int,
                     consider_mpcadd: bool=True,
                     stop_on_failure: bool=True) -> List[List[int]]:
    r"""
    Get the MPC/MPCADD IDs.

    Parameters
    ----------
    mpc_id : int
        the MPC id
    consider_mpcadd : bool
        MPCADDs should not be considered when referenced from an MPCADD
        from a case control, True should be used.
    stop_on_failure : bool; default=True
        errors if parsing something new

    Returns
    -------
    lines : List[[independent, dependent]]
        independent : int
           the independent node id
        dependent : int
           the dependent node id

    I      I
      \   /
    I---D---I

    """
    lines = []
    mpcs = model.get_reduced_mpcs(
        mpc_id, consider_mpcadd=consider_mpcadd,
        stop_on_failure=stop_on_failure)

    # dependent, independent
    for card in mpcs:
        if card.type == 'MPC':
            nids = card.node_ids
            nid0 = nids[0]
            #component0 = card.components[0]
            #enforced0 = card.coefficients[0]
            #card.constraints[1:]
            for nid, coefficient in zip(nids[1:], card.coefficients[1:]):
                if coefficient != 0.0:
                    lines.append([nid0, nid])
        else:
            msg = 'get_MPCx_node_ids doesnt support %r' % card.type
            if stop_on_failure:
                raise RuntimeError(msg)
            model.log.warning(msg)
    return lines

def get_mpc_node_ids_c1(model: BDF, mpc_id: int,
                        consider_mpcadd: bool=True,
                        stop_on_failure: bool=True) -> Tuple[Dict[str, List[int]],
                                                             Dict[str, List[int]]]:
    r"""
    Get the MPC/MPCADD IDs.

    Parameters
    ----------
    mpc_id : int
        the MPC id
    consider_mpcadd : bool
        MPCADDs should not be considered when referenced from an MPCADD
        from a case control, True should be used.
    stop_on_failure : bool; default=True
        errors if parsing something new

    Returns
    -------
    independent_node_ids_c1 : Dict[component] = node_ids
        component : str
            the DOF to constrain
        node_ids : List[int]
            the constrained node ids
    dependent_node_ids_c1 : Dict[component] = node_ids
        component : str
            the DOF to constrain
        node_ids : List[int]
            the constrained node ids

    I      I
      \   /
    I---D---I

    """
    if not isinstance(mpc_id, integer_types):
        msg = 'mpc_id must be an integer; type=%s, mpc_id=\n%r' % (type(mpc_id), mpc_id)
        raise TypeError(msg)

    mpcs = model.get_reduced_mpcs(
        mpc_id, consider_mpcadd=consider_mpcadd,
        stop_on_failure=stop_on_failure)

    # dependent, independent
    independent_node_ids_c1 = defaultdict(list)
    dependent_node_ids_c1 = defaultdict(list)
    for card in mpcs:
        if card.type == 'MPC':
            nids = card.node_ids
            nid0 = nids[0]
            #component0 = card.components[0]
            #coefficient0 = card.coefficients[0]
            #card.constraints[1:]
            dofs = card.components
            for dof in dofs:
                independent_node_ids_c1[dof].append(nid0)
            for nid, coefficient in zip(nids[1:], card.coefficients[1:]):
                if coefficient != 0.0:
                    for dof in dofs:
                        dependent_node_ids_c1[dof].append(nid)
        else:
            msg = 'get_MPCx_node_ids_c1 doesnt support %r' % card.type
            if stop_on_failure:
                raise RuntimeError(msg)
            model.log.warning(msg)
    return dict(independent_node_ids_c1), dict(dependent_node_ids_c1)

def get_rigid_elements_with_node_ids(model: BDF, node_ids):
    """
    Gets the series of rigid elements that use specific nodes

    Parameters
    ----------
    node_ids : List[int]
        the node ids to check

    Returns
    -------
    rbes : List[int]
        the set of self.rigid_elements

    """
    try:
        nids = set(node_ids)
    except TypeError:
        print(node_ids)
        raise
    rbes = []
    for eid, rigid_element in model.rigid_elements.items():
        if rigid_element.type in ['RBE3', 'RBE2', 'RBE1', 'RBAR', 'RSPLINE', 'RROD', 'RBAR1']:
            independent_nodes = set(rigid_element.independent_nodes)
            dependent_nodes = set(rigid_element.dependent_nodes)
            rbe_nids = independent_nodes | dependent_nodes
            if nids.intersection(rbe_nids):
                rbes.append(eid)
        elif rigid_element.type == 'RSSCON':
            msg = 'skipping card in get_rigid_elements_with_node_ids\n%s' % str(rigid_element)
            model.log.warning(msg)
        else:
            raise RuntimeError(rigid_element.type)
    return rbes

def get_dependent_nid_to_components(model: BDF, mpc_id=None, stop_on_failure=True):
    """
    Gets a dictionary of the dependent node/components.

    Parameters
    ----------
    mpc_id : int; default=None -> no MPCs are checked
        TODO: add
    stop_on_failure : bool; default=True
        errors if parsing something new

    Returns
    -------
    dependent_nid_to_components : dict[node_id] : components
        node_id : int
            the node_id
        components : str
            the DOFs that are linked

    Nastran can either define a load/motion at a given node.
    SPCs define constraints that may not have loads/motions.

    MPCs and rigid elements define independent and dependent nodes on
    specific DOFs.
      - independent nodes : loads/motions may be defined
      - dependent nodes : loads/motions may not be defined

    """
    dependent_nid_to_components = {}

    if mpc_id is not None:
        mpcs = get_mpcs(model, mpc_id)
        for mpc in mpcs:
            if mpc.type == 'MPC':
                for nid, component in zip(mpc.node_ids, mpc.components):
                    dependent_nid_to_components[nid] = component
            else:
                raise NotImplementedError(mpc)

    for unused_eid, rigid_element in model.rigid_elements.items():
        if rigid_element.type == 'RBE2':
            dependent_nodes = set(rigid_element.dependent_nodes)
            components = rigid_element.cm
            for nid in dependent_nodes:
                dependent_nid_to_components[nid] = components
        elif rigid_element.type == 'RBE3':
            dependent_nid_to_components[rigid_element.ref_grid_id] = rigid_element.refc
            for gmi, cmi in zip(rigid_element.Gmi_node_ids, rigid_element.Cmi):
                dependent_nid_to_components[gmi] = cmi
        #if rigid_element.type in ['RBE3', 'RBE2', 'RBE1', 'RBAR']:
            ##independent_nodes = set(rigid_element.independent_nodes)
            #dependent_nodes = set(rigid_element.dependent_nodes)
            #rbe_nids = independent_nodes | dependent_nodes
            #if nids.intersection(rbe_nids):
                #rbes.append(eid)
        #elif rigid_element == 'RSPLINE':
        elif rigid_element.type == 'RBAR':
            nodes = [rigid_element.ga, rigid_element.gb]
            components = [rigid_element.cma, rigid_element.cmb]
            for nid, componentsi in zip(nodes, components):
                dependent_nid_to_components[nid] = componentsi
        elif rigid_element.type == 'RBAR1':
            for componentsi in rigid_element.cb:
                dependent_nid_to_components[rigid_element.gb] = componentsi
        elif rigid_element.type == 'RBE1':
            # +------+-----+-----+-----+-------+-----+-----+-----+
            # |   1  |  2  |  3  |  4  |   5   |  6  |  7  |  8  |
            # +======+=====+=====+=====+=======+=====+=====+=====+
            # | RBE1 | EID | GN1 | CN1 |  GN2  | CN2 | GN3 | CN3 |
            # |      |     | GN4 | CN4 |  GN5  | CN5 | GN6 | CN6 |
            # |      | UM  | GM1 | CM1 |  GM2  | CM2 | GM3 | CM3 |
            # |      | GM4 | CM4 | etc | ALPHA |     |     |     |
            # +------+-----+-----+-----+-------+-----+-----+-----+
            # | RBE1 | 59  | 59  | 123 |  60   | 456 |     |     |
            # |      | UM  | 61  | 246 |       |     |     |     |
            # +------+-----+-----+-----+-------+-----+-----+-----+
            # dependent=m (independent=n)
            for nid, componentsi in zip(rigid_element.Gmi_node_ids, rigid_element.Cmi):
                dependent_nid_to_components[nid] = componentsi
            #dependent = elem.dependent_nodes
            #independent = elem.independent_nodes
            #assert len(dependent) == 1, dependent
            #assert len(independent) == 1, independent
            #lines_rigid.append([dependent[0], independent[0]])
        elif rigid_element.type == 'RROD':
            components = [rigid_element.cma, rigid_element.cmb]
            if rigid_element.cma is not None:
                nid = rigid_element.nodes[0]
                for component in rigid_element.cma:
                    dependent_nid_to_components[nid] = component

            if rigid_element.cmb is not None:
                nid = rigid_element.nodes[1]
                for component in rigid_element.cmb:
                    dependent_nid_to_components[nid] = component
        elif rigid_element.type == 'RSPLINE':
            #independent_nid = rigid_element.independent_nid
            for nid, component in zip(rigid_element.dependent_nids, rigid_element.dependent_components):
                if component is None:
                    continue
                dependent_nid_to_components[nid] = component
        elif rigid_element.type == 'RSSCON':
            msg = 'skipping card in get_dependent_nid_to_components\n%s' % str(rigid_element)
            model.log.warning(msg)
        else:
            raise RuntimeError(rigid_element.type)
    return dependent_nid_to_components

def get_lines_rigid(model: BDF) -> Any:
    """
    GUI helper function

    dependent = (lines[:, 0])
    independent = np.unique(lines[:, 1])

    """
    lines_rigid = []
    for eid, elem in model.rigid_elements.items():
        if elem.type == 'RBE3':
            if elem.Gmi != []:
                # UM are dependent
                msg = 'UM is not supported; RBE3 eid=%s Gmi=%s' % (elem.eid, elem.Gmi)
                raise RuntimeError(msg)
            #list_fields = ['RBE3', elem.eid, None, elem.ref_grid_id, elem.refc]
            n1 = elem.ref_grid_id
            assert isinstance(n1, integer_types), 'RBE3 eid=%s ref_grid_id=%s' % (elem.eid, n1)
            for (_weight, ci, Gij) in zip(elem.weights, elem.comps, elem.Gijs):
                Giji = elem._node_ids(nodes=Gij, allow_empty_nodes=True)
                # list_fields += [wt, ci] + Giji
                for n2 in Giji:
                    assert isinstance(n2, integer_types), 'RBE3 eid=%s Giji=%s' % (elem.eid, Giji)
                    lines_rigid.append([n1, n2])
        elif elem.type == 'RBE2':
            #list_fields = ['RBE2', elem.eid, elem.Gn(), elem.cm
                           #] + elem.Gmi_node_ids + [elem.alpha]
            n2 = elem.Gn() # independent
            nids1 = elem.Gmi_node_ids # dependent
            for n1 in nids1:
                lines_rigid.append([n1, n2])
        elif elem.type in ['RBAR', 'RBAR1', 'RROD']: ## TODO: these aren't quite right
            dependent = elem.Ga()
            independent = elem.Gb()
            lines_rigid.append([dependent, independent])
        elif elem.type == 'RBE1':
            # +------+-----+-----+-----+-------+-----+-----+-----+
            # |   1  |  2  |  3  |  4  |   5   |  6  |  7  |  8  |
            # +======+=====+=====+=====+=======+=====+=====+=====+
            # | RBE1 | EID | GN1 | CN1 |  GN2  | CN2 | GN3 | CN3 |
            # |      |     | GN4 | CN4 | GN5   | CN5 | GN6 | CN6 |
            # |      | UM  | GM1 | CM1 |  GM2  | CM2 | GM3 | CM3 |
            # |      | GM4 | CM4 | etc | ALPHA |     |     |     |
            # +------+-----+-----+-----+-------+-----+-----+-----+
            # | RBE1 | 59  | 59  | 123 |  60   | 456 |     |     |
            # |      | UM  | 61  | 246 |       |     |     |     |
            # +------+-----+-----+-----+-------+-----+-----+-----+
            dependent = elem.dependent_nodes
            independent = elem.independent_nodes
            #assert len(dependent) == 1, dependent
            #assert len(independent) == 1, independent
            if len(independent) != 1 or len(dependent) != 1:
                msg = 'skipping card because len(independent) != 1 or len(dependent) != 1\n'
                msg += '  independent = %s\n'  % independent
                msg += '  dependent = %s\n'  % dependent
                msg += str(elem)
                model.log.error(msg)
                continue
            lines_rigid.append([dependent[0], independent[0]])
        elif elem.type == 'RSPLINE':
            independent_nid = elem.independent_nid
            for dependent_nid in np.unique(elem.dependent_nids):
                lines_rigid.append([dependent_nid, independent_nid])
        elif elem.type == 'RSSCON':
            model.log.warning('skipping card in _get_rigid\n%s' % str(elem))
        else:
            print(str(elem))
            raise NotImplementedError(elem.type)
    return lines_rigid

def get_mpcs(model, mpc_id: int, consider_mpcadd: bool=True,
             stop_on_failure: bool=True) -> Tuple[List[int], List[str]]:
    """
    Gets the MPCs in a semi-usable form.

    Parameters
    ----------
    mpc_id : int
        the desired MPC ID
    stop_on_failure : bool; default=True
        errors if parsing something new

    Returns
    -------
    nids : List[int]
        the constrained nodes
    comps : List[str]
        the components that are constrained on each node

    Considers:
      - MPC
      - MPCADD

    """
    mpcs = model.get_reduced_mpcs(
        mpc_id, consider_mpcadd=consider_mpcadd,
        stop_on_failure=stop_on_failure)
    nids = []
    comps = []
    for mpc in mpcs:
        if mpc.type == 'MPC':
            for nid, comp, unused_coefficient in zip(mpc.nodes, mpc.components, mpc.coefficients):
                nids.append(nid)
                comps.append(comp)
        else:
            if stop_on_failure:
                model.log.error('not considering:\n%s' % str(mpc))
                raise NotImplementedError(mpc)
            model.log.warning('not considering:\n%s' % str(mpc))
    return nids, comps
