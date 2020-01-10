# coding: utf-8
"""
This file defines:
  - model, nid_offset, eid_offset = bdf_mirror(bdf_filename, plane='xz')
  - model, nid_offset, eid_offset = write_bdf_symmetric(
        bdf_filename, out_filename=None, encoding=None,
        size=8, is_double=False,
        enddata=None, close=True, plane='xz')

"""
from typing import Set, Tuple, Union, Optional
import numpy as np

from pyNastran.bdf.cards.coordinate_systems import CORD1R, CORD1C, CORD1S, CORD2R, CORD2C, CORD2S
from pyNastran.bdf.cards.loads.static_loads import (
    FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
    PLOAD, PLOAD2, PLOAD4)
from pyNastran.bdf.cards.thermal.loads import TEMP, QVOL, QBDY1, QBDY2, QBDY3, QHBDY
from pyNastran.bdf.cards.aero.aero import CAERO1, SPLINE1
from pyNastran.bdf.cards.bdf_sets import SET1 #, SET3
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils.internal_utils import get_bdf_model


def bdf_mirror_plane(bdf_filename: Union[str, BDF], plane, mirror_model=None,
                     log=None, debug: bool=True, use_nid_offset: bool=True):
    """mirrors a model about an arbitrary plane"""
    model = get_bdf_model(bdf_filename, xref=True, log=log, debug=debug)
    if mirror_model is None:
        mirror_model = BDF(debug=debug, log=log, mode='msc')

    nid_offset, plane = _mirror_nodes_plane(model, mirror_model, plane,
                                            use_nid_offset=use_nid_offset)
    eid_offset = _mirror_elements(model, mirror_model, nid_offset, use_eid_offset=True)
    #_mirror_loads(model, nid_offset, eid_offset)
    return model, mirror_model, nid_offset, eid_offset


def bdf_mirror(bdf_filename: Union[str, BDF],
               plane: str='xz', log=None, debug: bool=True):
    """
    Mirrors the model about the symmetry plane

    Parameters
    ----------
    bdf_filename : str / BDF()
        str : the bdf filename
        BDF : the BDF model object
    plane : str; {'xy', 'yz', 'xz'}; default='xz'
        the plane to mirror about
        xz : +y/-y
        yz : +x/-x
        xy : +z/-z

    Returns
    -------
    model : BDF()
        BDF : the BDF model object
    nid_offset : int
        the offset node id
    eid_offset : int
        the offset element id

    """
    model = get_bdf_model(bdf_filename, xref=True, log=log, debug=debug)
    mirror_model = model
    nid_offset, plane = _mirror_nodes(model, plane=plane)
    eid_offset = _mirror_elements(model, mirror_model, nid_offset, plane=plane)
    _mirror_loads(model, nid_offset, eid_offset)
    _mirror_aero(model, nid_offset, plane=plane)
    return model, nid_offset, eid_offset

def write_bdf_symmetric(bdf_filename: Union[str, BDF],
                        out_filename=None, encoding=None,
                        size: int=8, is_double: bool=False,
                        enddata: Optional[bool]=None, close: bool=True,
                        plane: str='xz', log=None):
    """
    Mirrors the model about the symmetry plane

    Parameters
    ----------
    bdf_filename : str / BDF()
        str : the bdf filename
        BDF : the BDF model object
    out_filename : varies; default=None
        str        - the name to call the output bdf
        file       - a file object
        StringIO() - a StringIO object
        None       - pops a dialog
    encoding : str; default=None -> system specified encoding
        the unicode encoding
        latin1, and utf8 are generally good options
    size : int; {8, 16}
        the field size
    is_double : bool; default=False
        False : small field
        True : large field
    enddata : bool; default=None
        bool - enable/disable writing ENDDATA
        None - depends on input BDF
    close : bool; default=True
        should the output file be closed
    plane : str; {'xy', 'yz', 'xz'}; default='xz'
        the plane to mirror about
        xz : +y/-y
        yz : +x/-x
        xy : +z/-z

    Returns
    -------
    model : BDF()
        BDF : the BDF model object
    nid_offset : int
        the offset node id
    eid_offset : int
        the offset element id

    Notes
    -----
    Updates the BDF object to be symmetric
      - see bdf_mirror if you don't want to write the model
    Doesn't equivalence nodes on the centerline.

    Considers
     - nodes : GRID
     - elements, rigid_elements, mass_elements : see ``_mirror_elements``
     - loads : see ``_mirror_loads``
     - aero cards : see ``_mirror_aero``

    """
    model, nid_offset, eid_offset = bdf_mirror(bdf_filename, plane=plane, log=None)
    model.write_bdf(out_filename=out_filename, encoding=encoding,
                    size=size, is_double=is_double,
                    interspersed=False, enddata=enddata, close=close)
    return model, nid_offset, eid_offset

def _mirror_nodes(model: BDF, plane: str='xz'):
    """
    Mirrors the GRIDs

    .. warning:: doesn't consider coordinate systems;
                  it could, but you'd need 20 new coordinate systems
    .. warning:: doesn't mirror SPOINTs, EPOINTs

    """
    iy, plane = _plane_to_iy(plane)

    nid_offset = max(model.nodes) if model.nodes else 0
    if model.spoints:
        nid_offset = max(max(model.spoints), nid_offset)
    if model.epoints:
        nid_offset = max(max(model.epoints), nid_offset)

    if model.nodes:
        for (nid, node) in sorted(model.nodes.items()):
            xyz = node.get_position()
            nid2 = nid + nid_offset
            xyz2 = xyz.copy()
            xyz2[iy] *= -1.
            model.add_grid(nid2, xyz2, cp=0, cd=node.cd, ps=node.ps, seid=node.seid)
    if model.spoints:
        for nid in list(model.spoints):
            nid2 = nid + nid_offset
            model.add_spoint(nid2)
    if model.epoints:
        for nid in list(model.epoints):
            nid2 = nid + nid_offset
            model.add_spoint(nid2)
    #if model.epoints:

    return nid_offset, plane

def _mirror_nodes_plane(model: BDF, mirror_model: BDF, plane,
                        use_nid_offset: bool=True) -> Tuple[int, str]:
    """
    Mirrors the GRIDs about an arbitrary plane

    Parameters
    ----------
    model : BDF
        ???
    mirror_model : BDF
        ???
    plane : str
        ???
    use_nid_offset : bool
        ???

    Returns
    -------
    nid_offset : int
        the node id offset
    plane : str
        the sorted plane; ZX -> xz

    .. warning:: doesn't consider coordinate systems;
                  it could, but you'd need 20 new coordinate systems
    .. warning:: doesn't mirror SPOINTs, EPOINTs

    https://mathinsight.org/distance_point_plane

    """
    nid_offset = 0

    if model.nodes:
        all_nodes, xyz_cid0 = model.get_xyz_in_coord_no_xref(
            cid=0, fdtype='float64', sort_ids=True)
        unused_cid = max(model.coords) + 1
        origin = plane[0, :]

        # just a triangle's normal vector
        n1 = plane[0, :]
        n2 = plane[1, :]
        n3 = plane[2, :]
        normal = np.cross(n2 - n1, n3 - n1)
        normal /= np.linalg.norm(normal)

        #cord2r = model.add_cord2r(cid, plane[0, :], plane[1, :], plane[2, :])
        #del model.coords[cid]
        #print(cord2r)

        #origin = cord2r.origin
        #normal = cord2r.beta()[2, :]
        #print('normal =', normal)
        vector = xyz_cid0 - origin
        assert xyz_cid0.shape == vector.shape, 'xyz_cid0.shape=%s; vector.shape=%s' % (xyz_cid0.shape, vector.shape)
        v_dot_n = vector * normal[np.newaxis, :]
        assert v_dot_n.shape == vector.shape, 'v_dot_n.shape=%s; vector.shape=%s' % (v_dot_n.shape, vector.shape)
        distance = np.linalg.norm(v_dot_n, axis=1)
        assert v_dot_n.shape[0] == len(distance), 'v_dot_n.shape=%s; distance.shape=%s' % (v_dot_n.shape, distance.shape)

        # we're some distance from the plane, but we don't know the
        # direction, so we take the max distance from the plane and
        # project it in the +normal direction and then check that
        # distance in comparison to the known distance
        #
        max_distance = distance.max()
        imax = np.where(distance == max_distance)[0][0]
        distance0 = distance[imax]
        xyz0 = xyz_cid0[imax, :] + distance0 * normal
        v0_dot_n = xyz0 * normal
        distance_plus = np.linalg.norm(v0_dot_n)

        if distance_plus > 1.1*distance0:
            xyz_cid0_2 = xyz_cid0 - 2 * distance[:, np.newaxis] * normal[np.newaxis, :]
        else:
            xyz_cid0_2 = xyz_cid0 + 2 * distance[:, np.newaxis] * normal[np.newaxis, :]

        if use_nid_offset:
            nid_offset = max(all_nodes)
        for nid, xyz2 in zip(all_nodes, xyz_cid0_2):
            node = model.nodes[nid]
            nid2 = nid + nid_offset
            mirror_model.add_grid(nid2, xyz2, cp=0, cd=node.cd, ps=node.ps, seid=node.seid)
    return nid_offset, plane

def _plane_to_iy(plane: str) -> Tuple[int, str]:
    """gets the index fo the mirror plane"""
    plane = plane.strip().lower()
    plane_sorted =  ''.join(sorted(set(plane)))
    if plane_sorted == 'yz':
        iy = 0
    if plane_sorted == 'xz':
        iy = 1
    elif plane_sorted == 'xy':
        iy = 2
    else:  # pragma: no cover
        raise NotImplementedError("plane=%r and must be 'yz', 'xz', or 'xy'." % plane)
    return iy, plane_sorted

def _mirror_elements(model: BDF, mirror_model: BDF,
                     nid_offset: int, use_eid_offset: bool=True, plane: str='xz') -> None:
    """
    Mirrors the elements

    Parameters
    ----------
    model : BDF
        the base model
    mirror_model : BDF
        the mirrored model
    nid_offset : int
        the node id offset
    use_eid_offset : bool; default=True
        what is this for???

    elements:
       0d : CELAS1, CELAS2, CELAS3, CELAS4, CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5
            CFAST, CBUSH, CBUSH1D
       1d:  CROD, CONROD, CTUBE, CBAR, CBEAM, CBEAM3
       2d : CTRIA3, CQUAD4, CTRIA6, CQUAD8, CQUAD, CTRIAR, CQUADR
       3d : ???
       missing : CVISC, CTRIAX, CTRIAX6, CQUADX, CQUADX8, CCONEAX
    rigid_elements:
       loaded: RBE2, RBE3, RBAR
       missing: RBAR1
    mass_elements:
       loaded: CONM2
       missing CONM1, CMASS1, CMASS2, CMASS3, CMASS4
    plotels:
       loaded: PLOTEL
       missing: ???

    Notes
    -----
    Doesn't handle CBAR/CBEAM offsets
    Doesn't handle CBEAM SPOINTs
    Do I need to invert the solids?

    """
    # why eid_offset not always calculated?
    eid_max_elements = 0
    eid_max_masses = 0
    eid_max_rigid = 0
    eid_max_plotels = 0
    if use_eid_offset:
        if model.elements:
            eid_max_elements = max(model.elements.keys())
        if model.masses:
            eid_max_masses = max(model.masses.keys())
        if model.rigid_elements:
            eid_max_rigid = max(model.rigid_elements.keys())
        if model.plotels:
            eid_max_plotels = max(model.plotels.keys())
    eid_offset = max(eid_max_elements, eid_max_masses, eid_max_rigid, eid_max_plotels)
    cid_offset = 1 if len(model.coords) == 1 else max(model.coords.keys())

    if model.elements:
        __mirror_elements(model, mirror_model, nid_offset, eid_offset, cid_offset,
                          plane=plane)

    if model.masses:
        __mirror_masses(model, mirror_model, nid_offset, eid_offset)

    if model.rigid_elements:
        __mirror_rigid_elements(model, mirror_model, nid_offset, eid_offset)

    if model.plotels:
        for eid, element in sorted(model.plotels.items()):
            eid_mirror = eid + eid_offset
            nodes = element.node_ids
            if element.type == 'PLOTEL':
                nodes = [node_id + nid_offset for node_id in nodes]
                mirror_model.add_plotel(eid_mirror, nodes)
            else:  # pragma: no cover
                mirror_model.log.warning('skipping:\n%s' % str(element))

    return eid_offset

def __mirror_elements(model: BDF, mirror_model: BDF,
                      nid_offset: int, eid_offset: int, cid_offset: int, plane: str='xz') -> None:
    """mirrors model.elements"""
    shells = {'CTRIA3', 'CQUAD4', 'CTRIAR', 'CQUADR'}
    shell_nones = {'CTRIA6', 'CQUAD8', 'CQUAD', }
    rods = {'CROD', 'CONROD', 'CTUBE'}
    spring_dampers = {
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
        'CVISC',
    }
    solids = {'CHEXA', 'CPENTA', 'CTETRA', 'CYPRAM'}
    generic_types = {
        'CSHEAR',
        'CTRIAX6',
        'CTRAX3', 'CTRAX6',
        'CQUADX4', 'CQUADX8',
        'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8', }
    generic_types_none = {
        'CHBDYP', 'CRAC2D', 'CRAC3D',
        'CQUADX',
        'CTRIAX',
        # acoustic
        'CHACAB'}

    def _set_nodes(element, nodes):
        try:
            element.nodes = nodes
        except AttributeError:
            print(element.get_stats())
            print(nodes)
            raise

    update_mcids = False
    element_cids = set([])
    etypes_skipped = set([])
    for eid, element in sorted(model.elements.items()):
        etype = element.type
        if etype in ['CHBDYG', 'CHBDYE']:
            continue

        nodes1 = element.node_ids
        #try:
            #nodes = [node_id + nid_offset for node_id in nodes]
        #except TypeError:
            #msg = 'cannot mirror %r eid=%s because None exists in nodes=%s' % (
                #element.type, eid, nodes)
            #model.log.warning(msg)
            #continue

        eid_mirror = eid + eid_offset
        fields = element.repr_fields()
        fields[1] = eid_mirror
        mirror_model.add_card(fields, etype)
        element2 = mirror_model.elements[eid_mirror]

        if etype in shells:
            nodes2 = [node_id + nid_offset for node_id in nodes1]
            element2.nodes = nodes2
            element.flip_normal() # nodes = nodes[::-1]
            if isinstance(element.theta_mcid, int) and update_mcids:
                element_cids.add(element.theta_mcid)
                element2.theta_mcid += cid_offset
        elif etype in solids:
            # what about inverting solids?
            nodes2 = [node_id + nid_offset if node_id is not None else None
                     for node_id in nodes1]
            etypes_skipped.add(etype)
            element2.nodes2 = nodes2
            element2.cross_reference(model)
            vol = element2.Volume()
            assert vol >= 0., vol
        elif etype in shell_nones:
            nodes2 = [node_id + nid_offset if node_id is not None else None
                     for node_id in nodes1]
            element2.nodes = nodes2
            element.flip_normal() # nodes = nodes[::-1]
            if isinstance(element.theta_mcid, int) and update_mcids:
                element_cids.add(element.theta_mcid)
                element2.theta_mcid += cid_offset
        elif etype in rods:
            nodes = [node_id + nid_offset for node_id in nodes1]
            element2.nodes = nodes
        elif etype in ['CBAR', 'CBEAM']:
            nodes2 = [node_id + nid_offset for node_id in nodes1]
            element2.nodes = nodes2
            g0 = element2.g0 + nid_offset if element2.g0 is not None else None
            element2.g0 = g0
        elif etype == 'CBEAM3':
            element2.ga = nodes1[0] + nid_offset
            element2.gb = nodes1[1] + nid_offset
            g0 = element2.g0 + nid_offset if element2.g0 is not None else None
            element2.g0 = g0

        elif etype == 'CGAP':
            #nodes = [node_id + nid_offset if node_id is not None else None
                     #for node_id in nodes]
            #_set_nodes(element2, nodes)
            ga = element2.ga + nid_offset if element2.ga is not None else None
            gb = element2.gb + nid_offset if element2.gb is not None else None
            g0 = element2.g0 + nid_offset if element2.g0 is not None else None
            element2.ga = ga
            element2.gb = gb
            element2.g0 = g0
            if isinstance(element.cid, int) and update_mcids:
                element_cids.add(element.cid)
                element2.cid += cid_offset
        elif etype == 'CBUSH':
            nodes2 = [node_id + nid_offset for node_id in nodes1]
            g0 = element2.g0 + nid_offset if element2.g0 is not None else None
            element2.g0 = g0
            element2.nodes = nodes2
            #  TODO: ocid
            if isinstance(element.cid, int) and update_mcids:
                element_cids.add(element.cid)
                element2.cid += cid_offset
        elif etype == 'CBUSH1D':
            ga = element2.ga + nid_offset if element2.ga is not None else None
            gb = element2.gb + nid_offset if element2.gb is not None else None
            element2.ga = ga
            element2.gb = gb
            if isinstance(element.cid, int) and update_mcids:
                element_cids.add(element.cid)
                element2.cid += cid_offset
        elif etype == 'CFAST':
            ga = element2.ga + nid_offset if element2.ga is not None else None
            gb = element2.gb + nid_offset if element2.gb is not None else None
            gs = element2.gs + nid_offset if element2.gs is not None else None
            element2.ga = ga
            element2.gb = gb
            element2.gs = gs

        elif etype == 'CCONEAX':
            pass
        elif etype in spring_dampers:
            nodes2 = [node_id + nid_offset if node_id is not None else None for node_id in nodes1]
            _set_nodes(element2, nodes2)
            #print(nodes)
            #element2.nodes = nodes
        elif etype == 'GENEL':
            element2.ul = element2.ul + nid_offset
            element2.ud = element2.ud + nid_offset
        elif etype in generic_types:
            try:
                nodes2 = [node_id + nid_offset for node_id in nodes1]
            except TypeError:  # pragma: no cover
                print(element.get_stats())
                raise
            _set_nodes(element2, nodes2)
        elif etype in generic_types_none:
            nodes2 = [node_id + nid_offset if node_id is not None else None for node_id in nodes1]
            _set_nodes(element2, nodes2)
        else:
            etypes_skipped.add(etype)
            try:
                nodes2 = [node_id + nid_offset for node_id in nodes1]
                element2.nodes = nodes2
            except (AttributeError, TypeError):  # pragma: no cover
                print(element.get_stats())
                raise
    if etypes_skipped:
        etypes = list(etypes_skipped)
        etypes.sort()
        model.log.warning(f'Verify element_types={etypes}')
    if element_cids:
        cids = list(element_cids)
        cids.sort()
        model.log.warning(f'verify material coordinate systems %s' % cids)
        _asymmetrically_mirror_coords(model, element_cids, cid_offset, plane)
    return

def __mirror_masses(model, mirror_model, nid_offset, eid_offset):
    """mirrors model.masses"""
    for eid, element in sorted(model.masses.items()):
        eid_mirror = eid + eid_offset

        if element.type == 'CONM2':
            old_nid = element.nid
            new_nid = old_nid + nid_offset
            mass = element.mass
            cid = element.cid
            X = element.X
            I = element.I
            mirror_model.add_conm2(eid_mirror, new_nid, mass, cid=cid, X=X, I=I, comment='')
        elif element.type == 'CONM1':
            new_nid = element.node_ids[0] + nid_offset
            mass_matrix = element.mass_matrix
            cid = element.Cid()
            mirror_model.add_conm1(eid_mirror, new_nid, mass_matrix, cid=cid, comment='')
        elif element.type == 'CMASS1':
            old_node_ids = [element.G1(), element.G2()]
            new_nids = [nid + nid_offset if nid is not None else None
                        for nid in old_node_ids]
            mirror_model.add_cmass1(eid_mirror, element.Pid(), new_nids,
                                    c1=element.c1, c2=element.c2, comment='')
        elif element.type == 'CMASS2':
            old_node_ids = [element.G1(), element.G2()]
            new_nids = [nid + nid_offset if nid is not None else None
                        for nid in old_node_ids]
            mirror_model.add_cmass2(eid_mirror, element.mass, new_nids,
                                    c1=element.c1, c2=element.c2, comment='')
        elif element.type == 'CMASS3':
            # this uses SPOINTs, so we don't need as fancy checking as CMASS1 and CMASS2
            new_nids = [nid + nid_offset if nid is not None else None
                        for nid in element.node_ids]
            mirror_model.add_cmass3(eid_mirror, element.Pid(), new_nids, comment='')
        elif element.type == 'CMASS4':
            # this uses SPOINTs, so we don't need as fancy checking as CMASS1 and CMASS2
            new_nids = [nid + nid_offset if nid is not None else None
                        for nid in element.node_ids]
            mirror_model.add_cmass4(eid_mirror, element.mass, new_nids, comment='')
        else:  # pragma: no cover
            mirror_model.log.warning('skipping mass element:\n%s' % str(element))
    del eid_mirror

def __mirror_rigid_elements(model: BDF, mirror_model: BDF,
                            nid_offset: int, eid_offset: int) -> None:
    """mirrors model.rigid_elements"""
    for eid, rigid_element in sorted(model.rigid_elements.items()):
        eid_mirror = eid + eid_offset
        if rigid_element.type == 'RBE2':
            Gmi_node_ids = rigid_element.Gmi_node_ids
            Gn = rigid_element.Gn()
            Gijs = None
            ref_grid_id = None
        elif rigid_element.type == 'RBE3':
            Gmi_node_ids = rigid_element.Gmi_node_ids
            Gijs = rigid_element.Gijs
            ref_grid_id = rigid_element.ref_grid_id
            Gn = None
        elif rigid_element.type == 'RROD':
            node_ids_mirror = [node_id + nid_offset for node_id in rigid_element.nodes]
            mirror_model.add_rrod(
                eid_mirror, node_ids_mirror,
                cma=rigid_element.cma, cmb=rigid_element.cmb,
                alpha=rigid_element.alpha, comment='')
            continue
        elif rigid_element.type == 'RBAR':
            node_ids_mirror = [node_id + nid_offset for node_id in rigid_element.nodes]
            mirror_model.add_rbar(
                eid_mirror, node_ids_mirror,
                rigid_element.cna, rigid_element.cnb,
                rigid_element.cma, rigid_element.cmb,
                alpha=rigid_element.alpha, comment='')
            continue
        elif rigid_element.type == 'RBAR1':
            node_ids = [node_id + nid_offset for node_id in
                        [rigid_element.ga, rigid_element.gb]]
            mirror_model.add_rbar1(
                eid_mirror, node_ids, rigid_element.cb, alpha=rigid_element.alpha, comment='')
            continue
        elif rigid_element.type == 'RBE1':
            Gni = [node_id + nid_offset for node_id in rigid_element.Gni]
            Gmi = [node_id + nid_offset for node_id in rigid_element.Gmi]
            mirror_model.add_rbe1(
                eid_mirror,
                Gni, rigid_element.Cni,
                Gmi, rigid_element.Cmi,
                alpha=rigid_element.alpha, comment='')
            continue
        else:
            model.log.warning('_write_elements_symmetric: %s not implemented' % rigid_element.type)
            continue
            #raise NotImplementedError(msg)

        #if rigid_element.type in ['RBE2', 'RBE3']:
        Gmi_node_ids_mirror = [node_id + nid_offset for node_id in Gmi_node_ids]
        if Gn:
            Gn_mirror = Gn + nid_offset
        if Gijs:
            Gijs_mirror = [[node_id + nid_offset for node_id in nodes] for nodes in Gijs]
        if ref_grid_id:
            ref_grid_id_mirror = ref_grid_id + nid_offset

        if rigid_element.type == 'RBE2':
            mirror_model.add_rbe2(eid_mirror, Gn_mirror, rigid_element.cm,
                                  Gmi_node_ids_mirror)
        elif rigid_element.type == 'RBE3':
            mirror_model.add_rbe3(
                eid_mirror, ref_grid_id_mirror, rigid_element.refc, rigid_element.weights,
                rigid_element.comps, Gijs_mirror
            )
        else:  # pragma: no cover
            mirror_model.log.warning('skipping:\n%s' % str(rigid_element))


def _mirror_loads(model: BDF, nid_offset: int=0, eid_offset: int=0) -> None:
    """
    Mirrors the loads.  A mirrored force acts in the same direction.

    Considers:
     - PLOAD4
        - no coordinate systems (assumes cid=0)
     - FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2
     - PLOAD, PLOAD2
     - TEMP, QVOL, QHBDY, QBDY1, QBDY2, QBDY3

    """
    # unmirrorable loads
    skip_loads = {
        'FORCEAX', 'TEMPAX', 'PRESAX',
    }
    for unused_load_id, loads in model.loads.items():
        for load in loads:
            loads_new = []
            load_type = load.type
            if load_type == 'PLOAD4':
                g1 = None
                g34 = None
                if load.g1 is not None:
                    g1 = load.g1 + nid_offset
                if load.g34 is not None:
                    g34 = load.g34 + nid_offset

                eids = [eid + eid_offset for eid in load.eids]
                load = PLOAD4(
                    load.sid, eids, load.pressures, g1, g34,
                    cid=load.cid, nvector=load.nvector,
                    surf_or_line=load.surf_or_line,
                    line_load_dir=load.line_load_dir, comment='')
                loads_new.append(load)

            elif load_type == 'FORCE':
                load = FORCE(load.sid, load.node + nid_offset, load.mag, load.xyz,
                              cid=load.cid, comment='')
                loads_new.append(load)
            elif load_type == 'FORCE1':
                load = FORCE1(load.sid, load.node + nid_offset, load.mag,
                              load.g1 + nid_offset, load.g2 + nid_offset, comment='')
                loads_new.append(load)
            elif load_type == 'FORCE2':
                load = FORCE2(load.sid, load.node + nid_offset, load.mag,
                              load.g1 + nid_offset, load.g2 + nid_offset,
                              load.g3 + nid_offset, load.g4 + nid_offset, comment='')
                loads_new.append(load)

            elif load_type == 'MOMENT':
                load = MOMENT(load.sid, load.node + nid_offset, load.mag, load.xyz,
                              cid=load.cid, comment='')
                loads_new.append(load)
            elif load_type == 'MOMENT1':
                load = MOMENT1(load.sid, load.node + nid_offset, load.mag,
                               load.g1 + nid_offset, load.g2 + nid_offset, comment='')
                loads_new.append(load)
            elif load_type == 'MOMENT2':
                load = MOMENT2(load.sid, load.node + nid_offset, load.mag,
                               load.g1 + nid_offset, load.g2 + nid_offset,
                               load.g3 + nid_offset, load.g4 + nid_offset, comment='')
                loads_new.append(load)

            elif load_type == 'PLOAD':
                nodes = [nid + nid_offset for nid in load.nodes]
                load = PLOAD(load.sid, load.pressure, nodes, comment='')
                loads_new.append(load)
            elif load_type == 'PLOAD2':
                eids = [eid + eid_offset for eid in load.eids]
                load = PLOAD2(load.sid, load.pressure, eids, comment='')
                loads_new.append(load)

            elif load_type == 'QVOL':
                elements = [eid + eid_offset for eid in load.elements]
                load = QVOL(load.sid, load.qvol, nid_offset + load.control_point, elements)
                loads_new.append(load)
            elif load_type == 'QHBDY':
                grids = [nid + nid_offset for nid in load.grids]
                load = QHBDY(load.sid, load.flag, load.q0, grids, af=load.af)
                loads_new.append(load)
            elif load_type == 'QBDY1':
                eids = [eid + eid_offset for eid in load.eids]
                load = QBDY1(load.sid, load.qflux, eids)
                loads_new.append(load)
            elif load_type == 'QBDY2':
                load = QBDY2(load.sid, load.eid + eid_offset, load.qfluxs, comment='')
                loads_new.append(load)
            elif load_type == 'QBDY3':
                eids = [eid + eid_offset for eid in load.eids]
                load = QBDY3(load.sid, load.q0, load.cntrlnd + nid_offset, eids)
                loads_new.append(load)

            elif load_type == 'TEMP':
                temperatures = {}
                for nid, temp in load.temperatures.items():
                    temperatures[nid + nid_offset] = temp
                load = TEMP(load.sid, temperatures)
                loads_new.append(load)
            elif load_type == 'GRAV':
                pass
            elif load_type in skip_loads:
                continue
            else:  # pragma: no cover
                model.log.warning('skipping:\n%s' % load.rstrip())
        if loads_new:
            loads += loads_new

def _mirror_aero(model: BDF, nid_offset: int, plane: str='xz') -> None:
    """
    Mirrors the aero cards

    Considers:
     - AEROS
      - doesn't consider sideslip
     - CAERO1
      - doesn't consider sideslip
      - considers Cp
      - considers lchord/lspan/nchord/nspan
     - SPLINE1
       - handle boxes
     - SET1
       - handle nodes
     - AELIST
       - handle boxes
     - AESURF
       - only supports names of length 7 or less (appends an M to the name)
       - handles AELIST
       - doesn't handle coords well
       - doesn't handle second AESURF

    Doesnt consider:
     - AERO
     - AEFORCE
     - AEPRES
     - CAERO2/3/4/5
     - PAERO1/2/3/4/5
     - AESURFS

    """
    is_aero = False
    aero_cids_set = set([])
    if model.aeros is not None:
        is_aero = True
        aeros = model.aeros
        # The ACSID must be a rectangular coordinate system.
        # Flow is in the positive x-direction (T1).
        if aeros.acsid > 0:
            model.log.error('Sideslip coordinate system (ACSID) mirroring is not supported')

        # REFB should be full span, even on half-span models
        # REFS should be half area on half-span models
        aeros.sref *= 2.
        if plane == 'xz':
            aeros.sym_xz = 0
        elif plane == 'yz':
            aeros.sym_yz = 0
        else:
            model.log.error('not mirroring plane %r; only xz, yz' % plane)

    caero_id_offset = 0
    if len(model.caeros):
        is_aero = True
        caero_id_max = max(model.caero_ids)
        caero_id_offset = np.max(model.caeros[caero_id_max].box_ids.flat)

        caeros = []
        for unused_caero_id, caero in model.caeros.items():
            if caero.type == 'CAERO1':
                # the AEFACTs are assumed to be the same on the left and right side
                # I think the spanwise direction will be fine
                lchord = caero.lchord
                nchord = caero.nchord
                lspan = caero.lspan
                nspan = caero.nspan
                p1, p4 = caero.get_leading_edge_points()
                p1 = p1.copy()
                p4 = p4.copy()
                x12 = caero.x12
                x43 = caero.x43
                if plane == 'xz': # flip the y
                    p1[1] *= -1.
                    p4[1] *= -1.
                elif  plane == 'xy':
                    p1[2] *= -1.
                    p4[2] *= -1.
                else:  # pragma: no cover
                    raise NotImplementedError('plane=%r not supported in CAERO1' % plane)
                eid2 = caero.eid + caero_id_offset
                caero_new = CAERO1(eid2, caero.pid, caero.igroup,
                                   p1, x12, p4, x43,
                                   cp=0,
                                   nspan=nspan, lspan=lspan,
                                   nchord=nchord, lchord=lchord,
                                   comment='')

                # we flip the normal so if we ever use W2GJ it's going to be consistent
                caero_new.flip_normal()
                caeros.append(caero_new)
            else:  # pragma: no cover
                model.log.error('skipping (only supports CAERO1):\n%s' % caero.rstrip())

        for caero in caeros:
            model._add_caero_object(caero)

    nsplines = len(model.splines)
    sets_max = max(model.sets) if len(model.sets) else 0
    if caero_id_offset == 0 and nsplines:
        model.log.error("cant mirror splines because CAEROs don't exist...")
    elif nsplines and sets_max == 0:
        model.log.error("cant mirror splines because SET1/3 don't exist...")
    elif nsplines:
        is_aero = True
        splines = []
        spline_sets_to_duplicate = []
        spline_max = max(model.splines)
        for unused_spline_id, spline in model.splines.items():
            if spline.type == 'SPLINE1':
                eid = spline.eid + spline_max
                caero = spline.caero + caero_id_offset
                method = spline.method
                usage = spline.usage
                box1 = spline.box1 + caero_id_offset
                box2 = spline.box2 + caero_id_offset
                setg = spline.setg + sets_max
                dz = spline.dz
                melements = spline.melements
                nelements = spline.nelements
                spline_new = SPLINE1(eid, caero, box1, box2, setg, dz=dz,
                                     method=method, usage=usage,
                                     nelements=nelements, melements=melements, comment='')
                splines.append(spline_new)
                spline_sets_to_duplicate.append(spline.setg)
            else:  # pragma: no cover
                model.log.error('skipping (only support SPLINE1):\n%s' % spline.rstrip())

        #print("spline_sets_to_duplicate =", spline_sets_to_duplicate)
        msg = ', which is required to mirror:\n%s' % spline.rstrip()

        sets_to_add = []
        for set_id in spline_sets_to_duplicate:
            set_card = model.Set(set_id, msg=msg)
            if set_card.type == 'SET1':
                sid = set_card.sid + sets_max
                ids = [nid + nid_offset for nid in set_card.ids]
                is_skin = set_card.is_skin
                set_card = SET1(sid, ids, is_skin=is_skin, comment='')
                sets_to_add.append(set_card)
            else:  # pragma: no cover
                model.log.error('skipping (only support SET1):\n%s' % set_card.rstrip())

        for spline in splines:
            model._add_spline_object(spline)
        for set_card in sets_to_add:
            model._add_set_object(set_card)

    aelist_id_offset = 0
    if len(model.aelists):
        is_aero = True
        aelist_id_offset = max(model.aelists)

    cid_offset = max(model.coords) if len(model.coords) > 1 else 0
    if len(model.aesurf):
        is_aero = True
        aesurf_id_offset = max(model.aesurf)
        for aelist_id, aelist in sorted(model.aelists.items()):
            aelist_id_new = aelist_id + aelist_id_offset
            elements = [eid + caero_id_offset for eid in aelist.elements]
            model.add_aelist(aelist_id_new, elements, comment='')

        for aesurf_id, aesurf in sorted(model.aesurf.items()):
            # TODO: doesn't handle cid2/aelist2
            aesurf_id_new = aesurf_id + aesurf_id_offset
            label = aesurf.label + 'M'
            cid1 = aesurf.cid1 + cid_offset
            alid1 = aesurf.alid1 + aelist_id_offset
            if cid_offset > 0:
                aero_cids_set.add(aesurf.cid1)

            #cid2 = None
            #alid2 = None
            if aesurf.cid2:
                model.log.warning("skipping aesurf='{aesurf.label}' second cid/aelist")
                # combine this into the first coordinate system
                #
                #  don't mirror the coordinate system because it's antisymmetric?
                #cid2 = aesurf.cid1 + cid_offset * 2
                # how is the second coord and aelist made?
                #aero_cids_set.add(aesurf.cid1)
                #alid2 = aesurf.alid1 + aelist_id_offset * 2

            model.add_aesurf(aesurf_id_new, label, cid1, alid1,
                             cid2=None, alid2=None,
                             eff=aesurf.eff, ldw=aesurf.ldw,
                             crefc=aesurf.crefc, crefs=aesurf.crefs,
                             pllim=aesurf.pllim, pulim=aesurf.pulim,
                             hmllim=aesurf.hmllim, hmulim=aesurf.hmulim,
                             tqllim=aesurf.tqllim, tqulim=aesurf.tqulim, comment='')

    if is_aero:
        _asymmetrically_mirror_coords(model, aero_cids_set, cid_offset, plane=plane)
    model.pop_parse_errors()

def _asymmetrically_mirror_coords2(model: BDF,
                                  cids_nominal_set: Set[int],
                                  cid_offset: int, plane: str='xz') -> None:
    """
    We'll invert i, but not j, which will invert k.

    This will allow us to mirror material coordinate systems (defined in x).
    Additionally, we can mirror CGAP coordinate systems (defined in x) as well
    as some arbitrary y-direction.  We'll try to make y look mirrored.
    """
    return

def _asymmetrically_mirror_coords(model: BDF,
                                  cids_nominal_set: Set[int],
                                  cid_offset: int, plane: str='xz') -> None:
    """we'll leave i the same, flip j, and invert k"""
    # doesn't handle CORD1x

    coord_map = {
        'CORD1R': CORD1R,
        'CORD1C': CORD1C,
        'CORD1S': CORD1S,

        'CORD2R': CORD2R,
        'CORD2C': CORD2C,
        'CORD2S': CORD2S,
    }
    cids = list(cids_nominal_set)
    cids.sort()
    model.log.debug(f'asymmetrically mirroring coordinate systems')
    for cid in cids:

        coord = model.Coord(cid)
        if cid == 0:
            continue
        cid_new = cid + cid_offset
        model.log.debug(f'  cid={cid} -> {cid_new}')

        if coord.type in {'CORD2R', 'CORD2C', 'CORD2S'}:
            i, j, unused_k = coord.beta().copy()
            origin = coord.origin.copy()
            if plane == 'xz':
                origin[1] *= -1
                j[1] *= -1
            elif plane == 'xy':
                origin[2] *= -1
                j[2] *= -1
            else:
                model.log.warning('skipping coord_id=%s' % coord.cid)
                return
            #k = np.cross(i, j)
            coord_obj = coord_map[coord.type]  # CORD2R/C/S
            coord_new = coord_obj.add_ijk(cid_new, origin=origin, i=i, j=j, k=None,
                                          rid=0, comment='')
        else:
            model.log.warning('skipping coord_id=%s' % coord.cid)
            continue
        model.coords[cid_new] = coord_new
    return
