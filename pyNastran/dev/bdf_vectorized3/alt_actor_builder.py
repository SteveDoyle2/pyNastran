"""Creates a series of Nastran-themed vtkActors"""
from __future__ import annotations
from itertools import count
from typing import TYPE_CHECKING

import numpy as np
from pyNastran.gui.vtk_common_core import vtkPoints
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid, vtkCellArray, vtkVertex

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.gui.utils.vtk.vectorized_geometry import build_vtk_geometry
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_type,
    create_vtk_cells_of_constant_element_types,
    numpy_to_vtk_points)
from pyNastran.gui.qt_files.colors import (
    RED_FLOAT, BLUE_FLOAT,
    LIGHT_GREEN_FLOAT, PINK_FLOAT, GREEN_FLOAT, #PURPLE_FLOAT,
    YELLOW_FLOAT, # ORANGE_FLOAT,
)
from pyNastran.bdf.cards.aero.utils import (
    points_elements_from_quad_points, # create_axisymmetric_body,
)
from pyNastran.femutils.utils import hstack_unique

#from pyNastran.gui.gui_objects.gui_result import GuiResult# , NormalResult
#from pyNastran.gui.gui_objects.displacements import ForceTableResults, ElementalTableResults

from pyNastran.dev.bdf_vectorized3.cards.constraints import spc_cards_to_nid_dof
from pyNastran.dev.bdf_vectorized3.cards.elements.shell_coords import (
    get_shell_element_coordinate_system, get_shell_material_coordinate_system)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import (
        #SPC, SPC1,
        AECOMP, AECOMPL, SET1, RBE2, RBE3 #, GRID
    )
    from pyNastran.gui.main_window import MainWindow
    from pyNastran.gui.gui_objects.settings import NastranSettings # Settings,
    from .nastran_io3 import Nastran3 as NastranIO


def create_alt_spcs(gui: MainWindow,
                    model: BDF,
                    grid_id: np.ndarray,
                    xyz_cid0: np.ndarray) -> None:
    """
    Creates SPC/SPC1 constraints:
     - Single DOF Set
       - All SPCs id=1; dof=123456
     - Multiple DOF Sets
       - All SPCs id=2
       - SPCs id=2; dof=123
       - SPCs id=2; dof=456

    Creates SPCADD constraints:
     - Single DOF Set
       - All SPCADD id=1; dof=123456
     - Multiple DOF Sets
       - All SPCADD id=2
       - SPCADD id=2; dof=123
       - SPCADD id=2; dof=456

    TODO: no support for enforced
    """
    spcadd = model.spcadd
    cards_all = (model.spc, model.spc1)
    cards = [card for card in cards_all if len(card)]
    ncards = sum([len(card) for card in cards])

    nastran_settings: NastranSettings = gui.settings.nastran_settings
    if ncards == 0 or not nastran_settings.is_constraints:  # and len(spcadd) == 0:
        return

    spc_ids_list = []
    for card in cards:
        spc_ids_list.append(card.spc_id)
    spc_ids = np.hstack(spc_ids_list)
    uspc_id = np.unique(spc_ids)

    for spc_id in uspc_id:
        is_failed, nids, comp = spc_cards_to_nid_dof(spc_id, cards)
        if is_failed:
            continue
        unids_all = np.unique(nids)
        ucomp_all = np.unique(comp)

        inid = np.searchsorted(grid_id, unids_all)
        xyz_cid0s = xyz_cid0[inid, :]
        if len(ucomp_all) == 1:
            # single dof; don't need "ALL SPCs" cause this
            # case all has the same DOFs
            ucomp = ucomp_all[0]
            name = f'All SPCs id={spc_id:g}; dof={ucomp}'
            _build_dots(gui, name, xyz_cid0s)
        else:
            # all nodes
            name = f'All SPCs id={spc_id:g}'
            _build_dots(gui, name, xyz_cid0s)

            # SPCs by DOF
            for ucompi in ucomp_all:
                name = f'SPCs id={spc_id:g}; DOF={ucompi}'
                icomp = np.where(comp == ucompi)[0]
                unids = nids[icomp]
                inid = np.searchsorted(grid_id, unids)
                xyz_cid0s = xyz_cid0[inid, :]
                _build_dots(gui, name, xyz_cid0s)

    if len(spcadd):
        reduced_spc_dict = spcadd.get_reduced_spcs()
        spcadd_node_component_dict = spcadd.get_reduced_node_component(reduced_spc_dict)
        #get_reduced_spcs
        #spcs_dict = model.spcadd.get_spcs_by_spc_id()
        for spcadd_id, (nids, comp) in spcadd_node_component_dict.items():
            unids_all = np.unique(nids)
            ucomp_all = np.unique(comp)

            inid = np.searchsorted(grid_id, unids_all)
            xyz_cid0s = xyz_cid0[inid, :]
            if len(ucomp_all) == 1:
                # single dof; don't need "ALL SPCs" cause this
                # case all has the same DOFs
                ucomp = ucomp_all[0]
                name = f'All SPCADD id={spcadd_id:g}; dof={ucomp}'
                _build_dots(gui, name, xyz_cid0s)
            else:
                # all nodes
                name = f'All SPCADD id={spcadd_id:g}'
                _build_dots(gui, name, xyz_cid0s)

                # SPCs by DOF
                for ucompi in ucomp_all:
                    name = f'SPCADD id={spcadd_id:g}; DOF={ucompi}'
                    icomp = np.where(comp == ucompi)[0]
                    unids = nids[icomp]
                    inid = np.searchsorted(grid_id, unids)
                    xyz_cid0s = xyz_cid0[inid, :]
                    _build_dots(gui, name, xyz_cid0s)
        #x = 1
    pass

def create_alt_conm2_grids(gui: MainWindow,
                           model: BDF,
                           grid_id: np.ndarray,
                           xyz_cid0: np.ndarray) -> None:
    """
    Creates:
     - global CG actor
     - idividual CONM2 actors

    """
    nastran_settings: NastranSettings = gui.settings.nastran_settings
    if not nastran_settings.is_mass:
        return
    log = model.log
    try:
        mass_total, cg, inertia = model.inertia_sum()
    except Exception as e:
        log.error(f'cannot build mass actors\n{str(e)}')
        #raise
        mass_total = 0.

    if mass_total != 0.0:
        name = f'CG mass={mass_total:g} xyz=[{cg[0]:g}, {cg[0]:g}, {cg[2]:g}]'
        log.info(f'mass = {mass_total:g}')
        log.info(f'cg = {cg}')
        #print('inertia =', inertia)
        _build_dot(gui, name, cg)

    element = model.conm2
    if element.n == 0:
        return

    try:
        mass = element.mass()
        centroid = element.centroid()
    except AssertionError:
        return

    # total - cg & all
    conm2_mass_total = mass.sum()
    cg_mass_total = np.full(3, np.nan, dtype=mass.dtype)
    if conm2_mass_total != 0:
        cg_mass_total = (mass[:, np.newaxis] * centroid).sum(axis=0) / conm2_mass_total
        assert len(cg_mass_total) == 3, cg_mass_total
        name = f'All CONM2s CG mass={conm2_mass_total:g} xyz=[{cg_mass_total[0]:g}, {cg_mass_total[1]:g}, {cg_mass_total[2]:g}]'
        _build_dot(gui, name, cg_mass_total)

        name = f'All CONM2s mass={conm2_mass_total:g}'
        _build_dots(gui, name, centroid)

    # individual
    for eid, nid, massi, xyzi in zip(element.element_id, element.node_id, mass, centroid):
        name = f'CONM2 {eid} nid={nid} mass={massi:g}  xyz=[{xyzi[0]:g}, {xyzi[1]:g}, {xyzi[2]:g}]'
        _build_dot(gui, name, xyzi)

def create_alt_rbe3_grids(gui: MainWindow,
                          model: BDF,
                          grid_id: np.ndarray,
                          xyz_cid0: np.ndarray) -> None:
    """
    Creates global actors for:
      - RBE3 independents
      - RBE3 reference dependents
      - RBE3 UM dependents
      - RBE3 lines

    TODO: add individual support

    """
    elem = model.rbe3
    nastran_settings: NastranSettings = gui.settings.nastran_settings
    if elem.n == 0 or not nastran_settings.is_rbe:
        return

    i_independent = np.searchsorted(grid_id, elem.independent_nodes)
    xyz_independents = xyz_cid0[i_independent, :]

    name = f'RBE3 independents'
    _build_dots(gui, name, xyz_independents, color=BLUE_FLOAT)

    name = f'RBE3 Reference dependents'
    inid_all_dependent = np.searchsorted(grid_id, elem.ref_grid)
    try:
        ref_xyz_dependents = xyz_cid0[inid_all_dependent, :]
    except IndexError:
        return
    _build_dots(gui, name, ref_xyz_dependents, color=RED_FLOAT)

    all_dep_nodes = elem.dependent_nodes
    if len(all_dep_nodes):
        name = f'RBE3 UM dependents'
        inid_all_dependent = np.searchsorted(grid_id, all_dep_nodes)
        all_xyz_dependents = xyz_cid0[inid_all_dependent, :]
        _build_dots(gui, name, all_xyz_dependents, color=RED_FLOAT)

    #----------------------
    ref_igrids = np.searchsorted(grid_id, elem.ref_grid)
    independent_inodes = np.searchsorted(grid_id, elem.independent_nodes)
    dependent_inodes = np.searchsorted(grid_id, elem.dependent_nodes)
    iweight0_for_grids = 0
    #iGij0 = 0
    igrid0 = 0
    for eid, ref_igrid_dep, refc, nweight, iweight, \
            idependent, ref_xyz_dependent in zip(elem.element_id, ref_igrids, elem.ref_component,
                                                 elem.nweight, elem.iweight,
                                                 elem.idependent, ref_xyz_dependents):
        #idim0, idim1 = idim
        #dep_nodes = dep_node_str[idim0:idim1].tolist()  # Gmi

        # independent
        iweight0, iweight1 = iweight
        #weights = self.weight[iweight0:iweight1]
        comps = elem.independent_dofs[iweight0:iweight1]

        name = f'RBE3_eid={eid} ddof={refc} idof={ref_igrid_dep}'
        Gijs = []
        for i in range(nweight):
            ngrid = elem.ngrid_per_weight[iweight0_for_grids+i]
            igrid1 = igrid0 + ngrid
            Gij = independent_inodes[igrid0:igrid1].tolist()
            Gijs.append(Gij)
            igrid0 = igrid1

        nindependent = len(Gijs)
        ref_igridsi = np.ones(nindependent) * ref_igrid_dep
        lines = np.column_stack([Gijs, ref_igridsi])

        iweight0_for_grids += nweight
        #Gij = self.independent_nodes[iweight0:iweight1]
        #igrid = make_idim(len(Gij), ngrid_per_weight)

        # dependent
        idependent0, idependent1 = idependent
        Gmi = dependent_inodes[idependent0:idependent1]
        i_independent = np.searchsorted(grid_id, Gmi)
        xyz_independentsi = xyz_cid0[i_independent, :]

        #Cmi = elem.dependent_dofs[idependent0:idependent1]
        assert len(ref_xyz_dependent) == 3, ref_xyz_dependent

        ndependent = idependent1 - idependent0
        if ndependent:
            ref_xyz_dependentsi = np.array([ref_xyz_dependent] * ndependent)
            assert xyz_independentsi.shape == ref_xyz_dependentsi.shape, (xyz_independentsi.shape, ref_xyz_dependentsi.shape)
            xyz_lines = np.vstack([xyz_independentsi, ref_xyz_dependentsi])
            _build_lines(gui, name, xyz_lines, node_lines, color=LIGHT_GREEN_FLOAT)

    #_build_rbe3_vtk_lines(gui, elem, xyz_independents, ref_xyz_dependents)
    #for eid, ref_dep_node, ref_dep_xyz, dep_component, (iweight0, iweight1) in zip(elem.element_id,
                                                                                   #elem.ref_grid, ref_xyz_dependents,
                                                                                   #elem.ref_component,
                                                                                   #elem.iweight):
        ##ind_nodes = elem.independent_nodes[iweight0:iweight1]
        ##inid_independent = np.searchsorted(grid_id, ind_nodes)
        ##xyz_independents = xyz_cid0[inid_independent, :]

        #name = f'RBE3 {eid} ref dependent (dof={ind_component})'
        #_build_dot(gui, name, ref_dep_xyz, color=BLUE_FLOAT)

        #name = f'RBE3 {eid} dependent'
        #_build_dots(gui, name, xyz_dependents, color=RED_FLOAT)

def _build_rbe3_vtk_lines(gui, elem: RBE3,
                          xyz_independents: np.ndarray,
                          ref_xyz_dependents: np.ndarray) -> None:
    """not quite right..."""
    pass
    #lines: list[tuple[int, int]] = []
    #i0 = len(elem.element_id)
    #assert i0 == len(ref_xyz_dependents)
    #xyz_lines = np.vstack([ref_xyz_dependents,
                           #xyz_independents])

    #iweight0_for_grids = 0
    #igrid0 = 0
    #for ielem, eid, ref_dep_node, ref_xyz, dep_component, nweight, (iweight0, iweight1) in zip(
        #count(),
        #elem.element_id,
        #elem.ref_grid, ref_xyz_dependents,
        #elem.ref_component,
        #elem.nweight, elem.iweight):

        #Gijs = []
        #for i in range(nweight):
            #ngrid = elem.ngrid_per_weight[iweight0_for_grids+i]
            #igrid1 = igrid0 + ngrid
            #Gij = elem.independent_nodes[igrid0:igrid1]
            ##ind_nodesi = elem.independent_nodes[iweight0:iweight1]
            #Gijs.append(Gij)
            #igrid0 = igrid1
        #iweight0_for_grids += nweight

        #ind_nodes = np.hstack(Gijs)
        #print('ind_nodes =', eid, iweight0_for_grids, ind_nodes)
        #nind_nodes = len(ind_nodes)

        ## initialize to the independent node
        #nodes = np.full((nind_nodes, 2), i, dtype='int32')
        #n2 = np.arange(i0, i0+nind_nodes, dtype='int32')
        #nodes[:, 1] = n2

        ## we leave off the +1 cause we're starting at
        #i0 += nind_nodes
        #lines.append(nodes)
    #assert i0 == len(xyz_lines), f'i0={i0} xyz_lines.shape={xyz_lines.shape} nref={len(ref_xyz_dependents)} nind={len(xyz_independents)}'
    #node_lines = np.vstack(lines)
    #name = f'RBE3 lines'
    #_build_lines(gui, name, xyz_lines, node_lines, color=LIGHT_GREEN_FLOAT)

def create_alt_axes(self: NastranIO,
                    gui: MainWindow,
                    model: BDF,
                    grid_id: np.ndarray,
                    xyz_cid0: np.ndarray):
    """
    Creates:
     - shell element coordinate systems
     - shell material coordinate systems

    Creates orientation vectors for:
     - CBAR
     - CBEAM
     - CBUSH
     - CGAP

    """
    assert gui is not None
    if not hasattr(gui, 'bar_eids') or gui.bar_eids is None:
        gui.bar_eids = {}
    if not hasattr(gui, 'bar_lines') or gui.bar_lines is None:
        gui.bar_lines = {}
    #from pyNastran.utils import object_attributes
    #print('gui attrs', object_attributes(gui))
    #print('gui type = ', type(gui))
    try:
        _create_shell_axes(self, gui, model, grid_id, xyz_cid0)
        _create_alt_axes(self, gui, model, grid_id, xyz_cid0, model.cbar, 'CBAR')
        _create_alt_axes(self, gui, model, grid_id, xyz_cid0, model.cbeam, 'CBEAM')
        _create_alt_axes(self, gui, model, grid_id, xyz_cid0, model.cbush, 'CBUSH')
        _create_alt_axes(self, gui, model, grid_id, xyz_cid0, model.cgap, 'CGAP')
    except IndexError:
        pass

def _create_shell_axes(self: NastranIO,
                       gui: MainWindow,
                       model: BDF,
                       grid_id: np.ndarray,
                       xyz_cid0: np.ndarray) -> None:
    """
    Creates:
     - shell element coordinate systems
     - shell material coordinate systems
    """
    nastran_settings: NastranSettings = gui.settings.nastran_settings
    if not nastran_settings.is_shell_mcids:
        return

    name_func = [
        ('eleement', get_shell_material_coordinate_system),
        ('material', get_shell_element_coordinate_system),
    ]
    for name, func in name_func:
        element_id, length, centroid, ielement, jelement = func(model)
        if len(element_id) == 0:
            return

        imax = ~np.isnan(ielement.max(axis=1))
        jmax = ~np.isnan(jelement.max(axis=1))
        assert len(imax) == len(length)

        p1 = centroid
        p2 = centroid + ielement * length[:, np.newaxis] / 2
        p3 = centroid + jelement * length[:, np.newaxis] / 2
        name_i = f'shell {name}_coord i axis'
        name_j = f'shell {name}_coord j axis'

        scale = 1.0
        gui.create_alternate_vtk_grid(
            name_i, color=RED_FLOAT, line_width=5, opacity=1.,
            point_size=5, representation='bar', bar_scale=scale, is_visible=False)
        gui.create_alternate_vtk_grid(
            name_j, color=GREEN_FLOAT, line_width=5, opacity=1.,
            point_size=5, representation='bar', bar_scale=scale, is_visible=False)

        if imax.sum():
            xyz_lines = np.column_stack([p1, p2])[imax, :]
            _add_nastran_bar_vectors_to_grid(gui, name_i, xyz_lines, element_id[imax])
        else:
            model.log.error('nan shell i')

        if jmax.sum():
            xyz_lines = np.column_stack([p1, p3])[jmax, :]
            _add_nastran_bar_vectors_to_grid(gui, name_j, xyz_lines, element_id[jmax])
        else:
            model.log.error('nan shell j')

def _create_alt_axes(self: NastranIO,
                     gui: MainWindow,
                     model: BDF,
                     grid_id: np.ndarray,
                     xyz_cid0: np.ndarray,
                     elem, card_name: str) -> None:
    """
    Creates orientation vectors for:
     - CBAR
     - CBEAM
     - CBUSH
     - CGAP

    Parameters
    ----------
    self : NastranIO
        do I need this?

    """
    assert gui is not None
    if elem.n == 0:
        return

    xyz1, xyz2 = elem.get_xyz()
    #neids = xyz1.shape[0]
    centroid = (xyz1 + xyz2) / 2.
    assert centroid.shape[0] == elem.nodes.shape[0]
    assert not np.isnan(np.max(xyz1)), xyz1
    assert not np.isnan(np.max(xyz2)), xyz2

    try:
        v, ihat, jhat, khat, wa, wb = elem.get_axes(xyz1, xyz2)
    except ValueError:
        self.log.error(f'failed getting {elem.type} orientation vectors')
        return
    length = elem.length()
    del v, ihat, wa, wb

    if elem.type in {'CBUSH', 'CGAP'} and np.isnan(length.max()):
        try:
            mean_length = np.nanmean(length)
        except:
            ## TODO: should be based on the mean edge length of the model
            mean_length = 1.0
            raise
        inan = np.isnan(length)
        length[inan] = mean_length
        del mean_length, inan

    p1 = centroid
    p2 = centroid + jhat * length[:, np.newaxis] / 2
    p3 = centroid + khat * length[:, np.newaxis] / 2
    #i1 = np.arange(neids)
    #i2 = i1 + neids
    #node_lines = np.column_stack([i1, i2])

    name_nodes = f'{card_name} nodes'
    name_bar_y = f'{card_name} y axis'
    name_bar_z = f'{card_name} z axis'

    scale = 1.0
    nids = np.unique(elem.nodes.ravel())
    inid = np.searchsorted(grid_id, nids)
    all_xyz = xyz_cid0[inid, :]
    _build_dots(gui, name_nodes, all_xyz, color=RED_FLOAT)

    eids = elem.element_id
    xyz_lines = np.column_stack([p1, p2])
    if len(xyz_lines):
        gui.create_alternate_vtk_grid(
            name_bar_y, color=GREEN_FLOAT, line_width=5, opacity=1.,
            point_size=5, representation='bar', bar_scale=scale, is_visible=False)
        _add_nastran_bar_vectors_to_grid(gui, name_bar_y, xyz_lines, eids)

    xyz_lines = np.column_stack([p1, p3])
    if len(xyz_lines):
        gui.create_alternate_vtk_grid(
            name_bar_z, color=BLUE_FLOAT, line_width=5, opacity=1.,
            point_size=5, representation='bar', bar_scale=scale, is_visible=False)
        _add_nastran_bar_vectors_to_grid(gui, name_bar_z, xyz_lines, eids)

    #gui._add_nastran_lines_xyz_to_grid(name_bar_y, xyz_lines, eids)
    #gui.bar_lines[name] = bar_lines
    #eids = elem.element_id

def _add_nastran_bar_vectors_to_grid(gui: MainWindow,
                                     name: str,
                                     lines: np.ndarray,
                                     eids: np.ndarray) -> None:
    """Creates the bar orientation vector lines

    Parameters
    ----------
    name: str
        the name of the actor
    lines : (nlines, 6) float array
        [x1, y1, z1, x2, y2, z2], # line1
        [x1, y1, z1, x2, y2, z2], # line2
    eids : (n,) int array
        informational

    """
    assert name is not None
    nlines = len(lines)
    nnodes = nlines * 2

    assert name != 'Bar Nodes', name
    grid = gui.alt_grids[name]

    bar_eids_array = np.asarray(eids, dtype='int32')
    bar_lines_array = np.asarray(lines, dtype='float32').reshape(nlines, 6)
    #for line in bar_lines_array[:10, :]:
        #print(line)
    #bar_eids = gui.bar_eids
    #bar_lines = gui.bar_lines

    gui.bar_eids[name] = bar_eids_array
    gui.bar_lines[name] = bar_lines_array

    nodes = bar_lines_array.reshape(nlines * 2, 3)
    points = numpy_to_vtk_points(nodes)
    elements = np.arange(0, nnodes, dtype='int32').reshape(nlines, 2)

    etype = 3 # vtkLine().GetCellType()
    create_vtk_cells_of_constant_element_type(grid, elements, etype)
    grid.SetPoints(points)

def create_alt_rbe2_grids(gui: MainWindow,
                          model: BDF,
                          grid_id: np.ndarray,
                          xyz_cid0: np.ndarray):
    """
    Creates global & individual actors for:
      - RBE2 independents
      - RBE2 dependents
      - RBE2 lines

    """
    elem = model.rbe2
    nastran_settings: NastranSettings = gui.settings.nastran_settings
    if elem.n == 0 or not nastran_settings.is_rbe:
        return

    i_independent = np.searchsorted(grid_id, elem.independent_node)
    try:
        xyz_independents = xyz_cid0[i_independent, :]
    except IndexError:
        return

    name = f'RBE2 independents'
    _build_dots(gui, name, xyz_independents, color=BLUE_FLOAT)

    name = f'RBE2 dependents'
    all_dep_nodes = elem.dependent_nodes
    inid_all_dependent = np.searchsorted(grid_id, all_dep_nodes)
    all_xyz_dependents = xyz_cid0[inid_all_dependent, :]
    _build_dots(gui, name, all_xyz_dependents, color=RED_FLOAT)

    _build_rbe2_vtk_lines(gui, elem, xyz_independents, all_xyz_dependents)

    for eid, ind_node, ind_component, idim, xyz_independent in zip(elem.element_id,
                                                                   elem.independent_node,
                                                                   elem.independent_dof,
                                                                   elem.idim, xyz_independents):
        idim0, idim1 = idim
        dep_nodes = elem.dependent_nodes[idim0:idim1]  # Gmi
        inid_dependent = np.searchsorted(grid_id, dep_nodes)
        xyz_dependents = xyz_cid0[inid_dependent, :]

        name = f'RBE2 {eid} independent (dof={ind_component})'
        _build_dot(gui, name, xyz_independent, color=BLUE_FLOAT)

        name = f'RBE2 {eid} dependent'
        _build_dots(gui, name, xyz_dependents, color=RED_FLOAT)


def _build_rbe2_vtk_lines(gui, elem: RBE2,
                          xyz_independents: np.ndarray,
                          xyz_dependents: np.ndarray):
    lines = []
    i0 = len(elem.element_id)
    xyz_lines = np.vstack([xyz_independents,
                           xyz_dependents])
    for i, eid, ind_node, idim, xyz_independent in zip(
            count(), elem.element_id, elem.independent_node, elem.idim, xyz_independents):
        idim0, idim1 = idim
        dep_nodes = elem.dependent_nodes[idim0:idim1]
        ndep_nodes = len(dep_nodes)

        # initialize to the independent node
        nodes = np.full((ndep_nodes, 2), i, dtype='int32')
        n2 = np.arange(i0, i0+ndep_nodes, dtype='int32')
        nodes[:, 1] = n2

        # we leave off the +1 cause we're starting at
        i0 += ndep_nodes
        lines.append(nodes)
    node_lines = np.vstack(lines)
    name = f'RBE2 lines'
    _build_lines(gui, name, xyz_lines, node_lines, color=LIGHT_GREEN_FLOAT)


def create_plotels(gui: MainWindow,
                   model: BDF,
                   grid_id: np.ndarray,
                   xyz_cid0: np.ndarray) -> None:
    """creates the PLOTELs"""
    nastran_settings: NastranSettings = gui.settings.nastran_settings
    if not nastran_settings.is_plotel:
        return
    card_types = [
        (model.plotel, 'line'),
        (model.plotel3, 'tri3'),
        (model.plotel4, 'quad4'),
        (model.plotel6, 'tri6'),
        (model.plotel8, 'quad8'),

        (model.plottet, 'tet4'),
        (model.plothex, 'hexa8'),
        (model.plotpyr, 'pyram5'),
        (model.plotpyr, 'penta6'),
    ]
    etypes = []
    elements_list = []
    for (card, etype) in card_types:
        if len(card) == 0:
            continue
        inodes = np.searchsorted(grid_id, card.nodes)
        elements_list.append(inodes)
        etypes.append(etype)
    if len(etypes) == 0:
        return

    stacked_ids = np.hstack([inode.ravel() for inode in elements_list])
    uids = np.unique(stacked_ids)
    xyzs = xyz_cid0[uids, :]
    color = nastran_settings.plotel_color
    name = 'plotel'
    gui.create_alternate_vtk_grid(
        name, color=color, line_width=2, opacity=0.8,
        representation='point', is_visible=False, is_pickable=False)

    ugrid = gui.alt_grids[name]
    create_vtk_cells_of_constant_element_types(ugrid, elements_list, etypes)
    vtk_points = ugrid.GetPoints()
    vtk_points = numpy_to_vtk_points(xyzs, points=vtk_points)
    ugrid.SetPoints(vtk_points)

def create_monpnt(gui: MainWindow,
                  model: BDF,
                  grid_id: np.ndarray,
                  xyz_cid0: np.ndarray):
    nastran_settings: NastranSettings = gui.settings.nastran_settings
    if not nastran_settings.is_aero:
        return
    create_monpnt1(gui, model, grid_id, xyz_cid0)
    create_caero(gui, model, grid_id, xyz_cid0)

def create_caero(gui: MainWindow,
                 model: BDF,
                 grid_id: np.ndarray,
                 xyz_cid0: np.ndarray):
    caero1 = model.caero1
    caero2 = model.caero2
    n = len(caero1)
    assert len(caero2) == 0, caero2
    if n == 0:
        #print('no caero1')
        return

    log = model.log
    if model.aeros is not None:
        local_coord_id = model.aeros.acsid
    elif model.aero is not None:
        local_coord_id = model.aero.acsid
    else:
        local_coord_id = 0
        log.error('missing AERO/AEROS; assuming AERO/AEROS acid=0')
    #log.info(f'local_coord_id={local_coord_id}')

    fdtype = caero1.p1.dtype
    p1_global = np.zeros((n, 3), dtype=fdtype)
    p4_global = np.zeros((n, 3), dtype=fdtype)
    for cp in np.unique(caero1.cp):
        #log.info(f'CAERO: cp={cp}')
        icp = np.where(cp == caero1.cp)[0]
        p1_global[icp, :] = model.coord.transform_local_xyz_to_global(caero1.p1[icp, :], cp)
        p4_global[icp, :] = model.coord.transform_local_xyz_to_global(caero1.p4[icp, :], cp)

    #62 + 28 = 90
    #if local_coord_id == 0 and 0:
        #p2_global = p1_global.copy()
        #p3_global = p4_global.copy()
        #p2_global[:, 0] += caero1.x12
        #p3_global[:, 0] += caero1.x43
    #else:
    dxyz12 = np.zeros((n, 3), dtype=fdtype)
    dxyz43 = np.zeros((n, 3), dtype=fdtype)
    dxyz12[:, 0] = caero1.x12
    dxyz43[:, 0] = caero1.x43
    coord = model.coord.slice_card_by_id(local_coord_id)
    #print(coord.i)
    #dxyz12_global = coord.transform_local_xyz_to_global(dxyz12, local_coord_id)
    #dxyz43_global = coord.transform_local_xyz_to_global(dxyz43, local_coord_id)
    dxyz12_global = coord.transform_force_local_to_global(dxyz12, local_coord_id)
    dxyz43_global = coord.transform_force_local_to_global(dxyz43, local_coord_id)
    #print('p1 =', caero.p1)
    #print('p1_global =', p1_global)
    #print('dxyz12 =', dxyz12)
    #print('dxyz12_global =', dxyz12_global)
    #print('-------------')
    p2_global = p1_global + dxyz12_global
    p3_global = p4_global + dxyz43_global

    aefact = model.aefact
    points = []
    elements = []
    ipoint = 0
    ipoint_aerobox = 0
    xyz_total = np.zeros((4 * n, 3), dtype='float64')
    elements_total = np.zeros((n, 4), dtype='int32')
    element_id = []
    for ielement, eid, p1, p2, p3, p4, nspan, nchord, lspan, lchord in zip(
        count(), caero1.element_id, p1_global, p2_global, p3_global, p4_global,
        caero1.nspan, caero1.nchord, caero1.lspan, caero1.lchord):
        xyz_total[ipoint, :] = p1
        xyz_total[ipoint+1, :] = p2
        xyz_total[ipoint+2, :] = p3
        xyz_total[ipoint+3, :] = p4
        elements_total[ielement] = [ipoint, ipoint + 1, ipoint + 2, ipoint + 3]
        ipoint += 4

        #p1, p2, p3, p4 = self.get_points()
        #x, y = self.xy
        if nchord == 0:
            lchord_ref = aefact.slice_card_by_id(lchord)
            #log.warning(f'CAERO1: eid={eid} nchord=0')
            #raise RuntimeError(f'CAERO1: eid={eid} nchord=0')
            #continue
            x = lchord_ref.fractions
            nchord = len(x) - 1  # elements, not nodes
        else:
            x = np.linspace(0., 1., nchord + 1)

        if nspan == 0:
            lspan_ref = aefact.slice_card_by_id(lchord)
            #log.warning(f'CAERO1: eid={eid} nspan=0')
            #raise RuntimeError(f'CAERO1: eid={eid} nspan=0')
            #continue
            y = lspan_ref.fractions
            nspan = len(y) - 1  # elements, not nodes
        else:
            y = np.linspace(0., 1., nspan + 1)

        if nchord < 1 or nspan < 1:
            msg = (f'CAERO1 eid={eid} '
                   f'nchord={nchord} nspan={nspan} '
                   f'lchord={lchord} lspan={lspan}')
            raise RuntimeError(msg)

        nelementsi = nchord * nspan
        element_id.append(np.arange(eid, eid+nelementsi))
        # We're reordering the points so we get the node ids and element ids
        # to be consistent with Nastran.  This is only useful if you're plotting
        # aero panel forces
        #
        # this gives us chordwise panels and chordwise nodes
        pointsi, elementsi = points_elements_from_quad_points(
            p1, p4, p3, p2, y, x, dtype='int32')
        points.append(pointsi)
        elements.append(elementsi + ipoint_aerobox)
        npoints = len(pointsi)
        ipoint_aerobox += npoints

    if len(xyz_total) == 0:
        log.warning(f'CAERO1: no panels/aeroboxes')
        return

    name = 'caero'
    _build_quads(gui, name,
                 xyz_total, elements_total,
                 line_width=3, color=YELLOW_FLOAT,
                 is_visible=True, representation='wire')

    if len(points) == 0:
        log.warning(f'CAERO1: no aeroboxes?')
        return
    aero_element_ids = np.hstack(element_id)
    aero_xyz = np.vstack(points)
    aero_elements = np.vstack(elements)
    #print(element_id)
    #print(aero_elements)
    assert len(aero_element_ids) == len(aero_elements), 'missing caero elements...'

    #counts = np.bincount(aero_element_ids)
    assert len(aero_element_ids) == len(np.unique(aero_element_ids)), 'overwrote caero elements...'

    name = 'caero_boxes'
    _build_quads(gui, name,
                 aero_xyz, aero_elements,
                 line_width=5, color=YELLOW_FLOAT,
                 representation='wire',
                 is_visible=False)
    create_aesurf(gui, model,
                  aero_xyz,
                  aero_elements,
                  aero_element_ids)

def create_aesurf(gui: MainWindow,
                  model: BDF,
                  aero_xyz: np.ndarray,
                  aero_elements: np.ndarray,
                  aero_element_ids: np.ndarray):
    """
    Creates the AESURF actors
     - flap_control_surface cid=100'
     - aileron_control_surface cid=200'
     - caero_control_surfaces

    Supports:
     - aelist_id2

    TODO: add zfighting_offset support
    """
    aelist = model.aelist
    aelist_ids = aelist.aelist_id
    all_aero_element_ids = []

    aesurf = model.aesurf
    #print(aesurf.write(size=8, is_double=False))
    for i in range(aesurf.n):
        aesurf_id = aesurf.aesurf_id[i]
        label = aesurf.label[i]
        local_aelist_ids = aesurf.aelist_id[i, :]
        #print(aesurf)
        assert isinstance(aesurf_id, integer_types), (aesurf_id, type(aesurf_id))
        assert isinstance(label, str), label

        jj = []
        for iaelist_id, aelist_id in enumerate(local_aelist_ids):
            assert isinstance(aelist_id, integer_types), (aelist_id, type(aelist_id))
            if aelist_id == 0:
                continue
            #if iaelist_id == 1: # and aelist_id is None:
                #continue
            aelisti = aelist.slice_card_by_aelist_id(aelist_id)
            aelist_idi = aelisti.aelist_id
            assert aelist_idi == aelist_id, f'AESURF label={label!r} aelist_id{local_aelist_ids+1} not found'

            aero_element_idsj = aelisti.elements
            j = np.searchsorted(aero_element_ids, aero_element_idsj)
            if not np.array_equal(aero_element_ids[j], aero_element_idsj):
                msg = (
                    f'aero_element_ids={aero_element_ids}\n'
                    f'aero_element_ids[j]={aero_element_ids[j]}\n'
                    f'aero_element_idsj={aero_element_idsj}')
                raise RuntimeError(msg)
            jj.append(j)

        j = hstack_unique(jj, unique=True)
        #aesurf_id = integer(card, 1, 'aesid')
        #label = string(card, 2, 'label')

        #cid1 = integer(card, 3, 'cid1')
        #alid1 = integer(card, 4, 'alid1')

        #cid2 = integer_or_blank(card, 5, 'cid2')
        #alid2 = integer_or_blank(card, 6, 'alid2')
        #name = f'AESURF label={label!r} cid={cid1}'
        #print(name)
        cid1 = aesurf.coord_id[i, 0]
        cs_name = f'{label}_control_surface cid={cid1}'
        if len(jj) == 2:
            cid2 = aesurf.coord_id[i, 1]
            cs_name += f' cid2={cid2}'
        aero_elementsj = aero_elements[j, :]
        _build_quads(gui, cs_name,
                     aero_xyz, aero_elementsj,
                     line_width=5, color=PINK_FLOAT,
                     is_visible=True, representation='surface',
                     opacity=0.5,
                     )
        all_aero_element_ids.append(aero_element_idsj)

    if len(all_aero_element_ids) >= 2:
        # only add this actor if there are multiple control surfaces
        cs_name = 'caero_control_surfaces'
        all_aero_element_ids = np.unique(np.hstack(all_aero_element_ids))
        k = np.searchsorted(aero_element_ids, all_aero_element_ids)
        aero_elementsk = aero_elements[k, :]

        _build_quads(gui, cs_name,
                     aero_xyz, aero_elementsk,
                     line_width=5, color=PINK_FLOAT,
                     is_visible=True, representation='surface',
                     opacity=0.5,
                     )


def create_monpnt1(gui: MainWindow,
                  model: BDF,
                  grid_id: np.ndarray,
                  xyz_cid0: np.ndarray) -> None:
    """
    Creates MONPNT1 actors
     - points
     - xyz summation

    """
    monpnt = model.monpnt1
    ncards = len(monpnt)
    if ncards == 0:
        return
    xyz_globals = np.empty((ncards, 3), dtype=xyz_cid0.dtype)
    log = model.log

    #print('monpnt1', monpnt)
    ucps = np.unique(monpnt.cp)
    monpnt_xyz = monpnt.xyz
    for local_coord_id in ucps:
        i = np.where(local_coord_id == monpnt.cp)[0]
        xyz_cid = monpnt_xyz[i, :]
        xyz_global = model.coord.transform_local_xyz_to_global(xyz_cid, local_coord_id)
        xyz_globals[i, :] = xyz_global
    del i, xyz_global

    aecomp = model.aecomp
    aecompl = model.aecompl
    all_aecomp_names = aecomp.name.tolist()
    all_aecompl_names = aecompl.name.tolist()

    for name, label, cp, aecomp_name, xyz_global in zip(monpnt.name, monpnt.label, monpnt.cp, monpnt.comp, xyz_globals):
        actor_name = f'MONPNT1: {name} xyz; cp={cp}'
        _build_dot(gui, actor_name, xyz_global, color=BLUE_FLOAT, point_size=5)

        if aecomp_name in all_aecomp_names:
            iaecomp = all_aecomp_names.index(aecomp_name)
            list_type = aecomp.list_type[iaecomp]
            name = str(name)
            label = str(label)
            aecomp_name = str(aecomp_name)

            if list_type == 'SET1':
                _create_monpnt_aecomp_set1(
                    gui, model, iaecomp,
                    grid_id, xyz_cid0,
                    name, label, aecomp_name, xyz_global)
            elif list_type == 'AELIST':
                log.warning(f'skipping MONPNT1: {name!r} label={label!r} '
                            f'because AECOMP={aecomp_name!r} has a list_type=AELIST') # \nkeys={keys}
                #_create_monpnt_aecomp_aelist(
                    #gui, model, iaecomp,
                    #aero_xyz, aero_element_ids, aero_elements,
                    #name, label, aecomp_name, xyz_global)
            elif list_type == 'CAERO':
                log.warning(f'skipping MONPNT1: {name!r} label={label!r} because '
                            f'AECOMP={aecomp_name!r} has a list_type=CAERO') # \nkeys={keys}
            else:  # pragma: no cover
                log.warning(f'skipping MONPNT1: {name!r} label={label!r} because '
                            f'AECOMP={aecomp_name!r} has a list_type={list_type!r}') # \nkeys={keys}
                raise RuntimeError(f'skipping MONPNT1: {name!r} label={label!r} '
                                   f'because AECOMP={aecomp_name!r} has a list_type={list_type!r}') # \nkeys={keys}
        elif aecomp_name in all_aecompl_names:
            aecomp_names = get_aecomp_from_aecompl(
                aecomp_name,
                aecomp, aecompl,
                all_aecomp_names, all_aecompl_names)

            i = np.searchsorted(all_aecomp_names, aecomp_names)
            assert np.array_equal(aecomp.name[i], aecomp_names)
            list_types = np.unique(aecomp.list_type[i])
            assert len(list_types) == 1, list_types
            list_type = list_types[0]
            if list_type == 'SET1':
                #_create_monpnt_aecomp_set1(
                    #gui, model, iaecomp,
                    #grid_id, xyz_cid0,
                    #name, label, aecomp_name, xyz_global)
                log.warning(f'skipping MONPNT1: {name!r} label={label!r} '
                            f'because AECOMPs={aecomp_names} has a list_type=SET1') # \nkeys={keys}
            elif list_type == 'AELIST':
                log.warning(f'skipping MONPNT1: {name!r} label={label!r} '
                            f'because AECOMPs={aecomp_names} has a list_type=AELIST') # \nkeys={keys}
            elif list_type == 'CAERO':
                log.warning(f'skipping MONPNT1: {name!r} label={label!r} '
                            f'because AECOMPs={aecomp_names} has a list_type=CAERO') # \nkeys={keys}
            else:  # pragma: no cover
                raise RuntimeError(f'skipping MONPNT1: {name!r} label={label!r} '
                                   f'because AECOMPs={aecomp_names} has a list_type={list_type!r}') # \nkeys={keys}

        else:  # pragma: no cover
            log.warning(f'skipping MONPNT1: {name!r} label={label!r} '
                        f'because AECOMP/L={aecomp_name!r} does not exist') # \nkeys={keys}
            continue
    return

def get_aecomp_from_aecompl(aecomp_name: str,
                            aecomp: AECOMP,
                            aecompl: AECOMPL,
                            all_aecomp_names: list[str],
                            all_aecompl_names: list[str]) -> list[str]:
    used_aecompl_names = set([aecomp_name])

    iaecompl = all_aecompl_names.index(aecomp_name)
    aecompl_ilabel = aecompl.ilabel
    ilabel0, ilabel1 = aecompl_ilabel[iaecompl]
    labels = aecompl.labels[ilabel0:ilabel1].tolist()
    labels = set(labels)

    #nlabels = len(labels)
    #nlabels_old = nlabels + 1
    aecomp_names = []
    #aecompl_names =
    while len(labels):
        aecompls = []
        for label in labels:
            if label in all_aecomp_names:
                aecomp_names.append(label)
            elif label in all_aecompl_names:
                assert label not in used_aecompl_names, f'AECOMPL label={label!r} was already used; circular reference'
                aecompls.append(label)
            else:  # pragma: no cover
                raise RuntimeError(label)

        all_labels = []
        for name2 in aecompls:
            iaecompl = all_aecompl_names.index(name2)
            ilabel0, ilabel1 = aecompl_ilabel[iaecompl]
            labelsi = aecompl.labels[ilabel0:ilabel1].tolist()
            all_labels.extend(labelsi)
        labels = list(set(all_labels))
    aecomp_names = list(set(aecomp_names))
    return aecomp_names

def _create_monpnt_aecomp_set1(gui: MainWindow,
                               model: BDF, iaecomp: int,
                               all_nids: np.ndarray,
                               xyz_cid0: np.ndarray,
                               name: str, label: str,
                               comp: str, xyz_global: np.ndarray) -> None:
    """
    Creates MONPNT1 actors
     - points
     - xyz summation
    """
    if len(model.set1) == 0:
        model.log.warning(f'skipping MONPNT1={name!r} because AECOMP={comp!r} and no SET1 points exist') # \nkeys={keys}
        return
    nids = get_aecomp_set1_nids_by_index(iaecomp, model.aecomp, model.set1)
    inids = np.searchsorted(all_nids, nids)
    xyzs = xyz_cid0[inids, :]

    #------------------------------------------------------------
    name = f'MONPNT1: {name} GRIDs'
    _build_dots(gui, name, xyzs, color=RED_FLOAT, point_size=5)

def get_aecomp_set1_nids_by_index(iaecomp: int,
                                  aecomp: AECOMP,
                                  set1: SET1) -> np.ndarray:
    ilist0, ilist1 = aecomp.ilist[iaecomp, :]
    set1_ids = aecomp.lists[ilist0:ilist1]
    #set1 = set1.slice_card_by_set_id(set1_ids)
    set1 = set1.slice_card_by_id(set1_ids)
    nids = np.unique(set1.ids)
    return nids

def _create_monpnt_aecomp_aelist(gui: MainWindow,
                                 model: BDF, iaecomp: int,
                                 aero_xyz: np.ndarray,
                                 aero_element_ids: np.ndarray,
                                 aero_elements: np.ndarray,
                                 name: str, label: str,
                                 comp: str, xyz_global: np.ndarray) -> None:
    """
    Creates MONPNT1 actors
     - aelist panels
     - xyz summation

    """
    aecomp = model.aecomp

    ilist0, ilist1 = aecomp.ilist[iaecomp, :]
    aelist_ids = aecomp.lists[ilist0:ilist1]
    aelist = model.aelist.slice_card_by_aelist_id(aelist_ids)

    elements = np.unique(aelist.elements)
    ieids = np.searchsorted(aero_element_ids, elements)
    aero_elementsi = aero_elements[ieids, :]

    #------------------------------------------------------------
    name = f'MONPNT1: {name} AELIST'
    _build_quads(gui, name,
                 aero_xyz, aero_elementsi,
                 line_width=3, color=RED_FLOAT,
                 is_visible=True, representation='wire')

#def _create_monpnt_aecomp_caero(gui: MainWindow,
                                #model: BDF, iaecomp: int,
                                #all_nids: np.ndarray,
                                #xyz_cid0: np.ndarray,
                                #name: str, label: str, cp: int,
                                #comp: str, xyz_global: np.ndarray):
    #aecomp = model.aecomp
    #nids = []

    #ilist0, ilist1 = aecomp.ilist[iaecomp, :]
    #caero_ids = aecomp.lists[ilist0:ilist1]
    #caero = model.caero1.slice_card_by_set_id(caero_ids)

    #inids = np.searchsorted(all_nids, nids)
    #xyzs = xyz_cid0[inids, :]

    ##------------------------------------------------------------
    #name = f'MONPNT1: {name} GRIDs; cid={cp}'
    #_build_dots(gui, name, xyzs, color=RED_FLOAT, point_size=5)

    #name = f'MONPNT1: {name} xyz'
    #_build_dot(gui, name, xyz_global, color=BLUE_FLOAT, point_size=5)


def _build_dot(gui: MainWindow, name: str, xyzi: np.ndarray,
               point_size: int=3, color=RED_FLOAT,
               is_visible: bool=False) -> None:
    """creates a vtkActor with 1 point"""
    j = 0
    gui.create_alternate_vtk_grid(
        name, color=color, point_size=point_size, opacity=1.0,
        representation='point', is_visible=is_visible, is_pickable=False)
    alt_grid = gui.alt_grids[name]

    alt_grid.Allocate(1, 10)
    #sphere_size = self._get_sphere_size(dim_max)

    points = vtkPoints()
    points.SetNumberOfPoints(1)
    points.InsertPoint(0, *xyzi)

    elem = vtkVertex()
    point_ids = elem.GetPointIds()
    point_ids.SetId(0, j)
    alt_grid.InsertNextCell(elem.GetCellType(), point_ids)

    alt_grid.SetPoints(points)


def _build_dots(gui: MainWindow, name: str, xyzs: np.ndarray,
                point_size: int=3, color=RED_FLOAT,
                is_visible: bool=False) -> None:
    """creates a vtkActor with N points"""
    assert len(xyzs.shape) == 2, xyzs.shape
    gui.create_alternate_vtk_grid(
        name, color=color, point_size=point_size, opacity=1.0,
        representation='point', is_visible=is_visible, is_pickable=False)
    alt_grid = gui.alt_grids[name]

    nelement = xyzs.shape[0]
    #ids = np.ones((nnode, 2), dtype='int32')

    nnodes = np.ones((nelement, 1), dtype='int64')
    nodes_index = np.arange(nelement).reshape((nelement, 1)).astype('int64')

    cell_type_point = 1
    _build_vtk_data_from_dnode(
        alt_grid, xyzs,
        nnodes, nodes_index,
        nelement, cell_type_point, dnode=1)
    return

def _build_lines(gui: MainWindow, name: str,
                 xyzs: np.ndarray,
                 nodes_index: np.ndarray,
                 line_width: int=3, color=RED_FLOAT,
                 representation: str = 'wire',
                 is_visible: bool=True) -> None:
    """Creates a vtkActor with N lines"""
    assert len(xyzs.shape) == 2, xyzs.shape
    assert len(nodes_index.shape) == 2, nodes_index.shape
    nelement = nodes_index.shape[0]

    gui.create_alternate_vtk_grid(
        name, color=color, line_width=line_width, opacity=1.0,
        representation=representation, is_visible=is_visible, is_pickable=False)
    alt_grid = gui.alt_grids[name]

    nnodes = np.ones((nelement, 1), dtype='int64') * 2
    #nodes_index = np.arange(nelement).reshape((nelement, 2)).astype('int64')
    #ids = np.ones((nnode, 2), dtype='int32')

    cell_type_line = 3
    _build_vtk_data_from_dnode(
        alt_grid, xyzs,
        nnodes, nodes_index,
        nelement, cell_type_line, dnode=4)
    return

def _build_quads(gui: MainWindow, name: str,
                 xyzs: np.ndarray,
                 nodes_index: np.ndarray,
                 line_width: int=3, color=RED_FLOAT,
                 opacity: float=1.0,
                 representation: str='wire+surf',
                 is_visible: bool=True) -> None:
    """Creates a vtkActor with N quads"""
    assert len(xyzs.shape) == 2, xyzs.shape
    assert len(nodes_index.shape) == 2, nodes_index.shape
    assert nodes_index.shape[1] == 4, nodes_index.shape
    nelement = nodes_index.shape[0]

    gui.create_alternate_vtk_grid(
        name, color=color, line_width=line_width, opacity=opacity,
        representation=representation, is_visible=is_visible, is_pickable=False)
    alt_grid = gui.alt_grids[name]

    nnodes = np.ones((nelement, 1), dtype='int64') * 4
    #nodes_index = np.arange(nelement).reshape((nelement, 2)).astype('int64')
    #ids = np.ones((nnode, 2), dtype='int32')

    cell_type_quad4 = 9
    _build_vtk_data_from_dnode(
        alt_grid, xyzs,
        nnodes, nodes_index,
        nelement, cell_type_quad4, dnode=2)
    return

def _build_vtk_data_from_dnode(alt_grid: vtkUnstructuredGrid,
                               xyz: np.ndarray,
                               nnodes: np.ndarray,
                               nodes_index: np.ndarray,
                               nelement: int, cell_type: int, dnode: int) -> None:
    """
    Creates a vtkUnstructuredGrid with a given cell_type for nelements.
    This helps to avoid the complexity of vtk requiring the number of
    nodes listed for each element and makes vectorization easy.
    """
    points = numpy_to_vtk_points(xyz)

    cell_type_ = np.ones(nelement, dtype='int64') * cell_type
    n_nodes = np.hstack([nnodes, nodes_index]).ravel()

    # (nnodes+1) = 4+1 = 5
    # [0, 5, 10, 15, 20, ... (nelements-1)*5]
    #
    # for 2 CQUAD4s elements (4 nodes; 5 columns including the node count of 4)
    # [0, 5]
    # should be length nelement
    cell_offset = np.arange(0, nelement * (dnode + 1), dnode + 1)

    nelement_total = nelement
    build_vtk_geometry(
        nelement_total, alt_grid,
        n_nodes, cell_type_, cell_offset)
    alt_grid.SetPoints(points)
