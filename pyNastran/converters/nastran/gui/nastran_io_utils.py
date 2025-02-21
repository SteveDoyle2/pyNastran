from __future__ import annotations
import os
import sys
import traceback
from collections import defaultdict
from typing import Any, Optional, cast, TYPE_CHECKING
import numpy as np
from numpy.linalg import norm

from pyNastran.converters.neu.neu import read_neu
from pyNastran.gui.vtk_interface import (
    vtkVertex, vtkLine,
    vtkTriangle, vtkQuad, vtkTetra, vtkWedge, vtkHexahedron,
    vtkQuadraticTriangle, vtkQuadraticQuad, vtkQuadraticTetra,
    vtkQuadraticWedge, vtkQuadraticHexahedron,
    vtkPyramid, vtkQuadraticPyramid,
    vtkQuadraticEdge, vtkBiQuadraticQuad,
    vtkUnstructuredGrid,
)

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.femutils.nan import (
    isfinite, isfinite_and_greater_than, isfinite_and_nonzero,
    isgreater_int)

# bdf
from pyNastran.bdf.bdf import (
    BDF,
    #CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
    CTRSHL,
    CTRAX3, CTRIAX6, CTRIAX, #CTRAX6,
    CQUADX4, CQUADX8, CQUADX,
    #CONM2,
    PCOMP, PCOMPG, PCOMPS, PCOMPLS,
    # nastran95
    CQUAD1,
    SET1, MONPNT1, CORD2R, AECOMP)

from pyNastran.bdf.cards.elements.shell import (
    CQUAD4, CQUAD8, CQUAD, CQUADR, CSHEAR,
    CTRIA3, CTRIA6, CTRIAR,
    CTRIA3, CTRIA6, CTRIAR,
    CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8,
    CPLSTS3, CPLSTS4, CPLSTS6, CPLSTS8,
)
from pyNastran.bdf.cards.elements.solid import (
    CTETRA4, CTETRA10, CPENTA6, CPENTA15,
    CHEXA8, CHEXA20, CIHEX1, CIHEX2, CHEXA1, CHEXA2,
    CPYRAM5, CPYRAM13,
)
from pyNastran.bdf.mesh_utils.delete_bad_elements import (
    tri_quality, quad_quality, get_min_max_theta)
from pyNastran.op2.op2_geom import OP2Geom

#
from .utils import (
    build_offset_normals_dims, # build_map_centroidal_result,
    #get_nastran_gui_layer_word, # check_for_missing_control_surface_boxes,
    get_elements_nelements_unvectorized, get_shell_material_coord,
    #make_nid_map, store_warning,
)

# gui
from pyNastran.gui.utils.vtk.vtk_utils import (
    numpy_to_vtk_points, create_vtk_cells_of_constant_element_type)
from pyNastran.gui.qt_files.colors import (
    RED_FLOAT, BLUE_FLOAT)
from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.gui.gui_objects.displacements import ElementalTableResults # ForceTableResults,

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.gui.gui_objects.settings import Settings, NastranSettings
    from pyNastran.gui.main_window import MainWindow


SIDE_MAP = {}
SIDE_MAP['CHEXA'] = {
    1 : [4, 3, 2, 1],
    2 : [1, 2, 6, 5],
    3 : [2, 3, 7, 6],
    4 : [3, 4, 8, 7],
    5 : [4, 1, 5, 8],
    6 : [5, 6, 7, 8],
}

NO_THETA = [
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    'CBAR', 'CBEAM', 'CBEAM3', 'CBEND',
    'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CVISC',
    'CONROD', 'CROD', 'CTUBE', 'PLOTEL',
    'CHBDYP', 'GENEL',
]

def map_elements1_quality_helper(self,
                                 model: BDF,
                                 xyz_cid0: np.ndarray,
                                 nid_cp_cd: np.ndarray,
                                 nid_map: dict[int, int],
                                 j: int):
    # fancy way to do an inplace "eid_to_nid_map = {}"
    #eid_to_nid_map.clear()
    eid_to_nid_map = {}
    self.eid_to_nid_map = eid_to_nid_map


    eid_map = self.gui.eid_map
    # --------------------------

    nids = nid_cp_cd[:, 0]
    #sphere_size = self._get_sphere_size(dim_max)

    # :param i: the element id in grid
    # :param j: the element id in grid2
    i = 0

    #nids = self.eid_to_nid_map[eid]

    # the list of all pids
    #pids = []

    # pid = pids_dict[eid]
    pids_dict = {}

    elements, nelements, superelements = get_elements_nelements_unvectorized(model)

    pids = np.zeros(nelements, 'int32')
    material_coord = np.full(nelements, -1, dtype='int32')
    material_theta = np.full(nelements, np.nan, dtype='float32')
    min_interior_angle = np.zeros(nelements, 'float32')
    max_interior_angle = np.zeros(nelements, 'float32')
    dideal_theta = np.zeros(nelements, 'float32')
    max_skew_angle = np.zeros(nelements, 'float32')
    max_warp_angle = np.zeros(nelements, 'float32')
    max_aspect_ratio = np.zeros(nelements, 'float32')
    area = np.zeros(nelements, 'float32')
    area_ratio = np.zeros(nelements, 'float32')
    taper_ratio = np.zeros(nelements, 'float32')
    min_edge_length = np.zeros(nelements, 'float32')

    # pids_good = []
    # pids_to_keep = []
    # pids_btm = []
    # pids_to_drop = []

    # 3
    # | \
    # |   \
    # |     \
    # 1------2


    # these normals point inwards
    #      4
    #    / | \
    #   /  |  \
    #  3-------2
    #   \  |   /
    #    \ | /
    #      1
    _ctetra_faces = (
        (0, 1, 2),  # (1, 2, 3),
        (0, 3, 1),  # (1, 4, 2),
        (0, 3, 2),  # (1, 3, 4),
        (1, 3, 2),  # (2, 4, 3),
    )

    # these normals point inwards
    #
    #
    #
    #
    #        /4-----3
    #       /       /
    #      /  5    /
    #    /    \   /
    #   /      \ /
    # 1---------2
    _cpyram_faces = (
        (0, 1, 2, 3),  # (1, 2, 3, 4),
        (1, 4, 2),  # (2, 5, 3),
        (2, 4, 3),  # (3, 5, 4),
        (0, 3, 4),  # (1, 4, 5),
        (0, 4, 1),  # (1, 5, 2),
    )

    # these normals point inwards
    #       /6
    #     /  | \
    #   /    |   \
    # 3\     |     \
    # |  \   /4-----5
    # |    \/       /
    # |   /  \     /
    # |  /    \   /
    # | /      \ /
    # 1---------2
    _cpenta_faces = (
        (0, 2, 1),  # (1, 3, 2),
        (3, 4, 5),  # (4, 5, 6),

        (0, 1, 4, 3),  # (1, 2, 5, 4), # bottom
        (1, 2, 5, 4),  # (2, 3, 6, 5), # right
        (0, 3, 5, 2),  # (1, 4, 6, 3), # left
    )

    # these normals point inwards
    #      8----7
    #     /|   /|
    #    / |  / |
    #   /  5-/--6
    # 4-----3   /
    # |  /  |  /
    # | /   | /
    # 1-----2
    _chexa_faces = (
        (4, 5, 6, 7),  # (5, 6, 7, 8),
        (0, 3, 2, 1),  # (1, 4, 3, 2),
        (1, 2, 6, 5),  # (2, 3, 7, 6),
        (2, 3, 7, 6),  # (3, 4, 8, 7),
        (0, 4, 7, 3),  # (1, 5, 8, 4),
        (0, 6, 5, 4),  # (1, 7, 6, 5),
    )
    nid_to_pid_map = defaultdict(list)
    pid = 0

    log = self.log
    grid = self.gui.grid
    self._build_plotels(model)

    #print("map_elements...")
    for (eid, element) in sorted(elements.items()):
        eid_map[eid] = i
        if i % 5000 == 0 and i > 0:
            print('  map_elements = %i' % i)
        etype = element.type
        # if element.Pid() >= 82:
            # continue
        # if element.Pid() in pids_to_drop:
            # continue
        # if element.Pid() not in pids_to_keep:
            # continue
        # if element.pid.type == 'PSOLID':
            # continue

        pid = np.nan
        dideal_thetai = np.nan
        min_thetai = np.nan
        max_thetai = np.nan
        #max_thetai = np.nan
        max_skew = np.nan
        #max_warp = np.nan
        max_warp = np.nan
        aspect_ratio = np.nan
        areai = np.nan
        area_ratioi = np.nan
        taper_ratioi = np.nan
        min_edge_lengthi = np.nan

        if isinstance(element, (CTRIA3, CTRIAR, CTRAX3, CPLSTN3)):
            if isinstance(element, (CTRIA3, CTRIAR)):
                mcid, theta = get_shell_material_coord(element)
                material_coord[i] = mcid
                material_theta[i] = theta
            elem = vtkTriangle()
            node_ids = element.node_ids
            pid = element.Pid()
            eid_to_nid_map[eid] = node_ids
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)  # or blank?

            n1, n2, n3 = [nid_map[nid] for nid in node_ids]
            p1 = xyz_cid0[n1, :]
            p2 = xyz_cid0[n2, :]
            p3 = xyz_cid0[n3, :]
            out = tri_quality(p1, p2, p3)
            (areai, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out

            point_ids = elem.GetPointIds()
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        elif isinstance(element, (CTRIA6, CPLSTN6, CTRIAX)):
            # the CTRIAX is a standard 6-noded element
            if isinstance(element, CTRIA6):
                mcid, theta = get_shell_material_coord(element)
                material_coord[i] = mcid
                material_theta[i] = theta
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids:
                elem = vtkQuadraticTriangle()
                point_ids = elem.GetPointIds()
                point_ids.SetId(3, nid_map[node_ids[3]])
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkTriangle()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:3]

            n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
            p1 = xyz_cid0[n1, :]
            p2 = xyz_cid0[n2, :]
            p3 = xyz_cid0[n3, :]
            out = tri_quality(p1, p2, p3)
            (areai, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        elif isinstance(element, (CTRIAX6, CTRSHL)):
            # the CTRIAX6 is not a standard second-order triangle
            #
            # 5
            # |\
            # |  \
            # 6    4
            # |     \
            # |       \
            # 1----2----3
            #
            #material_coord[i] = element.theta # TODO: no mcid
            # midside nodes are required, nodes out of order
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            if None not in node_ids:
                elem = vtkQuadraticTriangle()
                point_ids = elem.GetPointIds()
                point_ids.SetId(3, nid_map[node_ids[1]])
                point_ids.SetId(4, nid_map[node_ids[3]])
                point_ids.SetId(5, nid_map[node_ids[5]])
            else:
                elem = vtkTriangle()
                point_ids = elem.GetPointIds()

            n1 = nid_map[node_ids[0]]
            n2 = nid_map[node_ids[2]]
            n3 = nid_map[node_ids[4]]
            p1 = xyz_cid0[n1, :]
            p2 = xyz_cid0[n2, :]
            p3 = xyz_cid0[n3, :]
            out = tri_quality(p1, p2, p3)
            (areai, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            eid_to_nid_map[eid] = [node_ids[0], node_ids[2], node_ids[4]]
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, (CQUAD4, CSHEAR, CQUADR, CPLSTN4, CQUADX4, CQUAD1)):
            if isinstance(element, (CQUAD4, CQUADR, CQUAD1)):
                mcid, theta = get_shell_material_coord(element)
                material_coord[i] = mcid
                material_theta[i] = theta
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)  # or blank?
            eid_to_nid_map[eid] = node_ids

            try:
                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids]
            except KeyError:  # pragma: no cover
                print("node_ids =", node_ids)
                print(str(element))
                #print('nid_map = %s' % nid_map)
                raise
                #continue
            p1 = xyz_cid0[n1, :]
            p2 = xyz_cid0[n2, :]
            p3 = xyz_cid0[n3, :]
            p4 = xyz_cid0[n4, :]
            out = quad_quality(element, p1, p2, p3, p4)
            (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warp) = out

            elem = vtkQuad()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            point_ids.SetId(3, n4)
            grid.InsertNextCell(9, point_ids)

        elif isinstance(element, (CQUAD8, CPLSTN8, CQUADX8)):
            if isinstance(element, CQUAD8):
                mcid, theta = get_shell_material_coord(element)
                material_coord[i] = mcid
                material_theta[i] = theta
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
            p1 = xyz_cid0[n1, :]
            p2 = xyz_cid0[n2, :]
            p3 = xyz_cid0[n3, :]
            p4 = xyz_cid0[n4, :]
            out = quad_quality(element, p1, p2, p3, p4)
            (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warp) = out
            if None not in node_ids:
                elem = vtkQuadraticQuad()
                point_ids = elem.GetPointIds()
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkQuad()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:4]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            point_ids.SetId(3, n4)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, (CQUAD, CQUADX)):
            # CQUAD, CQUADX are 9 noded quads
            mcid, theta = get_shell_material_coord(element)
            material_coord[i] = mcid
            material_theta[i] = theta

            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
            p1 = xyz_cid0[n1, :]
            p2 = xyz_cid0[n2, :]
            p3 = xyz_cid0[n3, :]
            p4 = xyz_cid0[n4, :]
            out = quad_quality(element, p1, p2, p3, p4)
            (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warp) = out
            if None not in node_ids:
                elem = vtkBiQuadraticQuad()
                point_ids = elem.GetPointIds()
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                point_ids.SetId(8, nid_map[node_ids[8]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkQuad()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:4]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            point_ids.SetId(3, n4)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, CTETRA4):
            elem = vtkTetra()
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids[:4]
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            grid.InsertNextCell(10, point_ids)
            #elem_nid_map = {nid:nid_map[nid] for nid in node_ids[:4]}
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _ctetra_faces, node_ids[:4], nid_map, xyz_cid0)

        elif isinstance(element, CTETRA10):
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids:
                elem = vtkQuadraticTetra()
                point_ids = elem.GetPointIds()
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkTetra()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:4]
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _ctetra_faces, node_ids[:4], nid_map, xyz_cid0)

        elif isinstance(element, CPENTA6):
            elem = vtkWedge()
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids[:6]
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            grid.InsertNextCell(13, point_ids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _cpenta_faces, node_ids[:6], nid_map, xyz_cid0)

        elif isinstance(element, CPENTA15):
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids:
                elem = vtkQuadraticWedge()
                point_ids = elem.GetPointIds()
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                point_ids.SetId(10, nid_map[node_ids[10]])
                point_ids.SetId(11, nid_map[node_ids[11]])
                point_ids.SetId(12, nid_map[node_ids[12]])
                point_ids.SetId(13, nid_map[node_ids[13]])
                point_ids.SetId(14, nid_map[node_ids[14]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkWedge()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:6]
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _cpenta_faces, node_ids[:6], nid_map, xyz_cid0)

        elif isinstance(element, (CHEXA8, CIHEX1, CHEXA1)):
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids[:8]
            elem = vtkHexahedron()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            point_ids.SetId(6, nid_map[node_ids[6]])
            point_ids.SetId(7, nid_map[node_ids[7]])
            grid.InsertNextCell(12, point_ids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _chexa_faces, node_ids[:8], nid_map, xyz_cid0)

        elif isinstance(element, (CHEXA20, CIHEX2)):
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids:
                elem = vtkQuadraticHexahedron()
                point_ids = elem.GetPointIds()
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                point_ids.SetId(10, nid_map[node_ids[10]])
                point_ids.SetId(11, nid_map[node_ids[11]])

                # these two blocks are flipped
                point_ids.SetId(12, nid_map[node_ids[16]])
                point_ids.SetId(13, nid_map[node_ids[17]])
                point_ids.SetId(14, nid_map[node_ids[18]])
                point_ids.SetId(15, nid_map[node_ids[19]])

                point_ids.SetId(16, nid_map[node_ids[12]])
                point_ids.SetId(17, nid_map[node_ids[13]])
                point_ids.SetId(18, nid_map[node_ids[14]])
                point_ids.SetId(19, nid_map[node_ids[15]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkHexahedron()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:8]

            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            point_ids.SetId(6, nid_map[node_ids[6]])
            point_ids.SetId(7, nid_map[node_ids[7]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _chexa_faces, node_ids[:8], nid_map, xyz_cid0)

        elif isinstance(element, CPYRAM5):
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids[:5]
            elem = vtkPyramid()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            # etype = 14
            grid.InsertNextCell(elem.GetCellType(), point_ids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _cpyram_faces, node_ids[:5], nid_map, xyz_cid0)
        elif isinstance(element, CPYRAM13):
            node_ids = element.node_ids
            pid = element.Pid()
            if None not in node_ids:
                #print(' node_ids =', node_ids)
                elem = vtkQuadraticPyramid()
                point_ids = elem.GetPointIds()
                # etype = 27
                point_ids.SetId(5, nid_map[node_ids[5]])
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                point_ids.SetId(10, nid_map[node_ids[10]])
                point_ids.SetId(11, nid_map[node_ids[11]])
                point_ids.SetId(12, nid_map[node_ids[12]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkPyramid()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:5]
            #print('*node_ids =', node_ids[:5])


            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)
            min_thetai, max_thetai, dideal_thetai, min_edge_lengthi = get_min_max_theta(
                _cpyram_faces, node_ids[:5], nid_map, xyz_cid0)

        elif etype in ('CBUSH', 'CBUSH1D', 'CFAST',
                       'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                       'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                       'CVISC', 'CGAP'):

            # TODO: verify
            # CBUSH, CBUSH1D, CFAST, CELAS1, CELAS3
            # CDAMP1, CDAMP3, CDAMP4, CDAMP5, CVISC
            if hasattr(element, 'pid'):
                pid = element.pid
            else:
                # CELAS2, CELAS4?
                pid = 0

            node_ids = element.node_ids
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            if node_ids[0] is None and node_ids[1] is None: # CELAS2
                log.warning('removing CELASx eid=%i -> no node %s' % (eid, node_ids[0]))
                del self.eid_map[eid]
                continue
            if None in node_ids:  # used to be 0...
                if node_ids[0] is None:
                    slot = 1
                elif node_ids[1] is None:
                    slot = 0
                #print('node_ids=%s slot=%s' % (str(node_ids), slot))
                eid_to_nid_map[eid] = node_ids[slot]
                nid = node_ids[slot]
                if nid not in nid_map:
                    # SPOINT
                    log.warning('removing CELASx eid=%i -> SPOINT %i' % (eid, nid))
                    continue

                #c = nid_map[nid]

                #if 1:
                elem = vtkVertex()
                point_ids = elem.GetPointIds()
                point_ids.SetId(0, j)
                #else:
                    #elem = vtkSphere()
                    #elem = vtkSphereSource()
                    #if d == 0.:
                    #d = sphere_size
                    #elem.SetRadius(sphere_size)
            else:
                # 2 points
                #d = norm(element.nodes[0].get_position() - element.nodes[1].get_position())
                eid_to_nid_map[eid] = node_ids
                elem = vtkLine()
                point_ids = elem.GetPointIds()
                try:
                    point_ids.SetId(0, nid_map[node_ids[0]])
                    point_ids.SetId(1, nid_map[node_ids[1]])
                except KeyError:
                    print("node_ids =", node_ids)
                    print(str(element))
                    continue

            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif etype in ('CBAR', 'CBEAM', 'CROD', 'CONROD', 'CTUBE'):
            if etype == 'CONROD':
                pid = 0
                areai = element.Area()
            else:
                pid = element.Pid()
                try:
                    areai = element.pid_ref.Area()
                except AttributeError:
                    print(element)
                    areai = -1.

            node_ids = element.node_ids
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)

            # 2 points
            #min_edge_lengthi = norm(element.nodes_ref[0].get_position() -
                                    #element.nodes_ref[1].get_position())
            try:
                n1, n2 = np.searchsorted(nids, element.nodes)
            except Exception:
                print(element.get_stats())
                n1i, n2i = element.nodes
                print('nids =', nids)
                assert n1i in nids, 'n1=%s could not be found' % n1i
                assert n2i in nids, 'n2=%s could not be found' % n2i
                raise
            xyz1 = xyz_cid0[n1, :]
            xyz2 = xyz_cid0[n2, :]
            min_edge_lengthi = norm(xyz2 - xyz1)
            eid_to_nid_map[eid] = node_ids
            elem = vtkLine()
            try:
                n1, n2 = [nid_map[nid] for nid in node_ids]
            except KeyError:  # pragma: no cover
                print("node_ids =", node_ids)
                print(str(element))
                print('nid_map = %s' % nid_map)
                raise
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif etype == 'CBEND':
            pid = element.Pid()
            node_ids = element.node_ids

            # 2 points
            n1, n2 = np.searchsorted(nids, element.nodes)
            xyz1 = xyz_cid0[n1, :]
            xyz2 = xyz_cid0[n2, :]
            #min_edge_lengthi = norm(element.nodes_ref[0].get_position() -
                                    #element.nodes_ref[1].get_position())

            g0 = element.g0 #_vector
            if not isinstance(g0, integer_types):
                log.warning('removing\n%s' % (element))
                log.warning('removing CBEND/g0 eid=%s; %s' % (eid, element.type))
                del self.eid_map[eid]
                continue
                # msg = 'CBEND: g0 must be an integer; g0=%s x=%s\n%s' % (
                #     g0, element.x, element)
                # raise NotImplementedError(msg)

            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids
            # only supports g0 as an integer
            elem = vtkQuadraticEdge()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[g0])
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif etype == 'CHBDYG':
            node_ids = element.node_ids
            pid = 0
            #pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            if element.surface_type in ('AREA4', 'AREA8'):
                eid_to_nid_map[eid] = node_ids[:4]

                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                p3 = xyz_cid0[n3, :]
                p4 = xyz_cid0[n4, :]
                out = quad_quality(element, p1, p2, p3, p4)
                (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warp) = out
                if element.surface_type == 'AREA4' or None in node_ids:
                    elem = vtkQuad()
                    point_ids = elem.GetPointIds()
                else:
                    elem = vtkQuadraticQuad()
                    point_ids = elem.GetPointIds()
                    point_ids.SetId(4, nid_map[node_ids[4]])
                    point_ids.SetId(5, nid_map[node_ids[5]])
                    point_ids.SetId(6, nid_map[node_ids[6]])
                    point_ids.SetId(7, nid_map[node_ids[7]])

                point_ids.SetId(0, n1)
                point_ids.SetId(1, n2)
                point_ids.SetId(2, n3)
                point_ids.SetId(3, n4)
                grid.InsertNextCell(elem.GetCellType(), point_ids)
            elif element.surface_type in ['AREA3', 'AREA6']:
                eid_to_nid_map[eid] = node_ids[:3]
                if element.Type == 'AREA3' or None in node_ids:
                    elem = vtkTriangle()
                    point_ids = elem.GetPointIds()
                else:
                    elem = vtkQuadraticTriangle()
                    point_ids = elem.GetPointIds()
                    point_ids.SetId(3, nid_map[node_ids[3]])
                    point_ids.SetId(4, nid_map[node_ids[4]])
                    point_ids.SetId(5, nid_map[node_ids[5]])

                n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                p3 = xyz_cid0[n3, :]
                out = tri_quality(p1, p2, p3)
                (areai, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out
                point_ids.SetId(0, n1)
                point_ids.SetId(1, n2)
                point_ids.SetId(2, n3)
                grid.InsertNextCell(elem.GetCellType(), point_ids)
            else:
                #print('removing\n%s' % (element))
                log.warning('removing eid=%s; %s' % (eid, element.type))
                del self.eid_map[eid]
                self.gui.log_info("skipping %s" % element.type)
                continue
        elif etype == 'CHBDYP':
            #|    1   |    2    |    3    |   4  |    5   |    6   |  7 |  8 |  9 |
            #| CHBDYP |   EID   |   PID   | TYPE | IVIEWF | IVIEWB | G1 | G2 | G0 |
            #|        | RADMIDF | RADMIDB | GMID |   CE   |   E1   | E2 | E3 |    |
            pid = 0 # element.pid
            node_ids = element.node_ids
            if element.Type == 'LINE':
                n1, n2 = [nid_map[nid] for nid in node_ids[:2]]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                elem = vtkLine()
                point_ids = elem.GetPointIds()
                point_ids.SetId(0, n1)
                point_ids.SetId(1, n2)
            else:
                msg = 'element_solid:\n%s' % (str(element_solid))
                msg += 'mapped_inids = %s\n' % mapped_inids
                msg += 'side_inids = %s\n' % side_inids
                msg += 'nodes = %s\n' % nodes
                #msg += 'side_nodes = %s\n' % side_nodes
                raise NotImplementedError(msg)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif etype == 'CHBDYE':
            #|   1    |  2  |   3  |  4   |   5    |    6   |    7    |    8    |
            #| CHBDYE | EID | EID2 | SIDE | IVIEWF | IVIEWB | RADMIDF | RADMIDB |
            eid_solid = element.eid2
            side = element.side
            element_solid = model.elements[eid_solid]

            try:
                mapped_inids = SIDE_MAP[element_solid.type][side]
            except KeyError:  # pragma: no cover
                log.warning('removing\n%s' % (element))
                log.warning('removing eid=%s; %s' % (eid, element.type))
                del self.eid_map[eid]
                self.gui.log_info("skipping %s" % element.type)
                continue
            side_inids = [nid - 1 for nid in mapped_inids]
            nodes = element_solid.node_ids

            pid = 0
            unused_nnodes = len(side_inids)
            node_ids = [nodes[inid] for inid in side_inids]
            #inids = np.searchsorted(all_nids, node_ids)

            #if len(side_inids) == 2:
                #n1, n2 = [nid_map[nid] for nid in node_ids[:2]]
                #p1 = xyz_cid0[n1, :]
                #p2 = xyz_cid0[n2, :]
                #elem = vtkLine()
                #point_ids = elem.GetPointIds()
                #point_ids.SetId(0, n1)
                #point_ids.SetId(1, n2)
            if len(side_inids) == 3:
                n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                p3 = xyz_cid0[n3, :]
                out = tri_quality(p1, p2, p3)
                (areai, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi) = out

                elem = vtkTriangle()
                point_ids = elem.GetPointIds()
                point_ids.SetId(0, n1)
                point_ids.SetId(1, n2)
                point_ids.SetId(2, n3)
            elif len(side_inids) == 4:
                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
                p1 = xyz_cid0[n1, :]
                p2 = xyz_cid0[n2, :]
                p3 = xyz_cid0[n3, :]
                p4 = xyz_cid0[n4, :]
                out = quad_quality(element, p1, p2, p3, p4)
                (areai, taper_ratioi, area_ratioi, max_skew, aspect_ratio,
                 min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warp) = out

                elem = vtkQuad()
                point_ids = elem.GetPointIds()
                point_ids.SetId(0, n1)
                point_ids.SetId(1, n2)
                point_ids.SetId(2, n3)
                point_ids.SetId(3, n4)
            else:
                msg = 'element_solid:\n%s' % (str(element_solid))
                msg += 'mapped_inids = %s\n' % mapped_inids
                msg += 'side_inids = %s\n' % side_inids
                msg += 'nodes = %s\n' % nodes
                #msg += 'side_nodes = %s\n' % side_nodes
                raise NotImplementedError(msg)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif etype == 'GENEL':
            genel_nids = []
            if len(element.ul_nodes):
                genel_nids.append(element.ul_nodes)
            if len(element.ud_nodes):
                genel_nids.append(element.ud_nodes)
            node_ids = np.unique(np.hstack(genel_nids))
            node_ids = node_ids[:2]
            del genel_nids

            elem = vtkLine()
            try:
                n1, n2 = [nid_map[nid] for nid in node_ids]
            except KeyError:  # pragma: no cover
                print("node_ids =", node_ids)
                print(str(element))
                print('nid_map = %s' % nid_map)
                raise
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

            #areai = np.nan
            pid = 0
            #cell_type = cell_type_line
            #inids = np.searchsorted(all_nids, nids)
            #p1, p2 = xyz_cid0[inids, :]
            #min_edge_lengthi = norm(p2 - p1)
            #nnodes = len(nids)
            #dim = 1
        else:
            log.warning('removing\n%s' % (element))
            log.warning('removing eid=%s; %s' % (eid, element.type))
            del self.eid_map[eid]
            self.gui.log_info("skipping %s" % element.type)
            continue
        # what about MPCs, RBE2s (rigid elements)?
        #   are they plotted as elements?
        #   and thus do they need a property?

        if pid is None:
            # CONROD
            #print(element)
            #pids[i] = 0
            #pids_dict[eid] = 0
            pass
        else:
            pids[i] = pid
            pids_dict[eid] = pid

        if np.isnan(max_thetai) and etype not in NO_THETA:
            print('eid=%s theta=%s...setting to 360. deg' % (eid, max_thetai))
            print(element.rstrip())
            if isinstance(element.nodes[0], integer_types):
                print('  nodes = %s' % element.nodes)
            else:
                for node in element.nodes:
                    print(str(node).rstrip())
            max_thetai = 2 * np.pi
        #print(eid, min_thetai, max_thetai, '\n', element)

        min_interior_angle[i] = min_thetai
        max_interior_angle[i] = max_thetai
        dideal_theta[i] = dideal_thetai
        max_skew_angle[i] = max_skew
        max_warp_angle[i] = max_warp
        max_aspect_ratio[i] = aspect_ratio
        area[i] = areai
        area_ratio[i] = area_ratioi
        taper_ratio[i] = taper_ratioi
        min_edge_length[i] = min_edge_lengthi
        i += 1
    #assert len(self.eid_map) > 0, self.eid_map
    #print('mapped elements')

    nelements = i
    self.gui.nelements = nelements
    #print('nelements=%s pids=%s' % (nelements, list(pids)))
    pids = pids[:nelements]

    out = (
        nid_to_pid_map, xyz_cid0, superelements, pids, nelements,
        material_coord, material_theta,
        area, min_interior_angle, max_interior_angle, max_aspect_ratio,
        max_skew_angle, taper_ratio, dideal_theta,
        area_ratio, min_edge_length, max_warp_angle,
    )
    return out

def map_elements1_no_quality_helper(self,
                                    xyz_cid0: np.ndarray,
                                    nid_cp_cd: np.ndarray,
                                    model: BDF,
                                    nid_map: dict[int, int],
                                    j: int):
    #eid_to_nid_map: dict[int, list[int]],

    nids = nid_cp_cd[:, 0]
    #sphere_size = self._get_sphere_size(dim_max)

    # :param i: the element id in grid
    # :param j: the element id in grid2
    i = 0

    #nids = self.eid_to_nid_map[eid]
    self.eid_to_nid_map = {}

    # the list of all pids
    #pids = []

    # pid = pids_dict[eid]
    pids_dict = {}

    elements, nelements, superelements = get_elements_nelements_unvectorized(model)

    pids = np.zeros(nelements, 'int32')
    material_coord = np.full(nelements, -1, dtype='int32')
    material_theta = np.full(nelements, np.nan, dtype='float32')

    # pids_good = []
    # pids_to_keep = []
    # pids_btm = []
    # pids_to_drop = []

    # 3
    # | \
    # |   \
    # |     \
    # 1------2


    # these normals point inwards
    #      4
    #    / | \
    #   /  |  \
    #  3-------2
    #   \  |   /
    #    \ | /
    #      1
    #_ctetra_faces = (
        #(0, 1, 2), # (1, 2, 3),
        #(0, 3, 1), # (1, 4, 2),
        #(0, 3, 2), # (1, 3, 4),
        #(1, 3, 2), # (2, 4, 3),
    #)

    # these normals point inwards
    #
    #
    #
    #
    #        /4-----3
    #       /       /
    #      /  5    /
    #    /    \   /
    #   /      \ /
    # 1---------2
    #_cpyram_faces = (
        #(0, 1, 2, 3), # (1, 2, 3, 4),
        #(1, 4, 2), # (2, 5, 3),
        #(2, 4, 3), # (3, 5, 4),
        #(0, 3, 4), # (1, 4, 5),
        #(0, 4, 1), # (1, 5, 2),
    #)

    # these normals point inwards
    #       /6
    #     /  | \
    #   /    |   \
    # 3\     |     \
    # |  \   /4-----5
    # |    \/       /
    # |   /  \     /
    # |  /    \   /
    # | /      \ /
    # 1---------2
    #_cpenta_faces = (
        #(0, 2, 1), # (1, 3, 2),
        #(3, 4, 5), # (4, 5, 6),

        #(0, 1, 4, 3), # (1, 2, 5, 4), # bottom
        #(1, 2, 5, 4), # (2, 3, 6, 5), # right
        #(0, 3, 5, 2), # (1, 4, 6, 3), # left
    #)

    # these normals point inwards
    #      8----7
    #     /|   /|
    #    / |  / |
    #   /  5-/--6
    # 4-----3   /
    # |  /  |  /
    # | /   | /
    # 1-----2
    #_chexa_faces = (
        #(4, 5, 6, 7), # (5, 6, 7, 8),
        #(0, 3, 2, 1), # (1, 4, 3, 2),
        #(1, 2, 6, 5), # (2, 3, 7, 6),
        #(2, 3, 7, 6), # (3, 4, 8, 7),
        #(0, 4, 7, 3), # (1, 5, 8, 4),
        #(0, 6, 5, 4), # (1, 7, 6, 5),
    #)
    line_type = 3 # vtkLine().GetCellType()

    nid_to_pid_map = defaultdict(list)
    pid = 0

    log = self.log
    grid = self.gui.grid

    #print("map_elements...")
    eid_to_nid_map = self.eid_to_nid_map
    eid_map = self.gui.eid_map
    for (eid, element) in sorted(elements.items()):
        eid_map[eid] = i
        if i % 5000 == 0 and i > 0:
            print('  map_elements (no quality) = %i' % i)
        etype = element.type
        # if element.Pid() >= 82:
            # continue
        # if element.Pid() in pids_to_drop:
            # continue
        # if element.Pid() not in pids_to_keep:
            # continue
        # if element.pid.type == 'PSOLID':
            # continue

        pid = np.nan

        if isinstance(element, (CTRIA3, CTRIAR, CTRAX3, CPLSTN3, CPLSTS3)):
            if isinstance(element, (CTRIA3, CTRIAR)):
                mcid, theta = get_shell_material_coord(element)
                material_coord[i] = mcid
                material_theta[i] = theta
            elem = vtkTriangle()
            node_ids = element.node_ids
            pid = element.Pid()
            eid_to_nid_map[eid] = node_ids
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            n1, n2, n3 = [nid_map[nid] for nid in node_ids]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]

            point_ids = elem.GetPointIds()
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        elif isinstance(element, (CTRIA6, CPLSTN6, CPLSTS6, CTRIAX)):
            # the CTRIAX is a standard 6-noded element
            if isinstance(element, CTRIA6):
                mcid, theta = get_shell_material_coord(element)
                material_coord[i] = mcid
                material_theta[i] = theta
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids:
                elem = vtkQuadraticTriangle()
                point_ids = elem.GetPointIds()
                point_ids.SetId(3, nid_map[node_ids[3]])
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkTriangle()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:3]

            n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        elif isinstance(element, CTRIAX6):
            # the CTRIAX6 is not a standard second-order triangle
            #
            # 5
            # |\
            # |  \
            # 6    4
            # |     \
            # |       \
            # 1----2----3
            #
            #material_coord[i] = element.theta # TODO: no mcid
            # midside nodes are required, nodes out of order
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            if None not in node_ids:
                elem = vtkQuadraticTriangle()
                point_ids = elem.GetPointIds()
                point_ids.SetId(3, nid_map[node_ids[1]])
                point_ids.SetId(4, nid_map[node_ids[3]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                eid_to_nid_map[eid] = [node_ids[0], node_ids[2], node_ids[4],
                                       node_ids[1], node_ids[3], node_ids[5]]
            else:
                elem = vtkTriangle()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = [node_ids[0], node_ids[2], node_ids[4]]

            n1 = nid_map[node_ids[0]]
            n2 = nid_map[node_ids[2]]
            n3 = nid_map[node_ids[4]]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, CTRSHL):  # nastran95
            # the CTRIAX6 is not a standard second-order triangle
            #
            # 5
            # |\
            # |  \
            # 6    4
            # |     \
            # |       \
            # 1----2----3
            #
            #material_coord[i] = element.theta # TODO: no mcid
            # midside nodes are required, nodes out of order
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids and 0:
                elem = vtkQuadraticTriangle()
                point_ids = elem.GetPointIds()
                point_ids.SetId(3, nid_map[node_ids[1]])
                point_ids.SetId(4, nid_map[node_ids[3]])
                point_ids.SetId(5, nid_map[node_ids[5]])
            else:
                elem = vtkTriangle()
                point_ids = elem.GetPointIds()

            n1 = nid_map[node_ids[0]]
            n2 = nid_map[node_ids[2]]
            n3 = nid_map[node_ids[4]]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            eid_to_nid_map[eid] = [node_ids[0], node_ids[2], node_ids[4]]

            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, (CQUAD4, CSHEAR, CQUADR, CPLSTN4, CPLSTS4, CQUADX4, CQUAD1)):
            if isinstance(element, (CQUAD4, CQUADR, CQUAD1)):
                mcid, theta = get_shell_material_coord(element)
                material_coord[i] = mcid
                material_theta[i] = theta
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids

            try:
                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids]
            except KeyError:  # pragma: no cover
                print("node_ids =", node_ids)
                print(str(element))
                #print('nid_map = %s' % nid_map)
                raise
                #continue
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            #p4 = xyz_cid0[n4, :]

            elem = vtkQuad()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            point_ids.SetId(3, n4)
            grid.InsertNextCell(9, point_ids)

        elif isinstance(element, (CQUAD8, CPLSTN8, CPLSTS8, CQUADX8)):
            if isinstance(element, CQUAD8):
                mcid, theta = get_shell_material_coord(element)
                material_coord[i] = mcid
                material_theta[i] = theta
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            #p4 = xyz_cid0[n4, :]
            if None not in node_ids:
                elem = vtkQuadraticQuad()
                point_ids = elem.GetPointIds()
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkQuad()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:4]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            point_ids.SetId(3, n4)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, (CQUAD, CQUADX)):
            # CQUAD, CQUADX are 9 noded quads
            mcid, theta = get_shell_material_coord(element)
            material_coord[i] = mcid
            material_theta[i] = theta

            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            #p4 = xyz_cid0[n4, :]
            if None not in node_ids:
                elem = vtkBiQuadraticQuad()
                point_ids = elem.GetPointIds()
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                point_ids.SetId(8, nid_map[node_ids[8]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkQuad()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:4]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            point_ids.SetId(3, n4)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, CTETRA4):
            elem = vtkTetra()
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids[:4]
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            grid.InsertNextCell(10, point_ids)
            #elem_nid_map = {nid:nid_map[nid] for nid in node_ids[:4]}

        elif isinstance(element, CTETRA10):
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids:
                elem = vtkQuadraticTetra()
                point_ids = elem.GetPointIds()
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkTetra()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:4]
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, CPENTA6):
            elem = vtkWedge()
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids[:6]
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            grid.InsertNextCell(13, point_ids)

        elif isinstance(element, CPENTA15):
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids:
                elem = vtkQuadraticWedge()
                point_ids = elem.GetPointIds()
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                point_ids.SetId(10, nid_map[node_ids[10]])
                point_ids.SetId(11, nid_map[node_ids[11]])
                point_ids.SetId(12, nid_map[node_ids[12]])
                point_ids.SetId(13, nid_map[node_ids[13]])
                point_ids.SetId(14, nid_map[node_ids[14]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkWedge()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:6]
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, (CHEXA8, CIHEX1)):
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids[:8]
            elem = vtkHexahedron()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            point_ids.SetId(6, nid_map[node_ids[6]])
            point_ids.SetId(7, nid_map[node_ids[7]])
            grid.InsertNextCell(12, point_ids)

        elif isinstance(element, (CHEXA20, CIHEX2)):
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids:
                elem = vtkQuadraticHexahedron()
                point_ids = elem.GetPointIds()
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                point_ids.SetId(10, nid_map[node_ids[10]])
                point_ids.SetId(11, nid_map[node_ids[11]])

                # these two blocks are flipped
                point_ids.SetId(12, nid_map[node_ids[16]])
                point_ids.SetId(13, nid_map[node_ids[17]])
                point_ids.SetId(14, nid_map[node_ids[18]])
                point_ids.SetId(15, nid_map[node_ids[19]])

                point_ids.SetId(16, nid_map[node_ids[12]])
                point_ids.SetId(17, nid_map[node_ids[13]])
                point_ids.SetId(18, nid_map[node_ids[14]])
                point_ids.SetId(19, nid_map[node_ids[15]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkHexahedron()
                eid_to_nid_map[eid] = node_ids[:8]

            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            point_ids.SetId(6, nid_map[node_ids[6]])
            point_ids.SetId(7, nid_map[node_ids[7]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, CPYRAM5):
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids[:5]
            elem = vtkPyramid()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            # etype = 14
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        elif isinstance(element, CPYRAM13):
            node_ids = element.node_ids
            pid = element.Pid()
            if None not in node_ids:
                elem = vtkQuadraticPyramid()
                point_ids = elem.GetPointIds()
                #etype = 27
                _nids = [nid_map[node_ids[i]] for i in range(13)]
                point_ids.SetId(0, _nids[0])
                point_ids.SetId(1, _nids[1])
                point_ids.SetId(2, _nids[2])
                point_ids.SetId(3, _nids[3])
                point_ids.SetId(4, _nids[4])

                point_ids.SetId(5, _nids[5])
                point_ids.SetId(6, _nids[6])
                point_ids.SetId(7, _nids[7])
                point_ids.SetId(8, _nids[8])

                point_ids.SetId(9, _nids[9])
                point_ids.SetId(10, _nids[10])
                point_ids.SetId(11, _nids[11])
                point_ids.SetId(12, _nids[12])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkPyramid()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:5]
                point_ids.SetId(0, nid_map[node_ids[0]])
                point_ids.SetId(1, nid_map[node_ids[1]])
                point_ids.SetId(2, nid_map[node_ids[2]])
                point_ids.SetId(3, nid_map[node_ids[3]])
                point_ids.SetId(4, nid_map[node_ids[4]])
            #print('*node_ids =', node_ids[:5])


            #if min(node_ids) > 0:
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif etype in {'CBUSH', 'CBUSH1D', 'CFAST',
                       'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                       'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                       'CVISC', 'CGAP'}:

            # TODO: verify
            # CBUSH, CBUSH1D, CFAST, CELAS1, CELAS3
            # CDAMP1, CDAMP3, CDAMP4, CDAMP5, CVISC
            if hasattr(element, 'pid'):
                pid = element.pid
            else:
                # CELAS2, CELAS4?
                pid = 0

            node_ids = element.node_ids
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            if node_ids[0] is None and node_ids[1] is None: # CELAS2
                log.warning('removing CELASx eid=%i -> no node %s' % (eid, node_ids[0]))
                del self.eid_map[eid]
                continue
            if None in node_ids:  # used to be 0...
                if node_ids[0] is None:
                    slot = 1
                elif node_ids[1] is None:
                    slot = 0
                #print('node_ids=%s slot=%s' % (str(node_ids), slot))
                eid_to_nid_map[eid] = node_ids[slot]
                nid = node_ids[slot]
                if nid not in nid_map:
                    # SPOINT
                    log.warning('removing CELASx eid=%i -> SPOINT %i' % (eid, nid))
                    continue

                #c = nid_map[nid]

                #if 1:
                #print(str(element))
                elem = vtkVertex()
                point_ids = elem.GetPointIds()
                point_ids.SetId(0, j)
                #else:
                    #elem = vtkSphere()
                    #elem = vtkSphereSource()
                    #if d == 0.:
                    #d = sphere_size
                    #elem.SetRadius(sphere_size)
                grid.InsertNextCell(elem.GetCellType(), point_ids)
            else:
                # 2 points
                #d = norm(element.nodes[0].get_position() - element.nodes[1].get_position())
                eid_to_nid_map[eid] = node_ids
                elem = vtkLine()
                point_ids = elem.GetPointIds()
                try:
                    point_ids.SetId(0, nid_map[node_ids[0]])
                    point_ids.SetId(1, nid_map[node_ids[1]])
                except KeyError:
                    print("node_ids =", node_ids)
                    print(str(element))
                    continue
                grid.InsertNextCell(line_type, point_ids)

        elif etype in ('CBAR', 'CBEAM', 'CROD', 'CONROD', 'CTUBE'):
            if etype == 'CONROD':
                pid = 0
                #areai = element.Area()
            else:
                pid = element.Pid()
                #try:
                    #areai = element.pid_ref.Area()
                #except Exception:
                    #print(element)
                    #raise

            node_ids = element.node_ids
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)

            # 2 points
            n1, n2 = np.searchsorted(nids, element.nodes)
            #xyz1 = xyz_cid0[n1, :]
            #xyz2 = xyz_cid0[n2, :]
            eid_to_nid_map[eid] = node_ids
            elem = vtkLine()
            try:
                n1, n2 = [nid_map[nid] for nid in node_ids]
            except KeyError:  # pragma: no cover
                print("node_ids =", node_ids)
                print(str(element))
                print('nid_map = %s' % nid_map)
                raise
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            grid.InsertNextCell(line_type, point_ids)

        elif etype == 'CBEND':
            pid = element.Pid()
            node_ids = element.node_ids
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)

            # 2 points
            n1, n2 = np.searchsorted(nids, element.nodes)
            #xyz1 = xyz_cid0[n1, :]
            #xyz2 = xyz_cid0[n2, :]
            eid_to_nid_map[eid] = node_ids

            if 0:
                g0 = element.g0 #_vector
                if not isinstance(g0, integer_types):
                    msg = 'CBEND: g0 must be an integer; g0=%s x=%s\n%s' % (
                        g0, element.x, element)
                    raise NotImplementedError(msg)
                # only supports g0 as an integer
                elem = vtkQuadraticEdge()
                point_ids = elem.GetPointIds()
                point_ids.SetId(2, nid_map[g0])
            else:
                elem = vtkLine()
                point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif etype == 'CHBDYG':
            node_ids = element.node_ids
            pid = 0
            #pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            if element.surface_type in ['AREA4', 'AREA8']:
                eid_to_nid_map[eid] = node_ids[:4]

                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
                #p1 = xyz_cid0[n1, :]
                #p2 = xyz_cid0[n2, :]
                #p3 = xyz_cid0[n3, :]
                #p4 = xyz_cid0[n4, :]
                if element.surface_type == 'AREA4' or None in node_ids:
                    elem = vtkQuad()
                    point_ids = elem.GetPointIds()
                else:
                    elem = vtkQuadraticQuad()
                    point_ids = elem.GetPointIds()
                    point_ids.SetId(4, nid_map[node_ids[4]])
                    point_ids.SetId(5, nid_map[node_ids[5]])
                    point_ids.SetId(6, nid_map[node_ids[6]])
                    point_ids.SetId(7, nid_map[node_ids[7]])

                point_ids.SetId(0, n1)
                point_ids.SetId(1, n2)
                point_ids.SetId(2, n3)
                point_ids.SetId(3, n4)
                grid.InsertNextCell(elem.GetCellType(), point_ids)
            elif element.surface_type in ['AREA3', 'AREA6']:
                eid_to_nid_map[eid] = node_ids[:3]
                if element.surface_type == 'AREA3' or None in node_ids:
                    elem = vtkTriangle()
                    point_ids = elem.GetPointIds()
                else:
                    elem = vtkQuadraticTriangle()
                    point_ids = elem.GetPointIds()
                    point_ids.SetId(3, nid_map[node_ids[3]])
                    point_ids.SetId(4, nid_map[node_ids[4]])
                    point_ids.SetId(5, nid_map[node_ids[5]])

                n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
                #p1 = xyz_cid0[n1, :]
                #p2 = xyz_cid0[n2, :]
                #p3 = xyz_cid0[n3, :]
                point_ids.SetId(0, n1)
                point_ids.SetId(1, n2)
                point_ids.SetId(2, n3)
                grid.InsertNextCell(elem.GetCellType(), point_ids)
            else:
                #print('removing\n%s' % (element))
                self.log.warning('removing eid=%s; %s' % (eid, element.type))
                del eid_map[eid]
                self.gui.log_info("skipping %s" % element.type)
                continue
        #elif etype == 'CBYDYP':
        elif etype == 'CHBDYE':
            eid_solid = element.eid2
            side = element.side
            element_solid = model.elements[eid_solid]

            try:
                mapped_inids = SIDE_MAP[element_solid.type][side]
            except KeyError:  # pragma: no cover
                log.warning('removing\n%s' % (element))
                log.warning('removing eid=%s; %s' % (eid, element.type))
                del self.eid_map[eid]
                self.gui.log_info("skipping %s" % element.type)
                continue
            side_inids = [nid - 1 for nid in mapped_inids]
            nodes = element_solid.node_ids

            pid = 0
            unused_nnodes = len(side_inids)
            node_ids = [nodes[inid] for inid in side_inids]
            #inids = np.searchsorted(all_nids, node_ids)

            if len(side_inids) == 3:
                n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
                #p1 = xyz_cid0[n1, :]
                #p2 = xyz_cid0[n2, :]
                #p3 = xyz_cid0[n3, :]

                elem = vtkTriangle()
                point_ids = elem.GetPointIds()
                point_ids.SetId(0, n1)
                point_ids.SetId(1, n2)
                point_ids.SetId(2, n3)
            elif len(side_inids) == 4:
                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
                #p1 = xyz_cid0[n1, :]
                #p2 = xyz_cid0[n2, :]
                #p3 = xyz_cid0[n3, :]
                #p4 = xyz_cid0[n4, :]

                elem = vtkQuad()
                point_ids = elem.GetPointIds()
                point_ids.SetId(0, n1)
                point_ids.SetId(1, n2)
                point_ids.SetId(2, n3)
                point_ids.SetId(3, n4)
            else:
                msg = 'element_solid:\n%s' % (str(element_solid))
                msg += 'mapped_inids = %s\n' % mapped_inids
                msg += 'side_inids = %s\n' % side_inids
                msg += 'nodes = %s\n' % nodes
                #msg += 'side_nodes = %s\n' % side_nodes
                raise NotImplementedError(msg)
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        elif etype == 'GENEL':
            node_ids = element.node_ids
            pid = 0
            elem = vtkLine()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
        elif isinstance(element, CHEXA1):
            node_ids = element.node_ids
            pid = 0
            #mid = element.Mid()
            _set_nid_to_pid_map(nid_to_pid_map, pid, node_ids)
            eid_to_nid_map[eid] = node_ids[:8]
            elem = vtkHexahedron()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            point_ids.SetId(6, nid_map[node_ids[6]])
            point_ids.SetId(7, nid_map[node_ids[7]])
            grid.InsertNextCell(12, point_ids)
        elif isinstance(element, CHEXA2):
            node_ids = element.node_ids
            pid = 0
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids:
                elem = vtkQuadraticHexahedron()
                point_ids = elem.GetPointIds()
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                point_ids.SetId(10, nid_map[node_ids[10]])
                point_ids.SetId(11, nid_map[node_ids[11]])

                # these two blocks are flipped
                point_ids.SetId(12, nid_map[node_ids[16]])
                point_ids.SetId(13, nid_map[node_ids[17]])
                point_ids.SetId(14, nid_map[node_ids[18]])
                point_ids.SetId(15, nid_map[node_ids[19]])

                point_ids.SetId(16, nid_map[node_ids[12]])
                point_ids.SetId(17, nid_map[node_ids[13]])
                point_ids.SetId(18, nid_map[node_ids[14]])
                point_ids.SetId(19, nid_map[node_ids[15]])
                eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkHexahedron()
                point_ids = elem.GetPointIds()
                eid_to_nid_map[eid] = node_ids[:8]

            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            point_ids.SetId(6, nid_map[node_ids[6]])
            point_ids.SetId(7, nid_map[node_ids[7]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        else:
            log.warning('removing\n%s' % (element))
            log.warning('removing eid=%s; %s' % (eid, element.type))
            del self.eid_map[eid]
            self.gui.log_info("skipping %s" % element.type)
            continue
        # what about MPCs, RBE2s (rigid elements)?
        #   are they plotted as elements?
        #   and thus do they need a property?

        if pid is None:
            # CONROD
            #print(element)
            #pids[i] = 0
            #pids_dict[eid] = 0
            pass
        else:
            pids[i] = pid
            pids_dict[eid] = pid

        #print(eid, min_thetai, max_thetai, '\n', element)
        i += 1
    #assert len(self.eid_map) > 0, self.eid_map
    #print('mapped elements')

    nelements = i
    #print('nelements=%s pids=%s' % (nelements, list(pids)))
    pids = pids[:nelements]

    out = (
        nid_to_pid_map, xyz_cid0, superelements, pids, nelements,
        material_coord, material_theta,
    )
    return out

def create_ugrid_from_elements(gui: MainWindow,
                               grid: vtkUnstructuredGrid,
                               elements: dict[int, Any],
                               xyz_cid0: np.ndarray,
                               nid_cp_cd: np.ndarray,
                               #model: BDF,
                               nid_map: dict[int, int],
                               log: SimpleLogger):
    nids = nid_cp_cd[:, 0]
    # :param i: the element id in grid
    # :param j: the element id in grid2
    i = 0

    #nids = self.eid_to_nid_map[eid]
    #self.eid_to_nid_map = {}

    # the list of all pids
    #pids = []

    #nelements = len(elements)
    line_type = 3
    nid_to_pid_map = defaultdict(list)
    pid = 0

    #print("map_elements...")
    #eid_to_nid_map = self.eid_to_nid_map
    #eid_map = self.gui.eid_map
    for (eid, element) in sorted(elements.items()):
        #eid_map[eid] = i
        if i % 5000 == 0 and i > 0:
            print('  map_elements (no quality) = %i' % i)
        etype = element.type
        # if element.Pid() >= 82:
            # continue
        # if element.Pid() in pids_to_drop:
            # continue
        # if element.Pid() not in pids_to_keep:
            # continue
        # if element.pid.type == 'PSOLID':
            # continue

        pid = np.nan

        if isinstance(element, (CTRIA3, CTRIAR, CTRAX3, CPLSTN3, CPLSTS3)):
            elem = vtkTriangle()
            node_ids = element.node_ids
            pid = element.Pid()
            #eid_to_nid_map[eid] = node_ids
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            n1, n2, n3 = [nid_map[nid] for nid in node_ids]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]

            point_ids = elem.GetPointIds()
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        elif isinstance(element, (CTRIA6, CPLSTN6, CPLSTS6, CTRIAX)):
            # the CTRIAX is a standard 6-noded element
            node_ids = element.node_ids
            pid = element.Pid()
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)
            if None not in node_ids:
                elem = vtkQuadraticTriangle()
                point_ids = elem.GetPointIds()
                point_ids.SetId(3, nid_map[node_ids[3]])
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                #eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkTriangle()
                point_ids = elem.GetPointIds()
                #eid_to_nid_map[eid] = node_ids[:3]

            n1, n2, n3 = [nid_map[nid] for nid in node_ids[:3]]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        elif isinstance(element, CTRIAX6):
            # the CTRIAX6 is not a standard second-order triangle
            #
            # 5
            # |\
            # |  \
            # 6    4
            # |     \
            # |       \
            # 1----2----3
            #
            #material_coord[i] = element.theta # TODO: no mcid
            # midside nodes are required, nodes out of order
            node_ids = element.node_ids
            #pid = element.Pid()
            #_set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            if None not in node_ids:
                elem = vtkQuadraticTriangle()
                point_ids = elem.GetPointIds()
                point_ids.SetId(3, nid_map[node_ids[1]])
                point_ids.SetId(4, nid_map[node_ids[3]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                #eid_to_nid_map[eid] = [node_ids[0], node_ids[2], node_ids[4],
                                       #node_ids[1], node_ids[3], node_ids[5]]
            else:
                elem = vtkTriangle()
                point_ids = elem.GetPointIds()
                #eid_to_nid_map[eid] = [node_ids[0], node_ids[2], node_ids[4]]

            n1 = nid_map[node_ids[0]]
            n2 = nid_map[node_ids[2]]
            n3 = nid_map[node_ids[4]]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        elif isinstance(element, (CQUAD4, CSHEAR, CQUADR, CPLSTN4, CPLSTS4, CQUADX4, CQUAD1)):
            node_ids = element.node_ids
            #eid_to_nid_map[eid] = node_ids

            try:
                n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids]
            except KeyError:  # pragma: no cover
                print("node_ids =", node_ids)
                print(str(element))
                #print('nid_map = %s' % nid_map)
                raise
                #continue
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            #p4 = xyz_cid0[n4, :]

            elem = vtkQuad()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            point_ids.SetId(3, n4)
            grid.InsertNextCell(9, point_ids)

        elif isinstance(element, (CQUAD8, CPLSTN8, CPLSTS8, CQUADX8)):
            node_ids = element.node_ids
            #pid = element.Pid()
            #_set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            #p4 = xyz_cid0[n4, :]
            if None not in node_ids:
                elem = vtkQuadraticQuad()
                point_ids = elem.GetPointIds()
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                #eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkQuad()
                point_ids = elem.GetPointIds()
                #eid_to_nid_map[eid] = node_ids[:4]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            point_ids.SetId(3, n4)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, (CQUAD, CQUADX)):
            # CQUAD, CQUADX are 9 noded quads
            node_ids = element.node_ids

            n1, n2, n3, n4 = [nid_map[nid] for nid in node_ids[:4]]
            #p1 = xyz_cid0[n1, :]
            #p2 = xyz_cid0[n2, :]
            #p3 = xyz_cid0[n3, :]
            #p4 = xyz_cid0[n4, :]
            if None not in node_ids:
                elem = vtkBiQuadraticQuad()
                point_ids = elem.GetPointIds()
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                point_ids.SetId(8, nid_map[node_ids[8]])
                #eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkQuad()
                point_ids = elem.GetPointIds()
                #eid_to_nid_map[eid] = node_ids[:4]
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            point_ids.SetId(2, n3)
            point_ids.SetId(3, n4)
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, CTETRA4):
            elem = vtkTetra()
            node_ids = element.node_ids
            #eid_to_nid_map[eid] = node_ids[:4]
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            grid.InsertNextCell(10, point_ids)
            #elem_nid_map = {nid:nid_map[nid] for nid in node_ids[:4]}

        elif isinstance(element, CTETRA10):
            node_ids = element.node_ids
            if None not in node_ids:
                elem = vtkQuadraticTetra()
                point_ids = elem.GetPointIds()
                point_ids.SetId(4, nid_map[node_ids[4]])
                point_ids.SetId(5, nid_map[node_ids[5]])
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                #eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkTetra()
                point_ids = elem.GetPointIds()
                #eid_to_nid_map[eid] = node_ids[:4]
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, CPENTA6):
            elem = vtkWedge()
            node_ids = element.node_ids
            #eid_to_nid_map[eid] = node_ids[:6]
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            grid.InsertNextCell(13, point_ids)

        elif isinstance(element, CPENTA15):
            node_ids = element.node_ids
            if None not in node_ids:
                elem = vtkQuadraticWedge()
                point_ids = elem.GetPointIds()
                point_ids.SetId(6, nid_map[node_ids[6]])
                point_ids.SetId(7, nid_map[node_ids[7]])
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                point_ids.SetId(10, nid_map[node_ids[10]])
                point_ids.SetId(11, nid_map[node_ids[11]])
                point_ids.SetId(12, nid_map[node_ids[12]])
                point_ids.SetId(13, nid_map[node_ids[13]])
                point_ids.SetId(14, nid_map[node_ids[14]])
                #eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkWedge()
                point_ids = elem.GetPointIds()
                #eid_to_nid_map[eid] = node_ids[:6]
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, CHEXA8):
            node_ids = element.node_ids
            #eid_to_nid_map[eid] = node_ids[:8]
            elem = vtkHexahedron()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            point_ids.SetId(6, nid_map[node_ids[6]])
            point_ids.SetId(7, nid_map[node_ids[7]])
            grid.InsertNextCell(12, point_ids)

        elif isinstance(element, CHEXA20):
            node_ids = element.node_ids
            if None not in node_ids:
                elem = vtkQuadraticHexahedron()
                point_ids = elem.GetPointIds()
                point_ids.SetId(8, nid_map[node_ids[8]])
                point_ids.SetId(9, nid_map[node_ids[9]])
                point_ids.SetId(10, nid_map[node_ids[10]])
                point_ids.SetId(11, nid_map[node_ids[11]])

                # these two blocks are flipped
                point_ids.SetId(12, nid_map[node_ids[16]])
                point_ids.SetId(13, nid_map[node_ids[17]])
                point_ids.SetId(14, nid_map[node_ids[18]])
                point_ids.SetId(15, nid_map[node_ids[19]])

                point_ids.SetId(16, nid_map[node_ids[12]])
                point_ids.SetId(17, nid_map[node_ids[13]])
                point_ids.SetId(18, nid_map[node_ids[14]])
                point_ids.SetId(19, nid_map[node_ids[15]])
                #eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkHexahedron()
                #eid_to_nid_map[eid] = node_ids[:8]

            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            point_ids.SetId(5, nid_map[node_ids[5]])
            point_ids.SetId(6, nid_map[node_ids[6]])
            point_ids.SetId(7, nid_map[node_ids[7]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif isinstance(element, CPYRAM5):
            node_ids = element.node_ids
            #eid_to_nid_map[eid] = node_ids[:5]
            elem = vtkPyramid()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            point_ids.SetId(2, nid_map[node_ids[2]])
            point_ids.SetId(3, nid_map[node_ids[3]])
            point_ids.SetId(4, nid_map[node_ids[4]])
            # etype = 14
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        elif isinstance(element, CPYRAM13):
            node_ids = element.node_ids
            if None not in node_ids:
                elem = vtkQuadraticPyramid()
                point_ids = elem.GetPointIds()
                #etype = 27
                _nids = [nid_map[node_ids[i]] for i in range(13)]
                point_ids.SetId(0, _nids[0])
                point_ids.SetId(1, _nids[1])
                point_ids.SetId(2, _nids[2])
                point_ids.SetId(3, _nids[3])
                point_ids.SetId(4, _nids[4])

                point_ids.SetId(5, _nids[5])
                point_ids.SetId(6, _nids[6])
                point_ids.SetId(7, _nids[7])
                point_ids.SetId(8, _nids[8])

                point_ids.SetId(9, _nids[9])
                point_ids.SetId(10, _nids[10])
                point_ids.SetId(11, _nids[11])
                point_ids.SetId(12, _nids[12])
                #eid_to_nid_map[eid] = node_ids
            else:
                elem = vtkPyramid()
                point_ids = elem.GetPointIds()
                #eid_to_nid_map[eid] = node_ids[:5]
                point_ids.SetId(0, nid_map[node_ids[0]])
                point_ids.SetId(1, nid_map[node_ids[1]])
                point_ids.SetId(2, nid_map[node_ids[2]])
                point_ids.SetId(3, nid_map[node_ids[3]])
                point_ids.SetId(4, nid_map[node_ids[4]])
            #print('*node_ids =', node_ids[:5])

            #if min(node_ids) > 0:
            grid.InsertNextCell(elem.GetCellType(), point_ids)

        elif etype in {'CBUSH', 'CBUSH1D', 'CFAST',
                       'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                       'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                       'CVISC', 'CGAP'}:

            node_ids = element.node_ids
            _set_nid_to_pid_map_or_blank(nid_to_pid_map, pid, node_ids)

            if node_ids[0] is None and node_ids[1] is None: # CELAS2
                log.warning('removing CELASx eid=%i -> no node %s' % (eid, node_ids[0]))
                #del self.eid_map[eid]
                continue
            if None in node_ids:  # used to be 0...
                if node_ids[0] is None:
                    slot = 1
                elif node_ids[1] is None:
                    slot = 0
                #print('node_ids=%s slot=%s' % (str(node_ids), slot))
                #eid_to_nid_map[eid] = node_ids[slot]
                nid = node_ids[slot]
                if nid not in nid_map:
                    # SPOINT
                    log.warning('removing CELASx eid=%i -> SPOINT %i' % (eid, nid))
                    continue

                #c = nid_map[nid]

                #if 1:
                #print(str(element))
                elem = vtkVertex()
                point_ids = elem.GetPointIds()
                point_ids.SetId(0, j)
                #else:
                    #elem = vtkSphere()
                    #elem = vtkSphereSource()
                    #if d == 0.:
                    #d = sphere_size
                    #elem.SetRadius(sphere_size)
                grid.InsertNextCell(elem.GetCellType(), point_ids)
            else:
                # 2 points
                #d = norm(element.nodes[0].get_position() - element.nodes[1].get_position())
                #eid_to_nid_map[eid] = node_ids
                elem = vtkLine()
                point_ids = elem.GetPointIds()
                try:
                    point_ids.SetId(0, nid_map[node_ids[0]])
                    point_ids.SetId(1, nid_map[node_ids[1]])
                except KeyError:
                    print("node_ids =", node_ids)
                    print(str(element))
                    continue
                grid.InsertNextCell(line_type, point_ids)

        elif etype in ('CBAR', 'CBEAM', 'CROD', 'CONROD', 'CTUBE'):
            node_ids = element.node_ids

            # 2 points
            n1, n2 = np.searchsorted(nids, element.nodes)
            #xyz1 = xyz_cid0[n1, :]
            #xyz2 = xyz_cid0[n2, :]
            #eid_to_nid_map[eid] = node_ids
            elem = vtkLine()
            try:
                n1, n2 = [nid_map[nid] for nid in node_ids]
            except KeyError:  # pragma: no cover
                print("node_ids =", node_ids)
                print(str(element))
                print('nid_map = %s' % nid_map)
                raise
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, n1)
            point_ids.SetId(1, n2)
            grid.InsertNextCell(line_type, point_ids)

        elif etype == 'CBEND':
            node_ids = element.node_ids

            # 2 points
            n1, n2 = np.searchsorted(nids, element.nodes)
            #xyz1 = xyz_cid0[n1, :]
            #xyz2 = xyz_cid0[n2, :]
            #eid_to_nid_map[eid] = node_ids

            point_ids.SetId(0, nid_map[node_ids[0]])
            point_ids.SetId(1, nid_map[node_ids[1]])
            grid.InsertNextCell(elem.GetCellType(), point_ids)
        else:
            log.warning('removing\n%s' % (element))
            log.warning('removing eid=%s; %s' % (eid, element.type))
            #del self.eid_map[eid]
            gui.log_info("skipping %s" % element.type)
            continue
        # what about MPCs, RBE2s (rigid elements)?
        #   are they plotted as elements?
        #   and thus do they need a property?

        #print(eid, min_thetai, max_thetai, '\n', element)
        i += 1
    #assert len(self.eid_map) > 0, self.eid_map
    #print('mapped elements')

    #nelements = i
    #print('nelements=%s pids=%s' % (nelements, list(pids)))

    out = grid
    #(
        #nid_to_pid_map, xyz_cid0, superelements, pids, nelements,
        #material_coord, material_theta,
    #)
    return out

def create_monpnt(gui: MainWindow,
                  model: BDF,
                  xyz_cid0: np.ndarray,
                  nid_cp_cd: np.ndarray):
    cell_type_point = 1  # vtkVertex().GetCellType()

    all_nids = nid_cp_cd[:, 0]
    log = model.log
    for monpnt in model.monitor_points:
        monpnt = cast(MONPNT1, monpnt)
        coord: CORD2R = model.coords[monpnt.cp]
        xyz_global = coord.transform_node_to_global(monpnt.xyz).reshape(1, 3)
        label = monpnt.label
        try:
            aecomp: AECOMP = model.aecomps[monpnt.comp]
        except KeyError:
            key = monpnt.comp
            keys = list(model.aecomps.keys())
            log.warning(f'skipping:\n{monpnt}\nbecause AECOMP/L={key} does not exist\nkeys={keys}')
            continue

        if aecomp.type == 'AECOMP':
            if aecomp.list_type == 'SET1':
                #set1_ids = aecomp.lists
                nids = []
                for set1_id in aecomp.lists:
                    set1: SET1 = model.sets[set1_id]
                    nids += set1.ids

                nids = np.unique(nids)
                inids = np.searchsorted(all_nids, nids)
                xyz = xyz_cid0[inids, :]

                name = f'MONPNT1: {monpnt.name} GRIDs; cid={monpnt.cp}'
                #------------------------------------------------------------
                gui.create_alternate_vtk_grid(
                    name, color=RED_FLOAT, point_size=5, opacity=1.0,
                    representation='point', is_visible=False, is_pickable=False)
                grid = gui.alt_grids[name]

                npoints = len(xyz)
                points = numpy_to_vtk_points(xyz, points=None, dtype='<f', deep=1)
                grid.SetPoints(points)
                elements = np.arange(npoints).reshape(npoints, 1)
                create_vtk_cells_of_constant_element_type(grid, elements, cell_type_point)
                #------------------------------------------------------------
                name = f'MONPNT1: {monpnt.name} xyz'
                gui.create_alternate_vtk_grid(
                    name, color=BLUE_FLOAT, point_size=5, opacity=1.0,
                    representation='point', is_visible=False, is_pickable=False)
                grid = gui.alt_grids[name]

                points = numpy_to_vtk_points(xyz_global, points=None, dtype='<f', deep=1)
                grid.SetPoints(points)
                elements = np.array([[0]], dtype='int32')
                create_vtk_cells_of_constant_element_type(grid, elements, cell_type_point)
        elif aecomp.type == 'AECOMPL':
            log.warning(f'skipping:\n{monpnt}\nbecause AECOMPL is not supported{aecomp}')
        else:  # pragma: no cover
            raise NotImplementedError(aecomp)


def get_results_to_exclude(nastran_settings: NastranSettings) -> set[str]:
    exclude_results = set([])
    if not nastran_settings.eigenvector:
        exclude_results.add('eigenvectors')
    if not nastran_settings.displacement:
        exclude_results.add('displacements')
    if not nastran_settings.velocity:
        exclude_results.add('velocities')
    if not nastran_settings.acceleration:
        exclude_results.add('accelerations')
    if not nastran_settings.temperature:
        exclude_results.add('temperatures')

    if not nastran_settings.spc_force:
        exclude_results.add('spc_forces')
    if not nastran_settings.mpc_force:
        exclude_results.add('mpc_forces')

    plate_etypes = [
        'ctria3', 'ctria6', 'ctriar',
        'cquad4', 'cquad8', 'cquadr',
    ]
    rod_types = ['crod', 'ctube', 'conrod']
    bar_types = ['cbar']
    beam_types = ['cbeam']
    spring_types = ['celas1', 'celas2', 'celas3', 'celas4']

    #---------------------------------------------------------------------------
    flag = 'force'
    if not nastran_settings.force:
        exclude_results.add(flag)

    if not nastran_settings.spring_force:
        names = {f'force.{name}_{flag}' for name in spring_types}
        exclude_results.add(names)
    #if not nastran_settings.rod_force:
        #names = {f'force.{name}_{flag}' for name in rod_types}
        #exclude_results.add(names)
    if not nastran_settings.bar_force:
        names = {f'force.{name}_{flag}' for name in bar_types}
        exclude_results.add(names)
    if not nastran_settings.beam_force:
        names = {f'force.{name}_{flag}' for name in beam_types}
        exclude_results.add(names)

    if not nastran_settings.plate_force:
        names = {f'force.{name}_plate_{flag}' for name in plate_etypes}
        exclude_results.add(names)

    if not nastran_settings.gap_force:
        exclude_results.add(['force.cgap_force'])
    if not nastran_settings.cbush_force:
        exclude_results.add(['force.cbush_force'])

    #---------------------------------------------------------------------------
    flag = 'stress'
    if not nastran_settings.stress:
        exclude_results.add(flag)
    if not nastran_settings.composite_plate_stress:
        names = {f'{flag}.{name}_composite_plate_{flag}' for name in plate_etypes}
        exclude_results.add(names)
    if not nastran_settings.plate_stress:
        names = {f'{flag}.{name}_plate_{flag}' for name in plate_etypes}
        exclude_results.add(names)

    if not nastran_settings.spring_stress:
        names = {f'{flag}.{name}_{flag}' for name in spring_types}
        exclude_results.add(names)
    if not nastran_settings.rod_stress:
        names = {f'{flag}.{name}_{flag}' for name in rod_types}
        exclude_results.add(names)
    if not nastran_settings.bar_stress:
        names = {f'{flag}.{name}_{flag}' for name in bar_types}
        exclude_results.add(names)
    if not nastran_settings.beam_stress:
        names = {f'{flag}.{name}_{flag}' for name in beam_types}
        exclude_results.add(names)

    #---------------------------------------------------------------------------
    flag = 'strain'
    if not nastran_settings.strain:
        exclude_results.add(flag)
    if not nastran_settings.composite_plate_strain:
        names = {f'{flag}.{name}_composite_plate_{flag}' for name in plate_etypes}
        exclude_results.add(names)
    if not nastran_settings.plate_strain:
        names = {f'{flag}.{name}_plate_{flag}' for name in plate_etypes}
        exclude_results.add(names)

    if not nastran_settings.spring_strain:
        names = {f'{flag}.{name}_{flag}' for name in spring_types}
        exclude_results.add(names)
    if not nastran_settings.rod_strain:
        names = {f'{flag}.{name}_{flag}' for name in rod_types}
        exclude_results.add(names)
    if not nastran_settings.bar_strain:
        names = {f'{flag}.{name}_{flag}' for name in bar_types}
        exclude_results.add(names)
    if not nastran_settings.beam_strain:
        names = {f'{flag}.{name}_{flag}' for name in beam_types}
        exclude_results.add(names)
    #---------------------------------------------------------------------------

    if not nastran_settings.strain_energy:
        exclude_results.add('strain_energy*')
    if not nastran_settings.grid_point_force:
        exclude_results.add('grid_point_forces')
    return exclude_results

def get_pcomp_nplies(properties: dict[int, PCOMP | PCOMPG | PCOMPS | PCOMPLS],
                     property_ids_pcomp: list[int]) -> int:
    """
    layer 0 will be defined as the total, so:
    thickness   -> total thickness
    material_id -> null because it's meaningless
    theta       -> null because it's meaningless

    If we have no property_ids, then it's just 0.

    """
    if len(property_ids_pcomp) == 0:
        return 0

    npliesi = 0
    pcomp_nplies = 0
    for pid in property_ids_pcomp:
        try:
            prop = properties[pid]
        except KeyError:
            continue
        pcomp_nplies = max(pcomp_nplies, prop.nplies + 1)
    npliesi = max(npliesi, pcomp_nplies)
    return npliesi

def build_superelement_model(model: BDF, cid: int=0,
                             fdtype: str='float32'):
    models = {0 : model}
    models.update(model.superelement_models)
    #nmodels = len(models)

    xyz_cid0 = {}
    nid_cp_cd = {}
    icd_transform = {}
    #nid_map = {}
    #inode = 0

    for superelement_tuple, modeli in models.items():
        if isinstance(superelement_tuple, int):
            super_id = superelement_tuple
        else:
            super_id = superelement_tuple[1]
        out = modeli.get_displacement_index_xyz_cp_cd(
            fdtype=fdtype, idtype='int32', sort_ids=True)
        icd_transformi, icp_transformi, xyz_cpi, nid_cp_cdi = out
        icd_transform[super_id] = icd_transformi

        xyz_cid0i = modeli.transform_xyzcp_to_xyz_cid(
            xyz_cpi, nid_cp_cdi[:, 0], icp_transformi, cid=cid,
            in_place=False)

        if super_id in model.seloc and super_id: # in model.initial_superelement_models and 0:
            # TODO: when should seloc get applied?
            #       during superelement creation or now?
            #       I'm going with superelement creation...
            #       I think we need to update the node locations for the superelements
            #       that exist before mirroring
            seloc = model.seloc[super_id]
            xyz_cid0i = seloc.transform(model, xyz_cid0i)

        #print('model.spoints =', model.spoints)
        #import json
        #for spoint_id, spoint in model.spoints.items():
            #if spoint.comment: # or spoint._comment?
                #print('SPOINT comment=%r _comment=%r' % (spoint.comment, spoint._comment))
                #comment_lower = spoint.comment.lower()
                #print('comment_lower = %r' % comment_lower)
                ## pyNastran: SPOINT={'id':10, 'xyz':[10.,10.,10.]}
                #if 'pynastran' in comment_lower and 'spoint' in comment_lower:
                    #dict_str = jsonify(comment_lower)
                    #print('dict_str = %r' % dict_str)
                    #dicti = json.loads(dict_str)
                    #print(dicti)
        #for epoint_id, epoint in model.epoints.items():
            #if epoints.comment:
                #print('EPOINT comment=%r _comment=%r' % (spoint.comment, spoint._comment))
        #sys.stdout.flush()

        #------------------------------
        nid_cp_cd[super_id] = nid_cp_cdi
        xyz_cid0[super_id] = xyz_cid0i
    return xyz_cid0, nid_cp_cd, icd_transform

def build_normals_quality(settings: Settings,
                          model: BDF, eid_map, nelements: int, cases, form0, icase: int,
                          xyz_cid0: np.ndarray,
                          material_coord, material_theta,
                          min_interior_angle, max_interior_angle, dideal_theta,
                          area, max_skew_angle, taper_ratio,
                          max_warp_angle, area_ratio, min_edge_length, max_aspect_ratio,
                          make_offset_normals_dim: bool=True,
                          make_xyz: bool=False, make_nnodes_result: bool=False,
                          is_testing: bool=False) -> tuple[int, Any]:
    """
    Creates some nastran specific results

    creates:
     - ElementDim
     - Normal X/Y/Z
     - NNodes/Elem
     - Area
     - Min/Max Interior Angle
     - Skew Angle
     - Taper Ratio
     - Area Ratio
     - MaterialCoord
     - MaterialTheta
    """
    log = model.log
    nastran_settings: NastranSettings = settings.nastran_settings
    colormap = settings.colormap
    #ielement = 0
    #nelements = self.element_ids.shape[0]

    normals = None
    offset = None
    xoffset = None
    yoffset = None
    zoffset = None
    element_dim = None
    nnodes_array = None
    if make_offset_normals_dim:
        try:
            out = build_offset_normals_dims(model, eid_map, nelements)
            normals, offset, xoffset, yoffset, zoffset, element_dim, nnodes_array = out
        except KeyError as error:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
            log.error(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
            log.error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)
            log.error(str(error))
            make_offset_normals_dim = False

    # if not a flat plate
    #if min(nxs) == max(nxs) and min(nxs) != 0.0:
    #is_element_dim = element_dim is not None and np.max(element_dim) != np.min(element_dim)
    is_element_dim = element_dim is not None
    if is_element_dim and isfinite_and_greater_than(element_dim, -1.0):
        eid_dim_res = GuiResult(0, header='ElementDim', title='ElementDim',
                                location='centroid', scalar=element_dim, mask_value=-1)
        cases[icase] = (eid_dim_res, (0, 'ElementDim'))

    #is_shell = normals is not None and np.abs(normals).max() > 0.  # NaN -> 2.0
    is_shell = normals is not None and isfinite(normals)  # using NaNs

    # we have to add the 2nd/3rd lines to make sure bars are getting into this check
    is_solid = (
        isfinite_and_nonzero(min_interior_angle) and
        isfinite_and_nonzero(max_interior_angle)
    )

    #print('is_shell=%s is_solid=%s' % (is_shell, is_solid))
    is_element_quality = area is not None and nastran_settings.is_element_quality
    if is_shell:
        if make_offset_normals_dim:
            nx_res = GuiResult(
                0, header='NormalX', title='NormalX',
                location='centroid', scalar=normals[:, 0], data_format='%.2f')
            ny_res = GuiResult(
                0, header='NormalY', title='NormalY',
                location='centroid', scalar=normals[:, 1], data_format='%.2f')
            nz_res = GuiResult(
                0, header='NormalZ', title='NormalZ',
                location='centroid', scalar=normals[:, 2], data_format='%.2f')
            nxyz_res = NormalResult(0, 'Normals', 'Normals',
                                    nlabels=2, labelsize=5, ncolors=2,
                                    colormap=colormap, data_format='%.1f',
                                    uname='NormalResult')

        if is_element_quality:
            area_res = GuiResult(0, header='Area', title='Area',
                                 location='centroid', scalar=area)
            min_edge_length_res = GuiResult(
                0, header='Min Edge Length', title='Min Edge Length',
                location='centroid', scalar=min_edge_length)

            min_theta_res = GuiResult(
                0, header='Min Interior Angle', title='Min Interior Angle',
                location='centroid', scalar=np.degrees(min_interior_angle))
            max_theta_res = GuiResult(
                0, header='Max Interior Angle', title='Max Interior Angle',
                location='centroid', scalar=np.degrees(max_interior_angle))
            dideal_theta_res = GuiResult(
                0, header='Delta Ideal Angle', title='Delta Ideal Angle',
                location='centroid', scalar=np.degrees(dideal_theta))

            skew = np.degrees(max_skew_angle)
            skew_res = GuiResult(
                0, header='Max Skew Angle', title='MaxSkewAngle',
                location='centroid', scalar=skew)
            aspect_res = GuiResult(
                0, header='Aspect Ratio', title='AspectRatio',
                location='centroid', scalar=max_aspect_ratio)

        form_checks = []
        form0.append(('Element Checks', None, form_checks))
        if is_element_dim:
            form_checks.append(('ElementDim', icase, []))

        if make_offset_normals_dim and make_nnodes_result:
            nnodes_res = GuiResult(
                0, header='NNodes/Elem', title='NNodes/Elem',
                location='centroid', scalar=nnodes_array)
            form_checks.append(('NNodes', icase + 1, []))
            cases[icase + 1] = (nnodes_res, (0, 'NNodes'))
            icase += 1

        if make_offset_normals_dim:
            # 0 is element_dim
            cases[icase + 1] = (nx_res, (0, 'NormalX'))
            cases[icase + 2] = (ny_res, (0, 'NormalY'))
            cases[icase + 3] = (nz_res, (0, 'NormalZ'))
            cases[icase + 4] = (nxyz_res, (0, 'Normal'))

            form_checks.append(('NormalX', icase + 1, []))
            form_checks.append(('NormalY', icase + 2, []))
            form_checks.append(('NormalZ', icase + 3, []))
            form_checks.append(('Normal', icase + 4, []))
            icase += 5

        if is_element_quality:
            cases[icase] = (area_res, (0, 'Area'))
            cases[icase + 1] = (min_edge_length_res, (0, 'Min Edge Length'))
            cases[icase + 2] = (min_theta_res, (0, 'Min Interior Angle'))
            cases[icase + 3] = (max_theta_res, (0, 'Max Interior Angle'))
            cases[icase + 4] = (dideal_theta_res, (0, 'Delta Ideal Angle'))
            cases[icase + 5] = (skew_res, (0, 'Max Skew Angle'))
            cases[icase + 6] = (aspect_res, (0, 'Aspect Ratio'))

            form_checks.append(('Area', icase, []))
            form_checks.append(('Min Edge Length', icase + 1, []))
            form_checks.append(('Min Interior Angle', icase + 2, []))
            form_checks.append(('Max Interior Angle', icase + 3, []))
            form_checks.append(('Delta Ideal Angle', icase + 4, []))
            form_checks.append(('Max Skew Angle', icase + 5, []))
            form_checks.append(('Aspect Ratio', icase + 6, []))
            icase += 7

            if np.any(np.isfinite(area_ratio)) and np.nanmax(area_ratio) > 1.:
                arearatio_res = GuiResult(
                    0, header='Area Ratio', title='Area Ratio',
                    location='centroid', scalar=area_ratio)
                cases[icase] = (arearatio_res, (0, 'Area Ratio'))
                form_checks.append(('Area Ratio', icase, []))
                icase += 1

            if np.any(np.isfinite(taper_ratio)) and np.nanmax(taper_ratio) > 1.:
                taperratio_res = GuiResult(
                    0, header='Taper Ratio', title='Taper Ratio',
                    location='centroid', scalar=taper_ratio)
                cases[icase] = (taperratio_res, (0, 'Taper Ratio'))
                form_checks.append(('Taper Ratio', icase, []))
                icase += 1

            if isfinite_and_nonzero(max_warp_angle):
                warp_res = GuiResult(
                    0, header='Max Warp Angle', title='MaxWarpAngle',
                    location='centroid', scalar=np.degrees(max_warp_angle))
                cases[icase] = (warp_res, (0, 'Max Warp Angle'))
                form_checks.append(('Max Warp Angle', icase, []))
                icase += 1

            #if (np.abs(xoffset).max() > 0.0 or np.abs(yoffset).max() > 0.0 or
                #np.abs(zoffset).max() > 0.0):
            #if isfinite(max_warp_angle):

            # offsets
            if make_offset_normals_dim and np.any(np.isfinite(xoffset)):
                offset_res = GuiResult(
                    0, header='Offset', title='Offset',
                    location='centroid', scalar=offset, data_format='%g',
                    scale_type='length')
                offset_x_res = GuiResult(
                    0, header='OffsetX', title='OffsetX',
                    location='centroid', scalar=xoffset, data_format='%g',
                    scale_type='length')
                offset_y_res = GuiResult(
                    0, header='OffsetY', title='OffsetY',
                    location='centroid', scalar=yoffset, data_format='%g',
                    scale_type='length')
                offset_z_res = GuiResult(
                    0, header='OffsetZ', title='OffsetZ',
                    location='centroid', scalar=zoffset, data_format='%g',
                    scale_type='length')

                cases[icase] = (offset_res, (0, 'Offset'))
                cases[icase + 1] = (offset_x_res, (0, 'OffsetX'))
                cases[icase + 2] = (offset_y_res, (0, 'OffsetY'))
                cases[icase + 3] = (offset_z_res, (0, 'OffsetZ'))

                form_checks.append(('Offset', icase, []))
                form_checks.append(('OffsetX', icase + 1, []))
                form_checks.append(('OffsetY', icase + 2, []))
                form_checks.append(('OffsetZ', icase + 3, []))
                icase += 4

        if 0:  # pragma: no cover
            xyz_offset = np.vstack([xoffset, yoffset, zoffset]).T
            titles = ['Offset XYZ']
            headers = titles
            assert xyz_offset.shape[1] == 3, xyz_offset.shape
            assert xyz_offset.shape[0] == len(offset)
            scales = [1.0]
            subcase_id = 0
            #methods = ['magnitude', 'x', 'y', 'z']
            offset_xyz_res = ElementalTableResults(
                subcase_id, titles, headers, xyz_offset, offset, scales,
                #methods,
            )
            offset_xyz_res.save_defaults()
            cases[icase] = (offset_z_res, (0, 'OffsetZ'))
            form_checks.append(('OffsetXYZ', icase, []))
            icase += 1

        if make_xyz or is_testing:
            x_res = GuiResult(
                0, header='X', title='X',
                location='node', scalar=xyz_cid0[:, 0], data_format='%g',
                scale_type='length')
            y_res = GuiResult(
                0, header='Y', title='Y',
                location='node', scalar=xyz_cid0[:, 1], data_format='%g',
                scale_type='length')
            z_res = GuiResult(
                0, header='Z', title='Z',
                location='node', scalar=xyz_cid0[:, 2], data_format='%g',
                scale_type='length')
            cases[icase] = (x_res, (0, 'X'))
            cases[icase + 1] = (y_res, (0, 'Y'))
            cases[icase + 2] = (z_res, (0, 'Z'))
            form_checks.append(('X', icase + 0, []))
            form_checks.append(('Y', icase + 1, []))
            form_checks.append(('Z', icase + 2, []))
            icase += 3

    elif is_solid:
        # only solid elements
        form_checks = []
        form0.append(('Element Checks', None, form_checks))

        if is_element_dim:
            form_checks.append(('ElementDim', icase, []))
            icase += 1

        if is_element_quality:
            min_edge_length_res = GuiResult(
                0, header='Min Edge Length', title='Min Edge Length',
                location='centroid', scalar=min_edge_length)
            min_theta_res = GuiResult(
                0, header='Min Interior Angle', title='Min Interior Angle',
                location='centroid', scalar=np.degrees(min_interior_angle))
            max_theta_res = GuiResult(
                0, header='Max Interior Angle', title='Max Interior Angle',
                location='centroid', scalar=np.degrees(max_interior_angle))
            #skew = 90. - np.degrees(max_skew_angle)
            #skew_res = GuiResult(0, header='Max Skew Angle', title='MaxSkewAngle',
                                    #location='centroid', scalar=skew)

            form_checks.append(('Min Edge Length', icase, []))
            form_checks.append(('Min Interior Angle', icase + 1, []))
            form_checks.append(('Max Interior Angle', icase + 2, []))
            #form_checks.append(('Max Skew Angle', icase + 3, []))
            cases[icase] = (min_edge_length_res, (0, 'Min Edge Length'))
            cases[icase + 1] = (min_theta_res, (0, 'Min Interior Angle'))
            cases[icase + 2] = (max_theta_res, (0, 'Max Interior Angle'))
            #cases[icase + 3] = (skew_res, (0, 'Max Skew Angle'))
            icase += 3

    else:
        form0.append(('ElementDim', icase, []))
        icase += 1

    if isgreater_int(material_coord, -1):
        material_coord_res = GuiResult(
            0, header='MaterialCoord', title='MaterialCoord',
            location='centroid',
            scalar=material_coord, mask_value=-1, data_format='%i')
        cases[icase] = (material_coord_res, (0, 'MaterialCoord'))
        form0.append(('MaterialCoord', icase, []))
        icase += 1
    if isfinite(material_theta):
        material_theta_res = GuiResult(
            0, header='MaterialTheta', title='MaterialTheta',
            location='centroid',
            scalar=material_theta, data_format='%.3f')
        cases[icase] = (material_theta_res, (0, 'MaterialTheta'))
        form0.append(('MaterialTheta', icase, []))
        icase += 1
    return icase, normals

def _set_nid_to_pid_map(nid_to_pid_map: dict[int, list[int]],
                        pid: int,
                        node_ids: list[int]) -> None:
    for nid in node_ids:
        nid_to_pid_map[nid].append(pid)

def _set_nid_to_pid_map_or_blank(nid_to_pid_map: dict[int, list[int]],
                            pid: int,
                            node_ids: list[Optional[int]]) -> None:
    for nid in node_ids:
        if nid is not None:
            nid_to_pid_map[nid].append(pid)

def get_caero_control_surface_grid(grid: vtkUnstructuredGrid,
                                   box_id_to_caero_element_map: dict[int, int],
                                   caero_points: np.ndarray,
                                   boxes_to_show: list[int],
                                   log):
    j = 0
    areas = []
    centroids = []
    all_points = []
    plot_elements = []
    assert isinstance(boxes_to_show, list), type(boxes_to_show)

    vtk_type = 9 # vtkQuad
    for box_id in boxes_to_show:
        elementi = box_id_to_caero_element_map[box_id]
        pointsi = caero_points[elementi]
        p1, p2, p3, p4 = pointsi
        area = np.linalg.norm(np.cross(p3 - p1, p4 - p2)) / 2.
        if area == 0.0:
            log.warning(f'box_id={box_id:d} has 0 area')
            continue
        #centroid = pointsi.sum(axis=0) / 4.
        centroid = (p1 + p2 + p3 + p4) / 4.
        #assert len(centroid) == 3, centroid

        elem = vtkQuad()
        point_ids = elem.GetPointIds()
        point_ids.SetId(0, j)
        point_ids.SetId(1, j + 1)
        point_ids.SetId(2, j + 2)
        point_ids.SetId(3, j + 3)
        grid.InsertNextCell(vtk_type, point_ids)
        plot_elements.append(j + elementi)
        all_points.append(pointsi)
        centroids.append(centroid)
        areas.append(area)
        j += 4
    elements = np.asarray(plot_elements, dtype='int32')
    return all_points, elements, centroids, areas

def get_model_unvectorized(log: SimpleLogger,
                           bdf_filename: str | BDF,
                           xref_loads: bool=True, is_h5py: bool=True):
    """Loads the BDF/OP2 geometry"""
    ext = '.bdf'
    if isinstance(bdf_filename, str):
        ext = os.path.splitext(bdf_filename)[1].lower()
    elif isinstance(bdf_filename, BDF):
        model = bdf_filename
        xref_nodes = True
        return model, xref_nodes

    punch = None
    if ext == '.pch':
        punch = True

    if ext == '.neu':
        neu_model = read_neu(bdf_filename, log=log)
        model = neu_model.model
        #xref_nodes = True
        #return model, xref_nodes
    elif ext == '.op2':
        model = OP2Geom(make_geom=True, debug=False, log=log,
                        debug_file=None)
        model.clear_results()
        model.IS_TESTING = False
        model.read_op2(op2_filename=bdf_filename)
    elif ext == '.h5' and is_h5py:
        model = BDF(log=log, debug=True)
        model.load_hdf5_filename(bdf_filename)
        model.validate()
    elif ext == '.obj':
        model = BDF(log=log, debug=True)
        model.load(obj_filename=bdf_filename)
    else:  # read the bdf/punch
        model = BDF(log=log, debug=True)

        skip_cards = ['DMI']
        model.disable_cards(skip_cards)

        model.is_strict_card_parser = False
        #model.set_error_storage(nparse_errors=0,
        #                        stop_on_parsing_error=True,
        #                        nxref_errors=0,
        #                        stop_on_xref_error=True)
        model.read_bdf(bdf_filename,
                       punch=punch, xref=False,
                       validate=True)
        #print('done with read_bdf')
        #xref_loads = False
    #xref_aero = len(model.caeros) > 0

    xref_nodes = True
    #model.cross_reference()
    model.safe_cross_reference(
        xref=True,
        xref_nodes=xref_nodes,
        xref_elements=True,
        xref_nodes_with_elements=False,
        xref_properties=True,
        xref_masses=True,
        xref_materials=False,
        xref_loads=xref_loads,
        xref_constraints=False,
        xref_optimization=False,
        xref_aero=True,
        xref_sets=False,
        create_superelement_geometry=True,
    )
    return model, xref_nodes
