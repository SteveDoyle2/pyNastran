"""Defines the GUI IO file for Fluent."""
from __future__ import annotations
import os
from typing import Any, TYPE_CHECKING
import numpy as np

try:
    import tables
except ImportError:
    print('pytables was not found; no h5 support.  Run ">>> pip install tables"\n'
          'Do you have h5py installed?  That can cause conflicts.')
    raise

try:
    from tables import open_file, Group, Node, File
except ImportError:
    print('pytables was not found; no h5 support.  Run ">>> pip install tables"\n'
          'Do you have h5py installed?  That can cause conflicts.')
    print(f'tables.__path__ = {tables.__path__}')
    raise

import vtkmodules

from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult
from pyNastran.gui.utils.vtk.vtk_utils import (
    create_vtk_cells_of_constant_element_types, numpy_to_vtk_points)
from pyNastran.gui.utils.vtk.vectorized_geometry import (
    create_offset_arrays)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow
    from vtk import vtkUnstructuredGrid
    from pyNastran.gui.gui_objects.settings import (
        OtherSettings)


class VtkIO:
    def __init__(self, gui: MainWindow):
        self.gui = gui

    def get_vtk_wildcard_geometry_results_functions(self):
        data = ('Vtk CFD',
                'VTK CFD (*.h5)', self.load_vtk_geometry,
                None, None)
        return data

    def load_vtk_geometry(self, h5_filename: str,
                          name: str='main',
                          plot: bool=True):
        gui = self.gui
        model_name = name
        skip_reading = gui._remove_old_geometry(h5_filename)
        if skip_reading:
            return

        log = gui.log
        other_settings: OtherSettings = gui.settings.other_settings

        assert os.path.exists(h5_filename), h5_filename
        with open_file(h5_filename, mode="r", title="", root_uep="/", filters=None) as h5_file:
            nodes, tris, nodal_results, cell_results, centroid, normal = read_vtk_geometry(h5_file)

        node_id = np.unique(tris.ravel())
        assert len(nodes) > 0
        assert len(tris) > 0

        #regions_to_remove = [] #3
        #regions_to_include = [7, 4]
        #regions_to_include = []
        # regions_to_remove = other_settings.cart3d_fluent_remove
        # regions_to_include = other_settings.cart3d_fluent_include

        # out = model.get_area_centroid_normal(tris, quads)
        # (tri_area, quad_area,
        #  tri_centroid, quad_centroid,
        #  tri_normal, quad_normal) = out
        # normal = np.vstack([quad_normal, tri_normal])
        # centroid = np.vstack([quad_centroid, tri_centroid])
        # area = np.hstack([quad_area, tri_area])

        # nnodes_array = np.hstack([quad_area*0+4, tri_area*0+3])
        # del quad_centroid, tri_centroid, quad_area, tri_area

        assert nodes is not None
        nodes = gui.scale_length(nodes)
        # units_pressure_in = other_settings.units_model_in[-1]
        # units_pressure_out = other_settings.units_pressure

        nelement = len(tris)

        nnodes = len(nodes)
        gui.nnodes = nnodes
        gui.nelements = nelement

        gui.log.info(f'nnodes={gui.nnodes:d} nelements={gui.nelements:d}')
        ugrid = gui.grid
        ugrid.Allocate(gui.nelements, 1000)

        points = numpy_to_vtk_points(nodes)
        ugrid.SetPoints(points)
        log.info('created vtk points')

        xmax, ymax, zmax = nodes.max(axis=0)
        xmin, ymin, zmin = nodes.min(axis=0)
        gui.log_info("xmin=%s xmax=%s dx=%s" % (xmin, xmax, xmax-xmin))
        gui.log_info("ymin=%s ymax=%s dy=%s" % (ymin, ymax, ymax-ymin))
        gui.log_info("zmin=%s zmax=%s dz=%s" % (zmin, zmax, zmax-zmin))
        dim_max = max(xmax-xmin, ymax-ymin, zmax-zmin)

        quads = np.zeros((0, 4))
        _create_elements(ugrid, node_id, tris, quads)
        log.info('created vtk elements')

        gui.nid_map = {}
        gui.create_global_axes(dim_max)

        ugrid.Modified()

        # loadSTLResults - regions/loads
        gui.scalar_bar_actor.VisibilityOff()
        gui.scalar_bar_actor.Modified()

        cases = {}
        gui.isubcase_name_map = {}
        idi = 1

        element_id = np.arange(1, nelement+1)
        gui.isubcase_name_map[idi] = ('Vtk', '')
        form, cases = _fill_case(
            cases, idi, node_id, element_id,
            nodal_results, cell_results, normal)

        gui.node_ids = node_id
        gui.element_ids = element_id
        #log.debug(f'running _finish_results_io2')
        gui._finish_results_io2(model_name, form, cases)
        #log.info(f'finished')

def _create_elements(ugrid: vtkUnstructuredGrid,
                     node_id: np.ndarray,
                     tris: np.ndarray,
                     quads: np.ndarray) -> None:
    cell_type_list = []
    cell_offset_list = []

    nquad = len(quads)
    ntri = len(tris)
    cell_offset0 = 0
    n_nodes_list = []

    assert tris.shape[1] == 3, tris.shape
    assert quads.shape[1] == 4, quads.shape
    tri_nodes = tris
    quad_nodes = quads
    assert tri_nodes.shape[1] == 3, tri_nodes.shape
    assert quad_nodes.shape[1] == 4, quad_nodes.shape
    assert node_id.min() >= 0, node_id.min()
    all_nodes = np.unique(np.hstack([tri_nodes.ravel(), quad_nodes.ravel()]))
    assert all_nodes.min() >= 0, all_nodes.min()
    missing_nodes = np.setdiff1d(all_nodes, node_id)
    assert len(missing_nodes) == 0, missing_nodes

    elements_list = []
    etypes = []
    if nquad:
        #elem.GetCellType() = 9  # vtkQuad
        iquad_nodes = np.searchsorted(node_id, quad_nodes)
        elements_list.append(iquad_nodes)
        etypes.append('quad4')
        cell_type = 9
        dnode = 4
        cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
            node_id, quad_nodes,
            nquad, cell_type, cell_offset0, dnode)
        n_nodes_list.append(n_nodesi.ravel())
        cell_type_list.append(cell_typei)
        cell_offset_list.append(cell_offseti)

    if ntri:
        itri_nodes = np.searchsorted(node_id, tri_nodes)
        etypes.append('tri3')
        elements_list.append(itri_nodes)
        #elem.GetCellType() = 5  # vtkTriangle
        cell_type = 5
        dnode = 3
        cell_offset0, n_nodesi, cell_typei, cell_offseti = create_offset_arrays(
            node_id, tri_nodes,
            ntri, cell_type, cell_offset0, dnode)
        n_nodes_list.append(n_nodesi.ravel())
        cell_type_list.append(cell_typei)
        cell_offset_list.append(cell_offseti)

    create_vtk_cells_of_constant_element_types(ugrid, elements_list, etypes)
    return


def _fill_case(cases: dict[int, Any],
               idi: int,
               node_id: np.ndarray,
               element_id: np.ndarray,
               nodal_results: dict[str, np.ndarray],
               cell_results: dict[str, np.ndarray],
               normal: np.ndarray,
               ) -> tuple[list[tuple[str, int, list]],
                          dict[int, Any]]:
    """adds the sidebar results"""
    title = ''
    # reorg the ids
    #element_ids = np.unique(np.hstack([tris[:, 0], quads[:, 0]]))
    #colormap = 'jet'
    icase = 0
    itime = 0
    colormap = 'jet'

    nnode = len(node_id)
    nelement = len(element_id)

    nid_res = GuiResult(idi, header='NodeID', title='NodeID',
                        location='node', scalar=node_id)
    eid_res = GuiResult(idi, header='ElementID', title='ElementID',
                        location='centroid', scalar=element_id)

    nxyz_res = NormalResult(0, 'Normals', 'Normals',
                            nlabels=2, labelsize=5, ncolors=2,
                            #colormap=colormap,
                            data_format='%.1f',
                            uname='NormalResult')
    nx_res = GuiResult(idi, 'normal_x', 'NormalX', 'centroid', normal[:, 0],
                       data_format='%.3f', colormap=colormap, uname='NormalX')
    ny_res = GuiResult(idi, 'normal_y', 'NormalY', 'centroid', normal[:, 1],
                       data_format='%.3f', colormap=colormap, uname='NormalY')
    nz_res = GuiResult(idi, 'normal_z', 'NormalZ', 'centroid', normal[:, 2],
                       data_format='%.3f', colormap=colormap, uname='NormalZ')

    cases[icase] = (nid_res, (itime, 'NodeID'))
    cases[icase + 1] = (eid_res, (itime, 'ElementID'))
    cases[icase + 2] = (nxyz_res, (itime, 'Normal'))
    cases[icase + 3] = (nx_res, (itime, 'NormalX'))
    cases[icase + 4] = (ny_res, (itime, 'NormalY'))
    cases[icase + 5] = (nz_res, (itime, 'NormalZ'))
    form = [
        ('NodeID', icase, []),
        ('ElementID', icase + 1, []),
        ('Normal', icase + 2, []),
        ('NormalX', icase + 3, []),
        ('NormalY', icase + 4, []),
        ('NormalZ', icase + 5, []),
        # ('Results', icase + 5, []),
    ]
    icase += 6

    for name, cell_result in cell_results.items():
        if name == 'RegionId':
            cell_result = cell_result.astype('int32')
        assert len(cell_result) == nelement
        cell_res = GuiResult(idi, header=name, title=name,
                             location='centroid', scalar=cell_result)
        cases[icase] = (cell_res, (itime, title))
        form.append((name, icase, []))
        icase += 1

    for name, node_result in nodal_results.items():
        assert len(node_result) == nnode
        if name == 'vtkOriginalPointIds':
            continue

        if name == 'WallShearStress':
            assert node_result.shape == (nnode, 3), (name, node_result.shape)
            for i, comp in enumerate(['x', 'y', 'z']):
                name_comp = f'{name}_{comp}'
                node_resulti = node_result[:, i]
                node_res = GuiResult(idi, header=name_comp, title=name_comp,
                                     location='node', scalar=node_resulti)
                cases[icase] = (node_res, (itime, title))
                form.append((name, icase, []))
                icase += 1
        else:
            assert node_result.shape == (nnode,), (name, node_result.shape)
            node_res = GuiResult(idi, header=name, title=name,
                                 location='node', scalar=node_result)
            cases[icase] = (node_res, (itime, title))
            form.append((name, icase, []))
            icase += 1
    return form, cases

def read_vtk_geometry(h5_file: File) -> tuple:
    cell_results = {}
    nodal_results = {}

    xyz = h5_file.get_node('/points').read()
    tris = h5_file.get_node('/faces').read()

    cell_data = h5_file.get_node('/cell_data')
    for h5_element in cell_data._f_iter_nodes():
        name = h5_element.name
        data = h5_element.read()
        cell_results[name] = data

    point_data = h5_file.get_node('/point_data')
    for h5_element in point_data._f_iter_nodes():
        name = h5_element.name
        data = h5_element.read()
        nodal_results[name] = data

    itri = tris
    ntri = len(tris)
    nelement = ntri
    assert nelement > 0
    if ntri > 0:
        n1 = itri[:, 0]
        n2 = itri[:, 1]
        n3 = itri[:, 2]
        xyz1 = xyz[n1, :]
        xyz2 = xyz[n2, :]
        xyz3 = xyz[n3, :]
        centroid = (xyz1 + xyz2 + xyz3) / 3
        assert xyz1.shape == xyz2.shape
        assert xyz1.shape == xyz3.shape

        normal = np.cross(xyz2 - xyz1, xyz3 - xyz1, axis=1)
        assert normal.shape == (nelement, 3), normal.shape
        normi = np.linalg.norm(normal, axis=1)
        assert len(normi) == len(xyz1)
        normal /= normi[:, np.newaxis]
    else:
        centroid = None
        normal = None
    return xyz, tris, nodal_results, cell_results, centroid, normal
