from __future__ import print_function, absolute_import

"""
Quick and easy way to view model.  Not a full blown gui.
"""

from typing import List, Dict
from six import iteritems, itervalues, iterkeys
import numpy as np

from pyNastran.bdf.bdf import BDF


class VTKData(object):
    def __init__(self, ugrid, nid_list, nid_dict, eid_list, eid_dict):
        import vtk
        self.ugrid = ugrid  # type: vtk.vtkUnstructuredGrid
        self.nid_list = nid_list  # type: List[int]
        self.nid_dict = nid_dict  # type: Dict[int, int]
        self.eid_list = eid_list  # type: List[int]
        self.eid_dict = eid_dict  # type: Dict[int, int]

    def visualize(self):
        import vtk

        ugrid = self.ugrid

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(ugrid)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        renderer = vtk.vtkRenderer()
        rw = vtk.vtkRenderWindow()
        rw.AddRenderer(renderer)

        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(rw)
        renderer.AddActor(actor)

        interactor_style = vtk.vtkInteractorStyleTrackballCamera()

        iren.SetInteractorStyle(interactor_style)

        lut = vtk.vtkLookupTable()
        lut.SetHueRange(2 / 3, 0)
        lut.Build()

        mapper.SetLookupTable(lut)
        actor.GetProperty().SetRepresentationToWireframe()

        iren.Initialize()
        rw.Render()
        iren.Start()


_nastran_to_vtk = None


def _basic_elem(elem, cells, nid_dict):
    nids = elem.nodes
    _cells = [nid_dict[nid] for nid in nids]
    cells.append(len(nids))
    cells.extend(_cells)


def _rbe2(elem, cells, nid_dict):
    independent = elem.gn
    dependent = elem.gmi

    nids = []

    for i in range(len(dependent)):
        nids.append(nid_dict[independent])
        nids.append(nid_dict[dependent[i]])

    cells.append(len(nids))
    cells.extend(nids)


def _rbe3(elem, cells, nid_dict):
    independent = elem.refgrid
    dependent = elem.Gijs

    nids = []

    for i in range(len(dependent)):
        nids.append(nid_dict[independent])
        nids.append(nid_dict[dependent[i]])

    cells.append(len(nids))
    cells.extend(nids)


def _get_nastran_to_vtk():
    import vtk

    global _nastran_to_vtk

    if _nastran_to_vtk is not None:
        return _nastran_to_vtk

    _nastran_to_vtk = {
        'CQUAD4': (vtk.VTK_QUAD, _basic_elem),
        'CTRIA3': (vtk.VTK_TRIANGLE, _basic_elem),
        'CBEAM': (vtk.VTK_LINE, _basic_elem),
        'CBUSH': (vtk.VTK_LINE, _basic_elem),
        'CBAR': (vtk.VTK_LINE, _basic_elem),
        'RBE2': (vtk.VTK_POLY_LINE, _rbe2),
        'RBE3': (vtk.VTK_POLY_LINE, _rbe3)
    }

    return _nastran_to_vtk


def to_vtk(bdf):
    # type: (BDF) -> VTKData
    import vtk
    # TODO: use pynastran numpy_to_vtk
    from vtk.util.numpy_support import numpy_to_vtk

    nid_list = []
    nid_dict = {}

    node_pos = np.empty((len(bdf.nodes), 3), dtype=np.float64)

    i = 0
    for node in itervalues(bdf.nodes):
        node_pos[i] = node.get_position()
        nid = node.nid
        nid_list.append(nid)
        nid_dict[nid] = i
        i += 1

    _points = vtk.vtkPoints()
    _points.SetData(numpy_to_vtk(node_pos))

    points = vtk.vtkPoints()
    points.DeepCopy(_points)

    cells = []
    cell_types = []
    cell_count = 0

    eid_list = []
    eid_dict = {}

    _nastran_to_vtk = _get_nastran_to_vtk()

    _element_codes = []

    for elem in itervalues(bdf.elements):
        elem_type = elem.type
        cell_type, add_method = _nastran_to_vtk.get(elem_type, (None, None))

        if cell_type is None:
            continue

        cell_types.append(cell_type)
        _element_codes.append(element_codes[elem_type])

        add_method(elem, cells, nid_dict)

        eid_list.append(elem.eid)
        eid_dict[elem.eid] = cell_count

        cell_count += 1

    cells = np.array(cells, dtype=np.int64)

    id_array = vtk.vtkIdTypeArray()
    id_array.SetVoidArray(cells, len(cells), 1)

    vtk_cells = vtk.vtkCellArray()
    vtk_cells.SetCells(cell_count, id_array)

    cell_types = np.array(cell_types, 'B')
    vtk_cell_types = numpy_to_vtk(cell_types)

    cell_locations = np.array([i for i in range(cell_count)])

    vtk_cell_locations = numpy_to_vtk(cell_locations, deep=1, array_type=vtk.VTK_ID_TYPE)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    ugrid.SetCells(vtk_cell_types, vtk_cell_locations, vtk_cells)

    _element_codes = np.array(_element_codes, dtype=np.int32)
    ecode_array = vtk.vtkIntArray()
    ecode_array.SetVoidArray(_element_codes, len(_element_codes), 1)
    ecode_array.SetName('element_codes')

    ugrid.GetCellData().AddArray(ecode_array)

    copy = vtk.vtkUnstructuredGrid()
    copy.DeepCopy(ugrid)

    vtk_data = VTKData(copy, nid_list, nid_dict, eid_list, eid_dict)

    return vtk_data


element_codes = {
    # grids and other points
    'GRID': -1,
    # various types of mpc's
    'RBAR': -101,
    'RBAR1': -102,
    'RBAX3D': -103,
    'RBE1': -104,
    'RBE2': -105,
    'RBE3': -106,
    'RBE3U': -107,
    'RBJOINT': -108,
    'RCONN': -109,
    # coordinate systems
    'CORD1C': -201,
    'CORD1R': -202,
    'CORD1S': -203,
    'CORD2C': -204,
    'CORD2R': -205,
    'CORD2S': -206,
    'CORD3G': -207,
    'CORD3R': -208,
    # elements, defined by nastran item codes
    'CAXIF2': 47,
    'CAXIF3': 48,
    'CAXIF4': 49,
    'CAXISYM': 241,
    'CBAR': 34,
    'CBEAM': 2,
    'CBEAM3': 184,
    'CBEND': 69,
    'CBUSH': 102,
    'CBUSH1D': 40,
    'CCONEAX': 35,
    'CDUM3': 55,
    'CDUM4': 56,
    'CDUM5': 57,
    'CDUM6': 58,
    'CDUM7': 59,
    'CDUM8': 60,
    'CDUM9': 61,
    'CELAS1': 11,
    'CELAS2': 12,
    'CELAS3': 13,
    'CGAP': 86,
    'CHEXA': 67,
    'CHEXAFD': 202,
    'CIFHEX': 65,
    'CIFPENT': 66,
    'CIFQDX': 73,
    'CIFQUAD': 63,
    'CONROD': 10,
    'CPENTA': 68,
    'CPENTAFD': 204,
    'CQUAD4': 33,
    'CQUAD8': 64,
    'CQUADFD': 201,
    'CQUADR': 82,
    'CQUADX': 18,
    'CQUADXFD': 214,
    'CROD': 1,
    'CSHEAR': 4,
    'CSLOT3': 50,
    'CSLOT4': 51,
    'CTETRA': 39,
    'CTETRAFD': 205,
    'CTRIA3': 74,
    'CTRIA6': 75,
    'CTRIAFD': 206,
    'CTRIAR': 70,
    'CTRIAX': 17,
    'CTRIAX6': 53,
    'CTRIAXFD': 212,
    'CTUBE': 3,
    'CWELDP': 118,
    'CWELDC': 117,
    'CWELD': 200,
    'VUHEXA': 145,
    'VUPENTA': 146,
    'VUTETRA': 147,
}
