import numpy as np

from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points
from pyNastran.gui.vtk_common_core import vtkPoints
from pyNastran.gui.vtk_interface import (
    vtkUnstructuredGrid, vtkVertex, vtkLine, vtkTriangle, vtkQuad,
)

#from vtkmodules.vtkCommonDataModel import vtkCellData, vtkPointData

def add_user_geometry(alt_grid: vtkUnstructuredGrid,
                      geom_grid: vtkUnstructuredGrid,
                      xyz: np.ndarray,
                      nid_map: dict[int, int],
                      nnodes: int,
                      bars: np.ndarray,
                      tris: np.ndarray,
                      quads: np.ndarray,
                      nelements: int, nbars: int,
                      ntris: int, nquads: int) -> vtkPoints:
    """helper method for ``_add_user_geometry``"""
    # set points
    points = numpy_to_vtk_points(xyz, dtype='<f')

    if nelements > 0:
        for i in range(nnodes):
            elem = vtkVertex()
            elem.GetPointIds().SetId(0, i)
            alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            geom_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
    else:
        for i in range(nnodes):
            elem = vtkVertex()
            elem.GetPointIds().SetId(0, i)
            alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

    if nbars:
        for i, bar in enumerate(bars[:, 1:]):
            g1 = nid_map[bar[0]]
            g2 = nid_map[bar[1]]
            elem = vtkLine()
            elem.GetPointIds().SetId(0, g1)
            elem.GetPointIds().SetId(1, g2)
            geom_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

    if ntris:
        for i, tri in enumerate(tris[:, 1:]):
            g1 = nid_map[tri[0]]
            g2 = nid_map[tri[1]]
            g3 = nid_map[tri[2]]
            elem = vtkTriangle()
            elem.GetPointIds().SetId(0, g1)
            elem.GetPointIds().SetId(1, g2)
            elem.GetPointIds().SetId(2, g3)
            geom_grid.InsertNextCell(5, elem.GetPointIds())

    if nquads:
        for i, quad in enumerate(quads[:, 1:]):
            g1 = nid_map[quad[0]]
            g2 = nid_map[quad[1]]
            g3 = nid_map[quad[2]]
            g4 = nid_map[quad[3]]
            elem = vtkQuad()
            point_ids = elem.GetPointIds()
            point_ids.SetId(0, g1)
            point_ids.SetId(1, g2)
            point_ids.SetId(2, g3)
            point_ids.SetId(3, g4)
            geom_grid.InsertNextCell(9, elem.GetPointIds())

    alt_grid.SetPoints(points)
    if nelements > 0:
        geom_grid.SetPoints(points)
    return points


