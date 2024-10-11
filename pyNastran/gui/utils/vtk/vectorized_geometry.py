
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.gui.vtk_common_core import vtkUnsignedCharArray, vtkPoints, VTK_ID_TYPE
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid, vtkCellArray, vtkVertex
from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk, numpy_to_vtkIdTypeArray
if TYPE_CHECKING:
    from pyNastran.gui.vtk_interface import vtkUnstructuredGrid



def build_vtk_geometry(nelement: int,
                       ugrid: vtkUnstructuredGrid,
                       n_nodes: np.ndarray,
                       cell_type: np.ndarray,
                       cell_offset: np.ndarray) -> None:
    """
    # The cells must be int64 because numpy_to_vtkIdTypeArray requires that.
    # I think it depends on what vtkIdTypeArray() was built with.
    # nnodes_tetra, (nodes_tetra1)
    # nnodes_hexa, (nodes_hexa1)
    cells = np.array([
        4, 0, 1, 2, 3, # tetra
        8, 4, 5, 6, 7, 8, 9, 10, 11 # hex
    ], dtype='int64')

    # The offsets for the cells (i.e., the indices where the cells start)
    # one for each element
    cell_offsets = np.array([0, 5], dtype='int32')

    """
    deep = 1
    #print('cells =', n_nodes)
    #print('cell_type =', cell_type)
    #print('cell_offset =', cell_offset)
    cells = n_nodes
    cell_id_type = numpy_to_vtkIdTypeArray(cells, deep=1)
    vtk_cell = vtkCellArray()
    vtk_cell.SetCells(nelement, cell_id_type)

    # Cell types
    vtk_cell_type = numpy_to_vtk(
        cell_type, deep=deep,
        array_type=vtkUnsignedCharArray().GetDataType())

    vtk_cell_offset = numpy_to_vtk(cell_offset, deep=1,
                                   array_type=VTK_ID_TYPE)

    #ugrid = vtkUnstructuredGrid()
    ugrid.SetCells(vtk_cell_type, vtk_cell_offset, vtk_cell)
    #settings = {}

def create_offset_arrays(all_grid_id: np.ndarray,
                         element_nodes: np.ndarray,
                         nelement: int,
                         cell_type: int,
                         cell_offset0: int,
                         dnode: int) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    nnodesi = np.ones((nelement, 1), dtype='int32') * dnode
    nodes_indexi = np.searchsorted(all_grid_id, element_nodes)
    #n_nodesi = np.hstack([nnodesi, nodes_indexi])
    print(nnodesi.shape, nodes_indexi.shape)
    n_nodesi = np.column_stack([nnodesi, nodes_indexi])

    cell_offseti = _cell_offset(cell_offset0, nelement, dnode)

    #nodes_indexi = np.searchsorted(all_grid_id, element_nodes)
    #nnodesi = np.ones((nelement, 1), dtype='int32') * dnode
    cell_typei = np.ones(nelement, dtype='int32') * cell_type

    cell_offset0 += nelement * (dnode + 1)
    return cell_offset0, n_nodesi, cell_typei, cell_offseti

def _cell_offset(cell_offset0: int, nelement: int, dnode: int) -> np.ndarray:
    r"""
    (nnodes+1) = 4+1 = 5
    [0, 5, 10, 15, 20, ... (nelements-1)*5]

    for 2 CQUAD4s elements (4 nodes; 5 columns including the node count of 4)
    [0, 5]
    should be length nelement
    """
    cell_offseti = cell_offset0 + np.arange(0, nelement*(dnode +1), dnode + 1)
    assert len(cell_offseti) == nelement
    return cell_offseti
