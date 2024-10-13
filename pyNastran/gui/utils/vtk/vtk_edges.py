# coding: utf-8
from __future__ import annotations
from typing import TYPE_CHECKING

import vtkmodules.vtkFiltersCore
if hasattr(vtkmodules.vtkFiltersCore, 'vtkExtractEdges'):
    from vtkmodules.vtkFiltersCore import vtkExtractEdges # 9.3.0+
else:
    from vtkmodules.vtkFiltersExtraction import vtkExtractEdges # 9.0.3 - 9.2.6

if TYPE_CHECKING:
    from pyNastran.gui.vtk_rendering_core import (
        # vtkRenderWindow, vtkRenderWindowInteractor,
        #vtkDataSetMapper,
        vtkColorTransferFunction,
        vtkPolyDataMapper)
    from vtkmodules.vtkRenderingLOD import vtkLODActor
    from vtkmodules.vtkRenderingCore import vtkColorTransferFunction
    from pyNastran.gui.vtk_interface import vtkUnstructuredGrid


def create_edges_from_grid(grid_selected: vtkUnstructuredGrid,
                           edge_mapper: vtkPolyDataMapper,
                           edge_actor: vtkLODActor,
                           color_function: vtkColorTransferFunction,
                           is_edges_visible: bool) -> None:
    """helper method to create vtk edges"""
    edges = vtkExtractEdges()
    edges.SetInputData(grid_selected)
    edge_mapper.SetInputConnection(edges.GetOutputPort())

    #is_edges_visible
    edge_actor.SetMapper(edge_mapper)
    edge_actor.GetProperty().SetColor(0., 0., 0.)
    edge_actor.SetPickable(0)

    edge_mapper.SetLookupTable(color_function)
    edge_mapper.SetResolveCoincidentTopologyToPolygonOffset()

    prop = edge_actor.GetProperty()
    prop.SetColor(0., 0., 0.)
    edge_actor.SetVisibility(is_edges_visible)
