"""
defines GuiAttributes, which defines Gui getter/setter methods
and is inherited from many GUI classes
"""
from __future__ import annotations
from pyNastran.gui.qt_files.colors import RED_FLOAT
#import os
#import sys
#import traceback
from typing import Optional, TYPE_CHECKING

import numpy as np
from vtkmodules.vtkRenderingCore import vtkActor, vtkDataSetMapper

from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
#from pyNastran.gui.gui_objects.settings import Settings


#from pyNastran.gui.utils.vtk.gui_utils import remove_actors_from_gui
from pyNastran.gui.utils.vtk.vtk_utils import (
    numpy_to_vtk_points, create_vtk_cells_of_constant_element_type)
from pyNastran.bdf.cards.aero.utils import points_elements_from_quad_points

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.main_window import MainWindow
    from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
    from pyNastran.gui.typing import ColorFloat

class VtkActorActions:
    def __init__(self, gui: MainWindow):
        self.gui = gui


    def set_line_grid(self, name: str,
                      nodes: np.ndarray,
                      elements: np.ndarray,
                      color: ColorFloat,
                      point_size: int=5, line_width: int=5,
                      opacity: float=1.,
                      representation: str='wire',
                      add: bool=True) -> Optional[vtkUnstructuredGrid]:
        """Makes a line grid"""
        etype = 3  # vtkLine().GetCellType()
        grid = self.create_grid_from_nodes_elements_etype(
            name, nodes, elements,  etype, color,
            point_size=point_size, line_width=line_width,
            opacity=opacity,
            representation=representation, add=add)
        return grid

    def create_grid_from_nodes_elements_etype(
        self, name: str,
        nodes: np.ndarray,
        elements: np.ndarray,
        etype: int,
        color: ColorFloat,
        point_size: int=5, line_width: int=5, opacity: float=1.,
        representation: str='wire',
        add: bool=True) -> Optional[vtkUnstructuredGrid]:
        """Makes a generic grid of constant type"""
        gui = self.gui
        gui.create_alternate_vtk_grid(
            name, color=color,
            point_size=point_size, line_width=line_width,
            opacity=opacity, representation=representation)

        nnodes = nodes.shape[0]
        nelements = elements.shape[0]
        if nnodes == 0 or nelements == 0:
            return None

        #print(f'adding quad_grid {name!r}; nnodes={nnodes:d} nelements={nelements:d}')
        assert isinstance(nodes, np.ndarray), type(nodes)

        points = numpy_to_vtk_points(nodes)
        grid = gui.alt_grids[name]
        grid.SetPoints(points)

        create_vtk_cells_of_constant_element_type(grid, elements, etype)
        if add:
            gui._add_alt_actors({name : gui.alt_grids[name]})

            #if name in self.geometry_actors:
        gui.geometry_actors[name].Modified()
        return grid

    def create_plane_actor_from_points(self, center: np.ndarray,
                                       i: np.ndarray,
                                       k: np.ndarray,
                                       dim_max: float,
                                       color: Optional[ColorFloat]=None,
                                       opacity: float=1.0,
                                       representation: str='surface',
                                       actor_name: str='plane') -> vtkActor:
        """
        This is used by the cutting plane tool and the shear/moment/torque tool.

            ^ k
            |
            |

           4+------+3
            |      |
            p1  c  p2
            |      |
           1+------+2 ----> i

        """
        if color is None:
            color = RED_FLOAT
        shift = 1.1
        dshift = (shift - 1) / 2.
        half_shift = 0.5 + dshift
        delta = half_shift * dim_max
        #dim_xy = shift * dim_max

        #n1 = 1 - dim_max * (dshift * i + half_shift * k)
        #n2 = n1 + shift * dim_max * i
        #n3 = n2 + shift * dim_max * k
        #n4 = n1 + shift * dim_max * k
        n1 = center - delta * i - delta * k
        n2 = center + delta * i - delta * k
        n3 = center + delta * i + delta * k
        n4 = center - delta * i + delta * k

        x = np.linspace(0., 1., num=10)
        y = x
        gui = self.gui
        if actor_name in gui.alt_grids:
            plane_actor = gui.plane_actor
            add = False
            #alt_grid =
            #plane_source = vtkPlaneSource()
            #self.rend.AddActor(plane_actor)
            #self.plane_actor = plane_actor
        else:
            add = True
            alt_grid = vtkUnstructuredGrid()
            gui.alt_grids[actor_name] = alt_grid

            mapper = vtkDataSetMapper()
            mapper.SetInputData(alt_grid)
            plane_actor = vtkActor()
            plane_actor.SetMapper(mapper)

            #plane_source = self.plane_source
            #plane_actor = self.plane_actor
            gui.plane_actor = plane_actor
            gui.rend.AddActor(plane_actor)

        nodes, elements = points_elements_from_quad_points(n1, n2, n3, n4, x, y)
        gui.set_quad_grid(actor_name, nodes, elements, color=color,
                          line_width=1, opacity=opacity, representation=representation,
                          add=add)
        #plane_actor.Modified()
        return plane_actor
