"""
these are functions not exposed in the GUI that are still useful
"""
from copy import deepcopy
import numpy as np
from numpy import issubdtype

from vtk import vtkLODActor, vtkDataSetMapper, vtkTypeFloat32Array
#from vtkmodules.vtkRenderingLOD import vtkLODActor

from pyNastran.gui.vtk_common_core import VTK_INT, VTK_FLOAT
from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk


def flip_actor_visibility(actor: vtkLODActor):
    """flips the visibility for an actor"""
    is_visible = actor.GetVisibility()
    actor.SetVisibility(not is_visible)
    actor.Modified()

def add_actors_to_gui(gui,
                      actors: list[vtkLODActor],
                      render: bool=True):
    """adds multiple vtk actors"""
    if not len(actors):
        return
    renderer = gui.rend
    for actor in actors:
        renderer.AddActor(actor)
    if render:
        renderer.Render()

def remove_actors_from_gui(gui,
                           actors: list[vtkLODActor],
                           render: bool=True,
                           force_render: bool=False):
    """removes multiple vtk actors"""
    if not (len(actors) or force_render):
        return
    renderer = gui.rend
    for actor in actors:
        renderer.RemoveActor(actor)
    if render:
        renderer.Render()

def set_vtk_fringe(grid_mapper: vtkDataSetMapper,
                   case_og: np.ndarray,
                   vector_size: int, phase: float) -> vtkTypeFloat32Array:
    """
    https://pyscience.wordpress.com/2014/09/06/numpy-to-vtk-converting-your-numpy-arrays-to-vtk-arrays-and-files/
    """
    case, vtk_data_type = set_grid_mapper(
        grid_mapper, case_og,
        vector_size, phase)
    #if 0: # nan testing
        #if case.dtype.name == 'float32':
            #case[50] = np.float32(1) / np.float32(0)
        #else:
            #case[50] = np.int32(1) / np.int32(0)

    if vector_size == 1:
        if case.flags.contiguous:
            case2 = case
        else:
            case2 = deepcopy(case)
        grid_result = numpy_to_vtk(
            num_array=case2,
            deep=True,
            array_type=vtk_data_type
        )
        #print('grid_result = %s' % grid_result)
        #print('max2 =', grid_result.GetRange())
    else:
        # vector_size=3
        if case.flags.contiguous:
            case2 = case
        else:
            case2 = deepcopy(case)
        grid_result = numpy_to_vtk(
            num_array=case2,
            deep=True,
            array_type=vtk_data_type
        )
    return grid_result

def set_grid_mapper(grid_mapper: vtkDataSetMapper,
                    case: np.ndarray,
                    vector_size: int,
                    phase: float) -> tuple[np.ndarray, int]:
    assert isinstance(case, np.ndarray), case
    if issubdtype(case.dtype, np.integer):
        vtk_data_type = VTK_INT
        grid_mapper.InterpolateScalarsBeforeMappingOn()
    elif issubdtype(case.dtype, np.floating):
        vtk_data_type = VTK_FLOAT
        grid_mapper.InterpolateScalarsBeforeMappingOff()
    elif case.dtype.name == 'complex64':
        if phase:
            phaser = np.radians(phase)
            case = (np.cos(phaser) * case.real + np.sin(phaser) * case.imag).real
        else:
            case = case.real
        vtk_data_type = VTK_FLOAT
        grid_mapper.InterpolateScalarsBeforeMappingOff()
    else:  # pragma: no cover
        raise NotImplementedError(case.dtype.type)
    return case, vtk_data_type
