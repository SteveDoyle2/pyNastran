"""
This example attempts to uses VTK efficiently to:
  - create an model (unstructured grid) with multiple element types
  - create arrows (glyphs) to represent some vector quantity (e.g., forces)
  - TODO: makes use of masking to remove small forces


Multiple challenges are addressed in this example:
  - how to add elements efficiently (even of the same type)
  - how to add multiple types efficiently
  - how to link an unstructured grid to a glyph (e.g., an arrow or a point)
    we avoid making two unstructured grids

Possible inefficienies:
  - could we use 1 actor? instead of using
    1 grid, but multiple mappers

There are flags that let you disable one or more features.

Tested on:
==========
Python 3.9, VTK 9.1
Python 3.10

Adapted from:
=============
http://docs.enthought.com/mayavi/mayavi/auto/example_unstructured_grid.html
http://vtk.1045678.n5.nabble.com/SetCells-for-multiple-cell-types-td5735172.html


To adapt:
=========
http://www.vtk.org/Wiki/VTK/Examples/Python/Visualization/ClampGlyphSizes
"""
import numpy as np
import vtk
from pyNastran.gui.vtk_interface import (
    vtkTetra, vtkHexahedron,
    vtkUnstructuredGrid, vtkCellArray)
from pyNastran.gui.vtk_renering_core import (
    vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor, vtkActor)

from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk, numpy_to_vtkIdTypeArray
from pyNastran.gui.utils.vtk.base_utils import VTK_VERSION


# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

make_glyphs = True
apply_color_to_glyph = True  # sets them as red otherwise
filter_small_forces = False # doesn't work


def mixed_type_unstructured_grid():
    """A slightly more complex example of how to generate an
    unstructured grid with different cell types.  Returns a created
    unstructured grid.
    """
    pts = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], # tetra
        [2, 0, 0], [3, 0, 0], [3, 1, 0], [2, 1, 0],
        [2, 0, 1], [3, 0, 1], [3, 1, 1], [2, 1, 1], # Hex
    ], dtype='float32')
    # shift the points so we can show both.
    pts[:, 1] += 2.0
    npoints = len(pts)
    forces = pts


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

    # add one element_type for each element
    tetra_type = vtkTetra().GetCellType() # VTK_TETRA == 10
    hex_type = vtkHexahedron().GetCellType() # VTK_HEXAHEDRON == 12
    cell_types = np.array([tetra_type, hex_type], dtype='int32')

    # Create the array of cells
    vtk_cells = vtkCellArray()
    vtk_cells_id_type = numpy_to_vtkIdTypeArray(cells, deep=1)

    # ncells = 2
    vtk_cells.SetCells(2, vtk_cells_id_type)

    # Now create the unstructured grid
    ug = vtkUnstructuredGrid()

    points_data = numpy_to_vtk(pts, deep=1)

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(npoints)
    points.SetData(points_data)
    ug.SetPoints(points)

    # Now just set the cell types and reuse the ug locations and cells

    #ug.SetCells(cell_types, cell_offsets, vtk_cell_array)
    ug.SetCells(
        numpy_to_vtk(
            cell_types,
            deep=1,
            array_type=vtk.vtkUnsignedCharArray().GetDataType(),
        ),
        numpy_to_vtk(
            cell_offsets,
            deep=1,
            array_type=vtk.vtkIdTypeArray().GetDataType()
        ),
        vtk_cells,
    )
    return ug, forces


def main():
    #print('VTK_VERSION = %r' % vtk.VTK_VERSION)
    ug, forces = mixed_type_unstructured_grid()
    forces_array = numpy_to_vtk(forces, deep=1)

    # get the magnitude of the force vectors
    scalars = np.linalg.norm(forces, axis=1)

    # scalars are defined from [0., 1.]
    scalars = scalars / scalars.max()
    print(scalars)
    #scalars[3] = 0.

    assert len(forces) == len(scalars)
    force_scalar_array = numpy_to_vtk(scalars, deep=1)

    grid_mapper = vtk.vtkDataSetMapper()
    #vtk_version = int(VTK_VERSION[0])
    grid_mapper.SetInputData(ug)

    if make_glyphs:
        glyphs = vtk.vtkGlyph3D()
        #if filter_small_forces:
            #glyphs.SetRange(0.5, 1.)

        glyphs.SetVectorModeToUseVector()
        if apply_color_to_glyph:
            glyphs.SetColorModeToColorByScale()
        #glyphs.SetColorModeToColorByScalar()
        #glyphs.SetColorModeToColorByVector()

        glyphs.ScalingOn()
        glyphs.ClampingOn()
        #glyphs.Update()

        glyphSource = vtk.vtkArrowSource()
        glyphSource.InvertOn()  # flip this arrow direction
        glyphs.SetInputData(ug)

        #glyphs.SetSource(glyphSource.GetOutput())
        glyphs.SetSourceConnection(glyphSource.GetOutputPort())
        #glyphs.SetScaleModeToDataScalingOff()
        #glyphs.SetScaleFactor(0.1)
        glyph_mapper = vtk.vtkPolyDataMapper()
        #glyph_mapper.SetInput(glyphs.GetOutput())
        glyph_mapper.SetInputConnection(glyphs.GetOutputPort())
        #glyph_mapper.SetVectors(forces_array)
        # You can set what arrays to use for the (scalars, vectors, normals, and colors) scalars
        # by using the SetInputArrayToProcess methods in vtkAlgorithm.
        # The first array is scalars, the next vectors, the next normals and finally color scalars.
        #glyphs.SetInputArrayToProcess()

        # Tell glyph which attribute arrays to use for what
        #glyph.SetInputArrayToProcess(0,0,0,0,'Elevation')           # scalars
        #glyph.SetInputArrayToProcess(1,0,0,0,'RTDataGradient')      # vectors
        #glyph.SetInputArrayToProcess(2,0,0,0,'nothing')             # normals
        #glyph.SetInputArrayToProcess(3,0,0,0,'RTData')              # colors

        #grid_mapper.SetInputData(ug)

        #if filter_small_forces:
            ## we need the data to be contiguous or VTK will fail
            ## np.copy (I think) always works, while deep=1 doesn't always work
            ##
            ## It seems like it makes more sense to just use numpy's copy method
            #show_data = np.copy(np.where(forces > 0.6)[0])
            #ids = numpy_to_vtkIdTypeArray(show_data, deep=0)

            #selection = vtk.vtkSelection()
            #selection_node = vtk.vtkSelectionNode()
            #extract_selection = vtk.vtkExtractSelection()

            #selection_node.SetFieldType(vtk.vtkSelectionNode.POINT)  # POINT/CELL
            #selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)

            #if vtk.VTK_MAJOR_VERSION <= 5:
                #extract_selection.SetInput(0, ug)
                #extract_selection.SetInput(1, selection)
            #else:
                #extract_selection.SetInputData(0, ug)
                #extract_selection.SetInputData(1, selection)

            ##ug.ShallowCopy(extract_selection.GetOutput())
            #selection.AddNode(selection_node)
            ##if 0:
                ##if flip_flag:
                    ##selection.RemoveAllNodes()
                    ##selection_node.SetSelectionList(ids)


        arrow_actor = vtk.vtkLODActor()
        arrow_actor.SetMapper(glyph_mapper)

        if not apply_color_to_glyph:
            prop = arrow_actor.GetProperty()
            prop.SetColor(1., 0., 0.)

    ug.GetPointData().SetScalars(force_scalar_array)
    if make_glyphs:
        ug.GetPointData().SetVectors(forces_array)
        if not apply_color_to_glyph:
            glyph_mapper.ScalarVisibilityOff()
    geom_actor = vtkActor()
    geom_actor.SetMapper(grid_mapper)


    # Setup renderer
    renderer = vtkRenderer()
    renderer.AddActor(geom_actor)
    if make_glyphs:
        renderer.AddActor(arrow_actor)
    renderer.ResetCamera()
    renderer.SetBackground(0.7, 0.8, 1.0)

    # Setup render window
    renderWindow = vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    # Setup render window
    renderWindow = vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    # Setup render window interactor
    renderWindowInteractor = vtkRenderWindowInteractor()
    style = vtk.vtkInteractorStyleImage()

    # Render and start interaction
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    renderWindowInteractor.Start()


if __name__ == '__main__':   # pragma: no cover
    main()
