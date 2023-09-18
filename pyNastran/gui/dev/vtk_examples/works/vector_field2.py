import vtk
import numpy as np
from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid

def main():
    grid = vtkUnstructuredGrid()
    grid_mapper = vtk.vtkDataSetMapper()
    grid_mapper.SetInputData(grid)
    #grid_mapper.SetInputData(grid)


    nodes = np.array([
        [0., 0., 0.],
        [1., 0., 0.],
        [1., 1., 0.],
        [0., 2., 1.],
    ], dtype='float32')

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(4)
    points_array = numpy_to_vtk(
        num_array=nodes,
        deep=True,
        array_type=vtk.VTK_FLOAT,
    )
    nelements = 1
    grid.Allocate(nelements, 1000)
    grid.SetPoints(points)
    elem = vtk.vtkQuad()
    pts = elem.GetPointIds()
    pts.SetId(0, 0)
    pts.SetId(1, 1)
    pts.SetId(2, 2)
    pts.SetId(3, 3)
    grid.InsertNextCell(elem.GetCellType(), pts)
    grid.Modified()

    forces = np.array([
        [0., 0.1, 0.],
        [0., 0., 0.],
        [0., 0., 0.],
        [0., 0., .3],
    ], dtype='float32')

    rend = vtk.vtkRenderer()
    if 1:
        maskPts = vtk.vtkMaskPoints()
        maskPts.SetInputData(grid)

        arrow = vtk.vtkArrowSource()
        arrow.SetTipResolution(16)
        arrow.SetTipLength(0.3)
        arrow.SetTipRadius(0.1)

        glyph = vtk.vtkGlyph3D()
        glyph.SetSourceConnection(arrow.GetOutputPort())
        glyph.SetInputConnection(maskPts.GetOutputPort())
        glyph.SetVectorModeToUseNormal()
        glyph.SetScaleFactor(1)
        glyph.SetColorModeToColorByVector()
        glyph.SetScaleModeToScaleByVector()
        glyph.OrientOn()
        glyph.Update()


        glyph_mapper = vtk.vtkPolyDataMapper()
        glyph_mapper.SetInputConnection(glyph.GetOutputPort())
        glyph_mapper.SetScalarModeToUsePointFieldData()
        glyph_mapper.SetColorModeToMapScalars()
        glyph_mapper.ScalarVisibilityOn()
        glyph_mapper.SelectColorArray('Elevation')
        # Colour by scalars.
        #glyph_mapper.SetScalarRange(scalarRangeElevation)

        glyph_actor = vtk.vtkActor()
        glyph_actor.SetMapper(glyph_mapper)
        glyph_actor.RotateX(-45)
        glyph_actor.RotateZ(45)
        rend.AddViewProp(glyph_actor)
        #rend.AddActor(glyph_actor)


    geom_actor = vtk.vtkActor()
    geom_actor.SetMapper(grid_mapper)
    # ------------------------------------------------------------
    # Create the RenderWindow, Renderer and Interactor
    # ------------------------------------------------------------
    renWin = vtk.vtkRenderWindow()
    iren = vtk.vtkRenderWindowInteractor()

    renWin.AddRenderer(rend)
    iren.SetRenderWindow(renWin)

    # add actors
    #rend.AddViewProp(geom_actor)
    #rend.AddViewProp(edgeActor)
    rend.AddActor(geom_actor)
    #rend.AddViewProp(glyph_actor)
    #rend.AddActor2D(scalarBar)

    rend.SetBackground(0.7, 0.8, 1.0)
    renWin.SetSize(800, 800)
    renWin.Render()
    iren.Start()

main()
