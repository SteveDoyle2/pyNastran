# coding: utf-8
from vtkmodules.vtkRenderingLOD import vtkLODActor
from vtkmodules.vtkFiltersSources import vtkArrowSource
from vtkmodules.vtkFiltersCore import vtkGlyph3D

from pyNastran.gui.vtk_rendering_core import vtkPolyDataMapper
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid

def build_glyph(grid: vtkUnstructuredGrid) -> tuple[vtkArrowSource, vtkGlyph3D,
                                                    vtkPolyDataMapper, vtkLODActor]:
    """builds the glyph actor"""
    glyphs = vtkGlyph3D()
    #if filter_small_forces:
        #glyphs.SetRange(0.5, 1.)

    glyphs.SetVectorModeToUseVector()
    #apply_color_to_glyph = False
    #if apply_color_to_glyph:
    #glyphs.SetScaleModeToScaleByScalar()
    glyphs.SetScaleModeToScaleByVector()
    glyphs.SetColorModeToColorByScale()
    #glyphs.SetColorModeToColorByScalar()  # super tiny
    #glyphs.SetColorModeToColorByVector()  # super tiny

    glyphs.ScalingOn()
    glyphs.ClampingOn()
    #glyphs.Update()

    glyph_source = vtkArrowSource()
    #glyph_source.InvertOn()  # flip this arrow direction
    glyphs.SetInputData(grid)


    glyphs.SetSourceConnection(glyph_source.GetOutputPort())
    #glyphs.SetScaleModeToDataScalingOff()
    #glyphs.SetScaleFactor(10.0)  # bwb
    #glyphs.SetScaleFactor(1.0)  # solid-bending
    glyph_mapper = vtkPolyDataMapper()
    glyph_mapper.SetInputConnection(glyphs.GetOutputPort())
    glyph_mapper.ScalarVisibilityOff()

    arrow_actor = vtkLODActor()
    arrow_actor.SetMapper(glyph_mapper)

    prop = arrow_actor.GetProperty()
    prop.SetColor(1., 0., 0.)
    #self.grid.GetPointData().SetActiveVectors(None)
    arrow_actor.SetVisibility(False)
    return glyph_source, glyphs, glyph_mapper, arrow_actor
