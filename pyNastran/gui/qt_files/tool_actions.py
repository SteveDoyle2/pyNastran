# coding: utf-8
from __future__ import annotations
import os
#import traceback
from typing import Union, Optional, TYPE_CHECKING

import numpy as np

#from vtk import (
    #vtkTransform,
    #vtkPostScriptWriter, vtkBMPWriter, vtkJPEGWriter, vtkTIFFWriter, vtkPNGWriter,
    #vtkRenderLargeImage,
    #vtkAxesActor,
    #vtkOrientationMarkerWidget,
    #vtkXMLUnstructuredGridWriter,
#)
from vtkmodules.vtkCommonDataModel import vtkCellData, vtkPointData
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkIOImage import vtkPostScriptWriter, vtkBMPWriter, vtkJPEGWriter, vtkTIFFWriter, vtkPNGWriter
from vtkmodules.vtkFiltersHybrid import vtkRenderLargeImage
from vtkmodules.vtkRenderingAnnotation import vtkAxesActor
from vtkmodules.vtkInteractionWidgets import vtkOrientationMarkerWidget
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridWriter

from pyNastran.gui.vtk_common_core import VTK_FONT_FILE
from pyNastran.gui.vtk_rendering_core import (
    vtkDataSetMapper, vtkPolyDataMapper,
    vtkCamera, vtkTextActor, vtkProp, vtkActor, vtkRenderer)

from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.utils.locale import func_str
from pyNastran.gui import font_file
from pyNastran.gui.utils.paths import remove_invalid_filename_characters
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.gui_objects.coord_properties import CoordProperties
from pyNastran.gui.utils.qt.dialogs import save_file_dialog
from pyNastran.gui.utils.vtk.vtk_utils import update_axis_text_size
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui import MainWindow
    from pyNastran.gui.gui_objects.settings import Settings


class ToolActions:
    def __init__(self, gui: MainWindow):
        self.gui = gui
        self.itext = 0

    #---------------------------------------------------------------------------
    def get_icases(self, icases: Union[int, list[int], None]=None) -> list[int]:
        gui = self.gui
        if icases is None:
            icases = gui.result_cases.keys()
        elif isinstance(icases, integer_types):
            icases = [icases]
        return icases

    def export_case_data(self, icases: Optional[list[int]]=None) -> None:
        """exports CSVs of the requested cases"""
        gui = self.gui
        log = self.gui.log
        icases2 = self.get_icases(icases)

        csv_filename = 'None'
        for icase in icases2:
            (obj, (i, name)) = gui.result_cases[icase]
            subcase_id = obj.subcase_id
            location = obj.get_location(i, name)

            #method = obj.get_methods(i, name)[0]
            fringe, vector = obj.get_fringe_vector_result(i, name)
            case = vector if vector is not None else fringe
            if case is None:
                continue # normals, grid point forces
            subtitle, label = gui.get_subtitle_label(subcase_id)
            label2 = obj.get_annotation(i, name)
            data_format = obj.get_data_format(i, name)
            unused_vector_size = obj.get_vector_size(i, name)
            title = obj.get_legend_title(i, name)
            print(subtitle, label, label2, location, name, title)

            word, eids_nids = gui.get_mapping_for_location(location)
            is_failed, csv_filename = _export_case(
                title, eids_nids, case, icase,
                word, label2, data_format, log)
            if not is_failed:
                print(csv_filename)

        if icases:
            gui.log_command(f'self.export_case_data(icases={icases})\n'
                            f'# -> {csv_filename}')

    #---------------------------------------------------------------------------
    def create_corner_axis(self) -> None:
        """creates the axes that sits in the corner"""
        if not self.gui.run_vtk:
            return
        axes = vtkAxesActor()
        corner_axis = vtkOrientationMarkerWidget()
        corner_axis.SetOrientationMarker(axes)
        corner_axis.SetInteractor(self.vtk_interactor)
        corner_axis.SetEnabled(1)
        corner_axis.InteractiveOff()
        self.gui.corner_axis = corner_axis

    def create_coordinate_system(self, coord_id: int, dim_max: float, label: str='',
                                 origin=None, matrix_3x3=None,
                                 coord_type: str='xyz') -> None:
        """
        Creates a coordinate system

        Parameters
        ----------
        coord_id : int
            the coordinate system id
        dim_max : float
            the max model dimension; 10% of the max will be used for the coord length
        label : str
            the coord id or other unique label (default is empty to indicate the global frame)
        origin : (3, ) ndarray/list/tuple
            the origin
        matrix_3x3 : (3, 3) ndarray
            a standard Nastran-style coordinate system
        coord_type : str
            a string of 'xyz', 'Rtz', 'Rtp' (xyz, cylindrical, spherical)
            that changes the axis names

        .. todo::  coord_type is not supported ('xyz' ONLY)
        .. todo::  Can only set one coordinate system

        .. seealso::
            http://en.wikipedia.org/wiki/Homogeneous_coordinates
            http://www3.cs.stonybrook.edu/~qin/courses/graphics/camera-coordinate-system.pdf
            http://www.vtk.org/doc/nightly/html/classvtkTransform.html#ad58b847446d791391e32441b98eff151

        """
        gui: MainWindow = self.gui
        settings: Settings = gui.settings
        settings.dim_max = dim_max
        coord_scale = settings.coord_scale * dim_max
        coord_text_scale = settings.coord_text_scale
        linewidth = settings.coord_linewidth

        create_actor = True
        if coord_id in gui.axes:
            axes = gui.axes[coord_id]
            create_actor = False
        else:
            axes = vtkAxesActor()
            axes.DragableOff()
            axes.PickableOff()

        transform = make_vtk_transform(origin, matrix_3x3)
        _set_base_axes(axes, transform, coord_type, label,
                       coord_scale, coord_text_scale, linewidth)

        gui.transform[coord_id] = transform
        gui.axes[coord_id] = axes

        is_visible = False
        if label == '':
            label = 'Global XYZ'
            is_visible = True
        else:
            label = f'Coord {label}'
        gui.geometry_properties[label] = CoordProperties(
            label, coord_type, is_visible, coord_scale)
        gui.geometry_actors[label] = axes
        if create_actor:
            self.rend.AddActor(axes)

    #---------------------------------------------------------------------------
    def create_text(self, position: np.ndarray,
                    label: str,
                    text_size: int=18):
        """creates the lower left text actors"""
        text_actor = vtkTextActor()

        text_actor.SetInput(label)
        text_prop = text_actor.GetTextProperty()
        set_vtk_property_to_unicode(text_prop, font_file)
        #text_prop.SetFontFamilyToArial()
        text_prop.SetFontSize(int(text_size))
        text_prop.SetColor(self.settings.corner_text_color)
        text_actor.SetDisplayPosition(*position)

        text_actor.VisibilityOff()

        # assign actor to the renderer
        self.rend.AddActor(text_actor)
        self.gui.corner_text_actors[self.itext] = text_actor
        self.itext += 1

    def update_corner_text_actors(self, location: str,
                                  subcase_id: int,
                                  imin: int, min_value: float,
                                  imax: int, max_value: float,
                                  subtitle: str='case=NA',
                                  label: str='NA') -> None:
        """
        Updates the corner text actors in the lower left

        Max:  1242.3
        Min:  0.
        Subcase: 1 Subtitle:
        Label: SUBCASE 1; Static

        """
        gui = self.gui
        min_msg, max_msg = self._get_corner_min_max_text(
            imin, min_value,
            imax, max_value, location)

        texts = [
            max_msg,
            min_msg,
            'Subcase: %s Subtitle: %s' % (subcase_id, subtitle),
        ]
        if label:
            texts.append('Label: %s' % label)

        ntext = len(texts)
        text_actors = gui.corner_text_actors
        for itext, text in enumerate(texts):
            text_actors[itext].SetInput(text)

        if ntext == 4:  # label
            text_actors[3].VisibilityOn()
        else:
            text_actors[3].VisibilityOff()

    def _get_corner_min_max_text(self,
                                 imin: int, min_value: float,
                                 imax: int, max_value: float,
                                 location: str) -> tuple[str, str]:
        if location == 'normal':
            min_msgi = str(min_value)
            max_msgi = str(max_value)
            return min_msgi, max_msgi

        gui = self.gui
        if location == 'node':
            nodes: np.ndarray = gui.node_ids
            min_msgi = f'Node: %d' % nodes[imin]
            max_msgi = f'Node: %d' % nodes[imax]
        elif location == 'centroid':
            elements: np.ndarray = gui.element_ids
            min_msgi = f'Element: %d' % elements[imin]
            max_msgi = f'Element: %d' % elements[imax]
        else:  # pragma: no cover
            raise NotImplementedError(location)


        if isinstance(max_value, integer_types):
            max_msg = 'Max:  %d' % max_value
            min_msg = 'Min:  %d' % min_value
        elif isinstance(max_value, str):
            max_msg = 'Max:  %s' % str(max_value)
            min_msg = 'Min:  %s' % str(min_value)
        elif (isinstance(max_value, float) or
              hasattr(max_value, 'dtype') and max_value.dtype.name in ['float32', 'float64']):
            max_msg = 'Max:  %s' % func_str(max_value)
            min_msg = 'Min:  %s' % func_str(min_value)

        elif hasattr(max_value, 'dtype') and max_value.dtype.name in ['complex64', 'complex128']:
            raise RuntimeError(f'{max_value.dtype.name} should be a magnitude')
            #max_msg = 'Max:  %g, %gj' % (max_value.real, max_value.imag)
            #min_msg = 'Min:  %g, %gj' % (min_value.real, min_value.imag)
        else:
            max_msg = 'Max:  %g' % max_value
            min_msg = 'Min:  %g' % min_value

        max_out = max_msg + '; %s' % max_msgi
        min_out = min_msg + '; %s' % min_msgi
        return min_out, max_out

    def turn_corner_text_off(self) -> None:
        """turns all the text actors off"""
        for text in self.gui.corner_text_actors.values():
            text.VisibilityOff()

    def turn_corner_text_on(self) -> None:
        """turns all the text actors on"""
        for text in self.gui.corner_text_actors.values():
            text.VisibilityOn()

    #---------------------------------------------------------------------------
    def on_take_screenshot(self, fname: Optional[str]=None,
                           magnify: Optional[int]=None,
                           show_msg: bool=True) -> None:
        """
        Take a screenshot of a current view and save as a file

        Parameters
        ----------
        fname : str; default=None
            None : pop open a window
            str : bypass the popup window
        magnify : int; default=None
            None : use self.settings.magnify
            int : resolution increase factor
        show_msg : bool; default=True
            log the command

        """
        fname, flt = self._get_screenshot_filename(fname)

        if not fname:
            return
        render_large = vtkRenderLargeImage()
        render_large.SetInput(self.rend)

        out = self._screenshot_setup(magnify, render_large)
        line_widths0, point_sizes0, coord_scale0, coord_text_scale0, line_width0, axes_actor, magnify = out

        nam, ext = os.path.splitext(fname)
        ext = ext.lower()
        for nam, exts, obj in (('PostScript', ['.ps'], vtkPostScriptWriter),
                               ("BMP", ['.bmp'], vtkBMPWriter),
                               ('JPG', ['.jpg', '.jpeg'], vtkJPEGWriter),
                               ("TIFF", ['.tif', '.tiff'], vtkTIFFWriter)):
            if flt == nam:
                fname = fname if ext in exts else fname + exts[0]
                writer = obj()
                break
        else:
            fname = fname if ext == '.png' else fname + '.png'
            writer = vtkPNGWriter()

        writer.SetInputConnection(render_large.GetOutputPort())
        writer.SetFileName(fname)
        writer.Write()

        #self.log_info("Saved screenshot: " + fname)
        if show_msg:
            self.gui.log_command(f'self.on_take_screenshot({fname!r}, magnify={magnify})')
        self._screenshot_teardown(line_widths0, point_sizes0,
                                  coord_scale0, coord_text_scale0, line_width0, axes_actor)

    def _get_screenshot_filename(self, fname: Optional[str]) -> tuple[str, str]:
        """helper method for ``on_take_screenshot``"""
        if fname not in {None, False}:
            unused_base, ext = os.path.splitext(os.path.basename(fname))
            if ext.lower() in ['png', 'jpg', 'jpeg', 'tif', 'tiff', 'bmp', 'ps']:
                flt = ext.lower()
            else:
                flt = 'png'
            return fname, flt

        #filt = ''
        #default_png_filename = ''

        title = ''
        gui = self.gui
        if gui.title is not None:
            title = self.gui.title

        #dirname_file = os.path.splitext(self.gui.infile_name)[0]
        #default_vtk_filename = f'{dirname_file}.vtu'
        base = ''
        if gui.out_filename is None:
            if self.gui.infile_name is not None:
                base, ext = os.path.splitext(os.path.basename(gui.infile_name))
                #default_png_filename = gui.infile_name
        else:
            base, ext = os.path.splitext(os.path.basename(gui.out_filename))

        default_png_filename = ''
        if base:
            if title:
                default_png_filename = f'{title}_{base}.png'
            else:
                default_png_filename = f'{base}.png'

        file_types = (
            'PNG Image *.png (*.png);; '
            'JPEG Image *.jpg *.jpeg (*.jpg, *.jpeg);; '
            'TIFF Image *.tif *.tiff (*.tif, *.tiff);; '
            'BMP Image *.bmp (*.bmp);; '
            'PostScript Document *.ps (*.ps)')

        title = 'Choose a filename and type'
        fname, flt = save_file_dialog(
            gui, title,
            default_png_filename, file_types)
        #fname, flt = getsavefilename(parent=gui, caption=title, basedir='',
                                     #filters=file_types, selectedfilter=filt,
                                     #options=None)
        if fname in {None, ''}:
            return '', ''
        self.gui.load_actions._set_last_dir(fname)
        #print("fname=%r" % fname)
        #print("flt=%r" % flt)
        return fname, flt

    def _screenshot_setup(self, magnify: Optional[int],
                          render_large: vtkRenderLargeImage) -> tuple[
            dict[str, int], dict[str, int],
            dict[str, float], dict[str, float],
            dict[str, int],
            vtkAxesActor, int]:
        """helper method for ``on_take_screenshot``"""
        settings: Settings = self.settings
        if magnify is None:
            magnify_min = 1
            magnify = settings.magnify if settings.magnify > magnify_min else magnify_min

        if not isinstance(magnify, integer_types):
            msg = 'magnify=%r type=%s' % (magnify, type(magnify))
            raise TypeError(msg)
        self.settings.update_corner_text_size(magnify=magnify)

        coord_scale0 = settings.coord_scale
        coord_text_scale0 = settings.coord_text_scale
        linewidth0 = settings.coord_linewidth
        settings.scale_coord(magnify, render=False)
        settings.update_coord_text_scale(coord_text_scale0*magnify, render=False)
        #settings.scale_coord_text(magnify, render=False)
        #settings.update_coord_scale(
            #coord_scale=coord_scale0*magnify,
            #coord_text_scale=coord_text_scale0*magnify,
            #linewidth=linewidth0*magnify,
            #render=False)
        render_large.SetMagnification(magnify)

        # multiply linewidth by magnify
        line_widths0 = {}
        point_sizes0 = {}
        for key, geom_actor in self.gui.geometry_actors.items():
            if isinstance(geom_actor, vtkActor):
                prop = geom_actor.GetProperty()
                line_width0 = prop.GetLineWidth()
                point_size0 = prop.GetPointSize()
                line_widths0[key] = line_width0
                point_sizes0[key] = point_size0
                line_width = line_width0 * magnify
                point_size = point_size0 * magnify
                prop.SetLineWidth(line_width)
                prop.SetPointSize(point_size)
                prop.Modified()
            elif isinstance(geom_actor, vtkAxesActor):
                pass
            else:
                raise NotImplementedError(geom_actor)

        # hide corner axis
        axes_actor = self.gui.corner_axis.GetOrientationMarker()
        axes_actor.SetVisibility(False)
        out = (
            line_widths0, point_sizes0, coord_scale0, coord_text_scale0,
            linewidth0, axes_actor, magnify,
        )
        return out

    def _screenshot_teardown(self, line_widths0: dict[str, int],
                             point_sizes0: dict[str, int],
                             coord_scale0: dict[str, float],
                             coord_text_scale0: float,
                             linewidth0: int,
                             axes_actor: vtkAxesActor) -> None:
        """helper method for ``on_take_screenshot``"""
        self.settings.update_corner_text_size(magnify=1.0)
        # show corner axes
        axes_actor.SetVisibility(True)

        # set linewidth back
        for key, geom_actor in self.gui.geometry_actors.items():
            if isinstance(geom_actor, vtkActor):
                prop = geom_actor.GetProperty()
                prop.SetLineWidth(line_widths0[key])
                prop.SetPointSize(point_sizes0[key])
                prop.Modified()
            elif isinstance(geom_actor, vtkAxesActor):
                pass
            else:
                raise NotImplementedError(geom_actor)
        self.settings.scale_coord(magnify=1.0, render=False)
        self.settings.update_coord_text_scale(coord_text_scale0, render=True)
        #self.settings.scale_coord_text(magnify=1.0, render=True)
        #self.settings.update_coord_scale(coord_scale=coord_scale0,
                                         #coord_text_scale=coord_text_scale0,
                                         #linewidth=linewidth0,
                                         #render=True)

    #---------------------------------------------------------------------------
    def on_save_vtk(self, vtk_filename: Optional[str]=None) -> bool:
        """
        The result of "Export VTK..."
        """
        is_failed = True
        gui = self.gui
        grid = gui.grid
        if grid is None:
            return is_failed

        if vtk_filename in {None, False}:
            title = 'Select the VTK file name for export'
            wildcard_delimited = 'VTK (*.vtu)'
            dirname_file = os.path.splitext(self.gui.infile_name)[0]
            default_vtk_filename = f'{dirname_file}.vtu'
            vtk_filename, wildcard = save_file_dialog(
                gui, title,
                default_vtk_filename, wildcard_delimited)
            #assert wildcard == 'VTK (*.vtu; *.vtk)', wildcard 'VTK (*.vtu; *.vtk)'
            if not vtk_filename:
                return is_failed
            self.gui.load_actions._set_last_dir(vtk_filename)

        vtk_ugrid = self._get_vtk_ugrid()
        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName(vtk_filename)
        writer.SetInputData(vtk_ugrid)
        writer.Write()

        is_failed = False
        return is_failed

    def _get_vtk_ugrid(self) -> vtkUnstructuredGrid:
        """gets the vtkUnstructuredGrid with the results loaded"""
        gui = self.gui
        log = gui.log
        if gui.format == 'nastran':
            from pyNastran.converters.nastran.nastran_to_vtk import save_nastran_results, add_vtk_array
            vtk_ugrid = gui.grid
            save_nastran_results(gui, vtk_ugrid)
        else:
            used_titles: set[str] = set()
            vtk_ugrid = vtkUnstructuredGrid()
            point_data: vtkPointData = vtk_ugrid.GetPointData()
            cell_data: vtkCellData = vtk_ugrid.GetCellData()
            for case in gui.result_cases:
                if case.is_complex:
                    log.warning(f'skipping format={self.format!r}, case {str(case)!r} because it is complex')
                    continue
                if not isinstance(case, GuiResult):
                    log.warning(f'skipping format={self.format!r}, case {str(case)!r} because it is not a GuiResult')
                    continue
                vtk_array = case.save_vtk_result(used_titles)
                add_vtk_array(case.location, point_data, cell_data, vtk_array)
        return vtk_ugrid

    #---------------------------------------------------------------------------
    def add_alt_geometry(self, grid: vtkUnstructuredGrid,
                         name: str,
                         color: Optional[list[float]]=None,
                         line_width: Optional[int]=None,
                         opacity: Optional[float]=None,
                         representation: Optional[str]=None) -> None:
        """NOTE: color, line_width, opacity are ignored if name already exists"""
        gui = self.gui
        has_geometry_actor = name in gui.geometry_actors

        is_pickable = gui.geometry_properties[name].is_pickable
        quad_mapper = vtkDataSetMapper()

        if has_geometry_actor:
            alt_geometry_actor = gui.geometry_actors[name]
            alt_geometry_actor.GetMapper().SetInputData(grid)
        else:
            quad_mapper.SetInputData(grid)
            alt_geometry_actor = vtkActor()
            if not is_pickable:
                alt_geometry_actor.PickableOff()
                alt_geometry_actor.DragableOff()

            alt_geometry_actor.SetMapper(quad_mapper)
            gui.geometry_actors[name] = alt_geometry_actor

        #geometryActor.AddPosition(2, 0, 2)
        if name in gui.geometry_properties:
            geom = gui.geometry_properties[name]
        else:
            geom = AltGeometry(self, name, color=color, line_width=line_width,
                               opacity=opacity, representation=representation)
            gui.geometry_properties[name] = geom

        color_float: tuple[float, float, float] = geom.color_float
        opacity = geom.opacity
        point_size = geom.point_size
        representation = geom.representation
        line_width = geom.line_width
        #print('color_2014[%s] = %s' % (name, str(color_float)))
        assert isinstance(color_float[0], float), color_float
        assert color_float[0] <= 1.0, color_float

        prop = alt_geometry_actor.GetProperty()
        #prop.SetInterpolationToFlat()    # 0
        #prop.SetInterpolationToGouraud() # 1
        #prop.SetInterpolationToPhong()   # 2
        prop.SetDiffuseColor(color_float)
        prop.SetOpacity(opacity)
        #prop.Update()

        #print('prop.GetInterpolation()', prop.GetInterpolation()) # 1

        if representation == 'point':
            prop.SetRepresentationToPoints()
            prop.RenderPointsAsSpheresOn()
            prop.SetLighting(False)
            #prop.SetInterpolationToFlat()
            prop.SetPointSize(point_size)
        elif representation in ['surface', 'toggle']:
            prop.SetRepresentationToSurface()
            prop.SetLineWidth(line_width)
        elif representation == 'wire':
            prop.SetRepresentationToWireframe()
            prop.SetLineWidth(line_width)

        if not has_geometry_actor:
            self.rend.AddActor(alt_geometry_actor)
        vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()

        if geom.is_visible:
            alt_geometry_actor.VisibilityOn()
        else:
            alt_geometry_actor.VisibilityOff()

        #print('current_actors = ', self.geometry_actors.keys())
        alt_geometry_actor.Modified()

    #---------------------------------------------------------------------------
    def GetCamera(self) -> vtkCamera:
        return self.rend.GetActiveCamera()

    @property
    def settings(self) -> Settings:
        return self.gui.settings

    @property
    def rend(self) -> vtkRenderer:
        return self.gui.rend

    @property
    def vtk_interactor(self):
        return self.gui.vtk_interactor


def set_vtk_property_to_unicode(prop: vtkProp, font_filei: str) -> None:
    prop.SetFontFile(font_filei)
    prop.SetFontFamily(VTK_FONT_FILE)

def make_vtk_transform(origin: Optional[np.ndarray],
                       matrix_3x3: Optional[np.ndarray]) -> vtkTransform:
    """makes a vtkTransform"""
    transform = vtkTransform()
    if origin is None and matrix_3x3 is None:
        pass
    elif origin is not None and matrix_3x3 is None:
        #print('origin%s = %s' % (label, str(origin)))
        transform.Translate(*origin)
    elif matrix_3x3 is not None:  # origin can be None
        xform = xform3_to_xform4(matrix_3x3, origin)
        transform.SetMatrix(xform.ravel())
    else:
        raise RuntimeError('unexpected coordinate system')
    return transform

def xform3_to_xform4(matrix_3x3: np.ndarray,
                     origin: Optional[np.ndarray]) -> np.ndarray:
    """
    creates a 4x4 transform
            [[xform_3x3] [origin]]
    xform = [[0        ]    1    ]

    Parameters
    ----------
    origin : (3, ) float ndarray/list/tuple
        the origin
    matrix_3x3 : (3, 3) float ndarray
        a standard Nastran-style coordinate system

    https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html
    """
    xform = np.eye(4, dtype='float32')
    xform[:3, :3] = matrix_3x3
    if origin is not None:
        xform[:3, 3] = origin
    return xform

def _set_base_axes(axes: vtkAxesActor,
                   transform: vtkTransform,
                   coord_type: str, label: str,
                   coord_scale: float,
                   coord_text_scale: float,
                   linewidth: int) -> None:
    """sets the names of the axes"""
    #axes.GetLength() # pi
    #axes.GetNormalizedShaftLength() # (0.8, 0.8, 0.8)
    #axes.GetNormalizedTipLength() # (0.2, 0.2, 0.2)
    #axes.GetOrigin() # (0., 0., 0.)
    #axes.GetScale() # (1., 1., 1.)
    #axes.GetShaftType() # 1
    #axes.GetTotalLength() # (1., 1., 1.)

    yactor = axes.GetYAxisCaptionActor2D()
    zactor = axes.GetZAxisCaptionActor2D()

    axes.SetUserTransform(transform)
    axes.SetTotalLength(coord_scale, coord_scale, coord_scale)
    if coord_type == 'xyz':
        if label:
            xlabel = 'x%s' % label
            ylabel = 'y%s' % label
            zlabel = 'z%s' % label
            axes.SetXAxisLabelText(xlabel)
            axes.SetYAxisLabelText(ylabel)
            axes.SetZAxisLabelText(zlabel)
        #else:
            #xlabel = 'x'
            #ylabel = 'y'
            #zlabel = 'z'
    else:
        if coord_type == 'Rtz':  # cylindrical
            y = 'θ'
            x = 'R'
            z = 'z'
            if font_file:
                yprop = yactor.GetCaptionTextProperty()
                set_vtk_property_to_unicode(yprop, font_file)
            else:
                y = 't'

        elif coord_type == 'Rtp':  # spherical
            x = 'R'
            y = 'θ'
            z = 'Φ'
            if font_file:
                yprop = yactor.GetCaptionTextProperty()
                zprop = zactor.GetCaptionTextProperty()
                set_vtk_property_to_unicode(yprop, font_file)
                set_vtk_property_to_unicode(zprop, font_file)
            else:
                y = 't'
                z = 'p'
        else:  # pragma: no cover
            raise RuntimeError('invalid axis type; coord_type=%r' % coord_type)

        xlabel = '%s%s' % (x, label)
        ylabel = '%s%s' % (y, label)
        zlabel = '%s%s' % (z, label)
        axes.SetXAxisLabelText(xlabel)
        axes.SetYAxisLabelText(ylabel)
        axes.SetZAxisLabelText(zlabel)

    update_axis_text_size(axes, coord_text_scale)

    xaxis = axes.GetXAxisShaftProperty()
    yaxis = axes.GetYAxisShaftProperty()
    zaxis = axes.GetZAxisShaftProperty()
    #lw = xaxis.GetLineWidth()  #  1.0
    xaxis.SetLineWidth(linewidth)
    yaxis.SetLineWidth(linewidth)
    zaxis.SetLineWidth(linewidth)


def _export_case(name: str,
                 eids_nids: np.ndarray,
                 case: np.ndarray,
                 icase: int,
                 word: str, label2: str,
                 data_format: str,
                 log: SimpleLogger) -> tuple[bool, str]:
    """
    Writes a csv of a gui result in a form that you can load back in.
    """
    # handles 'Force\nIsElementOn', which would mess up the csv_filename
    name = name.replace('\n', '_')

    is_failed = True
    # fixing cast int data
    header = '%s(%%i),"%s(%s)"' % (word, label2, data_format)
    if 'i' in data_format and isinstance(case.dtype, np.floating):
        header = '%s(%%i),"%s"' % (word, label2)

    fname = '%s_%s.csv' % (icase, remove_invalid_filename_characters(name))
    try:
        out_data = np.column_stack([eids_nids, case])
    except ValueError:
        log.error(f'nodes/elements.shape={str(eids_nids.shape)} and '
                  f'case.shape={str(case.shape)} are not the same '
                  f'for {fname!r}')
        return is_failed, ''
    try:
        np.savetxt(fname, out_data, delimiter=',', header=header, fmt=b'%s')
    except UnicodeEncodeError:  # pragma: no cover
        try:
            np.savetxt(fname, out_data, delimiter=',', header=header, fmt='%s')
        except UnicodeEncodeError:  # pragma: no cover
            header = 'word,unicode_strikes_again'
            np.savetxt(fname, out_data, delimiter=',', header=header, fmt='%s')
    is_failed = False
    return is_failed, fname
