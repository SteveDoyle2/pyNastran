# coding: utf-8
from __future__ import annotations
import os
import traceback
from typing import List, Optional, TYPE_CHECKING

import numpy as np
import vtk

from qtpy.compat import getsavefilename

from pyNastran.utils import check_path
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.gui import font_file
from pyNastran.gui.gui_objects.coord_properties import CoordProperties
from pyNastran.gui.utils.vtk.vtk_utils import numpy_to_vtk_points, update_axis_text_size
from pyNastran.gui.utils.load_results import load_user_geom
from pyNastran.gui.gui_objects.alt_geometry_storage import AltGeometry
from pyNastran.femutils.io import loadtxt_nice
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui_objects.settings import Settings


class ToolActions:
    def __init__(self, gui):
        self.gui = gui
        self.itext = 0

    #---------------------------------------------------------------------------
    def export_case_data(self, icases: Optional[List[int]]=None) -> None:
        """exports CSVs of the requested cases"""
        if icases is None:
            icases = self.gui.result_cases.keys()
        elif isinstance(icases, integer_types):
            icases = [icases]

        for icase in icases:
            (obj, (i, name)) = self.gui.result_cases[icase]
            subcase_id = obj.subcase_id
            location = obj.get_location(i, name)

            case = obj.get_result(i, name)
            if case is None:
                continue # normals
            subtitle, label = self.gui.get_subtitle_label(subcase_id)
            label2 = obj.get_header(i, name)
            data_format = obj.get_data_format(i, name)
            unused_vector_size = obj.get_vector_size(i, name)
            print(subtitle, label, label2, location, name)

            word, eids_nids = self.gui.get_mapping_for_location(location)

            # fixing cast int data
            header = '%s(%%i),%s(%s)' % (word, label2, data_format)
            if 'i' in data_format and isinstance(case.dtype, np.floating):
                header = '%s(%%i),%s' % (word, label2)

            fname = '%s_%s.csv' % (icase, _remove_invalid_filename_characters(name))
            out_data = np.column_stack([eids_nids, case])
            np.savetxt(fname, out_data, delimiter=',', header=header, fmt=b'%s')
        if icases:
            self.gui.log_command(f'self.export_case_data(icases={icases})')

    #---------------------------------------------------------------------------
    def create_corner_axis(self) -> None:
        """creates the axes that sits in the corner"""
        if not self.gui.run_vtk:
            return
        axes = vtk.vtkAxesActor()
        self.gui.corner_axis = vtk.vtkOrientationMarkerWidget()
        self.gui.corner_axis.SetOrientationMarker(axes)
        self.gui.corner_axis.SetInteractor(self.vtk_interactor)
        self.gui.corner_axis.SetEnabled(1)
        self.gui.corner_axis.InteractiveOff()

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
        self.settings.dim_max = dim_max
        coord_scale = self.settings.coord_scale * dim_max
        coord_text_scale = self.settings.coord_text_scale
        linewidth = self.settings.coord_linewidth

        create_actor = True
        if coord_id in self.gui.axes:
            axes = self.gui.axes[coord_id]
            create_actor = False
        else:
            axes = vtk.vtkAxesActor()
            axes.DragableOff()
            axes.PickableOff()

        transform = make_vtk_transform(origin, matrix_3x3)
        _set_base_axes(axes, transform, coord_type, label,
                       coord_scale, coord_text_scale, linewidth)

        self.gui.transform[coord_id] = transform
        self.gui.axes[coord_id] = axes

        is_visible = False
        if label == '':
            label = 'Global XYZ'
            is_visible = True
        else:
            label = 'Coord %s' % label
        self.gui.geometry_properties[label] = CoordProperties(label, coord_type, is_visible, coord_scale)
        self.gui.geometry_actors[label] = axes
        if create_actor:
            self.rend.AddActor(axes)

    #---------------------------------------------------------------------------
    def create_text(self, position, label, text_size: int=18):
        """creates the lower left text actors"""
        text_actor = vtk.vtkTextActor()

        text_actor.SetInput(label)
        text_prop = text_actor.GetTextProperty()
        set_vtk_property_to_unicode(text_prop, font_file)
        #text_prop.SetFontFamilyToArial()
        text_prop.SetFontSize(int(text_size))
        text_prop.SetColor(self.settings.text_color)
        text_actor.SetDisplayPosition(*position)

        text_actor.VisibilityOff()

        # assign actor to the renderer
        self.rend.AddActor(text_actor)
        self.gui.text_actors[self.itext] = text_actor
        self.itext += 1

    def update_text_actors(self, subcase_id, subtitle,
                           imin, min_value,
                           imax, max_value, label, location):
        """
        Updates the text actors in the lower left

        Max:  1242.3
        Min:  0.
        Subcase: 1 Subtitle:
        Label: SUBCASE 1; Static

        """
        if location == 'node':
            nodes = self.gui.node_ids
            min_msgi = f'Node: %i' % nodes[imin]
            max_msgi = f'Node: %i' % nodes[imax]
        elif location == 'centroid':
            elements = self.gui.element_ids
            min_msgi = f'Element: %i' % elements[imin]
            max_msgi = f'Element: %i' % elements[imax]
        else:
            raise NotImplementedError(location)


        if isinstance(max_value, integer_types):
            max_msg = 'Max:  %i' % max_value
            min_msg = 'Min:  %i' % min_value
        elif isinstance(max_value, str):
            max_msg = 'Max:  %s' % str(max_value)
            min_msg = 'Min:  %s' % str(min_value)
        elif (isinstance(max_value, float) or
              hasattr(max_value, 'dtype') and max_value.dtype.name in ['float32']):
            max_msg = 'Max:  %g' % max_value
            min_msg = 'Min:  %g' % min_value

        elif hasattr(max_value, 'dtype') and max_value.dtype.name in ['complex64']:
            raise RuntimeError(f'{max_value.dtype.name} should be a magnitude')
            #max_msg = 'Max:  %g, %gj' % (max_value.real, max_value.imag)
            #min_msg = 'Min:  %g, %gj' % (min_value.real, min_value.imag)
        else:
            max_msg = 'Max:  %g' % max_value
            min_msg = 'Min:  %g' % min_value

        texts = [
            max_msg + '; %s' % max_msgi,
            min_msg + '; %s' % min_msgi,
            'Subcase: %s Subtitle: %s' % (subcase_id, subtitle),
        ]
        if label:
            texts.append('Label: %s' % label)

        ntext = len(texts)
        for itext, text in enumerate(texts):
            self.gui.text_actors[itext].SetInput(text)

        if ntext == 4:  # label
            self.gui.text_actors[3].VisibilityOn()
        else:
            self.gui.text_actors[3].VisibilityOff()

    def turn_text_off(self):
        """turns all the text actors off"""
        for text in self.gui.text_actors.values():
            text.VisibilityOff()

    def turn_text_on(self):
        """turns all the text actors on"""
        for text in self.gui.text_actors.values():
            text.VisibilityOn()

    #---------------------------------------------------------------------------
    def on_take_screenshot(self, fname: Optional[str]=None,
                           magnify: Optional[int]=None,
                           show_msg: bool=True):
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
        render_large = vtk.vtkRenderLargeImage()
        render_large.SetInput(self.rend)

        out = self._screenshot_setup(magnify, render_large)
        line_widths0, point_sizes0, coord_scale0, coord_text_scale0, line_width0, axes_actor, magnify = out

        nam, ext = os.path.splitext(fname)
        ext = ext.lower()
        for nam, exts, obj in (('PostScript', ['.ps'], vtk.vtkPostScriptWriter),
                               ("BMP", ['.bmp'], vtk.vtkBMPWriter),
                               ('JPG', ['.jpg', '.jpeg'], vtk.vtkJPEGWriter),
                               ("TIFF", ['.tif', '.tiff'], vtk.vtkTIFFWriter)):
            if flt == nam:
                fname = fname if ext in exts else fname + exts[0]
                writer = obj()
                break
        else:
            fname = fname if ext == '.png' else fname + '.png'
            writer = vtk.vtkPNGWriter()

        writer.SetInputConnection(render_large.GetOutputPort())
        writer.SetFileName(fname)
        writer.Write()

        #self.log_info("Saved screenshot: " + fname)
        if show_msg:
            self.gui.log_command('on_take_screenshot(%r, magnify=%s)' % (fname, magnify))
        self._screenshot_teardown(line_widths0, point_sizes0,
                                  coord_scale0, coord_text_scale0, line_width0, axes_actor)

    def _get_screenshot_filename(self, fname: Optional[str]):
        """helper method for ``on_take_screenshot``"""
        if fname is None or fname is False:
            filt = ''
            default_filename = ''

            title = ''
            if self.gui.title is not None:
                title = self.gui.title

            if self.gui.out_filename is None:
                default_filename = ''
                if self.gui.infile_name is not None:
                    base, ext = os.path.splitext(os.path.basename(self.gui.infile_name))
                    default_filename = self.gui.infile_name
                    default_filename = base + '.png'
            else:
                base, ext = os.path.splitext(os.path.basename(self.gui.out_filename))
                default_filename = title + '_' + base + '.png'

            file_types = (
                'PNG Image *.png (*.png);; '
                'JPEG Image *.jpg *.jpeg (*.jpg, *.jpeg);; '
                'TIFF Image *.tif *.tiff (*.tif, *.tiff);; '
                'BMP Image *.bmp (*.bmp);; '
                'PostScript Document *.ps (*.ps)')

            title = 'Choose a filename and type'
            fname, flt = getsavefilename(parent=self.gui, caption=title, basedir='',
                                         filters=file_types, selectedfilter=filt,
                                         options=None)
            if fname in [None, '']:
                return None, None
            #print("fname=%r" % fname)
            #print("flt=%r" % flt)
        else:
            unused_base, ext = os.path.splitext(os.path.basename(fname))
            if ext.lower() in ['png', 'jpg', 'jpeg', 'tif', 'tiff', 'bmp', 'ps']:
                flt = ext.lower()
            else:
                flt = 'png'
        return fname, flt

    def _screenshot_setup(self, magnify: Optional[int], render_large):
        """helper method for ``on_take_screenshot``"""
        if magnify is None:
            magnify_min = 1
            magnify = self.settings.magnify if self.settings.magnify > magnify_min else magnify_min

        if not isinstance(magnify, integer_types):
            msg = 'magnify=%r type=%s' % (magnify, type(magnify))
            raise TypeError(msg)
        self.settings.update_text_size(magnify=magnify)

        coord_scale0 = self.settings.coord_scale
        coord_text_scale0 = self.settings.coord_text_scale
        linewidth0 = self.settings.coord_linewidth
        self.settings.scale_coord(magnify, render=False)
        self.settings.update_coord_text_scale(coord_text_scale0*magnify, render=False)
        #self.settings.scale_coord_text(magnify, render=False)
        #self.settings.update_coord_scale(
            #coord_scale=coord_scale0*magnify,
            #coord_text_scale=coord_text_scale0*magnify,
            #linewidth=linewidth0*magnify,
            #render=False)
        render_large.SetMagnification(magnify)

        # multiply linewidth by magnify
        line_widths0 = {}
        point_sizes0 = {}
        for key, geom_actor in self.gui.geometry_actors.items():
            if isinstance(geom_actor, vtk.vtkActor):
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
            elif isinstance(geom_actor, vtk.vtkAxesActor):
                pass
            else:
                raise NotImplementedError(geom_actor)

        # hide corner axis
        axes_actor = self.gui.corner_axis.GetOrientationMarker()
        axes_actor.SetVisibility(False)
        return line_widths0, point_sizes0, coord_scale0, coord_text_scale0, linewidth0, axes_actor, magnify

    def _screenshot_teardown(self, line_widths0, point_sizes0,
                             coord_scale0, coord_text_scale0, linewidth0, axes_actor):
        """helper method for ``on_take_screenshot``"""
        self.settings.update_text_size(magnify=1.0)
        # show corner axes
        axes_actor.SetVisibility(True)

        # set linewidth back
        for key, geom_actor in self.gui.geometry_actors.items():
            if isinstance(geom_actor, vtk.vtkActor):
                prop = geom_actor.GetProperty()
                prop.SetLineWidth(line_widths0[key])
                prop.SetPointSize(point_sizes0[key])
                prop.Modified()
            elif isinstance(geom_actor, vtk.vtkAxesActor):
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
    def on_load_user_geom(self, csv_filename=None, name=None, color=None):
        """
        Loads a User Geometry CSV File of the form:

        #    id  x    y    z
        GRID, 1, 0.2, 0.3, 0.3
        GRID, 2, 1.2, 0.3, 0.3
        GRID, 3, 2.2, 0.3, 0.3
        GRID, 4, 5.2, 0.3, 0.3
        grid, 5, 5.2, 1.3, 2.3  # case insensitive

        #    ID, nodes
        BAR,  1, 1, 2
        TRI,  2, 1, 2, 3
        # this is a comment

        QUAD, 3, 1, 5, 3, 4
        QUAD, 4, 1, 2, 3, 4  # this is after a blank line

        #RESULT,4,CENTROID,AREA(%f),PROPERTY_ID(%i)
        # in element id sorted order: value1, value2
        #1.0, 2.0 # bar
        #1.0, 2.0 # tri
        #1.0, 2.0 # quad
        #1.0, 2.0 # quad

        #RESULT,NODE,NODEX(%f),NODEY(%f),NODEZ(%f)
        # same difference

        #RESULT,VECTOR3,GEOM,DXYZ
        # 3xN

        Parameters
        ----------
        csv_filename : str (default=None -> load a dialog)
            the path to the user geometry CSV file
        name : str (default=None -> extract from fname)
            the name for the user points
        color : (float, float, float)
            RGB values as 0.0 <= rgb <= 1.0

        """
        if csv_filename in [None, False]:
            title = 'Load User Geometry'
            csv_filename = self.gui._create_load_file_dialog(
                self.gui.wildcard_delimited + ';;STL (*.stl)', title)[1]
            if not csv_filename:
                return

        if color is None:
            # we mod the num_user_points so we don't go outside the range
            icolor = self.gui.num_user_points % len(self.gui.color_order)
            color = self.gui.color_order[icolor]
        if name is None:
            name = os.path.basename(csv_filename).rsplit('.', 1)[0]

        self._add_user_geometry(csv_filename, name, color)
        self.gui.log_command('on_load_user_geom(%r, %r, %s)' % (
            csv_filename, name, str(color)))

    def _add_user_geometry(self, csv_filename, name, color):
        """
        helper method for ``on_load_user_geom``

        A custom geometry can be the pyNastran custom form or an STL

        """
        if name in self.gui.geometry_actors:
            msg = 'Name: %s is already in geometry_actors\nChoose a different name.' % name
            raise ValueError(msg)
        if len(name) == 0:
            msg = 'Invalid Name: name=%r' % name
            raise ValueError(msg)

        point_name = name + '_point'
        geom_name = name + '_geom'

        grid_ids, xyz, bars, tris, quads = load_user_geom(csv_filename, self.gui.log,
                                                          encoding='latin1')
        nbars = len(bars)
        ntris = len(tris)
        nquads = len(quads)
        nelements = nbars + ntris + nquads
        self.gui.create_alternate_vtk_grid(point_name, color=color, opacity=1.0,
                                           point_size=5, representation='point')

        if nelements > 0:
            nid_map = {}
            i = 0
            for nid in grid_ids:
                nid_map[nid] = i
                i += 1
            self.gui.create_alternate_vtk_grid(geom_name, color=color, opacity=1.0,
                                               line_width=5, representation='toggle')

        # allocate
        nnodes = len(grid_ids)
        #self.alt_grids[point_name].Allocate(npoints, 1000)
        #if nelements > 0:
            #self.alt_grids[geom_name].Allocate(npoints, 1000)

        alt_grid = self.gui.alt_grids[point_name]
        geom_grid = self.gui.alt_grids[geom_name]

        add_user_geometry(
            alt_grid, geom_grid,
            xyz, nid_map, nnodes,
            bars, tris, quads,
            nelements, nbars, ntris, nquads)

        # create actor/mapper
        self._add_alt_geometry(alt_grid, point_name)
        if nelements > 0:
            self._add_alt_geometry(geom_grid, geom_name)

        # set representation to points
        #self.geometry_properties[point_name].representation = 'point'
        #self.geometry_properties[geom_name].representation = 'toggle'
        #actor = self.geometry_actors[name]
        #prop = actor.GetProperty()
        #prop.SetRepresentationToPoints()
        #prop.SetPointSize(4)

    #---------------------------------------------------------------------------
    def on_load_csv_points(self, csv_filename=None, name=None, color=None):
        """
        Loads a User Points CSV File of the form:

        1.0, 2.0, 3.0
        1.5, 2.5, 3.5

        Parameters
        -----------
        csv_filename : str (default=None -> load a dialog)
            the path to the user points CSV file
        name : str (default=None -> extract from fname)
            the name for the user points
        color : (float, float, float)
            RGB values as 0.0 <= rgb <= 1.0

        .. note:: no header line is required
        .. note:: nodes are in the global frame

        .. todo:: support changing the name
        .. todo:: support changing the color
        .. todo:: support overwriting points

        """
        is_failed = True
        if csv_filename in [None, False]:
            title = 'Load User Points'
            csv_filename = self.gui._create_load_file_dialog(self.gui.wildcard_delimited, title)[1]
            if not csv_filename:
                return is_failed
        if color is None:
            # we mod the num_user_points so we don't go outside the range
            icolor = self.gui.num_user_points % len(self.gui.color_order)
            color = self.gui.color_order[icolor]
        if name is None:
            sline = os.path.basename(csv_filename).rsplit('.', 1)
            name = sline[0]

        is_failed = self._add_user_points_from_csv(csv_filename, name, color)
        if not is_failed:
            self.gui.num_user_points += 1
            self.gui.log_command('on_load_csv_points(%r, %r, %s)' % (
                csv_filename, name, str(color)))
        return is_failed

    def _add_user_points_from_csv(self, csv_points_filename: str, name: str,
                                  color, point_size: int=4):
        """
        Helper method for adding csv nodes to the gui

        Parameters
        ----------
        csv_points_filename : str
            CSV filename that defines one xyz point per line
        name : str
            name of the geometry actor
        color : List[float, float, float]
            RGB values; [0. to 1.]
        point_size : int; default=4
            the nominal point size

        """
        is_failed = True
        try:
            check_path(csv_points_filename, 'csv_points_filename')
            # read input file
            try:
                user_points = np.loadtxt(csv_points_filename, comments='#', delimiter=',')
            except ValueError:
                user_points = loadtxt_nice(csv_points_filename, comments='#', delimiter=',')
                # can't handle leading spaces?
                #raise
        except ValueError as error:
            #self.log_error(traceback.print_stack(f))
            self.gui.log_error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)
            self.gui.log_error(str(error))
            return is_failed

        self._add_user_points(user_points, name, color, csv_points_filename,
                              point_size=point_size)
        is_failed = False
        return False

    def _add_user_points(self, user_points, name: str, color,
                         csv_points_filename: str='', point_size: int=4):
        """
        Helper method for adding csv nodes to the gui

        Parameters
        ----------
        user_points : (n, 3) float ndarray
            the points to add
        name : str
            name of the geometry actor
        color : List[float, float, float]
            RGB values; [0. to 1.]
        point_size : int; default=4
            the nominal point size

        """
        if name in self.gui.geometry_actors:
            msg = 'Name: %s is already in geometry_actors\nChoose a different name.' % name
            raise ValueError(msg)
        if len(name) == 0:
            msg = 'Invalid Name: name=%r' % name
            raise ValueError(msg)

        # create grid
        self.gui.create_alternate_vtk_grid(name, color=color, line_width=5, opacity=1.0,
                                           point_size=point_size, representation='point')

        npoints = user_points.shape[0]
        if npoints == 0:
            raise RuntimeError('npoints=0 in %r' % csv_points_filename)
        if len(user_points.shape) == 1:
            user_points = user_points.reshape(1, npoints)

        # allocate grid
        alt_grid = self.gui.alt_grids[name]
        alt_grid.Allocate(npoints, 1000)

        # set points
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(npoints)

        for i, point in enumerate(user_points):
            points.InsertPoint(i, *point)
            elem = vtk.vtkVertex()
            elem.GetPointIds().SetId(0, i)
            alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
        alt_grid.SetPoints(points)

        # create actor/mapper
        self._add_alt_geometry(alt_grid, name)

        # set representation to points
        self.gui.geometry_properties[name].representation = 'point'
        actor = self.gui.geometry_actors[name]
        prop = actor.GetProperty()
        prop.SetRepresentationToPoints()
        prop.RenderPointsAsSpheresOn()
        prop.SetLighting(False)
        #prop.SetInterpolationToFlat()
        prop.SetPointSize(point_size)

    #---------------------------------------------------------------------------
    def _add_alt_geometry(self, grid, name: str, color=None, line_width=None,
                          opacity=None, representation=None):
        """NOTE: color, line_width, opacity are ignored if name already exists"""
        has_geometry_actor = name in self.gui.geometry_actors

        is_pickable = self.gui.geometry_properties[name].is_pickable
        quad_mapper = vtk.vtkDataSetMapper()

        if has_geometry_actor:
            alt_geometry_actor = self.gui.geometry_actors[name]
            alt_geometry_actor.GetMapper().SetInputData(grid)
        else:
            quad_mapper.SetInputData(grid)
            alt_geometry_actor = vtk.vtkActor()
            if not is_pickable:
                alt_geometry_actor.PickableOff()
                alt_geometry_actor.DragableOff()

            alt_geometry_actor.SetMapper(quad_mapper)
            self.gui.geometry_actors[name] = alt_geometry_actor

        #geometryActor.AddPosition(2, 0, 2)
        if name in self.gui.geometry_properties:
            geom = self.gui.geometry_properties[name]
        else:
            geom = AltGeometry(self, name, color=color, line_width=line_width,
                               opacity=opacity, representation=representation)
            self.gui.geometry_properties[name] = geom

        color = geom.color_float
        opacity = geom.opacity
        point_size = geom.point_size
        representation = geom.representation
        line_width = geom.line_width
        #print('color_2014[%s] = %s' % (name, str(color)))
        assert isinstance(color[0], float), color
        assert color[0] <= 1.0, color

        prop = alt_geometry_actor.GetProperty()
        #prop.SetInterpolationToFlat()    # 0
        #prop.SetInterpolationToGouraud() # 1
        #prop.SetInterpolationToPhong()   # 2
        prop.SetDiffuseColor(color)
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
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()

        if geom.is_visible:
            alt_geometry_actor.VisibilityOn()
        else:
            alt_geometry_actor.VisibilityOff()

        #print('current_actors = ', self.geometry_actors.keys())
        alt_geometry_actor.Modified()

    #---------------------------------------------------------------------------
    def GetCamera(self):
        return self.rend.GetActiveCamera()

    @property
    def settings(self) -> Settings:
        return self.gui.settings

    @property
    def rend(self):
        return self.gui.rend

    @property
    def vtk_interactor(self):
        return self.gui.vtk_interactor


def add_user_geometry(alt_grid, geom_grid,
                      xyz, nid_map, nnodes,
                      bars, tris, quads,
                      nelements, nbars, ntris, nquads):
    """helper method for ``_add_user_geometry``"""
    # set points
    points = numpy_to_vtk_points(xyz, dtype='<f')

    if nelements > 0:
        for i in range(nnodes):
            elem = vtk.vtkVertex()
            elem.GetPointIds().SetId(0, i)
            alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            geom_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
    else:
        for i in range(nnodes):
            elem = vtk.vtkVertex()
            elem.GetPointIds().SetId(0, i)
            alt_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

    if nbars:
        for i, bar in enumerate(bars[:, 1:]):
            g1 = nid_map[bar[0]]
            g2 = nid_map[bar[1]]
            elem = vtk.vtkLine()
            elem.GetPointIds().SetId(0, g1)
            elem.GetPointIds().SetId(1, g2)
            geom_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

    if ntris:
        for i, tri in enumerate(tris[:, 1:]):
            g1 = nid_map[tri[0]]
            g2 = nid_map[tri[1]]
            g3 = nid_map[tri[2]]
            elem = vtk.vtkTriangle()
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
            elem = vtk.vtkQuad()
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


def set_vtk_property_to_unicode(prop, font_filei):
    prop.SetFontFile(font_filei)
    prop.SetFontFamily(vtk.VTK_FONT_FILE)

def make_vtk_transform(origin, matrix_3x3):
    """makes a vtkTransform"""
    transform = vtk.vtkTransform()
    if origin is None and matrix_3x3 is None:
        pass
    elif origin is not None and matrix_3x3 is None:
        #print('origin%s = %s' % (label, str(origin)))
        transform.Translate(*origin)
    elif matrix_3x3 is not None:  # origin can be None
        xform = np.eye(4, dtype='float32')
        xform[:3, :3] = matrix_3x3
        if origin is not None:
            xform[:3, 3] = origin
        transform.SetMatrix(xform.ravel())
    else:
        raise RuntimeError('unexpected coordinate system')
    return transform


def _set_base_axes(axes: vtk.vtkAxesActor,
                   transform: vtk.vtkTransform,
                   coord_type: str, label: str,
                   coord_scale: float, coord_text_scale: float, linewidth: int) -> None:
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

def _remove_invalid_filename_characters(basename):
    """
    Helper method for exporting cases of 12*I/t^3.csv,
    which have invalid characters.

    Invalid for Windows
     < (less than)
     > (greater than)
     : (colon - sometimes works, but is actually NTFS Alternate Data Streams)
     " (double quote)
     / (forward slash)
     \ (backslash)
     | (vertical bar or pipe)
     ? (question mark)
     * (asterisk)

    Invalid for Linux
     / (forward slash)

    .. todo:: do a check for linux

    """
    invalid_chars = ':*?<>|/\\'
    for char in invalid_chars:
        basename = basename.replace(char, '')
    return basename
