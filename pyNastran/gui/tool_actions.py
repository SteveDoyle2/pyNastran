# coding: utf-8
from __future__ import print_function
import os

from six import iteritems, itervalues, string_types

import numpy as np
import vtk

from qtpy.compat import getsavefilename#, getopenfilename

from pyNastran.utils import integer_types
from pyNastran.gui.gui_objects.coord_properties import CoordProperties


class ToolActions(object):
    def __init__(self, gui):
        self.gui = gui
        self.itext = 0

    #---------------------------------------------------------------------------
    def export_case_data(self, icases=None):
        """exports CSVs of the requested cases"""
        if icases is None:
            icases = self.gui.result_cases.keys()
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

            fname = '%s_%s.csv' % (icase, name)
            out_data = np.column_stack([eids_nids, case])
            np.savetxt(fname, out_data, delimiter=',', header=header, fmt=b'%s')

    #---------------------------------------------------------------------------
    def create_corner_axis(self):
        """creates the axes that sits in the corner"""
        if not self.gui.run_vtk:
            return
        axes = vtk.vtkAxesActor()
        self.gui.corner_axis = vtk.vtkOrientationMarkerWidget()
        self.gui.corner_axis.SetOrientationMarker(axes)
        self.gui.corner_axis.SetInteractor(self.vtk_interactor)
        self.gui.corner_axis.SetEnabled(1)
        self.gui.corner_axis.InteractiveOff()

    def create_coordinate_system(self, coord_id, dim_max, label='', origin=None, matrix_3x3=None,
                                 coord_type='xyz'):
        """
        Creates a coordinate system

        Parameters
        ----------
        coord_id : float
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
        scale = self.settings.coord_scale * dim_max

        transform = make_vtk_transform(origin, matrix_3x3)

        create_actor = True
        if coord_id in self.gui.transform:
            axes = self.gui.axes[coord_id]
            create_actor = False
        else:
            axes = vtk.vtkAxesActor()
            axes.DragableOff()
            axes.PickableOff()

        #axes.GetLength() # pi
        #axes.GetNormalizedShaftLength() # (0.8, 0.8, 0.8)
        #axes.GetNormalizedTipLength() # (0.2, 0.2, 0.2)
        #axes.GetOrigin() # (0., 0., 0.)
        #axes.GetScale() # (1., 1., 1.)
        #axes.GetShaftType() # 1
        #axes.GetTotalLength() # (1., 1., 1.)

        axes.SetUserTransform(transform)
        axes.SetTotalLength(scale, scale, scale)
        if coord_type == 'xyz':
            if label:
                xlabel = u'x%s' % label
                ylabel = u'y%s' % label
                zlabel = u'z%s' % label
                axes.SetXAxisLabelText(xlabel)
                axes.SetYAxisLabelText(ylabel)
                axes.SetZAxisLabelText(zlabel)
        else:
            if coord_type == 'Rtz':  # cylindrical
                #x = u'R'
                #y = u'θ'
                #z = u'z'
                x = 'R'
                y = 't'
                z = 'z'

            elif coord_type == 'Rtp':  # spherical
                #x = u'R'
                #y = u'θ'
                #z = u'Φ'
                x = 'R'
                y = 't'
                z = 'p'
            else:
                raise RuntimeError('invalid axis type; coord_type=%r' % coord_type)

            xlabel = '%s%s' % (x, label)
            ylabel = '%s%s' % (y, label)
            zlabel = '%s%s' % (z, label)
            axes.SetXAxisLabelText(xlabel)
            axes.SetYAxisLabelText(ylabel)
            axes.SetZAxisLabelText(zlabel)

        self.gui.transform[coord_id] = transform
        self.gui.axes[coord_id] = axes

        is_visible = False
        if label == '':
            label = 'Global XYZ'
            is_visible = True
        else:
            label = 'Coord %s' % label
        self.gui.geometry_properties[label] = CoordProperties(label, coord_type, is_visible, scale)
        self.gui.geometry_actors[label] = axes
        if create_actor:
            self.rend.AddActor(axes)

    #---------------------------------------------------------------------------
    def create_text(self, position, label, text_size=18):
        """creates the lower left text actors"""
        text_actor = vtk.vtkTextActor()
        text_actor.SetInput(label)
        text_prop = text_actor.GetTextProperty()
        #text_prop.SetFontFamilyToArial()
        text_prop.SetFontSize(int(text_size))
        text_prop.SetColor(self.settings.text_color)
        text_actor.SetDisplayPosition(*position)

        text_actor.VisibilityOff()

        # assign actor to the renderer
        self.rend.AddActor(text_actor)
        self.gui.text_actors[self.itext] = text_actor
        self.itext += 1

    def update_text_actors(self, subcase_id, subtitle, min_value, max_value, label):
        """
        Updates the text actors in the lower left

        Max:  1242.3
        Min:  0.
        Subcase: 1 Subtitle:
        Label: SUBCASE 1; Static
        """
        if isinstance(max_value, integer_types):
            max_msg = 'Max:  %i' % max_value
            min_msg = 'Min:  %i' % min_value
        elif isinstance(max_value, string_types):
            max_msg = 'Max:  %s' % str(max_value)
            min_msg = 'Min:  %s' % str(min_value)
        else:
            max_msg = 'Max:  %g' % max_value
            min_msg = 'Min:  %g' % min_value
        self.gui.text_actors[0].SetInput(max_msg)
        self.gui.text_actors[1].SetInput(min_msg)
        self.gui.text_actors[2].SetInput('Subcase: %s Subtitle: %s' % (subcase_id, subtitle))

        if label:
            self.gui.text_actors[3].SetInput('Label: %s' % label)
            self.gui.text_actors[3].VisibilityOn()
        else:
            self.gui.text_actors[3].VisibilityOff()

    def turn_text_off(self):
        """turns all the text actors off"""
        for text in itervalues(self.gui.text_actors):
            text.VisibilityOff()

    def turn_text_on(self):
        """turns all the text actors on"""
        for text in itervalues(self.gui.text_actors):
            text.VisibilityOn()

    #---------------------------------------------------------------------------
    def on_take_screenshot(self, fname=None, magnify=None, show_msg=True):
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

        TODO: screenshot doesn't work well with the coordinate system text size
        """
        fname, flt = self._get_screenshot_filename(fname)

        if not fname:
            return
        render_large = vtk.vtkRenderLargeImage()
        render_large.SetInput(self.rend)

        out = self._screenshot_setup(magnify, render_large)
        line_widths0, point_sizes0, coord_scale0, axes_actor, magnify = out

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
        self._screenshot_teardown(line_widths0, point_sizes0, coord_scale0, axes_actor)

    def _get_screenshot_filename(self, fname):
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
                return
            #print("fname=%r" % fname)
            #print("flt=%r" % flt)
        else:
            base, ext = os.path.splitext(os.path.basename(fname))
            if ext.lower() in ['png', 'jpg', 'jpeg', 'tif', 'tiff', 'bmp', 'ps']:
                flt = ext.lower()
            else:
                flt = 'png'
        return fname, flt

    def _screenshot_setup(self, magnify, render_large):
        """helper method for ``on_take_screenshot``"""
        if magnify is None:
            magnify_min = 1
            magnify = self.settings.magnify if self.settings.magnify > magnify_min else magnify_min
        else:
            magnify = magnify

        if not isinstance(magnify, integer_types):
            msg = 'magnify=%r type=%s' % (magnify, type(magnify))
            raise TypeError(msg)
        self.settings.update_text_size(magnify=magnify)

        coord_scale0 = self.settings.coord_scale
        #coord_text_scale0 = self.settings.coord_text_scale
        self.settings.update_coord_scale(
            coord_scale=coord_scale0*magnify, render=False)
        render_large.SetMagnification(magnify)

        # multiply linewidth by magnify
        line_widths0 = {}
        point_sizes0 = {}
        for key, geom_actor in iteritems(self.gui.geometry_actors):
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
        return line_widths0, point_sizes0, coord_scale0, axes_actor, magnify

    def _screenshot_teardown(self, line_widths0, point_sizes0, coord_scale0, axes_actor):
        """helper method for ``on_take_screenshot``"""
        self.settings.update_text_size(magnify=1.0)

        # show corner axes
        axes_actor.SetVisibility(True)

        # set linewidth back
        for key, geom_actor in iteritems(self.gui.geometry_actors):
            if isinstance(geom_actor, vtk.vtkActor):
                prop = geom_actor.GetProperty()
                prop.SetLineWidth(line_widths0[key])
                prop.SetPointSize(point_sizes0[key])
                prop.Modified()
            elif isinstance(geom_actor, vtk.vtkAxesActor):
                pass
            else:
                raise NotImplementedError(geom_actor)
        self.settings.update_coord_scale(coord_scale=coord_scale0, render=True)

    #---------------------------------------------------------------------------
    # camera
    def update_camera(self, code):
        camera = self.GetCamera()
        #print("code =", code)
        if code == '+x':  # set x-axis
            # +z up
            # +y right
            # looking forward
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(1., 0., 0.)
        elif code == '-x':  # set x-axis
            # +z up
            # +y to the left (right wing)
            # looking aft
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(-1., 0., 0.)

        elif code == '+y':  # set y-axis
            # +z up
            # +x aft to left
            # view from right wing
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(0., 1., 0.)
        elif code == '-y':  # set y-axis
            # +z up
            # +x aft to right
            # view from left wing
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(0., -1., 0.)

        elif code == '+z':  # set z-axis
            # +x aft
            # +y up (right wing up)
            # top view
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 1., 0.)
            camera.SetPosition(0., 0., 1.)
        elif code == '-z':  # set z-axis
            # +x aft
            # -y down (left wing up)
            # bottom view
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., -1., 0.)
            camera.SetPosition(0., 0., -1.)
        else:
            self.gui.log_error('invalid camera code...%r' % code)
            return
        self._update_camera(camera)
        self.rend.ResetCamera()
        self.gui.log_command('update_camera(%r)' % code)

    def _update_camera(self, camera=None):
        if camera is None:
            camera = self.GetCamera()
        camera.Modified()
        self.vtk_interactor.Render()

    #---------------------------------------------------------------------------
    def GetCamera(self):
        return self.rend.GetActiveCamera()

    @property
    def settings(self):
        return self.gui.settings

    @property
    def rend(self):
        return self.gui.rend

    @property
    def vtk_interactor(self):
        return self.gui.vtk_interactor


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
