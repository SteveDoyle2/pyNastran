from __future__ import print_function
import os

from six import iteritems#, itervalues, string_types,

import vtk

from qtpy.compat import getsavefilename#, getopenfilename

from pyNastran.utils import integer_types


class ToolActions(object):
    def __init__(self, gui):
        self.gui = gui
        self.itext = 0

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

    def on_take_screenshot(self, fname=None, magnify=None, show_msg=True):
        """see gui for definition"""
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
            fname, flt = getsavefilename(parent=self, caption=title, basedir='',
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

        if not fname:
            return
        render_large = vtk.vtkRenderLargeImage()
        render_large.SetInput(self.rend)

        line_widths0, point_sizes0, coord_scale0, axes_actor, magnify = self._screenshot_setup(
            magnify, render_large)

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
        self.settings.update_coord_scale(coord_scale=coord_scale0*magnify, render=False)
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

    def rotate(self, rotate_deg, render=True):
        """see the gui"""
        camera = self.GetCamera()
        camera.Roll(-rotate_deg)
        camera.Modified()
        if render:
            self.vtk_interactor.Render()
        self.gui.log_command('rotate(%s)' % rotate_deg)

    def zoom(self, value):
        camera = self.GetCamera()
        camera.Zoom(value)
        camera.Modified()
        self.vtk_interactor.Render()
        self.gui.log_command('zoom(%s)' % value)

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
    def rend(self):
        return self.gui.rend

    @property
    def vtk_interactor(self):
        return self.gui.vtk_interactor

    @property
    def settings(self):
        return self.gui.settings
