"""interface to the ScalarBar"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from vtkmodules.vtkRenderingCore import vtkColorTransferFunction, vtkTextProperty
from vtkmodules.vtkRenderingAnnotation import vtkScalarBarActor

from pyNastran.gui.utils.colormaps import colormap_dict, RGB_MAPS, HSV_MAPS
from pyNastran.gui.qt_files.tool_actions import set_vtk_property_to_unicode
from pyNastran.gui import bi_font_file
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.gui.gui_objects.settings import Settings

NCOLORS_DEFAULT = 11

class ScalarBar:
    """defines the ScalarBar at the side of the vtk panel"""
    def set_visibility(self, is_visible: bool) -> None:
        """show/hide the scalar bar"""
        #print('is_visible=%s; is_shown=%s' % (is_visible, self.is_shown))
        if is_visible:
            self.VisibilityOn()
        else:
            self.VisibilityOff()

    def VisibilityOn(self):
        """shows the scalar bar"""
        if not self.is_shown:
            self.scalar_bar.VisibilityOn()
            self.scalar_bar.Modified()
            self.is_shown = True

    def VisibilityOff(self):
        """hides the scalar bar"""
        if self.is_shown:
            self.scalar_bar.VisibilityOff()
            self.scalar_bar.Modified()
            self.is_shown = False

    def __init__(self, settings: Settings):
        """creates the scalar bar"""
        self.settings = settings
        self.is_horizontal = False

        self.scalar_bar = vtkScalarBarActor()
        prop = self.scalar_bar.GetTitleTextProperty()
        set_vtk_property_to_unicode(prop, bi_font_file)

        # these don't work...you have to change the font file
        #prop.BoldOn()
        #prop.ItalicOn()
        #prop.ShadowOn()

        self.color_function = vtkColorTransferFunction()
        self.color_function.SetNanColor(1., 1., 1.)

        self.is_shown = True
        self.colormap = 'jet'
        self.colormap_space = ''  # hsv, rgb
        self.is_low_to_high = True
        self.min_value = 10.
        self.max_value = 20.
        self.nlabels = 11
        self.ncolors = 11
        self.data_format = '%i'

        #self.color_function.SetNanColor(0., 0., 0.)
        #self.color_function.SetColorSpaceToLab()
        #self.color_function.SetColorSpaceToRGB()
        #self.scalar_bar.SetDragable(True)
        #self.scalar_bar.SetPickable(True)

        # blue - low
        # red - high
        drange = [10., 20.]
        self.color_function.SetColorSpaceToHSV()
        self.color_function.HSVWrapOff()
        self.color_function.SetRange(*drange)
        self.color_function.AddRGBPoint(drange[0], 0.0, 0.0, 1.0)
        self.color_function.AddRGBPoint(drange[1], 1.0, 0.0, 0.0)

        self.scalar_bar.SetTitle('Title1')

        self.scalar_bar.SetLookupTable(self.color_function)
        #self.scalar_bar.SetNanColor(0., 0., 0.) # RGB color - black
        #self.scalar_bar.SetNanColor(1., 1., 1., 0.) # RGBA color - white

        # old
        #self.scalar_bar.SetHeight(0.9)
        #self.scalar_bar.SetWidth(0.20)  # the width is set first
        #self.scalar_bar.SetPosition(0.77, 0.1)

        # put the scalar bar at the right side
        self.scalar_bar.SetOrientationToVertical()
        width = 0.2
        height = 0.9
        x = 1 - 0.01 - width
        y = (1 - height) / 2.
        self.scalar_bar.SetPosition(x, y)

        # the width is set first
        # after the width is set, this is adjusted
        self.scalar_bar.SetHeight(height)
        self.scalar_bar.SetWidth(width)
        self.scalar_bar.SetPosition(x, y)

        prop_title = vtkTextProperty()
        prop_title.SetFontFamilyToArial()
        #prop_title.ItalicOff()
        prop_title.BoldOn()
        prop_title.ShadowOn()

        prop_label = vtkTextProperty()
        prop_label.BoldOff()
        prop_label.ShadowOn()

        #self.scalar_bar.SetTitleTextProperty(prop_title)
        #self.scalar_bar.SetLabelTextProperty(prop_label)
        self.scalar_bar.SetLabelFormat('%i')

        # allows 0-1 to be nice number when ranging values (gotta pick something)
        self.scalar_bar.SetNumberOfLabels(11)
        self.scalar_bar.SetMaximumNumberOfColors(11)

        #self.scalar_bar.VisibilityOff()  # first load -> scalar bar off
        #self.scalar_bar.ShadowOn()
        #self.scalar_bar.RepositionableOn()
        self.scalar_bar.VisibilityOff()

    def update_color_function(self,
                              min_value, max_value,
                              ncolors: int,
                              colormap: str='jet',
                              colormap_space: str='',
                              is_low_to_high: bool=True):
        """updates the color_function"""
        # jet - HSV :)
        # jet with RGB is red to blue (not a bad colormap, but not jet)

        # viridis and plasma look good as HSV
        # (not sure on the exact difference, but it probably should be
        #  RGB based on the others)
        #
        # viridis - RGB
        # plasma  - RGB
        # magma   - not HSV, RGB
        # inferno - not HSV, RGB
        ncolors = NCOLORS_DEFAULT if ncolors is None else ncolors
        colormap_space, colormap = get_colormap_space(colormap_space, colormap)

        update_color_function = (
            colormap_space != self.colormap_space or
            is_low_to_high != self.is_low_to_high or
            min_value != self.min_value or
            max_value != self.max_value or
            colormap != self.colormap # this iss last, so it's faster
        )
        if not update_color_function:
            return

        self.colormap = colormap
        self.colormap_space = colormap_space
        self.is_low_to_high = is_low_to_high

        self.color_function.RemoveAllPoints()

        if colormap_space == 'rgb':
            self.color_function.SetColorSpaceToRGB()
        elif colormap_space == 'hsv':
            self.color_function.SetColorSpaceToHSV()
        else:  # pragma: no cover
            raise NotImplementedError(colormap_space)

        if colormap == 'jet':
            if is_low_to_high:
                self.color_function.AddRGBPoint(min_value, 0.0, 0.0, 1.0)  # blue
                self.color_function.AddRGBPoint(max_value, 1.0, 0.0, 0.0)  # red
            else:
                self.color_function.AddRGBPoint(min_value, 1.0, 0.0, 0.0)  # red
                self.color_function.AddRGBPoint(max_value, 0.0, 0.0, 1.0)  # blue
        else:
            nvalues = ncolors
            #nvalues = ncolors + 1
            #nvalues = ncolors * 2
            if isinstance(colormap, str):
                colormap_array = colormap_dict[colormap]
                ntotal = len(colormap_array)

                i = np.round(np.linspace(0., ntotal-1, num=nvalues)).astype('int32')
                colormap_array2 = colormap_array[i, :]
            else:
                assert isinstance(colormap[0][0], float), colormap

            vals = np.linspace(min_value, max_value, num=nvalues)
            if is_low_to_high:
                vals = vals[::-1]

            # AddHSVPoint
            ## TODO: requires a vtk object
            #self.color_function.AddRGBPoints(vals, colormap_array)\
            assert len(vals) == len(colormap_array2)
            for val, (red, green, blue) in zip(vals, colormap_array2):
                self.color_function.AddRGBPoint(val, red, green, blue)

    def update_position(self, is_horizontal: bool=True) -> None:
        """updates if the scalar bar is horizontal/vertical"""
        update_position = is_horizontal is not self.is_horizontal
        if not update_position: # pragma: no cover
            return
        if is_horizontal:
            # put the scalar bar at the top
            self.scalar_bar.SetOrientationToHorizontal()
            width = 0.95
            height = 0.15
            x = (1 - width) / 2.
            y = 1 - 0.02 - height
        else:
            # put the scalar bar at the right side
            self.scalar_bar.SetOrientationToVertical()
            width = 0.2
            height = 0.9
            x = 1 - 0.01 - width
            y = (1 - height) / 2.

        self.is_horizontal = is_horizontal
        self.settings.is_horizontal_scalar_bar = is_horizontal
        self.scalar_bar.SetHeight(height)
        self.scalar_bar.SetWidth(width)
        self.scalar_bar.SetPosition(x, y)

    def update_title(self, title: str) -> None:
        """Updates the title.  Pads the text so text isn't huge."""
        nchars = len(title)
        if nchars > 10:
            padding = ''
        else:
            nspaces = (10 - nchars) // 2 + nchars % 2
            padding = nspaces * ' '
        self.scalar_bar.SetTitle('%s%s%s' % (padding, title, padding))

    def update_data_format(self,
                           min_value: float, max_value: float,
                           data_format: str,
                           nlabels=None, ncolors=None):
        """updates the data format and number of values/colors"""
        data_format_display = data_format
        if nlabels is None: # and labelsize is None:
            nvalues = NCOLORS_DEFAULT
            if data_format == '%i':
                data_format_display = '%.0f'
                nvalues = int(max_value - min_value) + 1

                # old
                if nvalues < 7:
                    nvalues = 7
                elif nvalues > 30:
                    nvalues = 11

                # new
                #if 0:
                    #text_prop = self.scalar_bar.GetLabelTextProperty()
                    ##font_size = text_prop.GetFontSize()
                    #nvalues_max = 11
                    #nvalues_min = 7
                    #if nvalues > nvalues_max:
                        #font_size = 4
                        #nvalues = nvalues_max
                        #font_size = text_prop.SetFontSize(font_size)
                        #text_prop.Modified()
                    #elif nvalues < nvalues_min:
                        #nvalues = nvalues_min
                        #font_size = 12
                        #font_size = text_prop.SetFontSize(font_size)
                        #text_prop.Modified()
        else:
            if data_format == '%i':
                data_format_display = '%.0f'
            nvalues = nlabels

        if ncolors is None:
            ncolors = nvalues

        assert data_format_display is not None, f'data_format is invalid = {data_format_display!r}'

        # the code explodes if these are too big
        if nvalues > 100:
            nvalues = 100
        if ncolors > 100:
            ncolors = 100

        if ncolors < 2 and (max_value - min_value) > 0: # data_format == '%i' and
            ncolors = 2

        update_data_format = (
            min_value != self.min_value or
            max_value != self.max_value or
            data_format != self.data_format or
            #nlabels != self.nlabels or
            ncolors != self.ncolors
        )
        if not update_data_format:
            return

        self.scalar_bar.SetLabelFormat(data_format_display)

        #print('ncolors=%s nvalues=%s' % (ncolors, nvalues))
        self.scalar_bar.SetNumberOfLabels(nvalues)
        self.scalar_bar.SetMaximumNumberOfColors(ncolors)

    def update(self, title: str,
               min_value: float,
               max_value: float,
               data_format: str,
               nlabels=None,
               labelsize=None,
               ncolors=None,
               colormap: str='jet',
               colormap_space: str='',  # hsv, rgb
               is_low_to_high: bool=True,
               is_horizontal: bool=True,
               is_shown: bool=True):
        """updates the scalar bar"""
        #self.set_visibility(is_shown)
        settings: Settings = self.settings

        nlabels = None if nlabels == -1 else nlabels
        labelsize = None if labelsize == -1 else labelsize
        ncolors = None if ncolors == -1 else ncolors
        colormap = settings.colormap if colormap == '' else colormap

        self.update_color_function(min_value, max_value, ncolors,
                                   colormap=colormap, colormap_space=colormap_space,
                                   is_low_to_high=is_low_to_high)

        self.update_position(is_horizontal=is_horizontal)

        #if 0:
            #self.color_function.SetRange(min_value, max_value)
            #self.color_function.Update()

            #scalar_range = self.grid.GetScalarRange()
            #print('scalar_range', scalar_range)
            #self.grid_mapper.SetScalarRange(scalar_range)
            #self.grid_mapper.SetScalarRange(min_value, max_value)
            #self.grid_mapper.SetScalarRange(max_value, min_value)
            #self.grid_mapper.Update()

        #self.scalar_bar.SetLookupTable(self.color_function)
        self.update_title(title)
        if data_format is not None:
            # None is RealGridPointForces
            self.update_data_format(min_value, max_value, data_format,
                                    nlabels=nlabels, ncolors=ncolors)

        self.min_value = min_value
        self.max_value = max_value

        self.set_visibility(is_shown)
        self.scalar_bar.Modified()

#def _is_int_result(data_format):
    #if 'i' in data_format:
        #return True
    #return False

def get_colormap_space(colormap_space: str,
                       colormap: str) -> tuple[str, str]:
    if colormap_space == '':
        if colormap in ['jet', 'jet2', 'blend', '', None] or colormap in HSV_MAPS:
            colormap_space = 'hsv'
            colormap = 'jet'
        elif colormap in RGB_MAPS:
            colormap_space = 'rgb'
        else:  # pragma: no cover
            raise NotImplementedError(colormap)
    return colormap_space, colormap
