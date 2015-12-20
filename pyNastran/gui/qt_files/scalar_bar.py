import vtk


class ScalarBar(object):

    def set_visibility(self, is_visible):
        if is_visible:
            self.VisibilityOn()
        else:
            self.VisibilityOff()

    def VisibilityOn(self):
        if not self.is_shown:
            self.scalar_bar.VisibilityOn()
            self.scalar_bar.Modified()
            self.is_shown = True

    def VisibilityOff(self):
        if self.is_shown:
            self.scalar_bar.VisibilityOff()
            self.scalar_bar.Modified()
            self.is_shown = False

    def __init__(self, is_horizontal):
        self.scalar_bar = vtk.vtkScalarBarActor()
        self.color_function = vtk.vtkColorTransferFunction()
        self.is_shown = True
        self.is_horizontal = False

        self.color_function.SetColorSpaceToHSV()
        self.color_function.HSVWrapOff()
        #self.color_function.SetNanColor(1., 1., 1., 0.)
        #self.color_function.SetColorSpaceToLab()
        #self.color_function.SetColorSpaceToRGB()
        #self.scalar_bar.SetDragable(True)
        #self.scalar_bar.SetPickable(True)

        drange = [10., 20.]
        self.color_function.SetRange(*drange)

        # blue - low
        # red - high
        self.color_function.AddRGBPoint(drange[0], 0.0, 0.0, 1.0)
        self.color_function.AddRGBPoint(drange[1], 1.0, 0.0, 0.0)

        self.scalar_bar.SetTitle("Title1")

        self.scalar_bar.SetLookupTable(self.color_function)
        #self.scalar_bar.SetNanColor(0., 0., 0.) # RGB color - black
        #self.scalar_bar.SetNanColor(1., 1., 1., 0.) # RGBA color - white

        # old
        #self.scalar_bar.SetHeight(0.9)
        #self.scalar_bar.SetWidth(0.20)  # the width is set first
        #self.scalar_bar.SetPosition(0.77, 0.1)
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
            self.scalar_bar.SetPosition(x, y)

        # the width is set first
        # after the width is set, this is adjusted
        self.scalar_bar.SetHeight(height)
        self.scalar_bar.SetWidth(width)
        self.scalar_bar.SetPosition(x, y)

        prop_title = vtk.vtkTextProperty()
        prop_title.SetFontFamilyToArial()
        #prop_title.ItalicOff()
        prop_title.BoldOn()
        prop_title.ShadowOn()

        prop_label = vtk.vtkTextProperty()
        prop_label.BoldOff()
        prop_label.ShadowOn()

        #self.scalar_bar.SetTitleTextProperty(prop_title)
        #self.scalar_bar.SetLabelTextProperty(prop_label)
        self.scalar_bar.SetLabelFormat("%i")

        # allows 0-1 to be nice number when ranging values (gotta pick something)
        self.scalar_bar.SetNumberOfLabels(11)
        self.scalar_bar.SetMaximumNumberOfColors(11)

        #self.scalar_bar.VisibilityOff()  # first load -> scalar bar off
        #self.scalar_bar.ShadowOn()
        #self.scalar_bar.RepositionableOn()
        self.scalar_bar.VisibilityOff()


    def update(self, title, min_value, max_value, norm_value,
               data_format, is_blue_to_red=True, is_horizontal=True,
               is_shown=True):
        self.color_function.RemoveAllPoints()

        if is_blue_to_red:
            self.color_function.AddRGBPoint(min_value, 0.0, 0.0, 1.0)  # blue
            self.color_function.AddRGBPoint(max_value, 1.0, 0.0, 0.0)  # red
        else:
            self.color_function.AddRGBPoint(min_value, 1.0, 0.0, 0.0)  # red
            self.color_function.AddRGBPoint(max_value, 0.0, 0.0, 1.0)  # blue

        if is_horizontal:
            # put the scalar bar at the top
            self.is_horizontal = True
            self.scalar_bar.SetOrientationToHorizontal()
            width = 0.95
            height = 0.15
            x = (1 - width) / 2.
            y = 1 - 0.02 - height
        else:
            # put the scalar bar at the right side
            self.is_horizontal = False
            self.scalar_bar.SetOrientationToVertical()
            width = 0.2
            height = 0.9
            x = 1 - 0.01 - width
            y = (1 - height) / 2.
        self.scalar_bar.SetHeight(height)
        self.scalar_bar.SetWidth(width)
        self.scalar_bar.SetPosition(x, y)

        if 0:
            self.color_function.SetRange(min_value, max_value)
            #self.color_function.Update()

            #scalar_range = self.grid.GetScalarRange()
            #print('scalar_range', scalar_range)
            #self.aQuadMapper.SetScalarRange(scalar_range)
            #self.aQuadMapper.SetScalarRange(min_value, max_value)
            #self.aQuadMapper.SetScalarRange(max_value, min_value)
            #self.aQuadMapper.Update()

        #self.scalar_bar.SetLookupTable(self.color_function)
        nchars = len(title)
        if nchars > 10:
            padding = ''
        else:
            nspaces = (10 - nchars) // 2 + nchars % 2
            padding = nspaces * ' '

        self.scalar_bar.SetTitle('%s%s%s' % (padding, title, padding))

        nvalues = 11
        data_format_display = data_format
        if data_format == '%i':
            data_format_display = '%.0f'
            nvalues = int(max_value - min_value) + 1

            # old
            if nvalues < 7:
                nvalues = 7
            elif nvalues > 30:
                nvalues = 11

            # new
            if 0:
                text_prop = self.scalar_bar.GetLabelTextProperty()
                #font_size = text_prop.GetFontSize()
                nvalues_max = 11
                nvalues_min = 7
                if nvalues > nvalues_max:
                    font_size = 4
                    nvalues = nvalues_max
                    font_size = text_prop.SetFontSize(font_size)
                    text_prop.Modified()
                elif nvalues < nvalues_min:
                    nvalues = nvalues_min
                    font_size = 12
                    font_size = text_prop.SetFontSize(font_size)
                    text_prop.Modified()

        assert data_format_display is not None, 'data_format is invalid = %r' % data_format_display
        self.scalar_bar.SetLabelFormat(data_format_display)

        #if title in ['ElementID', 'Eids', 'Region'] and norm_value < 11:
            #nvalues = int(max_value - min_value) + 1
            #print("need to adjust axes...max_value=%s" % max_value)

        #if self.nvalues is not None:
            #if not _is_int_result(data_format):
                ## don't change nvalues for int results
                #nvalues = self.nvalues

        self.scalar_bar.SetNumberOfLabels(nvalues)
        self.scalar_bar.SetMaximumNumberOfColors(nvalues)
        #is_shown = False
        #if is_shown:
            #self.scalar_bar.VisibilityOn()
        #else:
            #self.scalar_bar.VisibilityOff()
        self.set_visibility(is_shown)
        self.scalar_bar.Modified()

def _is_int_result(data_format):
    if 'i' in data_format:
        return True
    return False
