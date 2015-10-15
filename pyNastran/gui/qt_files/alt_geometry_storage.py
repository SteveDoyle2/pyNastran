from copy import deepcopy

class AltGeometry(object):

    def __init__(self, parent, name, color=None, line_width=1, opacity=0.0,
                 point_size=1, bar_scale=1.0, representation='main', is_visible=True):
        """
        Parameters
        ----------
        line_width : int
            the width of the line for 'surface' and 'main'
        color : [int, int, int]
            the RGB colors
        opacity : float
            0.0 -> solid
            1.0 -> transparent
        point_size : int
            the point size for 'point'
        bar_scale : float
            the scale for the CBAR / CBEAM elements
        representation : str
            main - change with main mesh
            wire - always wireframe
            point - always points
            surface - always surface

        """
        if line_width is None:
            line_width = 1
        if opacity is None:
            opacity = 1.0

        self.parent = parent
        self.name = name
        self._color = None
        if color is not None:
            assert color is not None, color
            self.color = color
        self.line_width = line_width
        self.point_size = point_size
        self._opacity = opacity
        self.bar_scale = bar_scale

        assert isinstance(is_visible, bool), is_visible
        self.is_visible = is_visible

        assert representation in ['main', 'toggle', 'wire', 'point', 'surface'], 'representation=%r' % representation
        self.representation = representation

    def __deepcopy__(self, memo):
        keys = ['name', '_color', 'line_width', 'point_size', '_opacity',
                '_representation', 'is_visible', 'bar_scale']
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k in keys:
            v = self.__dict__[k]
            setattr(result, k, deepcopy(v, memo))
        return result

    @property
    def opacity(self):
        """
        0 -> transparent
        1 -> solid
        """
        assert 0.0 <= self._opacity <= 1.0, self._opacity
        return self._opacity

    @opacity.setter
    def opacity(self, opacity):
        assert 0.0 <= opacity <= 1.0, opacity
        self._opacity = opacity

    @property
    def transparency(self):
        """
        0 -> solid
        1 -> transparent
        """
        assert 0.0 <= self._opacity <= 1.0, self._opacity
        return 1.0 - self._opacity

    @transparency.setter
    def transparency(self, transparency):
        assert 0.0 <= transparency <= 1.0, transparency
        self._opacity = 1.0 - transparency

    @property
    def color(self):
        if self._color is None:
            return (255, 0, 0)  # the default color; red
        return self._color

    @color.setter
    def color(self, color):
        assert len(color) == 3, color
        if isinstance(color[0], int):
            assert isinstance(color[0], int), color[0]
            assert isinstance(color[1], int), color[1]
            assert isinstance(color[2], int), color[2]
            self._color = tuple(color)
        else:
            assert isinstance(color[0], float), color[0]
            assert isinstance(color[1], float), color[1]
            assert isinstance(color[2], float), color[2]
            self._color = (int(color[0] * 255), int(color[1] * 255), int(color[2] * 255))

    @property
    def color_float(self):
        return tuple([i/255. for i in self._color])

    def set_color(self, color, mode='rgb'):
        assert mode == 'rgb', 'mode=%r' % mode
        self.color = color
        assert len(color) == 3, color
        #self.mode = 'rgb'

    @property
    def representation(self):
        """
        * main - main mesh
        * toggle - change with main mesh
        * wire - always wireframe
        * point - always points
        * surface - always surface
        """
        return self._representation

    @representation.setter
    def representation(self, representation):
        assert representation in ['main', 'toggle', 'wire', 'point', 'surface'], 'representation=%r is invalid' % representation
        self._representation = representation
