from typing import Union, Optional
from copy import deepcopy
from pyNastran.gui.typing import ColorFloat, ColorInt


class AltGeometry:
    representations = ['main', 'toggle', 'wire', 'point', 'surface',
                       'bar', 'wire+point', 'wire+surf']
    displays = ['Wireframe', 'Surface', 'point', None]

    def __init__(self, parent, name: str,
                 color: Optional[ColorInt]=None,
                 line_width: int=1,
                 opacity: float=0.0,
                 point_size: int=1,
                 bar_scale: float=1.0,
                 representation: str='main',
                 display=None,
                 is_visible: bool=True,
                 is_pickable: bool=False, label_actors=None,
                 visible_in_geometry_properties: bool=True,
                 ):
        """
        Creates an AltGeometry object

        Parameters
        ----------
        line_width : int
            the width of the line for 'surface' and 'main'
        color : ColorInt
            the RGB colors (0-255)
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
            bar - can use bar scale
            toggle - follow the main mesh
            wire+point - is this used???
            wire+surf - two options
        display : str
            only relevant to wire+surf
            the active state of the mesh
            'Surface', 'Wireframe', 'point'
        is_visible : bool; default=True
            is this actor currently visible
        is_pickable : bool; default=False
            can you pick a node/cell on this actor
        label_actors : list[annotation]; None -> []
            stores annotations (e.g., for a control surface)
        visible_in_geometry_properties : bool; default=True
            True: show up in ``Edit Geometry Properties`` menu
            False: don't show up

        """
        representation_map = {
            'main' : None,
            'wire' : 'Wireframe',
            'point' : 'point',
            'surface' : 'Surface',
            'wire+surf' : 'Surface',
            'wire+point' : 'Wireframe',
            'bar' : None,
            'toggle' : None,
        }
        if display is None:
            try:
                display = representation_map[representation]
            except KeyError:
                valid_keys = list(representation_map.keys())
                valid_keys.sort()
                raise RuntimeError('%r is not a valid representation\nvalid=[%s]' % (
                    representation, ', '.join(valid_keys)))

        if line_width is None:
            line_width = 1
        if opacity is None:
            opacity = 1.0
        if label_actors is None:
            label_actors = []

        self.parent = parent
        self.name = name
        self.display = display
        assert display in self.displays, 'dislay=%r displays=%s' % (display, self.displays)

        assert isinstance(name, str), 'name=%r' % name
        assert isinstance(label_actors, list), 'name=%r label_actors=%s' % (name, str(label_actors))
        self._color = None
        if color is not None:
            assert color is not None, color
            self.color = color
        self.line_width = line_width
        self.point_size = point_size
        self._opacity = opacity
        self.bar_scale = bar_scale
        self.label_actors = []

        assert isinstance(is_visible, bool), is_visible
        self.is_visible = is_visible
        self.is_pickable = is_pickable

        if representation not in self.representations:
            msg = 'representation=%r is invalid\nrepresentations=%r' % (
                representation, self.representations)
            raise RuntimeError(msg)
        self.representation = representation

        self.visible_in_geometry_properties = visible_in_geometry_properties

    def __deepcopy__(self, memo):
        """doesn't copy the label_actors to speed things up?"""
        keys = ['name', '_color', 'display', 'line_width', 'point_size', '_opacity',
                '_representation', 'is_visible', 'bar_scale', 'is_pickable',
                'visible_in_geometry_properties']
        cls = self.__class__
        result = cls.__new__(cls)
        idi = id(self)
        memo[idi] = result
        for key in keys:
            value = self.__dict__[key]
            setattr(result, key, deepcopy(value, memo))

        # gotta set something or the repr will crash
        #result.label_actors = [] #= memo['label_actors']
        result.label_actors = []
        return result

    @property
    def opacity(self) -> float:
        """
        0 -> transparent
        1 -> solid

        """
        assert 0.0 <= self._opacity <= 1.0, self._opacity
        return self._opacity

    @opacity.setter
    def opacity(self, opacity: float) -> None:
        assert 0.0 <= opacity <= 1.0, opacity
        self._opacity = opacity

    @property
    def transparency(self) -> float:
        """
        0 -> solid
        1 -> transparent

        """
        assert 0.0 <= self._opacity <= 1.0, self._opacity
        return 1.0 - self._opacity

    @transparency.setter
    def transparency(self, transparency: float) -> None:
        assert 0.0 <= transparency <= 1.0, transparency
        self._opacity = 1.0 - transparency

    @property
    def color(self) -> ColorInt:
        if self._color is None:
            return (255, 0, 0)  # the default color; red
        return self._color

    @color.setter
    def color(self, color: Union[ColorInt, ColorFloat]) -> None:
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
    def color_float(self) -> ColorFloat:
        return tuple([i/255. for i in self._color])
    @property
    def color_int(self) -> ColorInt:
        return self._color

    def set_color(self, color: Union[ColorInt, ColorFloat],
                  mode: str='rgb') -> None:
        assert mode == 'rgb', 'mode=%r' % mode
        self.color = color
        assert len(color) == 3, color
        #self.mode = 'rgb'

    @property
    def representation(self) -> str:
        """
        Gets the representation

        * main - main mesh
        * toggle - change with main mesh
        * wire - always wireframe
        * point - always points
        * surface - always surface
        * bar - this can use bar scale
        * wire+point - point (vertex) and wireframe allowed
        * wire+surf - the user can switch between surface and wireframe as a selection

        """
        return self._representation

    @representation.setter
    def representation(self, representation: str) -> None:
        """Sets the representation"""
        if representation not in self.representations:
            msg = 'representation=%r is invalid\nrepresentations=%r' % (
                representation, self.representations)
            raise RuntimeError(msg)
        self._representation = representation

    def __repr__(self):
        msg = ('AltGeometry(%r, color=%s, line_width=%s, opacity=%s,\n'
              ' point_size=%s, bar_scale=%s, representation=%r, display=%r, is_visible=%s,\n'
              'is_pickable=%s, label_actors=%s)' % (
                  self.name, str(self.color), self.line_width, self.opacity, self.point_size,
                  self.bar_scale, self.representation, self.display, self.is_visible,
                  self.is_pickable, self.label_actors))
        return msg
